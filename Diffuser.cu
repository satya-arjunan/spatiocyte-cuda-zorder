//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
//
//        This file is part of the Spatiocyte package
//
//        Copyright (C) 2006-2009 Keio University
//        Copyright (C) 2010-2014 RIKEN
//
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
//
//
// Spatiocyte is free software; you can redistribute it and/or
// modify it under the terms of the GNU General Public
// License as published by the Free Software Foundation; either
// version 2 of the License, or (at your option) any later version.
// 
// Spatiocyte is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
// See the GNU General Public License for more details.
// 
// You should have received a copy of the GNU General Public
// License along with Spatiocyte -- see the file COPYING.
// If not, write to the Free Software Foundation, Inc.,
// 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
// 
//END_HEADER
//
// written by Satya Arjunan <satya.arjunan@gmail.com>
//

#include <time.h>
#include <thrust/execution_policy.h>
#include <thrust/random.h>
#include <thrust/system/cuda/detail/bulk/bulk.hpp>
#include <curand_kernel.h>
#include <Diffuser.hpp>
#include <Compartment.hpp>
#include <Model.hpp>
#include <Reaction.hpp>
#include <random>

Diffuser::Diffuser(const double D, Species& species):
  D_(D),
  species_(species),
  compartment_(species_.get_compartment()),
  mols_(species_.get_mols()),
  blocks_(compartment_.get_model().get_blocks()),
  voxels_(compartment_.get_lattice().get_voxels()),
  offsets_(compartment_.get_offsets()),
  species_id_(species_.get_id()),
  vac_id_(species_.get_vac_id()),
  null_id_(species_.get_model().get_null_id()),
  seed_(0) {
}

void Diffuser::initialize() {
  Model& model(species_.get_model());
  stride_ = model.get_stride();
  id_stride_ = species_id_*stride_;
  is_reactive_.resize(model.get_species().size(), false);
  reactions_.resize(model.get_species().size(), NULL);
  substrate_mols_.resize(model.get_species().size(), NULL);
  product_mols_.resize(model.get_species().size(), NULL);
  reacteds_.resize(mols_.size()+1, 0);
 
  std::vector<Reaction*>& reactions(species_.get_reactions());
  for(unsigned i(0); i != reactions.size(); ++i) {
    std::vector<Species*>& substrates(reactions[i]->get_substrates());
    for(unsigned j(0); j != substrates.size(); ++j) {
      voxel_t reactant_id(substrates[j]->get_id());
      if(reactant_id != species_id_) {
        reactions_[reactant_id] = reactions[i];
        is_reactive_[reactant_id] = true;
        substrate_mols_[reactant_id] = thrust::raw_pointer_cast(substrates[j]->get_mols().data());
        product_mols_[reactant_id] = thrust::raw_pointer_cast(reactions[i]->get_products()[0]->get_mols().data());
      } 
    } 
  } 
}

double Diffuser::get_D() const {
  return D_;
}

/* kernel<<<3, 5>>>()
   <-block0-><-block1-><-block2->
   |0|1|2|3|4|0|1|2|3|4|0|1|2|3|4|
   gridDim = number of blocks in a grid = 3
   blockDim = number of threads per block = 5
   blockIdx = index of the block = [0,1,2]
   threadIdx = index of the thread in a block = [0,1,2,3,4]
   980 GTX: multiProcessorCount = 16
*/


__global__
void concurrent_walk(
    const unsigned mol_size_,
    const unsigned seed_,
    const voxel_t stride_,
    const voxel_t id_stride_,
    const voxel_t vac_id_,
    const voxel_t null_id_,
    umol_t* mols_,
    voxel_t* voxels_) {
  __shared__ int offsets_[48];
  if(threadIdx.x == 0) {
    //col=even, layer=even
    offsets_[0] = -1;
    offsets_[1] = 1;
    offsets_[2] = -NUM_ROW-1;
    offsets_[3] = -NUM_ROW;
    offsets_[4] = NUM_ROW-1;
    offsets_[5] = NUM_ROW;
    offsets_[6] = -NUM_COLROW-NUM_ROW;
    offsets_[7] = -NUM_COLROW-1;
    offsets_[8] = -NUM_COLROW;
    offsets_[9] = NUM_COLROW-NUM_ROW;
    offsets_[10] = NUM_COLROW-1;
    offsets_[11] = NUM_COLROW;

    //col=even, layer=odd +24 = %layer*24
    offsets_[24] = -1;
    offsets_[25] = 1;
    offsets_[26] = -NUM_ROW;
    offsets_[27] = -NUM_ROW+1;
    offsets_[28] = NUM_ROW;
    offsets_[29] = NUM_ROW+1;
    offsets_[30] = -NUM_COLROW;
    offsets_[31] = -NUM_COLROW+1;
    offsets_[32] = -NUM_COLROW+NUM_ROW;
    offsets_[33] = NUM_COLROW;
    offsets_[34] = NUM_COLROW+1;
    offsets_[35] = NUM_COLROW+NUM_ROW;

    //col=odd, layer=even +12 = %col*12
    offsets_[12] = -1;
    offsets_[13] = 1;
    offsets_[14] = -NUM_ROW;
    offsets_[15] = -NUM_ROW+1;
    offsets_[16] = NUM_ROW;
    offsets_[17] = NUM_ROW+1;
    offsets_[18] = -NUM_COLROW-NUM_ROW;
    offsets_[19] = -NUM_COLROW;
    offsets_[20] = -NUM_COLROW+1;
    offsets_[21] = NUM_COLROW-NUM_ROW;
    offsets_[22] = NUM_COLROW;
    offsets_[23] = NUM_COLROW+1;

    //col=odd, layer=odd +36 = %col*12 + %layer*24
    offsets_[36] = -1;
    offsets_[37] = 1;
    offsets_[38] = -NUM_ROW-1;
    offsets_[39] = -NUM_ROW;
    offsets_[40] = NUM_ROW-1;
    offsets_[41] = NUM_ROW;
    offsets_[42] = -NUM_COLROW-1;
    offsets_[43] = -NUM_COLROW; //a
    offsets_[44] = -NUM_COLROW+NUM_ROW;
    offsets_[45] = NUM_COLROW-1;
    offsets_[46] = NUM_COLROW;
    offsets_[47] = NUM_COLROW+NUM_ROW;
  }
  __syncthreads();
  //index is the unique global thread id (size: total_threads)
  unsigned index(blockIdx.x*blockDim.x + threadIdx.x);
  const unsigned total_threads(blockDim.x*gridDim.x);
  curandState local_state = curand_states[blockIdx.x][threadIdx.x];
  while(index < mol_size_) {
    const uint32_t rand32(curand(&local_state));
    uint16_t rand16((uint16_t)(rand32 & 0x0000FFFFuL));
    uint32_t rand(((uint32_t)rand16*12) >> 16);
    umol_t vdx(mols_[index]);
    bool odd_lay((vdx/NUM_COLROW)&1);
    bool odd_col((vdx%NUM_COLROW/NUM_ROW)&1);
    mol2_t val(mol2_t(vdx)+offsets_[rand+(24&(-odd_lay))+(12&(-odd_col))]);
    //Atomically put the current molecule id, index+id_stride_ at the target
    //voxel if it is vacant: 
    voxel_t tar_mol_id(atomicCAS(voxels_+val, vac_id_, index+id_stride_));
    //If not occupied, finalize walk:
    if(tar_mol_id == vac_id_) {
      voxels_[vdx] = vac_id_;
      mols_[index] = val;
    }
    index += total_threads;
    if(index < mol_size_) {
      rand16 = (uint16_t)(rand32 >> 16);
      rand = ((uint32_t)rand16*12) >> 16;
      vdx = mols_[index];
      odd_lay = (vdx/NUM_COLROW)&1;
      odd_col = (vdx%NUM_COLROW/NUM_ROW)&1;
      val = mol2_t(vdx)+offsets_[rand+(24&(-odd_lay))+(12&(-odd_col))];
      tar_mol_id = (atomicCAS(voxels_+val, vac_id_, index+id_stride_));
      //If not occupied, finalize walk:
      if(tar_mol_id == vac_id_) {
        voxels_[vdx] = vac_id_;
        mols_[index] = val;
      }
      index += total_threads;
    }
  }
  curand_states[blockIdx.x][threadIdx.x] = local_state;
}

void Diffuser::walk() {
  const size_t size(mols_.size());
  concurrent_walk<<<64, 256>>>(
      size,
      seed_,
      stride_,
      id_stride_,
      vac_id_,
      null_id_,
      thrust::raw_pointer_cast(&mols_[0]),
      thrust::raw_pointer_cast(&voxels_[0]));
  seed_ += size;
  cudaThreadSynchronize();
}

/*
//without voxels: 18.2 BUPS
__global__
void concurrent_walk(
    const unsigned mol_size_,
    const unsigned seed_,
    const voxel_t stride_,
    const voxel_t id_stride_,
    const voxel_t vac_id_,
    const voxel_t null_id_,
    umol_t* mols_,
    voxel_t* voxels_) {
  __shared__ int offsets_[48];
  if(threadIdx.x == 0) {
    //col=even, layer=even
    offsets_[0] = -1;
    offsets_[1] = 1;
    offsets_[2] = -NUM_ROW-1;
    offsets_[3] = -NUM_ROW;
    offsets_[4] = NUM_ROW-1;
    offsets_[5] = NUM_ROW;
    offsets_[6] = -NUM_COLROW-NUM_ROW;
    offsets_[7] = -NUM_COLROW-1;
    offsets_[8] = -NUM_COLROW;
    offsets_[9] = NUM_COLROW-NUM_ROW;
    offsets_[10] = NUM_COLROW-1;
    offsets_[11] = NUM_COLROW;

    //col=even, layer=odd +24 = %layer*24
    offsets_[24] = -1;
    offsets_[25] = 1;
    offsets_[26] = -NUM_ROW;
    offsets_[27] = -NUM_ROW+1;
    offsets_[28] = NUM_ROW;
    offsets_[29] = NUM_ROW+1;
    offsets_[30] = -NUM_COLROW;
    offsets_[31] = -NUM_COLROW+1;
    offsets_[32] = -NUM_COLROW+NUM_ROW;
    offsets_[33] = NUM_COLROW;
    offsets_[34] = NUM_COLROW+1;
    offsets_[35] = NUM_COLROW+NUM_ROW;

    //col=odd, layer=even +12 = %col*12
    offsets_[12] = -1;
    offsets_[13] = 1;
    offsets_[14] = -NUM_ROW;
    offsets_[15] = -NUM_ROW+1;
    offsets_[16] = NUM_ROW;
    offsets_[17] = NUM_ROW+1;
    offsets_[18] = -NUM_COLROW-NUM_ROW;
    offsets_[19] = -NUM_COLROW;
    offsets_[20] = -NUM_COLROW+1;
    offsets_[21] = NUM_COLROW-NUM_ROW;
    offsets_[22] = NUM_COLROW;
    offsets_[23] = NUM_COLROW+1;

    //col=odd, layer=odd +36 = %col*12 + %layer*24
    offsets_[36] = -1;
    offsets_[37] = 1;
    offsets_[38] = -NUM_ROW-1;
    offsets_[39] = -NUM_ROW;
    offsets_[40] = NUM_ROW-1;
    offsets_[41] = NUM_ROW;
    offsets_[42] = -NUM_COLROW-1;
    offsets_[43] = -NUM_COLROW; //a
    offsets_[44] = -NUM_COLROW+NUM_ROW;
    offsets_[45] = NUM_COLROW-1;
    offsets_[46] = NUM_COLROW;
    offsets_[47] = NUM_COLROW+NUM_ROW;
  }
  __syncthreads();
  //index is the unique global thread id (size: total_threads)
  unsigned index(blockIdx.x*blockDim.x + threadIdx.x);
  const unsigned total_threads(blockDim.x*gridDim.x);
  curandState local_state = curand_states[blockIdx.x][threadIdx.x];
  while(index < mol_size_) {
    const uint32_t rand32(curand(&local_state));
    uint16_t rand16((uint16_t)(rand32 & 0x0000FFFFuL));
    uint32_t rand(((uint32_t)rand16*12) >> 16);
    umol_t vdx(mols_[index]);
    bool odd_lay((vdx/NUM_COLROW)&1);
    bool odd_col((vdx%NUM_COLROW/NUM_ROW)&1);
    mol2_t val(mol2_t(vdx)+offsets_[rand+(24&(-odd_lay))+(12&(-odd_col))]);
    //Atomically put the current molecule id, index+id_stride_ at the target
    //voxel if it is vacant: 
    //voxel_t tar_mol_id(atomicCAS(voxels_+val, vac_id_, index+id_stride_));
    //If not occupied, finalize walk:
    //if(tar_mol_id == vac_id_) {
      //voxels_[vdx] = vac_id_;
      mols_[index] = val;
    //}
    index += total_threads;
    if(index < mol_size_) {
      rand16 = (uint16_t)(rand32 >> 16);
      rand = ((uint32_t)rand16*12) >> 16;
      vdx = mols_[index];
      odd_lay = (vdx/NUM_COLROW)&1;
      odd_col = (vdx%NUM_COLROW/NUM_ROW)&1;
      val = mol2_t(vdx)+offsets_[rand+(24&(-odd_lay))+(12&(-odd_col))];
      //tar_mol_id = (atomicCAS(voxels_+val, vac_id_, index+id_stride_));
      //If not occupied, finalize walk:
      //if(tar_mol_id == vac_id_) {
        //voxels_[vdx] = vac_id_;
        mols_[index] = val;
      //}
      index += total_threads;
    }
  }
  curand_states[blockIdx.x][threadIdx.x] = local_state;
}

void Diffuser::walk() {
  const size_t size(mols_.size());
  concurrent_walk<<<64, 256>>>(
      size,
      seed_,
      stride_,
      id_stride_,
      vac_id_,
      null_id_,
      thrust::raw_pointer_cast(&mols_[0]),
      thrust::raw_pointer_cast(&voxels_[0]));
  seed_ += size;
  cudaThreadSynchronize();
}
*/

/*
//with shared offsets and local curand: 0.472 BUPS
__global__
void concurrent_walk(
    const unsigned mol_size_,
    const unsigned seed_,
    const voxel_t stride_,
    const voxel_t id_stride_,
    const voxel_t vac_id_,
    const voxel_t null_id_,
    umol_t* mols_,
    voxel_t* voxels_) {
  __shared__ int offsets_[48];
  if(threadIdx.x == 0) {
    //col=even, layer=even
    offsets_[0] = -1;
    offsets_[1] = 1;
    offsets_[2] = -NUM_ROW-1;
    offsets_[3] = -NUM_ROW;
    offsets_[4] = NUM_ROW-1;
    offsets_[5] = NUM_ROW;
    offsets_[6] = -NUM_COLROW-NUM_ROW;
    offsets_[7] = -NUM_COLROW-1;
    offsets_[8] = -NUM_COLROW;
    offsets_[9] = NUM_COLROW-NUM_ROW;
    offsets_[10] = NUM_COLROW-1;
    offsets_[11] = NUM_COLROW;

    //col=even, layer=odd +24 = %layer*24
    offsets_[24] = -1;
    offsets_[25] = 1;
    offsets_[26] = -NUM_ROW;
    offsets_[27] = -NUM_ROW+1;
    offsets_[28] = NUM_ROW;
    offsets_[29] = NUM_ROW+1;
    offsets_[30] = -NUM_COLROW;
    offsets_[31] = -NUM_COLROW+1;
    offsets_[32] = -NUM_COLROW+NUM_ROW;
    offsets_[33] = NUM_COLROW;
    offsets_[34] = NUM_COLROW+1;
    offsets_[35] = NUM_COLROW+NUM_ROW;

    //col=odd, layer=even +12 = %col*12
    offsets_[12] = -1;
    offsets_[13] = 1;
    offsets_[14] = -NUM_ROW;
    offsets_[15] = -NUM_ROW+1;
    offsets_[16] = NUM_ROW;
    offsets_[17] = NUM_ROW+1;
    offsets_[18] = -NUM_COLROW-NUM_ROW;
    offsets_[19] = -NUM_COLROW;
    offsets_[20] = -NUM_COLROW+1;
    offsets_[21] = NUM_COLROW-NUM_ROW;
    offsets_[22] = NUM_COLROW;
    offsets_[23] = NUM_COLROW+1;

    //col=odd, layer=odd +36 = %col*12 + %layer*24
    offsets_[36] = -1;
    offsets_[37] = 1;
    offsets_[38] = -NUM_ROW-1;
    offsets_[39] = -NUM_ROW;
    offsets_[40] = NUM_ROW-1;
    offsets_[41] = NUM_ROW;
    offsets_[42] = -NUM_COLROW-1;
    offsets_[43] = -NUM_COLROW; //a
    offsets_[44] = -NUM_COLROW+NUM_ROW;
    offsets_[45] = NUM_COLROW-1;
    offsets_[46] = NUM_COLROW;
    offsets_[47] = NUM_COLROW+NUM_ROW;
  }
  __syncthreads();
  //index is the unique global thread id (size: total_threads)
  unsigned index(blockIdx.x*blockDim.x + threadIdx.x);
  const unsigned total_threads(blockDim.x*gridDim.x);
  curandState local_state = curand_states[blockIdx.x][threadIdx.x];
  while(index < mol_size_) {
    const uint32_t rand32(curand(&local_state));
    uint16_t rand16((uint16_t)(rand32 & 0x0000FFFFuL));
    uint32_t rand(((uint32_t)rand16*12) >> 16);
    umol_t vdx(mols_[index]);
    bool odd_lay((vdx/NUM_COLROW)&1);
    bool odd_col((vdx%NUM_COLROW/NUM_ROW)&1);
    mol2_t val(mol2_t(vdx)+offsets_[rand+(24&(-odd_lay))+(12&(-odd_col))]);
    //Atomically put the current molecule id, index+id_stride_ at the target
    //voxel if it is vacant: 
    voxel_t tar_mol_id(atomicCAS(voxels_+val, vac_id_, index+id_stride_));
    //If not occupied, finalize walk:
    if(tar_mol_id == vac_id_) {
      voxels_[vdx] = vac_id_;
      mols_[index] = val;
    }
    index += total_threads;
    if(index < mol_size_) {
      rand16 = (uint16_t)(rand32 >> 16);
      rand = ((uint32_t)rand16*12) >> 16;
      vdx = mols_[index];
      odd_lay = (vdx/NUM_COLROW)&1;
      odd_col = (vdx%NUM_COLROW/NUM_ROW)&1;
      val = mol2_t(vdx)+offsets_[rand+(24&(-odd_lay))+(12&(-odd_col))];
      tar_mol_id = (atomicCAS(voxels_+val, vac_id_, index+id_stride_));
      //If not occupied, finalize walk:
      if(tar_mol_id == vac_id_) {
        voxels_[vdx] = vac_id_;
        mols_[index] = val;
      }
      index += total_threads;
    }
  }
  curand_states[blockIdx.x][threadIdx.x] = local_state;
}

void Diffuser::walk() {
  const size_t size(mols_.size());
  concurrent_walk<<<64, 256>>>(
      size,
      seed_,
      stride_,
      id_stride_,
      vac_id_,
      null_id_,
      thrust::raw_pointer_cast(&mols_[0]),
      thrust::raw_pointer_cast(&voxels_[0]));
  seed_ += size;
  cudaThreadSynchronize();
}
*/

/*
//with shared offsets: 0.475 BUPS
__global__
void concurrent_walk(
    const unsigned mol_size_,
    const unsigned seed_,
    const voxel_t stride_,
    const voxel_t id_stride_,
    const voxel_t vac_id_,
    const voxel_t null_id_,
    umol_t* mols_,
    voxel_t* voxels_) {
  __shared__ int offsets_[48];
  if(threadIdx.x == 0) {
    //col=even, layer=even
    offsets_[0] = -1;
    offsets_[1] = 1;
    offsets_[2] = -NUM_ROW-1;
    offsets_[3] = -NUM_ROW;
    offsets_[4] = NUM_ROW-1;
    offsets_[5] = NUM_ROW;
    offsets_[6] = -NUM_COLROW-NUM_ROW;
    offsets_[7] = -NUM_COLROW-1;
    offsets_[8] = -NUM_COLROW;
    offsets_[9] = NUM_COLROW-NUM_ROW;
    offsets_[10] = NUM_COLROW-1;
    offsets_[11] = NUM_COLROW;

    //col=even, layer=odd +24 = %layer*24
    offsets_[24] = -1;
    offsets_[25] = 1;
    offsets_[26] = -NUM_ROW;
    offsets_[27] = -NUM_ROW+1;
    offsets_[28] = NUM_ROW;
    offsets_[29] = NUM_ROW+1;
    offsets_[30] = -NUM_COLROW;
    offsets_[31] = -NUM_COLROW+1;
    offsets_[32] = -NUM_COLROW+NUM_ROW;
    offsets_[33] = NUM_COLROW;
    offsets_[34] = NUM_COLROW+1;
    offsets_[35] = NUM_COLROW+NUM_ROW;

    //col=odd, layer=even +12 = %col*12
    offsets_[12] = -1;
    offsets_[13] = 1;
    offsets_[14] = -NUM_ROW;
    offsets_[15] = -NUM_ROW+1;
    offsets_[16] = NUM_ROW;
    offsets_[17] = NUM_ROW+1;
    offsets_[18] = -NUM_COLROW-NUM_ROW;
    offsets_[19] = -NUM_COLROW;
    offsets_[20] = -NUM_COLROW+1;
    offsets_[21] = NUM_COLROW-NUM_ROW;
    offsets_[22] = NUM_COLROW;
    offsets_[23] = NUM_COLROW+1;

    //col=odd, layer=odd +36 = %col*12 + %layer*24
    offsets_[36] = -1;
    offsets_[37] = 1;
    offsets_[38] = -NUM_ROW-1;
    offsets_[39] = -NUM_ROW;
    offsets_[40] = NUM_ROW-1;
    offsets_[41] = NUM_ROW;
    offsets_[42] = -NUM_COLROW-1;
    offsets_[43] = -NUM_COLROW; //a
    offsets_[44] = -NUM_COLROW+NUM_ROW;
    offsets_[45] = NUM_COLROW-1;
    offsets_[46] = NUM_COLROW;
    offsets_[47] = NUM_COLROW+NUM_ROW;
  }
  __syncthreads();
  //index is the unique global thread id (size: total_threads)
  unsigned index(blockIdx.x*blockDim.x + threadIdx.x);
  const unsigned total_threads(blockDim.x*gridDim.x);
  while(index < mol_size_) {
    const umol_t vdx(mols_[index]);
    thrust::default_random_engine rng;
    rng.discard(seed_+index);
    thrust::uniform_int_distribution<unsigned> u(0, 11);
    const unsigned rand(u(rng));
    const bool odd_lay((vdx/NUM_COLROW)&1);
    const bool odd_col((vdx%NUM_COLROW/NUM_ROW)&1);
    mol2_t val(mol2_t(vdx)+offsets_[rand+(24&(-odd_lay))+(12&(-odd_col))]);
    //Atomically put the current molecule id, index+id_stride_ at the target
    //voxel if it is vacant: 
    const voxel_t tar_mol_id(atomicCAS(voxels_+val, vac_id_, index+id_stride_));
    //If not occupied, finalize walk:
    if(tar_mol_id == vac_id_) {
      voxels_[vdx] = vac_id_;
      mols_[index] = val;
    }
    //Do nothing, stay at original position
    index += total_threads;
  }
  //__syncthreads();
}

void Diffuser::walk() {
  const size_t size(mols_.size());
  concurrent_walk<<<32, 512>>>(
      size,
      seed_,
      stride_,
      id_stride_,
      vac_id_,
      null_id_,
      thrust::raw_pointer_cast(&mols_[0]),
      thrust::raw_pointer_cast(&voxels_[0]));
  seed_ += size;
  cudaThreadSynchronize();
}
   */

/*
// with thread sync: 0.48 BUPS
__global__
void concurrent_walk(
    const unsigned mol_size_,
    const unsigned seed_,
    const voxel_t stride_,
    const voxel_t id_stride_,
    const voxel_t vac_id_,
    const voxel_t null_id_,
    umol_t* mols_,
    const mol_t* offsets_,
    voxel_t* voxels_) {
  //index is the unique global thread id (size: total_threads)
  unsigned index(blockIdx.x*blockDim.x + threadIdx.x);
  const unsigned total_threads(blockDim.x*gridDim.x);
  while(index < mol_size_) {
    const umol_t vdx(mols_[index]);
    thrust::default_random_engine rng;
    rng.discard(seed_+index);
    thrust::uniform_int_distribution<unsigned> u(0, 11);
    const unsigned rand(u(rng));
    const bool odd_lay((vdx/NUM_COLROW)&1);
    const bool odd_col((vdx%NUM_COLROW/NUM_ROW)&1);
    mol2_t val(mol2_t(vdx)+offsets_[rand+(24&(-odd_lay))+(12&(-odd_col))]);
    //Atomically put the current molecule id, index+id_stride_ at the target
    //voxel if it is vacant: 
    const voxel_t tar_mol_id(atomicCAS(voxels_+val, vac_id_, index+id_stride_));
    //If not occupied, finalize walk:
    if(tar_mol_id == vac_id_) {
      voxels_[vdx] = vac_id_;
      mols_[index] = val;
    }
    //Do nothing, stay at original position
    index += total_threads;
  }
  //__syncthreads();
}

void Diffuser::walk() {
  const size_t size(mols_.size());
  concurrent_walk<<<32, 512>>>(
      size,
      seed_,
      stride_,
      id_stride_,
      vac_id_,
      null_id_,
      thrust::raw_pointer_cast(&mols_[0]),
      thrust::raw_pointer_cast(&offsets_[0]),
      thrust::raw_pointer_cast(&voxels_[0]));
  seed_ += size;
  cudaThreadSynchronize();
}
*/


/* without offsets: 43.5 s
__device__
unsigned get_tar(
    const unsigned vdx,
    const unsigned nrand) {
  const bool odd_col((vdx%NUM_COLROW/NUM_ROW)&1);
  const bool odd_lay((vdx/NUM_COLROW)&1);
  switch(nrand)
    {
    case 1:
      return vdx+1;
    case 2:
      return vdx+(odd_col^odd_lay)-NUM_ROW-1 ;
    case 3:
      return vdx+(odd_col^odd_lay)-NUM_ROW;
    case 4:
      return vdx+(odd_col^odd_lay)+NUM_ROW-1;
    case 5:
      return vdx+(odd_col^odd_lay)+NUM_ROW;
    case 6:
      return vdx+NUM_ROW*(odd_lay-NUM_COL-1)-(odd_col&odd_lay);
    case 7:
      return vdx+!odd_col*(odd_lay-!odd_lay)-NUM_COLROW;
    case 8:
      return vdx+NUM_ROW*(odd_lay-NUM_COL)+(odd_col&!odd_lay);
    case 9:
      return vdx+NUM_ROW*(NUM_COL-!odd_lay)-(odd_col&odd_lay);
    case 10:
      return vdx+NUM_COLROW+!odd_col*(odd_lay-!odd_lay);
    case 11:
      return vdx+NUM_ROW*(NUM_COL+odd_lay)+(odd_col&!odd_lay);
    }
  return vdx-1;
}

__global__
void concurrent_walk(
    const unsigned mol_size_,
    const voxel_t stride_,
    const voxel_t id_stride_,
    const voxel_t vac_id_,
    const voxel_t null_id_,
    umol_t* mols_,
    voxel_t* voxels_) {
  //index is the unique global thread id (size: total_threads)
  unsigned index(blockIdx.x*blockDim.x + threadIdx.x);
  const unsigned total_threads(blockDim.x*gridDim.x);
  curandState local_state = curand_states[blockIdx.x][threadIdx.x];
  while(index < mol_size_) {
    const umol_t vdx(mols_[index]);
    float ranf(curand_uniform(&local_state)*11.999999);
    const unsigned rand((unsigned)truncf(ranf));
    mol2_t val(get_tar(vdx, rand));
    //Atomically put the current molecule id, index+id_stride_ at the target
    //voxel if it is vacant: 
    const voxel_t tar_mol_id(atomicCAS(voxels_+val, vac_id_, index+id_stride_));
    //If not occupied, finalize walk:
    if(tar_mol_id == vac_id_) {
      voxels_[vdx] = vac_id_;
      mols_[index] = val;
    }
    //Do nothing, stay at original position
    index += total_threads;
  }
  curand_states[blockIdx.x][threadIdx.x] = local_state;
}

void Diffuser::walk() {
  const size_t size(mols_.size());
  concurrent_walk<<<blocks_, 256>>>(
      size,
      stride_,
      id_stride_,
      vac_id_,
      null_id_,
      thrust::raw_pointer_cast(&mols_[0]),
      thrust::raw_pointer_cast(&voxels_[0]));
}
*/

/* with curand persistent local states: 42.8 s
__global__
void concurrent_walk(
    const unsigned mol_size_,
    const voxel_t stride_,
    const voxel_t id_stride_,
    const voxel_t vac_id_,
    const voxel_t null_id_,
    umol_t* mols_,
    const mol_t* offsets_,
    voxel_t* voxels_) {
  //index is the unique global thread id (size: total_threads)
  unsigned index(blockIdx.x*blockDim.x + threadIdx.x);
  const unsigned total_threads(blockDim.x*gridDim.x);
  curandState local_state = curand_states[blockIdx.x][threadIdx.x];
  while(index < mol_size_) {
    const umol_t vdx(mols_[index]);
    float ranf(curand_uniform(&local_state)*11.999999);
    const unsigned rand((unsigned)truncf(ranf));
    const bool odd_lay((vdx/NUM_COLROW)&1);
    const bool odd_col((vdx%NUM_COLROW/NUM_ROW)&1);
    mol2_t val(mol2_t(vdx)+offsets_[rand+(24&(-odd_lay))+(12&(-odd_col))]);
    //Atomically put the current molecule id, index+id_stride_ at the target
    //voxel if it is vacant: 
    const voxel_t tar_mol_id(atomicCAS(voxels_+val, vac_id_, index+id_stride_));
    //If not occupied, finalize walk:
    if(tar_mol_id == vac_id_) {
      voxels_[vdx] = vac_id_;
      mols_[index] = val;
    }
    //Do nothing, stay at original position
    index += total_threads;
  }
  curand_states[blockIdx.x][threadIdx.x] = local_state;
}

void Diffuser::walk() {
  const size_t size(mols_.size());
  concurrent_walk<<<blocks_, 256>>>(
      size,
      stride_,
      id_stride_,
      vac_id_,
      null_id_,
      thrust::raw_pointer_cast(&mols_[0]),
      thrust::raw_pointer_cast(&offsets_[0]),
      thrust::raw_pointer_cast(&voxels_[0]));
}
*/

/* with curand: 42.5 s
__global__
void concurrent_walk(
    const unsigned mol_size_,
    const unsigned seed_,
    const voxel_t stride_,
    const voxel_t id_stride_,
    const voxel_t vac_id_,
    const voxel_t null_id_,
    umol_t* mols_,
    const mol_t* offsets_,
    voxel_t* voxels_,
    const float* randoms_) {
  //index is the unique global thread id (size: total_threads)
  unsigned index(blockIdx.x*blockDim.x + threadIdx.x);
  const unsigned total_threads(blockDim.x*gridDim.x);
  while(index < mol_size_) {
    const umol_t vdx(mols_[index]);
    float ranf(randoms_[index]*11.999999);
    const unsigned rand((unsigned)truncf(ranf));
    const bool odd_lay((vdx/NUM_COLROW)&1);
    const bool odd_col((vdx%NUM_COLROW/NUM_ROW)&1);
    mol2_t val(mol2_t(vdx)+offsets_[rand+(24&(-odd_lay))+(12&(-odd_col))]);
    //Atomically put the current molecule id, index+id_stride_ at the target
    //voxel if it is vacant: 
    const voxel_t tar_mol_id(atomicCAS(voxels_+val, vac_id_, index+id_stride_));
    //If not occupied, finalize walk:
    if(tar_mol_id == vac_id_) {
      voxels_[vdx] = vac_id_;
      mols_[index] = val;
    }
    //Do nothing, stay at original position
    index += total_threads;
  }
  //__syncthreads();
}

void Diffuser::walk() {
  const size_t size(mols_.size());
  if(randoms_counter_ > compartment_.get_model().get_randoms_size()-size) {
    compartment_.get_model().generate_randoms();
    randoms_counter_ = 0;
  }
  concurrent_walk<<<blocks_, 512>>>(
      size,
      seed_,
      stride_,
      id_stride_,
      vac_id_,
      null_id_,
      thrust::raw_pointer_cast(&mols_[0]),
      thrust::raw_pointer_cast(&offsets_[0]),
      thrust::raw_pointer_cast(&voxels_[0]),
      thrust::raw_pointer_cast(&randoms_[randoms_counter_]));
  randoms_counter_ += size;
  //barrier cudaDeviceSynchronize() is not needed here since all work will be
  //queued in the stream sequentially by the CPU to be executed by GPU.
  //kernel1<<<X,Y>>>(...); // kernel start execution, CPU continues to next
                           // statement
  //kernel2<<<X,Y>>>(...); // kernel is placed in queue and will start after
                           // kernel1 finishes, CPU continues to next statement
  //cudaMemcpy(...); // CPU blocks until ememory is copied, memory copy starts
                     // only after kernel2 finishes
}
*/

/* without reactions: 38.5 s
__global__
void concurrent_walk(
    const unsigned mol_size_,
    const unsigned seed_,
    const voxel_t stride_,
    const voxel_t id_stride_,
    const voxel_t vac_id_,
    const voxel_t null_id_,
    umol_t* mols_,
    const mol_t* offsets_,
    voxel_t* voxels_) {
  //index is the unique global thread id (size: total_threads)
  unsigned index(blockIdx.x*blockDim.x + threadIdx.x);
  const unsigned total_threads(blockDim.x*gridDim.x);
  while(index < mol_size_) {
    const umol_t vdx(mols_[index]);
    thrust::default_random_engine rng;
    rng.discard(seed_+index);
    thrust::uniform_int_distribution<unsigned> u(0, 11);
    const unsigned rand(u(rng));
    const bool odd_lay((vdx/NUM_COLROW)&1);
    const bool odd_col((vdx%NUM_COLROW/NUM_ROW)&1);
    mol2_t val(mol2_t(vdx)+offsets_[rand+(24&(-odd_lay))+(12&(-odd_col))]);
    //Atomically put the current molecule id, index+id_stride_ at the target
    //voxel if it is vacant: 
    const voxel_t tar_mol_id(atomicCAS(voxels_+val, vac_id_, index+id_stride_));
    //If not occupied, finalize walk:
    if(tar_mol_id == vac_id_) {
      voxels_[vdx] = vac_id_;
      mols_[index] = val;
    }
    //Do nothing, stay at original position
    index += total_threads;
  }
  //__syncthreads();
}

void Diffuser::walk() {
  const size_t size(mols_.size());
  concurrent_walk<<<blocks_, 512>>>(
      size,
      seed_,
      stride_,
      id_stride_,
      vac_id_,
      null_id_,
      thrust::raw_pointer_cast(&mols_[0]),
      thrust::raw_pointer_cast(&offsets_[0]),
      thrust::raw_pointer_cast(&voxels_[0]));
  //barrier cudaDeviceSynchronize() is not needed here since all work will be
  //queued in the stream sequentially by the CPU to be executed by GPU.
  //kernel1<<<X,Y>>>(...); // kernel start execution, CPU continues to next
                           // statement
  //kernel2<<<X,Y>>>(...); // kernel is placed in queue and will start after
                           // kernel1 finishes, CPU continues to next statement
  //cudaMemcpy(...); // CPU blocks until ememory is copied, memory copy starts
                     // only after kernel2 finishes
  seed_ += size;
}
*/

/* with minimal number of blocks: 41.2 s
__global__
void concurrent_walk(
    const unsigned mol_size_,
    const unsigned seed_,
    const voxel_t stride_,
    const voxel_t id_stride_,
    const voxel_t vac_id_,
    umol_t* mols_,
    const mol_t* offsets_,
    voxel_t* voxels_) {
  //index is the unique id of all the threads from all blocks
  unsigned index(blockIdx.x*blockDim.x + threadIdx.x);
  const unsigned total_threads(blockDim.x*gridDim.x);
  while(index < mol_size_) {
    const umol_t vdx(mols_[index]);
    thrust::default_random_engine rng;
    rng.discard(seed_+index);
    thrust::uniform_int_distribution<unsigned> u(0, 11);
    const unsigned rand(u(rng));
    const bool odd_lay((vdx/NUM_COLROW)&1);
    const bool odd_col((vdx%NUM_COLROW/NUM_ROW)&1);
    mol2_t val(mol2_t(vdx)+offsets_[rand+(24&(-odd_lay))+(12&(-odd_col))]);
    //Atomically put the current molecule id, index+id_stride_ at the target
    //voxel if it is vacant: 
    const voxel_t tar_mol_id(atomicCAS(voxels_+val, vac_id_, index+id_stride_));
    //If not occupied, finalize walk:
    if(tar_mol_id == vac_id_) {
      voxels_[vdx] = vac_id_;
      mols_[index] = val;
    }
    //Do nothing, stay at original position
    index += total_threads;
  }
  __syncthreads();
}

void Diffuser::walk() {
  const size_t size(mols_.size());
  concurrent_walk<<<blocks_, 512>>>(
      size,
      seed_,
      stride_,
      id_stride_,
      vac_id_,
      thrust::raw_pointer_cast(&mols_[0]),
      thrust::raw_pointer_cast(&offsets_[0]),
      thrust::raw_pointer_cast(&voxels_[0]));
  //barrier cudaDeviceSynchronize() is not needed here since all work will be
  //queued in the stream sequentially by the CPU to be executed by GPU.
  //kernel1<<<X,Y>>>(...); // kernel start execution, CPU continues to next
                           // statement
  //kernel2<<<X,Y>>>(...); // kernel is placed in queue and will start after
                           // kernel1 finishes, CPU continues to next statement
  //cudaMemcpy(...); // CPU blocks until ememory is copied, memory copy starts
                     // only after kernel2 finishes
  seed_ += size;
}
*/

/* with cudaDeviceSynchronize barrier: 42.8 s
__global__
void concurrent_walk(
    const unsigned mol_size_,
    const unsigned seed_,
    const voxel_t stride_,
    const voxel_t id_stride_,
    const voxel_t vac_id_,
    umol_t* mols_,
    const mol_t* offsets_,
    voxel_t* voxels_) {
  const unsigned index(blockIdx.x*blockDim.x + threadIdx.x);
  if(index < mol_size_) {
    const umol_t vdx(mols_[index]);
    thrust::default_random_engine rng;
    rng.discard(seed_+index);
    thrust::uniform_int_distribution<unsigned> u(0, 11);
    const unsigned rand(u(rng));
    const bool odd_lay((vdx/NUM_COLROW)&1);
    const bool odd_col((vdx%NUM_COLROW/NUM_ROW)&1);
    mol2_t val(mol2_t(vdx)+offsets_[rand+(24&(-odd_lay))+(12&(-odd_col))]);
    //Atomically put the current molecule id, index+id_stride_ at the target
    //voxel if it is vacant: 
    const voxel_t tar_mol_id(atomicCAS(voxels_+val, vac_id_, index+id_stride_));
    //If not occupied, finalize walk:
    if(tar_mol_id == vac_id_) {
      voxels_[vdx] = vac_id_;
      mols_[index] = val;
    }
    //Do nothing, stay at original position
  }
}

void Diffuser::walk() {
  const size_t size(mols_.size());
  concurrent_walk<<<(size+511)/512, 512>>>(
      size,
      seed_,
      stride_,
      id_stride_,
      vac_id_,
      thrust::raw_pointer_cast(&mols_[0]),
      thrust::raw_pointer_cast(&offsets_[0]),
      thrust::raw_pointer_cast(&voxels_[0]));
  //barrier until all CUDA calls have completed:
  cudaDeviceSynchronize();
  seed_ += size;
}
*/

/* Without thrust: 40.8 s
__global__
void concurrent_walk(
    const unsigned mol_size_,
    const unsigned seed_,
    const voxel_t stride_,
    const voxel_t id_stride_,
    const voxel_t vac_id_,
    umol_t* mols_,
    const mol_t* offsets_,
    voxel_t* voxels_) {
  const unsigned index(blockIdx.x*blockDim.x + threadIdx.x);
  if(index < mol_size_) {
    const umol_t vdx(mols_[index]);
    thrust::default_random_engine rng;
    rng.discard(seed_+index);
    thrust::uniform_int_distribution<unsigned> u(0, 11);
    const unsigned rand(u(rng));
    const bool odd_lay((vdx/NUM_COLROW)&1);
    const bool odd_col((vdx%NUM_COLROW/NUM_ROW)&1);
    mol2_t val(mol2_t(vdx)+offsets_[rand+(24&(-odd_lay))+(12&(-odd_col))]);
    //Atomically put the current molecule id, index+id_stride_ at the target
    //voxel if it is vacant: 
    const voxel_t tar_mol_id(atomicCAS(voxels_+val, vac_id_, index+id_stride_));
    //If not occupied, finalize walk:
    if(tar_mol_id == vac_id_) {
      voxels_[vdx] = vac_id_;
      mols_[index] = val;
    }
    //Do nothing, stay at original position
  }
}

void Diffuser::walk() {
  const size_t size(mols_.size());
  concurrent_walk<<<(size+511)/512, 512>>>(
      size,
      seed_,
      stride_,
      id_stride_,
      vac_id_,
      thrust::raw_pointer_cast(&mols_[0]),
      thrust::raw_pointer_cast(&offsets_[0]),
      thrust::raw_pointer_cast(&voxels_[0]));
  seed_ += size;
}
*/


/* Verified random walk without any reaction checks: 41.9 s
struct generate {
  __host__ __device__ generate(
      const unsigned mol_size,
      const unsigned seed,
      const voxel_t stride,
      const voxel_t id_stride,
      const voxel_t vac_id,
      const mol_t* offsets,
      voxel_t* voxels):
    mol_size_(mol_size),
    seed_(seed),
    stride_(stride),
    id_stride_(id_stride),
    vac_id_(vac_id),
    offsets_(offsets),
    voxels_(voxels) {} 
  __device__ umol_t operator()(const unsigned index, const umol_t vdx) const {
    thrust::default_random_engine rng;
    rng.discard(seed_+index);
    thrust::uniform_int_distribution<unsigned> u(0, 11);
    const unsigned rand(u(rng));
    const bool odd_lay((vdx/NUM_COLROW)&1);
    const bool odd_col((vdx%NUM_COLROW/NUM_ROW)&1);
    mol2_t val(mol2_t(vdx)+offsets_[rand+(24&(-odd_lay))+(12&(-odd_col))]);
    //Atomically put the current molecule id, index+id_stride_ at the target
    //voxel if it is vacant: 
    const voxel_t tar_mol_id(atomicCAS(voxels_+val, vac_id_, index+id_stride_));
    //If not occupied, finalize walk:
    if(tar_mol_id == vac_id_) {
      voxels_[vdx] = vac_id_;
      return val;
    }
    //Stay at original position:
    return vdx;
  }
  const unsigned mol_size_;
  const unsigned seed_;
  const voxel_t stride_;
  const voxel_t id_stride_;
  const voxel_t vac_id_;
  const mol_t* offsets_;
  voxel_t* voxels_;
};

void Diffuser::walk() {
  const size_t size(mols_.size());
  thrust::transform(thrust::device, 
      thrust::counting_iterator<unsigned>(0),
      thrust::counting_iterator<unsigned>(size),
      mols_.begin(),
      mols_.begin(),
      generate(
        size,
        seed_,
        stride_,
        id_stride_,
        vac_id_,
        thrust::raw_pointer_cast(&offsets_[0]),
        thrust::raw_pointer_cast(&voxels_[0])));
  seed_ += size;
}
*/

/* Verified random walk without reaction: 42.0 s
struct generate {
  __host__ __device__ generate(
      const unsigned mol_size,
      const unsigned seed,
      const voxel_t stride,
      const voxel_t id_stride,
      const voxel_t vac_id,
      const bool* is_reactive,
      const mol_t* offsets,
      umol_t* reacteds,
      voxel_t* voxels):
    mol_size_(mol_size),
    seed_(seed),
    stride_(stride),
    id_stride_(id_stride),
    vac_id_(vac_id),
    is_reactive_(is_reactive),
    offsets_(offsets),
    reacteds_(reacteds),
    voxels_(voxels) {} 
  __device__ umol_t operator()(const unsigned index, const umol_t vdx) const {
    thrust::default_random_engine rng;
    rng.discard(seed_+index);
    thrust::uniform_int_distribution<unsigned> u(0, 11);
    const unsigned rand(u(rng));
    const bool odd_lay((vdx/NUM_COLROW)&1);
    const bool odd_col((vdx%NUM_COLROW/NUM_ROW)&1);
    mol2_t val(mol2_t(vdx)+offsets_[rand+(24&(-odd_lay))+(12&(-odd_col))]);
    //Atomically put the current molecule id, index+id_stride_ at the target
    //voxel if it is vacant: 
    const voxel_t tar_mol_id(atomicCAS(voxels_+val, vac_id_, index+id_stride_));
    //If not occupied, finalize walk:
    if(tar_mol_id == vac_id_) {
      voxels_[vdx] = vac_id_;
      return val;
    }
    //it is occupied, so check if it is reactive and add reacted:
    const voxel_t tar_id(tar_mol_id/stride_);
    if(is_reactive_[tar_id]) {
      //const unsigned old(atomicAdd(reacteds_+mol_size_, 1));
      reacteds_[index] = tar_mol_id;
    }
    //Stay at original position:
    return vdx;
  }
  const unsigned mol_size_;
  const unsigned seed_;
  const voxel_t stride_;
  const voxel_t id_stride_;
  const voxel_t vac_id_;
  const bool* is_reactive_;
  const mol_t* offsets_;
  umol_t* reacteds_;
  voxel_t* voxels_;
};

struct is_reacted {
  __device__ bool operator()(const umol_t reacted) {
    return reacted;
  }
};

struct react {
  __host__ __device__ react(
      const unsigned mol_size,
      umol_t* reacteds):
    mol_size_(mol_size),
    reacteds_(reacteds) {}
  __device__ umol_t operator()(const unsigned index, const umol_t mol) const {
      const unsigned old(atomicSub(reacteds_+mol_size_, 1));
    reacteds_[index] = 0;
    return mol;
  }
  const unsigned mol_size_;
  umol_t* reacteds_;
};

void Diffuser::walk() {
  const size_t size(mols_.size());
  reacteds_.resize(size+1);
  thrust::transform(thrust::device, 
      thrust::counting_iterator<unsigned>(0),
      thrust::counting_iterator<unsigned>(size),
      mols_.begin(),
      mols_.begin(),
      generate(
        size,
        seed_,
        stride_,
        id_stride_,
        vac_id_,
        thrust::raw_pointer_cast(&is_reactive_[0]),
        thrust::raw_pointer_cast(&offsets_[0]),
        thrust::raw_pointer_cast(&reacteds_[0]),
        thrust::raw_pointer_cast(&voxels_[0])));
  seed_ += size;
}
*/

/*
struct generate {
  __host__ __device__ generate(
      const unsigned mol_size,
      const unsigned seed,
      const voxel_t stride,
      const voxel_t id_stride,
      const voxel_t vac_id,
      const bool* is_reactive,
      const mol_t* offsets,
      umol_t* reacteds,
      voxel_t* voxels):
    mol_size_(mol_size),
    seed_(seed),
    stride_(stride),
    id_stride_(id_stride),
    vac_id_(vac_id),
    is_reactive_(is_reactive),
    offsets_(offsets),
    reacteds_(reacteds),
    voxels_(voxels) {} 
  __device__ umol_t operator()(const unsigned index, const umol_t vdx) const {
    thrust::default_random_engine rng;
    rng.discard(seed_+index);
    thrust::uniform_int_distribution<unsigned> u(0, 11);
    const unsigned rand(u(rng));
    const bool odd_lay((vdx/NUM_COLROW)&1);
    const bool odd_col((vdx%NUM_COLROW/NUM_ROW)&1);
    mol2_t val(mol2_t(vdx)+offsets_[rand+(24&(-odd_lay))+(12&(-odd_col))]);
    //Atomically put the current molecule id, index+id_stride_ at the target
    //voxel if it is vacant: 
    const voxel_t tar_mol_id(atomicCAS(voxels_+val, vac_id_, index+id_stride_));
    //If not occupied, finalize walk:
    if(tar_mol_id == vac_id_) {
      voxels_[vdx] = vac_id_;
      return val;
    }
    //it is occupied, so check if it is reactive and add reacted:
    const voxel_t tar_id(tar_mol_id/stride_);
    if(is_reactive_[tar_id]) {
      //const unsigned old(atomicAdd(reacteds_+mol_size_, 1));
      reacteds_[index] = tar_mol_id;
    }
    //Stay at original position:
    return vdx;
  }
  const unsigned mol_size_;
  const unsigned seed_;
  const voxel_t stride_;
  const voxel_t id_stride_;
  const voxel_t vac_id_;
  const bool* is_reactive_;
  const mol_t* offsets_;
  umol_t* reacteds_;
  voxel_t* voxels_;
};

struct is_reacted {
  __device__ bool operator()(const umol_t reacted) {
    return reacted;
  }
};

struct react {
  __host__ __device__ react(
      const unsigned mol_size,
      umol_t* reacteds):
    mol_size_(mol_size),
    reacteds_(reacteds) {}
  __device__ umol_t operator()(const unsigned index, const umol_t mol) const {
      const unsigned old(atomicSub(reacteds_+mol_size_, 1));
    reacteds_[index] = 0;
    return mol;
  }
  const unsigned mol_size_;
  umol_t* reacteds_;
};

void Diffuser::walk() {
  const size_t size(mols_.size());
  reacteds_.resize(size+1);
  thrust::transform(thrust::device, 
      thrust::counting_iterator<unsigned>(0),
      thrust::counting_iterator<unsigned>(size),
      mols_.begin(),
      mols_.begin(),
      generate(
        size,
        seed_,
        stride_,
        id_stride_,
        vac_id_,
        thrust::raw_pointer_cast(&is_reactive_[0]),
        thrust::raw_pointer_cast(&offsets_[0]),
        thrust::raw_pointer_cast(&reacteds_[0]),
        thrust::raw_pointer_cast(&voxels_[0])));
  thrust::transform_if(thrust::device,
      thrust::counting_iterator<unsigned>(0),
      thrust::counting_iterator<unsigned>(size),
      mols_.begin(),
      reacteds_.begin(),
      mols_.begin(),
      react(
        size,
        thrust::raw_pointer_cast(&reacteds_[0])),
      is_reacted());
  seed_ += size;
}
*/

/* Verified random walk with thrust::random : 44.1 s
struct generate {
  __host__ __device__ generate(
      const unsigned seed,
      const voxel_t stride,
      const voxel_t id_stride,
      const voxel_t vac_id,
      const bool* is_reactive,
      const mol_t* offsets,
      umol_t* reacteds,
      voxel_t* voxels):
    seed_(seed),
    stride_(stride),
    id_stride_(id_stride),
    vac_id_(vac_id),
    is_reactive_(is_reactive),
    offsets_(offsets),
    reacteds_(reacteds),
    voxels_(voxels) {} 
  __device__ umol_t operator()(const unsigned index, const umol_t vdx) const {
    thrust::default_random_engine rng;
    rng.discard(seed_+index);
    thrust::uniform_int_distribution<unsigned> u(0, 11);
    const unsigned rand(u(rng));
    const bool odd_lay((vdx/NUM_COLROW)&1);
    const bool odd_col((vdx%NUM_COLROW/NUM_ROW)&1);
    mol2_t val(mol2_t(vdx)+offsets_[rand+(24&(-odd_lay))+(12&(-odd_col))]);
    const voxel_t res(atomicCAS(voxels_+val, vac_id_, index+id_stride_));
    //If not occupied, walk:
    if(res == vac_id_) {
      voxels_[vdx] = vac_id_;
      reacteds_[index] = 0;
      return val;
    }
    //If occupied, check and add reacted:
    const voxel_t tar_id(res/stride_);
    if(is_reactive_[tar_id]) {
      reacteds_[index] = res;
    }
    else {
      reacteds_[index] = 0;
    }
    //Stay at original position:
    return vdx;
  }
  const unsigned seed_;
  const voxel_t stride_;
  const voxel_t id_stride_;
  const voxel_t vac_id_;
  const bool* is_reactive_;
  const mol_t* offsets_;
  umol_t* reacteds_;
  voxel_t* voxels_;
};

struct is_reacted {
  __device__ bool operator()(const umol_t reacted) {
    return reacted;
  }
};

struct react {
  __device__ umol_t operator()(const umol_t mol, const umol_t reacted) const {
    return mol;
  }
};

void Diffuser::walk() {
  const size_t size(mols_.size());
  reacteds_.resize(size);
  thrust::transform(thrust::device, 
      thrust::counting_iterator<unsigned>(0),
      thrust::counting_iterator<unsigned>(size),
      mols_.begin(),
      mols_.begin(),
      generate(
        seed_,
        stride_,
        id_stride_,
        vac_id_,
        thrust::raw_pointer_cast(&is_reactive_[0]),
        thrust::raw_pointer_cast(&offsets_[0]),
        thrust::raw_pointer_cast(&reacteds_[0]),
        thrust::raw_pointer_cast(&voxels_[0])));
  thrust::transform_if(thrust::device,
      mols_.begin(),
      mols_.end(),
      reacteds_.begin(),
      reacteds_.begin(),
      mols_.begin(),
      react(),
      is_reacted());
  seed_ += size;
}
*/

/*
// With pseudo reaction and non-overlap population: 43.2 s
struct is_reacted {
  __device__ bool operator()(const umol_t reacted) {
    return reacted;
  }
};

struct react {
  __device__ umol_t operator()(const umol_t mol, const umol_t reacted) const {
    return mol;
  }
};

void Diffuser::walk() {
  const size_t size(mols_.size());
  reacteds_.resize(size);
  thrust::transform(thrust::device, 
      thrust::counting_iterator<unsigned>(0),
      thrust::counting_iterator<unsigned>(size),
      mols_.begin(),
      mols_.begin(),
      generate(
        seed_,
        stride_,
        id_stride_,
        vac_id_,
        thrust::raw_pointer_cast(&is_reactive_[0]),
        thrust::raw_pointer_cast(&offsets_[0]),
        thrust::raw_pointer_cast(&reacteds_[0]),
        thrust::raw_pointer_cast(&voxels_[0])));
  thrust::transform_if(thrust::device,
      mols_.begin(),
      mols_.end(),
      reacteds_.begin(),
      reacteds_.begin(),
      mols_.begin(),
      react(),
      is_reacted());
  seed_ += size;
}
*/

/*
//With reaction check list: 41.5 s
struct generate {
  __host__ __device__ generate(
      const unsigned seed,
      const voxel_t stride,
      const voxel_t id_stride,
      const voxel_t vac_id,
      const bool* is_reactive,
      const mol_t* offsets,
      umol_t* reacteds,
      voxel_t* voxels):
    seed_(seed),
    stride_(stride),
    id_stride_(id_stride),
    vac_id_(vac_id),
    is_reactive_(is_reactive),
    offsets_(offsets),
    reacteds_(reacteds),
    voxels_(voxels) {} 
  __device__ umol_t operator()(const unsigned index, const umol_t vdx) const {
    curandState s;
    curand_init(seed_+index, 0, 0, &s);
    float ranf(curand_uniform(&s)*11.999999);
    const unsigned rand((unsigned)truncf(ranf));
    const bool odd_lay((vdx/NUM_COLROW)&1);
    const bool odd_col((vdx%NUM_COLROW/NUM_ROW)&1);
    mol2_t val(mol2_t(vdx)+offsets_[rand+(24&(-odd_lay))+(12&(-odd_col))]);
    const voxel_t res(atomicCAS(voxels_+val, vac_id_, index+id_stride_));
    //If not occupied, walk:
    if(res == vac_id_) {
      voxels_[vdx] = vac_id_;
      reacteds_[index] = 0;
      return val;
    }
    //If occupied, check and add reacted:
    const voxel_t tar_id(res/stride_);
    if(is_reactive_[tar_id]) {
      reacteds_[index] = res;
    }
    else {
      reacteds_[index] = 0;
    }
    //Stay at original position:
    return vdx;
  }
  const unsigned seed_;
  const voxel_t stride_;
  const voxel_t id_stride_;
  const voxel_t vac_id_;
  const bool* is_reactive_;
  const mol_t* offsets_;
  umol_t* reacteds_;
  voxel_t* voxels_;
};

void Diffuser::walk() {
  const size_t size(mols_.size());
  reacteds_.resize(size);
  thrust::transform(thrust::device, 
      thrust::counting_iterator<unsigned>(0),
      thrust::counting_iterator<unsigned>(size),
      mols_.begin(),
      mols_.begin(),
      generate(
        seed_,
        stride_,
        id_stride_,
        vac_id_,
        thrust::raw_pointer_cast(&is_reactive_[0]),
        thrust::raw_pointer_cast(&offsets_[0]),
        thrust::raw_pointer_cast(&reacteds_[0]),
        thrust::raw_pointer_cast(&voxels_[0])));
  seed_ += size;
}
*/

/*
//Use atomicCAS to avoid race condition: 39.1 s
struct generate {
  __host__ __device__ generate(const mol_t* _offsets, voxel_t* _voxels):
    offsets(_offsets), voxels(_voxels) {} 
  __device__ umol_t operator()(const unsigned n, const umol_t vdx) const {
    curandState s;
    curand_init(n, 0, 0, &s);
    float ranf(curand_uniform(&s)*11.999999);
    const unsigned rand((unsigned)truncf(ranf));
    const bool odd_lay((vdx/NUM_COLROW)&1);
    const bool odd_col((vdx%NUM_COLROW/NUM_ROW)&1);
    mol2_t val(mol2_t(vdx)+offsets[rand+(24&(-odd_lay))+(12&(-odd_col))]);
    const voxel_t res(atomicCAS(voxels+val, 0, 1));
    if(res == 0) {
      voxels[vdx] = 0;
      return val;
    }
    return vdx;
  }
  const mol_t* offsets;
  voxel_t* voxels;
};

void Diffuser::walk() {
  const size_t size(mols_.size());
  thrust::transform(thrust::device, 
      thrust::counting_iterator<unsigned>(seed_),
      thrust::counting_iterator<unsigned>(seed_+size),
      mols_.begin(),
      mols_.begin(),
      generate(thrust::raw_pointer_cast(&offsets_[0]), thrust::raw_pointer_cast(&voxels_[0])));
  seed_ += size;
}
*/

/*
//Perform walk in predicate, only single transformation function: 38.7 s
struct generate {
  __host__ __device__ generate(const mol_t* _offsets, voxel_t* _voxels):
    offsets(_offsets), voxels(_voxels) {} 
  __device__ umol_t operator()(const unsigned n, const umol_t vdx) const {
    curandState s;
    curand_init(n, 0, 0, &s);
    float ranf(curand_uniform(&s)*11.999999);
    const unsigned rand((unsigned)truncf(ranf));
    const bool odd_lay((vdx/NUM_COLROW)&1);
    const bool odd_col((vdx%NUM_COLROW/NUM_ROW)&1);
    mol2_t val(mol2_t(vdx)+offsets[rand+(24&(-odd_lay))+(12&(-odd_col))]);
    if(val < 0 || val > NUM_VOXEL) {
      val = vdx;
    }
    if(!voxels[val]) {
      voxels[val] = true;
      voxels[vdx] = false;
      return val;
    }
    else {
      return vdx;
    }
  }
  const mol_t* offsets;
  voxel_t* voxels;
};

void Diffuser::walk() {
  const size_t size(mols_.size());
  thrust::transform(thrust::device, 
      thrust::counting_iterator<unsigned>(seed_),
      thrust::counting_iterator<unsigned>(seed_+size),
      mols_.begin(),
      mols_.begin(),
      generate(thrust::raw_pointer_cast(&offsets_[0]), thrust::raw_pointer_cast(&voxels_[0])));
  seed_ += size;
}
*/

/*
//Simplified full collision with transform_if: very small performance
//improvement
struct generate {
  __host__ __device__ generate(const mol_t* _offsets):
    offsets(_offsets) {} 
  __device__ umol_t operator()(const unsigned n, const umol_t vdx) const {
    curandState s;
    curand_init(n, 0, 0, &s);
    float ranf(curand_uniform(&s)*11.999999);
    const unsigned rand((unsigned)truncf(ranf));
    const bool odd_lay((vdx/NUM_COLROW)&1);
    const bool odd_col((vdx%NUM_COLROW/NUM_ROW)&1);
    mol2_t val(mol2_t(vdx)+offsets[rand+(24&(-odd_lay))+(12&(-odd_col))]);
    if(val < 0 || val > NUM_VOXEL) {
      return vdx;
    }
    return val;
  }
  const mol_t* offsets;
};

struct is_occupied {
  __device__ bool operator()(const voxel_t voxel) {
    return ((bool)voxel);
  }
};

struct update {
  __device__ umol_t operator()(const umol_t mol) const {
    return mol;
  }
};

void Diffuser::walk() {
  const size_t size(mols_.size());
  tars_.resize(size);
  thrust::transform(thrust::device, 
      thrust::counting_iterator<unsigned>(seed_),
      thrust::counting_iterator<unsigned>(seed_+size),
      mols_.begin(),
      tars_.begin(),
      generate(thrust::raw_pointer_cast(&offsets_[0])));
  thrust::permutation_iterator<thrust::device_vector<voxel_t>::iterator,
    thrust::device_vector<umol_t>::iterator> stencil(voxels_.begin(),
        tars_.begin());
  thrust::transform_if(
      mols_.begin(),
      mols_.end(),
      stencil,
      tars_.begin(),
      update(),
      is_occupied());
  thrust::permutation_iterator<thrust::device_vector<voxel_t>::iterator,
    thrust::device_vector<umol_t>::iterator> vacants(voxels_.begin(),
        mols_.begin());
  thrust::fill_n(thrust::device, vacants, size, false);
  thrust::permutation_iterator<thrust::device_vector<voxel_t>::iterator,
    thrust::device_vector<umol_t>::iterator> occupieds(voxels_.begin(),
        tars_.begin());
  thrust::fill_n(thrust::device, occupieds, size, true);
  thrust::copy(tars_.begin(), tars_.end(), mols_.begin());
  seed_ += size;
}
*/


/*
//Full collision lattice with transform_if: same performance
struct is_vacant {
  __device__ bool operator()(const voxel_t voxel) {
    return (bool)voxel;
  }
};

struct update {
  __device__ umol_t operator()(const voxel_t tar, const umol_t mol) const {
    return mol;
  }
};

void Diffuser::walk() {
  const size_t size(mols_.size());
  tars_.resize(size);
  thrust::transform(thrust::device, 
      thrust::counting_iterator<unsigned>(seed_),
      thrust::counting_iterator<unsigned>(seed_+size),
      mols_.begin(),
      tars_.begin(),
      generate(thrust::raw_pointer_cast(&offsets_[0])));
  thrust::permutation_iterator<thrust::device_vector<voxel_t>::iterator,
    thrust::device_vector<umol_t>::iterator> stencil(lattice_.begin(),
        tars_.begin());
  thrust::transform_if(thrust::device, 
      tars_.begin(),
      tars_.end(),
      mols_.begin(),
      stencil,
      tars_.begin(),
      update(),
      is_vacant());
  thrust::permutation_iterator<thrust::device_vector<voxel_t>::iterator,
    thrust::device_vector<umol_t>::iterator> vacants(lattice_.begin(),
        mols_.begin());
  thrust::fill_n(thrust::device, vacants, size, true);
  thrust::permutation_iterator<thrust::device_vector<voxel_t>::iterator,
    thrust::device_vector<umol_t>::iterator> occupieds(lattice_.begin(),
        tars_.begin());
  thrust::fill_n(thrust::device, occupieds, size, false);
  thrust::copy(tars_.begin(), tars_.end(), mols_.begin());
  //thrust::copy(mols_.begin(), mols_.end(), box_mols_[0].begin());
  seed_ += size;
}
*/


/*
//Full collision check with lattice: 0.3 s (1,600,000 molecules, 42.5 s)
struct update {
  __host__ __device__ update(const voxel_t* _lattice):
    lattice(_lattice) {} 
  __device__ umol_t operator()(const umol_t tar, const umol_t mol) const {
    if(lattice[tar]) {
      return tar;
    }
    return mol;
  }
  const voxel_t* lattice;
};

void Diffuser::walk() {
  const size_t size(mols_.size());
  tars_.resize(size);
  thrust::transform(thrust::device, 
      thrust::counting_iterator<unsigned>(seed_),
      thrust::counting_iterator<unsigned>(seed_+size),
      mols_.begin(),
      tars_.begin(),
      generate(thrust::raw_pointer_cast(&offsets_[0])));
  thrust::transform(thrust::device, 
      tars_.begin(),
      tars_.end(),
      mols_.begin(),
      tars_.begin(),
      update(thrust::raw_pointer_cast(&lattice_[0])));
  thrust::permutation_iterator<thrust::device_vector<voxel_t>::iterator,
    thrust::device_vector<umol_t>::iterator> vacants(lattice_.begin(),
        mols_.begin());
  thrust::fill_n(thrust::device, vacants, mols_.size(), true);
  thrust::permutation_iterator<thrust::device_vector<voxel_t>::iterator,
    thrust::device_vector<umol_t>::iterator> occupieds(lattice_.begin(),
        tars_.begin());
  thrust::fill_n(thrust::device, occupieds, tars_.size(), false);
  thrust::copy(tars_.begin(), tars_.end(), mols_.begin());
  //thrust::copy(mols_.begin(), mols_.end(), box_mols_[0].begin());
  seed_ += mols_.size();
}
*/

/*
//Partial collision check with lattice: 0.2 s
struct update {
  __host__ __device__ update(bool* _lattice):
    lattice(_lattice) {} 
  __device__ umol_t operator()(const umol_t tar, const umol_t mol) const {
    if(lattice[tar]) {
      lattice[tar] = false;
      lattice[mol] = true;
      return tar;
    }
    return mol;
  }
  bool* lattice;
};

void Diffuser::walk() {
  const size_t size(mols_.size());
  tars_.resize(size);
  thrust::transform(thrust::device, 
      thrust::counting_iterator<unsigned>(seed_),
      thrust::counting_iterator<unsigned>(seed_+size),
      mols_.begin(),
      tars_.begin(),
      generate(thrust::raw_pointer_cast(&offsets_[0])));
  thrust::transform(thrust::device, 
      tars_.begin(),
      tars_.end(),
      mols_.begin(),
      mols_.begin(),
      update(thrust::raw_pointer_cast(&lattice_[0])));
  //thrust::copy(mols_.begin(), mols_.end(), box_mols_[0].begin());
  seed_ += mols_.size();
}
*/

/*
//Full collision check using molecule list only: 4.3s

struct update {
  __host__ __device__ update(const size_t _size, const umol_t* _mols):
    size(_size), mols(_mols) {} 
  __device__ umol_t operator()(const umol_t tar, const umol_t mol) const {
    for(unsigned i(0); i != size; ++i) {
      if(tar == mols[i]) {
        return mol;
      }
    }
    return tar;
  }
  const size_t size;
  const umol_t* mols;
};

void Diffuser::walk() {
  const size_t size(mols_.size());
  tars_.resize(size);
  thrust::transform(thrust::device, 
      thrust::counting_iterator<unsigned>(seed_),
      thrust::counting_iterator<unsigned>(seed_+size),
      mols_.begin(),
      tars_.begin(),
      generate(thrust::raw_pointer_cast(&offsets_[0])));
  thrust::transform(thrust::device, 
      tars_.begin(),
      tars_.end(),
      mols_.begin(),
      tars_.begin(),
      update(size, thrust::raw_pointer_cast(&mols_[0])));
  thrust::copy(tars_.begin(), tars_.end(), mols_.begin());
  //thrust::copy(mols_.begin(), mols_.end(), box_mols_[0].begin());
  seed_ += mols_.size();
}
*/


/*
struct generate {
  __host__ __device__ generate(mol_t* _offsets):
    offsets(_offsets) {} 
  __device__ umol_t operator()(const unsigned n, const umol_t vdx) const {
    thrust::default_random_engine rng(hash(n));
    thrust::uniform_int_distribution<int> uniform(0, 11);
    const bool odd_lay((vdx/NUM_COLROW)&1);
    const bool odd_col((vdx%NUM_COLROW/NUM_ROW)&1);
    mol2_t val(mol2_t(vdx)+offsets[uniform(rng)+(24&(-odd_lay))+
        (12&(-odd_col))]);
    if(val < 0 || val > NUM_VOXEL) {
      return vdx;
    }
    return val;
  }
  mol_t* offsets;
};
*/


/*
//Collisions not check, with intersection: 7s
void Diffuser::walk() {
  tars_.resize(mols_.size());
  thrust::transform(thrust::device, 
      thrust::counting_iterator<unsigned>(seed_),
      thrust::counting_iterator<unsigned>(seed_+mols_.size()),
      mols_.begin(),
      tars_.begin(),
      generate(thrust::raw_pointer_cast(&offsets_[0])));
  thrust::sort(thrust::device, tars_.begin(), tars_.end());
  collisions_.resize(mols_.size());
  thrust::set_intersection(thrust::device, mols_.begin(), mols_.end(),
      tars_.begin(), tars_.end(), collisions_.begin());
  //if(!collisions.size()) { 
    thrust::copy(tars_.begin(), tars_.end(), mols_.begin());
  //}
  //thrust::copy(mols_.begin(), mols_.end(), box_mols_[0].begin());
  seed_ += mols_.size();
}
*/

/*
//Sequential original 10.5 s
void Diffuser::walk(umol_t* mols, const unsigned size) {
  for (unsigned i(0); i != size; ++i) {
    umol_t tar(compartment_.get_tar(mols[i], rng_.Ran16_12()));
    for(unsigned j(0); j != size; ++j) {
      if(mols[j] == tar) {
        goto next;
      }
    }
    mols[i] = tar;
next:
    continue;
  }
}
*/
