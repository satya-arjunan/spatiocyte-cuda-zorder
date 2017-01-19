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

#include <cmath>
#include <cstring>
#include <Lattice.hpp>
#include <Compartment.hpp>

Lattice::Lattice(const Vector<unsigned>& dimensions):
  voxels_(NUM_VOXEL, 0),
  dimensions_(dimensions) {
}

void Lattice::initialize() {
  std::cout << "num_x:" << get_dimensions().x << " num_y:" << 
    get_dimensions().y << " num_z:" << get_dimensions().z  << " num_voxel:" <<
    get_num_voxel() << " memory:" << get_num_voxel()*sizeof(voxel_t)/
    (1024*1024.0) << " MB" << std::endl;
}

unsigned Lattice::get_num_voxel() const {
  return voxels_.size();
}

bool Lattice::is_mol_at_edge(const umol_t mol) const {
  const Vector<unsigned>& coord(mol_to_coord(mol));
  if(coord.x == 0 || coord.y == 0 || coord.z == 0 ||
     coord.x == dimensions_.x-1 ||
     coord.y == dimensions_.y-1 ||
     coord.z == dimensions_.z-1) {
    return true;
  }
  return false;
}

umol_t Lattice::coord_to_mol(const Vector<unsigned>& coord) const {
  return coord.x*dimensions_.y+coord.y+coord.z*dimensions_.x*dimensions_.y;
}

Vector<unsigned> Lattice::mol_to_coord(const umol_t mol) const {
  const unsigned xy(dimensions_.x*dimensions_.y);
  return Vector<unsigned>(mol%xy/dimensions_.y, mol%xy%dimensions_.y,
                          mol/xy);
}

const Vector<unsigned>& Lattice::get_dimensions() const {
  return dimensions_;
}

thrust::device_vector<voxel_t>& Lattice::get_voxels() {
  return voxels_;
}
