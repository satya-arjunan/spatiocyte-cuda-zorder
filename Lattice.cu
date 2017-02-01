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
  voxels_(rint(double(NUM_VOXEL)/WORD), 0),
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

//Lattice voxel original convention:
//x:cols
//y:rows
//z:layers
//mol = y + x*rows + z*cols*rows
//Zorder convention:
//encode[z + y*rows + x*cols*rows]
void Lattice::i2zorder_xyz(const umol_t mol, humol_t& zx, humol_t& zy,
    humol_t& zz) {
  const umol_t cols(dimensions_.x);
  const umol_t rows(dimensions_.y);
  const umol_t colrows(rows*cols);
  const humol_t x(mol%colrows/rows);
  const humol_t y(mol%colrows%rows);
  const humol_t z(mol/colrows);
  zx = z; 
  zy = x;
  zz = y;
}

umol_t Lattice::zorder_xyz2i(const humol_t zx, const humol_t zy,
   const humol_t zz) {
  return zy*dimensions_.y+zz+zx*dimensions_.x*dimensions_.y;
}

umol_t Lattice::split_3bits(const humol_t a) {
	umol_t x = a;
	x = x & 0x000003ff;
	x = (x | x << 16) & 0x30000ff;
	x = (x | x << 8)  & 0x0300f00f;
	x = (x | x << 4)  & 0x30c30c3;
	x = (x | x << 2)  & 0x9249249;
	return x;
}

umol_t Lattice::encode_zorder(const humol_t x, const humol_t y,
    const humol_t z){
	return split_3bits(x) |
    (split_3bits(y) << 1) |
    (split_3bits(z) << 2);
}

umol_t Lattice::i2z(const umol_t vdx) {
  humol_t x,y,z;
  i2zorder_xyz(vdx, x, y, z);
  umol_t val(encode_zorder(x, y, z));
  if(val > dimensions_.x*dimensions_.y*dimensions_.z) {
    std::cout << "val:" << val << " max:" << dimensions_.x*dimensions_.y*dimensions_.z << " x:" << x << " y:" << y << " z:" << z << std::endl;
  }
  return val;
}

humol_t Lattice::get_third_bits(const umol_t m) {
	umol_t x = m & 0x9249249;
	x = (x ^ (x >> 2)) & 0x30c30c3;
	x = (x ^ (x >> 4)) & 0x0300f00f;
	x = (x ^ (x >> 8)) & 0x30000ff;
	x = (x ^ (x >> 16)) & 0x000003ff;
	return static_cast<humol_t>(x);
}

void Lattice::decode_zorder(const umol_t m, humol_t& x, humol_t& y,
    humol_t& z){
	x = get_third_bits(m);
	y = get_third_bits(m >> 1);
	z = get_third_bits(m >> 2);
}

umol_t Lattice::z2i(const umol_t zval) {
  humol_t x,y,z;
  decode_zorder(zval, x, y, z);
  return zorder_xyz2i(x, y, z);
}

const Vector<unsigned>& Lattice::get_dimensions() const {
  return dimensions_;
}

thrust::device_vector<voxel_t>& Lattice::get_voxels() {
  return voxels_;
}
