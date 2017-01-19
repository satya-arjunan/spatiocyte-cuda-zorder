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


#ifndef __Lattice_hpp
#define __Lattice_hpp

#include <thrust/device_vector.h>
#include <Common.hpp>

class Lattice {
 public: 
  Lattice(const Vector<unsigned>&);
  ~Lattice() {}
  void initialize();
  unsigned get_num_voxel() const; 
  bool is_mol_at_edge(const umol_t) const;
  umol_t coord_to_mol(const Vector<unsigned>&) const;
  Vector<unsigned> mol_to_coord(const umol_t) const;
  const Vector<unsigned>& get_dimensions() const;
  thrust::device_vector<voxel_t>& get_voxels();
 private:
  const Vector<unsigned> dimensions_;
  thrust::device_vector<voxel_t> voxels_;
};

#endif /* __Lattice_hpp */

