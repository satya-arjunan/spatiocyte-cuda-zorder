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


#ifndef __Compartment_hpp
#define __Compartment_hpp

#include <Spatiocyte.hpp>
#include <Common.hpp>
#include <Species.hpp>
#include <Lattice.hpp>

#define HCP_X double(VOXEL_RADIUS*1.732050807568877)
#define HCP_Z double(VOXEL_RADIUS*1.632993161855452)
#define NUM_COL mol_t(LENGTH_X/HCP_X+3)
#define NUM_LAY mol_t(LENGTH_Z/HCP_Z+3)
#define NUM_ROW mol_t(LENGTH_Y/VOXEL_RADIUS/2+3)
#define NUM_COLROW mol_t(NUM_COL*NUM_ROW)
#define NUM_COLROWROW umol_t(umol_t(NUM_COLROW)*NUM_ROW)
#define NUM_VOXEL umol_t(umol_t(NUM_COLROW)*NUM_LAY)

class Compartment { 
 public: 
  Compartment(std::string, const double, const double, const double, Model&);
  ~Compartment() {}
  void initialize();
  Vector<double> get_center() const;
  const Vector<double>& get_dimensions() const;
  Species& get_surface_species();
  Species& get_volume_species();
  Model& get_model();
  const std::string& get_name() const;
  Lattice& get_lattice();
  thrust::device_vector<mol_t>& get_offsets();
 private:
  void set_volume_structure();
  void set_surface_structure();
  void populate_mol(const umol_t);
  void set_offsets();
 private:
  const std::string name_;
  Model& model_;
  thrust::device_vector<mol_t> offsets_;
  Lattice lattice_;
  const Vector<double> dimensions_;
  Species volume_species_;
  Species surface_species_;
};

#endif /* __Compartment_hpp */

