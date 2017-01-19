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


#ifndef __Model_hpp
#define __Model_hpp

#include <Common.hpp>
#include <Compartment.hpp>
#include <Stepper.hpp>
#include <thrust/device_vector.h>
#include <curand_kernel.h>
#include <curand.h>
#include <sstream>

extern __device__ curandState* curand_states[64];

class Model {
 public: 
  Model();
  ~Model();
  void initialize();
  void step(const unsigned);
  unsigned run(const double);
  unsigned push_species(Species&);
  Compartment& get_compartment();
  Stepper& get_stepper();
  std::vector<Species*>& get_species();
  voxel_t get_null_id() const;
  voxel_t get_stride() const;
  unsigned& get_blocks();
 private:
  void initialize_random_generator();
  std::vector<Species*> species_;
  const voxel_t null_id_;
  unsigned blocks_;
  Stepper stepper_;
  Compartment compartment_; //must declare this after initializing others
  voxel_t stride_;
};

#endif /* __Model_hpp */
