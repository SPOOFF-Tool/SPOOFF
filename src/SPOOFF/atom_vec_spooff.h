/* -*- c++ -*- ----------------------------------------------------------
 *
 * *** Smoothed Particle Open-source in-Operation Fuel Fragmentation ***
 *
 * The SPOOFF package is a modification of the MACHDYN package specifically 
 * designed for modelling fragmentation within nuclear fuels. This was 
 * created in collaboration betweein UKNNL and University of Bangor Nuclear 
 * Futures Centre under the OperaHPC project.
 * (2026) matthew.horton@uknnl.com
 * 
 * This file is modified from files part of the MACHDYN package for LAMMPS.
 * Copyright (2014) Georg C. Ganzenmueller, georg.ganzenmueller@emi.fhg.de
 * Fraunhofer Ernst-Mach Institute for High-Speed Dynamics, EMI,
 * Eckerstrasse 4, D-79104 Freiburg i.Br, Germany.
 *
 * ----------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifdef ATOM_CLASS
// clang-format off
AtomStyle(spooff,AtomVecSPOOFF);
// clang-format on
#else

#ifndef LMP_ATOM_VEC_SPOOFF_H
#define LMP_ATOM_VEC_SPOOFF_H

#include "atom_vec.h"

namespace LAMMPS_NS {

class AtomVecSPOOFF : virtual public AtomVec {
 public:
  AtomVecSPOOFF(class LAMMPS *);

  void grow_pointers() override;
  void force_clear(int, size_t) override;
  void create_atom_post(int) override;
  void data_atom_post(int) override;

 private:
  tagint *molecule;
  double *tsph, *rhosph, *pesph, *dpesph, **qsph, *cpsph, *ksph, *hgsph, *fracsph, *yeild, *esph, *desph, *vfrac, *rmass, *radius, *contact_radius;
  double *eff_plastic_strain, *eff_plastic_strain_rate, *damage, *starting_neighs, *youngs, *poissons, *linear_expansion;
  double **x0, **smd_data_9, **smd_stress, **vest;
};

}    // namespace LAMMPS_NS

#endif
#endif
