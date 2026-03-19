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

/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifdef COMPUTE_CLASS
// clang-format off
ComputeStyle(damage/local,ComputeDamageLocal);
// clang-format on
#else

#ifndef LMP_COMPUTE_DAMAGE_LOCAL_H
#define LMP_COMPUTE_DAMAGE_LOCAL_H

#include "compute.h"

namespace LAMMPS_NS {

class ComputeDamageLocal : public Compute {
 public:
  ComputeDamageLocal(class LAMMPS *, int, char **);
  ~ComputeDamageLocal() override;
  void init() override;
  void init_list(int, class NeighList *) override;
  void compute_local() override;
  double memory_usage() override;

 protected:
  int ifix_tlsph;
  class FixSPOOFF_TLSPH_ReferenceConfiguration *fix_tlsph_reference_configuration;

 private:
  int nvalues, kindflag, cutstyle;

  int nmax;
  double *vlocal;
  double **alocal;
  double *buf;

  class NeighList *list;

  int ncount;
  int **indices;

  int count_pairs(int, int);
  void reallocate(int);

  typedef void (ComputeDamageLocal::*FnPtrPack)(int);
  FnPtrPack *pack_choice;    // ptrs to pack functions

  void pack_patom1(int);
  void pack_patom2(int);
  void pack_damage(int);
  void pack_degredation(int);
  void pack_energy(int);

};

}    // namespace LAMMPS_NS

#endif
#endif
