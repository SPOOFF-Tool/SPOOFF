/* -*- c++ -*- ----------------------------------------------------------
 *
 *                    *** Smooth Mach Dynamics ***
 *
 * This file is part of the MACHDYN package for LAMMPS.
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

#ifdef COMPUTE_CLASS
// clang-format off
ComputeStyle(spooff/damage,ComputeSPOOFFDamage);
// clang-format on
#else

#ifndef LMP_COMPUTE_SPOOFF_DAMAGE_H
#define LMP_COMPUTE_SPOOFF_DAMAGE_H

#include "compute.h"
#include "fix.h"

namespace LAMMPS_NS {

class ComputeSPOOFFDamage : public Compute {
 public:
  ComputeSPOOFFDamage(class LAMMPS *, int, char **);
  ~ComputeSPOOFFDamage() override;
  void init() override;
  void compute_peratom() override;
  double memory_usage() override;

 protected:
  int ifix_tlsph;
  class FixSPOOFF_TLSPH_ReferenceConfiguration *fix_tlsph_reference_configuration;

 private:
  int nmax;
  double **damage_vector;
};

}    // namespace LAMMPS_NS

#endif
#endif
