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

#ifdef FIX_CLASS
// clang-format off
FixStyle(spooff/setconvection,FixSPOOFFSetConvection);
// clang-format on
#else

#ifndef LMP_FIX_SPOOFF_SETCONVECTION_H
#define LMP_FIX_SPOOFF_SETCONVECTION_H

#include "fix.h"

namespace LAMMPS_NS {

class FixSPOOFFSetConvection : public Fix {
 public:
  FixSPOOFFSetConvection(class LAMMPS *, int, char **);
  ~FixSPOOFFSetConvection() override;
  int setmask() override;
  void init() override;
  void setup(int) override;
  void min_setup(int) override;
  //void initial_integrate(int);
  void post_force(int) override;
  double compute_vector(int) override;
  double memory_usage() override;

 protected:
  int ifix_tlsph;
  class FixSPOOFF_TLSPH_ReferenceConfiguration *fix_tlsph_reference_configuration;

 private:
  double hvalue, T_inf;
  int varflag;
  char *estr;
  char *idregion;
  class Region *region;
  int evar, estyle;
  double foriginal[3], foriginal_all[3];
  int force_flag;

  int maxatom;
  double **sforce;
};

}    // namespace LAMMPS_NS

#endif
#endif
