/* -*- c++ -*- ----------------------------------------------------------
 *
 *                    *** Smooth Mach Dynamics ***
 *
 * This file is part of the MACHDYN package for LAMMPS.
 * Copyright (2014) Georg C. Ganzenmueller, georg.ganzenmueller@emi.fhg.de
 * Fraunhofer Ernst-Mach Institute for High-Speed Dynamics, EMI,
 * Eckerstrasse 4, D-79104 Freiburg i.Br, Germany.
 *
 * This file is based on the FixShearHistory class.
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

#ifdef FIX_CLASS
// clang-format off
FixStyle(SPOOFF_TLSPH_NEIGHBORS,FixSPOOFF_TLSPH_ReferenceConfiguration);
// clang-format on
#else

#ifndef LMP_FIX_SPOOFF_TLSPH_REFERENCE_H
#define LMP_FIX_SPOOFF_TLSPH_REFERENCE_H

#include "fix.h"

namespace LAMMPS_NS {

class FixSPOOFF_TLSPH_ReferenceConfiguration : public Fix {
  friend class Neighbor;
  friend class PairTlsph;
  friend class FixSPOOFFSetConvection;
  friend class ComputeSPOOFFDamage;
  friend class ComputeDamageLocal;

 public:
  FixSPOOFF_TLSPH_ReferenceConfiguration(class LAMMPS *, int, char **);
  ~FixSPOOFF_TLSPH_ReferenceConfiguration() override;

  int setmask() override;
  void init() override;
  void setup(int) override;
  void pre_exchange() override;
  int pack_forward_comm(int, int *, double *, int, int *) override;
  void unpack_forward_comm(int, int, double *) override;

  double memory_usage() override;
  void grow_arrays(int) override;
  void copy_arrays(int, int, int) override;
  int pack_exchange(int, double *) override;
  int unpack_exchange(int, double *) override;
  int pack_restart(int, double *) override;
  void unpack_restart(int, int) override;
  int size_restart(int) override;
  int maxsize_restart() override;

  bool crack_exclude(int i, int j);
  bool get_line_intersection(int i, int j);

 protected:
  int updateFlag;    // flag to update reference configuration
  int nmax;
  int maxpartner;
  int *npartner;       // Num of touching partners of each atom
  tagint **partner;    // global atom IDs for the partners
  tagint **partnerfirst;    // original global atom IDs for the partners
  float **wfd_list, **wf_list, **energy_per_bond;
  float **degradation_ij;    // per-pair interaction degradation status
  float **impact_velocity;    // per-pair interaction impact velocity storage

  class Pair *pair;
};

}    // namespace LAMMPS_NS

#endif
#endif
