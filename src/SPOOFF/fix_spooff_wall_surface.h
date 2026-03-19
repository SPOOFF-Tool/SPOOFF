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
FixStyle(spooff/wall_surface,FixSPOOFFWallSurface);
// clang-format on
#else

#ifndef LMP_FIX_SPOOFF_WALL_SURFACE_H
#define LMP_FIX_SPOOFF_WALL_SURFACE_H

#include "fix.h"

namespace LAMMPS_NS {

class FixSPOOFFWallSurface : public Fix {

 public:
  FixSPOOFFWallSurface(class LAMMPS *, int, char **);
  ~FixSPOOFFWallSurface() override;

  int setmask() override;
  void init() override;
  void setup(int) override;
  void min_setup(int) override;

  void read_triangles(int pass);

 private:
  int first;                    // flag for first time initialization
  double sublo[3], subhi[3];    // epsilon-extended proc sub-box for adding atoms;
  char *filename;
  int wall_particle_type;
  int wall_molecule_id;
};
}    // namespace LAMMPS_NS

#endif
#endif
