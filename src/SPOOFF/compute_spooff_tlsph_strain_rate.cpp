// clang-format off
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

#include "compute_spooff_tlsph_strain_rate.h"

#include "atom.h"
#include "comm.h"
#include "error.h"
#include "force.h"
#include "memory.h"
#include "modify.h"
#include "pair.h"
#include "update.h"

#include <cstring>
#include <Eigen/Eigen>          // IWYU pragma: export

using namespace Eigen;
using namespace LAMMPS_NS;


/* ---------------------------------------------------------------------- */

ComputeSPOOFFTLSPHStrainRate::ComputeSPOOFFTLSPHStrainRate(LAMMPS *lmp, int narg, char **arg) :
                Compute(lmp, narg, arg) {
        if (narg != 3)
                error->all(FLERR, "Illegal compute spooff/ulsph_strain_rate command");

        peratom_flag = 1;
        size_peratom_cols = 6;

        nmax = 0;
        strain_rate_array = nullptr;
}

/* ---------------------------------------------------------------------- */

ComputeSPOOFFTLSPHStrainRate::~ComputeSPOOFFTLSPHStrainRate() {
        memory->sfree(strain_rate_array);
}

/* ---------------------------------------------------------------------- */

void ComputeSPOOFFTLSPHStrainRate::init() {

        int count = 0;
        for (int i = 0; i < modify->ncompute; i++)
                if (strcmp(modify->compute[i]->style, "spooff/ulsph_strain_rate") == 0)
                        count++;
        if (count > 1 && comm->me == 0)
                error->warning(FLERR, "More than one compute spooff/ulsph_strain_rate");
}

/* ---------------------------------------------------------------------- */

void ComputeSPOOFFTLSPHStrainRate::compute_peratom() {
        invoked_peratom = update->ntimestep;

        // grow vector array if necessary

        if (atom->nmax > nmax) {
                memory->destroy(strain_rate_array);
                nmax = atom->nmax;
                memory->create(strain_rate_array, nmax, size_peratom_cols, "stresstensorVector");
                array_atom = strain_rate_array;
        }

        int itmp = 0;
        auto D = (Matrix3d *) force->pair->extract("spooff/tlsph/strain_rate_ptr", itmp);
        if (D == nullptr) {
                error->all(FLERR,
                                "compute spooff/tlsph_strain_rate could not access strain rate. Are the matching pair styles present?");
        }

        int nlocal = atom->nlocal;

        for (int i = 0; i < nlocal; i++) {

                strain_rate_array[i][0] = D[i](0, 0); // xx
                strain_rate_array[i][1] = D[i](1, 1); // yy
                strain_rate_array[i][2] = D[i](2, 2); // zz
                strain_rate_array[i][3] = D[i](0, 1); // xy
                strain_rate_array[i][4] = D[i](0, 2); // xz
                strain_rate_array[i][5] = D[i](1, 2); // yz
        }
}

/* ----------------------------------------------------------------------
 memory usage of local atom-based array
 ------------------------------------------------------------------------- */

double ComputeSPOOFFTLSPHStrainRate::memory_usage() {
        double bytes = (double)size_peratom_cols * nmax * sizeof(double);
        return bytes;
}
