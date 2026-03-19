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

#include "compute_spooff_tlsph_strain_cylindrical.h"

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
using namespace std;
using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

ComputeSPOOFFTLSPHStrainCylindrical::ComputeSPOOFFTLSPHStrainCylindrical(LAMMPS *lmp, int narg, char **arg) :
                Compute(lmp, narg, arg) {
        if (narg != 3)
                error->all(FLERR, "Illegal compute spooff/tlsph_strain_cylindrical command");

        peratom_flag = 1;
        size_peratom_cols = 6;

        nmax = 0;
        strain_cylindrical_Vector = nullptr;
}

/* ---------------------------------------------------------------------- */

ComputeSPOOFFTLSPHStrainCylindrical::~ComputeSPOOFFTLSPHStrainCylindrical() {
        memory->sfree(strain_cylindrical_Vector);
}

/* ---------------------------------------------------------------------- */

void ComputeSPOOFFTLSPHStrainCylindrical::init() {

        int count = 0;
        for (int i = 0; i < modify->ncompute; i++)
                if (strcmp(modify->compute[i]->style, "spooff/tlsph_strain_cylindrical") == 0)
                        count++;
        if (count > 1 && comm->me == 0)
                error->warning(FLERR, "More than one compute spooff/tlsph_strain_cylindrical");
}

/* ---------------------------------------------------------------------- */

void ComputeSPOOFFTLSPHStrainCylindrical::compute_peratom() {
        double **defgrad0 = atom->smd_data_9;

        invoked_peratom = update->ntimestep;

        // grow Vector array if necessary

        if (atom->nmax > nmax) {
                memory->destroy(strain_cylindrical_Vector);
                nmax = atom->nmax;
                memory->create(strain_cylindrical_Vector, nmax, size_peratom_cols, "strain_cylindrical_Vector");
                array_atom = strain_cylindrical_Vector;
        }

        // copy data to output array
        int itmp = 0;
        auto Fincr = (Matrix3d *) force->pair->extract("spooff/tlsph/Fincr_ptr", itmp);
        if (Fincr == nullptr) {
                error->all(FLERR, "compute spooff/tlsph_strain_cylindrical failed to access Fincr array");
        }

        int *mask = atom->mask;
        int nlocal = atom->nlocal;
        double **x = atom->x;
        double theta;
        Matrix3d E, eye, Ftotal, F0;
        eye.setIdentity();

        for (int i = 0; i < nlocal; i++) {
                if (mask[i] & groupbit) {

                        // old deformation gradient
                        F0(0, 0) = defgrad0[i][0];
                        F0(0, 1) = defgrad0[i][1];
                        F0(0, 2) = defgrad0[i][2];
                        F0(1, 0) = defgrad0[i][3];
                        F0(1, 1) = defgrad0[i][4];
                        F0(1, 2) = defgrad0[i][5];
                        F0(2, 0) = defgrad0[i][6];
                        F0(2, 1) = defgrad0[i][7];
                        F0(2, 2) = defgrad0[i][8];

                        // compute current total deformation gradient
                        Ftotal = F0 * Fincr[i]; // this is the total deformation gradient: reference deformation times incremental deformation

                        theta = atan2(x[i][1],x[i][0]);
                        E = 0.5 * (Ftotal.transpose() * Ftotal - eye); // Green-Lagrange strain
                        strain_cylindrical_Vector[i][0] = E(0, 0)*pow(cos(theta),2) + E(1, 1)*pow(sin(theta),2) + 2*E(0, 1)*sin(theta)*cos(theta); // rr
                        strain_cylindrical_Vector[i][1] = E(0, 0)*pow(sin(theta),2) + E(1, 1)*pow(cos(theta),2) - 2*E(0, 1)*sin(theta)*cos(theta); // tt
                        strain_cylindrical_Vector[i][2] = E(2, 2); // zz
                        strain_cylindrical_Vector[i][3] = (E(1, 1)-E(0, 0))*sin(theta)*cos(theta) + E(0, 1)*(pow(cos(theta),2)-pow(sin(theta),2)); // rt
                        strain_cylindrical_Vector[i][4] = E(0, 2)*cos(theta) + E(1, 2)*sin(theta); // rz
                        strain_cylindrical_Vector[i][5] = -E(0, 2)*sin(theta) + E(1, 2)*cos(theta); // tz
                } else {
                        for (int j = 0; j < size_peratom_cols; j++) {
                                strain_cylindrical_Vector[i][j] = 0.0;
                        }
                }
        }
}

/* ----------------------------------------------------------------------
 memory usage of local atom-based array
 ------------------------------------------------------------------------- */

double ComputeSPOOFFTLSPHStrainCylindrical::memory_usage() {
        double bytes = (double)size_peratom_cols * nmax * sizeof(double);
        return bytes;
}
