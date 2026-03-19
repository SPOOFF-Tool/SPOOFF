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

#include "pair_spooff_tlsph.h"

#include "fix_spooff_tlsph_reference_configuration.h"

#include "atom.h"
#include "comm.h"
#include "domain.h"
#include "error.h"
#include "fix.h"
#include "force.h"
#include "group.h"
#include "memory.h"
#include "modify.h"
#include "neighbor.h"
#include "spooff_kernels.h"
#include "spooff_material_models.h"
#include "spooff_math.h"
#include "update.h"

#include <cmath>
#include <cstring>
#include <Eigen/Eigen>

using namespace SPOOFF_Kernels;
using namespace Eigen;
using namespace LAMMPS_NS;
using namespace SPOOFF_Math;

#define JAUMANN false
#define DETF_MIN 0 // maximum compression deformation allow was 0.2
#define DETF_MAX 10.0 // maximum tension deformation allowed  was 2
#define TLSPH_DEBUG 0
#define PLASTIC_STRAIN_AVERAGE_WINDOW 100.0

/* ---------------------------------------------------------------------- */

PairTlsph::PairTlsph(LAMMPS *lmp) :
  Pair(lmp) {

  onerad_dynamic = onerad_frozen = maxrad_dynamic = maxrad_frozen = nullptr;

  failureModel = nullptr;
  strengthModel = eos = nullptr;
  contact_mechanical = contact_thermal = nullptr;

  nmax = 0; // make sure no atom on this proc such that initial memory allocation is correct
  Fdot = Fincr = K = PK1 = nullptr;
  R = FincrInv = W = D = nullptr;
  detF = nullptr;
  smoothVelDifference = nullptr;
  numNeighsRefConfig = nullptr;
  CauchyStress = nullptr;
  hourglass_error = nullptr;
  Lookup = nullptr;
  particle_dt = nullptr;

  updateFlag = 0;
  first = true;
  second = true;
  dtCFL = 0.0; // initialize dtCFL so it is set to safe value if extracted on zero-th timestep

  comm_forward = 25; // this pair style communicates 20 doubles to ghost atoms : PK1 tensor + F tensor + shepardWeight   //MH added cpsph,ksph,hgsph,fracsph to forward comm
  fix_tlsph_reference_configuration = nullptr;

  cut_comm = MAX(neighbor->cutneighmax, comm->cutghostuser); // cutoff radius within which ghost atoms are communicated.
}

/* ---------------------------------------------------------------------- */

PairTlsph::~PairTlsph() {

  if (allocated) {
    memory->destroy(setflag);
    memory->destroy(cutsq);
    memory->destroy(strengthModel);
    memory->destroy(contact_thermal);
    memory->destroy(contact_mechanical);
    memory->destroy(eos);
    memory->destroy(Lookup);

    delete[] onerad_dynamic;
    delete[] onerad_frozen;
    delete[] maxrad_dynamic;
    delete[] maxrad_frozen;

    delete[] Fdot;
    delete[] Fincr;
    delete[] K;
    delete[] detF;
    delete[] PK1;
    delete[] smoothVelDifference;
    delete[] R;
    delete[] FincrInv;
    delete[] W;
    delete[] D;
    delete[] numNeighsRefConfig;
    delete[] CauchyStress;
    delete[] hourglass_error;
    delete[] particle_dt;

    delete[] failureModel;
  }
}

/* ----------------------------------------------------------------------
 *
 * use half neighbor list to re-compute shape matrix
 *
 ---------------------------------------------------------------------- */

void PairTlsph::PreCompute() {
  tagint *mol = atom->molecule;
  double *vfrac = atom->vfrac;
  double *radius = atom->radius;
  double **x0 = atom->x0;
  double **x = atom->x;
  double **v = atom->vest; // extrapolated velocities corresponding to current positions
  double **vint = atom->v; // Velocity-Verlet algorithm velocities
  double *damage = atom->damage;
  double *cpsph = atom->cpsph;
  double *ksph = atom->ksph;
  double *hgsph = atom->hgsph;
  double *fracsph = atom->fracsph;
  double **qsph = atom->qsph;
  double *tsph = atom->tsph;
  double *esph = atom->esph;
  double *rmass = atom->rmass;
  tagint *tag = atom->tag;
  int *type = atom->type;
  int nlocal = atom->nlocal;
  int jnum, jj, i, j, itype, idim;
  

  tagint **partner = (dynamic_cast<FixSPOOFF_TLSPH_ReferenceConfiguration *>(modify->fix[ifix_tlsph]))->partner;
  int *npartner = (dynamic_cast<FixSPOOFF_TLSPH_ReferenceConfiguration *>(modify->fix[ifix_tlsph]))->npartner;
  float **wfd_list = (dynamic_cast<FixSPOOFF_TLSPH_ReferenceConfiguration *>(modify->fix[ifix_tlsph]))->wfd_list;
  float **wf_list = (dynamic_cast<FixSPOOFF_TLSPH_ReferenceConfiguration *>(modify->fix[ifix_tlsph]))->wf_list;
  float **degradation_ij = (dynamic_cast<FixSPOOFF_TLSPH_ReferenceConfiguration *>(modify->fix[ifix_tlsph]))->degradation_ij;
  double r0, r0Sq, wf, wfd, h, irad, voli, volj, scale, shepardWeight;
  Vector3d dx, dx0, dv, g;
  Matrix3d Ktmp, Ftmp, Fdottmp, L, U, eye;
  Vector3d vi, vj, vinti, vintj, xi, xj, x0i, x0j, dvint;
  int periodic = (domain->xperiodic || domain->yperiodic || domain->zperiodic);
  bool status;
  Matrix3d F0;

  eye.setIdentity();

  for (i = 0; i < nlocal; i++) {

    itype = type[i];
    if (setflag[itype][itype] == 1) {

    //MH I've added these so that they can be varied per atom either by weibul or by temperature dependence 
    //cpsph[i] = Lookup[HEAT_CAPACITY][itype];
    //ksph[i] = Lookup[THERMAL_CONDUCTIVITY][itype];
    //fracsph[i] = Lookup[FAILURE_MAX_PAIRWISE_STRAIN_THRESHOLD][itype];   //moved to initialisation function only

    //tsph[i] = esph[i]/(cpsph[i]*rmass[i]);

    qsph[i][0] = 0.0; //set up heatflux ready to be computed
    qsph[i][1] = 0.0;
    qsph[i][2] = 0.0;

    //--------------------------

      K[i].setZero();
      Fincr[i].setZero();
      Fdot[i].setZero();
      numNeighsRefConfig[i] = 0;
      smoothVelDifference[i].setZero();
      hourglass_error[i] = 0.0;

      if (mol[i] < 0) { // valid SPH particle have mol > 0
        continue;
      }

      // initialize aveage mass density
      h = 2.0 * radius[i];
      r0 = 0.0;
      spiky_kernel_and_derivative(h, r0, domain->dimension, wf, wfd);

      jnum = npartner[i];
      irad = radius[i];
      voli = vfrac[i];
      shepardWeight = wf * voli;

      // initialize Eigen data structures from LAMMPS data structures
      for (idim = 0; idim < 3; idim++) {
        xi(idim) = x[i][idim];
        x0i(idim) = x0[i][idim];
        vi(idim) = v[i][idim];
        vinti(idim) = vint[i][idim];
      }

      for (jj = 0; jj < jnum; jj++) {

        if (partner[i][jj] == 0)
          continue;
        j = atom->map(partner[i][jj]);
        if (j < 0) { //                 // check if lost a partner without first breaking bond
          partner[i][jj] = 0;
          continue;
        }

        if (mol[j] < 0) { // particle has failed. do not include it for computing any property
          continue;
        }

        if (mol[i] != mol[j]) {
          continue;
        }

        // initialize Eigen data structures from LAMMPS data structures
        for (idim = 0; idim < 3; idim++) {
          xj(idim) = x[j][idim];
          x0j(idim) = x0[j][idim];
          vj(idim) = v[j][idim];
          vintj(idim) = vint[j][idim];
        }
        dx0 = x0j - x0i;
        dx = xj - xi;

        if (periodic)
          domain->minimum_image(dx0(0), dx0(1), dx0(2));

        r0Sq = dx0.squaredNorm();
        h = irad + radius[j];

        r0 = sqrt(r0Sq);
        volj = vfrac[j];

        // distance vectors in current and reference configuration, velocity difference
        dv = vj - vi;
        dvint = vintj - vinti;

        // scale the interaction according to the damage variable
        scale = 1.0 - degradation_ij[i][jj];
        wf = wf_list[i][jj] * scale;
        wfd = wfd_list[i][jj] * scale;
        g = (wfd / r0) * dx0;

        

        /* build matrices */
        Ktmp = -g * dx0.transpose();
        Fdottmp = -dv * g.transpose();
        Ftmp = -(dx - dx0) * g.transpose();

        K[i] += volj * Ktmp;
        Fdot[i] += volj * Fdottmp;
        Fincr[i] += volj * Ftmp;
        shepardWeight += volj * wf;
        smoothVelDifference[i] += volj * wf * dvint;
        numNeighsRefConfig[i]++;
      } // end loop over j

      // normalize average velocity field around an integration point
      if (shepardWeight > 0.0) {
        smoothVelDifference[i] /= shepardWeight;
      } else {
        smoothVelDifference[i].setZero();
      }

      pseudo_inverse_SVD(K[i]);
      Fdot[i] *= K[i];
      Fincr[i] *= K[i];
      Fincr[i] += eye;

      if (JAUMANN) {
        R[i].setIdentity(); // for Jaumann stress rate, we do not need a subsequent rotation back into the reference configuration
      } else {
        status = PolDec(Fincr[i], R[i], U, false); // polar decomposition of the deformation gradient, F = R * U
        if (!status) {
          error->message(FLERR, "Polar decomposition of deformation gradient failed.\n");
          mol[i] = -1;
        } else {
          Fincr[i] = R[i] * U;
        }
      }

      if(Fincr[i]==eye) {
        detF[i] = 1;
      } else {
        detF[i] = Fincr[i].determinant();
      }
      
      FincrInv[i] = Fincr[i].inverse();

      // velocity gradient
      L = Fdot[i] * FincrInv[i];

      // symmetric (D) and asymmetric (W) parts of L
      D[i] = 0.5 * (L + L.transpose());
      W[i] = 0.5 * (L - L.transpose()); // spin tensor:: need this for Jaumann rate

      // unrotated rate-of-deformation tensor d, see right side of Pronto2d, eqn.(2.1.7)
      // convention: unrotated frame is that one, where the true rotation of an integration point has been subtracted.
      // stress in the unrotated frame of reference is denoted sigma (stress seen by an observer doing rigid body rotations along with the material)
      // stress in the true frame of reference (a stationary observer) is denoted by T, "true stress"
      D[i] = (R[i].transpose() * D[i] * R[i]).eval();

      // limit strain rate
      //double limit = 1.0e-3 * Lookup[SIGNAL_VELOCITY][itype] / radius[i];
      //D[i] = LimitEigenvalues(D[i], limit);

      /*
       * make sure F stays within some limits
       */

      if ((detF[i] < DETF_MIN) || (detF[i] > DETF_MAX) || (numNeighsRefConfig[i] == 0)) {
        utils::logmesg(lmp, "deleting particle [{}] because det(F)={}f is outside stable range"
                       " {} -- {} \n", tag[i], Fincr[i].determinant(), DETF_MIN, DETF_MAX);
                       std::string sep = "\n----------------------------------------\n";
                       std::cout << Fincr[i] << sep;
        utils::logmesg(lmp,"nn = {}, damage={}\n", numNeighsRefConfig[i], damage[i]);
        mol[i] = -1;
      }

      if (mol[i] < 0) {
        D[i].setZero();
        Fdot[i].setZero();
        Fincr[i].setIdentity();
        smoothVelDifference[i].setZero();
        detF[i] = 1.0;
        K[i].setIdentity();

        vint[i][0] = 0.0;
        vint[i][1] = 0.0;
        vint[i][2] = 0.0;
      }
    } // end loop over i
  } // end check setflag
}

/* ---------------------------------------------------------------------- */

void PairTlsph::compute(int eflag, int vflag) {
  double *ksph = atom->ksph;
  double *hgsph = atom->hgsph;
  double *cpsph = atom->cpsph;
  double *fracsph = atom->fracsph;
  double *esph = atom->esph;
  double *tsph = atom->tsph;
  double *rhosph = atom->rhosph;
  double *vfrac = atom->vfrac;
  double *rmass = atom->rmass;
  double *starting_neighs = atom->starting_neighs;
  double *youngs = atom->youngs;
  double *yeild = atom->yeild;
  double *poissons = atom->poissons;
  double *linear_expansion = atom->linear_expansion;
  int *type = atom->type;

  if (atom->nmax > nmax) {
    nmax = atom->nmax;
    delete[] Fdot;
    Fdot = new Matrix3d[nmax]; // memory usage: 9 doubles
    delete[] Fincr;
    Fincr = new Matrix3d[nmax]; // memory usage: 9 doubles
    delete[] K;
    K = new Matrix3d[nmax]; // memory usage: 9 doubles
    delete[] PK1;
    PK1 = new Matrix3d[nmax]; // memory usage: 9 doubles; total 5*9=45 doubles
    delete[] detF;
    detF = new double[nmax]; // memory usage: 1 double; total 46 doubles
    delete[] smoothVelDifference;
    smoothVelDifference = new Vector3d[nmax]; // memory usage: 3 doubles; total 49 doubles
    delete[] R;
    R = new Matrix3d[nmax]; // memory usage: 9 doubles; total 67 doubles
    delete[] FincrInv;
    FincrInv = new Matrix3d[nmax]; // memory usage: 9 doubles; total 85 doubles
    delete[] W;
    W = new Matrix3d[nmax]; // memory usage: 9 doubles; total 94 doubles
    delete[] D;
    D = new Matrix3d[nmax]; // memory usage: 9 doubles; total 103 doubles
    delete[] numNeighsRefConfig;
    numNeighsRefConfig = new int[nmax]; // memory usage: 1 int; total 108 doubles
    delete[] CauchyStress;
    CauchyStress = new Matrix3d[nmax]; // memory usage: 9 doubles; total 118 doubles
    delete[] hourglass_error;
    hourglass_error = new double[nmax];
    delete[] particle_dt;
    particle_dt = new double[nmax];
  }

  if (first) { // return on first call, because reference connectivity lists still needs to be built. Also zero quantities which are otherwise undefined.
    first = false;

    for (int i = 0; i < atom->nlocal; i++) {
      Fincr[i].setZero();
      detF[i] = 0.0;
      smoothVelDifference[i].setZero();
      D[i].setZero();
      numNeighsRefConfig[i] = 0;
      starting_neighs[i] = 1.0;  //protects for divide by zero erros
      CauchyStress[i].setZero();
      hourglass_error[i] = 0.0;
      particle_dt[i] = 0.0;
      //MH I've added these so that they can be varied per atom either by weibul or by temperature dependence 
      cpsph[i] = Lookup[HEAT_CAPACITY][type[i]];
      ksph[i] = Lookup[THERMAL_CONDUCTIVITY][type[i]];
      hgsph[i] = Lookup[HEAT_PRODUCED][type[i]];
      yeild[i] = Lookup[YIELD_STRESS][type[i]];
      if (failureModel[type[i]].failure_max_pairwise_strain) {
        fracsph[i] = Lookup[FAILURE_MAX_PAIRWISE_STRAIN_THRESHOLD][type[i]];
      }else if (failureModel[type[i]].failure_max_plastic_strain) {
        fracsph[i] = Lookup[FAILURE_MAX_PLASTIC_STRAIN_THRESHOLD][type[i]];
      }else if (failureModel[type[i]].failure_max_principal_stress) {
        fracsph[i] = Lookup[FAILURE_MAX_PRINCIPAL_STRESS_THRESHOLD][type[i]]/Lookup[YOUNGS_MODULUS][type[i]];
      }else if (failureModel[type[i]].failure_max_principal_strain) {
        fracsph[i] = Lookup[FAILURE_MAX_PRINCIPAL_STRAIN_THRESHOLD][type[i]];
      }
      tsph[i] = esph[i]/(cpsph[i]*rmass[i]);
      youngs[i] = Lookup[YOUNGS_MODULUS][type[i]];
      poissons[i] = Lookup[POISSON_RATIO][type[i]];
      linear_expansion[i] = Lookup[EOS_THERMAL_ALPHA][type[i]];
      rhosph[i]=rmass[i]/vfrac[i];
    }


    ///MH Add Weibull distrubution here

    return;
  }

  /*
   * calculate deformations and rate-of-deformations
   */
  PairTlsph::PreCompute();

  if (second) { // return on first call, because reference connectivity lists still needs to be built. Also zero quantities which are otherwise undefined.
    second = false;

    for (int i = 0; i < atom->nlocal; i++) {
      starting_neighs[i] = numNeighsRefConfig[i];
    }
  }

  /*
   * calculate stresses from constitutive models
   */
  PairTlsph::AssembleStress();

  /*
   * QUANTITIES ABOVE HAVE ONLY BEEN CALCULATED FOR NLOCAL PARTICLES.
   * NEED TO DO A FORWARD COMMUNICATION TO GHOST ATOMS NOW
   */
  comm->forward_comm(this);

  /*
   * compute forces between particles
   */
  updateFlag = 0;
  ComputeForces(eflag, vflag);
}

void PairTlsph::ComputeForces(int eflag, int vflag) {
  tagint *mol = atom->molecule;
  double **x = atom->x;
  double **v = atom->vest;
  double **x0 = atom->x0;
  double **f = atom->f;
  double *vfrac = atom->vfrac;
  double *desph = atom->desph;
  double *esph = atom->esph;
  double *pesph = atom->pesph;
  double *dpesph = atom->dpesph;
  double *tsph = atom->tsph;
  double **qsph = atom->qsph;
  double *ksph = atom->ksph;
  double *hgsph = atom->hgsph;
  double *cpsph = atom->cpsph;
  double *rhosph = atom->rhosph;
  double *fracsph = atom->fracsph;
  double *rmass = atom->rmass;
  double *radius = atom->radius;
  double *contact_radius = atom->contact_radius;
  double *damage = atom->damage;
  double *youngs = atom->youngs;
  double *poissons = atom->poissons;
  double *linear_expansion = atom->linear_expansion;
  double *plastic_strain = atom->eff_plastic_strain;
  int *type = atom->type;
  int nlocal = atom->nlocal;
  int i, j, jj, jnum, itype, idim;
  double r, hg_mag, wf, wfd, h, r0, r0Sq, voli, volj;
  double delVdotDelR, visc_magnitude, deltaE, deltaHE, TotaldeltaHE, deltaQ, mu_ij, hg_err, gamma_dot_dx, delta, scale;
  double strain1d, strain1d_max, softening_strain, shepardWeight;
  Vector3d fi, fj, dx0, dx, dv, f_stress, f_hg, dxp_i, dxp_j, gamma, g, gamma_i, gamma_j, x0i, x0j;
  Vector3d xi, xj, vi, vj, f_visc, sumForces, f_spring;
  int periodic = (domain->xperiodic || domain->yperiodic || domain->zperiodic);

  tagint **partner = (dynamic_cast<FixSPOOFF_TLSPH_ReferenceConfiguration *>(modify->fix[ifix_tlsph]))->partner;
  tagint **partnerfirst = (dynamic_cast<FixSPOOFF_TLSPH_ReferenceConfiguration *>(modify->fix[ifix_tlsph]))->partnerfirst;
  int *npartner = (dynamic_cast<FixSPOOFF_TLSPH_ReferenceConfiguration *>(modify->fix[ifix_tlsph]))->npartner;
  float **wfd_list = (dynamic_cast<FixSPOOFF_TLSPH_ReferenceConfiguration *>(modify->fix[ifix_tlsph]))->wfd_list;
  float **wf_list = (dynamic_cast<FixSPOOFF_TLSPH_ReferenceConfiguration *>(modify->fix[ifix_tlsph]))->wf_list;
  float **degradation_ij = (dynamic_cast<FixSPOOFF_TLSPH_ReferenceConfiguration *>(modify->fix[ifix_tlsph]))->degradation_ij;
  float **impact_velocity = (dynamic_cast<FixSPOOFF_TLSPH_ReferenceConfiguration *>(modify->fix[ifix_tlsph]))->impact_velocity;
  float **energy_per_bond = (dynamic_cast<FixSPOOFF_TLSPH_ReferenceConfiguration *>(modify->fix[ifix_tlsph]))->energy_per_bond;
  Matrix3d eye;
  eye.setIdentity();

  ev_init(eflag, vflag);

  /*
   * iterate over pairs of particles i, j and assign forces using PK1 stress tensor
   */

  //updateFlag = 0;
  hMin = 1.0e22;
  dtRelative = 1.0e22;

  for (i = 0; i < nlocal; i++) {

    if (mol[i] < 0) {
      continue; // Particle i is not a valid SPH particle (anymore). Skip all interactions with this particle.
    }

    if(false){
          printf("--------------------------");  //added for testing
    printf("i: %d \n",i); //added for testing
    printf("cpsph: %f \n",cpsph[i]); //added for testing
    printf("ksph: %f \n",ksph[i]); //added for testing
    printf("fracsph: %f \n",fracsph[i]); //added for testing
    printf("qsph: %f \n",qsph[i][0]);  //added for testing
    printf("tsph: %f \n",tsph[i]);  //added for testing
    }

    itype = type[i];
    jnum = npartner[i];
    voli = vfrac[i];

    // initialize average mass density
    h = 2.0 * radius[i];
    r = 0.0;
    TotaldeltaHE = 0.0;
    spiky_kernel_and_derivative(h, r, domain->dimension, wf, wfd);
    shepardWeight = wf * voli;

    for (idim = 0; idim < 3; idim++) {
      x0i(idim) = x0[i][idim];
      xi(idim) = x[i][idim];
      vi(idim) = v[i][idim];
    }

    for (jj = 0; jj < jnum; jj++) {
      if (partner[i][jj] == 0){
        if(contact_mechanical[itype]){ 
        if(partnerfirst[i][jj] != 0){
          //hertz contact force for disconected material
          j = atom->map(partnerfirst[i][jj]);
          if (mol[j] < 0) {
            continue; // Particle j is not a valid SPH particle (anymore). Skip all interactions with this particle.
          }

          if (mol[i] != mol[j]) {
            continue;
          }

          if (type[j] != itype)
            error->all(FLERR, "particle pair is not of same type!");
            
          double ri, rj, rsq, rcut, rcutSq, r_geom, delta, fpair, evdwl, scale;
          scale = Lookup[CONTACT_DISTANCE_SCALING][itype];
          for (idim = 0; idim < 3; idim++) {
          x0j(idim) = x0[j][idim];
          xj(idim) = x[j][idim];
          vj(idim) = v[j][idim];
          }
          dx = xj - xi;
          dv = vj - vi;
          r = dx.norm(); // current distance

          if(contact_mechanical[itype]){
          //check for contact
          ri = scale * contact_radius[i];
          rj = scale * contact_radius[j];
          rsq = dx.squaredNorm();
          rcut = ri + rj;
          rcutSq = rcut * rcut;
          if (rsq >= rcutSq) {
            impact_velocity[i][jj] = -1.0;
            continue;
          }
          // Hertzian short-range forces
                                delta = rcut - r; // overlap distance
                                r_geom = ri * rj / rcut;
                                double K_eff, mu_eff, shear, bulk, m_modulus, lame_lambda, E_reduced;
                                        shear = youngs[i] / (2.0 * (1.0 + poissons[i]));
                                        lame_lambda = youngs[i] * poissons[i]  / ((1.0 + poissons[i]) * (1.0 - 2.0 * poissons[i]));
                                        bulk = lame_lambda + 2.0 * shear / 3.0;
                                        E_reduced = Lookup[CONTACT_FORCE_SCALING][itype]*1.0/(((1-pow(poissons[i],2))/youngs[i])+((1-pow(poissons[j],2))/youngs[j]));
                                fpair = 1.066666667e0 * E_reduced * delta * sqrt(delta * r_geom); //  units: N
                                //damping
                                double damping_coeff, velocity, dt_crit;
                                velocity = dv.norm();
                                if(velocity==0.0){
                                  damping_coeff = 0.0;
                                }else{
                                  if(impact_velocity[i][jj]<0.0){
                                    impact_velocity[i][jj] = velocity;  
                                  }
                                  //damping_coeff = Lookup[CONTACT_DAMPING_SCALING][itype]*(8*E_reduced*(1-Lookup[CONTACT_RESTITUTION][itype]))/(5*Lookup[CONTACT_RESTITUTION][itype]*impact_velocity[i][jj]);
                                  damping_coeff = Lookup[CONTACT_DAMPING_SCALING][itype]*0.587*(1-Lookup[CONTACT_RESTITUTION][itype]);
                                }
                                if(dx.dot(dv)<0){
                                fpair += damping_coeff * 1.066666667e0 * sqrt(delta * r_geom)* velocity; //moving towards each other
                                }else{
                                fpair -= damping_coeff * 1.066666667e0 * sqrt(delta * r_geom)* velocity; //moving away from each other
                                }
                                if(fpair!=fpair){
                                  printf("damping_coeff: %f, delta: %f, r_geom: %f , velocity: %f \n",damping_coeff,delta,r_geom,velocity);
                                  printf("DAMPING_SCALING: %f , RESTITUTION: %f , E_reduced: %f  \n",Lookup[CONTACT_DAMPING_SCALING][itype],Lookup[CONTACT_RESTITUTION][itype],E_reduced);
                                  printf("impact_velocity: %f \n",impact_velocity[i][jj]);
                                  printf("------------------------------------------- \n");
                                }
                                fpair = MAX(fpair,0.0);
                                evdwl = fpair * 0.4e0 * delta; // GCG 25 April: this expression conserves total energy
                                if (r > 2.0e-16) {
                                        fpair /= r; // divide by r and multiply with non-normalized distance vector
                                } else {
                                        fpair = 0.0;
                                }
                                if (evflag) {
                                        ev_tally(i, j, nlocal, 0, evdwl, 0.0, fpair, dx[0], dx[1], dx[2]);
                                }
                                f[i][0] += dx[0] * fpair;
                                f[i][1] += dx[1] * fpair;
                                f[i][2] += dx[2] * fpair;

                                dt_crit = 3.14 * sqrt(0.5 * (rmass[i] + rmass[j]) / (abs(fpair) / delta));
                                if(dt_crit==dt_crit){
                                dtCFL = MIN(dtCFL, dt_crit);
                                }
                              }

                                

        }
      }
        continue;
      }
      j = atom->map(partner[i][jj]);
      if (j < 0) { //                 // check if lost a partner without first breaking bond
        partner[i][jj] = 0;
        continue;
      }

      if (mol[j] < 0) {
        continue; // Particle j is not a valid SPH particle (anymore). Skip all interactions with this particle.
      }

      if (mol[i] != mol[j]) {
        continue;
      }

      if (type[j] != itype)
        error->all(FLERR, "particle pair is not of same type!");

      for (idim = 0; idim < 3; idim++) {
        x0j(idim) = x0[j][idim];
        xj(idim) = x[j][idim];
        vj(idim) = v[j][idim];
      }

      if (periodic)
        domain->minimum_image(dx0(0), dx0(1), dx0(2));

      // check that distance between i and j (in the reference config) is less than cutoff
      dx0 = x0j - x0i;
      r0Sq = dx0.squaredNorm();
      h = radius[i] + radius[j];
      hMin = MIN(hMin, h);
      r0 = sqrt(r0Sq);
      volj = vfrac[j];

      // distance vectors in current and reference configuration, velocity difference
      dx = xj - xi;
      dv = vj - vi;
      r = dx.norm(); // current distance

      // scale the interaction according to the damage variable
      scale = 1.0 - degradation_ij[i][jj];
      wf = wf_list[i][jj] * scale;
      wfd = wfd_list[i][jj] * scale;

      g = (wfd / r0) * dx0; // uncorrected kernel gradient


      /*
       * force contribution -- note that the kernel gradient correction has been absorbed into PK1
       */

      f_stress = -voli * volj * (PK1[i] + PK1[j]) * g;
      

      /*
       * artificial viscosity
       */
      delVdotDelR = dx.dot(dv) / (r + 0.1 * h); // project relative velocity onto unit particle distance vector [m/s]
      LimitDoubleMagnitude(delVdotDelR, 0.01 * Lookup[SIGNAL_VELOCITY][itype]);
      mu_ij = h * delVdotDelR / (r + 0.1 * h); // units: [m * m/s / m = m/s]
      visc_magnitude = (-Lookup[VISCOSITY_Q1][itype] * Lookup[SIGNAL_VELOCITY][itype] * mu_ij
                        + Lookup[VISCOSITY_Q2][itype] * mu_ij * mu_ij) / Lookup[REFERENCE_DENSITY][itype]; // units: m^5/(s^2 kg))
      f_visc = rmass[i] * rmass[j] * visc_magnitude * wfd * dx / (r + 1.0e-2 * h); // units: kg^2 * m^5/(s^2 kg) * m^-4 = kg m / s^2 = N

      /*
       * hourglass deviation of particles i and j
       */

      gamma = 0.5 * (Fincr[i] + Fincr[j]) * dx0 - dx;
      hg_err = gamma.norm() / r0;
      hourglass_error[i] += volj * wf * hg_err;

      /* SPH-like hourglass formulation */

      if (MAX(plastic_strain[i], plastic_strain[j]) > 1.0e-3) {
        /*
         * viscous hourglass formulation for particles with plastic deformation
         */
        delta = gamma.dot(dx);
        if (delVdotDelR * delta < 0.0) {
          hg_err = MAX(hg_err, 0.05); // limit hg_err to avoid numerical instabilities
          hg_mag = -hg_err * Lookup[HOURGLASS_CONTROL_AMPLITUDE][itype] * Lookup[SIGNAL_VELOCITY][itype] * mu_ij
            / Lookup[REFERENCE_DENSITY][itype]; // this has units of pressure
        } else {
          hg_mag = 0.0;
        }
        f_hg = rmass[i] * rmass[j] * hg_mag * wfd * dx / (r + 1.0e-2 * h);

      } else {
        /*
         * stiffness hourglass formulation for particle in the elastic regime
         */

        gamma_dot_dx = gamma.dot(dx); // project hourglass error vector onto pair distance vector
        LimitDoubleMagnitude(gamma_dot_dx, 0.1 * r); // limit projected vector to avoid numerical instabilities
        delta = 0.5 * gamma_dot_dx / (r + 0.1 * h); // delta has dimensions of [m]
        hg_mag = Lookup[HOURGLASS_CONTROL_AMPLITUDE][itype] * delta / (r0Sq + 0.01 * h * h); // hg_mag has dimensions [m^(-1)]
        hg_mag *= -voli * volj * wf * youngs[i]; // hg_mag has dimensions [J*m^(-1)] = [N]
        f_hg = (hg_mag / (r + 0.01 * h)) * dx;
      }

      // scale hourglass force with damage
      f_hg *= (1.0 - damage[i]) * (1.0 - damage[j]);

      // sum stress, viscous, and hourglass forces
      sumForces = f_stress + f_visc + f_hg; // + f_spring;

      //potential energy rate -- project velocity onto force vector
      deltaE = 0.5 * sumForces.dot(dv);

      if((deltaE>0|deltaE<0)&false){
      printf("voli: %f, volj: %f, PK1: %f %f %f, PK1: %f %f %f, g: %f %f %f \n",voli,volj,PK1[i](0,0),PK1[i](1,1),PK1[i](2,2),PK1[j](0,0),PK1[j](1,1),PK1[j](2,2),g(0),g(1),g(2));
      printf("f_stress: %f %f %f, f_visc: %f %f %f, f_hg: %f %f %f \n",f_stress(0),f_stress(1),f_stress(2),f_visc(0),f_visc(1),f_visc(2),f_hg(0),f_hg(1),f_hg(2));
      printf("sumforces: %f %f %f, dv: %f %f %f \n",sumForces(0),sumForces(1),sumForces(2),dv(0),dv(1),dv(2));
      printf("deltaE: %f \n",deltaE);
      }
      // acount for heat transfer    
      if(Fincr[j]==eye) {
        detF[j] = 1;
      } else {
        detF[j] = Fincr[j].determinant();
      }  
      rhosph[i] = rmass[i] / (detF[i] * vfrac[i]);
      rhosph[j] = rmass[j] / (detF[j] * vfrac[j]);

      if(contact_thermal[itype]){
      deltaHE =  2.0*rmass[i] * rmass[j] / (rmass[i]+rmass[j]);
      deltaHE *=  (rhosph[i] + rhosph[j]) / (rhosph[i] * rhosph[j]);
      deltaHE *=  2/((1/ksph[i])+(1/ksph[j]));  //harmonic mean
      deltaHE *=  (cpsph[i] + cpsph[j])/(2*cpsph[i] * cpsph[j]);  //harmonic mean
      deltaHE *= (1/(rhosph[i])) * (esph[i] - esph[j]) * wfd/r;      


      // calculate heatflux
      deltaQ = rmass[i] * rmass[j] / (rmass[i]+rmass[j]);
      deltaQ *=  (rhosph[i] + rhosph[j]) / (rhosph[i] * rhosph[j]);
      deltaQ *= (-2/((1/ksph[i])+(1/ksph[j]))) * (tsph[i] - tsph[j]) * wfd/r; 

      if(false){
    printf("j: %d \n",i);  //added for testing
    printf("cpsph: %f \n",cpsph[j]);  //added for testing
    printf("ksph: %f \n",ksph[j]);  //added for testing
    printf("fracsph: %f \n",fracsph[j]);  //added for testing
    printf("qsph: %f \n",qsph[j][0]);  //added for testing
    printf("tsph: %f \n",tsph[j]);  //added for testing
    printf("deltaQ: %f \n",deltaQ);  //added for testing
    printf("deltaHE: %f \n",deltaHE);  //added for testing
      }

      dpesph[i] += deltaE; //potential energy stored seperate to thermal energy
      desph[i] += deltaHE;
      qsph[i][0] += deltaQ*dx[0]* K[i](0,0);
      qsph[i][1] += deltaQ*dx[1]* K[i](1,1);
      qsph[i][2] += deltaQ*dx[2]* K[i](2,2);
    }

      // apply forces to pair of particles
      f[i][0] += sumForces(0);
      f[i][1] += sumForces(1);
      f[i][2] += sumForces(2);


      // tally atomistic stress tensor
      if (evflag) {
        ev_tally_xyz(i, j, nlocal, 0, 0.0, 0.0, sumForces(0), sumForces(1), sumForces(2), dx(0), dx(1), dx(2));
      }

      shepardWeight += wf * volj;

      // check if a particle has moved too much w.r.t another particle
      if (r > r0) {
        if (update_method == UPDATE_CONSTANT_THRESHOLD) {
          if (r - r0 > update_threshold) {
            updateFlag = 1;
          }
        } else if (update_method == UPDATE_PAIRWISE_RATIO) {
          if ((r - r0) / h > update_threshold) {
            updateFlag = 1;
          }
        }
      }

      if (failureModel[itype].failure_max_pairwise_strain) {
        double r0alpha = r0 + r0*linear_expansion[i]*(tsph[i]-Lookup[EOS_THERMAL_T0][itype]);
        
        //strain1d = MAX(((r - r0alpha) / r0alpha),0.0); // - linear_expansion[i]*(tsph[i]-Lookup[EOS_THERMAL_T0][itype]);
        strain1d = MAX(((r - r0) / r0),0.0);
        if(false){
          printf("r: %f \n",r);  //added for testing
          printf("r0: %f \n",r0);  //added for testing
          printf("r0alpha: %f \n",r0alpha);  //added for testing
          printf("strain1d: %f \n",strain1d);  //added for testing
          printf("-------------------------------------------: \n");  //added for testing
        }
        //strain1d_max = Lookup[FAILURE_MAX_PAIRWISE_STRAIN_THRESHOLD][itype];
        strain1d_max = fracsph[i]; //MH changed to allow weibul per particle variation
        softening_strain = 2.0 * strain1d_max;

        if (strain1d > strain1d_max) {
          degradation_ij[i][jj] = (strain1d - strain1d_max) / softening_strain;
        } else {
          degradation_ij[i][jj] = 0.0;
        }

        if (degradation_ij[i][jj] >= 1.0) { // delete interaction if fully damaged
          //damage[i] += 0.0001;
          partner[i][jj] = 0;
        }
      }

      if (failureModel[itype].failure_energy_release_rate) {

        // integration approach
        energy_per_bond[i][jj] += update->dt * f_stress.dot(dv) / (voli * volj);
        double Vic = (2.0 / 3.0) * h * h * h; // interaction volume for 2d plane strain
        double critical_energy_per_bond = Lookup[CRITICAL_ENERGY_RELEASE_RATE][itype] / (2.0 * Vic);

        if (energy_per_bond[i][jj] > critical_energy_per_bond) {
          //degradation_ij[i][jj] = 1.0;
          partner[i][jj] = 0;
        }
      }

      if (failureModel[itype].integration_point_wise) {

        strain1d = (r - r0) / r0;

        if (strain1d > 0.0) {

          if ((damage[i] == 1.0) && (damage[j] == 1.0)) {
            // check if damage_onset is already defined
            if (energy_per_bond[i][jj] == 0.0) { // pair damage not defined yet
              energy_per_bond[i][jj] = strain1d;
            } else { // damage initiation strain already defined
              strain1d_max = energy_per_bond[i][jj];
              softening_strain = 2.0 * strain1d_max;

              if (strain1d > strain1d_max) {
                degradation_ij[i][jj] = (strain1d - strain1d_max) / softening_strain;
              } else {
                degradation_ij[i][jj] = 0.0;
              }
            }
          }

          if (degradation_ij[i][jj] >= 1.0) { // delete interaction if fully damaged
            //damage[i] += 0.0001;
            partner[i][jj] = 0;
          }

        } else {
          degradation_ij[i][jj] = 0.0;
        } // end failureModel[itype].integration_point_wise

      }




    } // end loop over jj neighbors of i


    //acount for heat generation
    if(contact_thermal[itype]){
      desph[i] += (hgsph[i]*rmass[i]);
    }

    //heat change due to stress thermo elastic:
      //desph[i] += -Lookup[EOS_THERMAL_ALPHA][itype]*Lookup[EOS_THERMAL_T0][itype]*(PK1[i](0, 0)  +  PK1[i](1, 1) + PK1[i](2, 2));
    

    // avoid division by zero and overflow
    if ((shepardWeight != 0.0) && (fabs(hourglass_error[i]) < 1.0e300)) {
      hourglass_error[i] /= shepardWeight;
    }

  } // end loop over i

  if (vflag_fdotr)
    virial_fdotr_compute();
}

/* ----------------------------------------------------------------------
 assemble unrotated stress tensor using deviatoric and pressure components.
 Convert to corotational Cauchy stress, then to PK1 stress and apply
 shape matrix correction
 ------------------------------------------------------------------------- */
void PairTlsph::AssembleStress() {
  tagint *mol = atom->molecule;
  double *eff_plastic_strain = atom->eff_plastic_strain;
  double *eff_plastic_strain_rate = atom->eff_plastic_strain_rate;
  double **tlsph_stress = atom->smd_stress;
  int *type = atom->type;
  double *radius = atom->radius;
  double *damage = atom->damage;
  double *rmass = atom->rmass;
  double *vfrac = atom->vfrac;
  double *esph = atom->esph;
  double *pesph = atom->pesph;
  double *tsph = atom->tsph;
  double *youngs = atom->youngs;
  double *poissons = atom->poissons;
  double pInitial, d_iso, pFinal, p_rate, plastic_strain_increment;
  int i, itype;
  int nlocal = atom->nlocal;
  double dt = update->dt;
  double M_eff, p_wave_speed, mass_specific_energy, vol_specific_energy, rho;
  Matrix3d sigma_rate, eye, sigmaInitial, sigmaFinal, T, T_damaged, Jaumann_rate, sigma_rate_check;
  Matrix3d d_dev, sigmaInitial_dev, sigmaFinal_dev, sigma_dev_rate, strain;
  Vector3d x0i, xi, xp;

  eye.setIdentity();
  dtCFL = 1.0e22;
  pFinal = 0.0;

  for (i = 0; i < nlocal; i++) {
    particle_dt[i] = 0.0;

    itype = type[i];
    if (setflag[itype][itype] == 1) {
      if (mol[i] > 0) { // only do the following if particle has not failed -- mol < 0 means particle has failed

        /*
         * initial stress state: given by the unrotateted Cauchy stress.
         * Assemble Eigen 3d matrix from stored stress state
         */
        sigmaInitial(0, 0) = tlsph_stress[i][0];
        sigmaInitial(0, 1) = tlsph_stress[i][1];
        sigmaInitial(0, 2) = tlsph_stress[i][2];
        sigmaInitial(1, 1) = tlsph_stress[i][3];
        sigmaInitial(1, 2) = tlsph_stress[i][4];
        sigmaInitial(2, 2) = tlsph_stress[i][5];
        sigmaInitial(1, 0) = sigmaInitial(0, 1);
        sigmaInitial(2, 0) = sigmaInitial(0, 2);
        sigmaInitial(2, 1) = sigmaInitial(1, 2);

        //cout << "this is sigma initial" << endl << sigmaInitial << endl;

        pInitial = sigmaInitial.trace() / 3.0; // isotropic part of initial stress
        sigmaInitial_dev = Deviator(sigmaInitial);
        d_iso = D[i].trace(); // volumetric part of stretch rate
        d_dev = Deviator(D[i]); // deviatoric part of stretch rate
        strain = 0.5 * (Fincr[i].transpose() * Fincr[i] - eye);
        mass_specific_energy = (esph[i] + pesph[i]) / rmass[i]; // energy per unit mass
        rho = rmass[i] / (detF[i] * vfrac[i]);
        vol_specific_energy = mass_specific_energy * rho; // energy per current volume

        /*
         * pressure: compute pressure rate p_rate and final pressure pFinal
         */

        ComputePressure(i, rho, mass_specific_energy, tsph[i], vol_specific_energy, pInitial, d_iso, pFinal, p_rate);

        /*
         * material strength
         */

        //cout << "this is the strain deviator rate" << endl << d_dev << endl;
        ComputeStressDeviator(i, sigmaInitial_dev, d_dev, sigmaFinal_dev, sigma_dev_rate, plastic_strain_increment);
        //cout << "this is the stress deviator rate" << endl << sigma_dev_rate << endl;

        // keep a rolling average of the plastic strain rate over the last 100 or so timesteps
        eff_plastic_strain[i] += plastic_strain_increment;

        // compute a characteristic time over which to average the plastic strain
        double tav = 1000 * radius[i] / (Lookup[SIGNAL_VELOCITY][itype]);
        eff_plastic_strain_rate[i] -= eff_plastic_strain_rate[i] * dt / tav;
        eff_plastic_strain_rate[i] += plastic_strain_increment / tav;
        eff_plastic_strain_rate[i] = MAX(0.0, eff_plastic_strain_rate[i]);

        /*
         *  assemble total stress from pressure and deviatoric stress
         */
        sigmaFinal = pFinal * eye + sigmaFinal_dev; // this is the stress that is kept

        if (JAUMANN) {
          /*
           * sigma is already the co-rotated Cauchy stress.
           * The stress rate, however, needs to be made objective.
           */

          if (dt > 1.0e-16) {
            sigma_rate = (1.0 / dt) * (sigmaFinal - sigmaInitial);
          } else {
            sigma_rate.setZero();
          }

          Jaumann_rate = sigma_rate + W[i] * sigmaInitial + sigmaInitial * W[i].transpose();
          sigmaFinal = sigmaInitial + dt * Jaumann_rate;
          T = sigmaFinal;
        } else {
          /*
           * sigma is the unrotated stress.
           * need to do forward rotation of the unrotated stress sigma to the current configuration
           */
          T = R[i] * sigmaFinal * R[i].transpose();
        }

        /*
         * store unrotated stress in atom vector
         * symmetry is exploited
         */
        tlsph_stress[i][0] = sigmaFinal(0, 0);
        tlsph_stress[i][1] = sigmaFinal(0, 1);
        tlsph_stress[i][2] = sigmaFinal(0, 2);
        tlsph_stress[i][3] = sigmaFinal(1, 1);
        tlsph_stress[i][4] = sigmaFinal(1, 2);
        tlsph_stress[i][5] = sigmaFinal(2, 2);

        /*
         *  Damage due to failure criteria.
         */

        if (failureModel[itype].integration_point_wise) {
          ComputeDamage(i, strain, T, T_damaged);
          //T = T_damaged; Do not do this, it is undefined as of now
        }

        // store rotated, "true" Cauchy stress
        CauchyStress[i] = T;

        /*
         * We have the corotational Cauchy stress.
         * Convert to PK1. Note that reference configuration used for computing the forces is linked via
         * the incremental deformation gradient, not the full deformation gradient.
         */
        PK1[i] = detF[i] * T * FincrInv[i].transpose();

        /*
         * pre-multiply stress tensor with shape matrix to save computation in force loop
         */
        PK1[i] = PK1[i] * K[i];

        /*
         * compute stable time step according to Pronto 2d
         */

        Matrix3d deltaSigma;
        deltaSigma = sigmaFinal - sigmaInitial;
        p_rate = deltaSigma.trace() / (3.0 * dt + 1.0e-16);
        sigma_dev_rate = Deviator(deltaSigma) / (dt + 1.0e-16);

        double K_eff, mu_eff, shear, bulk, m_modulus, lame_lambda;
         shear = youngs[i] / (2.0 * (1.0 + poissons[i]));
         lame_lambda = youngs[i] * poissons[i]  / ((1.0 + poissons[i]) * (1.0 - 2.0 * poissons[i]));
         m_modulus = lame_lambda + 2.0 * shear;
         bulk = lame_lambda + 2.0 * shear / 3.0;

        effective_longitudinal_modulus(itype, shear, bulk, m_modulus, dt, d_iso, p_rate, d_dev, sigma_dev_rate, damage[i], K_eff, mu_eff, M_eff);
        p_wave_speed = sqrt(M_eff / rho);

        if (mol[i] < 0) {
          error->one(FLERR, "this should not happen");
        }

        particle_dt[i] = 2.0 * radius[i] / p_wave_speed;
        dtCFL = MIN(dtCFL, particle_dt[i]);

      } else { // end if mol > 0
        PK1[i].setZero();
        K[i].setIdentity();
        CauchyStress[i].setZero();
        sigma_rate.setZero();
        tlsph_stress[i][0] = 0.0;
        tlsph_stress[i][1] = 0.0;
        tlsph_stress[i][2] = 0.0;
        tlsph_stress[i][3] = 0.0;
        tlsph_stress[i][4] = 0.0;
        tlsph_stress[i][5] = 0.0;
      } // end  if mol > 0
    } // end setflag
  } // end for
}

/* ----------------------------------------------------------------------
 allocate all arrays
 ------------------------------------------------------------------------- */

void PairTlsph::allocate() {
  allocated = 1;
  int n = atom->ntypes;

  memory->create(setflag, n + 1, n + 1, "pair:setflag");
  for (int i = 1; i <= n; i++)
    for (int j = i; j <= n; j++)
      setflag[i][j] = 0;

  memory->create(strengthModel, n + 1, "pair:strengthmodel");
  memory->create(eos, n + 1, "pair:eosmodel");
  memory->create(contact_mechanical, n + 1, "pair:contact_mechanical");
  memory->create(contact_thermal, n + 1, "pair:contact_thermal");
  failureModel = new failure_types[n + 1];
  memory->create(Lookup, MAX_KEY_VALUE, n + 1, "pair:LookupTable");

  memory->create(cutsq, n + 1, n + 1, "pair:cutsq"); // always needs to be allocated, even with granular neighborlist

  onerad_dynamic = new double[n + 1];
  onerad_frozen = new double[n + 1];
  maxrad_dynamic = new double[n + 1];
  maxrad_frozen = new double[n + 1];

}

/* ----------------------------------------------------------------------
 global settings
 ------------------------------------------------------------------------- */

void PairTlsph::settings(int narg, char **arg) {

  if (comm->me == 0)
    utils::logmesg(lmp,"\n>>========>>========>>========>>========>>========>>========>>========>>========>>========>>========>>========>>========\n"
                   "TLSPH settings\n");

  /*
   * default value for update_threshold for updates of reference configuration:
   * The maximum relative displacement which is tracked by the construction of LAMMPS' neighborlists
   * is the following.
   */

  cut_comm = MAX(neighbor->cutneighmax, comm->cutghostuser); // cutoff radius within which ghost atoms are communicated.
  update_threshold = cut_comm;
  update_method = UPDATE_NONE;

  int iarg = 0;

  while (true) {

    if (iarg >= narg) {
      break;
    }

    if (strcmp(arg[iarg], "*UPDATE_CONSTANT") == 0) {
      iarg++;
      if (iarg == narg) {
        error->all(FLERR, "expected number following *UPDATE_CONSTANT keyword");
      }

      update_method = UPDATE_CONSTANT_THRESHOLD;
      update_threshold = utils::numeric(FLERR, arg[iarg],false,lmp);

    } else if (strcmp(arg[iarg], "*UPDATE_PAIRWISE") == 0) {
      iarg++;
      if (iarg == narg) {
        error->all(FLERR, "expected number following *UPDATE_PAIRWISE keyword");
      }

      update_method = UPDATE_PAIRWISE_RATIO;
      update_threshold = utils::numeric(FLERR, arg[iarg],false,lmp);

    } else {
      error->all(FLERR, "Illegal keyword for spooff/integrate_tlsph: {}\n", arg[iarg]);
    }
    iarg++;
  }

  if ((update_threshold > cut_comm) && (update_method == UPDATE_CONSTANT_THRESHOLD)) {
    if (comm->me == 0) {
      utils::logmesg(lmp, "\n                ***** WARNING ***\n");
      utils::logmesg(lmp, "requested reference configuration update threshold is {} length units\n", update_threshold);
      utils::logmesg(lmp, "This value exceeds the maximum value {} beyond which TLSPH displacements can be tracked at current settings.\n",cut_comm);
      utils::logmesg(lmp, "Expect loss of neighbors!\n");
    }
  }

  if (comm->me == 0) {
    if (update_method == UPDATE_CONSTANT_THRESHOLD) {
      utils::logmesg(lmp, "... will update reference configuration if magnitude of relative displacement exceeds {} length units\n",
             update_threshold);
    } else if (update_method == UPDATE_PAIRWISE_RATIO) {
      utils::logmesg(lmp, "... will update reference configuration if ratio pairwise distance / smoothing length  exceeds {}\n",
             update_threshold);
    } else if (update_method == UPDATE_NONE) {
      utils::logmesg(lmp, "... will never update reference configuration\n");
    }
    utils::logmesg(lmp,">>========>>========>>========>>========>>========>>========>>========>>========>>========>>========>>========>>========\n");
  }
}

/* ----------------------------------------------------------------------
 set coeffs for one or more type pairs
 ------------------------------------------------------------------------- */

void PairTlsph::coeff(int narg, char **arg) {
  int ioffset, iarg, iNextKwd, itype;
  std::string s, t;

  if (narg < 3)
    error->all(FLERR, "number of arguments for pair tlsph is too small!");

  if (!allocated)
    allocate();

  /*
   * check that TLSPH parameters are given only in i,i form
   */
  if (utils::inumeric(FLERR, arg[0], false, lmp) != utils::inumeric(FLERR, arg[1], false, lmp))
    error->all(FLERR, "TLSPH coefficients can only be specified between particles of same type!");

  itype = utils::inumeric(FLERR, arg[0],false,lmp);

// set all eos, strength and failure models to inactive by default
  eos[itype] = EOS_NONE;
  strengthModel[itype] = STRENGTH_NONE;
  contact_thermal[itype] = true;
  contact_mechanical[itype] = false;

  if (comm->me == 0)
    utils::logmesg(lmp,"\n>>========>>========>>========>>========>>========>>========>>========>>========>>========>>========>>========>>========\n"
                   "SPOOFF / TLSPH PROPERTIES OF PARTICLE TYPE {}:\n", itype);

  /*
   * read parameters which are common -- regardless of material / eos model
   */

  ioffset = 2;
  if (strcmp(arg[ioffset], "*COMMON") != 0)
    error->all(FLERR, "common keyword missing!");

  t = std::string("*");
  iNextKwd = -1;
  for (iarg = ioffset + 1; iarg < narg; iarg++) {
    s = std::string(arg[iarg]);
    if (s.compare(0, t.length(), t) == 0) {
      iNextKwd = iarg;
      break;
    }
  }

  if (iNextKwd < 0) error->all(FLERR, "no *KEYWORD terminates *COMMON");
  if (iNextKwd - ioffset != 9 + 1)
    error->all(FLERR, "expected 8 arguments following *COMMON but got {}\n", iNextKwd - ioffset - 1);

  Lookup[REFERENCE_DENSITY][itype] = utils::numeric(FLERR, arg[ioffset + 1],false,lmp);
  Lookup[YOUNGS_MODULUS][itype] = utils::numeric(FLERR, arg[ioffset + 2],false,lmp);
  Lookup[POISSON_RATIO][itype] = utils::numeric(FLERR, arg[ioffset + 3],false,lmp);
  Lookup[VISCOSITY_Q1][itype] = utils::numeric(FLERR, arg[ioffset + 4],false,lmp);
  Lookup[VISCOSITY_Q2][itype] = utils::numeric(FLERR, arg[ioffset + 5],false,lmp);
  Lookup[HOURGLASS_CONTROL_AMPLITUDE][itype] = utils::numeric(FLERR, arg[ioffset + 6],false,lmp);
  Lookup[HEAT_CAPACITY][itype] = utils::numeric(FLERR, arg[ioffset + 7],false,lmp);
  Lookup[THERMAL_CONDUCTIVITY][itype] = utils::numeric(FLERR, arg[ioffset + 8],false,lmp);
  Lookup[HEAT_PRODUCED][itype] = utils::numeric(FLERR, arg[ioffset + 9],false,lmp);

  Lookup[LAME_LAMBDA][itype] = Lookup[YOUNGS_MODULUS][itype] * Lookup[POISSON_RATIO][itype]
    / ((1.0 + Lookup[POISSON_RATIO][itype]) * (1.0 - 2.0 * Lookup[POISSON_RATIO][itype]));
  Lookup[SHEAR_MODULUS][itype] = Lookup[YOUNGS_MODULUS][itype] / (2.0 * (1.0 + Lookup[POISSON_RATIO][itype]));
  Lookup[M_MODULUS][itype] = Lookup[LAME_LAMBDA][itype] + 2.0 * Lookup[SHEAR_MODULUS][itype];
  Lookup[SIGNAL_VELOCITY][itype] = sqrt(
    (Lookup[LAME_LAMBDA][itype] + 2.0 * Lookup[SHEAR_MODULUS][itype]) / Lookup[REFERENCE_DENSITY][itype]);
  Lookup[BULK_MODULUS][itype] = Lookup[LAME_LAMBDA][itype] + 2.0 * Lookup[SHEAR_MODULUS][itype] / 3.0;

  if (comm->me == 0) {
    utils::logmesg(lmp, "\nmaterial unspecific properties for SPOOFF/TLSPH definition of particle type {}:\n", itype);
    utils::logmesg(lmp, "{:60} : {}\n", "reference density", Lookup[REFERENCE_DENSITY][itype]);
    utils::logmesg(lmp, "{:60} : {}\n", "Young's modulus", Lookup[YOUNGS_MODULUS][itype]);
    utils::logmesg(lmp, "{:60} : {}\n", "Poisson ratio", Lookup[POISSON_RATIO][itype]);
    utils::logmesg(lmp, "{:60} : {}\n", "linear viscosity coefficient", Lookup[VISCOSITY_Q1][itype]);
    utils::logmesg(lmp, "{:60} : {}\n", "quadratic viscosity coefficient", Lookup[VISCOSITY_Q2][itype]);
    utils::logmesg(lmp, "{:60} : {}\n", "hourglass control coefficient", Lookup[HOURGLASS_CONTROL_AMPLITUDE][itype]);
    utils::logmesg(lmp, "{:60} : {}\n", "heat capacity [energy / (mass * temperature)]", Lookup[HEAT_CAPACITY][itype]);
    utils::logmesg(lmp, "{:60} : {}\n", "thermal conductivity ", Lookup[THERMAL_CONDUCTIVITY][itype]);
    utils::logmesg(lmp, "{:60} : {}\n", "heat produced per unit mass ", Lookup[HEAT_PRODUCED][itype]);
    utils::logmesg(lmp, "{:60} : {}\n", "Lame constant lambda", Lookup[LAME_LAMBDA][itype]);
    utils::logmesg(lmp, "{:60} : {}\n", "shear modulus", Lookup[SHEAR_MODULUS][itype]);
    utils::logmesg(lmp, "{:60} : {}\n", "bulk modulus", Lookup[BULK_MODULUS][itype]);
    utils::logmesg(lmp, "{:60} : {}\n", "signal velocity", Lookup[SIGNAL_VELOCITY][itype]);
  }

  /*
   * read following material cards
   */

  eos[itype] = EOS_NONE;
  strengthModel[itype] = STRENGTH_NONE;
  contact_thermal[itype] = true;
  contact_mechanical[itype] = false;

  while (true) {
    if (strcmp(arg[iNextKwd], "*END") == 0) {
      if (comm->me == 0)
        utils::logmesg(lmp,"found *END keyword"
                       "\n>>========>>========>>========>>========>>========>>========>>========>>========>>========>>========>>========>>========\n\n");
      break;
    }

    /*
     * Linear Elasticity model based on deformation gradient
     */
    ioffset = iNextKwd;
    if (strcmp(arg[ioffset], "*LINEAR_DEFGRAD") == 0) {
      strengthModel[itype] = LINEAR_DEFGRAD;

      if (comm->me == 0) utils::logmesg(lmp, "reading *LINEAR_DEFGRAD\n");

      t = std::string("*");
      iNextKwd = -1;
      for (iarg = ioffset + 1; iarg < narg; iarg++) {
        s = std::string(arg[iarg]);
        if (s.compare(0, t.length(), t) == 0) {
          iNextKwd = iarg;
          break;
        }
      }

      if (iNextKwd < 0)
        error->all(FLERR, "no *KEYWORD terminates *LINEAR_DEFGRAD");

      if (iNextKwd - ioffset != 1)
        error->all(FLERR, "expected 0 arguments following *LINEAR_DEFGRAD but got {}\n", iNextKwd - ioffset - 1);

      if (comm->me == 0) utils::logmesg(lmp, "\nLinear Elasticity model based on deformation gradient\n");

    } else if (strcmp(arg[ioffset], "*STRENGTH_LINEAR") == 0) {

      /*
       * Linear Elasticity strength only model based on strain rate
       */

      strengthModel[itype] = STRENGTH_LINEAR;
      if (comm->me == 0) utils::logmesg(lmp,"reading *STRENGTH_LINEAR\n");

      t = std::string("*");
      iNextKwd = -1;
      for (iarg = ioffset + 1; iarg < narg; iarg++) {
        s = std::string(arg[iarg]);
        if (s.compare(0, t.length(), t) == 0) {
          iNextKwd = iarg;
          break;
        }
      }

      if (iNextKwd < 0)
        error->all(FLERR, "no *KEYWORD terminates *STRENGTH_LINEAR");

      if (iNextKwd - ioffset != 1)
        error->all(FLERR, "expected 0 arguments following *STRENGTH_LINEAR but got {}\n", iNextKwd - ioffset - 1);

      if (comm->me == 0) utils::logmesg(lmp, "Linear Elasticity strength based on strain rate\n");

    } // end Linear Elasticity strength only model based on strain rate

    else if (strcmp(arg[ioffset], "*STRENGTH_LINEAR_PLASTIC") == 0) {

      /*
       * Linear Elastic / perfectly plastic strength only model based on strain rate
       */

      strengthModel[itype] = STRENGTH_LINEAR_PLASTIC;
      if (comm->me == 0) utils::logmesg(lmp,"reading *STRENGTH_LINEAR_PLASTIC\n");

      t = std::string("*");
      iNextKwd = -1;
      for (iarg = ioffset + 1; iarg < narg; iarg++) {
        s = std::string(arg[iarg]);
        if (s.compare(0, t.length(), t) == 0) {
          iNextKwd = iarg;
          break;
        }
      }

      if (iNextKwd < 0)
        error->all(FLERR, "no *KEYWORD terminates *STRENGTH_LINEAR_PLASTIC");

      if (iNextKwd - ioffset != 2 + 1)
        error->all(FLERR, "expected 2 arguments following *STRENGTH_LINEAR_PLASTIC but got {}\n", iNextKwd - ioffset - 1);

      Lookup[YIELD_STRESS][itype] = utils::numeric(FLERR, arg[ioffset + 1],false,lmp);
      strainfail = false;
      Lookup[HARDENING_PARAMETER][itype] = utils::numeric(FLERR, arg[ioffset + 2],false,lmp);

      if (comm->me == 0) {
        utils::logmesg(lmp, "Linear elastic / perfectly plastic strength based on strain rate");
        utils::logmesg(lmp, "{:60} : {}\n", "Young's modulus", Lookup[YOUNGS_MODULUS][itype]);
        utils::logmesg(lmp, "{:60} : {}\n", "Poisson ratio", Lookup[POISSON_RATIO][itype]);
        utils::logmesg(lmp, "{:60} : {}\n", "shear modulus", Lookup[SHEAR_MODULUS][itype]);
        utils::logmesg(lmp, "{:60} : {}\n", "constant yield stress", Lookup[YIELD_STRESS][itype]);
        utils::logmesg(lmp, "{:60} : {}\n", "constant hardening parameter", Lookup[HARDENING_PARAMETER][itype]);
      }
    } // end Linear Elastic / perfectly plastic strength only model based on strain rate

    else if (strcmp(arg[ioffset], "*JOHNSON_COOK") == 0) {

      /*
       * JOHNSON - COOK
       */

      strengthModel[itype] = STRENGTH_JOHNSON_COOK;
      if (comm->me == 0) utils::logmesg(lmp, "reading *JOHNSON_COOK\n");

      t = std::string("*");
      iNextKwd = -1;
      for (iarg = ioffset + 1; iarg < narg; iarg++) {
        s = std::string(arg[iarg]);
        if (s.compare(0, t.length(), t) == 0) {
          iNextKwd = iarg;
          break;
        }
      }

      if (iNextKwd < 0)
        error->all(FLERR, "no *KEYWORD terminates *JOHNSON_COOK");

      if (iNextKwd - ioffset != 8 + 1)
        error->all(FLERR, "expected 8 arguments following *JOHNSON_COOK but got {}\n", iNextKwd - ioffset - 1);

      Lookup[JC_A][itype] = utils::numeric(FLERR, arg[ioffset + 1],false,lmp);
      Lookup[JC_B][itype] = utils::numeric(FLERR, arg[ioffset + 2],false,lmp);
      Lookup[JC_a][itype] = utils::numeric(FLERR, arg[ioffset + 3],false,lmp);
      Lookup[JC_C][itype] = utils::numeric(FLERR, arg[ioffset + 4],false,lmp);
      Lookup[JC_epdot0][itype] = utils::numeric(FLERR, arg[ioffset + 5],false,lmp);
      Lookup[JC_T0][itype] = utils::numeric(FLERR, arg[ioffset + 6],false,lmp);
      Lookup[JC_Tmelt][itype] = utils::numeric(FLERR, arg[ioffset + 7],false,lmp);
      Lookup[JC_M][itype] = utils::numeric(FLERR, arg[ioffset + 8],false,lmp);

      if (comm->me == 0) {
        utils::logmesg(lmp, "Johnson Cook material strength model\n");
        utils::logmesg(lmp, "{:60} : {}\n", "A: initial yield stress", Lookup[JC_A][itype]);
        utils::logmesg(lmp, "{:60} : {}\n", "B : proportionality factor for plastic strain dependency", Lookup[JC_B][itype]);
        utils::logmesg(lmp, "{:60} : {}\n", "a : exponent for plastic strain dependency", Lookup[JC_a][itype]);
        utils::logmesg(lmp, "{:60} : {}\n", "C : proportionality factor for logarithmic plastic strain rate dependency",Lookup[JC_C][itype]);
        utils::logmesg(lmp, "{:60} : {}\n", "epdot0 : dimensionality factor for plastic strain rate dependency", Lookup[JC_epdot0][itype]);
        utils::logmesg(lmp, "{:60} : {}\n", "T0 : reference (room) temperature", Lookup[JC_T0][itype]);
        utils::logmesg(lmp, "{:60} : {}\n", "Tmelt : melting temperature", Lookup[JC_Tmelt][itype]);
        utils::logmesg(lmp, "{:60} : {}\n", "M : exponent for temperature dependency", Lookup[JC_M][itype]);
      }

    } else if (strcmp(arg[ioffset], "*EOS_NONE") == 0) {

      /*
       * no eos
       */

      eos[itype] = EOS_NONE;
      if (comm->me == 0) utils::logmesg(lmp, "reading *EOS_NONE\n");

      t = std::string("*");
      iNextKwd = -1;
      for (iarg = ioffset + 1; iarg < narg; iarg++) {
        s = std::string(arg[iarg]);
        if (s.compare(0, t.length(), t) == 0) {
          iNextKwd = iarg;
          break;
        }
      }

      if (iNextKwd < 0)
        error->all(FLERR, "no *KEYWORD terminates *EOS_NONE");

      if (iNextKwd - ioffset != 1)
        error->all(FLERR, "expected 0 arguments following *EOS_NONE but got {}\n", iNextKwd - ioffset - 1);

      if (comm->me == 0) utils::logmesg(lmp, "\nno EOS selected\n");

    } else if (strcmp(arg[ioffset], "*EOS_LINEAR") == 0) {

      /*
       * linear eos
       */

      eos[itype] = EOS_LINEAR;
      if (comm->me == 0) utils::logmesg(lmp, "reading *EOS_LINEAR\n");

      t = std::string("*");
      iNextKwd = -1;
      for (iarg = ioffset + 1; iarg < narg; iarg++) {
        s = std::string(arg[iarg]);
        if (s.compare(0, t.length(), t) == 0) {
          iNextKwd = iarg;
          break;
        }
      }

      if (iNextKwd < 0)
        error->all(FLERR, "no *KEYWORD terminates *EOS_LINEAR");

      if (iNextKwd - ioffset != 1)
        error->all(FLERR, "expected 0 arguments following *EOS_LINEAR but got {}\n", iNextKwd - ioffset - 1);

      if (comm->me == 0) {
        utils::logmesg(lmp, "\nlinear EOS based on strain rate\n");
        utils::logmesg(lmp, "{:60} : {}\n", "bulk modulus", Lookup[BULK_MODULUS][itype]);
      }
    } // end linear eos

    else if (strcmp(arg[ioffset], "*EOS_THERMAL") == 0) {

      /*
       * linear with thermal expansion eos
       */

      eos[itype] = EOS_THERMAL;
      if (comm->me == 0) utils::logmesg(lmp, "reading *EOS_THERMAL\n");

      t = std::string("*");
      iNextKwd = -1;
      for (iarg = ioffset + 1; iarg < narg; iarg++) {
        s = std::string(arg[iarg]);
        if (s.compare(0, t.length(), t) == 0) {
          iNextKwd = iarg;
          break;
        }
      }

      if (iNextKwd < 0)
        error->all(FLERR, "no *KEYWORD terminates *EOS_THERMAL");

      if (iNextKwd - ioffset != 2 + 1)
        error->all(FLERR, "expected 2 arguments following *EOS_THERMAL but got {}\n", iNextKwd - ioffset - 1);
        
      Lookup[EOS_THERMAL_T0][itype] = utils::numeric(FLERR, arg[ioffset + 1],false,lmp);
      Lookup[EOS_THERMAL_ALPHA][itype] = utils::numeric(FLERR, arg[ioffset + 2],false,lmp);

      if (comm->me == 0) {
        utils::logmesg(lmp, "\nlinear EOS based on strain rate with thermal expansion\n");
        utils::logmesg(lmp, "{:60} : {}\n", "bulk modulus", Lookup[BULK_MODULUS][itype]);
        utils::logmesg(lmp, "{:60} : {}\n", "heat capacity", Lookup[HEAT_CAPACITY][itype]);
        utils::logmesg(lmp, "{:60} : {}\n", "reference temperature", Lookup[EOS_THERMAL_T0][itype]);
        utils::logmesg(lmp, "{:60} : {}\n", "thermal expansion coefficient", Lookup[EOS_THERMAL_ALPHA][itype]);
      }
    } // end linear with thermal expansion eos

    else if (strcmp(arg[ioffset], "*EOS_SHOCK") == 0) {

      /*
       * shock eos
       */

      eos[itype] = EOS_SHOCK;
      if (comm->me == 0) utils::logmesg(lmp, "reading *EOS_SHOCK\n");

      t = std::string("*");
      iNextKwd = -1;
      for (iarg = ioffset + 1; iarg < narg; iarg++) {
        s = std::string(arg[iarg]);
        if (s.compare(0, t.length(), t) == 0) {
          iNextKwd = iarg;
          break;
        }
      }

      if (iNextKwd < 0)
        error->all(FLERR, "no *KEYWORD terminates *EOS_SHOCK");

      if (iNextKwd - ioffset != 3 + 1)
        error->all(FLERR, "expected 3 arguments (c0, S, Gamma) following *EOS_SHOCK but got {}\n", iNextKwd - ioffset - 1);

      Lookup[EOS_SHOCK_C0][itype] = utils::numeric(FLERR, arg[ioffset + 1],false,lmp);
      Lookup[EOS_SHOCK_S][itype] = utils::numeric(FLERR, arg[ioffset + 2],false,lmp);
      Lookup[EOS_SHOCK_GAMMA][itype] = utils::numeric(FLERR, arg[ioffset + 3],false,lmp);
      if (comm->me == 0) {
        utils::logmesg(lmp, "\nshock EOS based on strain rate\n");
        utils::logmesg(lmp, "{:60} : {}\n", "reference speed of sound", Lookup[EOS_SHOCK_C0][itype]);
        utils::logmesg(lmp, "{:60} : {}\n", "Hugoniot parameter S", Lookup[EOS_SHOCK_S][itype]);
        utils::logmesg(lmp, "{:60} : {}\n", "Grueneisen Gamma", Lookup[EOS_SHOCK_GAMMA][itype]);
      }
    } // end shock eos

    else if (strcmp(arg[ioffset], "*EOS_POLYNOMIAL") == 0) {
      /*
       * polynomial eos
       */

      eos[itype] = EOS_POLYNOMIAL;
      if (comm->me == 0) utils::logmesg(lmp, "reading *EOS_POLYNOMIAL\n");

      t = std::string("*");
      iNextKwd = -1;
      for (iarg = ioffset + 1; iarg < narg; iarg++) {
        s = std::string(arg[iarg]);
        if (s.compare(0, t.length(), t) == 0) {
          iNextKwd = iarg;
          break;
        }
      }

      if (iNextKwd < 0)
        error->all(FLERR, "no *KEYWORD terminates *EOS_POLYNOMIAL");

      if (iNextKwd - ioffset != 7 + 1)
        error->all(FLERR, "expected 7 arguments following *EOS_POLYNOMIAL but got {}\n", iNextKwd - ioffset - 1);

      Lookup[EOS_POLYNOMIAL_C0][itype] = utils::numeric(FLERR, arg[ioffset + 1],false,lmp);
      Lookup[EOS_POLYNOMIAL_C1][itype] = utils::numeric(FLERR, arg[ioffset + 2],false,lmp);
      Lookup[EOS_POLYNOMIAL_C2][itype] = utils::numeric(FLERR, arg[ioffset + 3],false,lmp);
      Lookup[EOS_POLYNOMIAL_C3][itype] = utils::numeric(FLERR, arg[ioffset + 4],false,lmp);
      Lookup[EOS_POLYNOMIAL_C4][itype] = utils::numeric(FLERR, arg[ioffset + 5],false,lmp);
      Lookup[EOS_POLYNOMIAL_C5][itype] = utils::numeric(FLERR, arg[ioffset + 6],false,lmp);
      Lookup[EOS_POLYNOMIAL_C6][itype] = utils::numeric(FLERR, arg[ioffset + 7],false,lmp);
      if (comm->me == 0) {
        utils::logmesg(lmp, "\npolynomial EOS based on strain rate\n");
        utils::logmesg(lmp, "{:60} : {}\n", "parameter c0", Lookup[EOS_POLYNOMIAL_C0][itype]);
        utils::logmesg(lmp, "{:60} : {}\n", "parameter c1", Lookup[EOS_POLYNOMIAL_C1][itype]);
        utils::logmesg(lmp, "{:60} : {}\n", "parameter c2", Lookup[EOS_POLYNOMIAL_C2][itype]);
        utils::logmesg(lmp, "{:60} : {}\n", "parameter c3", Lookup[EOS_POLYNOMIAL_C3][itype]);
        utils::logmesg(lmp, "{:60} : {}\n", "parameter c4", Lookup[EOS_POLYNOMIAL_C4][itype]);
        utils::logmesg(lmp, "{:60} : {}\n", "parameter c5", Lookup[EOS_POLYNOMIAL_C5][itype]);
        utils::logmesg(lmp, "{:60} : {}\n", "parameter c6", Lookup[EOS_POLYNOMIAL_C6][itype]);
      }
    } // end polynomial eos

    else if (strcmp(arg[ioffset], "*FAILURE_MAX_PLASTIC_STRAIN") == 0) {

      /*
       * maximum plastic strain failure criterion
       */

      if (comm->me == 0) utils::logmesg(lmp, "reading *FAILURE_MAX_PLASTIC_SRTRAIN\n");

      t = std::string("*");
      iNextKwd = -1;
      for (iarg = ioffset + 1; iarg < narg; iarg++) {
        s = std::string(arg[iarg]);
        if (s.compare(0, t.length(), t) == 0) {
          iNextKwd = iarg;
          break;
        }
      }

      if (iNextKwd < 0) error->all(FLERR, "no *KEYWORD terminates *FAILURE_MAX_PLASTIC_STRAIN");
      if (iNextKwd - ioffset != 1 + 1)
        error->all(FLERR, "expected 1 arguments following *FAILURE_MAX_PLASTIC_STRAIN but got {}\n", iNextKwd - ioffset - 1);

      failureModel[itype].failure_max_plastic_strain = true;
      failureModel[itype].integration_point_wise = true;
      Lookup[FAILURE_MAX_PLASTIC_STRAIN_THRESHOLD][itype] = utils::numeric(FLERR, arg[ioffset + 1],false,lmp);

      if (comm->me == 0) {
        utils::logmesg(lmp, "\nmaximum plastic strain failure criterion\n");
        utils::logmesg(lmp, "{:60} : {}\n", "failure occurs when plastic strain reaches limit",
                       Lookup[FAILURE_MAX_PLASTIC_STRAIN_THRESHOLD][itype]);
      }
    } // end maximum plastic strain failure criterion
    else if (strcmp(arg[ioffset], "*FAILURE_MAX_PAIRWISE_STRAIN") == 0) {

      /*
       * failure criterion based on maximum strain between a pair of TLSPH particles.
       */

      if (comm->me == 0) utils::logmesg(lmp, "reading *FAILURE_MAX_PAIRWISE_STRAIN\n");

      if (update_method != UPDATE_NONE) {
        error->all(FLERR, "cannot use *FAILURE_MAX_PAIRWISE_STRAIN with updated Total-Lagrangian formalism");
      }

      t = std::string("*");
      iNextKwd = -1;
      for (iarg = ioffset + 1; iarg < narg; iarg++) {
        s = std::string(arg[iarg]);
        if (s.compare(0, t.length(), t) == 0) {
          iNextKwd = iarg;
          break;
        }
      }

      if (iNextKwd < 0)
        error->all(FLERR, "no *KEYWORD terminates *FAILURE_MAX_PAIRWISE_STRAIN");

      if (iNextKwd - ioffset != 1 + 1)
        error->all(FLERR, "expected 1 arguments following *FAILURE_MAX_PAIRWISE_STRAIN but got {}\n", iNextKwd - ioffset - 1);

      failureModel[itype].failure_max_pairwise_strain = true;
      failureModel[itype].integration_point_wise = true;
      Lookup[FAILURE_MAX_PAIRWISE_STRAIN_THRESHOLD][itype] = utils::numeric(FLERR, arg[ioffset + 1],false,lmp);

      if (comm->me == 0) {
        utils::logmesg(lmp, "\nmaximum pairwise strain failure criterion\n");
        utils::logmesg(lmp, "{:60} : {}\n", "failure occurs when pairwise strain reaches limit",
               Lookup[FAILURE_MAX_PAIRWISE_STRAIN_THRESHOLD][itype]);
      }
    } // end pair based maximum strain failure criterion
    else if (strcmp(arg[ioffset], "*FAILURE_MAX_PRINCIPAL_STRAIN") == 0) {
      //error->all(FLERR, "this failure model is currently unsupported");

      /*
       * maximum principal strain failure criterion
       */
      if (comm->me == 0) utils::logmesg(lmp, "reading *FAILURE_MAX_PRINCIPAL_STRAIN\n");

      t = std::string("*");
      iNextKwd = -1;
      for (iarg = ioffset + 1; iarg < narg; iarg++) {
        s = std::string(arg[iarg]);
        if (s.compare(0, t.length(), t) == 0) {
          iNextKwd = iarg;
          break;
        }
      }

      if (iNextKwd < 0)
        error->all(FLERR, "no *KEYWORD terminates *FAILURE_MAX_PRINCIPAL_STRAIN");

      if (iNextKwd - ioffset != 1 + 1)
        error->all(FLERR, "expected 1 arguments following *FAILURE_MAX_PRINCIPAL_STRAIN but got {}\n", iNextKwd - ioffset - 1);

      failureModel[itype].failure_max_principal_strain = true;
      failureModel[itype].integration_point_wise = true;
      Lookup[FAILURE_MAX_PRINCIPAL_STRAIN_THRESHOLD][itype] = utils::numeric(FLERR, arg[ioffset + 1],false,lmp);

      if (comm->me == 0) {
        utils::logmesg(lmp, "\nmaximum principal strain failure criterion\n");
        utils::logmesg(lmp, "{:60} : {}\n", "failure occurs when principal strain reaches limit",
               Lookup[FAILURE_MAX_PRINCIPAL_STRAIN_THRESHOLD][itype]);
      }
    } // end maximum principal strain failure criterion
    else if (strcmp(arg[ioffset], "*FAILURE_JOHNSON_COOK") == 0) {
      error->all(FLERR, "this failure model is currently unsupported");
      if (comm->me == 0) utils::logmesg(lmp, "reading *FAILURE_JOHNSON_COOK\n");

      t = std::string("*");
      iNextKwd = -1;
      for (iarg = ioffset + 1; iarg < narg; iarg++) {
        s = std::string(arg[iarg]);
        if (s.compare(0, t.length(), t) == 0) {
          iNextKwd = iarg;
          break;
        }
      }

      if (iNextKwd < 0)
        error->all(FLERR, "no *KEYWORD terminates *FAILURE_JOHNSON_COOK");

      if (iNextKwd - ioffset != 5 + 1)
        error->all(FLERR, "expected 5 arguments following *FAILURE_JOHNSON_COOK but got {}\n", iNextKwd - ioffset - 1);

      failureModel[itype].failure_johnson_cook = true;
      failureModel[itype].integration_point_wise = true;

      Lookup[FAILURE_JC_D1][itype] = utils::numeric(FLERR, arg[ioffset + 1],false,lmp);
      Lookup[FAILURE_JC_D2][itype] = utils::numeric(FLERR, arg[ioffset + 2],false,lmp);
      Lookup[FAILURE_JC_D3][itype] = utils::numeric(FLERR, arg[ioffset + 3],false,lmp);
      Lookup[FAILURE_JC_D4][itype] = utils::numeric(FLERR, arg[ioffset + 4],false,lmp);
      Lookup[FAILURE_JC_EPDOT0][itype] = utils::numeric(FLERR, arg[ioffset + 5],false,lmp);

      if (comm->me == 0) {
        utils::logmesg(lmp, "\nJohnson-Cook failure criterion\n");
        utils::logmesg(lmp, "{:60} : {}\n", "parameter d1", Lookup[FAILURE_JC_D1][itype]);
        utils::logmesg(lmp, "{:60} : {}\n", "parameter d2", Lookup[FAILURE_JC_D2][itype]);
        utils::logmesg(lmp, "{:60} : {}\n", "parameter d3", Lookup[FAILURE_JC_D3][itype]);
        utils::logmesg(lmp, "{:60} : {}\n", "parameter d4", Lookup[FAILURE_JC_D4][itype]);
        utils::logmesg(lmp, "{:60} : {}\n", "reference plastic strain rate", Lookup[FAILURE_JC_EPDOT0][itype]);
      }

    } else if (strcmp(arg[ioffset], "*FAILURE_MAX_PRINCIPAL_STRESS") == 0) {
      //error->all(FLERR, "this failure model is currently unsupported");

      /*
       * maximum principal stress failure criterion
       */

      if (comm->me == 0) utils::logmesg(lmp, "reading *FAILURE_MAX_PRINCIPAL_STRESS\n");

      t = std::string("*");
      iNextKwd = -1;
      for (iarg = ioffset + 1; iarg < narg; iarg++) {
        s = std::string(arg[iarg]);
        if (s.compare(0, t.length(), t) == 0) {
          iNextKwd = iarg;
          break;
        }
      }

      if (iNextKwd < 0)
        error->all(FLERR, "no *KEYWORD terminates *FAILURE_MAX_PRINCIPAL_STRESS");

      if (iNextKwd - ioffset != 1 + 1)
        error->all(FLERR, "expected 1 arguments following *FAILURE_MAX_PRINCIPAL_STRESS but got {}\n", iNextKwd - ioffset - 1);

      failureModel[itype].failure_max_principal_stress = true;
      failureModel[itype].integration_point_wise = true;
      Lookup[FAILURE_MAX_PRINCIPAL_STRESS_THRESHOLD][itype] = utils::numeric(FLERR, arg[ioffset + 1],false,lmp);

      if (comm->me == 0) {
        utils::logmesg(lmp, "\nmaximum principal stress failure criterion\n");
        utils::logmesg(lmp, "{:60} : {}\n", "failure occurs when principal stress reaches limit",
               Lookup[FAILURE_MAX_PRINCIPAL_STRESS_THRESHOLD][itype]);
      }
    } // end maximum principal stress failure criterion

    else if (strcmp(arg[ioffset], "*FAILURE_ENERGY_RELEASE_RATE") == 0) {
      if (comm->me == 0) utils::logmesg(lmp, "reading *FAILURE_ENERGY_RELEASE_RATE\n");

      t = std::string("*");
      iNextKwd = -1;
      for (iarg = ioffset + 1; iarg < narg; iarg++) {
        s = std::string(arg[iarg]);
        if (s.compare(0, t.length(), t) == 0) {
          iNextKwd = iarg;
          break;
        }
      }

      if (iNextKwd < 0)
        error->all(FLERR, "no *KEYWORD terminates *FAILURE_ENERGY_RELEASE_RATE");

      if (iNextKwd - ioffset != 1 + 1)
        error->all(FLERR, "expected 1 arguments following *FAILURE_ENERGY_RELEASE_RATE but got {}\n", iNextKwd - ioffset - 1);

      failureModel[itype].failure_energy_release_rate = true;
      Lookup[CRITICAL_ENERGY_RELEASE_RATE][itype] = utils::numeric(FLERR, arg[ioffset + 1],false,lmp);

      if (comm->me == 0) {
        utils::logmesg(lmp,"\ncritical energy release rate failure criterion\n");
        utils::logmesg(lmp, "{:60} : {}\n", "failure occurs when energy release rate reaches limit",
               Lookup[CRITICAL_ENERGY_RELEASE_RATE][itype]);
      }
    } // end energy release rate failure criterion

    else if (strcmp(arg[ioffset], "*THERMAL_OFF") == 0) {
      if (comm->me == 0) utils::logmesg(lmp, "reading *THERMAL_OFF\n");

      t = std::string("*");
      iNextKwd = -1;
      for (iarg = ioffset + 1; iarg < narg; iarg++) {
        s = std::string(arg[iarg]);
        if (s.compare(0, t.length(), t) == 0) {
          iNextKwd = iarg;
          break;
        }
      }

      if (iNextKwd < 0)
        error->all(FLERR, "no *KEYWORD terminates *THERMAL_OFF");

      if (iNextKwd - ioffset != 1 + 0)
        error->all(FLERR, "expected 0 arguments following *THERMAL_OFF but got {}\n", iNextKwd - ioffset - 1);

      contact_thermal[itype] = false;
      //Lookup[CRITICAL_ENERGY_RELEASE_RATE][itype] = utils::numeric(FLERR, arg[ioffset + 1],false,lmp);

      if (comm->me == 0) {
        utils::logmesg(lmp,"\n thermal contact across undamaged pairs no longer active\n");
        //utils::logmesg(lmp, "{:60} : {}\n", "failure occurs when energy release rate reaches limit",
               //Lookup[CRITICAL_ENERGY_RELEASE_RATE][itype]);
      }
    } // end contact thermal criteria

    else if (strcmp(arg[ioffset], "*CONTACT_MECHANICAL") == 0) {
      if (comm->me == 0) utils::logmesg(lmp, "reading *CONTACT_MECHANICAL\n");

      t = std::string("*");
      iNextKwd = -1;
      for (iarg = ioffset + 1; iarg < narg; iarg++) {
        s = std::string(arg[iarg]);
        if (s.compare(0, t.length(), t) == 0) {
          iNextKwd = iarg;
          break;
        }
      }

      if (iNextKwd < 0)
        error->all(FLERR, "no *KEYWORD terminates *CONTACT_MECHANICAL");

      if (iNextKwd - ioffset != 1 + 4)
        error->all(FLERR, "expected 4 arguments (DISTANCE_SCALING, FORCE_SCALING, DAMPING_SCALING, RESTITUTION) following *CONTACT_MECHANICAL but got {}\n", iNextKwd - ioffset - 1);

      contact_mechanical[itype] = true;
      Lookup[CONTACT_DISTANCE_SCALING][itype] = utils::numeric(FLERR, arg[ioffset + 1],false,lmp);
      Lookup[CONTACT_FORCE_SCALING][itype] = utils::numeric(FLERR, arg[ioffset + 2],false,lmp);
      Lookup[CONTACT_DAMPING_SCALING][itype] = utils::numeric(FLERR, arg[ioffset + 3],false,lmp);
      Lookup[CONTACT_RESTITUTION][itype] = utils::numeric(FLERR, arg[ioffset + 4],false,lmp);

      if (comm->me == 0) {
        utils::logmesg(lmp,"\n mechanical contact across damaged pairs active\n");
        utils::logmesg(lmp, "{:60} : {}\n", "contact distance scaleing used is",
               Lookup[CONTACT_DISTANCE_SCALING][itype]);
        utils::logmesg(lmp, "{:60} : {}\n", "contact force scaleing used is",
               Lookup[CONTACT_FORCE_SCALING][itype]);
        utils::logmesg(lmp, "{:60} : {}\n", "contact force damping scaleing used is",
               Lookup[CONTACT_DAMPING_SCALING][itype]);
        utils::logmesg(lmp, "{:60} : {}\n", "contact force ceoficient of restitution is",
               Lookup[CONTACT_RESTITUTION][itype]);
      }
    } // end contact mechanical criteria

    else error->all(FLERR, "unknown *KEYWORD: {}", arg[ioffset]);
  }
  setflag[itype][itype] = 1;
}

/* ----------------------------------------------------------------------
 init for one type pair i,j and corresponding j,i
 ------------------------------------------------------------------------- */

double PairTlsph::init_one(int i, int j) {

  if (!allocated)
    allocate();

  if (setflag[i][j] == 0)
    error->all(FLERR, "All pair coeffs are not set");

  if (force->newton == 1)
    error->all(FLERR, "Pair style tlsph requires newton off");

// cutoff = sum of max I,J radii for
// dynamic/dynamic & dynamic/frozen interactions, but not frozen/frozen

  double cutoff = maxrad_dynamic[i] + maxrad_dynamic[j];
  cutoff = MAX(cutoff, maxrad_frozen[i] + maxrad_dynamic[j]);
  cutoff = MAX(cutoff, maxrad_dynamic[i] + maxrad_frozen[j]);
  return cutoff;
}

/* ----------------------------------------------------------------------
 init specific to this pair style
 ------------------------------------------------------------------------- */

void PairTlsph::init_style() {
  int i;

  if (force->newton_pair == 1) {
    error->all(FLERR, "Pair style tlsph requires newton pair off");
  }

// request a granular neighbor list
  neighbor->add_request(this, NeighConst::REQ_SIZE);

// set maxrad_dynamic and maxrad_frozen for each type
// include future Fix pour particles as dynamic

  for (i = 1; i <= atom->ntypes; i++)
    onerad_dynamic[i] = onerad_frozen[i] = 0.0;

  double *radius = atom->radius;
  int *type = atom->type;
  int nlocal = atom->nlocal;

  for (i = 0; i < nlocal; i++)
    onerad_dynamic[type[i]] = MAX(onerad_dynamic[type[i]], radius[i]);

  MPI_Allreduce(&onerad_dynamic[1], &maxrad_dynamic[1], atom->ntypes, MPI_DOUBLE, MPI_MAX, world);
  MPI_Allreduce(&onerad_frozen[1], &maxrad_frozen[1], atom->ntypes, MPI_DOUBLE, MPI_MAX, world);

// if first init, create Fix needed for storing reference configuration neighbors

  int igroup = group->find("tlsph");
  if (igroup == -1)
    error->all(FLERR, "Pair style tlsph requires its particles to be part of a group named tlsph. This group does not exist.");

  if (fix_tlsph_reference_configuration == nullptr) {
    auto fixarg = new char*[3];
    fixarg[0] = (char *) "SPOOFF_TLSPH_NEIGHBORS";
    fixarg[1] = (char *) "tlsph";
    fixarg[2] = (char *) "SPOOFF_TLSPH_NEIGHBORS";
    modify->add_fix(3, fixarg);
    delete[] fixarg;
    fix_tlsph_reference_configuration = dynamic_cast<FixSPOOFF_TLSPH_ReferenceConfiguration *>(modify->fix[modify->nfix - 1]);
    fix_tlsph_reference_configuration->pair = this;
  }

// find associated SPOOFF_TLSPH_NEIGHBORS fix that must exist
// could have changed locations in fix list since created

  ifix_tlsph = -1;
  for (int i = 0; i < modify->nfix; i++)
    if (strcmp(modify->fix[i]->style, "SPOOFF_TLSPH_NEIGHBORS") == 0)
      ifix_tlsph = i;
  if (ifix_tlsph == -1)
    error->all(FLERR, "Fix SPOOFF_TLSPH_NEIGHBORS does not exist");

}

/* ----------------------------------------------------------------------
 neighbor callback to inform pair style of neighbor list to use
 optional granular history list
 ------------------------------------------------------------------------- */

void PairTlsph::init_list(int id, class NeighList *ptr) {
  if (id == 0) list = ptr;
}

/* ----------------------------------------------------------------------
 memory usage of local atom-based arrays
 ------------------------------------------------------------------------- */

double PairTlsph::memory_usage() {
  return 118.0 * nmax * sizeof(double); 
}

/* ----------------------------------------------------------------------
 extract method to provide access to this class' data structures
 ------------------------------------------------------------------------- */

void *PairTlsph::extract(const char *str, int &/*i*/) {
  if (strcmp(str, "spooff/tlsph/Fincr_ptr") == 0) {
    return (void *) Fincr;
  } else if (strcmp(str, "spooff/tlsph/detF_ptr") == 0) {
    return (void *) detF;
  } else if (strcmp(str, "spooff/tlsph/PK1_ptr") == 0) {
    return (void *) PK1;
  } else if (strcmp(str, "spooff/tlsph/smoothVel_ptr") == 0) {
    return (void *) smoothVelDifference;
  } else if (strcmp(str, "spooff/tlsph/numNeighsRefConfig_ptr") == 0) {
    return (void *) numNeighsRefConfig;
  } else if (strcmp(str, "spooff/tlsph/stressTensor_ptr") == 0) {
    return (void *) CauchyStress;
  } else if (strcmp(str, "spooff/tlsph/updateFlag_ptr") == 0) {
    return (void *) &updateFlag;
  } else if (strcmp(str, "spooff/tlsph/strain_rate_ptr") == 0) {
    return (void *) D;
  } else if (strcmp(str, "spooff/tlsph/hMin_ptr") == 0) {
    return (void *) &hMin;
  } else if (strcmp(str, "spooff/tlsph/dtCFL_ptr") == 0) {
    return (void *) &dtCFL;
  } else if (strcmp(str, "spooff/tlsph/dtRelative_ptr") == 0) {
    return (void *) &dtRelative;
  } else if (strcmp(str, "spooff/tlsph/hourglass_error_ptr") == 0) {
    return (void *) hourglass_error;
  } else if (strcmp(str, "spooff/tlsph/particle_dt_ptr") == 0) {
    return (void *) particle_dt;
  } else if (strcmp(str, "spooff/tlsph/rotation_ptr") == 0) {
    return (void *) R;
  }

  return nullptr;
}

/* ---------------------------------------------------------------------- */

int PairTlsph::pack_forward_comm(int n, int *list, double *buf, int /*pbc_flag*/, int * /*pbc*/) {
  int i, j, m;
  tagint *mol = atom->molecule;
  double *damage = atom->damage;
  double *cpsph = atom->cpsph;
  double *ksph = atom->ksph;
  double *fracsph = atom->fracsph;
  double *eff_plastic_strain = atom->eff_plastic_strain;
  double *eff_plastic_strain_rate = atom->eff_plastic_strain_rate;
//assume hgsph not needed in forward com as not used in pair calculations.
  m = 0;
  for (i = 0; i < n; i++) {
    j = list[i];
    buf[m++] = PK1[j](0, 0); // PK1 is not symmetric
    buf[m++] = PK1[j](0, 1);
    buf[m++] = PK1[j](0, 2);
    buf[m++] = PK1[j](1, 0);
    buf[m++] = PK1[j](1, 1);
    buf[m++] = PK1[j](1, 2);
    buf[m++] = PK1[j](2, 0);
    buf[m++] = PK1[j](2, 1);
    buf[m++] = PK1[j](2, 2); // 9

    buf[m++] = Fincr[j](0, 0); // Fincr is not symmetric
    buf[m++] = Fincr[j](0, 1);
    buf[m++] = Fincr[j](0, 2);
    buf[m++] = Fincr[j](1, 0);
    buf[m++] = Fincr[j](1, 1);
    buf[m++] = Fincr[j](1, 2);
    buf[m++] = Fincr[j](2, 0);
    buf[m++] = Fincr[j](2, 1);
    buf[m++] = Fincr[j](2, 2); // 9 + 9 = 18

    buf[m++] = mol[j]; //19
    buf[m++] = damage[j]; //20
    buf[m++] = eff_plastic_strain[j]; //21
    buf[m++] = eff_plastic_strain_rate[j]; //22
    buf[m++] = cpsph[j]; //23
    buf[m++] = ksph[j]; //24
    buf[m++] = fracsph[j]; //25



  }
  return m;
}

/* ---------------------------------------------------------------------- */

void PairTlsph::unpack_forward_comm(int n, int first, double *buf) {
  int i, m, last;
  tagint *mol = atom->molecule;
  double *damage = atom->damage;
  double *cpsph = atom->cpsph;
  double *ksph = atom->ksph;
  double *fracsph = atom->fracsph;
  double *eff_plastic_strain = atom->eff_plastic_strain;
  double *eff_plastic_strain_rate = atom->eff_plastic_strain_rate;

  m = 0;
  last = first + n;
  for (i = first; i < last; i++) {

    PK1[i](0, 0) = buf[m++]; // PK1 is not symmetric
    PK1[i](0, 1) = buf[m++];
    PK1[i](0, 2) = buf[m++];
    PK1[i](1, 0) = buf[m++];
    PK1[i](1, 1) = buf[m++];
    PK1[i](1, 2) = buf[m++];
    PK1[i](2, 0) = buf[m++];
    PK1[i](2, 1) = buf[m++];
    PK1[i](2, 2) = buf[m++];

    Fincr[i](0, 0) = buf[m++];
    Fincr[i](0, 1) = buf[m++];
    Fincr[i](0, 2) = buf[m++];
    Fincr[i](1, 0) = buf[m++];
    Fincr[i](1, 1) = buf[m++];
    Fincr[i](1, 2) = buf[m++];
    Fincr[i](2, 0) = buf[m++];
    Fincr[i](2, 1) = buf[m++];
    Fincr[i](2, 2) = buf[m++];

    mol[i] = static_cast<int>(buf[m++]);
    damage[i] = buf[m++];
    eff_plastic_strain[i] = buf[m++]; //22
    eff_plastic_strain_rate[i] = buf[m++]; //23
    cpsph[i] = buf[m++]; //24
    ksph[i] = buf[m++];  //25
    fracsph[i] = buf[m++]; //26
  }
}

/* ----------------------------------------------------------------------
 compute effective P-wave speed
 determined by longitudinal modulus
 ------------------------------------------------------------------------- */

void PairTlsph::effective_longitudinal_modulus(const int itype, const double shear, const double bulk, const double m_modulus, const double dt, const double d_iso, const double p_rate,
                                               const Matrix3d& d_dev, const Matrix3d& sigma_dev_rate, const double /*damage*/, double &K_eff, double &mu_eff, double &M_eff) {
  double M0; // initial longitudinal modulus
  double shear_rate_sq;

  M0 = m_modulus;

  if (dt * d_iso > 1.0e-6) {
    K_eff = p_rate / d_iso;
    if (K_eff < 0.0) { // it is possible for K_eff to become negative due to strain softening
      K_eff = bulk;
    }
  } else {
    K_eff = bulk;
  }

  if (domain->dimension == 3) {
// Calculate 2 mu by looking at ratio shear stress / shear strain. Use numerical softening to avoid divide-by-zero.
    mu_eff = 0.5
      * (sigma_dev_rate(0, 1) / (d_dev(0, 1) + 1.0e-16) + sigma_dev_rate(0, 2) / (d_dev(0, 2) + 1.0e-16)
         + sigma_dev_rate(1, 2) / (d_dev(1, 2) + 1.0e-16));

// Calculate magnitude of deviatoric strain rate. This is used for deciding if shear modulus should be computed from current rate or be taken as the initial value.
    shear_rate_sq = d_dev(0, 1) * d_dev(0, 1) + d_dev(0, 2) * d_dev(0, 2) + d_dev(1, 2) * d_dev(1, 2);
  } else {
    mu_eff = 0.5 * (sigma_dev_rate(0, 1) / (d_dev(0, 1) + 1.0e-16));
    shear_rate_sq = d_dev(0, 1) * d_dev(0, 1);
  }

  if (dt * dt * shear_rate_sq < 1.0e-8) {
    mu_eff = shear;
  }

  if (mu_eff < shear) { // it is possible for mu_eff to become negative due to strain softening
    mu_eff = shear;
  }

  if (mu_eff < 0.0) {
    error->one(FLERR, "mu_eff = {}, tau={}, gamma={}", mu_eff, sigma_dev_rate(0, 1), d_dev(0, 1));

  }

  M_eff = (K_eff + 4.0 * mu_eff / 3.0); // effective dilational modulus, see Pronto 2d eqn 3.4.8

  if (M_eff < M0) { // do not allow effective dilatational modulus to decrease beyond its initial value
    M_eff = M0;
  }
}

/* ----------------------------------------------------------------------
 compute pressure. Called from AssembleStress().
 ------------------------------------------------------------------------- */
void PairTlsph::ComputePressure(const int i, const double rho, const double mass_specific_energy, const double temperature, const double vol_specific_energy,
                                const double pInitial, const double d_iso, double &pFinal, double &p_rate) {
  int *type = atom->type;
  double dt = update->dt;
  double *youngs = atom->youngs;
  double *poissons = atom->poissons;
  double *linear_expansion = atom->linear_expansion;

  int itype;

  double shear, bulk, m_modulus, lame_lambda;
         shear = youngs[i] / (2.0 * (1.0 + poissons[i]));
         lame_lambda = youngs[i] * poissons[i]  / ((1.0 + poissons[i]) * (1.0 - 2.0 * poissons[i]));
         m_modulus = lame_lambda + 2.0 * shear;
         bulk = lame_lambda + 2.0 * shear / 3.0;

  itype = type[i];

  switch (eos[itype]) {
  case EOS_LINEAR:
    LinearEOS(bulk, pInitial, d_iso, dt, pFinal, p_rate);
    break;
  case EOS_THERMAL:
    ThermalEOS(rho, Lookup[REFERENCE_DENSITY][itype], mass_specific_energy, temperature, Lookup[EOS_THERMAL_T0][itype], bulk, linear_expansion[i], Lookup[HEAT_CAPACITY][itype], pInitial, dt,
                pFinal, p_rate);
    break;
  case EOS_NONE:
    pFinal = 0.0;
    p_rate = 0.0;
    break;
  case EOS_SHOCK:
//  rho,  rho0,  e,  e0,  c0,  S,  Gamma,  pInitial,  dt,  &pFinal,  &p_rate);
    ShockEOS(rho, Lookup[REFERENCE_DENSITY][itype], mass_specific_energy, 0.0, Lookup[EOS_SHOCK_C0][itype],
             Lookup[EOS_SHOCK_S][itype], Lookup[EOS_SHOCK_GAMMA][itype], pInitial, dt, pFinal, p_rate);
    break;
  case EOS_POLYNOMIAL:
    polynomialEOS(rho, Lookup[REFERENCE_DENSITY][itype], vol_specific_energy, Lookup[EOS_POLYNOMIAL_C0][itype],
                  Lookup[EOS_POLYNOMIAL_C1][itype], Lookup[EOS_POLYNOMIAL_C2][itype], Lookup[EOS_POLYNOMIAL_C3][itype],
                  Lookup[EOS_POLYNOMIAL_C4][itype], Lookup[EOS_POLYNOMIAL_C5][itype], Lookup[EOS_POLYNOMIAL_C6][itype], pInitial, dt,
                  pFinal, p_rate);

    break;
  default:
    error->one(FLERR, "unknown EOS.");
    break;
  }
}

/* ----------------------------------------------------------------------
 Compute stress deviator. Called from AssembleStress().
 ------------------------------------------------------------------------- */
void PairTlsph::ComputeStressDeviator(const int i, const Matrix3d& sigmaInitial_dev, const Matrix3d& d_dev, Matrix3d &sigmaFinal_dev,
                                      Matrix3d &sigma_dev_rate, double &plastic_strain_increment) {
  double *eff_plastic_strain = atom->eff_plastic_strain;
  double *eff_plastic_strain_rate = atom->eff_plastic_strain_rate;
  int *type = atom->type;
  double *rmass = atom->rmass;
  double *esph = atom->esph;
  double *fracsph = atom->fracsph;
  double *yeild = atom->yeild;
  double *pesph = atom->pesph;
  double *youngs = atom->youngs;
  double *poissons = atom->poissons;
  double dt = update->dt;
  double yieldStress;
  int itype;

  double mass_specific_energy = (esph[i]) / rmass[i]; // internal energy per unit mass (not included elastic potential energy)
  plastic_strain_increment = 0.0;
  itype = type[i];

  double shear, bulk, m_modulus, lame_lambda;
         shear = youngs[i] / (2.0 * (1.0 + poissons[i]));
         lame_lambda = youngs[i] * poissons[i]  / ((1.0 + poissons[i]) * (1.0 - 2.0 * poissons[i]));
         m_modulus = lame_lambda + 2.0 * shear;
         bulk = lame_lambda + 2.0 * shear / 3.0;

  switch (strengthModel[itype]) {
  case STRENGTH_LINEAR:

    sigma_dev_rate = 2.0 * shear * d_dev;
    sigmaFinal_dev = sigmaInitial_dev + dt * sigma_dev_rate;

    break;
  case LINEAR_DEFGRAD:
//LinearStrengthDefgrad(Lookup[LAME_LAMBDA][itype], Lookup[SHEAR_MODULUS][itype], Fincr[i], &sigmaFinal_dev);
//eff_plastic_strain[i] = 0.0;
//p_rate = pInitial - sigmaFinal_dev.trace() / 3.0;
//sigma_dev_rate = sigmaInitial_dev - Deviator(sigmaFinal_dev);
    error->one(FLERR, "LINEAR_DEFGRAD is only for debugging purposes and currently deactivated.");
    R[i].setIdentity();
    break;
  case STRENGTH_LINEAR_PLASTIC:

    yieldStress = yeild[i] + Lookup[HARDENING_PARAMETER][itype] * eff_plastic_strain[i];
    LinearPlasticStrength(shear, yieldStress, sigmaInitial_dev, d_dev, dt, sigmaFinal_dev,
                          sigma_dev_rate, plastic_strain_increment);
    break;
  case STRENGTH_JOHNSON_COOK:
    JohnsonCookStrength(shear, Lookup[HEAT_CAPACITY][itype], mass_specific_energy, Lookup[JC_A][itype],
                        Lookup[JC_B][itype], Lookup[JC_a][itype], Lookup[JC_C][itype], Lookup[JC_epdot0][itype], Lookup[JC_T0][itype],
                        Lookup[JC_Tmelt][itype], Lookup[JC_M][itype], dt, eff_plastic_strain[i], eff_plastic_strain_rate[i],
                        sigmaInitial_dev, d_dev, sigmaFinal_dev, sigma_dev_rate, plastic_strain_increment);
    break;
  case STRENGTH_NONE:
    sigmaFinal_dev.setZero();
    sigma_dev_rate.setZero();
    break;
  default:
    error->one(FLERR, "unknown strength model.");
    break;
  }

}

/* ----------------------------------------------------------------------
 Compute damage. Called from AssembleStress().
 ------------------------------------------------------------------------- */
void PairTlsph::ComputeDamage(const int i, const Matrix3d& strain, const Matrix3d& stress, Matrix3d &/*stress_damaged*/) {
  double *eff_plastic_strain = atom->eff_plastic_strain;
  double *eff_plastic_strain_rate = atom->eff_plastic_strain_rate;
  double *radius = atom->radius;
  double *damage = atom->damage;
  double *fracsph = atom->fracsph;
  double *youngs = atom->youngs;
  int *type = atom->type;
  int itype = type[i];
  double jc_failure_strain;
  Matrix3d eye, stress_deviator;

  eye.setIdentity();
  stress_deviator = Deviator(stress);
  double pressure = -stress.trace() / 3.0;

  if (failureModel[itype].failure_max_principal_stress) {
    //error->one(FLERR, "not yet implemented");
    /*
     * maximum stress failure criterion:
     */
    if(IsotropicMaxStressDamage(stress, fracsph[i]*youngs[i])){
      damage[i] = 1.0;
    }else{
      damage[i] = 0.0;
    }
  } else if (failureModel[itype].failure_max_principal_strain) {
    //error->one(FLERR, "not yet implemented");
    /*
     * maximum strain failure criterion:
     */
    ;
    if(IsotropicMaxStrainDamage(strain, Lookup[FAILURE_MAX_PRINCIPAL_STRAIN_THRESHOLD][itype])){
      damage[i] = 1.0;
    }else{
      damage[i] = 0.0;
    }
  } else if (failureModel[itype].failure_max_plastic_strain) {
    if (eff_plastic_strain[i] >= fracsph[i]) { //Changed by MH
      damage[i] = 1.0;
    }else{
      damage[i] = 0.0;
    }
  } else if (failureModel[itype].failure_johnson_cook) {

    jc_failure_strain = JohnsonCookFailureStrain(pressure, stress_deviator, Lookup[FAILURE_JC_D1][itype],
                                                 Lookup[FAILURE_JC_D2][itype], Lookup[FAILURE_JC_D3][itype], Lookup[FAILURE_JC_D4][itype],
                                                 Lookup[FAILURE_JC_EPDOT0][itype], eff_plastic_strain_rate[i]);

    if (eff_plastic_strain[i] >= jc_failure_strain) {
      double damage_rate = Lookup[SIGNAL_VELOCITY][itype] / (100.0 * radius[i]);
      damage[i] += damage_rate * update->dt;
    }
  }
}

