/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under 
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing author: Mario Orsi (U Southampton), orsimario@gmail.com
------------------------------------------------------------------------- */

#include "math.h"
#include "stdlib.h"
#include "bond_dipole.h"
#include "atom.h"
#include "neighbor.h"
#include "domain.h"
#include "comm.h"
#include "force.h"
#include "memory.h"
#include "error.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

BondDipole::BondDipole(LAMMPS *lmp) : Bond(lmp) {}

/* ---------------------------------------------------------------------- */

BondDipole::~BondDipole()
{
  if (allocated) {
    memory->destroy(setflag);
    memory->destroy(k);
    memory->destroy(gamma0);
  }
}

/* ---------------------------------------------------------------------- */

void BondDipole::compute(int eflag, int vflag)
/* ----------------------------------------------------------------------
   This function implements an intramolecular pair interaction between a
   dipolar atom 'iDip' and a reference atom 'iRef'.
   A torque is applied on 'iDip' to restrain the dipole orientation along
   the direction defined by the vector 'x[iRef] - x[iDip]'.
------------------------------------------------------------------------- */
{
  int iDip,iRef,n,type;
  double delx,dely,delz,ebond,tbond;
  double r,dr,cosGamma,deltaGamma,kdg,rmu;

  ebond = 0.0;
  if (eflag || vflag) ev_setup(eflag,vflag);
  else evflag = 0;

  double **x = atom->x;
  double **mu = atom->mu;
  double **torque = atom->torque;
  int **bondlist = neighbor->bondlist;
  int nbondlist = neighbor->nbondlist;
  int nlocal = atom->nlocal;
  int newton_bond = force->newton_bond;

  for (n = 0; n < nbondlist; n++) {
    iDip = bondlist[n][0];
    iRef = bondlist[n][1];
    type = bondlist[n][2];

    delx = x[iRef][0] - x[iDip][0]; // x21
    dely = x[iRef][1] - x[iDip][1]; // y21
    delz = x[iRef][2] - x[iDip][2]; // z21
    domain->minimum_image(delx,dely,delz);

    r = sqrt(delx*delx + dely*dely + delz*delz);

    rmu = r * mu[iDip][3];
    cosGamma = (mu[iDip][0]*delx+mu[iDip][1]*dely+mu[iDip][2]*delz) / rmu;
    deltaGamma = cosGamma - cos(gamma0[type]);
    kdg = k[type] * deltaGamma;

    if (eflag) ebond = kdg * deltaGamma; // energy  
      
    tbond = 2.0 * kdg / rmu; 
      
    torque[iDip][0] += tbond * (dely*mu[iDip][2] - delz*mu[iDip][1]);
    torque[iDip][1] += tbond * (delz*mu[iDip][0] - delx*mu[iDip][2]);
    torque[iDip][2] += tbond * (delx*mu[iDip][1] - dely*mu[iDip][0]);
    
    if (evflag) 
      ev_tally(iDip,iRef,nlocal,newton_bond,ebond,0.0,delx,dely,delz);
  }
}

/* ---------------------------------------------------------------------- */

void BondDipole::allocate()
{
  allocated = 1;
  int n = atom->nbondtypes;

  memory->create(k,n+1,"bond:k");
  memory->create(gamma0,n+1,"bond:gamma0");

  memory->create(setflag,n+1,"bond:setflag");
  for (int i = 1; i <= n; i++) setflag[i] = 0;
}

/* ----------------------------------------------------------------------
   set coeffs for one or more types
------------------------------------------------------------------------- */

void BondDipole::coeff(int narg, char **arg)
{
  if (narg != 3) error->all("Incorrect args for bond coefficients");
  if (!allocated) allocate();

  int ilo,ihi;
  force->bounds(arg[0],atom->nbondtypes,ilo,ihi);

  double k_one = force->numeric(arg[1]);
  double gamma0_one = force->numeric(arg[2]);

  int count = 0;
  const double PI = 3.14159265358979323846; 
  for (int i = ilo; i <= ihi; i++) {
    k[i] = k_one;
    gamma0[i] = gamma0_one * PI / 180.0; // convert from degrees to radians
    setflag[i] = 1;
    count++;
  }

  if (count == 0) error->all("Incorrect args for bond coefficients");
}

/* ----------------------------------------------------------------------
   return an equilibrium bond length (used by SHAKE - not relevant here) 
------------------------------------------------------------------------- */

double BondDipole::equilibrium_distance(int i)
{
  return 0.0;
}

/* ----------------------------------------------------------------------
   proc 0 writes out coeffs to restart file 
------------------------------------------------------------------------- */

void BondDipole::write_restart(FILE *fp)
{
  fwrite(&k[1],sizeof(double),atom->nbondtypes,fp);
  fwrite(&gamma0[1],sizeof(double),atom->nbondtypes,fp);
}

/* ----------------------------------------------------------------------
   proc 0 reads coeffs from restart file, bcasts them 
------------------------------------------------------------------------- */

void BondDipole::read_restart(FILE *fp)
{
  allocate();

  if (comm->me == 0) {
    fread(&k[1],sizeof(double),atom->nbondtypes,fp);
    fread(&gamma0[1],sizeof(double),atom->nbondtypes,fp);
  }
  MPI_Bcast(&k[1],atom->nbondtypes,MPI_DOUBLE,0,world);
  MPI_Bcast(&gamma0[1],atom->nbondtypes,MPI_DOUBLE,0,world);

  for (int i = 1; i <= atom->nbondtypes; i++) setflag[i] = 1;
}

/* ----------------------------------------------------------------------
   used by ComputeBondLocal - not relevant here
------------------------------------------------------------------------- */

double BondDipole::single(int type, double rsq, int i, int j)
{
  return 0.0;
}
