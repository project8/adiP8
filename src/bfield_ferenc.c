using namespace std;
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <cstdio>
#include <iostream>
#include "paramanage.h"
#include "bfield_ferenc.h"

////////////////////////////////////////////////////////////////////////
//
//  Field calculation
//
////////////////////////////////////////////////////////////////////////

void get_bfield_fg(double *posvec, double *bvec)
{
  double z, r, A, Bz, Br;
  double xpos, ypos, zpos;

  //  xpos = posvec[0]-parameter.mag_x_offset;
  //  ypos = posvec[1]-parameter.mag_y_offset;
  //  zpos = posvec[2]-parameter.mag_z_offset;
  xpos = posvec[0];
  ypos = posvec[1];
  zpos = posvec[2];

  z = xpos / 100.0;
  r = sqrt(ypos * ypos + zpos * zpos) / 100.0;
  // posvec comes in cm units -> times 100 to have meters

  magfield(z, r, &A, &Bz, &Br);
  // call magfield of Ferenc Glck with converted values

  bvec[0] = Bz;
  if (Bz < 1e-12 || r < 1e-12) {
    bvec[1] = bvec[2] = 0.;
  } else {
    bvec[2] = posvec[2] * Br / (r * 100.);
    bvec[1] = posvec[1] * Br / (r * 100.);
  }

  //cout << z << " " << r << " " << Bz << " "<< Br << endl << flush;
  //cout << posvec[0] << " " << posvec[1] << " " << posvec[2] << " " << bvec[0] << " " << bvec[1] << " " << bvec[2] << endl << flush;
}

///////////////////////////////////////////////////////////////////////


void magfield(double z, double r, double *A, double *Bz, double *Br)
// This subroutine computes the axial and radial magnetic
//    field components Bz and Br, and the vectorpotential A
//    in a point with cylindrical coordinates z and r.
// Method: Legendre polynomial expansion around the source point with index k.
//  Nspmax: we need Nspmax>=Nsp
//  nmax: maximal index of the source constants (maximum of n).
//  Nspmax and nmax are given by global #define commands.
//  Important (if magfield is used separated from magsource):
//    the same nmax number used for magsource
//    should also be used for magfield !!!
{
  const int Nspmax = 20001;
  const int nmax = 200;

  FILE *fp;
  int kloop, kx, n, nmaxtest, k, kk;
  static double z0[Nspmax + 1], rocen[Nspmax + 1], B[Nspmax + 1][nmax + 1];
  static double c1[nmax + 1], c2[nmax + 1], c3[nmax + 1], c4[nmax + 1];
  static double c5[nmax + 1], c6[nmax + 1];
  static int Nsp, klast;
  static int iff = 0;
  double ro, u, delz, s, rcmin;
  double P[nmax + 1], Pp[nmax + 1], rc, rcn;
  int iA, iBz, iBr;
  double Aplus[nmax + 1], Bzplus[nmax + 1], Brplus[nmax + 1];
  double A1, Bz1, Br1;
  char mag_filename[255];       // filename for magsource data

  strcpy(mag_filename, parameter.filename);
  strcat(mag_filename, ".magsource"); // data input from file with extension .ini
// Input from file magsource.dat:
  if (iff == 0) {
    fp = fopen(mag_filename, "r");
    if (fp ==NULL) {
      cout << "Can't open magsource file " << mag_filename << endl;
    }
    fscanf(fp, "%i %i", &Nsp, &nmaxtest);
    if (Nsp > Nspmax || nmaxtest != nmax) {
      printf("Message from subroutine magfield: Nsp > Nspmax or different nmax values used in magsource and magfield !!! Computation is  stopped !!! \n\n");
      exit(0);
    }
    for (kloop = 1; kloop <= Nsp; kloop++) {
      fscanf(fp, "%i", &kx);
      fscanf(fp, "%le %le", &z0[kloop], &rocen[kloop]);
      for (n = 0; n <= nmax; n++) {
        fscanf(fp, "%le", &B[kloop][n]);
      }
    }
    fclose(fp);
// Initialization of c1,c2,c3,c4,c5,c6 vectors:
    for (n = 2; n <= nmax; n++) {
      c1[n] = (2. * n - 1.) / (1. * n);
      c2[n] = (n - 1.) / (1. * n);
      c3[n] = (2. * n - 1.) / (1. * (n - 1.));
      c4[n] = (1. * n) / (1. * (n - 1.));
      c5[n] = (1.) / (n * (n + 1.));
      c6[n] = (1.) / (n + 1.);
    }
// The best source point is searched here (with minimal
//    convergence ratio rc):
    rcmin = 1.e20;
    for (k = 1; k <= Nsp; k++) {
      delz = z - z0[k];
      ro = sqrt(r * r + delz * delz);
      rc = ro / rocen[k];
      if (rc < rcmin) {
        rcmin = rc;
        klast = k;
      }
    }
// End of source point searching
    iff = 1;
  }
//
// The best source point is searched here
//   (starting from the last source point)
  k = klast;
  delz = z - z0[k];
  ro = sqrt(r * r + delz * delz);
  rcmin = ro / rocen[k];
  kk = k + 1;
  if (kk <= Nsp) {
    delz = z - z0[kk];
    ro = sqrt(r * r + delz * delz);
    rc = ro / rocen[kk];
    if (rc < rcmin) {
      rcmin = rc;
      k = kk;
    }
  }
  kk = klast - 1;
  if (kk >= 1) {
    delz = z - z0[kk];
    ro = sqrt(r * r + delz * delz);
    rc = ro / rocen[kk];
    if (rc < rcmin) {
      k = kk;
    }
  }
  klast = k;
// If rc>0.999: new searching:
  delz = z - z0[k];
  ro = sqrt(r * r + delz * delz);
  rc = ro / rocen[k];
  if (rc > 0.999) {
    rcmin = 1.e20;
    for (k = 1; k <= Nsp; k++) {
      delz = z - z0[k];
      ro = sqrt(r * r + delz * delz);
      rc = ro / rocen[k];
      if (rc < rcmin) {
        rcmin = rc;
        klast = k;
      }
    }
    k = klast;
  }
// End of source point searching
//////////////////////////////////////
// If the field point is very close to the source point:
  if (r < 1.e-12 && fabs(z - z0[k]) < 1.e-12) {
    *A = 0.;
    *Bz = B[k][0];
    *Br = 0.;
    return;
  }
// ro,u,s,rc,rcn:
  delz = z - z0[k];
  ro = sqrt(r * r + delz * delz);
  u = delz / ro;
  s = r / ro;
  rc = ro / rocen[k];           // convergence ratio
  rcn = rc;
// First 2 terms of Legendre polynomial P and its derivative Pp (P-primed)
  P[0] = 1.;
  P[1] = u;
  Pp[0] = 0.;
  Pp[1] = 1.;
// If rc>0.999: computation is stopped
//   (the series is not convergent)
  if (rc > 0.999) {
    printf("Message from subroutine magfield: Convergence ratio ro/rocen is larger than 0.999 !!! Computation is  stopped !!! \n\n");
    exit(0);
  }
// First 2 terms of the series:
  *A = s * rocen[k] * B[k][0] / 2. * rc;
  *Bz = B[k][0] + B[k][1] * rc * u;
  *Br = -s * B[k][1] / 2. * rc;
//
  iA = 0;
  iBz = 0;
  iBr = 0;
  Aplus[1] = 1.e30;
  Bzplus[1] = 1.e30;
  Brplus[1] = 1.e30;
// We start here the series expansion:
  for (n = 2; n <= nmax - 1; n++) {
    rcn *= rc;
    P[n] = c1[n] * u * P[n - 1] - c2[n] * P[n - 2];
    Pp[n] = c3[n] * u * Pp[n - 1] - c4[n] * Pp[n - 2];
    Aplus[n] = rocen[k] * s * B[k][n - 1] * c5[n] * rcn * Pp[n];
    Bzplus[n] = B[k][n] * rcn * P[n];
    Brplus[n] = -s * B[k][n] * c6[n] * rcn * Pp[n];
    *A += Aplus[n];
    *Bz += Bzplus[n];
    *Br += Brplus[n];
    A1 = 1.e-15 * fabs(*A);
    Bz1 = 1.e-15 * fabs(*Bz);
    Br1 = 1.e-15 * fabs(*Br);
    if (n > 8) {
      if (fabs(Aplus[n]) < A1 && fabs(Aplus[n - 1]) < A1) {
        iA = 1;
      }
      if (fabs(Bzplus[n]) < Bz1 && fabs(Bzplus[n - 1]) < Bz1) {
        iBz = 1;
      }
      if (fabs(Brplus[n]) < Br1 && fabs(Brplus[n - 1]) < Br1) {
        iBr = 1;
      }
      if (fabs(*Br) < 1.e-12) {
        iBr = 1;
      }
    }
    if (iA * iBz * iBr == 1) {
      break;
    }
  }
  return;
}

////////////////////////////////////////////////////////
