using namespace std;

#include <cstdio>
//#include <iostream>
#include <cmath>
#include "eH2.h"
double Del(double E, double c);
double Dexc(double E, double c);
double sumexc(double K);
double sigmaBC(double E);
double sigmadiss10(double E);
double sigmadiss15(double E);
double Dinel(double E, double c);
void gensecelen(double E, double *W);
double sigmainel(double E);
void subrn(double *u, int len);
double random();
double lagrange(int n, double *xn, double *fn, double x);

//////////////////////////////////////////////////
///
///  The main scattering generation subroutine:
///
//////////////////////////////////////////////////
void eH2(double E, double s, double n, double *Eloss, double *theta, int *index) 
{
  
// This subroutine computes the electron - hydrogen molecule scattering
// probability, and it generates an
// elastic, excitation or ionization scattering event by this probability.
// Input parameters:
//    E: electron kinetic energy in eV
//    s: pathlength of the electron in m
//    n: number density of residual gas H2 in molecule/m^3
// Output parameters:
//    *Eloss: electron energy loss after scattering, in eV
//    *theta: change of polar scattering angle, in degrees
//    *index: this integer shows the scattering type
//          index=0: no scattering
//          index=1: elastic scattering
//          index=2: electronic excitation
//          index=3: ionization
//
// If the electron kinetic energy changes during the motion defined
// by the pathlength, an average kinetic energy should be used.
//
// It is assumed that the scattering probability P=sigma*s*n
// is much  smaller than 1
//   (sigma is the total cross section).
//
  double u, P, sigmel, sigmexc, sigmion, sigma, Pel, Pexc;
  *Eloss = 0.;
  *theta = 0.;
  *index = 0;
  sigmel = sigmael(E);
  sigmexc = sigmaexc(E);
  sigmion = sigmaion(E);
  sigma = sigmel + sigmexc + sigmion;  // total scattering cross section
  P = sigma * s * n;            // scattering probability
  u = random();
  if (u > P) {                   // no scattering in this case
    return;
  }
  
// Start of scattering generation:
  u = random();
  Pel = sigmel / sigma;
  Pexc = sigmexc / sigma;
  if (u < Pel) {                // elastic scattering
    randomel(E, Eloss, theta);
    *index = 1;
  } else if (u >= Pel && u < Pel + Pexc) {  // electronic excitation scattering
    randomexc(E, Eloss, theta);
    *index = 2;
  } else {                         // ionization scattering
    randomion(E, Eloss, theta);
    *index = 3;
  }
}

//////////////////////////////////////////////////
///
///  Elastic scattering:
///
//////////////////////////////////////////////////
double sigmael(double E) 
{
// This function computes the total elastic cross section of
// electron scatt. on molecular hydrogen.
// See: Liu, Phys. Rev. A35 (1987) 591,
//      Trajmar, Phys Reports 97 (1983) 221.
// E: incident electron energy in eV
// sigmael: cross section in m^2
  static double e[14] = { 0., 1.5, 5., 7., 10., 15., 20., 30., 60., 100., 150., 200., 300., 400. };
  static double s[14] = { 9.6, 13., 15., 12., 10., 7., 5.6, 3.3, 1.1, 0.9, 0.5, 0.36, 0.23, 0.15 };
  static double emass = 18780., a02 = 28.e-22;
  double gam, sigma = 0., T;
  int i;
  T = E / 27.2;
  if (E >= 400.) {
    gam = (emass + T) / emass;
    sigma = gam * gam * M_PI / (2. * T) * (4.2106 - 1. / T) * a02;
  } else {
    for (i = 0; i <= 12; i++) {
      if (E >= e[i] && E < e[i + 1]) {
        sigma = 1.e-20 * (s[i] + (s[i + 1] - s[i]) * (E - e[i]) / (e[i + 1] - e[i]));
      }
    }
  }
  return sigma;
}

//////////////////////////////////////////////////////////////////
void randomel(double E, double *Eloss, double *theta) 
{
  
// This subroutine generates  energy loss and polar scatt. angle according to
// electron elastic scattering in molecular hydrogen.
// Input:
//    E: incident electron energy in eV.
// Output:
//   *Eloss: energy loss in eV
//   *theta: change of polar angle in degrees
  static double emass = 18780.; // electron mass in atomic units
  static double clight = 137.;  // velocity of light in atomic units
  static double H2molmass = 69.e6;
  double T, c, b, u[3], G, a, gam, K2, Gmax;
  int i;
  if (E >= 250.) {
    Gmax = 1.e-19;
  } else if (E < 250. && E >= 150.) {
    Gmax = 2.5e-19;
  } else {
    Gmax = 1.e-18;
  }
  T = E / 27.2;
  gam = 1. + T / (clight * clight);      // relativistic correction factor
  b = 2. / (1. + gam) / T;
  for (i = 1; i < 5000; i++) {
    subrn(u, 2);
    c = 1. + b - b * (2. + b) / (b + 2. * u[1]);
    K2 = 2. * T * (1. + gam) * fabs(1. - c);     // momentum transfer squared
    a = (4. + K2) * (4. + K2) / (gam * gam);
    G = a * Del(E, c);
    if (G > Gmax * u[2]) {
      break;
    }
  }
  *theta = acos(c) * 180. / M_PI;
  *Eloss = 2. * emass / H2molmass * (1. - c) * E;
  return;
}


//////////////////////////////////////////////////
///
///  Electronic excitation scattering:
///
//////////////////////////////////////////////////
double sigmaexc(double E) 
{
  
// This function computes the electronic excitation cross section of
// electron scatt. on molecular hydrogen.
// E: incident electron energy in eV,
// sigmaexc: cross section in m^2
  static double a02 = 28.e-22, R = 13.6;
  double sigma;
  if (E < 9.8) {
    sigma = 1.e-40;
  } else if (E >= 9.8 && E <= 250.) {
    sigma = sigmaBC(E) + sigmadiss10(E) + sigmadiss15(E);
  } else {
    sigma = 4. * M_PI * a02 * R / E * (0.80 * log(E / R) + 0.28);
  }
  
//    sigma=sigmainel(E)-sigmaion(E);
      return sigma;
}

////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////
void randomexc(double E, double *Eloss, double *theta) 
{
  
// This subroutine generates  energy loss and polar scatt. angle according to
// electron excitation scattering in molecular hydrogen.
// Input:
//    E: incident electron energy in eV.
// Output:
//   *Eloss: energy loss in eV
//   *theta: change of polar angle in degrees
  static int iff = 0;
  static double sum[1001], fmax, Ecen = 12.6 / 27.21;
  double T, c = 0., u[3], K, xmin, ymin, ymax, x, y, fy, dy, pmax;
  double D, Dmax;
  int i, j, n = 0, N, v = 0;
  
// Energy values of the excited electronic states:
//  (from Mol. Phys. 41 (1980) 1501, in Hartree atomic units)
  static double En[7] = {12.73 / 27.2, 13.2 / 27.2, 14.77 / 27.2, 15.3 / 27.2, 14.93 / 27.2, 15.4 / 27.2, 13.06 / 27.2};
  
// Probability numbers of the electronic states:
//  (from testelectron7.c calculation )
  static double p[7] = {35.86, 40.05, 6.58, 2.26, 9.61, 4.08, 1.54};
  
// Energy values of the B vibrational states:
//   (from: Phys. Rev. A51 (1995) 3745 , in Hartree atomic units)
  static double EB[28] = { 0.411, 0.417, 0.423, 0.428, 0.434, 0.439, 0.444, 0.449, 0.454, 0.459, 0.464, 0.468, 0.473, 0.477, 0.481, 0.485, 0.489, 0.493, 0.496, 0.500, 0.503, 0.507, 0.510, 0.513, 0.516, 0.519, 0.521, 0.524 };
  
// Energy values of the C vibrational states:
//   (from: Phys. Rev. A51 (1995) 3745 , in Hartree atomic units)
  static double EC[14] = { 0.452, 0.462, 0.472, 0.481, 0.490, 0.498, 0.506, 0.513, 0.519, 0.525, 0.530, 0.534, 0.537, 0.539 };
  
// Franck-Condon factors of the B vibrational states:
//   (from: Phys. Rev. A51 (1995) 3745 )
  static double pB[28] = { 4.2e-3, 1.5e-2, 3.0e-2, 4.7e-2, 6.3e-2, 7.3e-2, 7.9e-2, 8.0e-2, 7.8e-2, 7.3e-2, 6.6e-2, 5.8e-2, 5.1e-2, 4.4e-2, 3.7e-2, 3.1e-2, 2.6e-2, 2.2e-2, 1.8e-2, 1.5e-2, 1.3e-2, 1.1e-2, 8.9e-3,
7.4e-3, 6.2e-3, 5.2e-3, 4.3e-3, 3.6e-3 };
  
// Franck-Condon factors of the C vibrational states:
//   (from: Phys. Rev. A51 (1995) 3745 )
  static double pC[14] = { 1.2e-1, 1.9e-1, 1.9e-1, 1.5e-1, 1.1e-1, 7.5e-2, 5.0e-2, 3.3e-2, 2.2e-2, 1.4e-2, 9.3e-3, 6.0e-3, 3.7e-3, 1.8e-3 };
  T = 20000. / 27.2;
  
//
  xmin = Ecen * Ecen / (2. * T);
  ymin = log(xmin);
  ymax = log(8. * T + xmin);
  dy = (ymax - ymin) / 1000.;
  
// Initialization of the sum[] vector, and fmax calculation:
  if (iff == 0) {
    fmax = 0;
    for (i = 0; i <= 1000; i++) {
      y = ymin + dy * i;
      K = exp(y / 2.);
      sum[i] = sumexc(K);
      if (sum[i] > fmax) {
        fmax = sum[i];
      }
    }
    fmax = 1.05 * fmax;
    iff = 1;
  }
  
//
//  Scattering angle *theta generation:
//
  T = E / 27.2;
  if (E >= 100.) {
    xmin = Ecen * Ecen / (2. * T);
    ymin = log(xmin);
    ymax = log(8. * T + xmin);
    dy = (ymax - ymin) / 1000.;
    
// Generation of y values with the Neumann acceptance-rejection method:
    for (j = 1; j < 5000; j++) {
      subrn(u, 2);
      y = ymin + (ymax - ymin) * u[1];
      K = exp(y / 2.);
      fy = sumexc(K);
      if (fmax * u[2] < fy) {
        break;
      }
    }
    
// Calculation of c=cos(theta) and theta:
    x = exp(y);
    c = 1. - (x - xmin) / (4. * T);
    *theta = acos(c) * 180. / M_PI;
  } else if (E <= 25.) {
    Dmax = 60.;
  } else if (E > 25. && E <= 35.) {
    Dmax = 95.;
  } else if (E > 35. && E <= 50.) {
    Dmax = 150.;
  } else {
    Dmax = 400.;
  }
  for (j = 1; j < 5000; j++) {
    subrn(u, 2);
    c = -1. + 2. * u[1];
    D = Dexc(E, c) * 1.e22;
    if (Dmax * u[2] < D) {
      break;
    }
  }
  *theta = acos(c) * 180. / M_PI;
  
// Energy loss *Eloss generation:
      
// First we generate the electronic state, using the Neumann
// acceptance-rejection method for discrete distribution:
  N = 7;                    // the number of electronic states in our calculation
  pmax = p[1];                  // the maximum of the p[] values
  for (j = 1; j < 5000; j++) {
    subrn(u, 2);
    n = (int) (N * u[1]);
    if (u[2] * pmax < p[n]) {
      break;
    }
  }
  if (n < 0) {
    n = 0;
  } 
  if (n > 6) {
    n = 6;
  }
  if (n > 1) {                  // Bp, Bpp, D, Dp, EF states
    *Eloss = En[n] * 27.2;
    return;
  }
  if (n == 0) {                 // B state; we generate now a vibrational state,
                                // using the Frank-Condon factors
    N = 28;                     // the number of B vibrational states in our calculation
    pmax = pB[7];               // maximum of the pB[] values
    for (j = 1; j < 5000; j++) {
      subrn(u, 2);
      v = (int) (N * u[1]);
      if (u[2] * pmax < pB[v]) {
        break;
      }
    }
    if (v < 0) {
      v = 0;
    }
    if (v > 27) {
      v = 27;
    }
    *Eloss = EB[v] * 27.2;
  }
  if (n == 1) {                 // C state; we generate now a vibrational state,
    // using the Franck-Condon factors
    N = 14;                     // the number of C vibrational states in our calculation
    pmax = pC[1];               // maximum of the pC[] values
    for (j = 1; j < 5000; j++) {
      subrn(u, 2);
      v = (int) (N * u[1]);
      if (u[2] * pmax < pC[v]) {
        break;
      }
    }
    if (v < 0) {
      v = 0;
    }
    if (v > 13) {
      v = 13;
    }
    *Eloss = EC[v] * 27.2;
  }
  return;
}


///////////////////////////////////////////////////////////////////
    
//////////////////////////////////////////////////
///
///  Ionization scattering:
///
//////////////////////////////////////////////////
double sigmaion(double E) 
{
  
// This function computes the total ionization cross section of
// electron scatt. on molecular hydrogen.
// E: incident electron energy in eV,
// sigmaion: total ionization cross section of
//   e+H2 --> e+e+H2^+  or  e+e+H^+ +H
// process in m^2.
//
// E<250 eV: Eq. 5 of J. Chem. Phys. 104 (1996) 2956
// E>250: sigma_i formula on page 107 in
//   Phys. Rev. A7 (1973) 103.
// Good agreement with measured results of
// PR A 54 (1996) 2146, and
//   Physica 31 (1965) 94.
//
  static double B = 15.43, U = 15.98, R = 13.6, a02 = 0.28e-20;
  double sigma, t, u, S, r, lnt;
  if (E < 16.) {
    sigma = 1.e-40;
  } else if (E >= 16. && E <= 250.) {
    t = E / B;
    u = U / B;
    r = R / B;
    S = 4. * M_PI * a02 * 2. * r * r;
    lnt = log(t);
    sigma = S / (t + u + 1.) * (lnt / 2. * (1. - 1. / (t * t)) + 1. - 1. / t - lnt / (t + 1.));
  } else {
    sigma = 4. * M_PI * a02 * R / E * (0.82 * log(E / R) + 1.3);
  }
  return sigma;
}


//////////////////////////////////////////////////////////////////
void randomion(double E, double *Eloss, double *theta) 
{
  
// This subroutine generates  energy loss and polar scatt. angle according to
// electron ionization scattering in molecular hydrogen.
// Input:
//    E: incident electron energy in eV.
// Output:
//   *Eloss: energy loss in eV
//   *theta: change of polar angle in degrees
// The kinetic energy of the secondary electron is: Eloss-15.4 eV
//
  static double Ei = 15.45 / 27.21;
  double c, b, u[3], K, xmin, ymin, ymax, x, y, T, G, W, Gmax;
  double q, h, F, Fmin, Fmax, Gp, Elmin, Elmax, qmin, qmax, El, wmax;
  double WcE, Jstarq, WcstarE, w, D2ion;
  int j;
  double K2, KK, fE, kej, ki, kf, Rex, arg, arctg;
  int i;
  double st1, st2;
  
//
// I. Generation of theta
// -----------------------
      Gmax = 1.e-20;
  if (E < 200.) {
    Gmax = 2.e-20;
  }
  T = E / 27.2;
  xmin = Ei * Ei / (2. * T);
  b = xmin / (4. * T);
  ymin = log(xmin);
  ymax = log(8. * T + xmin);
  
// Generation of y values with the Neumann acceptance-rejection method:
  for (j = 1; j < 5000; j++) {
    subrn(u, 2);
    y = ymin + (ymax - ymin) * u[1];
    K = exp(y / 2.);
    c = 1. + b - K * K / (4. * T);
    G = K * K * (Dinel(E, c) - Dexc(E, c));
    if (Gmax * u[2] < G) {
      break;
    }
  }
  
// y --> x --> c --> theta
  x = exp(y);
  c = 1. - (x - xmin) / (4. * T);
  *theta = acos(c) * 180. / M_PI;
  
//
// II. Generation of Eloss, for fixed theta
// ----------------------------------------
//
// For E<=100 eV we use subr. gensecelen
//   (in this case no correlation between theta and Eloss)
  if (E <= 100.) {
    gensecelen(E, &W);
    *Eloss = 15.45 + W;
    return;
  }
  
// For theta>=20 the free electron model is used
//   (with full correlation between theta and Eloss)
  if (*theta >= 20.) {
    *Eloss = E * (1. - c * c);
    return;
  }
  
// For E>100 eV and theta<20: analytical first Born approximation
//   formula of Bethe for H atom (with modification for H2)
//
// Calc. of wmax:
  if (*theta >= 0.7) {
    wmax = 1.1;
  } else if (*theta <= 0.7 && *theta > 0.2) {
    wmax = 2.;
  } else if (*theta <= 0.2 && *theta > 0.05) {
    wmax = 4.;
  } else {
    wmax = 8.;
  }
  
// We generate the q value according to the Jstarq pdf. We have to
// define the qmin and qmax limits for this generation:
  K = sqrt(4. * T * (1. - Ei / (2. * T) - sqrt(1. - Ei / T) * c));
  Elmin = Ei;
  Elmax = (E + 15.45) / 2. / 27.2;
  qmin = Elmin / K - K / 2.;
  qmax = Elmax / K - K / 2.;
  
//
  q = qmax;
  Fmax = 1. / 2. + 1. / M_PI * (q / (1. + q * q) + atan(q));
  q = qmin;
  Fmin = 1. / 2. + 1. / M_PI * (q / (1. + q * q) + atan(q));
  h = Fmax - Fmin;
  
// Generation of Eloss with the Neumann acceptance-rejection method:
  for (j = 1; j < 5000; j++) {
    
// Generation of q with inverse transform method
// (we use the Newton-Raphson method in order to solve the nonlinear eq.
// for the inversion) :
    subrn(u, 2);
    F = Fmin + h * u[1];
    y = 0.;
    for (i = 1; i <= 30; i++) {
      G = 1. / 2. + (y + sin(2. * y) / 2.) / M_PI;
      Gp = (1. + cos(2. * y)) / M_PI;
      y = y - (G - F) / Gp;
      if (fabs(G - F) < 1.e-8) {
        break;
      }
    }
    q = tan(y);
    
// We have the q value, so we can define El, and calculate the weight:
    El = q * K + K * K / 2.;
    
// First Born approximation formula of Bethe for e-H ionization:
    KK = K;
    ki = sqrt(2. * T);
    kf = sqrt(2. * (T - El));
    K2 = 4. * T * (1. - El / (2. * T) - sqrt(1. - El / T) * c);
    if (K2 < 1.e-9) {
      K2 = 1.e-9;
    }
    K = sqrt(K2);              // momentum transfer
    Rex = 1. - K * K / (kf * kf) + K2 * K2 / (kf * kf * kf * kf);
    kej = sqrt(2. * fabs(El - Ei) + 1.e-8);
    st1 = K2 - 2. * El + 2.;
    if (fabs(st1) < 1.e-9) {
      st1 = 1.e-9;
    }
    arg = 2. * kej / st1;
    if (arg >= 0.) {
      arctg = atan(arg);
    }
    
    else {
      arctg = atan(arg) + M_PI;
    }
    st1 = (K + kej) * (K + kej) + 1.;
    st2 = (K - kej) * (K - kej) + 1.;
    fE = 1024. * El * (K2 + 2. / 3. * El) / (st1 * st1 * st1 * st2 * st2 * st2) * exp(-2. / kej * arctg) / (1. - exp(-2. * M_PI / kej));
    D2ion = 2. * kf / ki * Rex / (El * K2) * fE;
    K = KK;
    
//
    WcE = D2ion;
    Jstarq = 16. / (3. * M_PI * (1. + q * q) * (1. + q * q));
    WcstarE = 4. / (K * K * K * K * K) * Jstarq;
    w = WcE / WcstarE;
    if (wmax * u[2] < w) {
      break;
    }
  }
  
//
  *Eloss = El * 27.2;
  
//
  return;
}


///////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////
//
// Additional subroutines
//
//////////////////////////////////////////////////////////////////
double Del(double E, double c) 
{
  
// This subroutine computes the differential cross section
// Del= d sigma/d Omega  of  elastic electron scattering
// on molecular hydrogen.
// See: Nishimura et al., J. Phys. Soc. Jpn. 54 (1985) 1757.
// Input:  E= electron kinetic energy in eV
//         c= cos(theta), where theta is the polar scatt. angle
// Del: in m^2/steradian
  static double a02 = 28.e-22;  // Bohr radius squared
  static double clight = 137.;  // velocity of light in atomic units
  static double Cel[50] = { -0.512, -0.512, -0.509, -0.505, -0.499, -0.491, -0.476, -0.473, -0.462, -0.452, -0.438, -0.422, -0.406, -0.388, -0.370, -0.352, -0.333, -0.314, -0.296, -0.277, -0.258, -0.239, -0.221, -0.202, -0.185, -0.167, -0.151, -0.135, -0.120, -0.105, -0.092, -0.070, -0.053, -0.039, -0.030, -0.024, -0.019, -0.016, -0.014, -0.013, -0.012, -0.009, -0.008, -0.006, -0.005, -0.004, -0.003, -0.002, -0.002, -0.001 };
  static double e[10] = { 0., 3., 6., 12., 20., 32., 55., 85., 150., 250. };
  static double t[10] = { 0., 10., 20., 30.,  40., 60., 80., 100., 140., 180. };
  static double D[9][10] = { {2.9, 2.7,  2.5,  2.1,  1.8,  1.2,    0.9,   1.,     1.6,   1.9}, 
                             {4.2, 3.6,  3.1,  2.5,  1.9,  1.1,    0.8,   0.9,    1.3,   1.4}, 
                             {6.,  4.4,  3.2,  2.3,  1.8,  1.1,    0.7,   0.54,   0.5,   0.6}, 
                             {6.,  4.1,  2.8,  1.9,  1.3,  0.6,    0.3,   0.17,   0.16,  0.23}, 
                             {4.9, 3.2,  2.,   1.2,  0.8,  0.3,    0.15,  0.09,   0.05,  0.05}, 
                             {5.2, 2.5,  1.2,  0.64, 0.36, 0.13,   0.05,  0.03,   0.016, 0.02}, 
                             {4.,  1.7,  0.7,  0.3,  0.16, 0.05,   0.02,  0.013,  0.01,  0.01},
                             {2.8, 1.1,  0.4,  0.15, 0.07, 0.02,   0.01,  0.007,  0.004, 0.003},
                             {1.2, 0.53, 0.2,  0.08, 0.03, 0.0074, 0.003, 0.0016, 0.001, 0.0008}
  };
  double T, K2, K, d, st1, st2, DH, gam, Delreturn = 0., CelK, Ki, theta;
  int i, j;
  T = E / 27.2;
  if (E >= 250.) {
    gam = 1. + T / (clight * clight);     // relativistic correction factor
    K2 = 2. * T * (1. + gam) * (1. - c);
    if (K2 < 0.) {
      K2 = 1.e-30;
    }
    K = sqrt(K2);
    if (K < 1.e-9) {
      K = 1.e-9;                // momentum transfer
    }
    d = 1.4009;                 // distance of protons in H2
    st1 = 8. + K2;
    st2 = 4. + K2;
    
// DH is the diff. cross section for elastic electron scatt.
// on atomic hydrogen within the first Born approximation :
    DH = 4. * st1 * st1 / (st2 * st2 * st2 * st2) * a02;
    
// CelK calculation with linear interpolation.
// CelK is the correction of the elastic electron
// scatt. on molecular hydrogen compared to the independent atom
// model.
    if (K < 3.) {
      i = (int) (K / 0.1);
      Ki = i * 0.1;
      CelK = Cel[i] + (K - Ki) / 0.1 * (Cel[i + 1] - Cel[i]);
    } else if (K >= 3. && K < 5.) {
      i = (int) (30 + (K - 3.) / 0.2);
      Ki = 3. + (i - 30) * 0.2;
      CelK = Cel[i] + (K - Ki) / 0.2 * (Cel[i + 1] - Cel[i]);
    } else if (K >= 5. && K < 9.49) {
      i = (int) (40 + (K - 5.) / 0.5);
      Ki = 5. + (i - 40) * 0.5;
      CelK = Cel[i] + (K - Ki) / 0.5 * (Cel[i + 1] - Cel[i]);
    } else {
      CelK = 0.;
    }
    Delreturn = 2. * gam * gam * DH * (1. + sin(K * d) / (K * d)) * (1. + CelK);
  } else {
    theta = acos(c) * 180. / M_PI;
    for (i = 0; i <= 8; i++) {
      if (E >= e[i] && E < e[i + 1]) {
        for (j = 0; j <= 8; j++) {
          if (theta >= t[j] && theta < t[j + 1]) {
            Delreturn = 1.e-20 * (D[i][j] + (D[i][j + 1] - D[i][j]) * (theta - t[j]) / (t[j + 1] - t[j]));
          }
        }
      }
    }
  }
  return Delreturn;
}


////////////////////////////////////////////////////////////////
double Dexc(double E, double c) 
{
  
// This subroutine computes the differential cross section
// Del= d sigma/d Omega  of excitation electron scattering
// on molecular hydrogen.
// Input:  E= electron kinetic energy in eV
//         c= cos(theta), where theta is the polar scatt. angle
// Dexc: in m^2/steradian
  double K2, K, sigma = 0., T, theta;
  static double a02 = 28.e-22; // Bohr radius squared
  static double EE = 12.6 / 27.2;
  static double e[5] = { 0., 25., 35., 50., 100. };
  static double t[9] = { 0., 10., 20., 30., 40., 60., 80., 100., 180. };
  static double D[4][9] = { {60., 43., 27., 18., 13., 8., 6., 6., 6.}, 
                            {95., 70., 21., 9., 6., 3., 2., 2., 2.,}, 
                            {150., 120., 32., 8., 3.7, 1.9, 1.2, 0.8, 0.8}, 
                            {400., 200., 12., 2., 1.4, 0.7, 0.3, 0.2, 0.2}
                          };
  int i, j;
  
//
  T = E / 27.2;
  if (E >= 100.) {
    K2 = 4. * T * (1. - EE / (2. * T) - sqrt(1. - EE / T) * c);
    if (K2 < 1.e-9) {
      K2 = 1.e-9;
    }
    K = sqrt(K2);              // momentum transfer
    sigma = 2. / K2 * sumexc(K) * a02;
  } else if (E <= 10.) {
    sigma = 0.;
  } else {
    theta = acos(c) * 180. / M_PI;
    for (i = 0; i <= 3; i++) {
      if (E >= e[i] && E < e[i + 1]) {
        for (j = 0; j <= 7; j++) {
          if (theta >= t[j] && theta < t[j + 1]) {
            sigma = 1.e-22 * (D[i][j] + (D[i][j + 1] - D[i][j]) * (theta - t[j]) / (t[j + 1] - t[j]));
          }
        }
      }
    }
  }
  return sigma;
}

double sumexc(double K) 
{
  static double Kvec[15] = { 0., 0.1, 0.2, 0.4, 0.6, 0.8, 1., 1.2, 1.5, 1.8, 2., 2.5, 3., 4., 5. };
  static double fvec[7][15] = { {2.907e-1, 2.845e-1, 2.665e-1, 2.072e-1, 1.389e-1,   // B
                                 8.238e-2, 4.454e-2, 2.269e-2, 7.789e-3, 2.619e-3, 1.273e-3, 2.218e-4, 4.372e-5, 2.889e-6, 4.247e-7},
                                {3.492e-1, 3.367e-1, 3.124e-1, 2.351e-1, 1.507e-1,   // C
                                 8.406e-2, 4.214e-2, 1.966e-2, 5.799e-3, 1.632e-3, 6.929e-4, 8.082e-5, 9.574e-6, 1.526e-7, 7.058e-9},
                                {6.112e-2, 5.945e-2, 5.830e-2, 5.072e-2, 3.821e-2,   // Bp
                                 2.579e-2, 1.567e-2, 8.737e-3, 3.305e-3, 1.191e-3, 6.011e-4, 1.132e-4, 2.362e-5, 1.603e-6, 2.215e-7},
                                {2.066e-2, 2.127e-2, 2.137e-2, 1.928e-2, 1.552e-2,   // Bpp
                                 1.108e-2, 7.058e-3, 4.069e-3, 1.590e-3, 5.900e-4, 3.046e-4, 6.142e-5, 1.369e-5, 9.650e-7, 1.244e-7},
                                {9.405e-2, 9.049e-2, 8.613e-2, 7.301e-2, 5.144e-2,   // D
                                 3.201e-2, 1.775e-2, 8.952e-3, 2.855e-3, 8.429e-4, 3.655e-4, 4.389e-5, 5.252e-6, 9.010e-8, 7.130e-9},
                                {4.273e-2, 3.862e-2, 3.985e-2, 3.362e-2, 2.486e-2,   // Dp
                                 1.612e-2, 9.309e-3, 4.856e-3, 1.602e-3, 4.811e-4, 2.096e-4, 2.498e-5, 2.905e-6, 5.077e-8, 6.583e-9},
                                {0.000e-3, 2.042e-3, 7.439e-3, 2.200e-2, 3.164e-2,      // EF
                                     3.161e-2, 2.486e-2, 1.664e-2, 7.562e-3, 3.044e-3, 1.608e-3, 3.225e-4, 7.120e-5, 6.290e-6, 1.066e-6}
                              };
  static double EeV[7] = { 12.73, 13.20, 14.77, 15.3, 14.93, 15.4, 13.06 };
  int n, j, jmin = 0, nmax;
  double En, f[7], x4[4], f4[4], sum;
  
//
  sum = 0.;
  nmax = 6;
  for (n = 0; n <= nmax; n++) {
    En = EeV[n] / 27.21;        // En is the excitation energy in Hartree atomic units
    if (K >= 5.) {
      f[n] = 0.;
    } else if (K >= 3. && K <= 4.) {
      f[n] = fvec[n][12] + (K - 3.) * (fvec[n][13] - fvec[n][12]);
    } else if (K >= 4. && K <= 5.) {
      f[n] = fvec[n][13] + (K - 4.) * (fvec[n][14] - fvec[n][13]);
    } else {
      for (j = 0; j < 14; j++) {
        if (K >= Kvec[j] && K <= Kvec[j + 1]) {
          jmin = j - 1;
        }
      }
      if (jmin < 0) {
        jmin = 0;
      }
      if (jmin > 11) {
        jmin = 11;
      }
      for (j = 0; j <= 3; j++) {
        x4[j] = Kvec[jmin + j];
        f4[j] = fvec[n][jmin + j];
      }
      f[n] = lagrange(4, x4, f4, K);
    }
    sum += f[n] / En;
  }
  return sum;
}


///////////////////////////////////////////////////////////////////
double sigmaBC(double E) 
{
  
// This function computes the sigmaexc electronic excitation
// cross section to the B and C states, with energy loss
// about 12.5 eV.
// E is incident electron energy in eV,
// sigmaexc in m^2
  static double aB[9] = { -4.2935194e2, 5.1122109e2, -2.8481279e2, 8.8310338e1, -1.6659591e1, 1.9579609, -1.4012824e-1, 5.5911348e-3, -9.5370103e-5 };
  static double aC[9] = { -8.1942684e2, 9.8705099e2, -5.3095543e2, 1.5917023e2, -2.9121036e1, 3.3321027, -2.3305961e-1, 9.1191781e-3, -1.5298950e-4 };
  double lnsigma, lnE, lnEn, sigmaB, Emin, sigma, sigmaC;
  int n;
  sigma = 0.;
  Emin = 12.5;
  lnE = log(E);
  lnEn = 1.;
  lnsigma = 0.;
  if (E < Emin) {
    sigmaB = 0.;
  } else {
    for (n = 0; n <= 8; n++) {
      lnsigma += aB[n] * lnEn;
      lnEn = lnEn * lnE;
    }
    sigmaB = exp(lnsigma);
  }
  sigma += sigmaB;
  
//  sigma=0.;
// C state:
  Emin = 15.8;
  lnE = log(E);
  lnEn = 1.;
  lnsigma = 0.;
  if (E < Emin) {
    sigmaC = 0.;
  } else {
    for (n = 0; n <= 8; n++) {
      lnsigma += aC[n] * lnEn;
      lnEn = lnEn * lnE;
    }
    sigmaC = exp(lnsigma);
  }
  sigma += sigmaC;
  return sigma * 1.e-4;
}


//////////////////////////////////////////////////////////////////
double sigmadiss10(double E) 
{
  
// This function computes the sigmadiss10 electronic
// dissociative excitation
// cross section, with energy loss
// about 10 eV.
// E is incident electron energy in eV,
// sigmadiss10 in m^2
  static double a[9] = { -2.297914361e5, 5.303988579e5, -5.316636672e5, 3.022690779e5, -1.066224144e5, 2.389841369e4, -3.324526406e3, 2.624761592e2, -9.006246604 };
  double lnsigma, lnE, lnEn, Emin, sigma;
  int n;
  
//  E is in eV
  sigma = 0.;
  Emin = 9.8;
  lnE = log(E);
  lnEn = 1.;
  lnsigma = 0.;
  if (E < Emin) {
    sigma = 0.;
  } else {
    for (n = 0; n <= 8; n++) {
      lnsigma += a[n] * lnEn;
      lnEn = lnEn * lnE;
    }
    sigma = exp(lnsigma);
  }
  return sigma * 1.e-4;
}


//////////////////////////////////////////////////////////////////
double sigmadiss15(double E) 
{
  
// This function computes the sigmadiss15 electronic
// dissociative excitation
// cross section, with energy loss
// about 15 eV.
// E is incident electron energy in eV,
// sigmadiss15 in m^2
  static double a[9] = { -1.157041752e3, 1.501936271e3, -8.6119387e2, 2.754926257e2, -5.380465012e1, 6.573972423, -4.912318139e-1, 2.054926773e-2, -3.689035889e-4 };
  double lnsigma, lnE, lnEn, Emin, sigma;
  int n;
  
//  E is in eV
  sigma = 0.;
  Emin = 16.5;
  lnE = log(E);
  lnEn = 1.;
  lnsigma = 0.;
  if (E < Emin) {
    sigma = 0.;
  } else {
    for (n = 0; n <= 8; n++) {
      lnsigma += a[n] * lnEn;
      lnEn = lnEn * lnE;
    }
    sigma = exp(lnsigma);
  }
  return sigma * 1.e-4;
}


///////////////////////////////////////////////////////////////////
double Dinel(double E, double c) 
{
  
// This subroutine computes the differential cross section
// Dinel= d sigma/d Omega  of  inelastic electron scattering
// on molecular hydrogen, within the first Born approximation.
// Input:  E= electron kinetic energy in eV
//         c= cos(theta), where theta is the polar scatt. angle
// Dinel: in m2/steradian
  static double a02 = 28.e-22;  // Bohr radius squared
  static double Cinel[50] = { -0.246, -0.244, -0.239, -0.234, -0.227, -0.219, -0.211, -0.201, -0.190, -0.179, -0.167, -0.155, -0.142, -0.130, -0.118, -0.107, -0.096, -0.085, -0.076, -0.067, -0.059, -0.051, -0.045, -0.039, -0.034, -0.029, -0.025, -0.022, -0.019, -0.016, -0.014, -0.010, -0.008, -0.006, -0.004, -0.003, -0.003, -0.002, -0.002, -0.001, -0.001, -0.001, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000 };
  static double Ei = 0.568;
  double T, K2, K, st1, F, DH, Dinelreturn, CinelK, Ki;
  int i;
  if (E < 16.) {
    return Dexc(E, c);
  }
  T = E / 27.2;
  K2 = 4. * T * (1. - Ei / (2. * T) - sqrt(1. - Ei / T) * c);
  if (K2 < 1.e-9) {
    K2 = 1.e-9;
  }
  K = sqrt(K2);                // momentum transfer
  st1 = 1. + K2 / 4.;
  F = 1. / (st1 * st1);        // scatt. formfactor of hydrogen atom
// DH is the diff. cross section for inelastic electron scatt.
// on atomic hydrogen within the first Born approximation :
  DH = 4. / (K2 * K2) * (1. - F * F) * a02;
  
// CinelK calculation with linear interpolation.
// CinelK is the correction of the inelastic electron
// scatt. on molecular hydrogen compared to the independent atom
// model.
  if (K < 3.) {
    i = (int) (K / 0.1);
    Ki = i * 0.1;
    CinelK = Cinel[i] + (K - Ki) / 0.1 * (Cinel[i + 1] - Cinel[i]);
  } else if (K >= 3. && K < 5.) {
    i = (int) (30 + (K - 3.) / 0.2);
    Ki = 3. + (i - 30) * 0.2;
    CinelK = Cinel[i] + (K - Ki) / 0.2 * (Cinel[i + 1] - Cinel[i]);
  } else if (K >= 5. && K < 9.49) {
    i = (int) (40 + (K - 5.) / 0.5);
    Ki = 5. + (i - 40) * 0.5;
    CinelK = Cinel[i] + (K - Ki) / 0.5 * (Cinel[i + 1] - Cinel[i]);
  } else {
    CinelK = 0.;
  }
  Dinelreturn = 2. * DH * (1. + CinelK);
  return Dinelreturn;
}


////////////////////////////////////////////////////////////////////
void gensecelen(double E, double *W) 
{
  
// This subroutine generates secondary electron energy W
// from ionization of incident electron energy E, by using
// the Lorentzian of Aseev  et al. (Eq. 8).
// E and W in eV.
  static double Ei = 15.45, eps2 = 14.3, b = 6.25;
  static double B;
  double C, A, eps, a, u, epsmax;
  static int iff = 0;
  if (iff == 0) {
    B = atan((Ei - eps2) / b);
    iff = 1;
  }
  epsmax = (E + Ei) / 2.;
  A = atan((epsmax - eps2) / b);
  C = b / (A - B);
  u = random();
  a = b / C * (u + C / b * B);
  eps = eps2 + b * tan(a);
  *W = eps - Ei;
  return;
}


///////////////////////////////////////////////////////////////////
double sigmainel(double E) 
{
  
// This function computes the total inelastic cross section of
// electron scatt. on molecular hydrogen,
//  in the first Born approximation.
// See: Liu, Phys. Rev. A35 (1987) 591.
// E: incident electron energy in eV,
// sigmainel: cross section in m^2
  static double Ei = 0.568;     // ionization energy of molecular
  // hydrogen in Hartree atomic units
  //  (15.45 eV)
  static double a02 = 28.e-22;
  double sigma, gamtot, T;
  T = E / 27.2;
  gamtot = 2. * (-7. / 4. + log(Ei / (2. * T)));
  sigma = 2. * M_PI / T * (1.5487 * log(2. * T) + 2.4036 + gamtot / (2. * T));
  sigma = sigma * a02;
  return sigma;
}


////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////
void subrn(double *u, int len) 
{
  
// This subroutine computes random numbers u[1],...,u[len]
// in the (0,1) interval. It uses the 0<IJKLRANDOM<900000000
// integer as initialization seed.
//  In the calling program the dimension
// of the u[] vector should be larger than len (the u[0] value is
// not used).
// For each IJKLRANDOM
// numbers the program computes completely independent random number
// sequences (see: F. James, Comp. Phys. Comm. 60 (1990) 329, sec. 3.3).
      //
      // remark by T. Thuemmler:
      // same random numbers appear each time one restarts the whole program
      //
  static long IJKLRANDOM = 100;
  static int iff = 0;
  static long ijkl, ij, kl, i, j, k, l, ii, jj, m, i97, j97, ivec;
  static float s, t, uu[98], c, cd, cm, uni;
  if (iff == 0) {
    ijkl = IJKLRANDOM;
    if (ijkl < 1 || ijkl >= 900000000) {
      ijkl = 1;
    }
    ij = ijkl / 30082;
    kl = ijkl - 30082 * ij;
    i = ((ij / 177) % 177) + 2;
    j = (ij % 177) + 2;
    k = ((kl / 169) % 178) + 1;
    l = kl % 169;
    for (ii = 1; ii <= 97; ii++) {
      s = 0;
      t = 0.5;
      for (jj = 1; jj <= 24; jj++) {
        m = (((i * j) % 179) * k) % 179;
        i = j;
        j = k;
        k = m;
        l = (53 * l + 1) % 169;
        if ((l * m) % 64 >= 32) {
          s = s + t;
        }
        t = 0.5 * t;
      }
      uu[ii] = s;
    }
    c = 362436. / 16777216.;
    cd = 7654321. / 16777216.;
    cm = 16777213. / 16777216.;
    i97 = 97;
    j97 = 33;
    iff = 1;
  }
  for (ivec = 1; ivec <= len; ivec++) {
    uni = uu[i97] - uu[j97];
    if (uni < 0.) {
      uni = uni + 1.;
    }
    uu[i97] = uni;
    i97 = i97 - 1;
    if (i97 == 0) {
      i97 = 97;
    }
    j97 = j97 - 1;
    if (j97 == 0) {
      j97 = 97;
    }
    c = c - cd;
    if (c < 0.) {
      c = c + cm;
    }
    uni = uni - c;
    if (uni < 0.) {
      uni = uni + 1.;
    }
    if (uni == 0.) {
      uni = uu[j97] * 0.59604644775391e-07;
      if (uni == 0.) {
        uni = 0.35527136788005e-14;
      }
    }
    u[ivec] = uni;
  }
  
      //  cout << endl<<  "random: " << u[1] << endl << flush;
  return;
}


////////////////////////////////////////////////////////////////
double random() 
{
  
// This function computes 1 random number in the (0,1) interval,
// using the subrn subroutine.
  double u[2];
  subrn(u, 1);
  return u[1];
}


///////////////////////////////////////////////////////////
double lagrange(int n, double *xn, double *fn, double x) 
{
  int i, j;
  double f, a[100], b[100], aa, bb;
  f = 0.;
  for (j = 0; j < n; j++) {
    for (i = 0; i < n; i++) {
      a[i] = x - xn[i];
      b[i] = xn[j] - xn[i];
    }
    a[j] = b[j] = aa = bb = 1.;
    for (i = 0; i < n; i++) {
      aa = aa * a[i];
      bb = bb * b[i];
    }
    f += fn[j] * aa / bb;
  }
  return f;
}


//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
