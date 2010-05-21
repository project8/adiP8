#include "radiation.h"
#include "frequency.h"
#include <iostream>
#include <fstream>
#include <cmath>
#include <stdio.h>

using namespace std;
// Adi2fft_v2 performs a calculation of particle position vs. time, specifically in terms of the x,y,z of the guiding center and (epar, eperp, phi) of the orbit.
// This code is meant to translate that into an electric field.


/* The electric field strength at a given distance should be a bunch of angle/position/distance-dependent prefactors (the power) times cos(theta) (the phase factor) and maybe also times the antenna area.   The subroutine antenna_at_infinity returns that prefactor for an antenna located at (0,0,zdet).  */
double antenna_at_infinity(double phi,double x,double y,double z,double eperp,double epar,double b)
{
  const double zdet = 0.5;
  const double me2 = 510998.0; //electron mass in eV
  double pitchangle = acos(eperp/(sqrt(eperp*eperp + epar*epar)));
  double d2_to_det = x*x + y*y + pow(z-zdet,2);
  double cangle_to_det = x/sqrt(d2_to_det);
  double gamma2 = (eperp*eperp+epar*epar+me2)/me2;
  return dpdd2(b,sqrt(1-1/gamma2),pitchangle,cangle_to_det)/d2_to_det;
}
