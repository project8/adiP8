#include "radiation.h"
#include "frequency.h"
#include "paramanage.h"
#include <iostream>
#include <fstream>
#include <cmath>
#include <stdio.h>

// adi2fft performs a calculation of particle position vs. time, specifically in terms of the x,y,z of the guiding center and (epar, eperp, phase) of the orbit.
// This code is meant to translate that into an electric field.
//and calculate the power transmitted in the + or - x-direction

using namespace std;
double eps = 8.85e-31;          //= e0 in C^2/fN/cm^2 

//structures store geometry information
struct t_tl_data tl_data;
double c = Clight * M2CM * 1e-6;	// in cm per us 

double coeff_of_t(double *efield, double *vel, int dir)
{
  //returns A C(t), the coefficient of each waveguide mode in +x-direction
  // in units of fN.cm/C
  double q = 1.60217646e-19;    //C 

  double coeff = q * tl_data.Zw / 2 * (dir * vel[0] * efield[0] + vel[1] * efield[1] + vel[2] * efield[2]);
  return coeff;
}

void init_data()
{
//default data
  tl_data.Zw = 3.77e+25;        //wave imp = Z0 for TEM modes in fN.cm.us/C^2 
  tl_data.R = 0;                //resistance
  tl_data.C = 0;                //capacitance
  tl_data.Zc = tl_data.Zw;      //characteristic imp
  tl_data.att = 0;              //attenuation coefficient
  cout << "Default Characteristic Impedance is wave impedance: " << tl_data.Zc * 1.0e-23 << " Ohms" << endl;
}

void init_tl_data(bool offset)
{
  //Parallel Wire TL
  tl_data.y1 = -0.092;          //cm, y position of 1st wire or strip
  tl_data.y2 = 0.092;           //cm, y position of 2nd wire or strip
  if (offset) {
    tl_data.y1 = -0.2;          //cm, y position of 1st wire
    tl_data.y2 = -0.0;          //cm, y position of 1st wire
    cout << "Offset ";
  }
  tl_data.rI = 0.005;           //wire radius (50 um) in cm, 38 awg
  double l = abs(tl_data.y2 - tl_data.y1) / 2;	//half wire separation in cm
  cout << "Parallel Wire transmission line with separation: " << 2 * l << " cm" << endl;
  double rm = 0.015;            //Ohms, resistivity/skin depth for 80K copper
  tl_data.R = rm / M_PI / tl_data.rI * (l / tl_data.rI / sqrt(pow(l / tl_data.rI, 2) - 1));	//in Ohms/cm 
  cout << " PW Resistance: " << tl_data.R << " Ohms/cm" << endl;
  tl_data.C = eps * M_PI / log(l / tl_data.rI + sqrt(pow(l / tl_data.rI, 2) - 1));	//in C^2/fN/cm^2
  cout << " PW Capacitance : " << tl_data.C * 1e+17 << " Farads/cm" << endl;
  tl_data.Zc = tl_data.Zw * eps / tl_data.C;	//in fN.cm.us/C^2
  cout << " PW Characteristic Impedance: " << tl_data.Zc * 1.0e-23 << " Ohms" << endl;
}

int get_tl_efield(double *p, double *efield)
{
  //function returns the electric field between two infinite parallel wires
  //efield normalized the jackson way with 1/cm units
  //wires extend in x-direction and are can be offset 
  double e_amp = 1 / 2.0 / M_PI * sqrt(tl_data.C / eps);

  //check for hitting wires 
  int status = 0;
  if ((pow((p[1] + tl_data.y1), 2) + pow(p[2], 2) < pow(tl_data.rI, 2)) || (pow((p[1] + tl_data.y2), 2) + pow(p[2], 2) < pow(tl_data.rI, 2))) {
    cout << "Problem!!! Electron hit a wire! " << endl;
    status = 1;
  }
  //calculate effective wire position for efield
  double l = (tl_data.y2 - tl_data.y1) / 2;	//half separation of wires in cm
  double a = sqrt(pow(l, 2) - pow(tl_data.rI, 2));	//effective wire half separation in cm
  double a1 = tl_data.y2 - l - a;	//effective position of first wire in cm
  double a2 = tl_data.y2 - l + a;	//effective position of second wire in cm

  double ra1sq = pow((p[1] + a1), 2) + pow(p[2], 2);
  double ra2sq = pow((p[1] + a2), 2) + pow(p[2], 2);
  efield[0] = 0;                //only true for TE or TEM modes
  efield[1] = e_amp * ((p[1] + a1) / ra1sq - (p[1] + a2) / ra2sq);
  efield[2] = e_amp * (p[2] / ra1sq - p[2] / ra2sq);
  return status;
}

void init_pp_data()
{
//Parallel Plate TL
  tl_data.y1 = -0.092;          //cm, y position of 1st wire or strip
  tl_data.y2 = 0.092;           //cm, y position of 2nd wire or strip
  tl_data.l = 0.37;             //cm, length of strips in z-direction
  tl_data.n = 1.8;              //correction factor for finite plate
  double g = abs(tl_data.y2 - tl_data.y1);	//plate separation in cm
  cout << "Parallel Plate transmission line with separation: " << g << " cm" << endl;
  double rm = 0.015;            //Ohms, resistivity/skin depth for 80K copper
  tl_data.R = rm / tl_data.l;   //in Ohms/cm 
  cout << " PP Resistance: " << tl_data.R << " Ohms/cm" << endl;
  tl_data.C = tl_data.n * eps * tl_data.l / g;	//in C^2/fN/cm^2
  cout << " PP Capacitance : " << tl_data.C * 1e+17 << " Farads/cm" << endl;
  tl_data.Zc = tl_data.Zw * eps / tl_data.C;	//in fN.cm.us/C^2
  cout << " PP Characteristic Impedance: " << tl_data.Zc * 1.0e-23 << " Ohms" << endl;
}

int get_pp_efield(double *p, double *efield)
{
  //function returns the electric field between two infinite parallel plates 
  //efield normalized the jackson way with 1/cm units
  //strips extend in x- and z- directions, so the field is in the y-direction
  double e_amp = sqrt(1 / (tl_data.n * tl_data.l * (tl_data.y2 - tl_data.y1)));
  //factor of sqrt(n=1.8) lower than ideal capacitor 
  int status = 0;

  if ((p[1] < tl_data.y1) || (p[1] > tl_data.y2)) {
    cout << "Problem!!! Electron hit a plate! " << endl;
    status = 1;
  }
  if ((p[2] < -tl_data.l / 2) || (p[2] > tl_data.l / 2)) {
    cout << "Problem!!! Electron outside plates! " << endl;
    status = 1;
  }
  efield[0] = 0;                //only true for TE or TEM modes
  efield[1] = e_amp;
  efield[2] = 0;
  return status;
}

void init_coax_data()
{
  //Coaxial Cable TL
  tl_data.rI = 0.005;           //wire radius (50 um) in cm, 38 awg
  tl_data.rO = 0.2;             //outer radius of coax, in cm
  cout << "Coaxial Cable transmission line with Radius: " << tl_data.rO << " cm" << endl;
  double rm = 0.015;            //Ohms, resistivity/skin depth for 80K copper
  tl_data.R = rm / 2 / M_PI * (1 / tl_data.rI + 1 / tl_data.rO);	//in Ohms/cm 
  cout << " Coaxial Cable Resistance: " << tl_data.R << " Ohms/cm" << endl;
  tl_data.C = eps * 2 * M_PI / log(tl_data.rO / tl_data.rI);	//in C^2/fN/cm^2
  cout << " Coax Capacitance : " << tl_data.C * 1e+17 << " Farads/cm" << endl;
  tl_data.Zc = tl_data.Zw * eps / tl_data.C;	//in fN.cm.us/C^2
  cout << " Coax Characteristic Impedance: " << tl_data.Zc * 1.0e-23 << " Ohms" << endl;
}

int get_coax_efield(double *p, double *efield)
{
  //function returns the electric field inside an infinite coaxial cable 
  //efield normalized the jackson way with 1/cm units
  //cable extend in x-direction 
  double e_amp = 1 / sqrt(2.0 * M_PI * log(tl_data.rO / tl_data.rI));
  int status = 0;

  double r = sqrt(pow(p[1], 2) + pow(p[2], 2));
  if (r < tl_data.rI || r > tl_data.rO) {
    cout << "Problem!!! Electron hit the cable! " << endl;
    status = 1;
  }
  efield[0] = 0;                //only true for TE or TEM modes
  efield[1] = e_amp * p[1] / r / r;
  efield[2] = e_amp * p[2] / r / r;
  return status;
}

void init_sq_wg_data(double k0)
{
//Square Waveguide, centered at origin
  tl_data.y1 = 1.0;             //cm, y-dir, largest dim of wg, >wavelength/2 
  tl_data.y2 = 0.3;             //cm, z-dir, smallest dim of wg, < wavelength/2
  //Warning!  Wave imp. for TE modes freq dep, not implemented properly
  //set impdence for TE modes in fN.cm.us/C^2
  //k0 = omega_cyclotron/c;
  //kj is cutoff frequency for that waveguide
  double z0 = 3.77e+25;         //in fN.cm.us/C^2 
  double kj = M_PI / tl_data.y1;
  tl_data.Zw = 0;
  if (k0 > kj) {
    tl_data.Zw = k0 * z0 / sqrt(pow(k0, 2) - pow(kj, 2));	//in fN.cm.us/C^2
  }
  cout << "Square Waveguide with largest dim: " << tl_data.y1 << " cm" << endl;
  cout << " SqWg Cutoff Freq: " << kj * c / 2 / M_PI << "MHz " << endl;
  cout << " SqWg Wave Impedance: " << tl_data.Zw * 1.0e-23 << " Ohms" << endl;
}

void set_sq_wg_Zw(double k0)
{
  double z0 = 3.77e+25;         //in fN.cm.us/C^2 
  double kj = M_PI / tl_data.y1;
  tl_data.Zw = 0;
  if (k0 > kj) {
    tl_data.Zw = k0 * z0 / sqrt(pow(k0, 2) - pow(kj, 2));	//in fN.cm.us/C^2
  }
}

int get_sq_wg_efield(double *pos, double *efield)
{
  //function returns the TE_10 efield inside an infinite square waveguide 
  //centered at y=0, z=0 (not with corner at origin as is standard)
  //efield normalized the jackson way with 1/cm units
  //waveguide extends in x-direction
  // wavelength of 27 GHz radiation is 1.1 cm
  double e_amp = sqrt(2 / tl_data.y1 / tl_data.y2);
  int status = 0;

  if ((pos[1] > tl_data.y1 / 2) || (pos[1] < -tl_data.y1 / 2)) {
    cout << "Problem!!! Electron hit a wall in y-dir! " << endl;
    status = 1;
  }
  if ((pos[2] > tl_data.y2 / 2) || (pos[2] < -tl_data.y2 / 2)) {
    cout << "Problem!!! Electron hit a wall in z-dir! " << endl;
    status = 1;
  }
  efield[0] = 0;                //is true for TE mode
  efield[1] = 0;
  efield[2] = e_amp * sin(M_PI * (pos[1] + tl_data.y1 / 2) / tl_data.y1);
  return status;
}

void init_circ_wg_data(double k0)
{
  //Circular Waveguide
  tl_data.rO = 1.0;             //cm, radius of wg, >wavelength/3.41 
  //Warning!  Wave imp. for TE modes freq dep, not implemented properly
  //set impdence for TE modes in fN.cm.us/C^2
  //k0 = omega_s/c;
  //kj is cutoff frequency for that waveguide
  double p1 = 1.841;            //1st zero of the derivate of bessel function
  double z0 = 3.77e+25;         //in fN.cm.us/C^2 
  double kj = p1 / tl_data.rO;
  tl_data.Zw = 0;
  if (k0 > kj) {
    tl_data.Zw = k0 * z0 / sqrt(pow(k0, 2) - pow(kj, 2));	//in fN.cm.us/C^2
  }
  cout << "Circ Waveguide with Radius: " << tl_data.rO << " cm" << endl;
  cout << " Circ WG Cutoff Freq: " << kj * c / 2 / M_PI << "MHz " << endl;
  cout << " Circ WG Wave Impedance: " << tl_data.Zw * 1.0e-23 << " Ohms" << endl;
}

void set_circ_wg_Zw(double k0)
{
  //Circular Waveguide
  double p1 = 1.841;            //1st zero of the derivate of bessel function
  double kj = p1 / tl_data.rO;
  double z0 = 3.77e+25;         //in fN.cm.us/C^2 
  tl_data.Zw = 0;
  if (k0 > kj) {
    tl_data.Zw = k0 * z0 / sqrt(pow(k0, 2) - pow(kj, 2));	//in fN.cm.us/C^2
  }
}

int get_circ_wg_efield(double phase, double *pos, double *efield)
{
  //returns the TE_11 efield inside an infinite circ. waveguide 
  //efield normalized the jackson way with 1/cm units
  //waveguide extends in x-direction
  // wavelength of 27 GHz radiation is 1.1 cm
  double p1 = 1.841;            //1st zero of the derivate of bessel function
  double radius = sqrt(pos[1] * pos[1] + pos[2] * pos[2]);
  //bessel functions expanded for small arguement
  //these are constant if the orbit is centered around wg axis
  double J1 = p1 * radius / tl_data.rO / 2 - pow(p1 * radius / tl_data.rO / 2, 3) / 2;	//this term cancels
  double Jp = p1 / tl_data.rO / 2 - p1 / tl_data.rO * 3 * pow(p1 * radius / tl_data.rO, 2) / 16;
  double e_amp = 5.574;
  int status = 0;

  if (radius > tl_data.rO) {
    cout << "Problem!!! Electron hit a wall! " << endl;
    status = 1;
  }
  efield[0] = 0;                //only true for TE mode
  efield[1] = e_amp * (J1 / radius + p1 / tl_data.rO * Jp) * sin(phase) * cos(phase);
  efield[2] = e_amp * (J1 / radius * sin(phase) * sin(phase)
                       - p1 / tl_data.rO * Jp * cos(phase) * cos(phase));

  return status;
}

int get_circ_cavity_efield(double phase, double *pos, double *efield)
{
  //DO NOT USE YET!  UNFINISHED AND UNTESTED
  //function returns the TE_111 efield inside an circular cavity 
  //efield normalized the jackson way with 1/cm units
  //cavity extends in x-direction
  // wavelength of 27 GHz radiation is 1.1 cm
  double a = 2.3;               //cm, radius of cavity, > wavelength/3.41
  double d = 6.97;              //cm, length of cavity, tuneable, d>2.03a
  double p1 = 1.841;            //1st zero of the derivate of bessel function
  double e_amp = 5.574;
  double radius = sqrt(pos[1] * pos[1] + pos[2] * pos[2]);
  //bessel functions expanded for small arguement
  //these are constant if the orbit is centered around cavity axis
  double J1 = p1 * radius / a / 2 - pow(p1 * radius / a / 2, 3) / 2;	//this term cancels
  double Jp = p1 / a / 2 - p1 / a * 3 * pow(p1 * radius / a, 2) / 16;
  int status = 0;

  if ((radius < a) || (pos[0] < -d / 2) || (pos[0] > d / 2)) {
    cout << "Problem!!! Electron hit a wall! " << endl;
    status = 1;
  }
  efield[0] = 0;                //only true for TE mode
  efield[1] = e_amp * (J1 / radius + p1 / a * Jp) * sin(phase) * cos(phase);
  efield[2] = e_amp * (J1 / radius * sin(phase) * sin(phase)
                       - p1 / a * Jp * cos(phase) * cos(phase));

  return status;
}

double antenna_at_infinity(double phase, double x, double y, double z, double eperp, double epar, double b)
{
// The electric field strength at a given distance should be a bunch of angle/position/distance-dependent prefactors (the power) times cos(theta) (the phase factor) and maybe also times the antenna area.   The subroutine antenna_at_infinity returns that prefactor for an antenna located at (0,0,zdet).  
//no longer maintained (-MLL)
  const double zdet = 0.5;
  const double me2 = 510998.0;  //electron mass in eV
  double pitchangle = acos(eperp / (sqrt(eperp * eperp + epar * epar)));
  double d2_to_det = x * x + y * y + pow(z - zdet, 2);
  double cangle_to_det = x / sqrt(d2_to_det);
  double gamma2 = (eperp * eperp + epar * epar + me2) / me2;
  return dpdd2(b, sqrt(1 - 1 / gamma2), pitchangle, cangle_to_det) / d2_to_det;
}
