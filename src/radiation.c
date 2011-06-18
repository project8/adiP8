#include "radiation.h"
#include "paramanage.h"
#include <iostream>
#include <fstream>
#include <cmath>
#include <stdio.h>

// adi2fft performs a calculation of particle position vs. time, 
//specifically in terms of the x,y,z of the guiding center plus
//cyclotron radius and phase of the orbit.
// This code is meant to translate that into an electric field.
//and calculate the power transmitted in the + or - x-direction
//attempting to use SI units

using namespace std;

//structures store geometry information
struct t_tl_data tl_data;

double coeff_of_t(double *efield, double *vel, int dir)
{
  //returns A C(t), the coefficient of each waveguide mode in +x-direction
  // in units of C*Ohm/s or volts
  return Echarge * tl_data.Zw / 2 * (dir * vel[0] * efield[0] + vel[1] * efield[1] + vel[2] * efield[2])/US2S;
}

void init_data()
{
//default data
  tl_data.Zw = Z0;        //wave imp = Z0 for TEM modes in Ohms 
  tl_data.R = 0;                //resistance
  tl_data.C = 0;                //capacitance
  tl_data.Zc = tl_data.Zw;      //characteristic imp
  tl_data.vg = Clight * M2CM / S2US;//group velocity is speed of light
  tl_data.att = 0;              //attenuation coefficient
  cout << "Default Wave Impedance : " << tl_data.Zw << " Ohms" << endl;
  cout << "Default Characteristic Impedance : " << tl_data.Zc << " Ohms" << endl;
  cout << "Default Group velocity is speed of light: " << tl_data.vg << " cm/us " << endl;
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
  tl_data.R = RM / M_PI / tl_data.rI * (l / tl_data.rI / sqrt(pow(l / tl_data.rI, 2) - 1));//in Ohms/cm 
  tl_data.C = EPS0 / M2CM * M_PI / log(l / tl_data.rI + sqrt(pow(l / tl_data.rI, 2) - 1));//in F/cm
  calculate_tl_parameters();
}

int get_tl_efield(double *p, double *efield)
{
  //function returns the electric field between two infinite parallel wires
  //efield normalized the jackson way with 1/cm units
  //wires extend in x-direction and are can be offset 
  double e_amp = 1 / 2.0 / M_PI * sqrt(tl_data.C / (EPS0 / M2CM));

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
  tl_data.R = RM / tl_data.l;   //in Ohms/cm 
  tl_data.C = tl_data.n * EPS0 / M2CM * tl_data.l / g;	//in F/cm
  calculate_tl_parameters();
}

int get_pp_efield(double *p, double *efield)
{
  //function returns the electric field between two parallel plates 
  //amplitude corrected for capacitance of finite width strips
  //not valid outside plates
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
  tl_data.R = RM / 2 / M_PI * (1 / tl_data.rI + 1 / tl_data.rO);	//in Ohms/cm 
  tl_data.C = EPS0 / M2CM * 2 * M_PI / log(tl_data.rO / tl_data.rI);	//in F/cm
  calculate_tl_parameters();
}

void calculate_tl_parameters() {
  cout << " Resistance: " << tl_data.R << " Ohms/cm" << endl;
  cout << " Capacitance : " << tl_data.C << " Farads/cm" << endl;
  tl_data.Zc = tl_data.Zw * EPS0 / M2CM / tl_data.C;	//in Ohms
  cout << " Characteristic Impedance: " << tl_data.Zc << " Ohms" << endl;
  
  double tand = 0.0012;         //Rogers Duroid tangent delta dimensionless
  double attenG = tl_data.C * tl_data.Zc * OMEGA0 * tand / 2;
  double attenR = tl_data.R / tl_data.Zc / 2;
  tl_data.att = attenR + attenG;
  cout << " Resistive atten coeff: " << attenR << " Neper per cm " << endl;
  cout << " Conductive atten coeff: " << attenG << " Neper per cm " << endl;
  
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
  //Square Waveguide, WR-42, centered at origin
  tl_data.y1 = 1.0668;             //cm, y-dir, largest dim of wg, >wavelength/2 
  tl_data.z1 = 0.4318;             //cm, z-dir, smallest dim of wg, < wavelength/2
  //Warning!  Wave imp. for TE modes freq dep, not implemented properly
  //set impdence for TE modes in Ohm
  //k0 = omega_cyclotron/c;
  //kj is cutoff frequency for that waveguide
  double kj = M_PI / tl_data.y1;
  double betaj = 0;
  tl_data.Zw = 0;
  if (k0 > kj) {
    betaj = sqrt(pow(k0, 2) - pow(kj, 2));
    tl_data.Zw = k0 * Z0 / betaj;	//in Ohm
    tl_data.vg = Clight * M2CM / S2US * betaj / k0;	//group velocity, cm/us
    tl_data.att = RM/tl_data.y1/tl_data.z1/betaj/k0/Z0*(2*tl_data.z1*kj*kj+tl_data.y1*k0*k0);
  }
  cout << "Square Waveguide with largest dim: " << tl_data.y1 << " cm" << endl;
  cout << " SqWg Cutoff Freq: " << kj * Clight * M2CM * 1e-9 / 2 / M_PI << " GHz " << endl;
  cout << " SqWg Wave Impedance: " << tl_data.Zw << " Ohms" << endl;
  cout << " SqWg Attenuation: " << tl_data.att << " Nepers per cm" << endl;
  cout << " SqWg Group Velocity: " << tl_data.vg << " cm/us" << endl;

}

int get_sq_wg_efield(double *pos, double *efield)
{
  //function returns the TE_10 efield inside an infinite square waveguide 
  //centered at y=0, z=0 (not with corner at origin as is standard)
  //efield normalized the jackson way with 1/cm units
  //waveguide extends in x-direction
  // wavelength of 27 GHz radiation is 1.1 cm
  double e_amp = sqrt(2 / tl_data.y1 / tl_data.z1);
  int status = 0;

  if ((pos[1] > tl_data.y1 / 2) || (pos[1] < -tl_data.y1 / 2)) {
    cout << "Problem!!! Electron hit a wall in y-dir! " << endl;
    status = 1;
  }
  if ((pos[2] > tl_data.z1 / 2) || (pos[2] < -tl_data.z1 / 2)) {
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
  //set impdence for TE modes in Ohm
  //k0 = omega_s/c;
  //kj is cutoff frequency for that waveguide
  double p1 = 1.841;            //1st zero of the derivate of bessel function
  double kj = p1 / tl_data.rO;
  tl_data.Zw = 0;
  if (k0 > kj) {
    tl_data.Zw = k0 * Z0 / sqrt(pow(k0, 2) - pow(kj, 2));	//in Ohms
  }
  cout << "Circ Waveguide with Radius: " << tl_data.rO << " cm" << endl;
  cout << " Circ WG Cutoff Freq: " << kj * Clight * M2CM * 1e-9 / 2 / M_PI << " GHz " << endl;
  cout << " Circ WG Wave Impedance: " << tl_data.Zw << " Ohms" << endl;
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

