#include "radiation.h"
#include "paramanage.h"
#include <iostream>
#include <fstream>
#include <cmath>
#include <stdio.h>
#include "TMath.h"
#include "TVector3.h"
// adi2fft performs a calculation of particle position vs. time, 
//specifically in terms of the x,y,z of the guiding center plus
//cyclotron radius and phase of the orbit.
// This code is meant to translate that into an electric field.
//and calculate the power transmitted in the + or - z-direction
//(note different coordinates from adipark!!)

using namespace std;

//structures store geometry information
struct t_tl_data tl_data;

double coeff_of_t(TVector3 &efield, TVector3 &vel, int dir)
{
  //returns A C(t), the coefficient of each waveguide mode in +x-direction
  // in units of C*Ohm/s or volts
  efield.SetZ( dir * efield.Z() );//only effect TM modes
  return Echarge * tl_data.Zw / 2 * ( efield * vel ) / US2S;
}

void init_data()
{
//default data
  tl_data.Zw = Z0;        //wave imp = Z0 for TEM modes in Ohms 
  tl_data.R = 0;                //resistance
  tl_data.C = 0;                //capacitance
  tl_data.Zc = tl_data.Zw;      //characteristic imp
  tl_data.vg = Clight * M2CM / S2US;//group velocity is speed of light
  tl_data.vp = Clight * M2CM / S2US;//phase velocity is speed of light
  tl_data.att = 0;              //attenuation coefficient
  cout << "Default Wave Impedance : " << tl_data.Zw << " Ohms" << endl;
  cout << "Default Characteristic Impedance : " << tl_data.Zc << " Ohms" << endl;
  cout << "Default Group and Phase velocities are the speed of light: ";
  cout  << tl_data.vg << " cm/us " << endl;
}

void init_tl_data(bool offset)
{
  //Parallel Wire TL
  tl_data.x1 = -0.092;          //cm, x position of 1st wire or strip
  tl_data.x2 = 0.092;           //cm, x position of 2nd wire or strip
  if (offset) {
    tl_data.x1 = -0.2;          //cm, x position of 1st wire
    tl_data.x2 = -0.0;          //cm, x position of 1st wire
    cout << "Offset ";
  }
  tl_data.rI = 0.005;           //wire radius (50 um) in cm, 38 awg
  double l = abs(tl_data.x2 - tl_data.x1) / 2;	//half wire separation in cm
  cout << endl << endl << "Parallel Copper Wire transmission line with separation: ";
  cout << 2 * l << " cm" << endl;
  tl_data.skinD = sqrt(2*CU_R/OMEGA0/MU0);//m
  double Rm = CU_R/tl_data.skinD;//Ohms, for copper at 27 GHz and 80 K
  tl_data.R = Rm / M_PI / tl_data.rI * (l / tl_data.rI / sqrt(pow(l / tl_data.rI, 2) - 1.));//in Ohms/cm 
  tl_data.C = EPS0 / M2CM * M_PI / log(l / tl_data.rI + sqrt(pow(l / tl_data.rI, 2) - 1.));//in F/cm
  calculate_tl_parameters();
  double R = 0.045;//cm, cyclotron radius
  cout << "Approximate cutoff freq of 1st mode: " << Clight * M2CM / 4 / l * 1e-9 << " GHz " << endl;
  print_tl_power(R);
}


void print_tl_power(double R)
{
  //returns the time-average power radiated into 
  //an infinite parallel wire transmission line
  // from an electron at the center of the line with:
  double l = abs(tl_data.x2 - tl_data.x1) / 2;	//half wire separation in cm
  double power = tl_data.C*M2CM/2/Clight * pow(Echarge*R*OMEGA0/(M_PI*EPS0*2*l),2);
  cout << "Power at center of TL: " << 1.e15 * power << " fW " << endl;
}

int get_tl_efield(TVector3 &p,TVector3 &efield)
{
  //function returns the electric field between two infinite parallel wires
  //efield normalized the jackson way with 1/cm units
  //wires extend in z-direction and are can be offset 
  double e_amp = 1. / 2.0 / M_PI * sqrt(tl_data.C / (EPS0 / M2CM));

  //calculate effective wire position for efield
  TVector3 x1(tl_data.x1, 0, 0);	//real position of first wire in cm
  TVector3 x2(tl_data.x2, 0, 0);	//real position of second wire in cm
  TVector3 l = x1 - x2;
  l *= 1./2;
  double a = sqrt(l.Mag2() - pow(tl_data.rI, 2));	//effective wire half separation in cm
  TVector3 a1 = x1 - l + a*l.Unit();	//effective position of first wire in cm
  TVector3 a2 = x1 - l - a*l.Unit();	//effective position of second wire in cm
 
  //vector from point to wires
  p.SetZ(0);
  TVector3 r1 = p + a1;
  TVector3 r2 = p + a2;
   
  //check for hitting wires 
  int status = 0;
  if ( ( (x1-r1).Mag() < tl_data.rI ) || ( (x2-r2).Mag() < tl_data.rI ) ) {
    cout << "Problem!!! Electron hit a wire! Radius " << p.Mag() << endl;
    status = 1;
  }
  //efield = e_amp * (r1 * 1/r1.Mag2() - r2 * 1/r2.Mag2());
  efield.SetX( e_amp * ( r1.X() / r1.Mag2() - r2.X() / r2.Mag2() ) );
  efield.SetY( e_amp * ( r1.Y() / r1.Mag2() - r2.Y() / r2.Mag2() ) );
  efield.SetZ( 0 );                //only true for TE or TEM modes
  return status;
}

void init_pp_data()
{
//Parallel Plate TL
  tl_data.x1 = -0.092;          //cm, x position of 1st wire or strip
  tl_data.x2 = 0.092;           //cm, x position of 2nd wire or strip
  tl_data.l = 0.37;             //cm, length of strips in y-direction
  tl_data.n = 1.8;              //correction factor for finite plate
  double g = abs(tl_data.x2 - tl_data.x1);	//plate separation in cm
  cout << endl << endl << "Parallel Copper Plate transmission line with separation: ";
  cout  << g << " cm" << endl;
  tl_data.skinD = sqrt(2.*CU_R/OMEGA0/MU0);//m
  double Rm = CU_R/tl_data.skinD;//Ohms, for copper at 27 GHz and 80 K
  tl_data.R = Rm / tl_data.l;   //in Ohms/cm 
  tl_data.C = tl_data.n * EPS0 / M2CM * tl_data.l / g;	//in F/cm
  calculate_tl_parameters();
}

int get_pp_efield(TVector3 &p, TVector3 &efield)
{
  //function returns the electric field between two parallel plates 
  //amplitude corrected for capacitance of finite width strips
  //not valid outside plates
  //efield normalized the jackson way with 1/cm units
  //strips extend in y- and z- directions, so the field is in the x-direction
  double e_amp = sqrt(1 / (tl_data.n * tl_data.l * (tl_data.x2 - tl_data.x1)));
  //factor of sqrt(n=1.8) lower than ideal capacitor 
  int status = 0;

  if ((p.X() < tl_data.x1) || (p.X() > tl_data.x2)) {
    cout << "Problem!!! Electron hit a plate! " << endl;
    status = 1;
  }
  if ((p.Y() < -tl_data.l / 2) || (p.Y() > tl_data.l / 2)) {
    cout << "Problem!!! Electron outside plates! " << endl;
    status = 1;
  }
  efield.SetX( e_amp );
  efield.SetY( 0 );
  efield.SetZ( 0 );                //only true for TE or TEM modes
  return status;
}

void init_coax_data()
{
  //Coaxial Cable TL
  tl_data.rI = 0.005;           //wire radius (50 um) in cm, 38 awg
  tl_data.rO = 0.2;             //outer radius of coax, in cm
  cout << endl << endl << "Copper Coaxial Cable transmission line with Radius: ";
  cout << tl_data.rO << " cm" << endl;
  tl_data.skinD = sqrt(2*CU_R/OMEGA0/MU0);//m
  double Rm = CU_R/tl_data.skinD;//Ohms, for copper at 27 GHz and 80 K
  tl_data.R = Rm / 2. / M_PI * (1. / tl_data.rI + 1. / tl_data.rO);	//in Ohms/cm 
  tl_data.C = EPS0 / M2CM * 2. * M_PI / log(tl_data.rO / tl_data.rI);	//in F/cm
  calculate_tl_parameters();
}

void calculate_tl_parameters() {
  cout << " Skin Depth " << tl_data.skinD*1e6 << " micron " << endl;
  cout << " Resistance: " << tl_data.R << " Ohms/cm" << endl;
  cout << " Capacitance : " << tl_data.C << " Farads/cm" << endl;
  tl_data.Zc = tl_data.Zw * EPS0 / M2CM / tl_data.C;	//in Ohms
  cout << " Characteristic Impedance: " << tl_data.Zc << " Ohms" << endl;
  
  double tand = 0.0012;         //Rogers Duroid tangent delta dimensionless
  double attenG = tl_data.C * tl_data.Zc * OMEGA0 * tand / 2.;
  double attenR = tl_data.R / tl_data.Zc / 2.;
  tl_data.att = attenR + attenG;
  cout << " Resistive atten coeff: " << attenR << " Neper per cm " << endl;
  cout << " Conductive atten coeff: " << attenG << " Neper per cm " << endl;
  
}

int get_coax_efield(TVector3 &p, TVector3 &efield)
{
  //function returns the electric field inside an infinite coaxial cable 
  //efield normalized the jackson way with 1/cm units
  //cable extend in z-direction 
  double e_amp = 1 / sqrt(2.0 * M_PI * log(tl_data.rO / tl_data.rI));
  int status = 0;

  p.SetZ(0);
  if (p.Mag() < tl_data.rI || p.Mag() > tl_data.rO) {
    cout << "Problem!!! Electron hit the cable! " << endl;
    status = 1;
  }
  efield.SetX( e_amp * p.X() / p.Mag2() );
  efield.SetY( e_amp * p.Y() / p.Mag2() );
  efield.SetZ( 0 );                //only true for TE or TEM modes
  return status;
}

void init_sq_wg_data(double k0)
{
  //Aluminum Square Waveguide, WR-42, centered at origin
  //Ensures only TE10 mode propagates
  tl_data.x1 = 1.0668;             //cm, x-dir, largest dim of wg, >wavelength/2 
  tl_data.y1 = 0.4318;             //cm, y-dir, smallest dim of wg, < wavelength/2
  //WR-34
  //tl_data.x1 = 0.8636;            //cm, x-dir, largest dim of wg, >wavelength/2 
  //tl_data.y1 = 0.4318;            //cm, y-dir, smallest dim of wg, < wavelength/2
  tl_data.skinD = sqrt(2.*AL_R/OMEGA0/MU0);//m
  double Rm = AL_R/tl_data.skinD;//Ohms, for Aluminum at 27 GHz and 80 K

  //set impdence for TE modes in Ohm
  //k0 = omega_cyclotron/c;//angular wavenumber of free-space cyclotron radiatin
  //k10 is angular wavenumber of cutoff frequency for TE10 mode
  double k10 = M_PI / tl_data.x1;
  double k11 = M_PI * sqrt(1 / pow(tl_data.x1,2) + 1/ pow(tl_data.y1,2));
  double k20 = M_PI * 2. /tl_data.x1;
  double k01 = M_PI /tl_data.y1;
  double betaj = 0;//angular wavenumber of propagating radiation
  tl_data.Zw = 0;
  if (k0 > k10) {
    betaj = sqrt(pow(k0, 2) - pow(k10, 2));
    tl_data.Zw = k0 * Z0 / betaj;	//in Ohm
    tl_data.vg = Clight * M2CM / S2US * betaj / k0;	//group velocity, cm/us
    tl_data.vp = Clight * M2CM / S2US * k0 / betaj;	//phase velocity, cm/us
    tl_data.att = Rm/tl_data.x1/tl_data.y1/betaj/k0/Z0*(2.*tl_data.y1*k10*k10+tl_data.x1*k0*k0);
  }
  cout << endl << endl << "Aluminum Square Waveguide with largest dim: ";
  cout << tl_data.x1 << " cm" << endl;
  cout << " SqWg TE10 Cutoff (Min) Freq: " << k10 / 2. / M_PI * Clight * M2CM * 1e-9 << " GHz " << endl;
  cout << " SqWg TE20 Cutoff (Min) Freq: " << k20 / 2. / M_PI * Clight * M2CM * 1e-9 << " GHz " << endl;
  cout << " SqWg TE01 Cutoff (Min) Freq: " << k01 / 2. / M_PI * Clight * M2CM * 1e-9 << " GHz " << endl;
  cout << " SqWg TE11 or TM11 Cutoff (Min) Freq: " << k11 / 2. / M_PI * Clight * M2CM * 1e-9 << " GHz " << endl;
  cout << " SqWg Wave Impedance: " << tl_data.Zw << " Ohms" << endl;
  cout << " Skin Depth " << tl_data.skinD*1e6 << " micron " << endl;
  cout << " SqWg Attenuation: " << tl_data.att << " Nepers per cm" << endl;
  cout << " SqWg Group Velocity: " << tl_data.vg << " cm/us" << endl;
  cout << " SqWg Phase Velocity: " << tl_data.vp << " cm/us" << endl;
  print_sq_wg_power(k0);
}

void print_sq_wg_power(double k0)
{
  //returns the time-average power radiated into the TE_10 
  //mode of an infinite square waveguide 
  // from an electron at the center of the guide with:
  double R = 0.045;//cm, cyclotron radius
  TVector3 pos(0,0,0);
  //k0 = omega_cyclotron/c;//angular wavenumber of free-space cyclotron radiatin
  //k10 is angular wavenumber of cutoff frequency for TE10 mode
  double k10 = M_PI / tl_data.x1;
  double betaj = sqrt(pow(k0, 2) - pow(k10, 2));//angular wavenumber of propagating radiation
  double c = Clight * M2CM ;//cm/s
  double x0 = pos.X() + tl_data.x1;
  double power = Z0 * pow(Echarge*R*c*sin(M_PI*x0/2),2) * pow(k0,3)/(4 * betaj * tl_data.y1 * tl_data.x1);
  cout << "Power at center of WG: " << 1e15 * power << " fW " << endl;
}

int get_sq_wg_efield(TVector3 &pos, TVector3 &efield)
{
  //function returns the TE_10 efield inside an infinite square waveguide 
  //centered at x=0, y=0 (not with corner at origin as is standard)
  //efield normalized the jackson way with 1/cm units
  //waveguide extends in z-direction
  // wavelength of 27 GHz radiation is 1.1 cm
  double e_amp = sqrt(2 / tl_data.x1 / tl_data.y1);
  int status = 0;

  if ((pos.X() > tl_data.x1 / 2.) || (pos.X() < -tl_data.x1 / 2.)) {
    cout << "Problem!!! Electron hit a wall in x-dir! " << endl;
    status = 1;
  }
  if ((pos.Y() > tl_data.y1 / 2.) || (pos.Y() < -tl_data.y1 / 2.)) {
    cout << "Problem!!! Electron hit a wall in y-dir! " << endl;
    status = 1;
  }
  efield.SetX( 0 );
  efield.SetY( e_amp * sin(M_PI * (pos.X() + tl_data.x1 / 2.) / tl_data.x1) );
  efield.SetZ( 0 );                //is true for TE mode
  return status;
}

void init_circ_wg_data(double k0)
{
  //Circular Waveguide with TE_11 mode propagating
  tl_data.rO = 0.38;             //cm, radius of wg, >wavelength/3.41 
  tl_data.skinD = sqrt(2.*AL_R/OMEGA0/MU0);//m
  double Rm = AL_R/tl_data.skinD;//Ohms, for Aluminum at 27 GHz and 80 K
  //Warning!  Wave imp. for TE modes freq dep, not implemented properly
  //set impdence for TE modes in Ohm
  //k0 = omega_s/c;//angular wavenumber of free-space cyclotron radiation
  //k11 is angular wavenumber for cutoff frequency for TE11
  double p11 = 1.841;            //1st zero of the derivate of bessel function
  double k11 = p11 / tl_data.rO;
  double p01 = 2.405;            //1st zero of the bessel function
  double k01 = p01 / tl_data.rO;
  double betaj = 0;//angular wavenumber of propagating radiation
  tl_data.Zw = 0;
  if (k0 > k11) {
    betaj = sqrt(pow(k0, 2) - pow(k11, 2));
    tl_data.Zw = k0 * Z0 / betaj;	//in Ohms
    tl_data.vg = Clight * M2CM / S2US * betaj / k0;	//group velocity, cm/us
    tl_data.vp = Clight * M2CM / S2US * k0 / betaj;	//phase velocity, cm/us
    tl_data.att = Rm/tl_data.rO/Z0*sqrt(1-pow(k11/k0,2))*(pow(k11/k0, 2)+1/(p11*p11-1));
  }
  cout << endl << endl << "Circ Waveguide with Radius: " << tl_data.rO << " cm" << endl;
  cout << " Circ WG TE11 Cutoff (Min) Freq: " << k11 / 2. / M_PI * Clight * M2CM * 1e-9 << " GHz " << endl;
  cout << " Circ WG TM01 Cutoff (Min) Freq: " << k01 / 2. / M_PI * Clight * M2CM * 1e-9 << " GHz " << endl;
  cout << " Circ WG Wave Impedance: " << tl_data.Zw << " Ohms" << endl;
  cout << " Skin Depth " << tl_data.skinD*1e6 << " micron " << endl;
  cout << " Circ WG Attenuation: " << tl_data.att << " Nepers per cm" << endl;
  cout << " Circ WG Group Velocity: " << tl_data.vg << " cm/us" << endl;
  cout << " Circ WG Phase Velocity: " << tl_data.vp << " cm/us" << endl;
  print_circ_wg_power(k0);
}

void print_circ_wg_power(double k0)
{
  //returns the time-average power radiated into the TE_11 
  //mode of an infinite circ. waveguide 
  // from an electron at the center of the guide with:
  double R = 0.045;//cm, cyclotron radius 
  TVector3 pos(0,0,0);//cm, guiding center position
  pos.SetZ(0);
  double radius = pos.Mag();
  double c = Clight * M2CM ;//cm/s
  double p11 = 1.841;            //1st zero of the derivate of bessel function
  //k0 = omega_cyclotron/c;//angular wavenumber of free-space cyclotron radiatin
  //k11 is angular wavenumber for cutoff frequency for TE11
  double k11 = p11 / tl_data.rO;
  double e_amp = 1.63303/tl_data.rO;//cm
  double Jp = TMath::BesselJ0(k11 * radius) - TMath::BesselJ1(k11 * radius) / k11 / radius;
  if (radius == 0) {
    Jp = 1.0/2; 
  }
  double power = tl_data.Zw / 4. * pow( Echarge * R * k0 * c * e_amp * Jp,2);
  cout << "Power at center of Circ WG: " << power*1e15 << " fW " << endl;

}

int get_circ_wg_efield_vert(TVector3 &pos, TVector3 &efield)
{
  //returns the TE_11 efield inside an infinite circ. waveguide 
  //efield normalized the jackson way with 1/cm units
  //waveguide extends in x-direction
  // wavelength of 27 GHz radiation is 1.1 cm
  double p11 = 1.841;            //1st zero of the derivate of bessel function
  //k11 is angular wavenumber for cutoff frequency for TE11
  double k11 = p11 / tl_data.rO;
  //convert position to cylindrical
  pos.SetZ(0);
  double radius = pos.Mag();
  double phi = pos.Phi();//azimuthal position
  //double phi = acos(pos.Y()/radius);//azimuthal position, see definition of phase
  double J1 = TMath::BesselJ1(k11 * radius);//this term cancels in dot product w/ vel
  double Jp = TMath::BesselJ0(k11 * radius) - TMath::BesselJ1(k11 * radius) / k11 / radius;
  double e_amp = 1.63303/tl_data.rO;//cm
  int status = 0;

  if (radius > tl_data.rO) {
    cout << "Problem!!! Electron hit a wall! At radius " << radius << endl;
    status = 1;
  }
  efield.SetX(e_amp * (J1 / k11 / radius * cos(phi) * sin(phi) - Jp * sin(phi) * cos(phi)));
  efield.SetY(e_amp * (J1 / k11 / radius * sin(phi) * sin(phi) + Jp * cos(phi) * cos(phi)));
  efield.SetZ(0);                //only true for TE mode
  if (radius == 0) {
    Jp = 1.0/2; 
    phi = 0;
    efield.SetX(e_amp * (1. / 2 * cos(phi) * sin(phi) - Jp * sin(phi) * cos(phi)));
    efield.SetY(e_amp * (1. / 2 * sin(phi) * sin(phi) + Jp * cos(phi) * cos(phi)));
  }
  return status;
}

int get_circ_wg_efield_horiz(TVector3 &pos, TVector3 &efield)
{
  //returns the TE_11 efield inside an infinite circ. waveguide 
  //efield normalized the jackson way with 1/cm units
  //waveguide extends in x-direction
  // wavelength of 27 GHz radiation is 1.1 cm
  double p11 = 1.841;            //1st zero of the derivate of bessel function
  //k11 is angular wavenumber for cutoff frequency for TE11
  double k11 = p11 / tl_data.rO;
  //convert position to cylindrical
  pos.SetZ(0);
  double radius = pos.Mag();
  double phi = pos.Phi();
  //double phi = acos(pos.Y()/radius);//azimuthal position, see definition of phase
  double J1 = TMath::BesselJ1(k11 * radius);//this term cancels in dot product w/ vel
  double Jp = TMath::BesselJ0(k11 * radius) - TMath::BesselJ1(k11 * radius) / k11 / radius;
  double e_amp = 1.63303/tl_data.rO;//cm
  int status = 0;

  if (radius > tl_data.rO) {
    cout << "Problem!!! Electron hit a wall! " << endl;
    status = 1;
  }
  efield.SetX(- e_amp * (J1 / k11 / radius * cos(phi) * cos(phi) + Jp * sin(phi) * sin(phi)));
  efield.SetY(e_amp * (-J1 / k11 / radius * cos(phi) * sin(phi) + Jp * cos(phi) * sin(phi)));
  efield.SetZ(0);                //only true for TE mode
  if (radius == 0) {
    Jp = 1.0/2; 
    phi = 0;
    efield.SetX(- e_amp * (1. / 2 * cos(phi) * cos(phi) + Jp * sin(phi) * sin(phi)));
    efield.SetY(e_amp * (- 1. / 2 * cos(phi) * sin(phi) + Jp * cos(phi) * sin(phi)));
  }
  cout << endl;
  
  return status;
}


