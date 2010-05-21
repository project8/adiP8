using namespace std;
#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <iostream>
#include "eH2.h"
#include "sim_scatter.h"
#include "paramanage.h"

double dx_per_path_init = 0.;
double dx_per_path_total = 0.;
double current_eloss = -1.;
double current_angle = -1.;
double current_sigma    = -1.;
double current_sigmaexc = -1.;
double current_sigmaion = -1.;
double current_sigmael  = -1.;
int current_process = -1;


double get_wq(double argument)
     // interface function to cross section module of
     // Ferenc Glueck.
     // Switch between different processes here!
{
  return 1e4*(sigmaion(argument)   // sigma for ionisation of H2
	      +sigmaexc(argument)  // sigma for excitation of H2
	      +sigmael(argument)); // sigma for elastic scattering
  // sigmaion return a m² value --> need to transform to cm²
}


double get_std_rand_1(void)
    // function gives a random number between 0 and 1
    // using only the C command
{
  int j;
  j=1+(int) (1e9*rand()/(RAND_MAX+1.0));
  
  return j/1e9;
}


int rand_seed=-1; // for random calculation
double get_random_1_cw(void)
     // function gives a random value between 0 and 1.
     // NO imput value needed.
     // numerical rec. version provided by Ch. Weinheimer
{
#define IA 16807
#define IM 2147483647
#define AM (1./IM)
#define IQ 127773
#define IR 2836
#define NTAB 32
#define NDIV (1+(IM-1)/NTAB)
#define RNMX (1.-1.2e-7)
  
  int j, k;
  static int iy=0, iv[NTAB];
  double temp;
  
  if (rand_seed<=0 || !iy){
    if (-rand_seed<1) rand_seed=1;
    else rand_seed=-rand_seed;
    for (j=NTAB+7; j>=0; j--){
      k=rand_seed/IQ;
      rand_seed=IA*(rand_seed-k*IQ)-IR*k;
      if (rand_seed<0) rand_seed+=IM;
      if (j<NTAB) iv[j]=rand_seed;
    }
    iy=iv[0];
  }
  k=rand_seed/IQ;
  rand_seed=IA*(rand_seed-k*IQ)-IR*k;
  if (rand_seed<0) rand_seed+=IM;
  j=iy/NDIV;
  iy=iv[j];
  iv[j]=rand_seed;
  if ((temp=AM*iy)>RNMX) return RNMX;
  else return temp;
  
#undef IA
#undef IM
#undef AM
#undef IQ
#undef IR
#undef NTAB
#undef NDIV
#undef RNMX
}


double get_residual_gas_density(void)
     // function to calculate the number of residual gas
     // particles depending on temperature and pressure
{
#define RGN_pressure parameter.residual_gas_pressure  // in mbar
#define RGN_temp     30.            // in Kelvin
#define RGN_k_Boltzmann 1.380658e-23 // Boltzmann's number
#define RGN_unit_conversion 1e-4     
  // converts mbar to N/m² and 1/m³ to 1/cm³

  if (RGN_pressure < 0.) return 0.;
  else return RGN_pressure * RGN_unit_conversion /
         (RGN_k_Boltzmann * RGN_temp);

#undef RGN_pressure 
#undef RGN_temp     
#undef RGN_atom_per_molecule 
#undef RGN_k_Boltzmann 
#undef RGN_unit_conversion     
}


void scatter_debug(void)
{
  int j;
  double ekin;
  
  dx_per_path_total = 2.;
  dx_per_path_init  = 1.; 
  ekin = 40000.;
  {
      for (j=1; j <= 50; j++)
	{
	  current_process = -1;
	  scatter_store_sigmas(ekin);
	  scatter_store_angle_eloss(ekin);
	  
	  //	  printf(" %f  %e  %e  %d\n",ekin,current_eloss,
	  //		 current_angle, current_process);
	}
    }

  //  printf(" %e  %e  %e\n",dx_per_path_total,dx_per_path_init,
  //	 scatter_get_ratio());
}


void scatter_init_rel_cross_section(void)
     // init vaiables for scattering detection
{
  //  scatter_test();

  dx_per_path_total = 0.;
  dx_per_path_init  = -log(get_random_1_cw());
  current_eloss     = -1.;
  current_angle     = -1.;
  current_sigma     = -1.;
  current_sigmaexc  = -1.;
  current_sigmaion  = -1.;
  current_sigmael   = -1.;
  current_process   = -1;
}


void scatter_add_step(double dx, double Ekin)
     // add dx/path fraction of current step to total fraction
{
  //  double meen_free_path = 1/(get_residual_gas_density()*get_wq(Ekin));
  if (parameter.residual_gas_pressure >= 0.)
    dx_per_path_total += get_residual_gas_density()*get_wq(Ekin)*dx;
}


int scatter_check_event(void)
     // check if dx/path fraction is reached already 
{
  if ((dx_per_path_total >= dx_per_path_init) && 
      (parameter.residual_gas_pressure >= 0.))
    return 1;
  else return 0;
}


void scatter_store_sigmas(double ekin)
     // remembers all sigma values in case of scatter event
{
  current_sigma     = -1.;
  current_sigmaexc  = -1.;
  current_sigmaion  = -1.;
  current_sigmael   = -1.;
  if (scatter_check_event())
    {
      current_sigma     = get_wq(ekin);
      current_sigmaexc  = sigmaexc(ekin)*1e4;
      current_sigmaion  = sigmaion(ekin)*1e4;
      current_sigmael   = sigmael(ekin)*1e4;
    }
}


int scatter_get_process(void)
     // returns the type of scattering process
     // 0 = elastic
     // 1 = ionization
     // 2 = excitation
     // -1 = not scattered
{
  double random;
  double edge1,edge2;
  
  if ((scatter_check_event()) && (current_process < 0.))
    {
      if (current_sigma >= 0.)
	{
	  random = get_random_1_cw();
	  
	  edge1 = current_sigmaion/current_sigma;
	  edge2 = (current_sigmaion+current_sigmaexc)/current_sigma;
	  
	  if (random < edge1) current_process = 1;
	  else if (random < edge2) current_process = 2;
	  else current_process = 0;

	  //	  printf(" %f  %f  %f  %d \n",edge1,edge2,random,current_process);
	}
    }
  return current_process;
}


void scatter_store_angle_eloss(double ekin)
     // returns the energy loss if particle has been scattered
     // -1 if not scattered
{
  if (scatter_check_event())
    {
      scatter_store_sigmas(ekin);
      switch(scatter_get_process())
	{
	case 0 :
	  randomel(ekin,&current_eloss,&current_angle);
	  break;
	case 1 :
	  randomion(ekin,&current_eloss,&current_angle);
	  break;
	case 2 :
	  randomexc(ekin,&current_eloss,&current_angle);
	  break;
	default :
	  printf("ERROR in scatter_get_eloss: process not set!\n");
	}
    }
}


double scatter_get_eloss(void)
     // returns precalculated energy loss
{
  return current_eloss;
}


double scatter_get_angle(void)
     // returns precalculated scattering angle
{
  return current_angle;
}


double scatter_get_ratio(void)
     // returns the ratio between current and destination value
{
  return dx_per_path_total/dx_per_path_init;
}


