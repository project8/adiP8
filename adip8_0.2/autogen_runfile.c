#include <iostream>
#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include <cmath>
#include "TRandom.h"

using namespace std;

Double_t dnde(Double_t ke)
{ // beta spectrum for mnu=0
  if (ke > 18575 || ke < 0) return 0;
  return sqrt(ke*ke+2*511000*ke)*(ke+511000)*(18575-ke)*(18575-ke);
}


int main(int argc, char* argv[])
{
  
  if (argc==1 || (argc > 1 && (strcmp(argv[1],"--help")==0 ||strcmp(argv[1],"-help")==0 || strcmp(argv[1],"-?")==0)))
    {
      cout << "Usage: autogen_runfile [nevents] [min_pa] [min_ene] [flatene]" << endl;
      cout << " pa in degrees (default 0); ene in eV (default 0)" << endl;
      cout << " flatene [0 (default) = tritium spectrum, 1=flat in e]" << endl; 
      exit(0);
    }

  cout << "repeat___x______y______z______ekin_____theta___phi____mass___charge\n"; 


  int nevents = atoi(argv[1]);
  double min_pa  =0;
  if (argc > 2) min_pa = atof(argv[2]);
  //if (argc > 2) min_pa = atof(argv[2])*3.14159265/180;
  double min_ene  =0;
  if (argc > 3) min_ene = atof(argv[3]);
  int flatene = 0;
  if (argc > 4) flatene = atoi(argv[4]);
    
  double specscale;
  if (min_ene > 3769.26) specscale = dnde(min_ene);
  else specscale = dnde(3769.26)*1.01;
  double y0=0.95, z0 =0.9;
  //double y0=0.95, z0 =0.9;
  double dy = 2*y0/nevents;
  double dz = 2*z0/(nevents-1);
  for (int i=0;i<=nevents;i++)
    {
    for (int j=0;j<nevents;j++)
      {
      double x=0;
      double y=y0-i*dy;
      double z=z0-j*dz;
      double ene, pa;

      if (flatene==0)
        {
        pa = acos(gRandom->Rndm()*cos(min_pa))*180/3.14159265;
	double ry;
	do {
	  ene= gRandom->Rndm()*(18575-min_ene)+min_ene;
	  ry = gRandom->Rndm()*specscale;
	} while (ry < dnde(ene))	 ;   
      }
      else 
	{
	  //ene = gRandom->Rndm()*(18575-min_ene)+min_ene;
	  ene = min_ene;
          pa = min_pa;
	}
      
      
      printf("%d %g %g %g %g %g %g %g %g\n",1,x,y,z,ene,pa, 0.,1.,1.);

    }
  }

}
