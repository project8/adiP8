using namespace std;
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <string>
#include <iostream>
#include "paramanage.h"
#include "sim_core.h"
#include "sim_help.h"
#include "sim_scatter.h"
#include "vector_tool.h"
#include "el_pa_tool.h"
#include "mag_pa_tool.h"


/******************** tracking_loop ***********************/
void tracking_loop(char* filename)
{
  // Loop for Particle Tracking of more than one Particle

  struct particle_data particle;    // structure containing particle data
  int particle_count = 0;           // counter for particles
  int repeat_number = 1;            // number of repeated simulations
  int repeated = 0;                 // repeater counter
  char track_filename[255];         // filename for tracking data
  char run_filename[255];         // filename for parameter list
  char charline[255];            // commentline in parameter file
  char set_filename[255];         // dilename for set results
  FILE* f_run;                    // parameter input file
  FILE* f_track;                  // output file
  FILE* f_set;                    //output of set results

  strcpy(run_filename,filename);
  strcat(run_filename,".run");   // data input from file with extension .run

  f_run = fopen(run_filename, "r");

  printf("Particle Tracking parameter file: %s\n",run_filename);

  if (f_run == (FILE *)0)                 // disk error handling
    {
      fprintf(stderr, "ERROR: Can't read parameter file %s\n",run_filename);
      exit(1);
    }
  else 
    {
      // use ini-parameters for initialization of particle data
      particle.calc_order = parameter.calc_order;
      particle.max_mirrors = parameter.max_mirrors;
      particle.max_tof = parameter.max_tof_in_sec;
      particle.e_para_min = parameter.e_para_min;
      particle.max_step_length = parameter.max_step_length;
      set_dipole(parameter.dipole_value); 
      // sets the dipole magnifying value for the efield module

      fscanf(f_run,"%s\n",charline);
      printf("MAX_MIRRORS       = %d \n",particle.max_mirrors);
      printf("MAX_TOF[s]        = %4.2f \n",particle.max_tof);
      printf("E_PARA_MIN        = %11.9f \n",particle.e_para_min);
      printf("MAX_STEP_LENGTH   = %2.3f \n",particle.max_step_length);
      printf("Dipol Value       = %4.2f V/m\n",get_dipole());
      
      printf("calculate       : ");
      if (particle.calc_order == 0) printf("all drifts\n");
      if (particle.calc_order == 1) printf("only ExB drift\n");
      if (particle.calc_order == 2) printf("only curvature and gradient drift\n");
      if (particle.calc_order > 2) printf("no drift\n");
    }
  printf("\n-----------------------------------------\n");

  while (fscanf(f_run,"%d %lf %lf %lf %lf %lf %lf %lf %lf\n",
		&repeat_number,
		&particle.start_pos[0],&particle.start_pos[1],&particle.start_pos[2],
		&particle.e_start,&particle.starting_theta,&particle.starting_phi,
		&particle.mass,&particle.charge)!=EOF)
    {// read data out of file until file ends

      particle_count++;
      particle.mass = particle.mass*M0Clight2;
      if (particle.e_start == 0.)     //< particle.e_para_min)
      {
	particle.e_start = particle.e_para_min;
	printf("warning: starting energie low, set to E_PARA_MIN for calculation!\n");
      }

      sprintf(set_filename,
	      "%s.set%d",filename,particle_count);        
      f_set = fopen(set_filename, "w");
      
      for (repeated = 1; repeated <= repeat_number; repeated++)
	{

	  printf("starting with run %d of run-set %d.\n\n",
		 repeated,particle_count);

	  sprintf(track_filename,
		  "%s.track%d",
		  filename,particle_count);        // data output to file with extension .track
	  f_track = fopen(track_filename, "w");
	  sim_help_init_scatter_file(track_filename);

	  printf("Particle Tracking output file: %s\n",track_filename);
	  printf("Simulating set %d of total %d.\n",repeated,repeat_number);
	  printf("Startparameters x=%.3f y=%.3f z=%.3f E=%.1f theta=%.1f\n",
		 particle.start_pos[0],particle.start_pos[1],
		 particle.start_pos[2],particle.e_start,particle.starting_theta);
	  
	  if (f_track == (FILE *)0)                 // disk error handling
	    {
	      fprintf(stderr, "ERROR: Can't write output file %s\n",track_filename);
	      exit(1);
	    }
	  
	  single_track(f_track, &particle); // start a single track simulation

	  fprintf(f_set,"%.3f  %.3f  %.3f  %.3f  %.3f  %.3f\n",
		  particle.position[0],particle.e_start,
		  particle.e_kin,particle.e_para,particle.e_perp,
		  asin(sqrt(particle.e_perp/particle.e_kin))*180/M_PI);

	  printf("\nfinished run %d of run-set %d.\n",
		 repeated,particle_count);
	  printf("-----------------------------------------\n");
	  
	  fclose(f_track);
	}
      fclose(f_set);
    }
  fclose(f_run);           // close all files
} // end of tracking_loop



int trapping_loop(char* filename)
                  //parameters of the stored electrons
{

  struct particle_data particle;    // structure containing particle data
  double position[3];
  double startx = 0.,stopx = 0.,stepx = 0.;
  double theta,energy;
  double phi;
  double y_max;
  double y_limit;
  int e_tag =0;
  double trapped = 0.;
  int trap_tag = 0;
  char tlog_filename[255];
  char tvol_filename[255];
  FILE *f_tlog;
  FILE *f_tvol;

  strcpy(tlog_filename,filename);
  strcat(tlog_filename,".tlog");

  strcpy(tvol_filename,filename);
  strcat(tvol_filename,".tvol");

  f_tvol=fopen(tvol_filename,"w");
  
  if (f_tvol == (FILE *)0)
    {
      fprintf(stderr, "ERROR: Can't open trapping-volume data file %s\n",
	      tvol_filename);
      exit(1);
    }
  else
    {
      f_tlog=fopen(tlog_filename,"w");
      
      if (f_tlog == (FILE *)0)
	{
	  fprintf(stderr, "ERROR: Can't open trapping-volume log file %s\n",
		  tlog_filename);
	  exit(1);
	}
      else {
	// disable any scattering or energy loss functions
	set_parameter("ENABLE_PARA_ENERGY_LOSS" , 0.);
	set_parameter("ENABLE_PERP_ENERGY_LOSS" , 0.);
	set_parameter("RESIDUAL_GAS_PRESSURE" , -1.);
	set_parameter("REL_START_ANGLE",1.);
	
	particle.max_mirrors = parameter.trap_max_mirrors;  
	// tracking stops after this count of mirrors

	particle.mass = M0Clight2*parameter.trap_mass;      // mass
	particle.charge = parameter.trap_charge;            // charge
	particle.max_step_length = parameter.trap_max_step_length; 
	// steplength

	particle.e_para_min = parameter.trap_e_para_min; 
	// minimum of parallel energy for mirror points

	particle.max_tof = parameter.trap_max_tof;
	// break after this tof

	particle.calc_order = parameter.trap_calc_order; 
	// order of drift calculations
	
	// x-direction loop
	y_limit = parameter.trap_stop_y;
	
	if (parameter.trap_neg_y == 1) 
	  cout << "Remark: negative volume included!" << endl;
	else
	  cout << "Remark: only positive volume included!" << endl;

	if ((parameter.trap_y_plane == 1) 
	    && (parameter.trap_z_plane == 2))
	  cout << "Remark: calculate in x-y-plane!" << endl;
	else 
	  cout << "Remark: calculate in x-z-plane!" << endl; 
	
	stepx = -sqrt(parameter.trap_step_x * parameter.trap_step_x);
	// generate absolute value of step in x direction.

	if (parameter.trap_start_x < parameter.trap_stop_x)
	  {  // if start left of stop exchange values
	    startx = parameter.trap_stop_x;
	    stopx = parameter.trap_start_x;
	  }
	else
	  {
	    startx = parameter.trap_start_x;
	    stopx = parameter.trap_stop_x;
	  }

	for (position[0]=startx; position[0]>=stopx; position[0] += stepx)
	  {
	    position[parameter.trap_z_plane]=0.;
	    y_max = 0.;
	    e_tag = 0;
	    do
	      {
		position[parameter.trap_y_plane] = y_max;
		epot3d(position,&e_tag);
		if (e_tag == 0) y_max = y_max + parameter.el_mm_unit/10.;
	      }
	    while ((y_max <= y_limit) && (e_tag == 0));
	    y_limit = y_max;
	    printf("x= %f    y_limit= %f\n",position[0],y_max);
	    position[parameter.trap_y_plane] = parameter.trap_start_y;

	    if (parameter.trap_neg_y == 1)
	      position[parameter.trap_y_plane] = -parameter.trap_stop_y;

	    e_tag = 0;
	    epot3d(position,&e_tag);
	    fprintf(f_tvol," \n");

	    // while (/*(e_tag == 0) &&*/ (position[y_plane] <= trap_stopy))
	    
	    while ((position[parameter.trap_y_plane] <= parameter.trap_stop_y))
	      // y-direction loop
	      {   
		for(energy = parameter.trap_energy_start;
		    energy >= parameter.trap_energy_end;
		    energy = energy*0.5)  // energy loop
		  {

		    fprintf(f_tvol,"%f %f %f ",
			    position[0],
			    position[parameter.trap_y_plane],
			    energy);

		    theta=parameter.trap_theta_start;
		    phi=parameter.trap_phi_start;
		    particle.starting_theta = theta;
		    particle.starting_phi = phi;
		    particle.e_start = energy; 
		    
		    if (particle.e_start == 0.)
		      particle.e_start=energy+particle.e_para_min;

		    particle.start_pos[0] = position[0];

		    particle.start_pos[parameter.trap_y_plane]
		      = position[parameter.trap_y_plane];

		    particle.start_pos[parameter.trap_z_plane]
		      = position[parameter.trap_z_plane];

		    trapped = 0.;
		    
		    do {//theta loop
			
			//*********phi loop to be placed here***********
			
			trap_tag = 0;

			// do the actual tracking simulation for one electron
			if (sqrt(position[parameter.trap_y_plane]*
				 position[parameter.trap_y_plane]) <= y_limit)
			  {
			    trap_tag = single_track_trap(&particle);
			    //trap_tag = 1; // for overriding single_track_trap
			  }

			switch(trap_tag)
			  {
			  case 2 :
			    fprintf(f_tlog,"ADIPARK broken loop at (x,y,z,E,a) %f %f %f %f %f\n",position[0],position[1],position[2],energy,theta);
			    break;
			  case 3 :
			    fprintf(f_tlog,"ADIPARK relativistic ERROR at (x,y,z,E,a) %f %f %f %f %f\n",position[0],position[1],position[2],energy,theta);
			    break;
			  }
			
			if (trap_tag == 1)
			  {
			    trapped++;
			    theta=theta-parameter.trap_theta_step;
			    particle.starting_theta = theta;
			    if (theta <= 1.) theta = -1;
			  }
			else theta= -1;
			
		    } while (theta> -1);

		    if (sqrt(position[parameter.trap_y_plane]
			     *position[parameter.trap_y_plane]) > y_limit)
		      trapped = 0.;
		    
		    fprintf(f_tvol,"%f\n",trapped);
		    fflush(f_tvol);
		  }
		position[parameter.trap_y_plane]=
		  position[parameter.trap_y_plane]+parameter.trap_step_y; 
		
		e_tag = 0;
		epot3d(position,&e_tag);		
	      }
	  }
	fclose(f_tlog);
      }
      fclose(f_tvol);
    }
  return 0;
}    



/******************** transmission_loop ***********************/
void transmission_loop(char* filename)
{
  // Loop for Particle Tracking of more than one Particle

  struct particle_data particle;    // structure containing particle data
  int particle_count = 0;           // counter for particles 
  char track_filename[255];         // filename for tracking data
  FILE* f_track;                    // output file
  //  double pos_spec_in[3],b_spec_in[3];
  double test_pos[3];
  double b_abs_spec_in, u_spec_in;   // spectrometer entrance
  int e_tag = 0;

  particle.calc_order = 3; // no drifts
  ENABLE_EPOT = 1; 
  set_dipole(0.);
  particle.max_mirrors = 2;
  particle.max_tof = 0.1;
  particle.e_para_min = 0.01;
  particle.max_step_length = 0.1;
  for (particle_count=0;particle_count<N_TRANSMISSION;particle_count++){
    particle.start_pos[0] = 0.0;
    particle.start_pos[1] = MAX_RADIUS/N_TRANSMISSION*particle_count;
    particle.start_pos[2] = 0.0;
    printf("particle.start.pos %f %f %f\n", particle.start_pos[0],particle.start_pos[1],particle.start_pos[2]);
    get_bfield(particle.start_pos, particle.b_start); // starting bvec vector at start position
    printf("b start.pos %f %f %f\n", particle.b_start[0],particle.b_start[1],particle.b_start[2]);
    printf("particle.start.pos %f %f %f\n", particle.start_pos[0],particle.start_pos[1],particle.start_pos[2]);
    particle.b_start_value = absvalue(particle.b_start);
                                                // and its value 
    particle.u_start = epot3d(particle.start_pos,&e_tag);// potential at start position
    printf("b start.pos %f, u start.pos %f\n",particle.b_start_value,particle.u_start);

    test_pos[0] = SPEC_IN;
    test_pos[1] = 0.0;
    test_pos[2] = 0.0;
    printf("u_pot at entrance spectrometer %f\n",epot3d(test_pos,&e_tag));
    test_pos[0] = SPEC_OUT;
    test_pos[1] = 0.0;
    test_pos[2] = 0.0;
    printf("u_pot at exit spectrometer %f\n",epot3d(test_pos,&e_tag));
    // spectrometer entrance = pinch magnet
    b_abs_spec_in = B_PINCH;
    u_spec_in = U_PINCH;   // potential at PINCH 
    particle.e_perp = (-particle.u_start + u_spec_in) * particle.b_start_value/b_abs_spec_in;
    particle.e_para = 1.0;
    particle.e_start = particle.e_perp+particle.e_para;
    particle.starting_theta = -180.0/M_PI*asin(sqrt(particle.e_perp/particle.e_start));  
    printf("particle.e_start %f, particle.starting_theta %f\n",particle.e_start, particle.starting_theta);
    //    particle.starting_theta = -89.0;
    particle.starting_phi = 0.0;
    particle.mass = M0Clight2;
    particle.charge = 1.0;
    printf("starting run %d.\n\n",particle_count);

    sprintf(track_filename,"%s.track%d",filename,particle_count);       
        // data output to file with extension .track

    f_track = fopen(track_filename, "w");
    printf("Particle Tracking output file: %s\n",track_filename);
    if (f_track == (FILE *)0){                // disk error handling
      fprintf(stderr, "ERROR: Can't write output file %s\n",track_filename);
      exit(1);
    }
    single_track(f_track, &particle); 
        // start a single track simulation
    printf("\nfinished run %d.\n",particle_count);
    printf("-----------------------------------------\n");
    fclose(f_track);
  }
} // end of transmission_loop


void scan_field_gradient(char* filename)
{
  int x,y,e_tag;
  int xmax = (int)parameter.spec_out;
  int ymax = (int)parameter.max_radius;
  int number = 0;
  double evec[3];
  double pos[3];
  double epot;
  char pot_filename[255];
  FILE* f_pot;                  // output file
  
  strcpy(pot_filename,filename);
  strcat(pot_filename,".pot");   // data input from file with extension .run
  
  f_pot = fopen(pot_filename, "w");
  
  printf("Potential and Field datafile: %s\n",pot_filename);
  
  for (x = 0; x <= xmax; x++)
    {
      for (y = 0; y <= ymax; y++)
	{
	  number++;
	  pos[0] = (double)x;
	  pos[1] = (double)y;
	  pos[2] = 0.;
	  e_tag = 0;
	  efield3d(pos,evec,&e_tag);
	  epot = epot3d(pos,&e_tag);

	  if (number % 100 == 0)
	    cout << "  writing: " 
		 << (int)(100.*(float)number/(float)(xmax*ymax)) 
		 << " %   \r" << flush;
	  
	  fprintf(f_pot,"%e %e %e %e %e %e %e %d\n",
		  pos[0], pos[1], pos[2],
		  epot, evec[0], evec[1], evec[2], e_tag); 
	}
      fprintf(f_pot,"\n");
    }
  cout << "    DONE.         " << endl << endl; 
  fclose(f_pot);
}

void scan_Bfield(char* filename)
{
  int x,y,number = 0;
  double bvec[3];
  double pos[3];
  char bvec_filename[255];
  FILE* f_bvec;                  // output file
  
  strcpy(bvec_filename,filename);
  strcat(bvec_filename,".bvec");   
  
  f_bvec = fopen(bvec_filename, "w");
  
  printf("Bfield datafile: %s\n",bvec_filename);
  
  for (x = 0; x <= (int)parameter.spec_out; x++)
    {
      for (y = 0; y <= (int)parameter.max_radius; y++)
	{
	  number++;
	  pos[0] = (double)x;
	  pos[1] = (double)y;
	  pos[2] = 0.;
	  get_bfield(pos, bvec);

	  if (number % 100 == 0)
	    cout << "  writing: " 
		 << (int)(100.*(float)number/
			  (float)(parameter.spec_out*parameter.max_radius)) 
		 << " %   \r" << flush;
	  
	  fprintf(f_bvec,"%e %e %e %e %e %e %e\n",
		  pos[0], pos[1], pos[2],
		  bvec[0], bvec[1], bvec[2],
		  sqrt(bvec[0]*bvec[0]+bvec[1]*bvec[1]+bvec[2]*bvec[2])); 
	}
      fprintf(f_bvec,"\n");
      fflush(f_bvec);
    }
  cout << "    DONE.         " << endl << endl; 
  fclose(f_bvec);
}




/******************** mc_tracking ***********************/
void mc_tracking(char* filename)
{
  // Loop for Particle Tracking of more than one Particle

  const double minr = 0.;
  const double maxr = 15.;      // in cm
  const double fluxtube = 192.; // in tesla*cm*cm

  struct particle_data particle;    // structure containing particle data
  double fluxhere = 0.;
  double minx, maxx;
  int particle_count = 0;           // counter for particles
  int repeat_number = 1;            // number of repeated simulations
  int repeated = 0;                 // repeater counter
  char track_filename[255];         // filename for tracking data
  char run_filename[255];         // filename for parameter list
  char charline[255];            // commentline in parameter file
  char set_filename[255];         // dilename for set results
  FILE* f_run;                    // parameter input file
  FILE* f_track;                  // output file
  FILE* f_set;                    //output of set results

  strcpy(run_filename,filename);
  strcat(run_filename,".run");   // data input from file with extension .run

  f_run = fopen(run_filename, "r");

  printf("\nDoing Monte Carlo tracking of trapped particles!\n\n");
  printf("Particle Tracking parameter file: %s\n",run_filename);

  if (f_run == (FILE *)0)                 // disk error handling
    {
      fprintf(stderr, "ERROR: Can't read parameter file %s\n",run_filename);
      exit(1);
    }
  else 
    {
      // use ini-parameters for initialization of particle data
      particle.calc_order = parameter.calc_order;
      particle.max_mirrors = parameter.max_mirrors;
      particle.max_tof = parameter.max_tof_in_sec;
      particle.e_para_min = parameter.e_para_min;
      particle.max_step_length = parameter.max_step_length;
      set_dipole(parameter.dipole_value);
      // sets the dipole magnifying value for the efield module

      set_parameter("SAVE_EVERY",1e12);
      // set write control to almost infinity, because
      // track files are not used here

      fscanf(f_run,"%s\n",charline);
      printf("MAX_MIRRORS       = %d \n",particle.max_mirrors);
      printf("MAX_TOF[s]        = %4.2f \n",particle.max_tof);
      printf("E_PARA_MIN        = %11.9f \n",particle.e_para_min);
      printf("MAX_STEP_LENGTH   = %2.3f \n",particle.max_step_length);
      printf("Dipol Value       = %4.2f V/m\n",get_dipole());
      
      printf("calculate       : ");
      if (particle.calc_order == 0) printf("all drifts\n");
      if (particle.calc_order == 1) printf("only ExB drift\n");
      if (particle.calc_order == 2) printf("only curvature and gradient drift\n");
      if (particle.calc_order > 2) printf("no drift\n");      
    }
  printf("\n-----------------------------------------\n");

  while (fscanf(f_run,"%d %lf %lf %lf %lf %lf %lf %lf %lf\n",
		&repeat_number,
		&particle.start_pos[0],&particle.start_pos[1],
		&particle.start_pos[2],&particle.e_start,
		&particle.starting_theta,&particle.starting_phi,
		&particle.mass,&particle.charge)!=EOF)
    {// read data out of file until file ends
      particle_count++;
      particle.mass = particle.mass*M0Clight2;

      minx = particle.start_pos[0];
      maxx = particle.start_pos[1];
      // use run-file x value as left xmargin
      // and run-file y value as right xmargin
      if (minx >= maxx)
	{
	  cout << endl << "ERROR in MC_Tracking: minx >= maxx, check run-file parameters." << endl << endl << flush;
	  exit (1);
	}
      
      if (particle.e_start == 0.)     //< particle.e_para_min)
	{
	  particle.e_start = particle.e_para_min;
	  printf("warning: starting energie low, set to E_PARA_MIN for calculation!\n");
	}

      sprintf(set_filename,
	      "%s.mc%d",filename,particle_count);        
      f_set = fopen(set_filename, "w");
      
      for (repeated = 1; repeated <= repeat_number; repeated++)
	{

	  printf("starting with run %d of run-set %d.\n\n",
		 repeated,particle_count);

	  sprintf(track_filename,
		  "%s.track%d",
		  filename,particle_count);  
	  // data output to .track file

	  f_track = fopen(track_filename, "w");
	  sim_help_init_scatter_file(track_filename);

	  printf("Particle Tracking output file: %s\n",track_filename);
	  printf("Simulating set %d of total %d.\n",repeated,repeat_number);

	  // now dice the startpos and startangle numbers
	  particle.start_pos[0] = minx+get_std_rand_1()*(maxx-minx);

	  do{
	    particle.start_pos[0] = minx + get_std_rand_1()*(maxx-minx);
	    // random between specin <= x <= specout
	    
	    particle.start_pos[1] = sqrt(get_std_rand_1()*maxr*maxr);
	    // random between 0 <= ran_r <= maxr;
	    
	    particle.start_pos[2] = 0.;
	    // no z koordinate needed

	    get_bfield(particle.start_pos,particle.b_start);
	    // get magnetic field at position
	    
	    fluxhere = M_PI*particle.start_pos[1]*particle.start_pos[1]
	      *absvalue(particle.b_start);
	    // calculate the magnetic flux with given values
	  } while ((fluxhere > fluxtube) || (particle.start_pos[1] < minr));
	  
	  particle.starting_theta = acos(get_std_rand_1())*180./M_PI;
	  // dice starting theta according to isotropic nature

	  printf("Startparameters x=%.3f y=%.3f z=%.3f E=%.1f theta=%.1f\n",
		 particle.start_pos[0],particle.start_pos[1],
		 particle.start_pos[2],particle.e_start,
		 particle.starting_theta);
	  
	  if (f_track == (FILE *)0)                 // disk error handling
	    {
	      fprintf(stderr, "ERROR: Can't write output file %s\n",track_filename);
	      exit(1);
	    }
	  
	  single_track(f_track, &particle); // start a single track simulation

	  if (particle.mirrors >= 2)
	    {
	      fprintf(f_set,"%10d %.3f %.3f %.3f %.3f %.3f %.3f %.3f %.3f %.3f\n",
		      repeated,particle.start_pos[0],
		      particle.start_pos[1],particle.starting_theta,
		      particle.position[0],particle.e_start,
		      particle.e_kin,particle.e_para,particle.e_perp,
		      asin(sqrt(particle.e_perp/particle.e_kin))*180/M_PI);
	      fflush(f_set);
	    }
	  else cout << "Particle NOT trapped! -> skipped" << endl;

	  printf("\nfinished run %d of run-set %d.\n",
		 repeated,particle_count);
	  printf("-----------------------------------------\n");
	  
	  fclose(f_track);
	}
      fclose(f_set);
    }
  fclose(f_run);           // close all files
} // end of mc_tracking
