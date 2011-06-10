using namespace std;
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <string>
#include <cstring>
#include <iostream>
#include "paramanage.h"
#include "sim_scatter.h"

static int sim_help_save_data_count = 0; // counter to count between data to save
static char scatlog_filename[255];     // filename for scatter log
static int scatlog_set = 0;           // is scatlog_filename set?

void do_debug(char* nr, struct particle_data* pdata)
{
  printf("%s: ",nr);
  printf("%f  %f  %f  %f",(*pdata).e_para,(*pdata).e_perp,
	 (*pdata).e_accel,(*pdata).e_syncro);
  printf("\n");
}


void sim_help_init_scatter_file(char* filename)
{
  FILE* f_log;                    // parameter output file
  
  strcpy(scatlog_filename,filename);
  strcat(scatlog_filename,".scatlog");
  // data input from file with additional extension .scatlog
  scatlog_set = 1;

  f_log = fopen(scatlog_filename, "w");
  fprintf(f_log,"process Eloss[eV] angle[degree] Ekin[eV] tof[us] Etot[eV] x[cm] r[cm]\n\n"); 
  fclose(f_log);  
}


void sim_help_store_scatter_data(struct particle_data* particle, double e_loss, double angle, int process)
     // stores data to the scatter log file "adi_scatter.log"
{
  FILE* f_log;                    // parameter output file

  if (scatlog_set) 
    {
      f_log = fopen(scatlog_filename, "a");
      
      fprintf(f_log,"%d %8.3f %6.2f %10.3f  %e  %8.3f",process,
	      e_loss,angle,(*particle).e_kin,
	      (*particle).time_of_flight,
	      e_loss-(*particle).e_pot);
      fprintf(f_log," %8.3f %8.3f",
	      (*particle).position[0],
	      sqrt((*particle).position[1]*(*particle).position[1]+
		   (*particle).position[2]*(*particle).position[2]));
      fprintf(f_log,"\n"); 
      
      fclose(f_log);  
    }
  else cout << endl << "ERROR: scatlog filename not set!" << endl << endl; 
}



/******************** show_run_file ***********************/
void sim_help_show_run_file(char* filename)
{
  // Loop for output of runfile entries

  char run_filename[255];         // filename for parameter list
  char charline[255];            // commentline in parameter file

  double start_pos[3];
  double e_start,starting_theta,starting_phi,mass,charge,phase;
  int repeat_number = 0;

  FILE* f_run;                    // parameter input file

  strcpy(run_filename,filename);
  strcat(run_filename,".run");   // data input from file with extension .run

  f_run = fopen(run_filename, "r");

  printf("\n-----------------General Simulation Parameters------------------\n");
  printf("Particle Tracking parameter file: %s\n",run_filename);

  if (f_run == (FILE *)0)                 // disk error handling
    {
      fprintf(stderr, "ERROR: Can't read parameter file %s\n",run_filename);
      exit(1);
    }
  else 
    {
      fscanf(f_run,"%s\n",charline);  // step over first line

      // use ini-parameters for initialization of particle data
      printf("\nMAX_MIRRORS       = %d \n",parameter.max_mirrors);
      printf("MAX_TOF[s]        = %4.2f \n",parameter.max_tof_in_sec);
      printf("E_PARA_MIN        = %11.9f \n",parameter.e_para_min);
      printf("MAX_STEP_LENGTH   = %2.3f \n",parameter.max_step_length);
      printf("Dipol Value       = %4.2f V/m\n",parameter.dipole_value);
      
      printf("Calculation order is %d = ",parameter.calc_order);
      if (parameter.calc_order == 0) printf("all drifts\n");
      if (parameter.calc_order == 1) printf("only ExB drift\n");
      if (parameter.calc_order == 2) printf("only curvature and gradient drift\n");
      if (parameter.calc_order > 2) printf("no drift\n");
    }
  printf("\nParameters for each run:\n");
  printf("   Multi     X         Y         Z      Ekin    theta  phi  mass  charge phase\n");

  while (fscanf(f_run,"%d %lf %lf %lf %lf %lf %lf %lf %lf %lf\n",
		&repeat_number,&start_pos[0],&start_pos[1],&start_pos[2],
		&e_start,&starting_theta,&starting_phi,
		&mass,&charge,&phase)!=EOF)
    {// read data out of RUN-file until EOF
      printf("%8d %9.3f %9.3f %9.3f %7.1f %4.1f %4.1f %3.1f %3.1f %4.1f\n",
	     repeat_number,start_pos[0],start_pos[1],start_pos[2],
	     e_start,starting_theta,starting_phi,mass,charge,phase);
    }
  printf("----------------------------------------------------------------\n");
  fclose(f_run);           // close all files
}

void sim_help_force_save_data()
{
  sim_help_save_data_count = parameter.save_every;
}

void sim_help_reset_save_data()
{
  sim_help_save_data_count = 0;
}

/******************** savedata ***********************/
void sim_help_save_particle_data(FILE* f_fly,struct particle_data* particle)
{	// save one datablock to harddisc
  if (((*particle).position[0] >= -(*particle).max_step_length) && 
      ((*particle).position[0] <= (*particle).max_step_length))
    sim_help_save_data_count = SAVE_EVERY; // force data saved around analysing plane
  if (((*particle).position[0] >= (SPEC_OUT-(*particle).max_step_length*10)) || 
      ((*particle).position[0] <= (SPEC_IN+(*particle).max_step_length*10)))
    sim_help_save_data_count = SAVE_EVERY; // force data saved at source and detector
  if (sim_help_save_data_count == SAVE_EVERY)
    {
      sim_help_save_data_count = 0;
      fprintf(f_fly,"%1.9f %1.9f %1.9f ",
	      (*particle).position[0],     // 1
	      (*particle).position[1],     // 2
	      (*particle).position[2]);    // 3
      fprintf(f_fly,"%1.9f ",
	      (*particle).cyclrad);        // 4
      fprintf(f_fly,"%1.9f %1.9f ",
	      (*particle).e_pot,           // 5
	      (*particle).e_kin);          // 6
      fprintf(f_fly,"%1.9f %1.9f ",
	      (*particle).e_para,          // 7
	      (*particle).e_perp);         // 8
      fprintf(f_fly,"%1.9e %1.9e ",
	      (*particle).e_accel_sum,    // 9
	      (*particle).e_syncro_sum);    // 10
      fprintf(f_fly,"%1.9f %1.15e %1.15e %1.12e %e ",
	      (*particle).b_value,         // 11
	      (*particle).phase,     // 12
	      (*particle).time_of_flight,  // 13
	      scatter_get_ratio(),         // 14
	      get_wq((*particle).e_kin));  // 15
       fprintf(f_fly,"%1.9f %1.9f %1.9f %1.9e ", 
 	      (*particle).v_para[0],       // 16 
 	      (*particle).v_para[1],       // 17 
 	      (*particle).v_para[2],      // 18 
 	      (*particle).omega);      // 19 
/*       fprintf(f_fly,"%1.9f %1.9f ",  */
/* 	      (*particle).v_para_value,    // 17  */
/* 	      (*particle).v_perp_value);   // 18  */
/*       fprintf(f_fly,"%1.9f %1.9f %1.9f ", */
/* 	      (*particle).exb_vel[0],      // 19 */
/* 	      (*particle).exb_vel[1],      // 20 */
/* 	      (*particle).exb_vel[2]);     // 21 */
/*       fprintf(f_fly,"%1.9f %1.9f %1.9f ", */
/* 	      (*particle).rxb_vel[0],      // 22 */
/* 	      (*particle).rxb_vel[1],      // 23 */
/* 	      (*particle).rxb_vel[2]);     // 24 */
/*       fprintf(f_fly,"%1.9f ", */
/* 	      (*particle).curv_rad);       // 25 */
/*       fprintf(f_fly,"%1.9f %1.9f %1.9f ", */
/* 	      (*particle).dbxb_vel[0],     // 26 */
/* 	      (*particle).dbxb_vel[1],     // 27 */
/* 	      (*particle).dbxb_vel[2]);    // 28 */
/*       fprintf(f_fly,"%1.9f %1.9f ", */
/* 	      (*particle).v_signum,        // 29 */
/* 	      (*particle).shrinkfactor);   // 30 */
       fprintf(f_fly,"\n");
    }
  else sim_help_save_data_count++;
}
