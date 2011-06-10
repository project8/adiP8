using namespace std;
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <string>
#include <iostream>
#include "paramanage.h"

/******************** init_parameters ***********************/
void init_parameters(void)
{
  int para_count;
  
  printf("Init all ADIPARK 4 Parameters.");

  for (para_count = 0;para_count <= max_parameter+1; para_count++)
    {
      parameter.inited[para_count] = 0;
    }

  para_count = 0;
  sprintf(parameter.name[++para_count],"EL_MM_PER_UNIT");
  sprintf(parameter.name[++para_count],"EL_X_OFFSET_IN_CM");
  sprintf(parameter.name[++para_count],"EL_Y_OFFSET_IN_CM");
  sprintf(parameter.name[++para_count],"EL_Z_OFFSET_IN_CM");
  sprintf(parameter.name[++para_count],"MAG_MM_PER_UNIT");
  sprintf(parameter.name[++para_count],"MAG_X_OFFSET_IN_CM");
  sprintf(parameter.name[++para_count],"MAG_Y_OFFSET_IN_CM");
  sprintf(parameter.name[++para_count],"MAG_Z_OFFSET_IN_CM");
  sprintf(parameter.name[++para_count],"SPEC_IN");
  sprintf(parameter.name[++para_count],"SPEC_OUT");
  sprintf(parameter.name[++para_count],"MAX_RADIUS");
  sprintf(parameter.name[++para_count],"SAVE_EVERY");
  sprintf(parameter.name[++para_count],"PULS_TIME"); 
  sprintf(parameter.name[++para_count],"MIN_SHRINK_FACTOR");
  sprintf(parameter.name[++para_count],"INTERPOL_SPLAT");  
  sprintf(parameter.name[++para_count],"MAX_LOOPS");
  sprintf(parameter.name[++para_count],"N_POT_ARRAY");
  sprintf(parameter.name[++para_count],"USE_MAG_PA");   
  sprintf(parameter.name[++para_count],"ENABLE_EPOT");   
  sprintf(parameter.name[++para_count],"TRAP_THETA_START");
  sprintf(parameter.name[++para_count],"TRAP_THETA_STEP");
  sprintf(parameter.name[++para_count],"TRAP_ENERGY_START");
  sprintf(parameter.name[++para_count],"TRAP_ENERGY_END");
  sprintf(parameter.name[++para_count],"TRAP_PHI_START");
  sprintf(parameter.name[++para_count],"ENABLE_NEG_Y");
  sprintf(parameter.name[++para_count],"SET_Y_PLANE");
  sprintf(parameter.name[++para_count],"SET_Z_PLANE");
  sprintf(parameter.name[++para_count],"TRAP_MAX_MIRRORS");
  sprintf(parameter.name[++para_count],"TRAP_MASS");
  sprintf(parameter.name[++para_count],"TRAP_CHARGE");
  sprintf(parameter.name[++para_count],"TRAP_MAX_STEP_LENGTH");
  sprintf(parameter.name[++para_count],"TRAP_E_PARA_MIN");
  sprintf(parameter.name[++para_count],"TRAP_MAX_TOF");
  sprintf(parameter.name[++para_count],"TRAP_CALC_ORDER");
  sprintf(parameter.name[++para_count],"TRANS_B_PINCH");
  sprintf(parameter.name[++para_count],"TRANS_U_PINCH");
  sprintf(parameter.name[++para_count],"TRANS_STEPS");
  sprintf(parameter.name[++para_count],"TRAP_START_X");
  sprintf(parameter.name[++para_count],"TRAP_STOP_X");
  sprintf(parameter.name[++para_count],"TRAP_STEP_X");
  sprintf(parameter.name[++para_count],"TRAP_START_Y");
  sprintf(parameter.name[++para_count],"TRAP_STOP_Y");
  sprintf(parameter.name[++para_count],"TRAP_STEP_Y");
  sprintf(parameter.name[++para_count],"ENABLE_PARA_ENERGY_LOSS");
  sprintf(parameter.name[++para_count],"ENABLE_PERP_ENERGY_LOSS");
  sprintf(parameter.name[++para_count],"RESIDUAL_GAS_PRESSURE");
  sprintf(parameter.name[++para_count],"CALC_ORDER");
  sprintf(parameter.name[++para_count],"DIPOLE_VALUE");
  sprintf(parameter.name[++para_count],"MAX_MIRRORS");
  sprintf(parameter.name[++para_count],"MAX_TOF_IN_SEC");
  sprintf(parameter.name[++para_count],"E_PARA_MIN");
  sprintf(parameter.name[++para_count],"MAX_STEP_LENGTH");
  sprintf(parameter.name[++para_count],"E_MIN_COOLING");
  sprintf(parameter.name[++para_count],"REL_START_ANGLE");

  sprintf(parameter.name[++para_count],"B_FIELD_BEN1");
  sprintf(parameter.name[++para_count],"B_FIELD_BEN2");
  sprintf(parameter.name[++para_count],"RAD_CALC_MODE");
  sprintf(parameter.name[++para_count],"RAD_SHIFT");
  sprintf(parameter.name[++para_count],"ANTENNA_TEMP");
  sprintf(parameter.name[++para_count],"ANTENNA_POS");
  sprintf(parameter.name[++para_count],"IMPEDANCE");
  sprintf(parameter.name[++para_count],"RAD_ATTEN");

  sprintf(parameter.name[++para_count],"FFT_ON");
  sprintf(parameter.name[++para_count],"FFT_RESAMPLE_TSTEP");
  sprintf(parameter.name[++para_count],"FFT_MAX_NPTS");

  parameter.el_mm_unit = 0.;
  parameter.el_x_offset = 0.;
  parameter.el_y_offset = 0.;
  parameter.el_z_offset = 0.;
  parameter.mag_mm_unit = 0.;
  parameter.mag_x_offset = 0.;
  parameter.mag_y_offset = 0.;
  parameter.mag_z_offset = 0.;
  parameter.spec_in = 0.;
  parameter.spec_out = 0.;
  parameter.max_radius = 0.;
  parameter.save_every = 0;
  parameter.puls_time = 0.;
  parameter.min_shrink_factor = 0.;
  parameter.interpol_splat = 0;
  parameter.max_loops = 0;
  parameter.n_pot_array = 0;
  parameter.use_mag_pa = 0;
  parameter.enable_epot = 0;
  parameter.enable_rel_start_angle = 0;
  parameter.trap_theta_start = 0.;
  parameter.trap_theta_step = 0.;
  parameter.trap_energy_start = 0.;
  parameter.trap_energy_end = 0.;
  parameter.trap_phi_start = 0.;
  parameter.trap_neg_y = 0;
  parameter.trap_y_plane = 0;
  parameter.trap_z_plane = 0;
  parameter.trap_max_mirrors = 0;  
  parameter.trap_mass = 0.;      
  parameter.trap_charge = 0.;
  parameter.trap_max_step_length = 0.;
  parameter.trap_e_para_min = 0.;
  parameter.trap_max_tof = 0.;
  parameter.trap_calc_order = 0;  
  parameter.trans_b_pinch = 0;  
  parameter.trans_u_pinch = 0;  
  parameter.trans_steps = 0;  
  parameter.trap_start_x = 0.;
  parameter.trap_stop_x = 0.;
  parameter.trap_step_x = 0;
  parameter.trap_start_y = 0.;
  parameter.trap_stop_y = 0.;
  parameter.trap_step_y = 0.;
  parameter.para_e_loss = 0;
  parameter.perp_e_loss = 0;
  parameter.residual_gas_pressure = -1.;  // default is off
  parameter.calc_order = 0;
  parameter.dipole_value = 0.;
  parameter.max_mirrors = 0;
  parameter.max_tof_in_sec = 0.;
  parameter.e_para_min = 0.;
  parameter.max_step_length = 0.;
  parameter.e_min_cooling = 0.;

  parameter.b_field_ben1 = 10.0;
  parameter.b_field_ben2 = 1.0;
  parameter.rad_calc_mode = 2;
  parameter.rad_shift = 1;
  parameter.antenna_temp = 0;
  parameter.antenna_pos = +1;
  parameter.impedance = +1;
  parameter.rad_atten = 1;

  parameter.fft_on = 1;
  parameter.fft_resample_tstep = 1e-6;
  parameter.fft_max_npts = 10000000;

  printf("    DONE.\n");
}

/******************** store_parameter ***********************/
int store_parameter(char *identifier,double value)
{
  int status;
  int number = 0;

  if (strcmp(identifier,parameter.name[++number]) == 0)  
    {parameter.el_mm_unit = value; parameter.inited[number] = 1; status = 1;}
  else if (strcmp(identifier,parameter.name[++number]) == 0)   
    {parameter.el_x_offset = value; parameter.inited[number] = 1; status = 1; }
  else if (strcmp(identifier,parameter.name[++number]) == 0)   
    {parameter.el_y_offset = value; parameter.inited[number] = 1;  status = 1;}
  else if (strcmp(identifier,parameter.name[++number]) == 0)   
    {parameter.el_z_offset = value; parameter.inited[number] = 1;  status = 1;}
  else if (strcmp(identifier,parameter.name[++number]) == 0)   
    {parameter.mag_mm_unit = value; parameter.inited[number] = 1;  status = 1;}
  else if (strcmp(identifier,parameter.name[++number]) == 0)   
    {parameter.mag_x_offset = value; parameter.inited[number] = 1;  status = 1;}
  else if (strcmp(identifier,parameter.name[++number]) == 0)   
    {parameter.mag_y_offset = value; parameter.inited[number] = 1;  status = 1;}
  else if (strcmp(identifier,parameter.name[++number]) == 0)   
    {parameter.mag_z_offset = value; parameter.inited[number] = 1;  status = 1;}
  else if (strcmp(identifier,parameter.name[++number]) == 0)   
    {parameter.spec_in = value; parameter.inited[number] = 1;  status = 1;}
  else if (strcmp(identifier,parameter.name[++number]) == 0)   
    {parameter.spec_out = value; parameter.inited[number] = 1;  status = 1;}
  else if (strcmp(identifier,parameter.name[++number]) == 0)   
    {parameter.max_radius = value; parameter.inited[number] = 1;  status = 1;}
  else if (strcmp(identifier,parameter.name[++number]) == 0)   
    {parameter.save_every = (int)value; parameter.inited[number] = 1;  status = 1;}
  else if (strcmp(identifier,parameter.name[++number]) == 0)   
    {parameter.puls_time = value; parameter.inited[number] = 1;  status = 1;}
  else if (strcmp(identifier,parameter.name[++number]) == 0)   
    {parameter.min_shrink_factor = value; parameter.inited[number] = 1;  status = 1;}
  else if (strcmp(identifier,parameter.name[++number]) == 0)   
    {parameter.interpol_splat = (int)value; parameter.inited[number] = 1;  status = 1;}
  else if (strcmp(identifier,parameter.name[++number]) == 0)   
    {parameter.max_loops = (int)value; parameter.inited[number] = 1;  status = 1;}
  else if (strcmp(identifier,parameter.name[++number]) == 0)   
    {parameter.n_pot_array = (int)value; parameter.inited[number] = 1;  status = 1;}
  else if (strcmp(identifier,parameter.name[++number]) == 0)   
    {parameter.use_mag_pa = (int)value; parameter.inited[number] = 1;  status = 1;}
  else if (strcmp(identifier,parameter.name[++number]) == 0)   
    {parameter.enable_epot = (int)value; parameter.inited[number] = 1;  status = 1;}
  else if (strcmp(identifier,parameter.name[++number]) == 0)   
    {parameter.trap_theta_start = value; parameter.inited[number] = 1;  status = 1;}
  else if (strcmp(identifier,parameter.name[++number]) == 0)   
    {parameter.trap_theta_step = value; parameter.inited[number] = 1;  status = 1;}
  else if (strcmp(identifier,parameter.name[++number]) == 0)   
    {parameter.trap_energy_start = value; parameter.inited[number] = 1;  status = 1;}
  else if (strcmp(identifier,parameter.name[++number]) == 0)   
    {parameter.trap_energy_end = value; parameter.inited[number] = 1;  status = 1;}
  else if (strcmp(identifier,parameter.name[++number]) == 0)   
    {parameter.trap_phi_start = value; parameter.inited[number] = 1;  status = 1;}
  else if (strcmp(identifier,parameter.name[++number]) == 0)   
    {parameter.trap_neg_y = (int)value; parameter.inited[number] = 1;  status = 1;}
  else if (strcmp(identifier,parameter.name[++number]) == 0)   
    {parameter.trap_y_plane = (int)value; parameter.inited[number] = 1;  status = 1;}
  else if (strcmp(identifier,parameter.name[++number]) == 0)   
    {parameter.trap_z_plane = (int)value; parameter.inited[number] = 1;  status = 1;}
  else if (strcmp(identifier,parameter.name[++number]) == 0)   
    {parameter.trap_max_mirrors = (int)value; parameter.inited[number] = 1;  status = 1;}
  else if (strcmp(identifier,parameter.name[++number]) == 0)   
    {parameter.trap_mass = value; parameter.inited[number] = 1;  status = 1;}
  else if (strcmp(identifier,parameter.name[++number]) == 0)   
    {parameter.trap_charge = value; parameter.inited[number] = 1;  status = 1;}
  else if (strcmp(identifier,parameter.name[++number]) == 0)   
    {parameter.trap_max_step_length = value; parameter.inited[number] = 1;  status = 1;}
  else if (strcmp(identifier,parameter.name[++number]) == 0)   
    {parameter.trap_e_para_min = value; parameter.inited[number] = 1;  status = 1;}
  else if (strcmp(identifier,parameter.name[++number]) == 0)   
    {parameter.trap_max_tof = value; parameter.inited[number] = 1;  status = 1;}
  else if (strcmp(identifier,parameter.name[++number]) == 0)   
    {parameter.trap_calc_order = (int)value; parameter.inited[number] = 1;  status = 1;}
  else if (strcmp(identifier,parameter.name[++number]) == 0)   
    {parameter.trans_b_pinch = (int)value; parameter.inited[number] = 1;  status = 1;}
  else if (strcmp(identifier,parameter.name[++number]) == 0)   
    {parameter.trans_u_pinch = (int)value; parameter.inited[number] = 1;  status = 1;}
  else if (strcmp(identifier,parameter.name[++number]) == 0)   
    {parameter.trans_steps = (int)value; parameter.inited[number] = 1;  status = 1;}
  else if (strcmp(identifier,parameter.name[++number]) == 0)   
    {parameter.trap_start_x = value; parameter.inited[number] = 1;  status = 1;}
  else if (strcmp(identifier,parameter.name[++number]) == 0)   
    {parameter.trap_stop_x = value; parameter.inited[number] = 1;  status = 1;}
  else if (strcmp(identifier,parameter.name[++number]) == 0)   
    {parameter.trap_step_x = value; parameter.inited[number] = 1;  status = 1;}
  else if (strcmp(identifier,parameter.name[++number]) == 0)   
    {parameter.trap_start_y = value; parameter.inited[number] = 1;  status = 1;}
  else if (strcmp(identifier,parameter.name[++number]) == 0)   
    {parameter.trap_stop_y = value; parameter.inited[number] = 1;  status = 1;}
  else if (strcmp(identifier,parameter.name[++number]) == 0)   
    {parameter.trap_step_y = value; parameter.inited[number] = 1;  status = 1;}
  else if (strcmp(identifier,parameter.name[++number]) == 0)   
    {parameter.para_e_loss = (int)value; parameter.inited[number] = 1;  status = 1;}
  else if (strcmp(identifier,parameter.name[++number]) == 0)   
    {parameter.perp_e_loss = (int)value; parameter.inited[number] = 1;  status = 1;}
  else if (strcmp(identifier,parameter.name[++number]) == 0)   
    {parameter.residual_gas_pressure = value; parameter.inited[number] = 1;  status = 1;}
  else if (strcmp(identifier,parameter.name[++number]) == 0)   
    {parameter.calc_order = (int)value; parameter.inited[number] = 1;  status = 1;}
  else if (strcmp(identifier,parameter.name[++number]) == 0)   
    {parameter.dipole_value = value; parameter.inited[number] = 1;  status = 1;}
  else if (strcmp(identifier,parameter.name[++number]) == 0)   
    {parameter.max_mirrors = (int)value; parameter.inited[number] = 1;  status = 1;}
  else if (strcmp(identifier,parameter.name[++number]) == 0)   
    {parameter.max_tof_in_sec = value; parameter.inited[number] = 1;  status = 1;}
  else if (strcmp(identifier,parameter.name[++number]) == 0)   
    {parameter.e_para_min = value; parameter.inited[number] = 1;  status = 1;}
  else if (strcmp(identifier,parameter.name[++number]) == 0)   
    {parameter.max_step_length = value; parameter.inited[number] = 1;  status = 1;}
  else if (strcmp(identifier,parameter.name[++number]) == 0)   
    {parameter.e_min_cooling = value; parameter.inited[number] = 1;  status = 1;} 
  else if (strcmp(identifier,parameter.name[++number]) == 0)   
    {parameter.enable_rel_start_angle = (int)value; parameter.inited[number] = 1;  status = 1;}
  else if (strcmp(identifier,parameter.name[++number]) == 0)   
    {parameter.b_field_ben1 = value; parameter.inited[number] = 1;  status = 1;}
  else if (strcmp(identifier,parameter.name[++number]) == 0)   
    {parameter.b_field_ben2 = value; parameter.inited[number] = 1;  status = 1;}

  else if (strcmp(identifier,parameter.name[++number]) == 0)   
    {parameter.rad_calc_mode = (int)value; parameter.inited[number] = 1;  status = 1;}

  else if (strcmp(identifier,parameter.name[++number]) == 0)   
    {parameter.rad_shift = (int)value; parameter.inited[number] = 1;  status = 1;}
  
  else if (strcmp(identifier,parameter.name[++number]) == 0)   
    {parameter.antenna_temp = value; parameter.inited[number] = 1;  status = 1;}

  else if (strcmp(identifier,parameter.name[++number]) == 0)   
    {parameter.antenna_pos = value; parameter.inited[number] = 1;  status = 1;}

  else if (strcmp(identifier,parameter.name[++number]) == 0)   
    {parameter.impedance = value; parameter.inited[number] = 1;  status = 1;}

  else if (strcmp(identifier,parameter.name[++number]) == 0)   
    {parameter.rad_atten = (int)value; parameter.inited[number] = 1;  status = 1;}
  
  else if (strcmp(identifier,parameter.name[++number]) == 0)   
    {parameter.fft_on= (int)value; parameter.inited[number] = 1;  status = 1;}

  else if (strcmp(identifier,parameter.name[++number]) == 0)   
    {parameter.fft_resample_tstep = value; parameter.inited[number] = 1;  status = 1;}

  else if (strcmp(identifier,parameter.name[++number]) == 0)   
    {parameter.fft_max_npts = (int)value; parameter.inited[number] = 1;  status = 1;}

  else 
    {
      printf("\nERROR: Identifier %s unknown!\n\n",identifier);
      status = 0;
    }
  return status;
}

/******************** check_parameter ***********************/
int check_parameters()
{
  int status;
  int paracount;

  printf("\nCheck input Parameters.");
  status = 1;
  printf("\n");
  for (paracount = 1;paracount <= max_parameter; paracount++)
    {
      if (parameter.inited[paracount] == 0)
	{
	  printf("ERROR: Missing Parameter %s!\n",parameter.name[paracount]);
	  status = 0;
	}
    }
  if (status == 1) printf("                                      DONE.\n");
  else printf("\n");
  return status;
}


/******************** set_parameter ***********************/
int set_parameter(char* identifier, double value)
{
  int status = store_parameter(identifier,value);

  if (status == 1)
    {
      printf("set parameter: %-25s = %17.8f\n",identifier,value);
    }
  return status;
}



/******************** load_init_data ***********************/
int load_init_data(char* filename)
{
  // main procedure to load all init data


  char ini_filename[255];    // filename for ini data
  char identifier[32];      // container for itentifier
  char remark[256];         // container for remarks
  double value;             // container for value of itentifier
  int f_status;
  FILE* f_ini;               // parameter input file

  strcpy(ini_filename,filename);
  strcat(ini_filename,".ini"); // data input from file with extension .ini
  f_ini = fopen(ini_filename, "r");

  printf("Read ADIPARK 4 Parameter File: %s\n\n",ini_filename);

  if (f_ini == (FILE *)0)                 // disk error handling
    {
      fprintf(stderr, "\nERROR: Can't read file %s\n\n",ini_filename);
      f_status = 1; //exit(1);
    }
  else 
    {
      f_status = 0;
      while (f_status!=EOF)
	{
	  f_status = fscanf(f_ini,"#define %s %lf",identifier,&value);
	  //	  printf("status1: %d\n",f_status);
	  if (f_status == 2)
	    {
	      if (store_parameter(identifier,value) == 1)
		{
		  printf("input value: %-25s = %17.8f   DONE.\n",identifier,value);
		}

	      f_status = fscanf(f_ini,"%[^\n]\n",remark);
	      //	      printf("status2: %d\n",f_status);
	      if (f_status == 0)
		{
		  f_status = fscanf(f_ini,"\n");
		}
	      //	      printf("status3: %d\n",f_status);
	      //	      printf("remark: %s\n",remark);
	    }
	  else {
	    f_status = fscanf(f_ini,"%[^\n]\n",remark);
	    if (f_status == 0)
	      {
		f_status = fscanf(f_ini,"\n");
	      }
	  }
	  //	  printf("status4: %d\n",f_status);
	}
      fclose(f_ini);
    }
  return check_parameters();
}

