#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <fstream>
#include "parameter.hpp"

using namespace std;

/***********Creator & Destructor*****************/
parameter::parameter():
  // adipark settings
  adip_run_mode(1),
  use_mag_pa(2),
  // tracking module
  max_mirrors(5000),
  max_tof_in_sec(0.000003),
  max_loops(1000000),
  max_step_length(0.001),
  save_every(1),
  para_e_loss(1),
  perp_e_loss(1),
  calc_order(2),
  enable_rel_start_angle(0),
  puls_time(1),
  min_shrink_factor(0.000000001),
  interpol_splat(1),
  residual_gas_pressure(-1),
  e_min_cooling(1),
  dipole_value(0),
  e_para_min(0.001),
  spec_in(-420),
  spec_out(420),
  max_radius(475),
  // magnetic potential array
  mag_mm_unit(5),
  mag_x_offset(0),
  mag_y_offset(0),
  mag_z_offset(0),
  n_pot_array(30125000),
  // parabolic magnetic field
  b_field_ben1(10),
  b_field_ben2(1),
  // electric fields
  enable_epot(0),
  el_mm_unit(1),
  el_x_offset(0),
  el_y_offset(0),
  el_z_offset(0),
  // transmission module
  trans_b_pinch(10),
  trans_u_pinch(0),
  trans_steps(10),
  // trapping module
  trap_start_x(0),
  trap_stop_x(-200),
  trap_step_x(-5),
  trap_start_y(0),
  trap_stop_y(45),
  trap_step_y(5),
  trap_energy_start(32),
  trap_energy_end(8),
  trap_theta_start(80),
  trap_theta_step(10),
  trap_phi_start(0),
  trap_neg_y(0),
  trap_y_plane(1),
  trap_z_plane(2),
  trap_max_mirrors(250),  
  trap_mass(1),
  trap_charge(1),
  trap_max_step_length(0.01),
  trap_e_para_min(0),
  trap_calc_order(0),
  trap_max_tof(0.000003),
  // adi2fft settings
  rad_calc_mode(2),
  rad_shift(1),
  antenna_temp(0),
  impedance(1),
  rad_atten(1),
  antenna_pos(500),
  fft_on(1),
  fft_resample_tstep(0.000003),
  fft_max_npts(20000000),
  leave_antenna(0),
  filter_lo(26.9E9),
  filter_sf(6E8)
{
}

parameter::~parameter()
{
}

/***********General Use**************************/
void parameter::restore_defaults()
{
  // adipark settings
  adip_run_mode = 1;
  use_mag_pa = 2;
  // tracking module
  max_mirrors = 5000;
  max_tof_in_sec = 0.000003;
  max_loops = 1000000;
  max_step_length = 0.001;
  save_every = 1;
  para_e_loss = 1;
  perp_e_loss = 1;
  calc_order = 2;
  enable_rel_start_angle = 0;
  puls_time = 1;
  min_shrink_factor = 0.000000001;
  interpol_splat = 1;
  residual_gas_pressure = -1;
  e_min_cooling = 1;
  dipole_value = 0;
  e_para_min = 0.001;
  spec_in = -420;
  spec_out = 420;
  max_radius = 475;
  // magnetic potential array
  mag_mm_unit = 5;
  mag_x_offset = 0;
  mag_y_offset = 0;
  mag_z_offset = 0;
  n_pot_array = 30125000;
  // parabolic magnetic field
  b_field_ben1 = 10;
  b_field_ben2 = 1;
  // electric fields
  enable_epot = 0;
  el_mm_unit = 1;
  el_x_offset = 0;
  el_y_offset = 0;
  el_z_offset = 0;
  // transmission module
  trans_b_pinch = 10;
  trans_u_pinch = 0;
  trans_steps = 10;
  // trapping module
  trap_start_x = 0;
  trap_stop_x = -200;
  trap_step_x = -5;
  trap_start_y = 0;
  trap_stop_y = 45;
  trap_step_y = 5;
  trap_energy_start = 32;
  trap_energy_end = 8;
  trap_theta_start = 80;
  trap_theta_step = 10;
  trap_phi_start = 0;
  trap_neg_y = 0;
  trap_y_plane = 1;
  trap_z_plane = 2;
  trap_max_mirrors = 250;  
  trap_mass = 1;
  trap_charge = 1;
  trap_max_step_length = 0.01;
  trap_e_para_min = 0;
  trap_calc_order = 0;
  trap_max_tof = 0.000003;
  // adi2fft settings
  rad_calc_mode = 2;
  rad_shift = 1;
  antenna_temp = 0;
  impedance = 1;
  rad_atten = 1;
  antenna_pos = 500;
  fft_on = 1;
  fft_resample_tstep = 0.000003;
  fft_max_npts = 20000000;
  leave_antenna = 0;
  filter_lo = 26.9E9;
  filter_sf = 6E8;
}

void parameter::parse_file(char inifile[255])
{
  FILE *infile;
  infile = fopen(inifile, "r");
  int stat;
  char identifier[32];
  char dummy[255];
  double value;
  stat = 0;
  while (stat != EOF) {
    stat = fscanf(infile, "#define %s %lf", identifier, &value);
    if (stat == 2) {
      if (strcmp(identifier,"ADIP_RUN_MODE") == 0) {
        set_adip_run_mode(int(value));
      } else if (strcmp(identifier,"USE_MAG_PA") == 0) {
        set_use_mag_pa(int(value));
      } else if (strcmp(identifier,"MAX_MIRRORS") == 0) {
        set_max_mirrors(int(value));
      } else if (strcmp(identifier,"MAX_TOF_IN_SEC") == 0) {
        set_max_tof_in_sec(double(value));
      } else if (strcmp(identifier,"MAX_LOOPS") == 0) {
        set_max_loops(int(value));
      } else if (strcmp(identifier,"MAX_STEP_LENGTH") == 0) {
        set_max_step_length(double(value));
      } else if (strcmp(identifier,"SAVE_EVERY") == 0) {
        set_save_every(int(value));
      } else if (strcmp(identifier,"ENABLE_PERP_ENERGY_LOSS") == 0) {
        set_perp_e_loss(int(value));
      } else if (strcmp(identifier,"ENABLE_PARA_ENERGY_LOSS") == 0) {
        set_para_e_loss(int(value));
      } else if (strcmp(identifier,"CALC_ORDER") == 0) {
        set_calc_order(int(value));
      } else if (strcmp(identifier,"REL_START_ANGLE") == 0) {
        set_rel_start_angle(int(value));
      } else if (strcmp(identifier,"PULS_TIME") == 0) {
        set_puls_time(double(value));
      } else if (strcmp(identifier,"MIN_SHRINK_FACTOR") == 0) {
        set_min_shrink_factor(double(value));
      } else if (strcmp(identifier,"INTERPOL_SPLAT") == 0) {
        set_interpol_splat(int(value));
      } else if (strcmp(identifier,"RESIDUAL_GAS_PRESSURE") == 0) {
        set_residual_gas_pressure(double(value));
      } else if (strcmp(identifier,"E_MIN_COOLING") == 0) {
        set_e_min_cooling(double(value));
      } else if (strcmp(identifier,"DIPOLE_VALUE") == 0) {
        set_dipole_value(double(value));
      } else if (strcmp(identifier,"E_PARA_MIN") == 0) {
        set_e_para_min(double(value));
      } else if (strcmp(identifier,"SPEC_IN") == 0) {
        set_spec_in(double(value));
      } else if (strcmp(identifier,"SPEC_OUT") == 0) {
        set_spec_out(double(value));
      } else if (strcmp(identifier,"MAX_RADIUS") == 0) {
        set_max_radius(double(value));
      } else if (strcmp(identifier,"MAG_MM_PER_UNIT") == 0) {
        set_mag_mm_unit(double(value));
      } else if (strcmp(identifier,"MAG_X_OFFSET_IN_CM") == 0) {
        set_mag_x_offset(double(value));
      } else if (strcmp(identifier,"MAG_Y_OFFSET_IN_CM") == 0) {
        set_mag_y_offset(double(value));
      } else if (strcmp(identifier,"MAG_Z_OFFSET_IN_CM") == 0) {
        set_mag_z_offset(double(value));
      } else if (strcmp(identifier,"N_POT_ARRAY") == 0) {
        set_n_pot_array(int(value));
      } else if (strcmp(identifier,"B_FIELD_BEN1") == 0) {
        set_b_field_ben1(double(value));
      } else if (strcmp(identifier,"B_FIELD_BEN2") == 0) {
        set_b_field_ben2(double(value));
      } else if (strcmp(identifier,"ENABLE_EPOT") == 0) {
        set_enable_epot(int(value));
      } else if (strcmp(identifier,"EL_MM_PER_UNIT") == 0) {
        set_el_mm_unit(double(value));
      } else if (strcmp(identifier,"EL_X_OFFSET_IN_CM") == 0) {
        set_el_x_offset(double(value));
      } else if (strcmp(identifier,"EL_Y_OFFSET_IN_CM") == 0) {
        set_el_y_offset(double(value));
      } else if (strcmp(identifier,"EL_Z_OFFSET_IN_CM") == 0) {
        set_el_z_offset(double(value));
      } else if (strcmp(identifier,"TRANS_B_PINCH") == 0) {
        set_trans_b_pinch(double(value));
      } else if (strcmp(identifier,"TRANS_U_PINCH") == 0) {
        set_trans_u_pinch(double(value));
      } else if (strcmp(identifier,"TRANS_STEPS") == 0) {
        set_trans_steps(int(value));
      } else if (strcmp(identifier,"TRAP_START_X") == 0) {
        set_trap_start_x(double(value));
      } else if (strcmp(identifier,"TRAP_STOP_X") == 0) {
        set_trap_stop_x(double(value));
      } else if (strcmp(identifier,"TRAP_STEP_X") == 0) {
        set_trap_step_x(double(value));
      } else if (strcmp(identifier,"TRAP_START_Y") == 0) {
        set_trap_start_y(double(value));
      } else if (strcmp(identifier,"TRAP_STOP_Y") == 0) {
        set_trap_stop_y(double(value));
      } else if (strcmp(identifier,"TRAP_STEP_Y") == 0) {
        set_trap_step_y(double(value));
      } else if (strcmp(identifier,"TRAP_ENERGY_START") == 0) {
        set_trap_energy_start(double(value));
      } else if (strcmp(identifier,"TRAP_ENERGY_END") == 0) {
        set_trap_energy_end(double(value));
      } else if (strcmp(identifier,"TRAP_THETA_START") == 0) {
        set_trap_theta_start(double(value));
      } else if (strcmp(identifier,"TRAP_THETA_STEP") == 0) {
        set_trap_theta_step(double(value));
      } else if (strcmp(identifier,"TRAP_PHI_START") == 0) {
        set_trap_phi_start(double(value));
      } else if (strcmp(identifier,"ENABLE_NEG_Y") == 0) {
        set_trap_neg_y(int(value));
      } else if (strcmp(identifier,"SET_Y_PLANE") == 0) {
        set_trap_y_plane(int(value));
      } else if (strcmp(identifier,"SET_Z_PLANE") == 0) {
        set_trap_z_plane(int(value));
      } else if (strcmp(identifier,"TRAP_MAX_MIRRORS") == 0) {
        set_trap_max_mirrors(int(value));
      } else if (strcmp(identifier,"TRAP_MASS") == 0) {
        set_trap_mass(double(value));
      } else if (strcmp(identifier,"TRAP_CHARGE") == 0) {
        set_trap_charge(double(value));
      } else if (strcmp(identifier,"TRAP_MAX_STEP_LENGTH") == 0) {
        set_trap_max_step_length(double(value));
      } else if (strcmp(identifier,"TRAP_E_PARA_MIN") == 0) {
        set_trap_e_para_min(double(value));
      } else if (strcmp(identifier,"TRAP_CALC_ORDER") == 0) {
        set_trap_calc_order(int(value));
      } else if (strcmp(identifier,"TRAP_MAX_TOF") == 0) {
        set_trap_max_tof(double(value));
      } else if (strcmp(identifier,"RAD_CALC_MODE") == 0) {
        set_rad_calc_mode(int(value));
      } else if (strcmp(identifier,"RAD_SHIFT") == 0) {
        set_rad_shift(int(value));
      } else if (strcmp(identifier,"ANTENNA_TEMP") == 0) {
        set_antenna_temp(double(value));
      } else if (strcmp(identifier,"IMPEDANCE") == 0) {
        set_impedance(double(value));
      } else if (strcmp(identifier,"RAD_ATTEN") == 0) {
        set_rad_atten(int(value));
      } else if (strcmp(identifier,"ANTENNA_POS") == 0) {
        set_antenna_pos(double(value));
      } else if (strcmp(identifier,"FFT_ON") == 0) {
        set_fft_on(int(value));
      } else if (strcmp(identifier,"FFT_RESAMPLE_TSTEP") == 0) {
        set_fft_resample_tstep(double(value));
      } else if (strcmp(identifier,"FFT_MAX_NPTS") == 0) {
        set_fft_max_npts(int(value));
      } else if (strcmp(identifier,"LEAVE_ANTENNA") == 0) {
        set_leave_antenna(int(value));
      } else if (strcmp(identifier,"FILTER_LO") == 0) {
        set_filter_lo(double(value));
      } else if (strcmp(identifier,"FILTER_SF") == 0) {
        set_filter_sf(double(value));
      } else {
        cout << "identifier not found" << endl;
        cin.ignore(1);
      }
    }
    stat = fscanf(infile, "%[^\n]\n", dummy);
    if (stat == 0) {
      stat = fscanf(infile, "\n");
    } else {
      stat = fscanf(infile, "%[^\n]\n", dummy);
      if (stat == 0) {
        stat = fscanf(infile, "\n");
      }
    }
  }
  fclose(infile);
}

/***********Get & Set****************************/
  // adipark settings*
int parameter::get_adip_run_mode()
/************************************************/
{
  return adip_run_mode;
}
void parameter::set_adip_run_mode()
{
  adip_run_mode = 1;
}
void parameter::set_adip_run_mode(int value)
{
  adip_run_mode = value;
}

/************************************************/
int parameter::get_use_mag_pa()
{
  return use_mag_pa;
}
void parameter::set_use_mag_pa()
{
  use_mag_pa = 2;
}
void parameter::set_use_mag_pa(int value)
{
  use_mag_pa = value;
}
  // tracking module
/************************************************/
int parameter::get_max_mirrors()
{
  return max_mirrors;
}
void parameter::set_max_mirrors()
{
  max_mirrors = 5000;
}
void parameter::set_max_mirrors(int value)
{
  max_mirrors = value;
}

/************************************************/
double parameter::get_max_tof_in_sec()
{
  return max_tof_in_sec;
}
void parameter::set_max_tof_in_sec()
{
  max_tof_in_sec = 0.000003;
}
void parameter::set_max_tof_in_sec(double value)
{
  max_tof_in_sec = value;
}

/************************************************/
int parameter::get_max_loops()
{
  return max_loops;
}
void parameter::set_max_loops()
{
  max_loops = 1000000;
}
void parameter::set_max_loops(int value)
{
  max_loops = value;
}

/************************************************/
double parameter::get_max_step_length()
{
  return max_step_length;
}
void parameter::set_max_step_length()
{
  max_step_length = 0.001;
}
void parameter::set_max_step_length(double value)
{
  max_step_length = value;
}

/************************************************/
int parameter::get_save_every()
{
  return save_every;
}
void parameter::set_save_every()
{
  save_every = 1;
}
void parameter::set_save_every(int value)
{
  save_every = value;
}

/************************************************/
int parameter::get_para_e_loss()
{
  return para_e_loss;
}
void parameter::set_para_e_loss()
{
  para_e_loss = 1;
}
void parameter::set_para_e_loss(int value)
{
  para_e_loss = value;
}

/************************************************/
int parameter::get_perp_e_loss()
{
  return perp_e_loss;
}
void parameter::set_perp_e_loss()
{
  perp_e_loss = 1;
}
void parameter::set_perp_e_loss(int value)
{
  perp_e_loss = value;
}

/************************************************/
int parameter::get_calc_order()
{
  return calc_order;
}
void parameter::set_calc_order()
{
  calc_order = 2;
}
void parameter::set_calc_order(int value)
{
  calc_order = value;
}

/************************************************/
int parameter::get_rel_start_angle()
{
  return enable_rel_start_angle;
}
void parameter::set_rel_start_angle()
{
  enable_rel_start_angle = 0;
}
void parameter::set_rel_start_angle(int value)
{
  enable_rel_start_angle = value;
}

/************************************************/
double parameter::get_puls_time()
{
  return puls_time;
}
void parameter::set_puls_time()
{
  puls_time = 1;
}
void parameter::set_puls_time(double value)
{
  puls_time = value;
}

/************************************************/
double parameter::get_min_shrink_factor()
{
  return min_shrink_factor;
}
void parameter::set_min_shrink_factor()
{
  min_shrink_factor = 0.000000001;
}
void parameter::set_min_shrink_factor(double value)
{
  min_shrink_factor = value;
}

/************************************************/
int parameter::get_interpol_splat()
{
  return interpol_splat;
}
void parameter::set_interpol_splat()
{
  interpol_splat = 1;
}
void parameter::set_interpol_splat(int value)
{
  interpol_splat = value;
}

/************************************************/
double parameter::get_residual_gas_pressure()
{
  return residual_gas_pressure;
}
void parameter::set_residual_gas_pressure()
{
  residual_gas_pressure = -1;
}
void parameter::set_residual_gas_pressure(double value)
{
  residual_gas_pressure = value;
}

/************************************************/
double parameter::get_e_min_cooling()
{
  return e_min_cooling;
}
void parameter::set_e_min_cooling()
{
  e_min_cooling = 1;
}
void parameter::set_e_min_cooling(double value)
{
  e_min_cooling = value;
}

/************************************************/
double parameter::get_dipole_value()
{
  return dipole_value;
}
void parameter::set_dipole_value()
{
  dipole_value = 0;
}
void parameter::set_dipole_value(double value)
{
  dipole_value = value;
}

/************************************************/
double parameter::get_e_para_min()
{
  return e_para_min;
}
void parameter::set_e_para_min()
{
  e_para_min = 0.001;
}
void parameter::set_e_para_min(double value)
{
  e_para_min = value;
}

/************************************************/
double parameter::get_spec_in()
{
  return spec_in;
}
void parameter::set_spec_in()
{
  spec_in = -420;
}
void parameter::set_spec_in(double value)
{
  spec_in = value;
}

/************************************************/
double parameter::get_spec_out()
{
  return spec_out;
}
void parameter::set_spec_out()
{
  spec_out = 420;
}
void parameter::set_spec_out(double value)
{
  spec_out = value;
}

/************************************************/
double parameter::get_max_radius()
{
  return max_radius;
}
void parameter::set_max_radius()
{
  max_radius = 475;
}
void parameter::set_max_radius(double value)
{
  max_radius = value;
}

  // magnetic potential array
/************************************************/
double parameter::get_mag_mm_unit()
{
  return mag_mm_unit;
}
void parameter::set_mag_mm_unit()
{
  mag_mm_unit = 5;
}
void parameter::set_mag_mm_unit(double value)
{
  mag_mm_unit = value;
}

/************************************************/
double parameter::get_mag_x_offset()
{
  return mag_x_offset;
}
void parameter::set_mag_x_offset()
{
  mag_x_offset = 0;
}
void parameter::set_mag_x_offset(double value)
{
  mag_x_offset = value;
}

/************************************************/
double parameter::get_mag_y_offset()
{
  return mag_y_offset;
}
void parameter::set_mag_y_offset()
{
  mag_y_offset = 0;
}
void parameter::set_mag_y_offset(double value)
{
  mag_y_offset = value;
}

/************************************************/
double parameter::get_mag_z_offset()
{
  return mag_z_offset;
}
void parameter::set_mag_z_offset()
{
  mag_z_offset = 0;
}
void parameter::set_mag_z_offset(double value)
{
  mag_z_offset = value;
}

/************************************************/
int parameter::get_n_pot_array()
{
  return n_pot_array;
}
void parameter::set_n_pot_array()
{
  n_pot_array = 30125000;
}
void parameter::set_n_pot_array(int value)
{
  n_pot_array = value;
}

  // parabolic magnetic field
/************************************************/
double parameter::get_b_field_ben1()
{
  return b_field_ben1;
}
void parameter::set_b_field_ben1()
{
  b_field_ben1 = 10;
}
void parameter::set_b_field_ben1(double value)
{
  b_field_ben1 = value;
}

/************************************************/
double parameter::get_b_field_ben2()
{
  return b_field_ben2;
}
void parameter::set_b_field_ben2()
{
  b_field_ben2 = 1;
}
void parameter::set_b_field_ben2(double value)
{
  b_field_ben2 = value;
}

  // electric fields
/************************************************/
int parameter::get_enable_epot()
{
  return enable_epot;
}
void parameter::set_enable_epot()
{
  enable_epot = 0;
}
void parameter::set_enable_epot(int value)
{
  enable_epot = value;
}

/************************************************/
double parameter::get_el_mm_unit()
{
  return el_mm_unit;
}
void parameter::set_el_mm_unit()
{
  el_mm_unit = 1;
}
void parameter::set_el_mm_unit(double value)
{
  el_mm_unit = value;
}

/************************************************/
double parameter::get_el_x_offset()
{
  return el_x_offset;
}
void parameter::set_el_x_offset()
{
  el_x_offset = 0;
}
void parameter::set_el_x_offset(double value)
{
  el_x_offset = value;
}

/************************************************/
double parameter::get_el_y_offset()
{
  return el_y_offset;
}
void parameter::set_el_y_offset()
{
  el_y_offset = 0;
}
void parameter::set_el_y_offset(double value)
{
  el_y_offset = value;
}

/************************************************/
double parameter::get_el_z_offset()
{
  return el_z_offset;
}
void parameter::set_el_z_offset()
{
  el_z_offset = 0;
}
void parameter::set_el_z_offset(double value)
{
  el_z_offset = value;
}

  // transmission module
/************************************************/
double parameter::get_trans_b_pinch()
{
  return trans_b_pinch;
}
void parameter::set_trans_b_pinch()
{
  trans_b_pinch = 10;
}
void parameter::set_trans_b_pinch(double value)
{
  trans_b_pinch = value;
}

/************************************************/
double parameter::get_trans_u_pinch()
{
  return trans_u_pinch;
}
void parameter::set_trans_u_pinch()
{
  trans_u_pinch = 0;
}
void parameter::set_trans_u_pinch(double value)
{
  trans_u_pinch = value;
}

/************************************************/
int parameter::get_trans_steps()
{
  return trans_steps;
}
void parameter::set_trans_steps()
{
  trans_steps = 10;
}
void parameter::set_trans_steps(int value)
{
  trans_steps = value;
}

  // trapping module
/************************************************/
double parameter::get_trap_start_x()
{
  return trap_start_x;
}
void parameter::set_trap_start_x()
{
  trap_start_x = 0;
}
void parameter::set_trap_start_x(double value)
{
  trap_start_x = value;
}

/************************************************/
double parameter::get_trap_stop_x()
{
  return trap_stop_x;
}
void parameter::set_trap_stop_x()
{
  trap_stop_x = -200;
}
void parameter::set_trap_stop_x(double value)
{
  trap_stop_x = value;
}

/************************************************/
double parameter::get_trap_step_x()
{
  return trap_step_x;
}
void parameter::set_trap_step_x()
{
  trap_step_x = -5;
}
void parameter::set_trap_step_x(double value)
{
  trap_step_x = value;
}

/************************************************/
double parameter::get_trap_start_y()
{
  return trap_start_y;
}
void parameter::set_trap_start_y()
{
  trap_start_y = 0;
}
void parameter::set_trap_start_y(double value)
{
  trap_start_y = value;
}

/************************************************/
double parameter::get_trap_stop_y()
{
  return trap_stop_y;
}
void parameter::set_trap_stop_y()
{
  trap_stop_y = 45;
}
void parameter::set_trap_stop_y(double value)
{
  trap_stop_y = value;
}

/************************************************/
double parameter::get_trap_step_y()
{
  return trap_step_y;
}
void parameter::set_trap_step_y()
{
  trap_step_y = 5;
}
void parameter::set_trap_step_y(double value)
{
  trap_step_y = value;
}

/************************************************/
double parameter::get_trap_energy_start()
{
  return trap_energy_start;
}
void parameter::set_trap_energy_start()
{
  trap_energy_start = 32;
}
void parameter::set_trap_energy_start(double value)
{
  trap_energy_start = value;
}

/************************************************/
double parameter::get_trap_energy_end()
{
  return trap_energy_end;
}
void parameter::set_trap_energy_end()
{
  trap_energy_end = 8;
}
void parameter::set_trap_energy_end(double value)
{
  trap_energy_end = value;
}

/************************************************/
double parameter::get_trap_theta_start()
{
  return trap_theta_start;
}
void parameter::set_trap_theta_start()
{
  trap_theta_start = 80;
}
void parameter::set_trap_theta_start(double value)
{
  trap_theta_start = value;
}

/************************************************/
double parameter::get_trap_theta_step()
{
  return trap_theta_step;
}
void parameter::set_trap_theta_step()
{
  trap_theta_step = 10;
}
void parameter::set_trap_theta_step(double value)
{
  trap_theta_step = value;
}

/************************************************/
double parameter::get_trap_phi_start()
{
  return trap_phi_start;
}
void parameter::set_trap_phi_start()
{
  trap_phi_start = 0;
}
void parameter::set_trap_phi_start(double value)
{
  trap_phi_start = value;
}

/************************************************/
int parameter::get_trap_neg_y()
{
  return trap_neg_y;
}
void parameter::set_trap_neg_y()
{
  trap_neg_y = 0;
}
void parameter::set_trap_neg_y(int value)
{
  trap_neg_y = value;
}

/************************************************/
int parameter::get_trap_y_plane()
{
  return trap_y_plane;
}
void parameter::set_trap_y_plane()
{
  trap_y_plane = 1;
}
void parameter::set_trap_y_plane(int value)
{
  trap_y_plane = value;
}

/************************************************/
int parameter::get_trap_z_plane()
{
  return trap_z_plane;
}
void parameter::set_trap_z_plane()
{
  trap_z_plane = 2;
}
void parameter::set_trap_z_plane(int value)
{
  trap_z_plane = value;
}

/************************************************/
int parameter::get_trap_max_mirrors()
{
  return trap_max_mirrors;
}
void parameter::set_trap_max_mirrors()
{
  trap_max_mirrors = 250;  
}
void parameter::set_trap_max_mirrors(int value)
{
  trap_max_mirrors = value;
}

/************************************************/
double parameter::get_trap_mass()
{
  return trap_mass;
}
void parameter::set_trap_mass()
{
  trap_mass = 1;
}
void parameter::set_trap_mass(double value)
{
  trap_mass = value;
}

/************************************************/
double parameter::get_trap_charge()
{
  return trap_charge;
}
void parameter::set_trap_charge()
{
  trap_charge = 1;
}
void parameter::set_trap_charge(double value)
{
  trap_charge = value;
}

/************************************************/
double parameter::get_trap_max_step_length()
{
  return trap_max_step_length;
}
void parameter::set_trap_max_step_length()
{
  trap_max_step_length = 0.01;
}
void parameter::set_trap_max_step_length(double value)
{
  trap_max_step_length = value;
}

/************************************************/
double parameter::get_trap_e_para_min()
{
  return trap_e_para_min;
}
void parameter::set_trap_e_para_min()
{
  trap_e_para_min = 0;
}
void parameter::set_trap_e_para_min(double value)
{
  trap_e_para_min = value;
}

/************************************************/
int parameter::get_trap_calc_order()
{
  return trap_calc_order;
}
void parameter::set_trap_calc_order()
{
  trap_calc_order = 0;
}
void parameter::set_trap_calc_order(int value)
{
  trap_calc_order = value;
}

/************************************************/
double parameter::get_trap_max_tof()
{
  return trap_max_tof;
}
void parameter::set_trap_max_tof()
{
  trap_max_tof = 0.000003;
}
void parameter::set_trap_max_tof(double value)
{
  trap_max_tof = value;
}

  // adi2fft settings
/************************************************/
int parameter::get_rad_calc_mode()
{
  return rad_calc_mode;
}
void parameter::set_rad_calc_mode()
{
  rad_calc_mode = 2;
}
void parameter::set_rad_calc_mode(int value)
{
  rad_calc_mode = value;
}

/************************************************/
int parameter::get_rad_shift()
{
  return rad_shift;
}
void parameter::set_rad_shift()
{
  rad_shift = 1;
}
void parameter::set_rad_shift(int value)
{
  rad_shift = value;
}

/************************************************/
double parameter::get_antenna_temp()
{
  return antenna_temp;
}
void parameter::set_antenna_temp()
{
  antenna_temp = 0;
}
void parameter::set_antenna_temp(double value)
{
  antenna_temp = value;
}

/************************************************/
double parameter::get_impedance()
{
  return impedance;
}
void parameter::set_impedance()
{
  impedance = 1;
}
void parameter::set_impedance(double value)
{
  impedance = value;
}

/************************************************/
int parameter::get_rad_atten()
{
  return rad_atten;
}
void parameter::set_rad_atten()
{
  rad_atten = 1;
}
void parameter::set_rad_atten(int value)
{
  rad_atten = value;
}

/************************************************/
double parameter::get_antenna_pos()
{
  return antenna_pos;
}
void parameter::set_antenna_pos()
{
  antenna_pos = 500;
}
void parameter::set_antenna_pos(double value)
{
  antenna_pos = value;
}

/************************************************/
int parameter::get_fft_on()
{
  return fft_on;
}
void parameter::set_fft_on()
{
  fft_on = 1;
}
void parameter::set_fft_on(int value)
{
  fft_on = value;
}

/************************************************/
double parameter::get_fft_resample_tstep()
{
  return fft_resample_tstep;
}
void parameter::set_fft_resample_tstep()
{
  fft_resample_tstep = 0.000003;
}
void parameter::set_fft_resample_tstep(double value)
{
  fft_resample_tstep = value;
}

/************************************************/
int parameter::get_fft_max_npts()
{
  return fft_max_npts;
}
void parameter::set_fft_max_npts()
{
  fft_max_npts = 20000000;
}
void parameter::set_fft_max_npts(int value)
{
  fft_max_npts = value;
}

/************************************************/
int parameter::get_leave_antenna()
{
  return leave_antenna;
}
void parameter::set_leave_antenna()
{
  leave_antenna = 0;
}
void parameter::set_leave_antenna(int value)
{
  leave_antenna = value;
}
