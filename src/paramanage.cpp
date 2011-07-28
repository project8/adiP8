#include "paramanage.hpp"

using namespace std;

/***********Creator & Destructor*****************/
paramanage::paramanage()
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
}

paramanage::~paramanage()
{
}

/***********General Use**************************/
void paramanage::restore_defaults()
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
}

/***********Get & Set****************************/
  // adipark settings*
int paramanage::get_adip_run_mode()
/************************************************/
{
  return adip_run_mode;
}
void paramanage::set_adip_run_mode()
{
  adip_run_mode = 1;
}
void paramanage::set_adip_run_mode(int value)
{
  adip_run_mode = value;
}

/************************************************/
int paramanage::get_use_mag_pa()
{
  return use_mag_pa;
}
void paramanage::set_use_mag_pa()
{
  use_mag_pa = 2;
}
void paramanage::set_use_mag_pa(int value)
{
  use_mag_pa = value;
}
  // tracking module
/************************************************/
int paramanage::get_max_mirrors()
{
  return max_mirrors;
}
void paramanage::set_max_mirrors()
{
  max_mirrors = 5000;
}
void paramanage::set_max_mirrors(int value)
{
  max_mirrors = value;
}

/************************************************/
double paramanage::get_max_tof_in_sec()
{
  return max_tof_in_sec;
}
void paramanage::set_max_tof_in_sec()
{
  max_tof_in_sec = 0.000003;
}
void paramanage::set_max_tof_in_sec(double value)
{
  max_tof_in_sec = value;
}

/************************************************/
int paramanage::get_max_loops()
{
  return max_loops;
}
void paramanage::set_max_loops()
{
  max_loops = 1000000;
}
void paramanage::set_max_loops(int value)
{
  max_loops = value;
}

/************************************************/
double paramanage::get_max_step_length()
{
  return max_step_length;
}
void paramanage::set_max_step_length()
{
  max_step_length = 0.001;
}
void paramanage::set_max_step_length(double value)
{
  max_step_length = value;
}

/************************************************/
int paramanage::get_save_every()
{
  return save_every;
}
void paramanage::set_save_every()
{
  save_every = 1;
}
void paramanage::set_save_every(int value)
{
  save_every = value;
}

/************************************************/
int paramanage::get_para_e_loss()
{
  return para_e_loss;
}
void paramanage::set_para_e_loss()
{
  para_e_loss = 1;
}
void paramanage::set_para_e_loss(int value)
{
  para_e_loss = value;
}

/************************************************/
int paramanage::get_perp_e_loss()
{
  return perp_e_loss;
}
void paramanage::set_perp_e_loss()
{
  perp_e_loss = 1;
}
void paramanage::set_perp_e_loss(int value)
{
  perp_e_loss = value;
}

/************************************************/
int paramanage::get_calc_order()
{
  return calc_order;
}
void paramanage::set_calc_order()
{
  calc_order = 2;
}
void paramanage::set_calc_order(int value)
{
  calc_order = value;
}

/************************************************/
int paramanage::get_rel_start_angle()
{
  return enable_rel_start_angle;
}
void paramanage::set_rel_start_angle()
{
  enable_rel_start_angle = 0;
}
void paramanage::set_rel_start_angle(int value)
{
  enable_rel_start_angle = value;
}

/************************************************/
double paramanage::get_puls_time()
{
  return puls_time;
}
void paramanage::set_puls_time()
{
  puls_time = 1;
}
void paramanage::set_puls_time(double value)
{
  puls_time = value;
}

/************************************************/
double paramanage::get_min_shrink_factor()
{
  return min_shrink_factor;
}
void paramanage::set_min_shrink_factor()
{
  min_shrink_factor = 0.000000001;
}
void paramanage::set_min_shrink_factor(double value)
{
  min_shrink_factor = value;
}

/************************************************/
int paramanage::get_interpol_splat()
{
  return interpol_splat;
}
void paramanage::set_interpol_splat()
{
  interpol_splat = 1;
}
void paramanage::set_interpol_splat(int value)
{
  interpol_splat = value;
}

/************************************************/
double paramanage::get_residual_gas_pressure()
{
  return residual_gas_pressure;
}
void paramanage::set_residual_gas_pressure()
{
  residual_gas_pressure = -1;
}
void paramanage::set_residual_gas_pressure(double value)
{
  residual_gas_pressure = value;
}

/************************************************/
double paramanage::get_e_min_cooling()
{
  return e_min_cooling;
}
void paramanage::set_e_min_cooling()
{
  e_min_cooling = 1;
}
void paramanage::set_e_min_cooling(double value)
{
  e_min_cooling = value;
}

/************************************************/
double paramanage::get_dipole_value()
{
  return dipole_value;
}
void paramanage::set_dipole_value()
{
  dipole_value = 0;
}
void paramanage::set_dipole_value(double value)
{
  dipole_value = value;
}

/************************************************/
double paramanage::get_e_para_min()
{
  return e_para_min;
}
void paramanage::set_e_para_min()
{
  e_para_min = 0.001;
}
void paramanage::set_e_para_min(double value)
{
  e_para_min = value;
}

/************************************************/
double paramanage::get_spec_in()
{
  return spec_in;
}
void paramanage::set_spec_in()
{
  spec_in = -420;
}
void paramanage::set_spec_in(double value)
{
  spec_in = value;
}

/************************************************/
double paramanage::get_spec_out()
{
  return spec_out;
}
void paramanage::set_spec_out()
{
  spec_out = 420;
}
void paramanage::set_spec_out(double value)
{
  spec_out = value;
}

/************************************************/
double paramanage::get_max_radius()
{
  return max_radius;
}
void paramanage::set_max_radius()
{
  max_radius = 475;
}
void paramanage::set_max_radius(double value)
{
  max_radius = value;
}

  // magnetic potential array
/************************************************/
double paramanage::get_mag_mm_unit()
{
  return mag_mm_unit;
}
void paramanage::set_mag_mm_unit()
{
  mag_mm_unit = 5;
}
void paramanage::set_mag_mm_unit(double value)
{
  mag_mm_unit = value;
}

/************************************************/
double paramanage::get_mag_x_offset()
{
  return mag_x_offset;
}
void paramanage::set_mag_x_offset()
{
  mag_x_offset = 0;
}
void paramanage::set_mag_x_offset(double value)
{
  mag_x_offset = value;
}

/************************************************/
double paramanage::get_mag_y_offset()
{
  return mag_y_offset;
}
void paramanage::set_mag_y_offset()
{
  mag_y_offset = 0;
}
void paramanage::set_mag_y_offset(double value)
{
  mag_y_offset = value;
}

/************************************************/
double paramanage::get_mag_z_offset()
{
  return mag_z_offset;
}
void paramanage::set_mag_z_offset()
{
  mag_z_offset = 0;
}
void paramanage::set_mag_z_offset(double value)
{
  mag_z_offset = value;
}

/************************************************/
int paramanage::get_n_pot_array()
{
  return n_pot_array;
}
void paramanage::set_n_pot_array()
{
  n_pot_array = 30125000;
}
void paramanage::set_n_pot_array(int value)
{
  n_pot_array = value;
}

  // parabolic magnetic field
/************************************************/
double paramanage::get_b_field_ben1()
{
  return b_field_ben1;
}
void paramanage::set_b_field_ben1()
{
  b_field_ben1 = 10;
}
void paramanage::set_b_field_ben1(double value)
{
  b_field_ben1 = value;
}

/************************************************/
double paramanage::get_b_field_ben2()
{
  return b_field_ben2;
}
void paramanage::set_b_field_ben2()
{
  b_field_ben2 = 1;
}
void paramanage::set_b_field_ben2(double value)
{
  b_field_ben2 = value;
}

  // electric fields
/************************************************/
int paramanage::get_enable_epot()
{
  return enable_epot;
}
void paramanage::set_enable_epot()
{
  enable_epot = 0;
}
void paramanage::set_enable_epot(int value)
{
  enable_epot = value;
}

/************************************************/
double paramanage::get_el_mm_unit()
{
  return el_mm_unit;
}
void paramanage::set_el_mm_unit()
{
  el_mm_unit = 1;
}
void paramanage::set_el_mm_unit(double value)
{
  el_mm_unit = value;
}

/************************************************/
double paramanage::get_el_x_offset()
{
  return el_x_offset;
}
void paramanage::set_el_x_offset()
{
  el_x_offset = 0;
}
void paramanage::set_el_x_offset(double value)
{
  el_x_offset = value;
}

/************************************************/
double paramanage::get_el_y_offset()
{
  return el_y_offset;
}
void paramanage::set_el_y_offset()
{
  el_y_offset = 0;
}
void paramanage::set_el_y_offset(double value)
{
  el_y_offset = value;
}

/************************************************/
double paramanage::get_el_z_offset()
{
  return el_z_offset;
}
void paramanage::set_el_z_offset()
{
  el_z_offset = 0;
}
void paramanage::set_el_z_offset(double value)
{
  el_z_offset = value;
}

  // transmission module
/************************************************/
double paramanage::get_trans_b_pinch()
{
  return trans_b_pinch;
}
void paramanage::set_trans_b_pinch()
{
  trans_b_pinch = 10;
}
void paramanage::set_trans_b_pinch(double value)
{
  trans_b_pinch = value;
}

/************************************************/
double paramanage::get_trans_u_pinch()
{
  return trans_u_pinch;
}
void paramanage::set_trans_u_pinch()
{
  trans_u_pinch = 0;
}
void paramanage::set_trans_u_pinch(double value)
{
  trans_u_pinch = value;
}

/************************************************/
int paramanage::get_trans_steps()
{
  return trans_steps;
}
void paramanage::set_trans_steps()
{
  trans_steps = 10;
}
void paramanage::set_trans_steps(int value)
{
  trans_steps = value;
}

  // trapping module
/************************************************/
double paramanage::get_trap_start_x()
{
  return trap_start_x;
}
void paramanage::set_trap_start_x()
{
  trap_start_x = 0;
}
void paramanage::set_trap_start_x(double value)
{
  trap_start_x = value;
}

/************************************************/
double paramanage::get_trap_stop_x()
{
  return trap_stop_x;
}
void paramanage::set_trap_stop_x()
{
  trap_stop_x = -200;
}
void paramanage::set_trap_stop_x(double value)
{
  trap_stop_x = value;
}

/************************************************/
double paramanage::get_trap_step_x()
{
  return trap_step_x;
}
void paramanage::set_trap_step_x()
{
  trap_step_x = -5;
}
void paramanage::set_trap_step_x(double value)
{
  trap_step_x = value;
}

/************************************************/
double paramanage::get_trap_start_y()
{
  return trap_start_y;
}
void paramanage::set_trap_start_y()
{
  trap_start_y = 0;
}
void paramanage::set_trap_start_y(double value)
{
  trap_start_y = value;
}

/************************************************/
double paramanage::get_trap_stop_y()
{
  return trap_stop_y;
}
void paramanage::set_trap_stop_y()
{
  trap_stop_y = 45;
}
void paramanage::set_trap_stop_y(double value)
{
  trap_stop_y = value;
}

/************************************************/
double paramanage::get_trap_step_y()
{
  return trap_step_y;
}
void paramanage::set_trap_step_y()
{
  trap_step_y = 5;
}
void paramanage::set_trap_step_y(double value)
{
  trap_step_y = value;
}

/************************************************/
double paramanage::get_trap_energy_start()
{
  return trap_energy_start;
}
void paramanage::set_trap_energy_start()
{
  trap_energy_start = 32;
}
void paramanage::set_trap_energy_start(double value)
{
  trap_energy_start = value;
}

/************************************************/
double paramanage::get_trap_energy_end()
{
  return trap_energy_end;
}
void paramanage::set_trap_energy_end()
{
  trap_energy_end = 8;
}
void paramanage::set_trap_energy_end(double value)
{
  trap_energy_end = value;
}

/************************************************/
double paramanage::get_trap_theta_start()
{
  return trap_theta_start;
}
void paramanage::set_trap_theta_start()
{
  trap_theta_start = 80;
}
void paramanage::set_trap_theta_start(double value)
{
  trap_theta_start = value;
}

/************************************************/
double paramanage::get_trap_theta_step()
{
  return trap_theta_step;
}
void paramanage::set_trap_theta_step()
{
  trap_theta_step = 10;
}
void paramanage::set_trap_theta_step(double value)
{
  trap_theta_step = value;
}

/************************************************/
double paramanage::get_trap_phi_start()
{
  return trap_phi_start;
}
void paramanage::set_trap_phi_start()
{
  trap_phi_start = 0;
}
void paramanage::set_trap_phi_start(double value)
{
  trap_phi_start = value;
}

/************************************************/
int paramanage::get_trap_neg_y()
{
  return trap_neg_y;
}
void paramanage::set_trap_neg_y()
{
  trap_neg_y = 0;
}
void paramanage::set_trap_neg_y(int value)
{
  trap_neg_y = value;
}

/************************************************/
int paramanage::get_trap_y_plane()
{
  return trap_y_plane;
}
void paramanage::set_trap_y_plane()
{
  trap_y_plane = 1;
}
void paramanage::set_trap_y_plane(int value)
{
  trap_y_plane = value;
}

/************************************************/
int paramanage::get_trap_z_plane()
{
  return trap_z_plane;
}
void paramanage::set_trap_z_plane()
{
  trap_z_plane = 2;
}
void paramanage::set_trap_z_plane(int value)
{
  trap_z_plane = value;
}

/************************************************/
int paramanage::get_trap_max_mirrors()
{
  return trap_max_mirrors;
}
void paramanage::set_trap_max_mirrors()
{
  trap_max_mirrors = 250;  
}
void paramanage::set_trap_max_mirrors(int value)
{
  trap_max_mirrors = value;
}

/************************************************/
double paramanage::get_trap_mass()
{
  return trap_mass;
}
void paramanage::set_trap_mass()
{
  trap_mass = 1;
}
void paramanage::set_trap_mass(double value)
{
  trap_mass = value;
}

/************************************************/
double paramanage::get_trap_charge()
{
  return trap_charge;
}
void paramanage::set_trap_charge()
{
  trap_charge = 1;
}
void paramanage::set_trap_charge(double value)
{
  trap_charge = value;
}

/************************************************/
double paramanage::get_trap_max_step_length()
{
  return trap_max_step_length;
}
void paramanage::set_trap_max_step_length()
{
  trap_max_step_length = 0.01;
}
void paramanage::set_trap_max_step_length(double value)
{
  trap_max_step_length = value;
}

/************************************************/
double paramanage::get_trap_e_para_min()
{
  return trap_e_para_min;
}
void paramanage::set_trap_e_para_min()
{
  trap_e_para_min = 0;
}
void paramanage::set_trap_e_para_min(double value)
{
  trap_e_para_min = value;
}

/************************************************/
int paramanage::get_trap_calc_order()
{
  return trap_calc_order;
}
void paramanage::set_trap_calc_order()
{
  trap_calc_order = 0;
}
void paramanage::set_trap_calc_order(int value)
{
  trap_calc_order = value;
}

/************************************************/
double paramanage::get_trap_max_tof()
{
  return trap_max_tof;
}
void paramanage::set_trap_max_tof()
{
  trap_max_tof = 0.000003;
}
void paramanage::set_trap_max_tof(double value)
{
  trap_max_tof = value;
}

  // adi2fft settings
/************************************************/
int paramanage::get_rad_calc_mode()
{
  return rad_calc_mode;
}
void paramanage::set_rad_calc_mode()
{
  rad_calc_mode = 2;
}
void paramanage::set_rad_calc_mode(int value)
{
  rad_calc_mode = value;
}

/************************************************/
int paramanage::get_rad_shift()
{
  return rad_shift;
}
void paramanage::set_rad_shift()
{
  rad_shift = 1;
}
void paramanage::set_rad_shift(int value)
{
  rad_shift = value;
}

/************************************************/
double paramanage::get_antenna_temp()
{
  return antenna_temp;
}
void paramanage::set_antenna_temp()
{
  antenna_temp = 0;
}
void paramanage::set_antenna_temp(double value)
{
  antenna_temp = value;
}

/************************************************/
double paramanage::get_impedance()
{
  return impedance;
}
void paramanage::set_impedance()
{
  impedance = 1;
}
void paramanage::set_impedance(double value)
{
  impedance = value;
}

/************************************************/
int paramanage::get_rad_atten()
{
  return rad_atten;
}
void paramanage::set_rad_atten()
{
  rad_atten = 1;
}
void paramanage::set_rad_atten(int value)
{
  rad_atten = value;
}

/************************************************/
double paramanage::get_antenna_pos()
{
  return antenna_pos;
}
void paramanage::set_antenna_pos()
{
  antenna_pos = 500;
}
void paramanage::set_antenna_pos(double value)
{
  antenna_pos = value;
}

/************************************************/
int paramanage::get_fft_on()
{
  return fft_on;
}
void paramanage::set_fft_on()
{
  fft_on = 1;
}
void paramanage::set_fft_on(int value)
{
  fft_on = value;
}

/************************************************/
double paramanage::get_fft_resample_tstep()
{
  return fft_resample_tstep;
}
void paramanage::set_fft_resample_tstep()
{
  fft_resample_tstep = 0.000003;
}
void paramanage::set_fft_resample_tstep(double value)
{
  fft_resample_tstep = value;
}

/************************************************/
int paramanage::get_fft_max_npts()
{
  return fft_max_npts;
}
void paramanage::set_fft_max_npts()
{
  fft_max_npts = 20000000;
}
void paramanage::set_fft_max_npts(int value)
{
  fft_max_npts = value;
}

/************************************************/
int paramanage::get_leave_antenna()
{
  return leave_antenna;
}
void paramanage::set_leave_antenna()
{
  leave_antenna = 0;
}
void paramanage::set_leave_antenna(int value)
{
  leave_antenna = value;
}
