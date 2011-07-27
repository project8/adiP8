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
void paramanage::set_max_mirrors()
{
  max_mirrors = 5000;
}

void paramanage::set_max_tof_in_sec()
{
  max_tof_in_sec = 0.000003;
}

void paramanage::set_max_loops()
{
  max_loops = 1000000;
}

void paramanage::set_max_step_length()
{
  max_step_length = 0.001;
}

void paramanage::set_save_every()
{
  save_every = 1;
}

void paramanage::set_para_e_loss()
{
  para_e_loss = 1;
}

void paramanage::set_perp_e_loss()
{
  perp_e_loss = 1;
}

void paramanage::set_calc_order()
{
  calc_order = 2;
}

void paramanage::set_rel_start_angle()
{
  enable_rel_start_angle = 0;
}

void paramanage::set_puls_time()
{
  puls_time = 1;
}

void paramanage::set_min_shrink_factor()
{
  min_shrink_factor = 0.000000001;
}

void paramanage::set_interpol_splat()
{
  interpol_splat = 1;
}

void paramanage::set_residual_gas_pressure()
{
  residual_gas_pressure = -1;
}

void paramanage::set_e_min_cooling()
{
  e_min_cooling = 1;
}

void paramanage::set_dipole_value()
{
  dipole_value = 0;
}

void paramanage::set_e_para_min()
{
  e_para_min = 0.001;
}

void paramanage::set_spec_in()
{
  spec_in = -420;
}

void paramanage::set_spec_out()
{
  spec_out = 420;
}

void paramanage::set_max_radius()
{
  max_radius = 475;
}

  // magnetic potential array
void paramanage::set_mag_mm_unit()
{
  mag_mm_unit = 5;
}

void paramanage::set_mag_x_offset()
{
  mag_x_offset = 0;
}

void paramanage::set_mag_y_offset()
{
  mag_y_offset = 0;
}

void paramanage::set_mag_z_offset()
{
  mag_z_offset = 0;
}

void paramanage::set_n_pot_array()
{
  n_pot_array = 30125000;
}

  // parabolic magnetic field
void paramanage::set_b_field_ben1()
{
  b_field_ben1 = 10;
}

void paramanage::set_b_field_ben2()
{
  b_field_ben2 = 1;
}

  // electric fields
void paramanage::set_enable_epot()
{
  enable_epot = 0;
}

void paramanage::set_el_mm_unit()
{
  el_mm_unit = 1;
}

void paramanage::set_el_x_offset()
{
  el_x_offset = 0;
}

void paramanage::set_el_y_offset()
{
  el_y_offset = 0;
}

void paramanage::set_el_z_offset()
{
  el_z_offset = 0;
}

  // transmission module
void paramanage::set_trans_b_pinch()
{
  trans_b_pinch = 10;
}

void paramanage::set_trans_u_pinch()
{
  trans_u_pinch = 0;
}

void paramanage::set_trans_steps()
{
  trans_steps = 10;
}

  // trapping module
void paramanage::set_trap_start_x()
{
  trap_start_x = 0;
}

void paramanage::set_trap_stop_x()
{
  trap_stop_x = -200;
}
void paramanage::set_trap_step_x()
{
  trap_step_x = -5;
}

void paramanage::set_trap_start_y()
{
  trap_start_y = 0;
}

void paramanage::set_trap_stop_y()
{
  trap_stop_y = 45;
}

void paramanage::set_trap_step_y()
{
  trap_step_y = 5;
}

void paramanage::set_trap_energy_start()
{
  trap_energy_start = 32;
}

void paramanage::set_trap_energy_end()
{
  trap_energy_end = 8;
}

void paramanage::set_trap_theta_start()
{
  trap_theta_start = 80;
}

void paramanage::set_trap_theta_step()
{
  trap_theta_step = 10;
}

void paramanage::set_trap_phi_start()
{
  trap_phi_start = 0;
}

void paramanage::set_trap_neg_y()
{
  trap_neg_y = 0;
}

void paramanage::set_trap_y_plane()
{
  trap_y_plane = 1;
}

void paramanage::set_trap_z_plane()
{
  trap_z_plane = 2;
}

void paramanage::set_trap_max_mirrors()
{
  trap_max_mirrors = 250;  
}

void paramanage::set_trap_mass()
{
  trap_mass = 1;
}

void paramanage::set_trap_charge()
{
  trap_charge = 1;
}

void paramanage::set_trap_max_step_length()
{
  trap_max_step_length = 0.01;
}

void paramanage::set_trap_e_para_min()
{
  trap_e_para_min = 0;
}

void paramanage::set_trap_calc_order()
{
  trap_calc_order = 0;
}

void paramanage::set_trap_max_tof()
{
  trap_max_tof = 0.000003;
}

  // adi2fft settings
void paramanage::set_rad_calc_mode()
{
  rad_calc_mode = 2;
}

void paramanage::set_rad_shift()
{
  rad_shift = 1;
}

void paramanage::set_antenna_temp()
{
  antenna_temp = 0;
}

void paramanage::set_impedance()
{
  impedance = 1;
}

void paramanage::set_rad_atten()
{
  rad_atten = 1;
}

void paramanage::set_antenna_pos()
{
  antenna_pos = 500;
}

void paramanage::set_fft_on()
{
  fft_on = 1;
}

void paramanage::set_fft_resample_tstep()
{
  fft_resample_tstep = 0.000003;
}

void paramanage::set_fft_max_npts()
{
  fft_max_npts = 20000000;
}

void paramanage::set_leave_antenna()
{
  leave_antenna = 0;
}
