#include "TObject.h"

class parameter : public TObject
{
private:
  // adipark settings
  int adip_run_mode;
  int use_mag_pa;
  // tracking module
  int max_mirrors;
  double max_tof_in_sec;
  int max_loops;
  double max_step_length;
  int save_every;
  int para_e_loss;
  int perp_e_loss;
  int calc_order;
  int enable_rel_start_angle;
  double puls_time;
  double min_shrink_factor;
  int interpol_splat;
  double residual_gas_pressure;
  double e_min_cooling;
  double dipole_value;
  double e_para_min;
  double spec_in;
  double spec_out;
  double max_radius;
  // magnetic potential array
  double mag_mm_unit;
  double mag_x_offset;
  double mag_y_offset;
  double mag_z_offset;
  int n_pot_array;
  // parabolic magnetic field
  double b_field_ben1;
  double b_field_ben2;
  // electric fields
  int enable_epot;
  double el_mm_unit;
  double el_x_offset;
  double el_y_offset;
  double el_z_offset;
  // transmission module
  double trans_b_pinch;
  double trans_u_pinch;
  int trans_steps;
  // trapping module
  double trap_start_x;
  double trap_stop_x;
  double trap_step_x;
  double trap_start_y;
  double trap_stop_y;
  double trap_step_y;
  double trap_energy_start;
  double trap_energy_end;
  double trap_theta_start;
  double trap_theta_step;
  double trap_phi_start;
  int trap_neg_y;
  int trap_y_plane;
  int trap_z_plane;
  int trap_max_mirrors;  
  double trap_mass;      
  double trap_charge;
  double trap_max_step_length;
  double trap_e_para_min;
  int trap_calc_order;  
  double trap_max_tof;
  // adi2fft settings
  int rad_calc_mode;
  int rad_shift;
  double antenna_temp;
  double impedance;
  int rad_atten;
  double antenna_pos;
  int fft_on;
  double fft_resample_tstep;
  int fft_max_npts;
  int leave_antenna;
  double filter_lo;
  double filter_sf;
public:
  /*********************************************/
  parameter();
  virtual ~parameter();

  /**********Other Use**************************/
  void restore_defaults();
  void parse_file(char inifile[255]);
  
  /***********Get & Set****************************/
  // adipark settings
 /*********************************************/
  int get_adip_run_mode();
  void set_adip_run_mode();
  void set_adip_run_mode(int);

 /*********************************************/
  int get_use_mag_pa();
  void set_use_mag_pa();
  void set_use_mag_pa(int);

  /*********************************************/
  // tracking module
  /*********************************************/
  int get_max_mirrors();
  void set_max_mirrors();
  void set_max_mirrors(int);
 
  /*********************************************/
  double get_max_tof_in_sec();
  void set_max_tof_in_sec();
  void set_max_tof_in_sec(double);

  /*********************************************/
  int get_max_loops();
  void set_max_loops();
  void set_max_loops(int);

  /*********************************************/
  double get_max_step_length();
  void set_max_step_length();
  void set_max_step_length(double);

  /*********************************************/
  int get_save_every();
  void set_save_every();
  void set_save_every(int);

  /*********************************************/
  int get_para_e_loss();
  void set_para_e_loss();
  void set_para_e_loss(int);

  /*********************************************/
  int get_perp_e_loss();
  void set_perp_e_loss();
  void set_perp_e_loss(int);

  /*********************************************/
  int get_calc_order();
  void set_calc_order();
  void set_calc_order(int);

  /*********************************************/
  int get_rel_start_angle();
  void set_rel_start_angle();
  void set_rel_start_angle(int);

  /*********************************************/
  double get_puls_time();
  void set_puls_time();
  void set_puls_time(double);
  
  /*********************************************/
  double get_min_shrink_factor();
  void set_min_shrink_factor();
  void set_min_shrink_factor(double);

  /*********************************************/
  int get_interpol_splat();
  void set_interpol_splat();
  void set_interpol_splat(int);

  /*********************************************/
  double get_residual_gas_pressure();
  void set_residual_gas_pressure();
  void set_residual_gas_pressure(double);

  /*********************************************/
  double get_e_min_cooling();
  void set_e_min_cooling();
  void set_e_min_cooling(double);

  /*********************************************/
  double get_dipole_value();
  void set_dipole_value();
  void set_dipole_value(double);

  /*********************************************/
  double get_e_para_min();
  void set_e_para_min();
  void set_e_para_min(double);

  /*********************************************/
  double get_spec_in();
  void set_spec_in();
  void set_spec_in(double);

  /*********************************************/
  double get_spec_out();
  void set_spec_out();
  void set_spec_out(double);

  /*********************************************/
  double get_max_radius();
  void set_max_radius();
  void set_max_radius(double);

  /*********************************************/
  // magnetic potential array
  /*********************************************/
  double get_mag_mm_unit();
  void set_mag_mm_unit();
  void set_mag_mm_unit(double);

  /*********************************************/
  double get_mag_x_offset();
  void set_mag_x_offset();
  void set_mag_x_offset(double);

  /*********************************************/
  double get_mag_y_offset();
  void set_mag_y_offset();
  void set_mag_y_offset(double);

  /*********************************************/
  double get_mag_z_offset();
  void set_mag_z_offset();
  void set_mag_z_offset(double);

  /*********************************************/
  int get_n_pot_array();
  void set_n_pot_array();
  void set_n_pot_array(int);

  /*********************************************/
  // parabolic magnetic field
  /*********************************************/
  double get_b_field_ben1();
  void set_b_field_ben1();
  void set_b_field_ben1(double);

  /*********************************************/
  double get_b_field_ben2();
  void set_b_field_ben2();
  void set_b_field_ben2(double);

  // electric fields
  int get_enable_epot();
  void set_enable_epot();
  void set_enable_epot(int);

  /*********************************************/
  double get_el_mm_unit();
  void set_el_mm_unit();
  void set_el_mm_unit(double);

  /*********************************************/
  double get_el_x_offset();
  void set_el_x_offset();
  void set_el_x_offset(double);

  /*********************************************/
  double get_el_y_offset();
  void set_el_y_offset();
  void set_el_y_offset(double);

  /*********************************************/
  double get_el_z_offset();
  void set_el_z_offset();
  void set_el_z_offset(double);

  /*********************************************/
  // transmission module
  /*********************************************/
  double get_trans_b_pinch();
  void set_trans_b_pinch();
  void set_trans_b_pinch(double);

  /*********************************************/
  double get_trans_u_pinch();
  void set_trans_u_pinch();
  void set_trans_u_pinch(double);

  /*********************************************/
  int get_trans_steps();
  void set_trans_steps();
  void set_trans_steps(int);

  /*********************************************/
  // trapping module
  /*********************************************/
  double get_trap_start_x();
  void set_trap_start_x();
  void set_trap_start_x(double);

  /*********************************************/
  double get_trap_stop_x();
  void set_trap_stop_x();
  void set_trap_stop_x(double);

  /*********************************************/
  double get_trap_step_x();
  void set_trap_step_x();
  void set_trap_step_x(double);

  /*********************************************/
  double get_trap_start_y();
  void set_trap_start_y();
  void set_trap_start_y(double);
  
  /*********************************************/
  double get_trap_stop_y();
  void set_trap_stop_y();
  void set_trap_stop_y(double);

  /*********************************************/
  double get_trap_step_y();
  void set_trap_step_y();
  void set_trap_step_y(double);

  /*********************************************/
  double get_trap_energy_start();
  void set_trap_energy_start();
  void set_trap_energy_start(double);

  /*********************************************/
  double get_trap_energy_end();
  void set_trap_energy_end();
  void set_trap_energy_end(double);
  
  /*********************************************/
  double get_trap_theta_start();
  void set_trap_theta_start();
  void set_trap_theta_start(double);

  /*********************************************/
  double get_trap_theta_step();
  void set_trap_theta_step();
  void set_trap_theta_step(double);

  /*********************************************/
  double get_trap_phi_start();
  void set_trap_phi_start();
  void set_trap_phi_start(double);

  /*********************************************/
  int get_trap_neg_y();
  void set_trap_neg_y();
  void set_trap_neg_y(int);

  /*********************************************/
  int get_trap_y_plane();
  void set_trap_y_plane();
  void set_trap_y_plane(int);

  /*********************************************/
  int get_trap_z_plane();
  void set_trap_z_plane();
  void set_trap_z_plane(int);
  
  /*********************************************/
  int get_trap_max_mirrors();
  void set_trap_max_mirrors();
  void set_trap_max_mirrors(int);

  /*********************************************/
  double get_trap_mass();
  void set_trap_mass();
  void set_trap_mass(double);

  /*********************************************/
  double get_trap_charge();
  void set_trap_charge();
  void set_trap_charge(double);

  /*********************************************/
  double get_trap_max_step_length();
  void set_trap_max_step_length();
  void set_trap_max_step_length(double);

  /*********************************************/
  double get_trap_e_para_min();
  void set_trap_e_para_min();
  void set_trap_e_para_min(double);

  /*********************************************/
  int get_trap_calc_order();
  void set_trap_calc_order();
  void set_trap_calc_order(int);

  /*********************************************/
  double get_trap_max_tof();
  void set_trap_max_tof();
  void set_trap_max_tof(double);

  /*********************************************/
  // adi2fft settings
  /*********************************************/
  int get_rad_calc_mode();
  void set_rad_calc_mode();
  void set_rad_calc_mode(int);

  /*********************************************/
  int get_rad_shift();
  void set_rad_shift();
  void set_rad_shift(int);

  /*********************************************/
  double get_antenna_temp();
  void set_antenna_temp();
  void set_antenna_temp(double);

  /*********************************************/
  double get_impedance();
  void set_impedance();
  void set_impedance(double);

  /*********************************************/
  int get_rad_atten();
  void set_rad_atten();
  void set_rad_atten(int);

  /*********************************************/
  double get_antenna_pos();
  void set_antenna_pos();
  void set_antenna_pos(double);

  /*********************************************/
  int get_fft_on();
  void set_fft_on();
  void set_fft_on(int);
  
  /*********************************************/
  double get_fft_resample_tstep();
  void set_fft_resample_tstep();
  void set_fft_resample_tstep(double);

  /*********************************************/
  int get_fft_max_npts();
  void set_fft_max_npts();
  void set_fft_max_npts(int);

  /*********************************************/
  int get_leave_antenna();
  void set_leave_antenna();
  void set_leave_antenna(int);

  /*********************************************/
  double get_filter_lo();
  void set_filter_lo();
  void set_filter_lo(double);

  /*********************************************/
  double get_filter_sf();
  void set_filter_sf();
  void set_filter_sf(double);

  /*********************************************/
  /*********************************************/
  ClassDef(parameter, 1);
};
