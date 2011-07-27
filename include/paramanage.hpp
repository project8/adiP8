class paramanage
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
public:
  /*********************************************/
  paramanage();
  ~paramanage();
  void restore_defaults();
  
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
  void set_puls_time();
  /*********************************************/
  void set_min_shrink_factor();
  /*********************************************/
  void set_interpol_splat();
  /*********************************************/
  void set_residual_gas_pressure();
  /*********************************************/
  void set_e_min_cooling();
  /*********************************************/
  void set_dipole_value();
  /*********************************************/
  void set_e_para_min();
  /*********************************************/
  void set_spec_in();
  /*********************************************/
  void set_spec_out();
  /*********************************************/
  void set_max_radius();
  /*********************************************/
  // magnetic potential array
  /*********************************************/
  void set_mag_mm_unit();
  /*********************************************/
  void set_mag_x_offset();
  /*********************************************/
  void set_mag_y_offset();
  /*********************************************/
  void set_mag_z_offset();
  /*********************************************/
  void set_n_pot_array();
  /*********************************************/
  // parabolic magnetic field
  /*********************************************/
  void set_b_field_ben1();
  /*********************************************/
  void set_b_field_ben2();
  // electric fields
  void set_enable_epot();
  /*********************************************/
  void set_el_mm_unit();
  /*********************************************/
  void set_el_x_offset();
  /*********************************************/
  void set_el_y_offset();
  /*********************************************/
  void set_el_z_offset();
  /*********************************************/
  // transmission module
  /*********************************************/
  void set_trans_b_pinch();
  /*********************************************/
  void set_trans_u_pinch();
  /*********************************************/
  void set_trans_steps();
  /*********************************************/
  // trapping module
  /*********************************************/
  void set_trap_start_x();
  /*********************************************/
  void set_trap_stop_x();
  /*********************************************/
  void set_trap_step_x();
  /*********************************************/
  void set_trap_start_y();
  /*********************************************/
  void set_trap_stop_y();
  /*********************************************/
  void set_trap_step_y();
  /*********************************************/
  void set_trap_energy_start();
  /*********************************************/
  void set_trap_energy_end();
  /*********************************************/
  void set_trap_theta_start();
  /*********************************************/
  void set_trap_theta_step();
  /*********************************************/
  void set_trap_phi_start();
  /*********************************************/
  void set_trap_neg_y();
  /*********************************************/
  void set_trap_y_plane();
  /*********************************************/
  void set_trap_z_plane();
  /*********************************************/
  void set_trap_max_mirrors();
  /*********************************************/
  void set_trap_mass();
  /*********************************************/
  void set_trap_charge();
  /*********************************************/
  void set_trap_max_step_length();
  /*********************************************/
  void set_trap_e_para_min();
  /*********************************************/
  void set_trap_calc_order();
  /*********************************************/
  void set_trap_max_tof();
  /*********************************************/
  // adi2fft settings
  /*********************************************/
  void set_rad_calc_mode();
  /*********************************************/
  void set_rad_shift();
  /*********************************************/
  void set_antenna_temp();
  /*********************************************/
  void set_impedance();
  /*********************************************/
  /*********************************************/
  void set_rad_atten();
  void set_antenna_pos();
  /*********************************************/
  void set_fft_on();
  /*********************************************/
  void set_fft_resample_tstep();
  /*********************************************/
  void set_fft_max_npts();
  /*********************************************/
  void set_leave_antenna();
};
