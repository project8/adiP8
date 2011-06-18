#define MM_PER_UNIT parameter.el_mm_unit
#define X_OFFSET_IN_CM parameter.el_x_offset
#define Y_OFFSET_IN_CM parameter.el_y_offset
#define Z_OFFSET_IN_CM parameter.el_z_offset
#define MAG_MM_PER_UNIT parameter.mag_mm_unit
#define MAG_X_OFFSET_IN_CM parameter.mag_x_offset
#define MAG_Y_OFFSET_IN_CM parameter.mag_y_offset
#define MAG_Z_OFFSET_IN_CM parameter.mag_z_offset
#define SPEC_IN parameter.spec_in
#define SPEC_OUT parameter.spec_out
#define N_ARRAY parameter.n_pot_array
#define MAX_RADIUS parameter.max_radius
#define MAX_LOOPS  parameter.max_loops
#define USE_MAG_PA parameter.use_mag_pa
#define ENABLE_EPOT parameter.enable_epot
#define SAVE_EVERY parameter.save_every
#define PULS_TIME  parameter.puls_time
#define MIN_SHRINK_FACTOR parameter.min_shrink_factor
#define INTERPOL_SPLAT parameter.interpol_splat
#define B_PINCH parameter.trans_b_pinch
#define U_PINCH parameter.trans_u_pinch
#define N_TRANSMISSION parameter.trans_steps
#define PARA_ENERGY_LOSS parameter.para_e_loss
#define PERP_ENERGY_LOSS parameter.perp_e_loss

#define Clight 299792458.0   // PDG 07 m/s
#define M0Clight2 510998.91  // PDG 07
#define MM2CM 0.1            // scale mm value to cm value
#define M2CM 100.            // scale m value to cm value
#define CM2MM 10.            // scale cm value to mm value
#define US2S 1.e-6            // scale us value to s value
#define S2US 1.e+6            // scale s value to us value
#define TESLA2GAUSS 10000.   // scale tesla to gauss
#define Echarge 1.602176487e-19    // charge of single electron in C, pdg 07
#define erg2eV  6.24150974e+11    // converts CGS unit "erg" to "eV"
#define SI2esE  3e9          // converts SI units to CGS unit "esE"
#define OMEGA0 1.758820150e+11 // nrel electron cyclotron frequency (rad/s) at B=1t,pdg 07
#define EPS0 8.854187817e-12  //permittivity of free space pdg 10, in F/m
#define Z0 376.7           //Impedance of free space in Ohms
#define CU_R 2.15e-9           //Resititivy of Cu @ 80K in Ohm-m 
#define RM 0.015             //Resititivy of Cu / skin depth @ OMEGA0, Ohms
#define K_BOL 1.3806488e-23     //Boltzam constant, J/K

#define display 0           // output switch mode: 0=nothing, 1=some, 2=all 

struct particle_data
{
  double start_pos[3];              // origin of particle
  double position[3];               // actual position
  double e_start;                   // kinetic energy at start_pos
  double b_start[3];                // magnetic vector bat start_pos
  double b_start_value;             // value of b_start
  double u_start;                   // el. potential at start_pos
  double sin2_alpha_start;          // sin² of starting alpha
  double starting_theta;            // 0 to 89 degree
  double starting_phi;              // 0 to 359 degree
  double phase;                     // radians
  double doppler_phase;              // radians, includes doppler shift
  double omega;                     // in degree per us
  double gamma_start;               // rel. gamma factor at start
  double b_vec[3];                  // actual magnetic vector
  double b_vec_old[3];              // mag. vector or prev. step used for curv. calc.
  double b_value;                   // value of b_vec
  double b_value_old;               // value of b_vec_old
  double v_para[3];                 // actual parallel velocity vector in m/s
  double v_para_value;              // value of it
  double v_signum;                  // director relative to magnetic field line
  double v_perp_start[3];           // perpendicular velocity vector at start_pos
  double v_perp_start_value;        // and its value
  double v_perp_value;              // actual perp velocity value
  double v_drift[3];                // velocity vector for all drift parts
  double exb_vel[3];                // vector for exb drift velocity
  double dbxb_vel[3];               // vector for gradient drift velocity
  double rxb_vel[3];                // vector for curvature drift velocity
  double curv_rad;                  // value of curvature radius
  double e_para;                    // parallel energy
  double e_para_min;                // lower energy limit
  double e_perp;                    // perpendicular energy only of magnetic motion, without ExB
  double e_kin;                     // kinetic energy
  double e_pot;                     // potential energy (el. field)
  double e_perp_all;                // kinetic energy in all perpendicular motions, with ExB
  double e_cycl;                    // kinetic energy only of cyclotron and gradientdrift motion
  double e_ExB;                     // kinetic energy in ExB motion
  double e_curv;                    // kinetic energy in curvation motion
  double e_grad;                    // kinetic energy in gradient drift motion
  double e_syncro;                  // perp energy loss due to syncrotron radiation for single step
  double e_syncro_sum;              // sum of all perp energy losses due to syncrotron radiation
  double delta_e_para;              // change of e_para, needed for e_accel
  double e_accel;                   // para energy loss due to acceleration
  double e_accel_sum;               // sum of all para energy losses
  double e_scatter_para;            // energy loss para part due to scattering
  double e_scatter_perp;            // and the same for perp part
  double cyclrad;                   // cyclotron radius in cm
  double mass;                      // mass in eV
  double charge;                    // charge in e
  double time_of_flight;            // time of flight
  double delta_tof;                 // time of flight for single step
  double max_tof;                   // upper limit of flight time
  double step_final[3];             // final step to go
  double max_step_length;           // upper limit for step length
  double shrinkfactor;              // factor to shrink step length under max_step_length
  int    mirrors;                   // counting number of mirror
  int    max_mirrors;               // upper limit for mirrors
  int    calc_order;                // tag to switch on and off drift calculation
                                    // 0=alldrift 1=exb only 2=RxB and gradBxB 3=nothing

};

#define max_parameter 65
 
struct t_parameter{
  int inited[max_parameter+1];
  char name[max_parameter+1][32];
  double el_mm_unit;
  double el_x_offset;
  double el_y_offset;
  double el_z_offset;
  double mag_mm_unit;
  double mag_x_offset;
  double mag_y_offset;
  double mag_z_offset;
  double spec_in;
  double spec_out;
  double max_radius;
  int save_every;
  double puls_time;
  double min_shrink_factor;
  int interpol_splat;
  int max_loops;
  int n_pot_array;
  int use_mag_pa;
  int enable_epot;
  int enable_rel_start_angle;
  double trap_theta_start;
  double trap_theta_step;
  double trap_energy_start;
  double trap_energy_end;
  double trap_phi_start;
  int trap_neg_y;
  int trap_y_plane;
  int trap_z_plane;
  int trap_max_mirrors;  
  double trap_mass;      
  double trap_charge;
  double trap_max_step_length;
  double trap_e_para_min;
  double trap_max_tof;
  int trap_calc_order;  
  double trans_b_pinch;
  double trans_u_pinch;
  int trans_steps;
  double trap_start_x;
  double trap_stop_x;
  double trap_step_x;
  double trap_start_y;
  double trap_stop_y;
  double trap_step_y;
  int para_e_loss;
  int perp_e_loss;
  int calc_order;
  double dipole_value;
  int max_mirrors;
  double max_tof_in_sec;
  double e_para_min;
  double max_step_length;
  double e_min_cooling;
  double residual_gas_pressure;

  double b_field_ben1;
  double b_field_ben2;

  int rad_calc_mode;
  int rad_shift;
  double antenna_temp;
  double antenna_pos;
  double impedance;
  int rad_atten;

  int fft_on;
  double fft_resample_tstep;
  int fft_max_npts;
  char filename[255];
};

extern struct t_parameter parameter;

void init_parameters();

int load_init_data(char* filename);

int set_parameter(const char* identifier, double value);
