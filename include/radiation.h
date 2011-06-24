struct t_tl_data
{
  //optional geometry parameters
  double y1;//cm, y position of 1st wire or strip
  double y2;//cm, y position of 2nd wire or strip
  double z1;//cm, z length of sq wg
  double rI;//wire radius (50 um) in cm, 38 awg
  double rO;//outer coax and circ wg radius in cm
  double l;//Length of strip in cm in z-dir
  double n;//correction factor for finite strips
  //parameters necessary for power calculation
  double Zw;//wave impedance =Z0 for TEM, k0/k*Z0 for TE 
  double Zc;//characteristic impedance, only relevant for TEM modes
  double skinD;//skin depth, m, freq. and temp. dependent
  double C;//capacitance per unit length, only relevant for TEM modes
  double R;//resistance per unit length, only relevant for TEM modes
  double att;//attenuation coefficient
  double vg;//group velocity
  double vp;//phase velocity
};
extern struct t_tl_data tl_data;

void init_data();
double coeff_of_t(double *efield, double *vel, int dir);
void calculate_tl_parameters();

void init_tl_data(bool offset);
int get_tl_efield(double *p, double *efield);

void init_pp_data();
int get_pp_efield(double *pos, double *efield);

void init_coax_data();
int get_coax_efield(double *pos, double *efield);

void init_sq_wg_data(double k0);
void print_sq_wg_power(double k0);
int get_sq_wg_efield(double *pos, double *efield);

void init_circ_wg_data(double k0);
void print_circ_wg_power(double k0);
int get_circ_wg_efield(double *pos, double *efield);

int get_circ_cavity_efield(double phase, double *pos, double *efield);
