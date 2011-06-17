struct t_tl_data
{
  double y1;//cm, y position of 1st wire or strip
  double y2;//cm, y position of 2nd wire or strip
  double rI;//wire radius (50 um) in cm, 38 awg
  double rO;//outer coax and circ wg radius in cm
  double l;//Length of strip in cm in z-dir
  double n;//correction factor for finite strips
  double Zw;//wave impedance =z0 for TEM, k0/k*z0 for TE 
  double C;//capacitance per unit length, only relevant for TEM modes
  double R;//resistance per unit length, only relevant for TEM modes
  double Zc;//characteristic impedance, only relevant for TEM modes
  double att;//attenuatio coefficient
};
extern struct t_tl_data tl_data;
extern double c;

void init_data();
double coeff_of_t(double *efield, double *vel, int dir);

void init_tl_data(bool offset);
int get_tl_efield(double *p, double *efield);
void init_pp_data();
int get_pp_efield(double *pos, double *efield);

void init_coax_data();
int get_coax_efield(double *pos, double *efield);

void init_sq_wg_data(double k0);
void set_sq_wg_Zw(double k0);
int get_sq_wg_efield(double *pos, double *efield);

void init_circ_wg_data(double k0);
void set_circ_wg_Zw(double k0);
int get_circ_wg_efield(double phase, double *pos, double *efield);

int get_circ_cavity_efield(double phase, double *pos, double *efield);
