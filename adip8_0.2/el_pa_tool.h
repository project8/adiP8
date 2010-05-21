void set_dipole(double dipole);  // sets dipole magnifying value
double get_dipole(void);
double get_max_voltage();
int get_symmetry();
double get_epot3D(double *position_cm,int *splat);
double get_epot2D(double *position_cm,int *splat);

void alloc_electric_arrays();
void free_electric_arrays();
int read_epot(char* pa_in_name);

double epot3d(double *position, int *e_tag);
// this function gives 3D potentials out of cylindric values from efield module

void efield3d_old(double *position, double *evec, int *e_tag);
// this function calcs the electric field vector E by gradient of potential

void efield3d(double *pos_cm, double *evec, int *e_tag);
// new efield3d with built in interpolation and direct access to potential array
