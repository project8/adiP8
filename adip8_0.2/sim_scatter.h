double get_wq(double argument);
     // function to calculate an approximated cross-section
     // for electron-H2 inelastic scattering.
     // input value is total kinetic energy in eV.

double get_std_rand_1(void);
    // function gives a random number between 0 and 1
    // using only the C command

double get_residual_gas_density(void);
     // function to calculate the number of residual gas
     // particles depending on temperature and pressure

void scatter_init_rel_cross_section(void);
     // init vaiables for scattering detection

void scatter_add_step(double dx, double Ekin);
     // add dx/path fraction of current step to total fraction

int scatter_check_event(void);
     // check if dx/path fraction is reached already 

void scatter_store_sigmas(double ekin);
     // remembers all sigma values in case of scatter event

void scatter_store_angle_eloss(double ekin);
     // remembers all angle and eloss values in case of scatter event

double scatter_get_eloss(void);
     // returns precalculated energy loss

double scatter_get_angle(void);
     // returns precalculated scattering angle

int scatter_get_process(void);
     // returns the type of scattering process
     // 0 = no ionization, 1 = ionization
     // -1 = not scattered

double scatter_get_ratio(void);
     // returns the ratio between current and destination value

void scatter_debug(void);


