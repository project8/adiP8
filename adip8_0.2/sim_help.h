void sim_help_show_run_file(char* filename);

void sim_help_save_particle_data(FILE* f_fly,struct particle_data* particle);

void sim_help_force_save_data(); 

void sim_help_reset_save_data();

void sim_help_init_scatter_file(char* filename);

void sim_help_store_scatter_data(struct particle_data* particle,
				 double e_loss, double angle, int process);

void do_debug(char* nr, struct particle_data* pdata);
