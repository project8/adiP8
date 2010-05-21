int signum(double x);
// return sign of a double number

double interpol_2dim(double d_z, double d_r,
                     double z_left,double z_right, double r_below, double r_above,
                     double y_left_below,  double y_right_below, 
                     double y_right_above, double y_left_above);

double interpol_3dim(double *position,
		     int *pos_left, int *pos_right,
		     double *oktett);
