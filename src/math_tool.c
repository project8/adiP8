using namespace std;
#include <cmath>
#include "math_tool.h"

/* ************************** signum ********************** */
int signum(double x)
{
  if (x >= 0) {
    return +1;                  // returns +1 if input >= zero
  } else {
    return -1;                  // returns -1 if input < zero
  }
}


/* ************************ interpol_2dim ********************** */
double interpol_2dim(double d_z, double d_r, double z_left, double z_right, double r_below, double r_above, double y_left_below, double y_right_below, double y_right_above, double y_left_above)
{
  double t, u;                  //see NUMERICAL RECIPES in C, chapter 3.6, page 123
  t = (d_z - z_left) / (z_right - z_left);	//relative distance in z
  u = (d_r - r_below) / (r_above - r_below);
  return (1 - t) * (1 - u) * y_left_below + t * (1 - u) * y_right_below + t * u * y_right_above + (1 - t) * u * y_left_above;
}


/* ************************ interpol_3dim ********************** */
double interpol_3dim(double *position, int *pos_left, int *pos_right, double *oktett)
{
  int comp;
  double rel_pos[3];            // relative position
  // see Numerical Recipes in C, chapter 3.6, page 123
  // and use 3 dimentions
  for (comp = 0; comp <= 2; comp++) {
    rel_pos[comp] = (position[comp] - (double) pos_left[comp]) / ((double) pos_right[comp] - (double) pos_left[comp]);
  }
  return
      (1. - rel_pos[0]) * (1. - rel_pos[1]) * (1. - rel_pos[2]) * oktett[0] +
      (rel_pos[0]) * (1. - rel_pos[1]) * (1. - rel_pos[2]) * oktett[1] +
      (rel_pos[0]) * (rel_pos[1]) * (1. - rel_pos[2]) * oktett[2] +
      (1. - rel_pos[0]) * (rel_pos[1]) * (1. - rel_pos[2]) * oktett[3] +
      (1. - rel_pos[0]) * (rel_pos[1]) * (rel_pos[2]) * oktett[4] +
      (rel_pos[0]) * (rel_pos[1]) * (rel_pos[2]) * oktett[5] + (rel_pos[0]) * (1. - rel_pos[1]) * (rel_pos[2]) * oktett[6] + (1. - rel_pos[0]) * (1. - rel_pos[1]) * (rel_pos[2]) * oktett[7];
}
