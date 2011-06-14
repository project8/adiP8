using namespace std;
#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <string>
#include <cstring>
#include "paramanage.h"
#include "matrix_tool.h"
#include "vector_tool.h"
#include "math_tool.h"
#include "el_pa_tool.h"

static int nx, ny, nz, sym;
static int max_nx = 20000;      // xline length to read entire lines in pa's
static void *electric_potential_array;
static double maxvoltage;
double dipol = 0.;

void set_dipole(double dipole)  // magnifying value for efield module
{
  dipol = dipole;
}

double get_dipole(void)
{
  return dipol;
}

double get_max_voltage()
{
  return maxvoltage;
}

int get_symmetry()
{
  return sym;
}

void alloc_electric_arrays()
{
  printf("\nAllocating memory for electric potential array.\n");
  electric_potential_array = malloc(N_ARRAY * sizeof(double));
  if (electric_potential_array == 0) {
    printf("\nERROR: not enough memory for N_POT_ARRAY of %d!\n\n", N_ARRAY);
  } else {
    printf("                                      DONE.\n");
  }
}

void free_electric_arrays()
{
  free(electric_potential_array);
  return;
}


int read_epot(char *el_filename)
{

  FILE *pa_in;                  // input file
  int x, y, z;
  int error = 0;

  double pot_xline[max_nx];
  char pa_in_name[255];

  struct header_3d {
    int mode;
    int symmetry;
    double max_voltage;
    int nx;
    int ny;
    int nz;
    int mirror;
  };

  struct header_3d h;

  strcpy(pa_in_name, el_filename);
  strcat(pa_in_name, ".pa0");   // data input from file with extension .pa0
  pa_in = fopen(pa_in_name, "r");  //* open electric PA-File in read mode */

  if (pa_in == (FILE *) 0) {     // disk error handling
    fprintf(stderr, "\nERROR: Can't read file %s\n\n", pa_in_name);
    error = 1;
  } else {
    printf("\nReading electric potential array: %s", pa_in_name);
    fread(&h, sizeof(h), 1, pa_in);  /* read the header out of the file */

    maxvoltage = h.max_voltage;
    sym = h.symmetry;
    printf("\nHeader informations:\n");
    printf("  mode %d\n", h.mode);
    printf("  symmetrie %d\n", sym);
    printf("  max_voltage %f\n", maxvoltage);
    printf("  Potentials at electrodes = potential + %f\n", 2 * h.max_voltage);
    printf("  mirror %d\n", h.mirror);
    printf("  nx, ny, nz : %d,%d,%d\n", h.nx, h.ny, h.nz);

    if (h.nx * h.ny * h.nz > N_ARRAY) {
      printf("\nERROR: not enough memory allocated! N_POT_ARRAY of %d needed!\n\n", h.nx * h.ny * h.nz);
      error = 1;
    } else {
      nx = h.nx;
      ny = h.ny;
      nz = h.nz;
      if (nx > max_nx) {
        printf("ERROR: read_epot: max_nx lower than nx!\n\n");
      } else {
        for (z = 0; z < h.nz; z++) {
          for (y = 0; y < h.ny; y++) {

            cout << "  reading: " << (int) (100. * (float) (y + 1) * (float) (z + 1) / (float) (ny) / (float) (nz)) << " %\r" << flush;

            fread(pot_xline, sizeof(double), nx, pa_in);

            for (x = 0; x < h.nx; x++) {
              //          read(pa_in, &pot, sizeof(double));          /* read potential of an (x/y/z)-Koord. */

              if (pot_xline[x] > h.max_voltage) {
                pot_xline[x] = pot_xline[x] - h.max_voltage;
              }
              write_matrix_value(electric_potential_array, h.nx, h.ny, h.nz, 0, x, y, z, pot_xline[x]);
            }
          }
        }
        printf("                                      DONE.\n");
      }
    }
    fclose(pa_in);
  }
  return error;
}


double get_epot3D(double *position_cm, int *splat)
{
  double oktett[8];             // potential values at all 8 corners
  double pos_backup[3];         // real position without any mirroring
  double position[3];           // position used for epot calculation
  double epot_interpol;         // interpolated potential value
  int pos_left[3];              // position of lower left corner in front - position of oktett[0]
  int pos_right[3];             // position of upper right corner in back - position of oktett[5]
  int electrode_count;          // counting points located in electrodes
  int comp;

  *splat = 0;
  vector_times_scalar(position_cm, 1., position);  // copy position vector
  position[0] = (position[0] - X_OFFSET_IN_CM) * CM2MM / MM_PER_UNIT;  // transform to simion grid
  position[1] = (position[1] - Y_OFFSET_IN_CM) * CM2MM / MM_PER_UNIT;
  position[2] = (position[2] - Z_OFFSET_IN_CM) * CM2MM / MM_PER_UNIT;
  vector_times_scalar(position, 1., pos_backup);  // backup position vector

  for (comp = 0; comp <= 2; comp++) {
    pos_left[comp] = (int) (position[comp] + 5000.) - 5000;  // corner 0
    // 5000 shifts all to positive than truncate and shift back
    // this is only needed when negative position
    if (pos_left[comp] < 0.) {
      pos_left[comp] = -pos_left[comp];  // use positive value
      pos_right[comp] = pos_left[comp] - 1;  // corner 5
    } else {
      pos_right[comp] = pos_left[comp] + 1;  // corner 5
    }
    if (position[comp] < 0.) {
      position[comp] = -position[comp];  // use positive values
    }
  }
  if ((position[0] > nx) || (position[1] > ny) || (position[2] > nz)) {
    return maxvoltage;          // watch for electrodes
  } else {                      // read potential in all 8 corners of cubus
    oktett[0] = read_matrix_value(electric_potential_array, nx, ny, nz, 0, pos_left[0], pos_left[1], pos_left[2]);
    oktett[1] = read_matrix_value(electric_potential_array, nx, ny, nz, 0, pos_right[0], pos_left[1], pos_left[2]);
    oktett[2] = read_matrix_value(electric_potential_array, nx, ny, nz, 0, pos_right[0], pos_right[1], pos_left[2]);
    oktett[3] = read_matrix_value(electric_potential_array, nx, ny, nz, 0, pos_left[0], pos_right[1], pos_left[2]);
    oktett[4] = read_matrix_value(electric_potential_array, nx, ny, nz, 0, pos_left[0], pos_right[1], pos_right[2]);
    oktett[5] = read_matrix_value(electric_potential_array, nx, ny, nz, 0, pos_right[0], pos_right[1], pos_right[2]);
    oktett[6] = read_matrix_value(electric_potential_array, nx, ny, nz, 0, pos_right[0], pos_left[1], pos_right[2]);
    oktett[7] = read_matrix_value(electric_potential_array, nx, ny, nz, 0, pos_left[0], pos_left[1], pos_right[2]);
    electrode_count = 0;
    for (comp = 0; comp <= 7; comp++) {  // patch epot when in electrode
      if (oktett[comp] >= maxvoltage / 2.) {
        oktett[comp] = oktett[comp] - maxvoltage;
        electrode_count++;
      }
    }
    for (comp = 0; comp <= 2; comp++) {  // redo corner 0 and 5 for interpolation
      pos_left[comp] = (int) (pos_backup[comp] + 5000.) - 5000;
      // now positive and negative values needed
      pos_right[comp] = pos_left[comp] + 1;
    }
    epot_interpol = interpol_3dim(pos_backup, pos_left, pos_right, oktett);
    // do interpolation
    //      if (electrode_count > 0) printf("electrodes: %d\n",electrode_count);
    if (electrode_count >= INTERPOL_SPLAT) {
      *splat = 1;
    }
    // if too many corners in electrode than interpolated point in electrode
    return epot_interpol;       // return epot value
  }
}

double get_epot2D(double *position_cm, int *splat)
{
  double oktett[8];             // potential values at all 8 corners
  double pos_backup[3];         // real position without any mirroring
  double position[3];           // position used for epot calculation
  double epot_interpol;         // interpolated potential value
  int pos_left[3];              // position of lower left corner in front - position of oktett[0]
  int pos_right[3];             // position of upper right corner in back - position of oktett[5]
  int electrode_count;          // counting points located in electrodes
  int comp;

  *splat = 0;
  position[0] = (position_cm[0] - X_OFFSET_IN_CM) * CM2MM / MM_PER_UNIT;  // transform to simion grid
  position[1] = (position_cm[1] - Y_OFFSET_IN_CM) * CM2MM / MM_PER_UNIT;
  position[2] = 0.;
  vector_times_scalar(position, 1., pos_backup);  // backup position vector

  for (comp = 0; comp <= 2; comp++) {
    pos_left[comp] = (int) (position[comp] + 5000.) - 5000;  // corner 0
    // 5000 shifts all to positive than truncate and shift back
    // this is only needed when negative position
    if (pos_left[comp] < 0.) {
      pos_left[comp] = -pos_left[comp];  // use positive value
      pos_right[comp] = pos_left[comp] - 1;  // corner 5
    } else {
      pos_right[comp] = pos_left[comp] + 1;  // corner 5
    }
    if (position[comp] < 0.) {
      position[comp] = -position[comp];  // use positive values
    }
  }
  if ((position[0] > nx) || (position[1] > ny) || (position[2] > nz) || (pos_left[0] > nx) || (pos_left[1] > ny) || (pos_left[2] > nz) || (pos_right[0] > nx) || (pos_right[1] > ny) || (pos_right[2] > nz)) {
    return maxvoltage;          // watch for electrodes
  } else {                      // read potential in all 8 corners of cubus
    oktett[0] = read_matrix_value(electric_potential_array, nx, ny, nz, 0, pos_left[0], pos_left[1], pos_left[2]);
    oktett[1] = read_matrix_value(electric_potential_array, nx, ny, nz, 0, pos_right[0], pos_left[1], pos_left[2]);
    oktett[2] = read_matrix_value(electric_potential_array, nx, ny, nz, 0, pos_right[0], pos_right[1], pos_left[2]);
    oktett[3] = read_matrix_value(electric_potential_array, nx, ny, nz, 0, pos_left[0], pos_right[1], pos_left[2]);

    oktett[4] = 0.;
    oktett[5] = 0.;
    oktett[6] = 0.;
    oktett[7] = 0.;

    electrode_count = 0;
    for (comp = 0; comp <= 7; comp++) { // patch epot when in electrode
      //          if (oktett[comp] >= maxvoltage)
      if (oktett[comp] > maxvoltage / 2.) {
        oktett[comp] = oktett[comp] - maxvoltage;
        electrode_count++;
      }
    }
    for (comp = 0; comp <= 2; comp++) {  // redo corner 0 and 5 for interpolation
      pos_left[comp] = (int) (pos_backup[comp] + 5000.) - 5000;
      // now positive and negative values needed
      pos_right[comp] = pos_left[comp] + 1;
    }
    epot_interpol = interpol_3dim(pos_backup, pos_left, pos_right, oktett);
    // do interpolation
    //      if (electrode_count > 0) printf("electrodes: %d\n",electrode_count);
    if (electrode_count >= INTERPOL_SPLAT) {
      *splat = 1;
    }
    // if too many corners in electrode than interpolated point in electrode
    return epot_interpol;       // return epot value
  }
}


/******************** epot3d ***********************/
double epot3d(double *position, int *e_tag)
{                               // this function gives 3 dim potentials out of cylindric values from efield module
  double cyl_pos[3];            // 2D coordinates
  double e_pot = 0.;            // electric potential
  int splat;

  if (ENABLE_EPOT == 1)         // switch on and off retarding potential
  {
    if (get_symmetry() == 0) {
      cyl_pos[0] = position[0]; // 2D z coordinate equal to 3D x coordinate
      cyl_pos[1] = sqrt(position[1] * position[1] + position[2] * position[2]);
      cyl_pos[2] = 0.;
      // read 2D flux radius by using distance between
      // x axis and position
      e_pot = get_epot2D(cyl_pos, &splat);
    } else {
      e_pot = get_epot3D(position, &splat);  // use 3D potential array
    }
    if (splat == 1) {            // tag for electrode detection
      *e_tag = *e_tag + 1;      // use e_tag to inform parent routine
      //          printf("WARNING: Potential was read in Electrode area! Possible Impact!\n");
    }
  }
  if (dipol != 0.) {
    e_pot = e_pot - dipol * (position[2] + 50.) / 100.;
  }
  // add a linear potential to force drift motion
  return e_pot;                 // return the value of epot
}


/******************** efield3d ***********************/
void efield3d_old(double *position, double *evec, int *e_tag)
{                               // this function calculates the electric field vector E by gradient of potential
  double pos1[3], pos3[3];
  double delta = MM_PER_UNIT;
  int komp, splat_tag;

  delta = delta * MM2CM;        // delta in mm must be used in cm
  if (delta == 0.) {
    printf("Error: MM_PER_UNIT constant not set!\n");
  }
  for (komp = 0; komp <= 2; komp++) {
    vector_times_scalar(position, 1., pos1);
    vector_times_scalar(position, 1., pos3);
    pos1[komp] = pos1[komp] - delta;
    pos3[komp] = pos3[komp] + delta;
    evec[komp] = -(epot3d(pos3, &splat_tag) - epot3d(pos1, &splat_tag))
        * M2CM / (2. * delta);
  }
  if (splat_tag > 0) {
    *e_tag = *e_tag + 1;
  }
}


void get_4corner_potentials(double *pos_left, double *pos_right, double *corner)
{
  corner[0] = read_matrix_value(electric_potential_array, nx, ny, nz, 0, (int) pos_left[0], (int) pos_left[1], (int) pos_left[2]);

  corner[1] = read_matrix_value(electric_potential_array, nx, ny, nz, 0, (int) pos_right[0], (int) pos_left[1], (int) pos_left[2]);

  corner[2] = read_matrix_value(electric_potential_array, nx, ny, nz, 0, (int) pos_right[0], (int) pos_right[1], (int) pos_left[2]);

  corner[3] = read_matrix_value(electric_potential_array, nx, ny, nz, 0, (int) pos_left[0], (int) pos_right[1], (int) pos_left[2]);
}


/******************** efield3d_new ***********************/
void efield3d(double *pos_cm, double *evec, int *e_tag)
{                               // this function calculates the electric field vector E by gradient of potential
  double azimuth;               // angle between 1 and 2 component
  double keep_pos[3];           // real position without any mirroring
  double grid_pos[3];           // position on grid used for epot calculation
  double gpos_ll[3];            // grid point upper left 
  double gpos_ur[3];            // grid point lower right 
  double corner_pot[8];         // potential at corners
  double delta = MM_PER_UNIT;   // 
  double upper_pot_derivative;
  double lower_pot_derivative;
  double left_pot_derivative;
  double right_pot_derivative;
  double deltax, deltay;        // for distance of pos to grid

  //  int electrode_count; // counting points located in electrodes
  int comp;                     // vector component counter
  int shift = 10000;            // position shift to get the right corner points

  if (ENABLE_EPOT == 1) {        // switch on and off retarding potential
    if (get_symmetry() == 0) {   // means cylindrical 2D symmetry
      azimuth = atan2(pos_cm[1], pos_cm[2]) + 3 * M_PI / 2.;
      // important to use atan2 because of different quadrants
      // see "man atan2" for more information
      // 3*M_PI/2 needed for transformation to rotation coordinates
      // -> only positive rotation possible
      // -> need offset of 90 degrees

      grid_pos[0] = (pos_cm[0] - X_OFFSET_IN_CM) * CM2MM / MM_PER_UNIT;
      grid_pos[1] = (sqrt(pos_cm[1] * pos_cm[1] + pos_cm[2] * pos_cm[2])
                     - Y_OFFSET_IN_CM) * CM2MM / MM_PER_UNIT;
      grid_pos[2] = 0.;         // transform position to 2D and to simion grid

      vector_times_scalar(grid_pos, 1., keep_pos);
      // backup grid position vector

      for (comp = 0; comp <= 2; comp++) {
        gpos_ll[comp] = (int) (grid_pos[comp] + (double) shift) - shift;
        // upper left corner
        // shift=10000 shifts all to positive than truncate and shift back
        // this is only needed when negative position accure to get the upper left corner
        if (gpos_ll[comp] < 0.) {
          gpos_ll[comp] = -gpos_ll[comp];  // use positive value
          gpos_ur[comp] = gpos_ll[comp] - 1;  // lower right corner for - gpos
        } else {
          gpos_ur[comp] = gpos_ll[comp] + 1;  // lower right corner for + gpos
        }
        if (grid_pos[comp] < 0.) {
          grid_pos[comp] = -grid_pos[comp];  // use positive values
        }
      }

      if ((grid_pos[0] > nx) || (grid_pos[1] > ny) || (grid_pos[2] > nz) || (gpos_ll[0] > nx) || (gpos_ll[1] > ny) || (gpos_ll[2] > nz) || (gpos_ur[0] > nx) || (gpos_ur[1] > ny) || (gpos_ur[2] > nz)) {
        e_tag++;                // splat on electrode
        vector_times_scalar(grid_pos, 0., evec);  // evec = 0
      } else {                   // read potential in 4 corners
        get_4corner_potentials(gpos_ll, gpos_ur, corner_pot);
      }

      delta = delta * MM2CM;    // delta in mm must be used in cm
      if (delta == 0.) {
        printf("Error: MM_PER_UNIT constant not set!\n");
      }

      upper_pot_derivative = -(corner_pot[2] - corner_pot[3]) * M2CM / delta;
      // potential derivative on top of rectangle

      lower_pot_derivative = -(corner_pot[1] - corner_pot[0]) * M2CM / delta;
      // potential derivative at bottom of rectangle

      left_pot_derivative = -(corner_pot[3] - corner_pot[0]) * M2CM / delta;
      // potential derivative left of rectangle

      right_pot_derivative = -(corner_pot[2] - corner_pot[1]) * M2CM / delta;
      // potential derivative right of rectangle

      deltax = grid_pos[0] - gpos_ll[0];
      deltay = grid_pos[1] - gpos_ll[1];

      evec[1] = left_pot_derivative * (1. - deltax)
          + right_pot_derivative * deltax;

      evec[0] = upper_pot_derivative * (1. - deltay)
          + lower_pot_derivative * deltay;

      evec[2] = 0.;

      double omega[3];
      double output[3];
      omega[0] = azimuth;
      omega[1] = 0.;
      omega[2] = 0.;
      rotate_vec(evec, omega, output);
      vector_times_scalar(output, 1., evec);

      //          evec[2] = sin(cylinder_angle)*evec[1]; // rotate back to cylinder angle
      //          evec[1] = cos(cylinder_angle)*evec[1]; // rotate back to cylinder angle
    } else {                     // if get_symmetry
      cout << endl << "ERROR: 3D maps not included yet!!!" << endl;
      cout << "Switching to old efield interpolation for 3D" << endl;
      cout << endl << endl << flush;
      efield3d_old(pos_cm, evec, e_tag);
    }
  } else                        // if not enabled
    vector_times_scalar(grid_pos, 0., evec);  // evec = 0
}
