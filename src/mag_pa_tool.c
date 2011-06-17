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
#include "bfield_ferenc.h"
#include "magfield3.h"

static int nx, ny, nz, sym;
static int max_nx = 20000;      // xline length to read entire lines in pa's
static void *mag_potential_array;

int get_magpa_sym()
{
  return sym;
}

void alloc_mag_arrays()
{
  printf("\nAllocating memory for magnetic potential array.\n");
  mag_potential_array = malloc(N_ARRAY * sizeof(double));
  if (mag_potential_array == 0) {
    printf("\nERROR: not enough memory for N_POT_ARRAY of %d!\n\n", N_ARRAY);
  } else {
    printf("                                      DONE.\n");
  }
}


void free_mag_arrays()
{
  free(mag_potential_array);
  return;
}


int read_mag_pa(char *mag_filename)
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

  strcpy(pa_in_name, mag_filename);
  strcat(pa_in_name, ".pa");    // data input from file with extension .pa

  pa_in = fopen(pa_in_name, "r");	//* open PA-File in read mode */

  if (pa_in == (FILE *) 0) {     // disk error handling
    fprintf(stderr, "\nERROR: Can't read file %s\n\n", pa_in_name);
    error = 1;
  } else {
    printf("\nReading magnetic potential array: %s", pa_in_name);
    fread(&h, sizeof(h), 1, pa_in);	/* read the header out of the file */

    sym = h.symmetry;
    printf("\nHeader information:\n");
    printf("  mode %d\n", h.mode);
    printf("  symmetrie %d\n", sym);
    printf("  max_voltage %f\n", h.max_voltage);
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
        printf("ERROR: read_mag_pa: max_nx lower than nx!\n\n");
      } else {
        for (z = 0; z < h.nz; z++) {
          for (y = 0; y < h.ny; y++) {

            cout << "  reading: " << (int) (100. * (float) (y + 1) * (float) (z + 1) / (float) (ny) / (float) (nz)) << " %\r" << flush;

            fread(pot_xline, sizeof(double), nx, pa_in);

            for (x = 0; x < h.nx; x++) {
              //                                                                                                                                                                                                                                                             read(pa_in, &pot, sizeof(double));                                                                                                                                                                                                                             /* read potential of an (x/y/z)-Koord. */

              if (pot_xline[x] > h.max_voltage) {
                pot_xline[x] = pot_xline[x] - h.max_voltage;
              }
              write_matrix_value(mag_potential_array, h.nx, h.ny, h.nz, 0, x, y, z, pot_xline[x]);
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


double get_magpot3D(double *position_cm)
{
  double oktett[8];             // potential values at all 8 corners
  double pos_backup[3];         // real position without any mirroring
  double position[3];           // position used for epot calculation
  double magpot_interpol;       // interpolated potential value
  int pos_left[3];              // position of lower left corner in front - position of oktett[0]
  int pos_right[3];             // position of upper right corner in back - position of oktett[5]
  int comp;

  vector_times_scalar(position_cm, 1., position);	// copy position vector
  position[0] = (position[0] - MAG_X_OFFSET_IN_CM) * CM2MM / MAG_MM_PER_UNIT;	// transform to simion grid
  position[1] = (position[1] - MAG_Y_OFFSET_IN_CM) * CM2MM / MAG_MM_PER_UNIT;
  position[2] = (position[2] - MAG_Z_OFFSET_IN_CM) * CM2MM / MAG_MM_PER_UNIT;
  vector_times_scalar(position, 1., pos_backup);	// backup position vector

  for (comp = 0; comp <= 2; comp++) {
    pos_left[comp] = (int) (position[comp] + 100000.) - 100000;	// corner 0
    // 10000 shifts all to positive than truncate and shift back
    // this is only needed when negative position
    if (pos_left[comp] < 0.) {
      pos_left[comp] = -pos_left[comp];	// use positive value
      pos_right[comp] = pos_left[comp] - 1;	// corner 5
    } else {
      pos_right[comp] = pos_left[comp] + 1;	// corner 5
    }
    if (position[comp] < 0.) {
      position[comp] = -position[comp];	// use positive values
    }
  }
  if ((position[0] > nx) || (position[1] > ny) || (position[2] > nz)) {
    printf("ERROR: in get_magpot3D: coordinates outside BFIELD_3D potential map!\n");

    printf("%f %f %f\n", position[0], position[1], position[2]);
    printf("%d %d %d\n", nx, ny, nz);
    return -1.;                 // watch for electrodes
  } else {                      // read potential in all 8 corners of cubus
    oktett[0] = read_matrix_value(mag_potential_array, nx, ny, nz, 0, pos_left[0], pos_left[1], pos_left[2]);
    oktett[1] = read_matrix_value(mag_potential_array, nx, ny, nz, 0, pos_right[0], pos_left[1], pos_left[2]);
    oktett[2] = read_matrix_value(mag_potential_array, nx, ny, nz, 0, pos_right[0], pos_right[1], pos_left[2]);
    oktett[3] = read_matrix_value(mag_potential_array, nx, ny, nz, 0, pos_left[0], pos_right[1], pos_left[2]);
    oktett[4] = read_matrix_value(mag_potential_array, nx, ny, nz, 0, pos_left[0], pos_right[1], pos_right[2]);
    oktett[5] = read_matrix_value(mag_potential_array, nx, ny, nz, 0, pos_right[0], pos_right[1], pos_right[2]);
    oktett[6] = read_matrix_value(mag_potential_array, nx, ny, nz, 0, pos_right[0], pos_left[1], pos_right[2]);
    oktett[7] = read_matrix_value(mag_potential_array, nx, ny, nz, 0, pos_left[0], pos_left[1], pos_right[2]);

    for (comp = 0; comp <= 2; comp++) {	// redo corner 0 and 5 for interpolation
      pos_left[comp] = (int) (pos_backup[comp] + 100000.) - 100000;
      // now positive and negative values needed
      pos_right[comp] = pos_left[comp] + 1;
    }
    magpot_interpol = interpol_3dim(pos_backup, pos_left, pos_right, oktett);
    // do interpolation
    return magpot_interpol;     // return epot value
  }
}

double get_magpot2D(double *position_cm)
{
  double oktett[8];             // potential values at all 8 corners
  double pos_backup[3];         // real position without any mirroring
  double position[3];           // position used for epot calculation
  double magpot_interpol;       // interpolated potential value
  int pos_left[3];              // position of lower left corner in front - position of oktett[0]
  int pos_right[3];             // position of upper right corner in back - position of oktett[5]
  int comp;

  vector_times_scalar(position_cm, 1., position);	// copy position vector
  position[0] = (position[0] - MAG_X_OFFSET_IN_CM) * CM2MM / MAG_MM_PER_UNIT;	// transform to simion grid
  position[1] = (position[1] - MAG_Y_OFFSET_IN_CM) * CM2MM / MAG_MM_PER_UNIT;
  position[2] = (position[2] - MAG_Z_OFFSET_IN_CM) * CM2MM / MAG_MM_PER_UNIT;
  position[1] = sqrt(position[1] * position[1] + position[2] * position[2]);
  position[2] = 0.;
  vector_times_scalar(position, 1., pos_backup);	// backup position vector

  for (comp = 0; comp <= 2; comp++) {
    pos_left[comp] = (int) (position[comp] + 100000.) - 100000;	// corner 0
    // 5000 shifts all to positive than truncate and shift back
    // this is only needed when negative position
    if (pos_left[comp] < 0.) {
      pos_left[comp] = -pos_left[comp];	// use positive value
      pos_right[comp] = pos_left[comp] - 1;	// corner 5
    } else {
      pos_right[comp] = pos_left[comp] + 1;	// corner 5
    }
    if (position[comp] < 0.) {
      position[comp] = -position[comp];	// use positive values
    }
  }
  if ((position[0] > nx) || (position[1] > ny) || (position[2] > nz)) {
    printf("ERROR: in get_magpot2D: coordinates outside BFIELD_3D potential map!\n");
    return -1.;                 // watch for electrodes
  } else {                      // read potential in all 8 corners of cubus
    oktett[0] = read_matrix_value(mag_potential_array, nx, ny, nz, 0, pos_left[0], pos_left[1], pos_left[2]);
    oktett[1] = read_matrix_value(mag_potential_array, nx, ny, nz, 0, pos_right[0], pos_left[1], pos_left[2]);
    oktett[2] = read_matrix_value(mag_potential_array, nx, ny, nz, 0, pos_right[0], pos_right[1], pos_left[2]);
    oktett[3] = read_matrix_value(mag_potential_array, nx, ny, nz, 0, pos_left[0], pos_right[1], pos_left[2]);
    oktett[4] = 0.;
    oktett[5] = 0.;
    oktett[6] = 0.;
    oktett[7] = 0.;

    for (comp = 0; comp <= 2; comp++) {	// redo corner 0 and 5 for interpolation
      pos_left[comp] = (int) (pos_backup[comp] + 100000.) - 100000;
      // now positive and negative values needed
      pos_right[comp] = pos_left[comp] + 1;
    }
    magpot_interpol = interpol_3dim(pos_backup, pos_left, pos_right, oktett);
    // do interpolation
    return magpot_interpol;     // return epot value
  }
}



void get_bfield_poles(double *position, double *bvec)
{
  double xscale = parameter.b_field_ben1;	// default 10.0 ;// in cm
  double bmax = parameter.b_field_ben2;	// default 1.0; //in tesla 
  bvec[0] = bmax * (1 + pow(position[0] / xscale, 2));
  bvec[1] = -bmax * position[1] * position[0] / (xscale * xscale);
  bvec[2] = -bmax * position[2] * position[0] / (xscale * xscale);
}


void get_bfield_std(double *position, double *bvec)
{
  double pos1[3], pos3[3];
  double delta = MAG_MM_PER_UNIT;
  int komp;

  delta = delta * MM2CM;        // delta in mm must be used in cm
  if (delta == 0.) {
    printf("Error: MAG_MM_PER_UNIT constant not set!\n");
  }

  if (get_magpa_sym() == 1) {
    for (komp = 0; komp <= 2; komp++) {
      vector_times_scalar(position, 1., pos1);
      vector_times_scalar(position, 1., pos3);
      pos1[komp] = pos1[komp] - delta;
      pos3[komp] = pos3[komp] + delta;
      bvec[komp] = (get_magpot3D(pos3) - get_magpot3D(pos1))
          * M2CM * MAG_MM_PER_UNIT * MM2CM / (2. * TESLA2GAUSS * delta);
    }
  } else {
    for (komp = 0; komp <= 2; komp++) {
      vector_times_scalar(position, 1., pos1);
      vector_times_scalar(position, 1., pos3);
      pos1[komp] = pos1[komp] - delta;
      pos3[komp] = pos3[komp] + delta;
      bvec[komp] = (get_magpot2D(pos3) - get_magpot2D(pos1))
          * M2CM * MAG_MM_PER_UNIT * MM2CM / (2. * TESLA2GAUSS * delta);
    }
  }
}



void get_bfield_homogen(double *bvec)
{
  bvec[0] = 1.0;
  bvec[1] = 0.;
  bvec[2] = 0.;
}


/************** b_toroid ******************/
// this function calculates a toroidal magnetic field
void b_toroid(double *position, double *bvec)
{
  double pos[3];                //position for calculation
  double radius_value;          //absolute value of magnetic field radius
  double radius_vec[3];         //magnetic field radius
  double r_0 = 100.;            //radius at b=b_0 in cm
  double z_vec[3] = { 0., 0., 1. };	//unit vector in z direction
  double b_0 = 1.;              //magnetic field in Tesla for r_0=100cm
  double b_value;               //magnetic field in Tesla
  double xprod[3];              //r x z

  vector_times_scalar(position, 1., pos);
  pos[2] = 0.;

  radius_value = absvalue(pos);
  vector_times_scalar(pos, (1. / radius_value), radius_vec);

  b_value = r_0 * b_0 / radius_value;

  cross_prod(radius_vec, z_vec, xprod);

  vector_times_scalar(xprod, b_0 * r_0 / radius_value, bvec);
}

void get_bfield_mag3(double *pos_ap, double *bvec_ap)
{
  //convert from cm to m, rotate x to z, use components 1-3
  double pos_m3[4];
  double bvec_m3[4];
  char mag_filename[255];       // filename for magsource data
  strcpy(mag_filename, parameter.filename);
  strcat(mag_filename, "_mag3parms"); // data input from soucepoints file
  magfield3_set_prefix(mag_filename);

  pos_m3[1] = pos_ap[1]/100.0;
  pos_m3[2] = pos_ap[2]/100.0;
  pos_m3[3] = pos_ap[0]/100.0;//symmetry axis of coils

  magfield3_field(pos_m3, bvec_m3);
  // call magfield of Ferenc Glck with converted values

  bvec_ap[1] = bvec_m3[1];
  bvec_ap[2] = bvec_m3[2];
  bvec_ap[0] = bvec_m3[3];//symmetry axis of coils
  //cout << "AP Pos " << pos_ap[0] << " " << pos_ap[1] << " " << pos_ap[2] << endl << flush;
  //cout << "AP BVec " << bvec_ap[0] << " " << bvec_ap[1] << " " << bvec_ap[2] << endl << flush;
}

void get_bfield(double *position, double *bvec)
{
  switch (USE_MAG_PA) {
    case 0:
      printf("ERROR: This version of ADIPARK does not support magnetic field calculation!\n");
    case 1:
      get_bfield_std(position, bvec);
      break;
    case 2:
      get_bfield_homogen(bvec);
      break;
    case 3:
      b_toroid(position, bvec);
      break;
    case 4:
      get_bfield_fg(position, bvec);
      break;
    case 5:
      get_bfield_poles(position, bvec);
      break;
    case 6:
      get_bfield_mag3(position, bvec);
      break;
    default:
      printf("ERROR: In mag_pa_reader/get_bfield: PA use not set properly!\n");
  }
}
