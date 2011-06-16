/* $Id: magfield3.h 235 2007-07-27 08:08:27Z s_voec01 $ */

#ifndef MAGFIELD3_H
#define MAGFIELD3_H

void magfield3_input_coils(const char *inputcoilfile);
void magfield3_test_coils();
void magfield3_source();
void magfield3_field(double *P, double *B);
void magfield3_field_elliptic(double *P, double *B);

/* Functions to modify the nmaxmag parameter */
void magfield3_set_nmaxmag(int n);
int magfield3_get_nmaxmag(void);

/* Functions to control the outout files */
void magfield3_set_magcoil(const char* filename);
void magfield3_set_magsource_central(const char* filename);
void magfield3_set_magsource_remote(const char* filename);
void magfield3_set_magsource_axisymm(const char* filename);
const char* magfield3_get_magcoil(void);
const char* magfield3_get_magsource_central(void);
const char* magfield3_get_magsource_remote(void);
const char* magfield3_get_magsource_axisymm(void);

/* A shortcut to set a prefix for all output files */
void magfield3_set_prefix(const char* prefix);

#endif
