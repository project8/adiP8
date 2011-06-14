/* $Id: magfield2.h 231 2007-07-26 11:43:35Z s_voec01 $ */

#ifndef MAGFIELD2_H_
#define MAGFIELD2_H_

void magfield2_source(double z0min, double z0max, double delz0, const char *inputcoil);

void magfield2_field(double z, double r, const char *inputcoil, int n, double *A,
               double *Bz, double *Br);

void magfield2_set_nmaxmag(int n);
int magfield2_get_nmaxmag(void);
void magfield2_set_magsource(const char* filename);
const char* magfield2_get_magsource();
static int Ncoil, Nspmag;
static double **coil = NULL;
#endif /*MAGFIELD2_H_*/
