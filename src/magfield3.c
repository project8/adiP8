/* $Id: magfield3.c 235 2007-07-27 08:08:27Z s_voec01 $ */
////////////////////////////////////////////////////////////////////////
//                                                                    //
//     magfield3.c :  magnetic field calculation                      //
//                    with 3-dimensional coil system                  //
//                                                                    //
////////////////////////////////////////////////////////////////////////

#include "magfield3.h"
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "array.h"

#define FMIN(a,b) ((a)<=(b)?(a):(b))
#define FMAX(a,b) ((a)>(b)?(a):(b))
#define FMIN3(a,b,c) (FMIN(a,b)<=(c)?(FMIN(a,b)):(c))
#define FMAX3(a,b,c) (FMAX(a,b)>(c)?(FMAX(a,b)):(c))
#define pow2(x) ((x)*(x))

// Part 1: input of 3-dimensional coil parameters, 
//         test of coil parameters,
//         calculation of source coefficients,
//         and magnetic field calculation of all coils.

// Part 2: calculation in local axisymmetric systems of the coils

static void magfield2_elliptic_1coil(int n, double *coilpar, double z, double r,
                                     double *Bz, double *Br);
static double funrocen(double *coilpar, double z0);
static void magsource2_central_1coil(double *coilpar, double z0, double rocen,
                                     double *Bcen1);
static void magfield2_central(int type, int s, double z0, double rocen,
                              double z, double r, double *Bz, double *Br,
                              double *rcp);
static double funrorem(double *coilpar, double z0);
static void magsource2_remote_1coil(double *coilpar, double z0, double rorem,
                                    double *Brem1);
static void magfield2_remote(int type, int s, double z0, double rorem, double z,
                             double r, double *Bz, double *Br, double *rcp);

// Part 3: calculation of  3-dimensional coils

static void magfield3_elliptic_1coil(int i, double *P, double *B);
static void magsource3_central();
static void magsource3_remote();
static void magfield3_1coil(int i, double *P, double *B);

// Part 4: calculation of the axisymmetric coils

static void magsource_axisymm();
static double funrocenaxisymm(double z0cen);
static void magfield_axisymm(double *P, double *B);

// Part 5: complete elliptic integral calculations
//           (according to Numerical Recipes)

static double RF_Carlson(double x, double y, double z);
static double RD_Carlson(double x, double y, double z);
static double RJ_Carlson(double x, double y, double z, double p);
static double RC_Carlson(double x, double y);

static int Ncoil = 0, Ncoilaxisymm = 0;
static double **coil = NULL;
static double **Bcen = NULL, **Brem = NULL;
static double **Bcenaxisymm = NULL;
static double *Bremaxisymm = NULL;
static int **indexaxisymm = NULL;
static int lastNsp;

static int nmaxmag = 200;

/* Output files */
static char *filename_magcoil = NULL;
static char *filename_magsource_central = NULL;
static char *filename_magsource_remote = NULL;
static char *filename_magsource_axisymm = NULL;

static void malloc_failed()
{
    fprintf(stderr, "Error: Not enough memory!\n");
    abort();
}

void allocate_arrays()
{
    int i;

    if (Ncoil > 0) {
        coil = (double**) malloc(sizeof(double *) * (Ncoil + 1));
        indexaxisymm = (int**) malloc(sizeof(int *) * (Ncoil + 1));
        Brem = (double**) malloc(sizeof(double) * (3 * Ncoil + 1));
        for (i = 1; i <= Ncoil; i++) {
            coil[i] = (double*) malloc(sizeof(double) * 14);
            indexaxisymm[i] = (int*) calloc(2, sizeof(int));
        }
        for (i = 1; i <= 3 * Ncoil; i++)
            Brem[i] = (double*) malloc(sizeof(double) * (nmaxmag + 1));
    }
    if (Bremaxisymm)
        free(Bremaxisymm);
    Bremaxisymm = (double*) malloc(sizeof(double) * (nmaxmag + 1));
}

void allocate_arraysNsp(int Nsp)
{
    int i;

    lastNsp = Nsp;
    
    printf("Nsp: %d\n", Nsp);

    if (Nsp > 0) {
        Bcenaxisymm = (double**) malloc(sizeof(double *) * (Nsp + 1));
        if (!Bcenaxisymm)
            malloc_failed();
        for (i = 1; i <= Nsp; i++) {
            Bcenaxisymm[i] = (double*) malloc(sizeof(double) * (nmaxmag + 1));
            if (!Bcenaxisymm[i])
                malloc_failed();
        }
    }
}

void free_arrays(int oldNcoil)
{
    int i;

    if (coil && indexaxisymm) {
        for (i = 0; i < oldNcoil; i++) {
            free(coil[i]);
            free(indexaxisymm[i]);
        }
        free(coil);
        free(indexaxisymm);
    }
}

void free_arraysNsp()
{
    int i;

    if (Bcenaxisymm) {
        for (i = 1; i <= lastNsp; i++)
            free(Bcenaxisymm[i]);
        free(Bcenaxisymm);
        Bcenaxisymm = NULL;
    }
}

// Part 1: input of 3-dimensional coil parameters, and
//         calculation of source coefficients

void magfield3_input_coils(const char *inputcoilfile)
{
// This function reads the coil parameters (current and geometry)
//  of a general 3-dimensional coil system,
// from the data file defined by the string inputcoilfile.
// The coils have rectangular shape cross section with various
//  symmetry axes.
// The data in the file defined by string inputcoilfile are:
//     First line: number of coils  (Ncoil).
//     Then there are Ncoil number of lines; each line contains:
//       cd xA yA zA xB yB zB Rmin Rmax n
//    cd: current density (A/m^2)
//    (the current density is positive if looking from the 
//     point A towards point B
//     the current flow direction is clockwise)
//    xA:  x component of endpoint A of coil axis (m) 
//    yA:  y component of endpoint A of coil axis (m) 
//    zA:  z component of endpoint A of coil axis (m)  
//    xB:  x component of endpoint B of coil axis (m)   
//    yB:  y component of endpoint B of coil axis (m)  
//    zB:  z component of endpoint B of coil axis (m)  
//    Rmin:  inner radius of coil   (m)
//    Rmax:  outer radius of coil  (m)
//    n: radial integration number of coil
// The coil parameters are written into the array coil[Ncoilmax+1][14].
//   i: coil index
//   coil[i][0]: current density (A/m^2)
//   coil[i][1]: x component of endpoint A of coil axis (m) 
//   coil[i][2]: y component of endpoint A of coil axis (m) 
//   coil[i][3]: z component of endpoint A of coil axis (m) 
//   coil[i][4]: x component of endpoint B of coil axis (m) 
//   coil[i][5]: y component of endpoint B of coil axis (m) 
//   coil[i][6]: z component of endpoint B of coil axis (m) 
//   coil[i][7]: inner radius of coil  (m)
//   coil[i][8]: outer radius of coil (m)
//   coil[i][9]: radial integration number of coil
//      (the integer is represented here by double).
//  SI units are used here!
// The user of the program can of course arbitrarily define the coil
// input file, and accordingly to use another (self-written) function
// that reads the coil parameters from the coil input file and
// computes the coil[i][0]-coil[i][9] coil parameters.
// It is important that the user should not change the definitions
// of the coil[i][0]-coil[i][9] parameters (otherwise one would have to
// rewrite the whole program package).
    int i, n;
    double cd, xA, yA, zA, xB, yB, zB, Rmin, Rmax;
    FILE *fp;

//
    fp = fopen(inputcoilfile, "r");
    if (!fp) {
        fprintf(stderr, "Message from function input_coils:\n");
        fprintf(stderr, "Cannot open the coil input file %s!\n", inputcoilfile);
        fprintf(stderr, "Program running is stopped !!!\n");
        exit(1);
    }
    fscanf(fp, "%i", &Ncoil);

    allocate_arrays();

    for (i = 1; i <= Ncoil; i++) {
        fscanf(fp, "%le %le %le %le %le %le %le %le %le %i ", &cd, &xA, &yA,
               &zA, &xB, &yB, &zB, &Rmin, &Rmax, &n);
// Coil parameter determination:
        coil[i][0] = cd;
        coil[i][1] = xA;
        coil[i][2] = yA;
        coil[i][3] = zA;
        coil[i][4] = xB;
        coil[i][5] = yB;
        coil[i][6] = zB;
        coil[i][7] = Rmin;
        coil[i][8] = Rmax;
        coil[i][9] = n;
    }
    fclose(fp);
}

void magfield3_test_coils()
{
// This function tests the coil parameters. It computes also the
// additional (10-13) coil parameters, and the indexaxisymm[ia][0],
// indexaxisymm[in][1] integers.
// The coil parameters are written into the data file 'magcoil.dat'.
// Additional (redundant) coil parameters: 
//     coil[i][10]: length of coil (distance of points A and B) (m)
//     coil[i][11],coil[i][12],coil[i][13]: unit vector in A --> B
//                                          direction (m)
// Axisymmetric coils (if Ncoilaxisymm>0):
//    index: ia=1,...,Ncoilaxisymm
//   The coil index i corresponding to ia:  i=indexaxisymm[ia][0]
// Non-axisymmetric coils (if Ncoilaxisymm<Ncoil):
//    index: in=1,...,Ncoil-Ncoilaxisymm
//   The coil index i corresponding to in:  i=indexaxisymm[in][1]
    int i, k, ia, in;
    double v[4], stock;
    FILE *fp;

// Test of coil parameters:     
    if (Ncoil < 1) {
        fprintf(stderr, "Message from function test_coils:\n"
                "Ncoil<1 !!!\n" "Program running is stopped !!!\n");
        exit(1);
    }
    for (i = 1; i <= Ncoil; i++) {
        if (coil[i][7] <= 1.e-9 || coil[i][8] <= 1.e-9 || coil[i][9] < 0.9) {
            fprintf(stderr, "Message from function test_coils:\n"
                    "non-positive coil parameters: Rmin, Rmax or n !!!\n"
                    "Program running is stopped !!! i= %9i\n", i);
            exit(1);
        }
        if (coil[i][7] >= coil[i][8]) {
            fprintf(stderr, "Message from function test_coils:\n"
                    "Rmin>=Rmax!!! Rmax is changed to Rmin+1.e-6!!! i= %9i\n",
                    i);
            coil[i][8] = coil[i][7] + 1.e-6;
        }
    }
//------------------------------------------------------ 
// Additional coil parameters, and test of coil length:
    for (i = 1; i <= Ncoil; i++) {
        for (k = 1; k <= 3; k++)
            v[k] = coil[i][k + 3] - coil[i][k];
        coil[i][10] = sqrt(v[1] * v[1] + v[2] * v[2] + v[3] * v[3]);
        if (coil[i][10] <= 1.e-12) {
            fprintf(stderr, "Message from function test_coils:\n"
                    "length of coil is too small !!!\n"
                    "Program running is stopped !!! i= %9i\n", i);
            exit(1);
        }
        for (k = 1; k <= 3; k++)
            coil[i][10 + k] = v[k] / coil[i][10];
    }
//------------------------------------------------------ 
// Calculation of coil axisymmetry indices 
//   indexaxisymm[ia][0], indexaxisymm[in][1]:  
    /*for(i=1;i<=Ncoil;i++)
       {
       indexaxisymm[i][0]=0;
       indexaxisymm[i][1]=0;
       } */
    Ncoilaxisymm = 0;
    ia = 0;
    in = 0;
    for (i = 1; i <= Ncoil; i++) {
        if (fabs(coil[i][1]) < 1.e-7 && fabs(coil[i][2]) < 1.e-7
            && fabs(coil[i][4]) < 1.e-7 && fabs(coil[i][5]) < 1.e-7) {
            ia += 1;
            indexaxisymm[ia][0] = i;
            coil[i][1] = coil[i][2] = coil[i][4] = coil[i][5] = 0.;
            Ncoilaxisymm += 1;
// A --> B direction is changed into +z direction:
            if (coil[i][3] > coil[i][6]) {
                stock = coil[i][3];
                coil[i][3] = coil[i][6];
                coil[i][6] = stock;
                coil[i][0] = -coil[i][0];
            }
            coil[i][10] = fabs(coil[i][6] - coil[i][3]);
            coil[i][11] = 0.;
            coil[i][12] = 0.;
            coil[i][13] = 1.;
            if (fabs(coil[i][3] - coil[i][6]) < 1.e-12) {
                puts("Warning message from function test_coils:");
                puts("Axisymmetric coil length is zero !!!");
                puts("Program running is stopped !!! ");
                exit(0);
            }
        } else {
            in += 1;
            indexaxisymm[in][1] = i;
        }
    }
// Output into data file 'magcoil.dat':
    /* fp=fopen("magcoil.dat","w"); */
    fp = fopen(magfield3_get_magcoil(), "w");
// We write the coil parameters:  
    fprintf(fp, "%9i %9i \n", Ncoil, Ncoilaxisymm);
    fprintf(fp, "%9i \n", nmaxmag);
    for (i = 1; i <= Ncoil; i++) {
        fprintf(fp, "%9i %9i %9i \n", i, indexaxisymm[i][0],
                indexaxisymm[i][1]);
        fprintf(fp, "%22.15e %22.15e %22.15e %22.15e  \n", coil[i][0],
                coil[i][1], coil[i][2], coil[i][3]);
        fprintf(fp, "%22.15e %22.15e %22.15e %22.15e %22.15e  \n", coil[i][4],
                coil[i][5], coil[i][6], coil[i][7], coil[i][8]);
        fprintf(fp, "%22.15e %22.15e %22.15e %22.15e %22.15e  \n", coil[i][9],
                coil[i][10], coil[i][11], coil[i][12], coil[i][13]);
    }
    fclose(fp);
}

void magfield3_source()
{
// It computes all the source coefficients that are needed
// for the central and remote Legendre polynomial calculations.
// The source coefficient data files 'magsource_central.dat',
// 'magsource_remote.dat' and 'magsource_axisymm.dat' are created.
// Central source coefficient calculation for all coils:
    magsource3_central();
// Remote source coefficient calculation for all coils:
    magsource3_remote();
// Source coefficient calculation for axisymmetric coils:
    if (Ncoilaxisymm > 0)
        magsource_axisymm();
}

void inputmagcoil()
{
// This function reads the coil parameters from the data file 'magcoil.dat'.
    int i, ix, nmaxmagtest;
    FILE *fp;

    free_arrays(Ncoil);
// Input from data file 'magcoil.dat':
    /* fp=fopen("magcoil.dat","r"); */
    fp = fopen(magfield3_get_magcoil(), "r");
    if (!fp) {
        fprintf(stderr, "Message from function inputmagcoil:\n");
        fprintf(stderr, "Cannot open the file '%s' !\n",
                magfield3_get_magcoil());
        fprintf(stderr, "Program running is stopped !!!\n");
        exit(1);
    }
// We read the coil parameters:  
    fscanf(fp, "%i %i \n", &Ncoil, &Ncoilaxisymm);
    if (Ncoil < 1) {
        fprintf(stderr, "Message from function inputmagcoil:\n"
                "Ncoil<1 !!!\n"
                "Computation is  stopped !!! Ncoil= %9i\n", Ncoil);
        exit(1);
    }
    fscanf(fp, "%i  \n", &nmaxmagtest);
    if (nmaxmag != nmaxmagtest) {
        fprintf(stderr, "Error message from function inputmagcoil:\n"
                "The value of nmaxmag has been changed !!!\n"
                "For the magnetic field calc. the same nmaxmag should be used as\n"
                "during the coil parameter and source coefficient calc. (old nmaxmag)!\n"
                "Computation is  stopped !!! old nmaxmag value= %9i\n",
                nmaxmagtest);
        exit(1);
    }

    allocate_arrays();

    for (i = 1; i <= Ncoil; i++) {
        fscanf(fp, "%i %i %i \n", &ix, &indexaxisymm[i][0],
               &indexaxisymm[i][1]);
        fscanf(fp, "%le %le %le %le  \n", &coil[i][0], &coil[i][1], &coil[i][2],
               &coil[i][3]);
        fscanf(fp, "%le %le %le %le %le  \n", &coil[i][4], &coil[i][5],
               &coil[i][6], &coil[i][7], &coil[i][8]);
        fscanf(fp, "%le %le %le %le %le  \n", &coil[i][9], &coil[i][10],
               &coil[i][11], &coil[i][12], &coil[i][13]);
    }
    fclose(fp);
}

void magfield3_field(double *P, double *B)
{
// This function computes the magnetic field components B[1],B[2],B[3]
// in a field point P[1],P[2],P[3], due to all coils,
// using central or remote Legendre polynomial expansion,
// or elliptic integrals.
// SI units are used (P[k] in m, B[k] in T, k=1,2,3)!
    int i, in, k;
    double Bi[4];
    static int iff = 0;

// At the first call of this function, the coil parameters from the
// data file 'magcoil.dat' are read.
    if (iff == 0) {
        inputmagcoil();
        iff = 1;
    }
// First we compute the axisymmetric coils (if Ncoilaxisymm>0).
    if (Ncoilaxisymm == 0)
        for (k = 1; k <= 3; k++)
            B[k] = 0.;
    else
        magfield_axisymm(P, B);
// Second, we compute the non-axisymmetric coils (if Ncoilaxisymm<Ncoil).
    if (Ncoilaxisymm == Ncoil)
        return;
    else {
        for (in = 1; in <= Ncoil - Ncoilaxisymm; in++) {
            i = indexaxisymm[in][1];
            magfield3_1coil(i, P, Bi);
            for (k = 1; k <= 3; k++)
                B[k] += Bi[k];
        }
    }
}

void magfield3_field_elliptic(double *P, double *B)
{
// This function computes the magnetic field components B[1],B[2],B[3]
// in a field point P[1],P[2],P[3], due to all coils,
// using elliptic integrals.
// SI units are used (P[k] in m, B[k] in T, k=1,2,3)!
    int i, k;
    double Bi[4];
    static int iff = 0;

// At the first call of this function, the coil parameters from the
// data file 'magcoil.dat' are read.
    if (iff == 0) {
        inputmagcoil();
        iff = 1;
    }
    for (k = 1; k <= 3; k++)
        B[k] = 0.;
    for (i = 1; i <= Ncoil; i++) {
        magfield3_elliptic_1coil(i, P, Bi);
        for (k = 1; k <= 3; k++)
            B[k] += Bi[k];
    }
}

// Part 2: axially symmetric coils

void magfield2_elliptic_1coil(int n, double *coilpar, double z, double r,
                              double *Bz, double *Br)
{
// This function computes the magnetic field components Bz and Br
// of an axially symmetric coil, with z axis as symmetry axis,
// in a fieldpoint with (z,r) cylindrical coordinates, using
// the first, second and third complete elliptic integrals.
// n: radial integration number
// coilpar[0]: current density of the coil.
// coilpar[1], coilpar[2]: Zmin and Zmax parameters of the coil.
// coilpar[3], coilpar[4]: Rmin and Rmax parameters of the coil. 
    double Zmin, Zmax, Rmin, Rmax, sigma;
    double R, delR[3], Z, delr2, sumr2, delz2, eta, d, K, EK, PIK, S;
    double xBz, xBr, sign, c, st, delRr;
    double Rlow[3], Rhigh[3];
    const double mu0 = 4. * M_PI * 1.e-7;
    int i, iR, iZ, M, m;
    double w[1001];
    const double w1[2] = { 0.5, 0.5 };
    const double w2[3] = { 1. / 3., 4. / 3., 1. / 3. };
    const double w4[5] =
        { 14. / 45., 64. / 45., 24. / 45., 64. / 45., 14. / 45. };
    const double w5[6] = { 0.3187500000000000e+00, 0.1376388888888889e+01,
        0.6555555555555556e+00, 0.1212500000000000e+01,
        0.9256944444444445e+00, 0.1011111111111111e+01
    };
    const double w9[10] = { 0.2803440531305107e0, 0.1648702325837748e1,
        -0.2027449845679092e0, 0.2797927414021179e1,
        -0.9761199294532843e0, 0.2556499393738999e1,
        0.1451083002645404e0, 0.1311227127425048e1,
        0.9324249063051143e0, 0.1006631393298060e1
    };
// Coil parameters:
    sigma = coilpar[0];
    Zmin = coilpar[1];
    Zmax = coilpar[2];
    Rmin = coilpar[3];
    Rmax = coilpar[4];
// Improvement of Zmin,Zmax,Rmin,Rmax,z,r values:
    if (Zmax < Zmin) {
        st = Zmax;
        Zmax = Zmin;
        Zmin = st;
    }
    if (Rmax < Rmin) {
        st = Rmax;
        Rmax = Rmin;
        Rmin = st;
    }
    if (Zmax - Zmin < 1.e-12)
        Zmax = Zmin + 1.e-12;
    if (Rmax - Rmin < 1.e-12)
        Rmax = Rmin + 1.e-12;
// Coil geometry test:
    if (Rmin <= 1.e-9) {
        fprintf(stderr, "Message from function magfield2_elliptic_1coil:\n"
                "Rmin<=1.e-9!!!\n" "Computation is stopped !!!\n");
        exit(1);
    }
    if (fabs(z - Zmin) < 1.e-8 && r >= Rmin - (Rmax - Rmin) * 1.e-8 &&
        r <= Rmax + (Rmax - Rmin) * 1.e-8)
        z = Zmin - 1.e-8;
    if (fabs(z - Zmax) < 1.e-8 && r >= Rmin - (Rmax - Rmin) * 1.e-8 &&
        r <= Rmax + (Rmax - Rmin) * 1.e-8)
        z = Zmax + 1.e-8;
// Improvement of n value:
    if (n < 1)
        n = 1;
    else if (n == 3)
        n = 2;
    else if (n > 3 && n <= 6)
        n = 4;
    else if (n > 6 && n < 12)
        n = 12;
    else if (n > 1000)
        n = 1000;
// Integration weight factors:
    if (n == 1)
        for (i = 0; i <= 1; i++)
            w[i] = w1[i];
    else if (n == 2)
        for (i = 0; i <= 2; i++)
            w[i] = w2[i];
    else if (n == 4)
        for (i = 0; i <= 4; i++)
            w[i] = w4[i];
    else if (n >= 12 && n < 20) {
        for (i = 0; i <= 5; i++)
            w[i] = w5[i];
        for (i = 6; i <= n - 6; i++)
            w[i] = 1.;
        for (i = n - 5; i <= n; i++)
            w[i] = w5[n - i];
    } else {
        for (i = 0; i <= 9; i++)
            w[i] = w9[i];
        for (i = 10; i <= n - 10; i++)
            w[i] = 1.;
        for (i = n - 9; i <= n; i++)
            w[i] = w9[n - i];
    }
//
    xBz = 0.;
    xBr = 0.;
// R-integration limits:
    if (z > Zmin && z < Zmax && r > Rmin && r < Rmax) {
        M = 2;
        Rlow[1] = Rmin;
        Rhigh[1] = r - (r - Rmin) * 1.e-12;
        Rlow[2] = r + (Rmax - r) * 1.e-12;
        Rhigh[2] = Rmax;
        delR[1] = (Rhigh[1] - Rlow[1]) / n;
        delR[2] = (Rhigh[2] - Rlow[2]) / n;
    } else {
        M = 1;
        Rlow[1] = Rmin;
        Rhigh[1] = Rmax;
        delR[1] = (Rhigh[1] - Rlow[1]) / n;
    }
// Integration:
    for (m = 1; m <= M; m++) {
        for (iR = 0; iR <= n; iR++) {
            R = Rlow[m] + delR[m] * iR;
            for (iZ = 1; iZ <= 2; iZ++) {
                if (iZ == 1) {
                    Z = Zmax;
                    sign = 1.;
                } else {
                    Z = Zmin;
                    sign = -1.;
                }
                delr2 = (r - R) * (r - R);
                delz2 = (z - Z) * (z - Z);
                sumr2 = (r + R) * (r + R);
                d = delr2 / sumr2;
                eta = (delr2 + delz2) / (sumr2 + delz2);
                S = sqrt(sumr2 + delz2);
                K = RF_Carlson(0., eta, 1.);
                EK = -1. / 3. * RD_Carlson(0., eta, 1.);
                delRr = R - r;
                if (d < 1.e-18) {
                    d = 1.e-18;
                    if (R > r)
                        delRr = (r + R) * 1.e-9;
                    else if (R < r)
                        delRr = -(r + R) * 1.e-9;
                }
                PIK = 1. / 3. * RJ_Carlson(0., eta, 1., d);
                xBz +=
                    -w[iR] * sign * R * (z - Z) / (S * (R + r)) * (K +
                                                                   delRr / (2. *
                                                                            R) *
                                                                   PIK * (1. -
                                                                          d)) *
                    delR[m];
                xBr += -w[iR] * sign * R / S * (2. * EK + K) * delR[m];
            }
        }
    }
    c = mu0 / M_PI * sigma;
    *Bz = c * xBz;
    *Br = c * xBr;
}

// Central convergence radius calculation:

double funrocen(double *coilpar, double z0)
// This function computes the central convergence radius
// funrocen=rocen for an axisymmetric coil, at the axis source point z0.
// rocen = minimal distance of the axis point z0 from the coil winding.
// The coil parameters are given by *coilpar:
// coilpar[0]: current density
// coilpar[1]: Zmin, coilpar[2]: Zmax,
// coilpar[3]: Rmin, coilpar[4]: Rmax.
{
    double ro, Zmin, Zmax, Rmin;

    Zmin = coilpar[1];
    Zmax = coilpar[2];
    Rmin = coilpar[3];
    if (z0 <= Zmin)
        ro = sqrt((z0 - Zmin) * (z0 - Zmin) + Rmin * Rmin);
    else if (z0 >= Zmax)
        ro = sqrt((z0 - Zmax) * (z0 - Zmax) + Rmin * Rmin);
    else
        ro = Rmin;
    return ro;
}

// Central source coefficient calculation for 1 coil:
//  (axisymm. case)

void magsource2_central_1coil(double *coilpar, double z0, double rocen,
                              double *Bcen1)
// This subroutine computes the magnetic field central source constants Bcen1[n]
//  (n=0,...,nmaxmag) of an axisymmetric coil, in the axis source point z0,
//   with central convergence radius rocen as input parameter.
// The coil parameters are given by *coilpar:
// coilpar[0]: current density
// coilpar[1]: Zmin, coilpar[2]: Zmax,
// coilpar[3]: Rmin, coilpar[4]: Rmax.
// nmaxmag is defined by #define command.
// Indices: m: numerical integration;
//          n: source constants
//   iz=1: Z=Zmax; iz=2: Z=Zmin
// Radial integration number: M; this depends on the radial thickness:
//  M is small if the(Rmax-Rmin)/Rmin ratio is small, and it is large
//   if this ratio is large.
{
    double b[nmaxmag + 1][3], bhat[nmaxmag + 1][1001], Pp[nmaxmag + 2];
    int m, n, iz, M;
    const double w9[10] = { 0.2803440531305107e0, 0.1648702325837748e1,
        -0.2027449845679092e0, 0.2797927414021179e1,
        -0.9761199294532843e0, 0.2556499393738999e1,
        0.1451083002645404e0, 0.1311227127425048e1,
        0.9324249063051143e0, 0.1006631393298060e1
    };
    const double mu0 = 4. * M_PI * 1.e-7;
    static double w[1001];
    double ro, Z, R, sigma, del, q, st;
    double Zmin, Zmax, Rmin, Rmax, z, u, constcen, sigmac, rcen, rcen1, ratio;

// Coil parameters:
    sigma = coilpar[0];
    Zmin = coilpar[1];
    Zmax = coilpar[2];
    Rmin = coilpar[3];
    Rmax = coilpar[4];
// Improvement of Zmin,Zmax,Rmin,Rmax values:
    if (Zmax < Zmin) {
        st = Zmax;
        Zmax = Zmin;
        Zmin = st;
    }
    if (Rmax < Rmin) {
        st = Rmax;
        Rmax = Rmin;
        Rmin = st;
    }
    if (Zmax - Zmin < 1.e-12)
        Zmax = Zmin + 1.e-12;
    if (Rmax - Rmin < 1.e-12)
        Rmax = Rmin + 1.e-12;
// Coil geometry test:
    if (Rmin <= 1.e-9) {
        fprintf(stderr, "Message from function magsource2_central_1coil:\n"
                "Rmin is not positive!!!\n" "Computation is stopped !!!\n");
        exit(1);
    }
// Radial integration number M:
    ratio = (Rmax - Rmin) / Rmin;
    if (ratio < 0.1)
        M = 30;
    else if (ratio >= 0.1 && ratio < 0.2)
        M = 50;
    else
        M = 60 * ratio / 0.2;
    if (M > 1000)
        M = 1000;
// Initialization of the integration weight factors:
    for (m = 0; m <= 9; m++)
        w[m] = w9[m];
    for (m = 10; m <= M - 10; m++)
        w[m] = 1.;
    for (m = M - 9; m <= M; m++)
        w[m] = w9[M - m];
// Initialization of Bcen1[n] and Pp[n]:
    for (n = 0; n <= nmaxmag; n++)
        Bcen1[n] = 0.;
    Pp[0] = 0.;
    Pp[1] = 1.;
// Bcen1[n] calculation
    del = (Rmax - Rmin) / M;
    sigmac = sigma;
// Integration loop:
    for (m = 0; m <= M; m++) {
        R = Rmin + del * m;
        for (iz = 1; iz <= 2; iz++) {
            if (iz == 1)
                Z = Zmax;
            else
                Z = Zmin;
            z = Z - z0;
            ro = sqrt(R * R + z * z);
            u = z / ro;
            for (n = 2; n <= nmaxmag + 1; n++)
                Pp[n] =
                    ((2 * n - 1) * u * Pp[n - 1] -
                     n * Pp[n - 2]) / (1. * (n - 1.));
            constcen = mu0 * sigmac / 2. * (1. - u * u) / rocen;
            rcen = rocen / ro;
            rcen1 = rcen;
            for (n = 0; n <= nmaxmag; n++) {
                b[n][iz] = constcen * rcen1 * Pp[n + 1];
                rcen1 *= rcen;
            }
        }
        Z = Zmax;
        z = Z - z0;
        ro = sqrt(R * R + z * z);
        bhat[0][m] = mu0 * sigmac / 2. * z / ro;
        Z = Zmin;
        z = Z - z0;
        ro = sqrt(R * R + z * z);
        bhat[0][m] -= mu0 * sigmac / 2. * z / ro;
        for (n = 1; n <= nmaxmag; n++)
            bhat[n][m] = -rocen / n * (b[n - 1][1] - b[n - 1][2]);
    }
// end of m-loop
    for (n = 0; n <= nmaxmag; n++) {
        q = 0.;
        for (m = 0; m <= M; m++)
            q += w[m] * bhat[n][m];
        q *= del;
        Bcen1[n] += q;
    }
}

// Magnetic field calculation by central Legendre-series expansion
//  (axisymmetric case)

void magfield2_central(int type, int s, double z0, double rocen,
                       double z, double r, double *Bz, double *Br, double *rcp)
{
// This function computes the magnetic field components Bz and Br,
// either of one axially symmetric coil (type=1), or of all the
// axisymmetric coils (type=0), 
// in a fieldpoint with (z,r)
// cylindrical coordinates, using the central Legendre-series expansion.
// It computes also the central convergence ratio  *rcp=rc=ro/rocen. 
// z0: source point.
// rocen: central convergence radius.
// type=0:
//     Bcenaxisymm[s][n], n=0,...,nmaxmag:  central source
//     coefficients for the axisymmetric coils.
//     s is the index of the central source points for axisymmetric coils. 
// type=1:
//     Bcen[s][n], s=1,...,Ncen; n=0,...,nmaxmag: central source
//     coefficients in the local system of the coils. 
//     The source index s is defined as s=smin[i]+j-1
//     (for coil i, source point j), where smin[i] is
//      the minimal s value corresponding to the coil i.
//  nmaxmag: maximal index of the source constants (maximum of n).
//  SI units are used !
// If the convergence of the Legendre-series expansion is too slow:
//     *rcp>0.999,   *Bz=1.e21,  *Br=1.e21
//  (in this case one should use another computation method !!!)
    /*static double c1[nmaxmag + 1], c2[nmaxmag + 1], c3[nmaxmag + 1],
        c4[nmaxmag + 1];
    static double c5[nmaxmag + 1], c6[nmaxmag + 1];*/
    static double *c1 = NULL, *c2 = NULL, *c3 = NULL, *c4 = NULL, *c5 = NULL;
    static double *c6 = NULL;
    static int iff = 0;
    int n;
    double ro, u, delz, sr;
    double rc, rcn;
    double Beps, Bdelta;
    double *P = NULL, *Pp = NULL, *Bzplus = NULL, *Brplus = NULL;
    
    P = (double*) malloc(sizeof(double) * (nmaxmag + 1));
    Pp = (double*) malloc(sizeof(double) * (nmaxmag + 1));
    Bzplus = (double*) malloc(sizeof(double) * (nmaxmag + 1));
    Brplus = (double*) malloc(sizeof(double) * (nmaxmag + 1));
//
// Initialization of c1,c2,c3,c4,c5,c6 vectors:
    if (iff == 0) {
        c1 = (double*) malloc(sizeof(double) * (nmaxmag + 1));
        c2 = (double*) malloc(sizeof(double) * (nmaxmag + 1));
        c3 = (double*) malloc(sizeof(double) * (nmaxmag + 1));
        c4 = (double*) malloc(sizeof(double) * (nmaxmag + 1));
        c5 = (double*) malloc(sizeof(double) * (nmaxmag + 1));
        c6 = (double*) malloc(sizeof(double) * (nmaxmag + 1));
        for (n = 2; n <= nmaxmag; n++) {
            c1[n] = (2. * n - 1.) / (1. * n);
            c2[n] = (n - 1.) / (1. * n);
            c3[n] = (2. * n - 1.) / (1. * (n - 1.));
            c4[n] = (1. * n) / (1. * (n - 1.));
            c5[n] = (1.) / (n * 1.);
            c6[n] = (1.) / (n + 1.);
        }
        iff = 1;
    }
// If the field point is very close to the source point:
    if (r < 1.e-12 && fabs(z - z0) < 1.e-12) {
        if (type == 1)
            *Bz = Bcen[s][0];
        else
            *Bz = Bcenaxisymm[s][0];
        *Br = 0.;
        *rcp = 0.;
        free(P);
        free(Pp);
        free(Bzplus);
        free(Brplus);
        return;
    }
// ro,u,s,rc,rcn:
    delz = z - z0;
    ro = sqrt(r * r + delz * delz);
    u = delz / ro;
    sr = r / ro;
    rc = ro / rocen;            // convergence ratio
    *rcp = rc;
    rcn = rc;
// If rc>=0.999: the Legendre polynomial series is too slow,
//     or not convergent !!!
    if (rc >= 0.999) {
        *Bz = 1.e21;
        *Br = 1.e21;
        free(P);
        free(Pp);
        free(Bzplus);
        free(Brplus);
        return;
    }
// First 2 terms of Legendre polynomial P and its derivative Pp (P-primed)
    P[0] = 1.;
    P[1] = u;
    Pp[0] = 0.;
    Pp[1] = 1.;
// First 2 terms of the series:
    if (type == 1) {
        *Bz = Bcen[s][0] + Bcen[s][1] * rc * u;
        *Br = -sr * Bcen[s][1] / 2. * rc;
    } else {
        *Bz = Bcenaxisymm[s][0] + Bcenaxisymm[s][1] * rc * u;
        *Br = -sr * Bcenaxisymm[s][1] / 2. * rc;
    }
//
// We start here the series expansion:
    for (n = 2; n <= nmaxmag - 1; n++) {
        rcn *= rc;
        P[n] = c1[n] * u * P[n - 1] - c2[n] * P[n - 2];
        Pp[n] = c3[n] * u * Pp[n - 1] - c4[n] * Pp[n - 2];
        if (type == 1) {
            Bzplus[n] = Bcen[s][n] * rcn * P[n];
            Brplus[n] = -sr * Bcen[s][n] * c6[n] * rcn * Pp[n];
        } else {
            Bzplus[n] = Bcenaxisymm[s][n] * rcn * P[n];
            Brplus[n] = -sr * Bcenaxisymm[s][n] * c6[n] * rcn * Pp[n];
        }
        *Bz += Bzplus[n];
        *Br += Brplus[n];
        if (n > 5) {
            Beps = 1.e-15 * (fabs(*Bz) + fabs(*Br));
            Bdelta = fabs(Bzplus[n]) + fabs(Brplus[n]) +
                fabs(Bzplus[n - 1]) + fabs(Brplus[n - 1]) +
                fabs(Bzplus[n - 2]) + fabs(Brplus[n - 2]) +
                fabs(Bzplus[n - 3]) + fabs(Brplus[n - 3]);
            if (Bdelta < Beps || Bdelta < 1.e-20)
                break;
        }
    }
    if (n >= nmaxmag - 1) {
        *rcp = 1.;
        *Bz = 1.e21;
        *Br = 1.e21;
    }
    free(P);
    free(Pp);
    free(Bzplus);
    free(Brplus);
}

// Remote convergence radius calculation:

double funrorem(double *coilpar, double z0)
// This function computes the remote convergence radius 
// funrorem=rorem for an axisymmetric coil, at the axis source point z0.
// rorem = maximum distance of the axis point z0 from the coil winding.
// The coil parameters are given by *coilpar:
// coilpar[0]: current density
// coilpar[1]: Zmin, coilpar[2]: Zmax,
// coilpar[3]: Rmin, coilpar[4]: Rmax.
{
    double ro, Zmin, Zmax, Rmax, Zbar;

    Zmin = coilpar[1];
    Zmax = coilpar[2];
    Rmax = coilpar[4];
    Zbar = (Zmin + Zmax) / 2.;
    if (z0 >= Zbar)
        ro = sqrt((z0 - Zmin) * (z0 - Zmin) + Rmax * Rmax);
    else
        ro = sqrt((z0 - Zmax) * (z0 - Zmax) + Rmax * Rmax);
    return ro;
}

// Remote source coefficient calculation for 1 coil:
//  (axisymm. case)

void magsource2_remote_1coil(double *coilpar, double z0, double rorem,
                             double *Brem1)
// This subroutine computes the magnetic field remote source constants Brem1[n]
//  (n=0,...,nmaxmag) of an axisymmetric coil, in the axis source point z0,
//   with remote convergence radius rorem as input parameter.
// The coil parameters are given by *coilpar:
// coilpar[0]: current density
// coilpar[1]: Zmin, coilpar[2]: Zmax,
// coilpar[3]: Rmin, coilpar[4]: Rmax.
// nmaxmag is defined by #define command.
// Indices: m: numerical integration;
//          n: source constants
//   iz=1: Z=Zmax; iz=2: Z=Zmin
// Radial integration number: M; this depends on the radial thickness:
//  M is small if the(Rmax-Rmin)/Rmin ratio is small, and it is large
//   if this ratio is large.
{
    double bs[nmaxmag + 2][3], bshat[nmaxmag + 1][1001], Pp[nmaxmag + 2];
    int m, n, iz, M;
    const double w9[10] = { 0.2803440531305107e0, 0.1648702325837748e1,
        -0.2027449845679092e0, 0.2797927414021179e1,
        -0.9761199294532843e0, 0.2556499393738999e1,
        0.1451083002645404e0, 0.1311227127425048e1,
        0.9324249063051143e0, 0.1006631393298060e1
    };
    const double mu0 = 4. * M_PI * 1.e-7;
    static double w[1001];
    double ro, Z, R, sigma, del, q, st;
    double Zmin, Zmax, Rmin, Rmax, z, u, constrem, sigmac, rrem, rrem1, ratio;

// Coil parameters:
    sigma = coilpar[0];
    Zmin = coilpar[1];
    Zmax = coilpar[2];
    Rmin = coilpar[3];
    Rmax = coilpar[4];
// Improvement of Zmin,Zmax,Rmin,Rmax values:
    if (Zmax < Zmin) {
        st = Zmax;
        Zmax = Zmin;
        Zmin = st;
    }
    if (Rmax < Rmin) {
        st = Rmax;
        Rmax = Rmin;
        Rmin = st;
    }
    if (Zmax - Zmin < 1.e-12)
        Zmax = Zmin + 1.e-12;
    if (Rmax - Rmin < 1.e-12)
        Rmax = Rmin + 1.e-12;
// Coil geometry test:
    if (Rmin <= 1.e-9) {
        fprintf(stderr, "Message from function magsource2_central_1coil:\n"
                "Rmin is not positive!!!\n" "Computation is stopped!!!\n");
        exit(1);
    }
// Radial integration number M:
    ratio = (Rmax - Rmin) / Rmin;
    if (ratio < 0.1)
        M = 30;
    else if (ratio >= 0.1 && ratio < 0.2)
        M = 50;
    else
        M = 60 * ratio / 0.2;
    if (M > 1000)
        M = 1000;
// Initialization of the integration weight factors:
    for (m = 0; m <= 9; m++)
        w[m] = w9[m];
    for (m = 10; m <= M - 10; m++)
        w[m] = 1.;
    for (m = M - 9; m <= M; m++)
        w[m] = w9[M - m];
// Initialization of Brem1[n] and Pp[n]:
    for (n = 0; n <= nmaxmag; n++)
        Brem1[n] = 0.;
    Pp[0] = 0.;
    Pp[1] = 1.;
// Brem1[n] calculation
    del = (Rmax - Rmin) / M;
    sigmac = sigma;
// Integration loop:
    for (m = 0; m <= M; m++) {
        R = Rmin + del * m;
        for (iz = 1; iz <= 2; iz++) {
            if (iz == 1)
                Z = Zmax;
            else
                Z = Zmin;
            z = Z - z0;
            ro = sqrt(R * R + z * z);
            u = z / ro;
            for (n = 2; n <= nmaxmag + 1; n++)
                Pp[n] =
                    ((2 * n - 1) * u * Pp[n - 1] -
                     n * Pp[n - 2]) / (1. * (n - 1.));
            constrem = mu0 * sigmac / 2. * (1. - u * u) / rorem;
            rrem = ro / rorem;
            rrem1 = rrem * rrem;
            for (n = 2; n <= nmaxmag + 1; n++) {
                bs[n][iz] = constrem * rrem1 * Pp[n - 1];
                rrem1 *= rrem;
            }
        }
        bshat[0][m] = bshat[1][m] = 0.;
        for (n = 2; n <= nmaxmag; n++)
            bshat[n][m] = rorem / (n + 1.) * (bs[n + 1][1] - bs[n + 1][2]);
    }
// end of m-loop
    for (n = 0; n <= nmaxmag; n++) {
        q = 0.;
        for (m = 0; m <= M; m++)
            q += w[m] * bshat[n][m];
        q *= del;
        Brem1[n] += q;
    }
}

// Magnetic field calculation by remote Legendre-series expansion
//  (axisymmetric case)

void magfield2_remote(int type, int s, double z0, double rorem,
                      double z, double r, double *Bz, double *Br, double *rcp)
{
// This function computes the magnetic field components Bz and Br
// either of one axially symmetric coil (type=1 or 2), or of all the
// axisymmetric coils (type=0), 
// in a fieldpoint with (z,r)
// cylindrical coordinates, using the remote Legendre-series expansion.
// It computes also the remote convergence ratio  *rcp=rc=rorem/ro. 
// z0: source point
// rorem: remote convergence radius.
// type=0:
//     Bremaxisymm[n], n=0,...,nmaxmag:  remote source
//     coefficients for the axisymmetric coils.
// type=1 or 2:
//    Brem[s][n], s=1,...,3*Ncoil; n=0,...,nmaxmag: remote source
//    coefficients. For each coil we have 3 source points (j=1,2,3),
//    the source index s is defined as s=3*(i-1)+j
//    (for coil i, source point j).
//  type=1:  source point in coil centre (j=1), 
//           only even terms are needed in the series
//  type=2:  source point not in coil centre (j=2 or 3),
//          both even and odd terms are needed in the series
//  nmaxmag: maximal index of the source constants (maximum of n).
//  SI units are used !
// If the convergence of the Legendre-series expansion is too slow:
//     *rcp>0.999,   *Bz=1.e21,  *Br=1.e21
//  (in this case one should use another computation method !!!)
    /*static double c1[nmaxmag + 1], c2[nmaxmag + 1], c3[nmaxmag + 1],
        c4[nmaxmag + 1];
    static double c5[nmaxmag + 1], c6[nmaxmag + 1];
    static double c7[nmaxmag + 1], c8[nmaxmag + 1], c9[nmaxmag + 1];
    static double c10[nmaxmag + 1], c11[nmaxmag + 1], c12[nmaxmag + 1];*/
    static double *c1 = NULL, *c2 = NULL, *c3 = NULL, *c4 = NULL, *c5 = NULL;
    static double *c6 = NULL, *c7 = NULL, *c8 = NULL, *c9 = NULL, *c10 = NULL;
    static double *c11 = NULL, *c12 = NULL;
    static int iff = 0;
    int n, m;
    double ro, u, delz, sr;
    double P[nmaxmag + 1], Pp[nmaxmag + 1], rc, rcn;
    double Bzplus[nmaxmag + 1], Brplus[nmaxmag + 1];
    double Beps, Bdelta;
    double u2, rc2, M, Mp, A, Ap, B, Bp, C, Cp;

//
// Initialization of c1,c2,c3,c4,c5,c6 vectors:
    if (iff == 0) {
        c1 = (double*) malloc(sizeof(double) * (nmaxmag + 1));
        c2 = (double*) malloc(sizeof(double) * (nmaxmag + 1));
        c3 = (double*) malloc(sizeof(double) * (nmaxmag + 1));
        c4 = (double*) malloc(sizeof(double) * (nmaxmag + 1));
        c5 = (double*) malloc(sizeof(double) * (nmaxmag + 1));
        c6 = (double*) malloc(sizeof(double) * (nmaxmag + 1));
        c7 = (double*) malloc(sizeof(double) * (nmaxmag + 1));
        c8 = (double*) malloc(sizeof(double) * (nmaxmag + 1));
        c9 = (double*) malloc(sizeof(double) * (nmaxmag + 1));
        c10 = (double*) malloc(sizeof(double) * (nmaxmag + 1));
        c11 = (double*) malloc(sizeof(double) * (nmaxmag + 1));
        c12 = (double*) malloc(sizeof(double) * (nmaxmag + 1));
        for (n = 2; n <= nmaxmag; n++) {
            c1[n] = (2. * n - 1.) / (1. * n);
            c2[n] = (n - 1.) / (1. * n);
            c3[n] = (2. * n - 1.) / (1. * (n - 1.));
            c4[n] = (1. * n) / (1. * (n - 1.));
            c5[n] = (1.) / (n * 1.);
            c6[n] = (1.) / (n + 1.);
            if (n >= 6) {
                m = n - 2;
                M = (m + 1.) * (m + 2.) * (2. * m - 1.);
                Mp = (m) * (m + 1.) * (2. * m - 1.);
                A = (2. * m - 1.) * (2. * m + 1.) * (2. * m + 3.);
                Ap = A;
                B = 2. * m * m * (2. * m + 3.) - 1.;
                Bp = 2. * m * (m + 2.) * (2. * m - 1.) - 3.;
                C = m * (m - 1.) * (2. * m + 3.);
                Cp = m * (m + 1.) * (2. * m + 3.);
                c7[n] = A / M;
                c8[n] = B / M;
                c9[n] = C / M;
                c10[n] = Ap / Mp;
                c11[n] = Bp / Mp;
                c12[n] = Cp / Mp;
            }
        }
        iff = 1;
    }
// ro,u,s,rc,rcn:
    delz = z - z0;
    ro = sqrt(r * r + delz * delz);
    if (ro < 1.e-9)
        ro = 1.e-9;
    u = delz / ro;
    sr = r / ro;
    rc = rorem / ro;            // convergence ratio
    *rcp = rc;
    rcn = rc * rc;
// If rc>=0.999: the Legendre polynomial series is too slow,
//     or not convergent !!!
    if (rc >= 0.999) {
        *Bz = 1.e21;
        *Br = 1.e21;
        return;
    }
// First 2 terms of Legendre polynomial P and its derivative Pp (P-primed)
    P[0] = 1.;
    P[1] = u;
    Pp[0] = 0.;
    Pp[1] = 1.;
//
// We start here the series expansion:
    *Bz = *Br = 0.;
    if (type == 1)
        goto label1;
// Both even and odd terms are needed in the series:
    for (n = 2; n <= nmaxmag - 1; n++) {
        P[n] = c1[n] * u * P[n - 1] - c2[n] * P[n - 2];
        Pp[n] = c3[n] * u * Pp[n - 1] - c4[n] * Pp[n - 2];
        rcn = rcn * rc;
        if (type == 0) {
            Bzplus[n] = Bremaxisymm[n] * rcn * P[n];
            Brplus[n] = sr * Bremaxisymm[n] * c5[n] * rcn * Pp[n];
        } else {
            Bzplus[n] = Brem[s][n] * rcn * P[n];
            Brplus[n] = sr * Brem[s][n] * c5[n] * rcn * Pp[n];
        }
        *Bz += Bzplus[n];
        *Br += Brplus[n];
        if (n > 5) {
            Beps = 1.e-15 * (fabs(*Bz) + fabs(*Br));
            Bdelta = fabs(Bzplus[n]) + fabs(Brplus[n]) +
                fabs(Bzplus[n - 1]) + fabs(Brplus[n - 1]) +
                fabs(Bzplus[n - 2]) + fabs(Brplus[n - 2]) +
                fabs(Bzplus[n - 3]) + fabs(Brplus[n - 3]);
            if (Bdelta < Beps || Bdelta < 1.e-20)
                break;
        }
    }
    goto label2;
  label1:;
// Only the even terms are needed in the series
//   (because for the middle source point in the coil centre
//     Brem[s][n]=0 for odd n)
    for (n = 2; n <= 4; n++) {
        P[n] = c1[n] * u * P[n - 1] - c2[n] * P[n - 2];
        Pp[n] = c3[n] * u * Pp[n - 1] - c4[n] * Pp[n - 2];
        rcn = rcn * rc;
        Bzplus[n] = Brem[s][n] * rcn * P[n];
        Brplus[n] = sr * Brem[s][n] * c5[n] * rcn * Pp[n];
        *Bz += Bzplus[n];
        *Br += Brplus[n];
    }
    u2 = u * u;
    rc2 = rc * rc;
    for (n = 6; n <= nmaxmag - 1; n += 2) {
        P[n] = (c7[n] * u2 - c8[n]) * P[n - 2] - c9[n] * P[n - 4];
        Pp[n] = (c10[n] * u2 - c11[n]) * Pp[n - 2] - c12[n] * Pp[n - 4];
        rcn = rcn * rc2;
        Bzplus[n] = Brem[s][n] * rcn * P[n];
        Brplus[n] = sr * Brem[s][n] * c5[n] * rcn * Pp[n];
        *Bz += Bzplus[n];
        *Br += Brplus[n];
        Beps = 1.e-15 * (fabs(*Bz) + fabs(*Br));
        Bdelta = fabs(Bzplus[n]) + fabs(Brplus[n]) +
            fabs(Bzplus[n - 2]) + fabs(Brplus[n - 2]);
        if (Bdelta < Beps || Bdelta < 1.e-20)
            break;
    }
  label2:;
    if (n >= nmaxmag - 2) {
        *rcp = 1.;
        *Bz = 1.e21;
        *Br = 1.e21;
    }
}

// Part 3: 3-dimensional coils 

void magfield3_elliptic_1coil(int i, double *P, double *B)
{
// This function computes the magnetic field components B[1],B[2],B[3]
// in a field point P[1],P[2],P[3], due to the 3-dimensional coil
// with index i, using the first, second and third complete elliptic
// integrals.  SI units are used (P[k] in m, B[k] in T, k=1,2,3)!
// The coil is defined by the coil parameters coil[i][j}, j=0,...,13.
    int n, k;
    double coilpar[5], z, r, Bz, Br, u[4], length, Ploc[4], Pr[4], w[4];

// We define local coordinate system: 
//    origo in endpoint A of the coil axis,
//    local z axis parallel to coil axis, from A to B  
// Unit vector u in local z axis direction:
    length = coil[i][10];
    for (k = 1; k <= 3; k++)
        u[k] = coil[i][10 + k];
// Coil parameters:
    coilpar[0] = coil[i][0];
    coilpar[1] = 0.;
    coilpar[2] = length;
    coilpar[3] = coil[i][7];
    coilpar[4] = coil[i][8];
    n = coil[i][9] + 0.1;
// Local z and r coordinates of the field point P:
    for (k = 1; k <= 3; k++)
        Ploc[k] = P[k] - coil[i][k];
    z = Ploc[1] * u[1] + Ploc[2] * u[2] + Ploc[3] * u[3];
    for (k = 1; k <= 3; k++)
        Pr[k] = Ploc[k] - z * u[k];
    r = sqrt(Pr[1] * Pr[1] + Pr[2] * Pr[2] + Pr[3] * Pr[3]);
// Bz and Br calculation:     
    magfield2_elliptic_1coil(n, coilpar, z, r, &Bz, &Br);
// B[k] calculation from Bz, Br: 
    if (r < 1.e-13)
        for (k = 1; k <= 3; k++)
            B[k] = Bz * u[k];
    else {
        for (k = 1; k <= 3; k++)
            w[k] = Pr[k] / r;
        for (k = 1; k <= 3; k++)
            B[k] = Bz * u[k] + Br * w[k];
    }
}

// Central source coefficient calculation for all coils:
//  (3-dim. case)

void magsource3_central()
// This function computes the Bcen[s][n] (s=1,...,Ncen,n=0,...,nmaxmag+1)  //   central source coefficients for all coils.
// The results are written into the data file 'magsource_central.dat':
//     Ncoil: number of coils
//     nmaxmag: number of source coefficients for a fixed coil and
//              source point is nmaxmag+1 (n=0,...,nmaxmag)
//     i: coil index
//     Nsp: number of source points for a fixed coil
//     smin: for a fixed coil (with index i) the source index s of
//            Bcen[s][n] starts at s=smin
//     j: source point index for a fixed coil with index i
//     z0: local z value of the source point, for a fixed coil with index i
//     rocen: central convergence radius for a fixed source point and coil
//     Bcen[s][n], n=0,...,nmaxmag: central source coefficients
//              (for a fixed source point and coil).
{
    int i, n, k, Nsp, j, smin, nn, Ncen;
    double z0, coilpar[5], u[4];
    double Bcen1[nmaxmag + 1], rocen;
    double Zcen, Rmin, Rmax, del, length, rorem, d;
    FILE *fp;

// Output to file magsource_central.dat:
    /* fp=fopen("magsource_central.dat","w"); */
    fp = fopen(magfield3_get_magsource_central(), "w");
// Coil index loop: 
    smin = 1;
    for (i = 1; i <= Ncoil; i++) {
// We define local coordinate system: 
//    origo in endpoint A of the coil axis,
//    local z axis parallel to coil axis, from A to B.  
//    Unit vector u in local z axis direction:
        length = coil[i][10];
        for (k = 1; k <= 3; k++)
            u[k] = coil[i][10 + k];
// coilpar calc.:     
        coilpar[0] = coil[i][0];
        coilpar[1] = 0.;
        coilpar[2] = length;
        coilpar[3] = coil[i][7];
        coilpar[4] = coil[i][8];
// Number of source points (Nsp): 
        Zcen = length / 2.;     // local z value of coil centre
        Rmin = coil[i][7];      // inner radius of coil
        Rmax = coil[i][8];      // outer radius of coil
        rorem = sqrt(Rmax * Rmax + Zcen * Zcen);
        del = Rmin * 0.3;
        d = 1.3 * rorem;
        nn = d / del + 1;
        d = nn * del;
        Nsp = 2 * nn + 1;
        fprintf(fp, "%9i %9i %9i  \n", i, Nsp, smin);
        smin += Nsp;
// Source point loop:      
        for (j = 1; j <= Nsp; j++) {
            z0 = Zcen - d + del * (j - 1.);
// B[n] calc.:  
            rocen = funrocen(coilpar, z0);
            magsource2_central_1coil(coilpar, z0, rocen, Bcen1);
// Output to file magsource_central.dat (j, z0, rocen, Bcen1[n])
            fprintf(fp, "%9i %22.15e %22.15e \n", j, z0, rocen);
            for (n = 0; n <= nmaxmag; n++) {
                if (fabs(Bcen1[n]) < 1.e-30)
                    Bcen1[n] = 0.;
                fprintf(fp, "%9i %22.15e  \n", n, Bcen1[n]);
            }
        }
    }
    fclose(fp);
    Ncen = smin - 1;
}

// Remote source coefficient calculation for all coils:
//  (3-dim. case)

void magsource3_remote()
// This function computes the Brem[s][n] (s=1,...,3*Ncoil,n=0,...,nmaxmag+1)  //   remote source coefficients for all coils.
// Number of source points for a fixed coil is 3.
// The results are written into the data file 'magsource_remote.dat':
//     Ncoil: number of coils
//     nmaxmag: number of source coefficients for a fixed coil and
//              source point is nmaxmag+1 (n=0,...,nmaxmag)
//     i: coil index
//     j: source point index for a fixed coil with index i (j=1,2,3)
//     z0: local z value of the source point, for a fixed coil with index i
//     rorem: remote convergence radius for a fixed source point and coil
//     Brem1[n], n=0,...,nmaxmag: remote source coefficients
//              (for a fixed source point and coil).
{
    int i, n, k, Nsp, j;
    double z0, coilpar[5], u[4];
    double Brem1[nmaxmag + 1];
    double Zcen, Rmax, length, rorem, d;
    FILE *fp;

// Output to file magsource_remote.dat:
    /* fp=fopen("magsource_remote.dat","w"); */
    fp = fopen(magfield3_get_magsource_remote(), "w");
// Coil index loop: 
    for (i = 1; i <= Ncoil; i++) {
// We define local coordinate system: 
//    origo in endpoint A of the coil axis,
//    local z axis parallel to coil axis, from A to B.  
//    Unit vector u in local z axis direction:
        length = coil[i][10];
        for (k = 1; k <= 3; k++)
            u[k] = coil[i][10 + k];
// coilpar calc.:     
        coilpar[0] = coil[i][0];
        coilpar[1] = 0.;
        coilpar[2] = length;
        coilpar[3] = coil[i][7];
        coilpar[4] = coil[i][8];
// Source points: 
        Zcen = length / 2.;     // local z value of coil centre
        Rmax = coil[i][8];      // outer radius of coil
        rorem = sqrt(Rmax * Rmax + Zcen * Zcen);
        d = 2. * rorem;
        Nsp = 3;
// Source point loop:      
        for (j = 1; j <= Nsp; j++) {
            if (j == 1)
                z0 = Zcen;
            else if (j == 2)
                z0 = Zcen - d;
            else if (j == 3)
                z0 = Zcen + d;
// Brem1[n] calc.:  
            rorem = funrorem(coilpar, z0);
            magsource2_remote_1coil(coilpar, z0, rorem, Brem1);
// Output to file magsource_remote.dat (j, z0, rorem, Brem1[n])
            fprintf(fp, "%9i %22.15e %22.15e \n", j, z0, rorem);
            for (n = 0; n <= nmaxmag; n++) {
                if (fabs(Brem1[n]) < 1.e-30)
                    Brem1[n] = 0.;
                fprintf(fp, "%9i %22.15e  \n", n, Brem1[n]);
            }
        }
    }
    fclose(fp);
}

void magfield3_1coil(int i, double *P, double *B)
{
// This function computes the magnetic field components B[1],B[2],B[3]
// in a field point P[1],P[2],P[3], due to the 3-dimensional coil
// with index i, using central or remote Legendre polynomial expansion,
// or elliptic integrals.
// SI units are used (P[k] in m, B[k] in T, k=1,2,3)!
// The coil is defined by the coil parameters coil[i][j}, j=0,...,13.
    FILE *fp;
    static int iff = 0;
    int ic, ix, j, jx, s, n, nx, k;
    static int *Nsp = NULL, *smin = NULL;
    static double *z0cen = NULL, *rocen = NULL;
    static double *z0rem = NULL, *rorem = NULL;
    double coilpar[5], z, r, Bz, Br, u[4], length, Ploc[4], Pr[4], w[4], delz;
    double rcp, ro2, rcmin2, rocen2, rorem2;
    int num, jbest = 0, jmin, jmax;
    static int *jlast = NULL;

// First call of this function: input of  source parameters 
//   from files magsource_central.dat and magsource_remote.dat.
    if (iff == 0) {
        Nsp = (int*) malloc(sizeof(int) * (Ncoil + 1));
        smin = (int*) malloc(sizeof(int) * (Ncoil + 1));
        z0rem = (double*) malloc(sizeof(double) * (3 * Ncoil + 1));
        rorem = (double*) malloc(sizeof(double) * (3 * Ncoil + 1));
        jlast = (int*) malloc(sizeof(int) * (Ncoil + 1));

//  Input from file magsource_central.dat:
        /* fp=fopen("magsource_central.dat","r"); */
        fp = fopen(magfield3_get_magsource_central(), "r");
        if (!fp) {
            fprintf(stderr, "Message from function magfield3_1coil:\n");
            fprintf(stderr, "Cannot open the file '%s'!\n",
                    magfield3_get_magsource_central());
            fprintf(stderr, "Program running is stopped !!!\n");
            exit(0);
        }
//  Central source parameters:
        for (ic = 1; ic <= Ncoil; ic++) {
            fscanf(fp, "%i %i %i  \n", &ix, &Nsp[ic], &smin[ic]);

            Bcen = (double**) realloc(Bcen, sizeof(double *) * (smin[ic] + Nsp[ic]));
            z0cen = (double*) realloc(z0cen, sizeof(double) * (smin[ic] + Nsp[ic]));
            rocen = (double*) realloc(rocen, sizeof(double) * (smin[ic] + Nsp[ic]));

            for (j = 1; j <= Nsp[ic]; j++) {
                s = smin[ic] + j - 1;
                Bcen[s] = (double*) malloc(sizeof(double) * (nmaxmag + 1));
                fscanf(fp, "%i %le %le \n", &jx, &z0cen[s], &rocen[s]);
                for (n = 0; n <= nmaxmag; n++)
                    fscanf(fp, "%i %le  \n", &nx, &Bcen[s][n]);
            }
        }
        fclose(fp);
//  Input from file magsource_remote.dat:
        /*fp=fopen("magsource_remote.dat","r"); */
        fp = fopen(magfield3_get_magsource_remote(), "r");
        if (!fp) {
            fprintf(stderr, "Message from function magfield3_1coil:\n");
            fprintf(stderr, "Cannot open the file '%s'!\n",
                    magfield3_get_magsource_remote());
            fprintf(stderr, "Program running is stopped!!! \n");
            exit(1);
        }
//  Remote source parameters:
        for (ic = 1; ic <= Ncoil; ic++) {
            for (j = 1; j <= 3; j++) {
                s = 3 * (ic - 1) + j;
                fscanf(fp, "%i %le %le \n", &jx, &z0rem[s], &rorem[s]);
                for (n = 0; n <= nmaxmag; n++)
                    fscanf(fp, "%i %le  \n", &nx, &Brem[s][n]);
            }
        }
        fclose(fp);
// Initialization of jlast[ic] values:
        for (ic = 1; ic <= Ncoil; ic++)
            jlast[ic] = 0;
        iff = 1;
    }
// End of input from data files magsource_central.dat and magsource_remote.dat.
//
// We start now the magnetic field calculation.
// First we define the local coordinate system of the coil, and the
//   z and r local coordinates of the field point P. 
// Local coordinate system: 
//    origo in endpoint A of the coil axis,
//    local z axis parallel to coil axis, from A to B  
// Unit vector u in local z axis direction:
    length = coil[i][10];
    for (k = 1; k <= 3; k++)
        u[k] = coil[i][10 + k];
// Local z and r coordinates of the field point P:
    for (k = 1; k <= 3; k++)
        Ploc[k] = P[k] - coil[i][k];
    z = Ploc[1] * u[1] + Ploc[2] * u[2] + Ploc[3] * u[3];
    for (k = 1; k <= 3; k++)
        Pr[k] = Ploc[k] - z * u[k];
    r = sqrt(Pr[1] * Pr[1] + Pr[2] * Pr[2] + Pr[3] * Pr[3]);
// Bz and Br calculation.
//  ---------------------------
// 1. Probably most of the coils are far from the field point, therefore
// we try first the remote series expansion
// (with source point in the middle of the coil) 
    s = 3 * (i - 1) + 1;
    delz = z - z0rem[s];
    ro2 = r * r + delz * delz;
    rorem2 = rorem[s] * rorem[s];
    if (ro2 > rorem2 * 1.4) {
        magfield2_remote(1, s, z0rem[s], rorem[s], z, r, &Bz, &Br, &rcp);
        if (rcp <= 0.999)
            goto labelend;
    }
// -----------------------------
// 2. The field point is close to the coil. We try next the central
// Legendre polynomial expansion. If some source point for this coil
// has been already used (jlast[i]>0), we search the central source
// point j close to the old source point jlast[i].
    if (jlast[i] == 0)
        goto label3;
// First we search for the best source point among the
// 5 closest points (including the source point jlast[i] itself):
    rcmin2 = 0.91;
    jmin = jlast[i] - 2;
    if (jmin < 1)
        jmin = 1;
    jmax = jlast[i] + 2;
    if (jmax > Nsp[i])
        jmax = Nsp[i];
    for (j = jmin; j <= jmax; j++) {
        s = smin[i] + j - 1;
        delz = z - z0cen[s];
        ro2 = r * r + delz * delz;
        rocen2 = rocen[s] * rocen[s];
        if (ro2 < rocen2 * rcmin2) {
            rcmin2 = ro2 / rocen2;
            jbest = j;
        }
    }
    if (rcmin2 < 0.9) {
        s = smin[i] + jbest - 1;
        magfield2_central(1, s, z0cen[s], rocen[s], z, r, &Bz, &Br, &rcp);
        if (rcp <= 0.999) {
            jlast[i] = jbest;
            goto labelend;
        }
    }
// Next, we search for the best source point among the
// 21 closest points (including the source point jlast[i] itself):
    rcmin2 = 0.981;
    jmin = jlast[i] - 10;
    if (jmin < 1)
        jmin = 1;
    jmax = jlast[i] + 10;
    if (jmax > Nsp[i])
        jmax = Nsp[i];
    for (j = jmin; j <= jmax; j++) {
        s = smin[i] + j - 1;
        delz = z - z0cen[s];
        ro2 = r * r + delz * delz;
        rocen2 = rocen[s] * rocen[s];
        if (ro2 < rocen2 * rcmin2) {
            rcmin2 = ro2 / rocen2;
            jbest = j;
        }
    }
    if (rcmin2 < 0.98) {
        s = smin[i] + jbest - 1;
        magfield2_central(1, s, z0cen[s], rocen[s], z, r, &Bz, &Br, &rcp);
        if (rcp <= 0.999) {
            jlast[i] = jbest;
            goto labelend;
        }
    }
// -----------------------------
  label3:;
// 3. We search now among all central source points.
    rcmin2 = 0.9981;
    for (j = 1; j <= Nsp[i]; j++) {
        s = smin[i] + j - 1;
        delz = z - z0cen[s];
        ro2 = r * r + delz * delz;
        rocen2 = rocen[s] * rocen[s];
        if (ro2 < rocen2 * rcmin2) {
            rcmin2 = ro2 / rocen2;
            jbest = j;
        }
    }
    if (rcmin2 < 0.998) {
        s = smin[i] + jbest - 1;
        magfield2_central(1, s, z0cen[s], rocen[s], z, r, &Bz, &Br, &rcp);
        if (rcp <= 0.999) {
            jlast[i] = jbest;
            goto labelend;
        }
    }
// -----------------------------
// 4. For this coil no central source point was found.  We try now the
//    3 remote source points.
    rcmin2 = 0.9981;
    for (j = 1; j <= 3; j++) {
        s = 3 * (i - 1) + j;
        delz = z - z0rem[s];
        ro2 = r * r + delz * delz;
        rorem2 = rorem[s] * rorem[s];
        if (ro2 * rcmin2 > rorem2) {
            rcmin2 = rorem2 / ro2;
            jbest = j;
        }
    }
    if (rcmin2 < 0.998) {
        j = jbest;
        s = 3 * (i - 1) + j;
        if (j == 1)
            magfield2_remote(1, s, z0rem[s], rorem[s], z, r, &Bz, &Br, &rcp);
        else
            magfield2_remote(2, s, z0rem[s], rorem[s], z, r, &Bz, &Br, &rcp);
        if (rcp <= 0.999)
            goto labelend;
    }
// -----------------------------
// 5. Unfortunately, no appropriate central or remote source point was
// found. We have to use elliptic integrals:      
    coilpar[0] = coil[i][0];
    coilpar[1] = 0.;
    coilpar[2] = length;
    coilpar[3] = coil[i][7];
    coilpar[4] = coil[i][8];
    num = coil[i][9] + 0.1;
    magfield2_elliptic_1coil(num, coilpar, z, r, &Bz, &Br);
// B[k],k=1,2,3 calculation from Bz, Br:
  labelend:;
    if (r < 1.e-16 || fabs(Br) < fabs(Bz) * 1.e-15)
        for (k = 1; k <= 3; k++)
            B[k] = Bz * u[k];
    else {
        for (k = 1; k <= 3; k++)
            w[k] = Pr[k] / r;
        for (k = 1; k <= 3; k++)
            B[k] = Bz * u[k] + Br * w[k];
    }
}

// Part 4: calculation of the axisymmetric coils     

void magsource_axisymm()
{
// This subroutine computes the magnetic field central and remote source
// coefficients for the axisymmetric coils.
// The results are written into the data file 'magsource_asisymm.dat':
    double z0rem, zmin, zmax, coilpar[5], rorem, rorem1, d;
    double *z0cen, *rocen, rocen1;
    double Bcen1[nmaxmag + 1], Brem1[nmaxmag + 1];
    int i, n, j, ia, Nsp;
    FILE *fp;
    Array *array_z0cen, *array_rocen;
    double z0cen_temp, rocen_temp;

//
// Remote source point of axisymmetric coils
//   (centre of axisymmetric coils): z0rem.
// zmin, zmax: minimal and maximal z values of the
// axisymmetric coils.
    zmin = 1.e9;
    zmax = -1.e9;
    for (ia = 1; ia <= Ncoilaxisymm; ia++) {
        i = indexaxisymm[ia][0];
        if (coil[i][3] < zmin)
            zmin = coil[i][3];
        if (coil[i][6] > zmax)
            zmax = coil[i][6];
    }
    z0rem = (zmax + zmin) / 2.;
// Remote convergence radius rorem of the axisymmetric coils, 
// with z0rem as remote source point: 
    rorem = 0.;
    for (ia = 1; ia <= Ncoilaxisymm; ia++) {
        i = indexaxisymm[ia][0];
        coilpar[0] = coil[i][0];
        coilpar[1] = coil[i][3];
        coilpar[2] = coil[i][6];
        coilpar[3] = coil[i][7];
        coilpar[4] = coil[i][8];
        rorem1 = funrorem(coilpar, z0rem);
        if (rorem1 > rorem)
            rorem = rorem1;
    }
// Central source points of the axisymmetric coils:
    d = 1.3 * rorem;
    Nsp = 0;
    z0cen_temp = z0rem - d;
    rocen_temp = funrocenaxisymm(z0cen_temp);
    array_z0cen = array_new();
    array_append(array_z0cen, 0);
    array_append(array_z0cen, z0cen_temp);
    array_rocen = array_new();
    array_append(array_rocen, 0);
    array_append(array_rocen, rocen_temp);
    j = 1;

    while (!Nsp) {
        if (z0cen_temp > z0rem + d)
            Nsp = j;
        rocen1 = funrocenaxisymm(z0cen_temp + rocen_temp * 0.2);
        z0cen_temp = z0cen_temp + FMIN(rocen_temp, rocen1) * 0.2;
        rocen_temp = funrocenaxisymm(z0cen_temp);
        array_append(array_z0cen, z0cen_temp);
        array_append(array_rocen, rocen_temp);
        j++;
    }
    z0cen = array_ptr(array_z0cen);
    rocen = array_ptr(array_rocen);
    array_free(array_z0cen);
    array_free(array_rocen);

// Output to file magsource_axisymm.dat:
    /*fp=fopen("magsource_axisymm.dat","w"); */
    fp = fopen(magfield3_get_magsource_axisymm(), "w");
    fprintf(fp, "%9i \n", Nsp);

    allocate_arraysNsp(Nsp);

// Source point loop:      
    for (j = 1; j <= Nsp; j++) {
// Bcenaxisymm[j][n] calc.:  
        for (n = 0; n <= nmaxmag; n++)
            Bcenaxisymm[j][n] = 0.;
        for (ia = 1; ia <= Ncoilaxisymm; ia++) {
            i = indexaxisymm[ia][0];
            coilpar[0] = coil[i][0];
            coilpar[1] = coil[i][3];
            coilpar[2] = coil[i][6];
            coilpar[3] = coil[i][7];
            coilpar[4] = coil[i][8];
            magsource2_central_1coil(coilpar, z0cen[j], rocen[j], Bcen1);
            for (n = 0; n <= nmaxmag; n++)
                Bcenaxisymm[j][n] += Bcen1[n];
        }
// Output to file magsource_axisymm.dat (j, z0cen, rocen, Bcenaxisymm[j][n])
        fprintf(fp, "%9i %22.15e %22.15e \n", j, z0cen[j], rocen[j]);
        for (n = 0; n <= nmaxmag; n++) {
            if (fabs(Bcenaxisymm[j][n]) < 1.e-30)
                Bcenaxisymm[j][n] = 0.;
            fprintf(fp, "%9i %22.15e  \n", n, Bcenaxisymm[j][n]);
        }
    }
// Calculation of remote source coefficients for 
// the axisymmetric coils:  
    for (n = 0; n <= nmaxmag; n++)
        Bremaxisymm[n] = 0.;
    for (ia = 1; ia <= Ncoilaxisymm; ia++) {
        i = indexaxisymm[ia][0];
        coilpar[0] = coil[i][0];
        coilpar[1] = coil[i][3];
        coilpar[2] = coil[i][6];
        coilpar[3] = coil[i][7];
        coilpar[4] = coil[i][8];
        magsource2_remote_1coil(coilpar, z0rem, rorem, Brem1);
        for (n = 0; n <= nmaxmag; n++)
            Bremaxisymm[n] += Brem1[n];
    }
// Output to file magsource_axisymm.dat (z0rem, rorem, Bremaxisymm[n])
    fprintf(fp, "%22.15e %22.15e \n", z0rem, rorem);
    for (n = 0; n <= nmaxmag; n++) {
        if (fabs(Bremaxisymm[n]) < 1.e-30)
            Bremaxisymm[n] = 0.;
        fprintf(fp, "%9i %22.15e  \n", n, Bremaxisymm[n]);
    }
    fclose(fp);

    free(z0cen);
    free(rocen);
}

double funrocenaxisymm(double z0cen)
{
// Computes the central convergence radius in an axis  point z0cen, for
//  all the axisymmetric coils.
    double rocen, coilpar[5], rocen1;
    int ia, i;

    rocen = 1.e9;
    for (ia = 1; ia <= Ncoilaxisymm; ia++) {
        i = indexaxisymm[ia][0];
        coilpar[0] = coil[i][0];
        coilpar[1] = coil[i][3];
        coilpar[2] = coil[i][6];
        coilpar[3] = coil[i][7];
        coilpar[4] = coil[i][8];
        rocen1 = funrocen(coilpar, z0cen);
        if (rocen1 < rocen)
            rocen = rocen1;
    }
    return rocen;
}

void magfield_axisymm(double *P, double *B)
{
// This function computes the magnetic field components B[1],B[2],B[3]
// in a field point P[1],P[2],P[3], due to the axisymmetric coils,
// using central or remote Legendre polynomial expansion,
// or elliptic integrals.
// This function should be called only if Ncoilaxisymm>0 and 
//  less equal to Ncoil!
// SI units are used (P[k] in m, B[k] in T, k=1,2,3)!
    FILE *fp;
    static int iff = 0;
    int i, j, jx, s, n, nx, jbest = 0, jmin, jmax, k;
    static int Nsp;
    static double *z0cen = NULL, *rocen = NULL;
    static double z0rem, rorem;
    double z, r, Bz, Br, delz, ro;
    double rcp, rc, Bi[4], ro2, rcmin2, rocen2, rc2, rc2old;
    int ia;
    static int jlast = 0;

// Test of Ncoilaxisymm:
    if (Ncoilaxisymm < 1 || Ncoilaxisymm > Ncoil) {
        fprintf(stderr, "Message from function magfield_axisymm:\n"
                "Ncoilaxisymm<1 or Ncoilaxisymm>Ncoil!!!\n"
                "Computation is stopped !!! Ncoilaxisymm= %9i\n", Ncoilaxisymm);
        exit(1);
    }
// First call of this function: input of source parameters 
//   from file magsource_axisymm.dat.
    if (iff == 0) {
//  Input from file magsource_axisymm.dat:
        /*fp=fopen("magsource_axisymm.dat","r"); */
        fp = fopen(magfield3_get_magsource_axisymm(), "r");
        if (!fp) {
            fprintf(stderr, "Message from function magfield_axisymm:\n");
            fprintf(stderr,
                    "Cannot open the file '%s'!\n",
                    magfield3_get_magsource_axisymm());
            fprintf(stderr, "Program running is stopped !!!\n");
            exit(1);
        }
//  Central source parameters:
        fscanf(fp, "%i  \n", &Nsp);
        z0cen = (double*) malloc(sizeof(double) * (Nsp + 1));
        rocen = (double*) malloc(sizeof(double) * (Nsp + 1));
        free_arraysNsp();
        allocate_arraysNsp(Nsp);
        for (j = 1; j <= Nsp; j++) {
            fscanf(fp, "%i %le %le \n", &jx, &z0cen[j], &rocen[j]);
            for (n = 0; n <= nmaxmag; n++)
                fscanf(fp, "%i %le  \n", &nx, &Bcenaxisymm[j][n]);
        }
        fscanf(fp, " %le %le \n", &z0rem, &rorem);
        for (n = 0; n <= nmaxmag; n++)
            fscanf(fp, "%i %le  \n", &nx, &Bremaxisymm[n]);
        fclose(fp);
        iff = 1;
    }
// End of input from data file magsource_axisymm.dat. 
//
// We start now the magnetic field calculation.
// z and r coordinates of the field point P:
    z = P[3];
    r = sqrt(P[1] * P[1] + P[2] * P[2]);
// Now we have to find the best source point. 
//
// 1. We try first the central Legendre polynomial expansion.
// If some central source point has been already used 
// (jlast>0), we search the central source
// point j close to the old source point jlast.
    if (jlast == 0)
        goto label2;
// We search for the best source point (local minimum
//   of rc2=ro2/rocen2) close to jlast
    rcmin2 = 1.e9;
    jmin = jlast - 1;
    if (jmin < 1)
        jmin = 1;
    jmax = jlast + 1;
    if (jmax > Nsp)
        jmax = Nsp;
    for (j = jmin; j <= jmax; j++) {
        delz = z - z0cen[j];
        ro2 = r * r + delz * delz;
        rocen2 = rocen[j] * rocen[j];
        if (ro2 < rocen2 * rcmin2) {
            rcmin2 = ro2 / rocen2;
            jbest = j;
        }
    }
    if (jbest == jlast + 1) {
        j = jlast;
        delz = z - z0cen[j];
        rc2old = (r * r + delz * delz) / (rocen[j] * rocen[j]);
        for (j = jbest; j <= Nsp; j++) {
            delz = z - z0cen[j];
            rc2 = (r * r + delz * delz) / (rocen[j] * rocen[j]);
            if (rc2 > rc2old) {
                jbest = j - 1;
                goto label1;
            } else
                rc2old = rc2;
        }
    } else if (jbest == jlast - 1) {
        j = jlast;
        delz = z - z0cen[j];
        rc2old = (r * r + delz * delz) / (rocen[j] * rocen[j]);
        for (j = jbest; j >= 1; j--) {
            delz = z - z0cen[j];
            rc2 = (r * r + delz * delz) / (rocen[j] * rocen[j]);
            if (rc2 > rc2old) {
                jbest = j + 1;
                goto label1;
            } else
                rc2old = rc2;
        }
    }
  label1:;
    s = jbest;
    //printf("%d\n", s);
    //printf("%f %f\n", z, r);
    //printf("%f\n", z0cen[s]);
    //printf("%f\n", rocen[s]);
    magfield2_central(0, s, z0cen[s], rocen[s], z, r, &Bz, &Br, &rcp);
    if (rcp <= 0.999) {
        jlast = jbest;
        goto labelend;
    }
// -----------------------------
// 2. We try the remote series expansion:
// (with remote source point in the middle of the axisymmetric coil system)
  label2:;
    delz = z - z0rem;
    ro = sqrt(r * r + delz * delz);
    rc = rorem / ro;
    if (rc < 0.9) {
        magfield2_remote(0, 1, z0rem, rorem, z, r, &Bz, &Br, &rcp);
        if (rcp <= 0.999)
            goto labelend;
    }
// 3. We search now among all central source points.
    rcmin2 = 0.9981;
    for (j = 1; j <= Nsp; j++) {
        delz = z - z0cen[j];
        ro2 = r * r + delz * delz;
        rocen2 = rocen[j] * rocen[j];
        if (ro2 < rocen2 * rcmin2) {
            rcmin2 = ro2 / rocen2;
            jbest = j;
        }
    }
    if (rcmin2 < 0.998) {
        s = jbest;
        magfield2_central(0, s, z0cen[s], rocen[s], z, r, &Bz, &Br, &rcp);
        if (rcp <= 0.999) {
            jlast = jbest;
            goto labelend;
        }
    }
// -----------------------------
// 4. We try again the remote series expansion:
// (with remote source point in the middle of the axisymmetric coil system)
    delz = z - z0rem;
    ro = sqrt(r * r + delz * delz);
    rc = rorem / ro;
    if (rc < 0.999) {
        magfield2_remote(0, 1, z0rem, rorem, z, r, &Bz, &Br, &rcp);
        if (rcp <= 0.999)
            goto labelend;
    }
// -----------------------------
// 5. Unfortunately, no appropriate central or remote source point was
// found for the axisymmetric system. We have to compute the magnetic field
// as a sum of the field from the axisymmetric coils
//  (ia=1,...,Ncoilaxisymm).      
    for (k = 1; k <= 3; k++)
        B[k] = 0.;
    for (ia = 1; ia <= Ncoilaxisymm; ia++) {
        i = indexaxisymm[ia][0];
        magfield3_1coil(i, P, Bi);
        for (k = 1; k <= 3; k++)
            B[k] += Bi[k];
    }
    return;
// -----------------------------
// B[k],k=1,2,3 calculation from Bz, Br:
  labelend:;
    B[3] = Bz;
    if (r < 1.e-16 || fabs(Br) < fabs(Bz) * 1.e-15)
        B[1] = B[2] = 0.;
    else {
        B[1] = P[1] / r * Br;
        B[2] = P[2] / r * Br;
    }
}

// Part 5: complete elliptic integral calculations
//           (according to Numerical Recipes)

double RF_Carlson(double x, double y, double z)
{
// This function computes Carlson's elliptic integral of the first kind:
// R_F(x,y,z). x, y, z must be nonnegative, and at most one can be zero
//  (see: Press et al., Numerical Recipes, Sec. 6.11).
    const double ERRTOL = 0.002, TINY = 1.e-38, BIG = 1.e38, C1 = 1. / 24., C2 =
        0.1, C3 = 3. / 44., C4 = 1. / 14., THIRD = 1. / 3.;
    double alamb, ave, delx, dely, delz, e2, e3, sqrtx, sqrty, sqrtz, xt, yt,
        zt;
    if (FMIN3(x, y, z) < 0. || FMIN3(x + y, x + z, y + z) < TINY
        || FMAX3(x, y, z) > BIG) {
        fprintf(stderr,
                "Message from function RF_Carlson: invalid arguments !!!\n");
        fprintf(stderr, "Program running is stopped !!!\n");
        exit(1);
    }
    xt = x;
    yt = y;
    zt = z;
    do {
        sqrtx = sqrt(xt);
        sqrty = sqrt(yt);
        sqrtz = sqrt(zt);
        alamb = sqrtx * (sqrty + sqrtz) + sqrty * sqrtz;
        xt = 0.25 * (xt + alamb);
        yt = 0.25 * (yt + alamb);
        zt = 0.25 * (zt + alamb);
        ave = THIRD * (xt + yt + zt);
        delx = (ave - xt) / ave;
        dely = (ave - yt) / ave;
        delz = (ave - zt) / ave;
    }
    while (FMAX3(fabs(delx), fabs(dely), fabs(delz)) > ERRTOL);
    e2 = delx * dely - delz * delz;
    e3 = delx * dely * delz;
    return (1. + (C1 * e2 - C2 - C3 * e3) * e2 + C4 * e3) / sqrt(ave);
}

double RD_Carlson(double x, double y, double z)
{
// This function computes Carlson's elliptic integral of the second kind:
// R_D(x,y,z). x and y must be nonnegative, and at most one can be zero.
// z must be positive
//  (see: Press et al., Numerical Recipes, Sec. 6.11).
    const double ERRTOL = 0.0015, TINY = 1.e-25, BIG = 1.e22, C1 =
        3. / 14., C2 = 1. / 6., C3 = 9. / 22., C4 = 3. / 26., C5 =
        0.25 * C3, C6 = 1.5 * C4;
    double alamb, ave, delx, dely, delz, ea, eb, ec, ed, ee, fac, sum, sqrtx,
        sqrty, sqrtz, xt, yt, zt;
    if (FMIN(x, y) < 0. || FMIN(x + y, z) < TINY || FMAX3(x, y, z) > BIG) {
        fprintf(stderr,
                "Message from function RD_Carlson: invalid arguments !!!\n");
        fprintf(stderr, "Program running is stopped !!!\n");
        exit(1);
    }
    xt = x;
    yt = y;
    zt = z;
    sum = 0.;
    fac = 1.;
    do {
        sqrtx = sqrt(xt);
        sqrty = sqrt(yt);
        sqrtz = sqrt(zt);
        alamb = sqrtx * (sqrty + sqrtz) + sqrty * sqrtz;
        sum += fac / (sqrtz * (zt + alamb));
        fac = 0.25 * fac;
        xt = 0.25 * (xt + alamb);
        yt = 0.25 * (yt + alamb);
        zt = 0.25 * (zt + alamb);
        ave = 0.2 * (xt + yt + 3. * zt);
        delx = (ave - xt) / ave;
        dely = (ave - yt) / ave;
        delz = (ave - zt) / ave;
    }
    while (FMAX3(fabs(delx), fabs(dely), fabs(delz)) > ERRTOL);
    ea = delx * dely;
    eb = delz * delz;
    ec = ea - eb;
    ed = ea - 6. * eb;
    ee = ed + ec + ec;
    return 3. * sum + fac * (1. + ed * (-C1 + C5 * ed - C6 * delz * ee) +
                             delz * (C2 * ee +
                                     delz * (-C3 * ec +
                                             delz * C4 * ea))) / (ave *
                                                                  sqrt(ave));
}

double RJ_Carlson(double x, double y, double z, double p)
{
// This function computes Carlson's elliptic integral of the third kind:
// R_J(x,y,z,p). x, y and z must be nonnegative, and at most one can be zero.
// p must be nonzero. If p<0, the Cauchy principal value is returned.
//  (see: Press et al., Numerical Recipes, Sec. 6.11).
    const double ERRTOL = 0.0015, TINY = 1.e-20, BIG = 1.e12, C1 =
        3. / 14., C2 = 1. / 3., C3 = 3. / 22., C4 = 3. / 26., C5 =
        0.75 * C3, C6 = 1.5 * C4, C7 = 0.5 * C2, C8 = 2. * C3;
    double RC_Carlson(double x, double y);
    double RF_Carlson(double x, double y, double z);
    double a = 0, alamb, alpha, ans, ave, b =
        0, beta, delp, delx, dely, delz, ea, eb, ec, ed, ee, fac, pt, rcx =
        0, rho, sum, sqrtx, sqrty, sqrtz, tau, xt, yt, zt;
    if (FMIN3(x, y, z) < 0.
        || FMIN(FMIN(x + y, x + z), FMIN(y + z, fabs(p))) < TINY
        || FMAX(FMAX(x, y), FMAX(z, fabs(p))) > BIG) {
        fprintf(stderr,
                "Message from function RJ_Carlson: invalid arguments !!!\n");
        fprintf(stderr, "Program running is stopped !!!\n");
        exit(0);
    }
    sum = 0.;
    fac = 1.;
    if (p > 0.) {
        xt = x;
        yt = y;
        zt = z;
        pt = p;
    } else {
        xt = FMIN3(x, y, z);
        zt = FMAX3(x, y, z);
        yt = x + y + z - xt - zt;
        a = 1. / (yt - p);
        b = a * (zt - yt) * (yt - xt);
        pt = yt + b;
        rho = xt * zt / yt;
        tau = p * pt / yt;
        rcx = RC_Carlson(rho, tau);
    }
    do {
        sqrtx = sqrt(xt);
        sqrty = sqrt(yt);
        sqrtz = sqrt(zt);
        alamb = sqrtx * (sqrty + sqrtz) + sqrty * sqrtz;
        alpha = pow2(pt * (sqrtx + sqrty + sqrtz) + sqrtx * sqrty * sqrtz);
        beta = pt * pow2(pt + alamb);
        sum += fac * RC_Carlson(alpha, beta);;
        fac = 0.25 * fac;
        xt = 0.25 * (xt + alamb);
        yt = 0.25 * (yt + alamb);
        zt = 0.25 * (zt + alamb);
        pt = 0.25 * (pt + alamb);
        ave = 0.2 * (xt + yt + zt + 2. * pt);
        delx = (ave - xt) / ave;
        dely = (ave - yt) / ave;
        delz = (ave - zt) / ave;
        delp = (ave - pt) / ave;
    }
    while (FMAX(FMAX(fabs(delx), fabs(dely)), FMAX(fabs(delz), fabs(delp)))
           > ERRTOL);
    ea = delx * (dely + delz) + dely * delz;
    eb = delx * dely * delz;
    ec = delp * delp;
    ed = ea - 3. * ec;
    ee = eb + 2. * delp * (ea - ec);
    ans =
        3. * sum + fac * (1. + ed * (-C1 + C5 * ed - C6 * ee) +
                          eb * (C7 + delp * (-C8 + delp * C4)) +
                          delp * ea * (C2 - delp * C3) -
                          C2 * delp * ec) / (ave * sqrt(ave));
    if (p < 0.)
        ans = a * (b * ans + 3. * (rcx - RF_Carlson(xt, yt, zt)));
    return ans;
}

double RC_Carlson(double x, double y)
{
// This function computes Carlson's degenerate elliptic integral:
// R_C(x,y). x must be nonnegative, and y must be nonzero.
// If y<0, the Cauchy principal value is returned.
//  (see: Press et al., Numerical Recipes, Sec. 6.11).
    const double ERRTOL = 0.001, TINY = 1.e-38, BIG = 1.e38, SQRTNY =
        1.e-19, TNBG = TINY * BIG, COMP1 = 2.236 / SQRTNY, COMP2 =
        TNBG * TNBG / 25., THIRD = 1. / 3., C1 = 0.3, C2 = 1. / 7., C3 =
        0.375, C4 = 9. / 22.;
    double alamb, ave, s, w, xt, yt;

    if (x < 0. || y == 0. || (x + fabs(y)) < TINY || x + fabs(y) > BIG ||
        (y < -COMP1 && x > 0. && x < COMP2)) {
        fprintf(stderr,
                "Message from function RC_Carlson: invalid arguments !!!\n");
        fprintf(stderr, "Program running is stopped !!!\n");
        exit(1);
    }
    if (y > 0.) {
        xt = x;
        yt = y;
        w = 1.;
    } else {
        xt = x - y;
        yt = -y;
        w = sqrt(x) / sqrt(xt);
    }
    do {
        alamb = 2. * sqrt(xt) * sqrt(yt) + yt;
        xt = 0.25 * (xt + alamb);
        yt = 0.25 * (yt + alamb);
        ave = THIRD * (xt + yt + yt);
        s = (yt - ave) / ave;
    }
    while (fabs(s) > ERRTOL);
    return w * (1. + s * s * (C1 + s * (C2 + s * (C3 + s * C4)))) / sqrt(ave);
}

void magfield3_set_prefix(const char *prefix)
{
    char *magcoil, *magsource_central, *magsource_remote, *magsource_axisymm;
    int baselen = strlen(prefix);

    magcoil = (char*) malloc(baselen + 13);
    magsource_central = (char*) malloc(baselen + 23);
    magsource_remote = (char*) malloc(baselen + 22);
    magsource_axisymm = (char*) malloc(baselen + 23);

    strcpy(magcoil, prefix);
    strcat(magcoil, "_magcoil.dat");
    strcpy(magsource_central, prefix);
    strcat(magsource_central, "_magsource_central.dat");
    strcpy(magsource_remote, prefix);
    strcat(magsource_remote, "_magsource_remote.dat");
    strcpy(magsource_axisymm, prefix);
    strcat(magsource_axisymm, "_magsource_axiymm.dat");

    filename_magcoil = magcoil;
    filename_magsource_central = magsource_central;
    filename_magsource_remote = magsource_remote;
    filename_magsource_axisymm = magsource_axisymm;
}

void magfield3_set_magcoil(const char *filename)
{
    if (filename_magcoil)
        free(filename_magcoil);
    filename_magcoil = strdup(filename);
}

void magfield3_set_magsource_central(const char *filename)
{
    if (filename_magsource_central)
        free(filename_magsource_central);
    filename_magsource_central = strdup(filename);
}

void magfield3_set_magsource_remote(const char *filename)
{
    if (filename_magsource_remote)
        free(filename_magsource_remote);
    filename_magsource_remote = strdup(filename);
}

void magfield3_set_magsource_axisymm(const char *filename)
{
    if (filename_magsource_axisymm)
        free(filename_magsource_axisymm);
    filename_magsource_axisymm = strdup(filename);
}

const char *magfield3_get_magcoil(void)
{
    if (filename_magcoil)
        return filename_magcoil;
    else
        return "magcoil.dat";
}

const char *magfield3_get_magsource_central(void)
{
    if (filename_magsource_central)
        return filename_magsource_central;
    else
        return "magsource_central.dat";
}

const char *magfield3_get_magsource_remote(void)
{
    if (filename_magsource_remote)
        return filename_magsource_remote;
    else
        return "magsource_remote.dat";
}

const char *magfield3_get_magsource_axisymm(void)
{
    if (filename_magsource_axisymm)
        return filename_magsource_axisymm;
    else
        return "magsource_axisymm.dat";
}

void magfield3_set_nmaxmag(int n)
{
    nmaxmag = n;
}

int magfield3_get_nmaxmag(void)
{
    return nmaxmag;
}

/*

The coil parameters:

double coil[Ncoilmax+1][14] matrix:
     maximal number of coils: Ncoilmax
     number of coils: Ncoil
     coil index : i
     coil[i][0] : current density (A/m^2)
     coil[i][1],coil[i][2],coil[i][3]: endpoint A of the coil axis (m)
     coil[i][4],coil[i][5],coil[i][6]: endpoint B of the coil axis (m)
     coil[i][7] : inner radius of the coil winding (m)
     coil[i][8] : outer radius of the coil winding (m)
     coil[i][9] : radial integration number (positive integer; it
                    is represented here by double)
     coil[i][10]: length of coil (distance of points A and B) (m)
     coil[i][11],coil[i][12],coil[i][13]: unit vector in A --> B
                                          direction (m)
     
Global integer parameters:

  Ncoil: number of coils
  Ncoilaxisymm: number of axisymmetric coils (it is between 0 and Ncoil)
  Ncoilmax: maximal number of coils (determined by #define)
  nmaxmag: maximal source coefficient index (determined by #define)
  Ncenmax: maximal number of central source points of all coils 
           (determined by #define)
  Ncenaxisymmmax: maximal number of central source points of the 
                  axisymmetric coil system (determined by #define)
  
There are 2 different types of coils: axisymmetric (with 
 xA=xB=yA=yB=0.), and not-axisymmetric. Axis of the axisymmetric coils: z.
 



*/
