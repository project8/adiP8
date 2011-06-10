/* $Id: magfield2.c 231 2007-07-26 11:43:35Z s_voec01 $ */
////////////////////////////////////////////////////////////////////////
//                                                                    //
//     magfield2.c : axially symmetric magnetic field calculation     //
//                                                                    //
////////////////////////////////////////////////////////////////////////

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "magfield2.h"

// changes made by K. Valerius, 19. Jan. 2006:
// 1) add check for file i/o success

#define FMIN(a,b) ((a)<=(b)?(a):(b))
#define FMAX(a,b) ((a)>(b)?(a):(b))
#define FMIN3(a,b,c) (FMIN(a,b)<=(c)?(FMIN(a,b)):(c))
#define FMAX3(a,b,c) (FMAX(a,b)>(c)?(FMAX(a,b)):(c))
#define pow2(x) ((x)*(x))

static void coilread(const char *inputcoil);
static void input_sourcepoints_mag(double z0min, double z0max, double delz0,
                            double **z0vec);
static void magfield2_elliptic(double z, double r, const char *inputcoil, int n, double *A,
                        double *Bz, double *Br);
static void magfield_elliptic_1coil(int n, double *coilpar, double z, double r,
                             double *A, double *Bz, double *Br);
static double RF_Carlson(double x, double y, double z);
static double RD_Carlson(double x, double y, double z);
static double RJ_Carlson(double x, double y, double z, double p);
static double RC_Carlson(double x, double y);

static int nmaxmag = 200;

static char *filename_magsource = "magsource.dat";

void magfield2_set_nmaxmag(int n)
{
    nmaxmag = n;
}

int magfield2_get_nmaxmag(void)
{
    return nmaxmag;
}

void magfield2_set_magsource(const char* filename)
{
    filename_magsource = strdup(filename);
}

const char* magfield2_get_magsource()
{
    return filename_magsource;
}

// Source coefficient calculation:

void magfield2_source(double z0min, double z0max, double delz0, const char *inputcoil)
// This subroutine computes the magnetic field source constants B[n]
//  (n=0,...,nmaxmag), in the axis source points z0=z0vec[k] (k=1,...,Nspmag).
//  *z0vec: the axis source point vector.
// It needs the global integers Ncoil, nmaxmag and Nspmag, and the matrix
//   coil[Ncoilmax+1][5].
//  Nspmag: number of z0 source points;
//  nmaxmag: maximal index of the source constants (maximum of n).
//  nmaxmag is given by #define command, Nspmag by global variable.
//
// Output into data file magsource.dat:
//    Nspmag, nmaxmag and k number-groups;
//    numbers in number-group k:
//              index k, z0[k], rocen, B[n] (n=0,...,nmaxmag)
{
    //double z0vec[Nspmagmax + 1];
    double *z0vec = NULL;
    double *B, **b, **bhat, *Pp;
// Indices: i: coil; k: source points; m: numerical integration;
//          n: source constants
//   iz=0: z=zmax; iz=1: z=zmin
    int i, k, m, n, iz;
    const int M = 30;
    const double w9[10] = { 0.2803440531305107e0, 0.1648702325837748e1,
        -0.2027449845679092e0, 0.2797927414021179e1,
        -0.9761199294532843e0, 0.2556499393738999e1,
        0.1451083002645404e0, 0.1311227127425048e1,
        0.9324249063051143e0, 0.1006631393298060e1
    };
    static double w[31];
    double z0, rocen, ro, Z, R, del;
    double q;
    double zmin, zmax, rmin, rmax, current, z, u, constcen, sigmac, rcen, rcen1;
    const double mu0 = 4. * M_PI * 1.e-7;
    FILE *fp;
    
    B = (double*) malloc(sizeof(double) * (nmaxmag + 1));
    b = (double**) malloc(sizeof(double*) * (nmaxmag + 1));
    bhat = (double**) malloc(sizeof(double*) * (nmaxmag + 1));
    Pp = (double*) malloc(sizeof(double) * (nmaxmag + 2));
    
    for (i = 0; i <= nmaxmag; i++) {
        b[i] = (double*) malloc(sizeof(double) * 3);
        bhat[i] = (double*) malloc(sizeof(double) * 31);
    }

// Reading the coil parameters:
    coilread(inputcoil);
// Calculation of source points:
    input_sourcepoints_mag(z0min, z0max, delz0, &z0vec);
// Output to file magsource.dat (Nspmag,nmaxmag)
    //fp=fopen("magsource.dat","w");
    //better check file opening success:
    if ((fp = fopen(filename_magsource, "w")) == NULL) {
        fprintf(stderr, "Error opening file '%s'for write access!\n",
                filename_magsource);
        exit(1);
    }
    fprintf(fp, "%9i %9i \n", Nspmag, nmaxmag);
// Initialization of the integration weight factors:
    for (m = 0; m <= 9; m++)
        w[m] = w9[m];
    for (m = 10; m <= M - 10; m++)
        w[m] = 1.;
    for (m = M - 9; m <= M; m++)
        w[m] = w9[M - m];
//
// Source point loop:
    for (k = 1; k <= Nspmag; k++) {
        z0 = z0vec[k];
// Initialization of B[n] and Pp[n]:
        for (n = 0; n <= nmaxmag; n++)
            B[n] = 0.;
        Pp[0] = 0.;
        Pp[1] = 1.;
// Calc. of rocen = minimum distance of the axis point z0 from the
//                  coils
        rocen = 1.e20;
        for (i = 1; i <= Ncoil; i++) {
            zmin = coil[i][1];
            zmax = coil[i][2];
            rmin = coil[i][3];
            rmax = coil[i][4];
// Coil geometry test:
            if (zmin >= zmax || rmin >= rmax || rmin <= 0.) {
                printf("Message from subroutine magsource:\
               zmin>=zmax or rmin>=rmax or rmin<=0.!!!\
               Computation is  stopped !!! \n");
                exit(0);
            }
            if (z0 <= zmin)
                ro = sqrt((z0 - zmin) * (z0 - zmin) + rmin * rmin);
            else if (z0 >= zmax)
                ro = sqrt((z0 - zmax) * (z0 - zmax) + rmin * rmin);
            else
                ro = rmin;
            if (ro < rocen)
                rocen = ro;
        }
// B[n] calculation
// Coil index loop:
        for (i = 1; i <= Ncoil; i++) {
            zmin = coil[i][1];
            zmax = coil[i][2];
            rmin = coil[i][3];
            rmax = coil[i][4];
            current = coil[i][0];
            del = (rmax - rmin) / M;
            sigmac = current / (zmax - zmin) / (rmax - rmin);
// Integration loop:
            for (m = 0; m <= M; m++) {
                R = rmin + del * m;
                for (iz = 1; iz <= 2; iz++) {
                    if (iz == 1)
                        Z = zmax;
                    else
                        Z = zmin;
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
                Z = zmax;
                z = Z - z0;
                ro = sqrt(R * R + z * z);
                bhat[0][m] = mu0 * sigmac / 2. * z / ro;
                Z = zmin;
                z = Z - z0;
                ro = sqrt(R * R + z * z);
                bhat[0][m] -= mu0 * sigmac / 2. * z / ro;
                for (n = 1; n <= nmaxmag; n++)
                    bhat[n][m] = -rocen / n * (b[n - 1][1] - b[n - 1][2]);
            }                   // end of m-loop
            for (n = 0; n <= nmaxmag; n++) {
                q = 0.;
                for (m = 0; m <= M; m++)
                    q += w[m] * bhat[n][m];
                q *= del;
                B[n] += q;
            }
        }                       // end of i-loop
// Output to file magsource.dat (k, z0, rocen, B[n])
        fprintf(fp, "%9i \n", k);
        fprintf(fp, "%22.14e %22.14e \n", z0, rocen);
        for (n = 0; n <= nmaxmag; n++) {
            if (fabs(B[n]) < 1.e-30)
                B[n] = 0.;
            fprintf(fp, "%22.14e  \n", B[n]);
        }
    }                           // end of k-loop
    fclose(fp);
    
    for (i = 0; i <= nmaxmag; i++) {
        free(b[i]);
        free(bhat[i]);
    }
    free(z0vec);
    free(B);
    free(b);
    free(Pp);
    free(bhat);
}

// Reading from file given by string inputcoil:

void coilread(const char *inputcoil)
{
// This function reads the coil parameters (current and geometry)
// from the data file defined by the string inputcoil.
// The coils have rectangular shape cross section in the (z,r)
//   meridian plane.
// The data in the file defined by string inputcoil are:
//  First line: number of coils  (Ncoil)
// Then there are Ncoil number of lines; each line contains:
//    zmid  rin  thick  length  current
//    zmid:  middle z value of the coil (zmid=(zmin+zmax)/2)  (m)
//    rin:  inner radius of coil (rin=rmin)  (m)
//    thick:  radial thickness of coil (thick=rmax-rmin)  (m)
//    length:  length of coil in z direction (length=zmax-zmin)  (m)
//    current: total coil current =  turn number * current in wire
//       (Amper)
// The coil parameters are written into the array coil[Ncoilmax+1][5].
//   i: coil index
//   coil[i][0]: total coil current =  turn number * current in wire (Amper)
//   coil[i][1] = zmin: minimal z value of the coil (meter)
//   coil[i][2] = zmax: maximal z value of the coil (meter)
//   coil[i][3] = rmin: minimal r value of the coil (meter)
//   coil[i][4] = rmax: maximal r value of the coil (meter)
//  SI units are used !
    int i;
    double zmid, rin, thick, length, current;
    FILE *fp;

//
    fp = fopen(inputcoil, "r");
    if (!fp) {
        puts("Message from function coilread:");
        puts("Cannot open the coil input file!");
        puts("Program running is stopped !!! ");
        exit(0);
    }
    
    if(coil) {
        for(i = 1; i <= Ncoil; i++)
            free(coil[i]);
        free(coil); 
    }
    
    fscanf(fp, "%i", &Ncoil);
    
    coil = (double**) malloc(sizeof(double) * (Ncoil + 1));

    printf("\nPrint current coil parameters list:");
    printf("\n===================================\n\n");
    printf
        ("coil #: zmid [m]\trin [m]\t\t thick [m]\tlength [m]\tampturns [A] \n");

    for (i = 1; i <= Ncoil; i++) {
        fscanf(fp, "%le %le %le %le %le ", &zmid, &rin, &thick, &length,
               &current);
        printf("#%i: \t %3.4f \t %3.3f \t\t %3.3f \t\t %3.3f \t\t %9.3f \n", i,
               zmid, rin, thick, length, current);

        if (rin <= 0. || thick <= 0. || length <= 0.) {
            puts("Warning message from function coilread:");
            puts("non-positive coil parameters are changed !!!");
        }
        if (rin <= 0.)
            rin = 1.e-12;
        if (thick <= 0.)
            thick = 1.e-12;
        if (length <= 0.)
            length = 1.e-12;
        coil[i] = (double*) malloc(sizeof(double) * 5);
        coil[i][0] = current;
        coil[i][1] = (zmid - length / 2.);
        coil[i][2] = (zmid + length / 2.);
        coil[i][3] = rin;
        coil[i][4] = (rin + thick);
    }
    fclose(fp);
}

// Source point calculation:

void input_sourcepoints_mag(double z0min, double z0max, double delz0,
                            double **z0vec)
{
// The z0vec[] source points are defined here.
// k: index of the source points
// We use here the simplest case: equidistant source point distribution.
// Input parameters:
//   z0min,z0max: minimal and maximal z values of the source points
//   delz0: distance between neighbouring source points
// Output:
//   z values of source points, z0vec[k], k=1,...,Nspmag
// z0vec: source point vector
// Nspmag: number of source points
// k: index of the source points
// SI units: meter !
    int k;

// Nspmag: number of source points
    Nspmag = (z0max - z0min) / delz0 + 1;
    *z0vec = (double*) malloc(sizeof(double) * (Nspmag + 1));
    for (k = 1; k <= Nspmag; k++)
        (*z0vec)[k] = z0min + delz0 * (k - 1);
}

// Magnetic field calculation by Legendre-polynomial expansion:

void magfield2_field(double z, double r, const char *inputcoil, int nn,
               double *A, double *Bz, double *Br)
// This subroutine computes the axial and radial magnetic
//    field components Bz and Br, and the vectorpotential A
//    in a point with cylindrical coordinates z and r.
// Method: Legendre polynomial expansion around the source point with index k.
//  Nspmagmax: we need Nspmagmax>=Nspmag
//  nmaxmag: maximal index of the source constants (maximum of n).
//  Nspmagmax and nmaxmag are given by global #define commands.
//  Important (if magfield2 is used separated from magsource):
//    the same nmaxmag number used for magsource
//    should also be used for magfield2!!!
//  SI units are used !
{
    FILE *fp;
    int kloop, kx, n, k, kk;
    static double *z0 = NULL, *rocen = NULL, **B = NULL; 
    static double *c1 = NULL, *c2 = NULL, *c3 = NULL, *c4 = NULL, *c5 = NULL;
    static double *c6 = NULL; 
    static int Nspmag, klast;
    static int iff = 0;
    static double rclimit;
    double ro, u, delz, s, rcmin;
    double *P = NULL, *Pp = NULL, rc, rcn;
    int iA, iBz, iBr;
    double *Aplus = NULL, *Bzplus = NULL, *Brplus = NULL;
    double A1, Bz1, Br1;

//
// Input from file magsource.dat:
    if (iff == 0) {
        //fp=fopen("magsource.dat","r");
        if ((fp = fopen(filename_magsource, "r")) == NULL) {
            fprintf(stderr, "Error opening file '%s' for read access\n",
                    filename_magsource);
            exit(1);
        }
        fscanf(fp, "%i %i", &Nspmag, &nmaxmag);

        z0 = (double*) malloc(sizeof(double) * (Nspmag + 1));
        rocen = (double*) malloc(sizeof(double) * (Nspmag + 1));
        B = (double**) malloc(sizeof(double*) * (Nspmag + 1));
        c1 = (double*) malloc(sizeof(double) * (nmaxmag + 1));
        c2 = (double*) malloc(sizeof(double) * (nmaxmag + 1));
        c3 = (double*) malloc(sizeof(double) * (nmaxmag + 1));
        c4 = (double*) malloc(sizeof(double) * (nmaxmag + 1));
        c5 = (double*) malloc(sizeof(double) * (nmaxmag + 1));
        c6 = (double*) malloc(sizeof(double) * (nmaxmag + 1));
        for (k = 1; k <= Nspmag; k++)
            B[k] = (double*) malloc(sizeof(double) * (nmaxmag + 1)); 
        for (kloop = 1; kloop <= Nspmag; kloop++) {
            fscanf(fp, "%i", &kx);
            fscanf(fp, "%le %le", &z0[kloop], &rocen[kloop]);
            for (n = 0; n <= nmaxmag; n++)
                fscanf(fp, "%le", &B[kloop][n]);
        }
        fclose(fp);
// Initialization of c1,c2,c3,c4,c5,c6 vectors:
        for (n = 2; n <= nmaxmag; n++) {
            c1[n] = (2. * n - 1.) / (1. * n);
            c2[n] = (n - 1.) / (1. * n);
            c3[n] = (2. * n - 1.) / (1. * (n - 1.));
            c4[n] = (1. * n) / (1. * (n - 1.));
            c5[n] = (1.) / (n * (n + 1.));
            c6[n] = (1.) / (n + 1.);
        }
// The best source point is searched here (with minimal
//    convergence ratio rc):
        rcmin = 1.e20;
        for (k = 1; k <= Nspmag; k++) {
            delz = z - z0[k];
            ro = sqrt(r * r + delz * delz);
            rc = ro / rocen[k];
            if (rc < rcmin) {
                rcmin = rc;
                klast = k;
            }
        }
// End of source point searching
// Calculation of rclimit: 
        if(nmaxmag < 100)
        {
            fprintf(stderr, "Message from subroutine magfield2:\n"
                            "nmaxmag<100\n"
                            "Computation is stopped!!!\n");
            exit(1);
        }
        else if(nmaxmag >= 100 && nmaxmag <= 150)
            rclimit = 0.9;
        else if(nmaxmag > 150 && nmaxmag <= 200)
            rclimit = 0.95;
        else if(nmaxmag > 200 && nmaxmag <= 350)
            rclimit = 0.97;
        else if(nmaxmag > 350 && nmaxmag <= 500)
            rclimit = 0.99;
        else if(nmaxmag > 500 && nmaxmag <= 650)
            rclimit = 0.993;
        else if(nmaxmag > 650 && nmaxmag <= 800)
            rclimit = 0.996;
        else 
            rclimit = 0.998;

        iff = 1;
    }
    P = (double*) malloc(sizeof(double) * (nmaxmag + 1));
    Pp = (double*) malloc(sizeof(double) * (nmaxmag + 1));
    Aplus = (double*) malloc(sizeof(double) * (nmaxmag + 1));
    Brplus = (double*) malloc(sizeof(double) * (nmaxmag + 1));
    Bzplus = (double*) malloc(sizeof(double) * (nmaxmag + 1));
//
// The best source point is searched here
//   (starting from the last source point)
    k = klast;
    delz = z - z0[k];
    ro = sqrt(r * r + delz * delz);
    rcmin = ro / rocen[k];
    kk = k + 1;
    if (kk <= Nspmag) {
        delz = z - z0[kk];
        ro = sqrt(r * r + delz * delz);
        rc = ro / rocen[kk];
        if (rc < rcmin) {
            rcmin = rc;
            k = kk;
        }
    }
    kk = klast - 1;
    if (kk >= 1) {
        delz = z - z0[kk];
        ro = sqrt(r * r + delz * delz);
        rc = ro / rocen[kk];
        if (rc < rcmin)
            k = kk;
    }
    klast = k;
// If rc>0.999: new searching:
    delz = z - z0[k];
    ro = sqrt(r * r + delz * delz);
    rc = ro / rocen[k];
    if (rc > rclimit) {
        rcmin = 1.e20;
        for (k = 1; k <= Nspmag; k++) {
            delz = z - z0[k];
            ro = sqrt(r * r + delz * delz);
            rc = ro / rocen[k];
            if (rc < rcmin) {
                rcmin = rc;
                klast = k;
            }
        }
        k = klast;
    }
// End of source point searching
//////////////////////////////////////
// If the field point is very close to the source point:
    if (r < 1.e-12 && fabs(z - z0[k]) < 1.e-12) {
        *A = 0.;
        *Bz = B[k][0];
        *Br = 0.;
        free(P);
        free(Pp);
        free(Aplus);
        free(Brplus);
        free(Bzplus);
        return;
    }
// ro,u,s,rc,rcn:
    delz = z - z0[k];
    ro = sqrt(r * r + delz * delz);
    u = delz / ro;
    s = r / ro;
    rc = ro / rocen[k];         // convergence ratio
    rcn = rc;
// If rc>rclimit: computation by elliptic integrals
//   (the Legendre polynomial series is not convergent)
    if (rc > rclimit) {
        magfield2_elliptic(z, r, inputcoil, nn, A, Bz, Br);
        free(P);
        free(Pp);
        free(Aplus);
        free(Brplus);
        free(Bzplus);
        return;
    }
// First 2 terms of Legendre polynomial P and its derivative Pp (P-primed)
    P[0] = 1.;
    P[1] = u;
    Pp[0] = 0.;
    Pp[1] = 1.;
// First 2 terms of the series:
    *A = s * rocen[k] * B[k][0] / 2. * rc;
    *Bz = B[k][0] + B[k][1] * rc * u;
    *Br = -s * B[k][1] / 2. * rc;
//
    iA = 0;
    iBz = 0;
    iBr = 0;
    Aplus[1] = 1.e30;
    Bzplus[1] = 1.e30;
    Brplus[1] = 1.e30;
// We start here the series expansion:
    for (n = 2; n <= nmaxmag - 1; n++) {
        rcn *= rc;
        P[n] = c1[n] * u * P[n - 1] - c2[n] * P[n - 2];
        Pp[n] = c3[n] * u * Pp[n - 1] - c4[n] * Pp[n - 2];
        Aplus[n] = rocen[k] * s * B[k][n - 1] * c5[n] * rcn * Pp[n];
        Bzplus[n] = B[k][n] * rcn * P[n];
        Brplus[n] = -s * B[k][n] * c6[n] * rcn * Pp[n];
        *A += Aplus[n];
        *Bz += Bzplus[n];
        *Br += Brplus[n];
        A1 = 1.e-15 * fabs(*A);
        Bz1 = 1.e-15 * fabs(*Bz);
        Br1 = 1.e-15 * fabs(*Br);
        if (n > 8) {
            if (fabs(Aplus[n]) < A1 && fabs(Aplus[n - 1]) < A1)
                iA = 1;
            if (fabs(Bzplus[n]) < Bz1 && fabs(Bzplus[n - 1]) < Bz1)
                iBz = 1;
            if (fabs(Brplus[n]) < Br1 && fabs(Brplus[n - 1]) < Br1)
                iBr = 1;
            if (fabs(*Br) < 1.e-12)
                iBr = 1;
        }
        if (iA * iBz * iBr == 1)
            break;
    }
    free(P);
    free(Pp);
    free(Aplus);
    free(Bzplus);
    free(Brplus);
}

// Magnetic field calculation by elliptic integrals:

void magfield2_elliptic(double z, double r, const char *inputcoil, int n,
                        double *A, double *Bz, double *Br)
// This subroutine computes the axial and radial magnetic
// field components Bz and Br, and the vectorpotential A
// in a point with cylindrical coordinates z and r, by using
// the first, second and third complete elliptic integrals.
// SI units are used !
{
    static int iff = 0;
    int i, k;
    double coilpar[5], Ai, Bzi, Bri;

// Reading the coil parameters:
    if (iff == 0) {
        coilread(inputcoil);
        iff = 1;
    }
// Summing the contributions of all coils:
    *A = 0.;
    *Bz = 0.;
    *Br = 0.;
    for (i = 1; i <= Ncoil; i++) {
        for (k = 0; k <= 4; k++)
            coilpar[k] = coil[i][k];
        magfield_elliptic_1coil(n, coilpar, z, r, &Ai, &Bzi, &Bri);
        *A += Ai;
        *Bz += Bzi;
        *Br += Bri;
    }
}

void magfield_elliptic_1coil(int n, double *coilpar, double z, double r,
                             double *A, double *Bz, double *Br)
{
    double Zmin, Zmax, Rmin, Rmax, current;
    double R, delR[3], Z, delr2, sumr2, delz2, eta, d, K, EK, PIK, S;
    double xA, xBz, xBr, sign, c, sigma, st, delRr;
    double Rlow[3], Rhigh[3];
    const double mu0 = 4. * M_PI * 1.e-7;
    int i, iR, iZ, M, m;
    double w[1001];
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
    current = coilpar[0];
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
    if (fabs(z - Zmin) < 1.e-8)
        z = Zmin - 1.e-8;
    if (fabs(z - Zmax) < 1.e-8)
        z = Zmax + 1.e-8;
// Current density:
    sigma = current / (Zmax - Zmin) / (Rmax - Rmin);
// Improvement of n value:
    if (n < 0)
        n = 0;
    if (n > 1000)
        n = 1000;
    if (n > 0 && n < 12)
        n = 12;
// Integration weight factors:
    if (n == 0)
        w[0] = 1.;
    if (n > 0 && n < 20) {
        for (i = 0; i <= 5; i++)
            w[i] = w5[i];
        for (i = 6; i <= n - 6; i++)
            w[i] = 1.;
        for (i = n - 5; i <= n; i++)
            w[i] = w5[n - i];
    } else if (n > 0 && n >= 20) {
        for (i = 0; i <= 9; i++)
            w[i] = w9[i];
        for (i = 10; i <= n - 10; i++)
            w[i] = 1.;
        for (i = n - 9; i <= n; i++)
            w[i] = w9[n - i];
    }
//
    xA = 0.;
    xBz = 0.;
    xBr = 0.;
// R-integration limits:
    if (z > Zmin && z < Zmax && r > Rmin && r < Rmax) {
        M = 2;
        Rlow[1] = Rmin;
        Rhigh[1] = r - (r - Rmin) * 1.e-12;
        Rlow[2] = r + (Rmax - r) * 1.e-12;
        Rhigh[2] = Rmax;
        if (n > 0) {
            delR[1] = (Rhigh[1] - Rlow[1]) / n;
            delR[2] = (Rhigh[2] - Rlow[2]) / n;
        } else {
            delR[1] = Rhigh[1] - Rlow[1];
            delR[2] = Rhigh[2] - Rlow[2];
        }
    } else {
        M = 1;
        Rlow[1] = Rmin;
        Rhigh[1] = Rmax;
        if (n > 0)
            delR[1] = (Rhigh[1] - Rlow[1]) / n;
        else
            delR[1] = Rhigh[1] - Rlow[1];
    }
// Integration:
    for (m = 1; m <= M; m++) {
        for (iR = 0; iR <= n; iR++) {
            if (n == 0)
                R = (Rlow[m] + Rhigh[m]) / 2.;
            else
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
                xA += w[iR] * sign * R * (z - Z) / S * (EK + d * PIK) * delR[m];
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
    *A = c * xA;
    *Bz = c * xBz;
    *Br = c * xBr;
}

double RF_Carlson(double x, double y, double z)
{
// This function computes Carlson's elliptic integral of the first kind:
// R_F(x,y,z). x, y, z must be nonnegative, and at most one can be zero.
    const double ERRTOL = 0.002, TINY = 1.e-38, BIG = 1.e38, C1 = 1. / 24., C2 =
        0.1, C3 = 3. / 44., C4 = 1. / 14., THIRD = 1. / 3.;
    double alamb, ave, delx, dely, delz, e2, e3, sqrtx, sqrty, sqrtz, xt, yt,
        zt;
    if (FMIN3(x, y, z) < 0. || FMIN3(x + y, x + z, y + z) < TINY
        || FMAX3(x, y, z) > BIG) {
        puts("Message from function RF_Carlson: invalid arguments !!!");
        puts("Program running is stopped !!! ");
        exit(0);
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
// z must be positive.
    const double ERRTOL = 0.0015, TINY = 1.e-25, BIG = 1.e22, C1 =
        3. / 14., C2 = 1. / 6., C3 = 9. / 22., C4 = 3. / 26., C5 =
        0.25 * C3, C6 = 1.5 * C4;
    double alamb, ave, delx, dely, delz, ea, eb, ec, ed, ee, fac, sum, sqrtx,
        sqrty, sqrtz, xt, yt, zt;
    if (FMIN(x, y) < 0. || FMIN(x + y, z) < TINY || FMAX3(x, y, z) > BIG) {
        puts("Message from function RD_Carlson: invalid arguments !!!");
        puts("Program running is stopped !!! ");
        exit(0);
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
// R_J(x,y,z,p). x,y and z must be nonnegative, and at most one can be zero.
// p must be nonzero. If p<0, the Cauchy principal value is returned.
    const double ERRTOL = 0.0015, TINY = 1.e-20, BIG = 1.e12, C1 =
        3. / 14., C2 = 1. / 3., C3 = 3. / 22., C4 = 3. / 26., C5 =
        0.75 * C3, C6 = 1.5 * C4, C7 = 0.5 * C2, C8 = 2. * C3;
    double RC_Carlson(double x, double y);
    double RF_Carlson(double x, double y, double z);
    double a = 0., alamb, alpha, ans, ave, b = 0., beta, delp, delx, dely, delz, ea, eb,
        ec, ed, ee, fac, pt, rcx = 0., rho, sum, sqrtx, sqrty, sqrtz, tau, xt, yt,
        zt;
    if (FMIN3(x, y, z) < 0.
        || FMIN(FMIN(x + y, x + z), FMIN(y + z, fabs(p))) < TINY
        || FMAX(FMAX(x, y), FMAX(z, fabs(p))) > BIG) {
        puts("Message from function RJ_Carlson: invalid arguments !!!");
        puts("Program running is stopped !!! ");
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
    const double ERRTOL = 0.001, TINY = 1.e-38, BIG = 1.e38, SQRTNY =
        1.e-19, TNBG = TINY * BIG, COMP1 = 2.236 / SQRTNY, COMP2 =
        TNBG * TNBG / 25., THIRD = 1. / 3., C1 = 0.3, C2 = 1. / 7., C3 =
        0.375, C4 = 9. / 22.;
    double alamb, ave, s, w, xt, yt;

    if (x < 0. || y == 0 || (x + fabs(y)) < TINY || x + fabs(y) > BIG ||
        (y < -COMP1 && x > 0. && x < COMP2)) {
        puts("Message from function RC_Carlson: invalid arguments !!!");
        puts("Program running is stopped !!! ");
        exit(0);
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
