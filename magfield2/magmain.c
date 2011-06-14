/* $Id: magmain.c 231 2007-07-26 11:43:35Z s_voec01 $ */

// magmain.c program

// changes made by K. Valerius, 19. Jan. 2006:
// 1) add request for file i/o success 

#define _GNU_SOURCE

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#include "magfield2.h"

void outputcoil()
{
    FILE *fp;
    extern int Ncoil;
    extern double **coil;
    double zmin, zmax, rmin, rmax;
    int i;

// Output to file 'coil.dat':
    //fp=fopen("coil.dat","w");
    if ((fp = fopen("coil.dat", "w")) == NULL) {
        fprintf(stderr, "Error opening file coil.dat \n");
        exit(1);
    }
    for (i = 1; i <= Ncoil; i++) {
        zmin = coil[i][1];
        zmax = coil[i][2];
        rmin = coil[i][3];
        rmax = coil[i][4];
        fprintf(fp, "%12.4f %12.4f %12.4f %12.4f  \n", zmin, zmax, rmin, rmax);
    }
    fclose(fp);
}

void test(char *inputcoil)
{
    double r, z, A, Bz, Br, B;
    int j, i, n = 20;
    double zmin = -0.05, zmax = 0.05, zstep = 0.00001;
    int imax = (zmax - zmin) / zstep;
    double rmin = 0.0000, rmax = 0.001, rstep = 0.0005;
    int jmax = (rmax - rmin) / rstep;
    FILE *file;
    char filename[20];
 
    for (j = 0; j <= jmax; j++) {
      sprintf(filename, "magnetfeld%d.dat", j);
      file = fopen(filename, "w");
      for (i = 0; i <= imax; i++) {
        z = zmin + i * zstep;
        r = rmin + j * rstep;

        if (i % 100)
            printf("%d\n", i);

        magfield2_field(z, r, inputcoil, n, &A, &Bz, &Br);
        B = sqrt(Bz * Bz + Br * Br);
        fprintf(file, "%f %f %f %f %f\n", r, z, B, Bz, Br);
      }
      fclose(file);
    }
    /*double z,r,A,Bz0,Bz,Br,B;
       int i,n;
       n=20;

       z=-17.45;
       r=0.0;
       magfield2(z,r,inputcoil,n,&A,&Bz,&Br);
       printf("z,r,Bz,Br=%19.11e  %19.11e %19.11e %19.11e \n \n",z,r,Bz,Br);
       z=-15.3;
       r=0.;
       magfield2(z,r,inputcoil,n,&A,&Bz,&Br);
       printf("z,r,Bz,Br=%19.11e  %19.11e %19.11e %19.11e \n \n",z,r,Bz,Br);
       z=-14.65;
       r=0.7;
       magfield2(z,r,inputcoil,n,&A,&Bz,&Br);
       printf("z,r,Bz,Br=%19.11e  %19.11e %19.11e %19.11e \n \n",z,r,Bz,Br);
       z=-14.65;
       r=1.5;
       magfield2(z,r,inputcoil,n,&A,&Bz,&Br);
       printf("z,r,Bz,Br=%19.11e  %19.11e %19.11e %19.11e \n \n",z,r,Bz,Br);
       //    magfield2_elliptic(z,r,inputcoil,n,&A,&Bz,&Br);
       //    printf("A,Bz,Br= %19.11e %19.11e %19.11e  \n",A,Bz,Br);
       z=0.;
       r=4.5;
       magfield2(z,r,inputcoil,n,&A,&Bz,&Br);
       printf("A,Bz,Br= %19.11e %19.11e %19.11e \n \n",A,Bz,Br);
       //    magfield2_elliptic(z,r,inputcoil,n,&A,&Bz,&Br);
       //    printf("A,Bz,Br= %19.11e %19.11e %19.11e  \n",A,Bz,Br);
       //  return;
       r=0.;
       for(i=-19*20;i<=18*20;i++)
       {
       z=i*0.05;
       magfield2(z,0.,inputcoil,n,&A,&Bz0,&Br);
       magfield2(z,4.5,inputcoil,n,&A,&Bz,&Br);
       printf("z,Bz0,Bz= %12.5f %14.5e %14.5e    \t\n",z,Bz0,Bz); 
       }
       r=0.301;
       for(i=294;i<=320;i++)
       {
       z=i*0.05;
       magfield2(z,r,inputcoil,n,&A,&Bz,&Br);
       B=sqrt(Bz*Bz+Br*Br);
       printf("z,r,B= %12.5f%12.5f  %17.5e    \t\n",z,r,B); 
       }
       r=0.32;
       for(i=294;i<=310;i++)
       {
       z=i*0.05;
       magfield2(z,r,inputcoil,n,&A,&Bz,&Br);
       B=sqrt(Bz*Bz+Br*Br);
       printf("z,r,B= %12.5f%12.5f  %17.5e    \t\n",z,r,B); 
       }
       return; */
}
char* strndupCus(char* str, int n) {
  // Stupid implementation of strndup since macos isn't born with
  // one.
  int len = strlen(str);
  if (len <= n)
    return strdup(str);
  char* result = (char*) malloc(n+1);
  int i;
  for (i = 0; i <= n; i++)
    result[i] = str[i];
  result[i] = '\0';
  return result;
}

int main(int argc, char **argv)
{
    double z0min, z0max, delz0;
    char *inputcoil;
    int c, run_test, len;

    z0min = -21.;
    z0max = 17.;
    delz0 = 0.04;
    run_test = 0;
    
    while((c = getopt(argc, argv, "m:M:d:tn:")) != -1) {
        switch(c) {
            case 'm':
            z0min = strtod(optarg, NULL);
            break;
            
            case 'M':
            z0max = strtod(optarg, NULL);
            break;
            
            case 'd':
            delz0 = strtod(optarg, NULL);
            break;
            
            case 't':
            run_test = 1;
            break;
            
            case 'n':
            magfield2_set_nmaxmag(strtol(optarg, NULL, 0));
            break;
        }
    }
    
    if (argc - optind < 1) {
        fprintf(stderr, "Syntax: magmain [options] <inputfile>\n");
        return 1;
    }
    inputcoil = argv[optind];
    
    len = strlen(inputcoil);
    if (len > 4 && !strcmp(inputcoil + len - 4, ".dat")) {
        char* prefix = strndupCus(inputcoil, len - 4);
        char* magsource = (char*) malloc(len + 15);
        
        strcpy(magsource, prefix);
        strcat(magsource, "_magsource.dat");
        magfield2_set_magsource(magsource);
        
        free(prefix);
        free(magsource);
    }
    
    magfield2_source(z0min, z0max, delz0, inputcoil);
    
    if(run_test)
        test(inputcoil);
    outputcoil();
    return 0;
}


