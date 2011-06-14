using namespace std;
#include <cstdio>

void write_matrix_value(void* matrix_ptr, int nx, int ny, int nz, int first, int x, int y, int z, double value)
{
  double* matrix;
  
  matrix = (double*)matrix_ptr;
  
  *(matrix + ((nz != 0)&&(nz != 1))*(z-first)*(ny*nx) + (y-first)*nx + (x-first)) = value;
}

double read_matrix_value(void* matrix_ptr, int nx, int ny, int nz, int first, int x, int y, int z)
{
  double* matrix;
  
  matrix = (double*)matrix_ptr;

  return ( *(matrix + ((nz != 0)&&(nz != 1))*(z-first)*(ny*nx) + (y-first)*nx + (x-first)) );
}
