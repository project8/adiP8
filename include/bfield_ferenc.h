////////////////////////////////////////////////////////////////////////
//
//  Field calculation
//
////////////////////////////////////////////////////////////////////////

void get_bfield_fg(double *posvec, double *bvec);

///////////////////////////////////////////////////////////////////////

void magfield(double z,double r,double *A,double *Bz,double *Br);
// This subroutine computes the axial and radial magnetic
//    field components Bz and Br, and the vectorpotential A
//    in a point with cylindrical coordinates z and r.
// Method: Legendre polynomial expansion around the source point with index k.
//  Nspmax: we need Nspmax>=Nsp
//  nmax: maximal index of the source constants (maximum of n).
//  Nspmax and nmax are given by global #define commands.
//  Important (if magfield is used separated from magsource):
//    the same nmax number used for magsource
//    should also be used for magfield !!!

