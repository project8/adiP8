double absvalue(double *vector);
// returns the absolute value of a vector

void cross_prod(double *vec1,double *vec2,double *vec3);
// returns in vec3 the cross product: vec3 = vec1 x vec2

void point_curvation(double *x1,double *x2,double *x3,double *radius);
// returns the curvation radius vector in radius[0,1,2]
// and the curvation value in radius[3]
// in point x2 of the curve x1[3]-x2[3]-x3[3]

void vector_curvation(double *p12,double *p23,double *radius);
// returns the curvation radius vector in radius[0,1,2]
// and the curvation value in radius[3]
// in destination of p12 of the curve p12[3]-p23[3]

double scalar_prod(double *vec1, double *vec2);
// returns the scalar product of vector1 and vector2

double angle(double *vec1, double *vec2);
// returns the angle between vec1 and vec2 in degree

void vector_times_scalar(double *vec, double value, double *result);
// product of vector "vec" and scalar "value" returned as vector "result"

void vector_sum(double *vec1, double *vec2, double *sum);
// returns the sum of vectors "vec1" and "vec2" in vector "sum"

void gradB(double delta,double b,double *b_delta, double *grB);
// returns in Vector grB the gradient of B-Field

void gradB_perp(double *Bvec,double *grB,double *grB_perp);
// returns in Vector grB_perp the perpendicular gradient of B-Field 

void spher2kart(double *vec1,double value, double theta,double phi);
// Converts an sperical value with direction by theta and phi in an kertesian vector vec1.

double angle_rad(double *vec1, double *vec2);
// returns the angle between vec1 and vec2 in radian

void rotate_vec(double *vec1,double *vec2, double *vec3);
//rotates vector vec1 (output vec3) in such a way, 
// that vec2 points in z-direction 
// see Bronstein: section 2.6.5.2.3, page 215,216 

void rotate_vec_on_z_axis(double *vec1,double *vec2, double *vec3);
// rotates vector vec1 (output vec3) in such a way, 
// that vec2 points in z-direction 
// see Bronstein: section 2.6.5.2.3, page 215,216 

void spherical_coordinate(double *vec,double* r, double* theta, double* phi);
// gives back from vector vec spherical coordinates 
// => inverse function to spher2kart
