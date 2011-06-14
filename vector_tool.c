using namespace std;
#include <cmath>
#include "math_tool.h"
#include "vector_tool.h"

/* ************************ absvalue ********************** */
double absvalue(double *vector)
{
  double value;

  value = sqrt(vector[0] * vector[0] + vector[1] * vector[1] + vector[2] * vector[2]);
  // returns the absolute value of input vector
  return value;
}

/* ************************ scalar_prod ********************** */
double scalar_prod(double *vec1, double *vec2)
{
  double value;

  value = vec1[0] * vec2[0] + vec1[1] * vec2[1] + vec1[2] * vec2[2];
  // returns the scalar product of two vectors 
  return value;
}

/* ************************ vector_times_scalar ********************** */
void vector_times_scalar(double *vec, double value, double *result)
{
  int index;

  for (index = 0; index <= 2; index++) { // multiply a vector with a scalar
    result[index] = vec[index] * value;
  }
}

/* ************************ vector_sum ********************** */
void vector_sum(double *vec1, double *vec2, double *sum)
{
  int index;

  for (index = 0; index <= 2; index++) { // returns the summe of two vectors
    sum[index] = vec1[index] + vec2[index];
  }
}

/* ************************ angle ********************** */
double angle(double *vec1, double *vec2)
{
  double value;

  value = (180. / M_PI) * acos((scalar_prod(vec1, vec2)) / (absvalue(vec1) * absvalue(vec2)));
  // returns the angle between two vectors in degree
  return value;
}

/* ************************ cross_prod ********************** */
void cross_prod(double *vec1, double *vec2, double *vec3)
{
  vec3[0] = (vec1[1] * vec2[2] - vec1[2] * vec2[1]);  // calculate the vector product
  vec3[1] = (vec1[2] * vec2[0] - vec1[0] * vec2[2]);
  vec3[2] = (vec1[0] * vec2[1] - vec1[1] * vec2[0]);
}

/* ************************ vec_rotate ********************** */
void spher2kart(double *vec1, double value, double theta, double phi)
{
  theta = M_PI * (theta) / 180.0;
  phi = M_PI * (phi) / 180.0;
  vec1[0] = value * cos(theta);
  vec1[1] = value * sin(theta) * cos(phi);
  vec1[2] = value * sin(theta) * sin(phi);
}

/* ************************ point_curvation ********************** */
void point_curvation(double *x1, double *x2, double *x3, double *radius)
{
  double p12[3], p23[3];
  int komp;

  for (komp = 0; komp <= 2; komp++) {
    p12[komp] = (x2[komp] - x1[komp]);  // calculate point to point vectors
    p23[komp] = (x3[komp] - x2[komp]);
  }

  vector_curvation(p12, p23, radius);
}

/* ************************ vector_curvation ********************** */
void vector_curvation(double *p12, double *p23, double *radius)
{
  double k0[3], p13[3], n[3], n1[3], n2[3];
  double a;
  int komp;

  for (komp = 0; komp <= 2; komp++) {
    p12[komp] = 0.5 * p12[komp];
    p23[komp] = 0.5 * p23[komp];
    // vectors to the centre 
    p13[komp] = p12[komp] + p23[komp];
    // vector from origin of p12 to destination of p23
  }

  cross_prod(p12, p23, n);      // calculate normal vector

  cross_prod(n, p12, n1);       // calculate perpendicular vector of
  // point to point vector 1 in normal plane

  cross_prod(n, p23, n2);       // calculate perpendicular vector of
  // point to point vector 2 in normal plane

  a = 0;
  if ((n1[0] * n2[1] - n1[1] * n2[0]) != 0) {
    a = (p13[0] * n2[1] - p13[1] * n2[0]) / (n1[0] * n2[1] - n1[1] * n2[0]);
  } else {
    if ((n1[0] * n2[2] - n1[2] * n2[0]) != 0) {
      a = (p13[0] * n2[2] - p13[2] * n2[0]) / (n1[0] * n2[2] - n1[2] * n2[0]);
    } else {
      if ((n1[1] * n2[2] - n1[2] * n2[1]) != 0) {
        a = (p13[1] * n2[2] - p13[2] * n2[1]) / (n1[1] * n2[2] - n1[2] * n2[1]);
      }
    }
  }
  // calculates the multiplication faktor of perpendicular vector 1 to
  // the crosspoint of both perpendicular vectors

  for (komp = 0; komp <= 2; komp++) {
    k0[komp] = p12[komp] + a * n1[komp];  // calculate origin of radius vector
    radius[komp] = 2. * p12[komp] - k0[komp];  // calculate radius vector
  }

  radius[3] = 1. / absvalue(radius);  // store curvation value in fourth component
}


void gradB(double delta, double b, double *b_delta, double *grB)
{
  grB[1] = (b_delta[1] - b) / delta;
}



/* ************************ gradB_perp ********************** */
void gradB_perp(double *Bvec, double *grB, double *grB_perp)
{
  double grB_para[3];

  vector_times_scalar(Bvec, -scalar_prod(Bvec, grB) / (absvalue(Bvec) * absvalue(Bvec)), grB_para);
  vector_sum(grB_para, grB, grB_perp);
}

/* ************************ angle_rad ********************** */
double angle_rad(double *vec1, double *vec2)
{
  double value;

  value = acos((scalar_prod(vec1, vec2)) / (absvalue(vec1) * absvalue(vec2)));

  return value;
}


/* ************************ rotate_vec ********************** */
void rotate_vec(double *vec1, double *omega, double *vec3)
// rotates vector vec1 (output vec3) around rotation axis omega 
// by angle |omega| 
// see Bronstein: section 2.6.5.2.3, page 215,216 
{
  double theta, stheta, ctheta, one_ctheta, alpha, beta, gamma;

  theta = absvalue(omega);
  ctheta = cos(theta);
  one_ctheta = 1.0 - cos(theta);
  stheta = sin(theta);
  alpha = omega[0] / absvalue(omega);
  beta = omega[1] / absvalue(omega);
  gamma = omega[2] / absvalue(omega);
  //  printf("alpha %f, beta %f gamma %f theta %f\n", alpha, beta, gamma, theta);
  vec3[0] = vec1[0] * (ctheta + alpha * alpha * one_ctheta) + vec1[1] * (gamma * stheta + alpha * beta * one_ctheta) + vec1[2] * (-beta * stheta + alpha * gamma * one_ctheta);
  vec3[1] = vec1[0] * (-gamma * stheta + alpha * beta * one_ctheta) + vec1[1] * (ctheta + beta * beta * one_ctheta) + vec1[2] * (alpha * stheta + beta * gamma * one_ctheta);
  vec3[2] = vec1[0] * (beta * stheta + gamma * alpha * one_ctheta) + vec1[1] * (-alpha * stheta + gamma * beta * one_ctheta) + vec1[2] * (ctheta + gamma * gamma * one_ctheta);
  return;
}


/* ************************ rotate_vec_on_z_axis ********************** */
void rotate_vec_on_z_axis(double *vec1, double *vec2, double *vec3)
// rotates vector vec1 (output vec3) in such a way, 
// that vec2 points in z-direction 
// see Bronstein: section 2.6.5.2.3, page 215,216 
{
  double vecz[3], rot_axis[3];
  double theta, stheta, ctheta, one_ctheta, alpha, beta, gamma, abs_rot_axis;
  vecz[0] = 0.0;
  vecz[1] = 0.0;
  vecz[2] = 1.0;
  theta = angle_rad(vec2, vecz);
  ctheta = cos(theta);
  one_ctheta = 1.0 - cos(theta);
  stheta = sin(theta);
  cross_prod(vecz, vec2, rot_axis);
  abs_rot_axis = absvalue(rot_axis);
  //  printf("rot_axis %f %f %f\n",rot_axis[0],rot_axis[1], rot_axis[2]);
  if (abs_rot_axis == 0.0) {
    vec3[0] = vec1[0];
    vec3[1] = vec1[1];
    vec3[2] = vec1[2];
  } else {
    alpha = rot_axis[0] / abs_rot_axis;
    beta = rot_axis[1] / abs_rot_axis;
    gamma = rot_axis[2] / abs_rot_axis;
    //printf("alpha %f, beta %f gamma %f theta %f\n", alpha, beta, gamma, theta);
    vec3[0] = vec2[0] * (ctheta + alpha * alpha * one_ctheta) + vec2[1] * (gamma * stheta + alpha * beta * one_ctheta) + vec2[2] * (-beta * stheta + alpha * gamma * one_ctheta);
    vec3[1] = vec2[0] * (-gamma * stheta + alpha * beta * one_ctheta) + vec2[1] * (ctheta + beta * beta * one_ctheta) + vec2[2] * (alpha * stheta + beta * gamma * one_ctheta);
    vec3[2] = vec2[0] * (beta * stheta + gamma * alpha * one_ctheta) + vec2[1] * (-alpha * stheta + gamma * beta * one_ctheta) + vec2[2] * (ctheta + gamma * gamma * one_ctheta);
    //    printf("rotated vec2: %f %f %f\n",vec3[0],vec3[1],vec3[2]);
    vec3[0] = vec1[0] * (ctheta + alpha * alpha * one_ctheta) + vec1[1] * (gamma * stheta + alpha * beta * one_ctheta) + vec1[2] * (-beta * stheta + alpha * gamma * one_ctheta);
    vec3[1] = vec1[0] * (-gamma * stheta + alpha * beta * one_ctheta) + vec1[1] * (ctheta + beta * beta * one_ctheta) + vec1[2] * (alpha * stheta + beta * gamma * one_ctheta);
    vec3[2] = vec1[0] * (beta * stheta + gamma * alpha * one_ctheta) + vec1[1] * (-alpha * stheta + gamma * beta * one_ctheta) + vec1[2] * (ctheta + gamma * gamma * one_ctheta);
  }
  return;
}

/* ************************ spherical_coordinate  ********************** */
void spherical_coordinate(double *vec, double *r, double *theta, double *phi)
// gives back from vector vec spherical coordinates 
{

  double vecz[3];

  vecz[0] = 0.0;
  vecz[1] = 0.0;
  vecz[2] = 1.0;
  //  printf("vec %f %f %f theta %f\n",vec[0],vec[1],vec[2],angle(vec,vecz));  

  *r = absvalue(vec);
  *theta = angle_rad(vec, vecz);
  if (vec[1] == 0.0) {
    *phi = signum(vec[0]) * M_PI / 2.0;
  } else {
    if (vec[1] > 0) {
      *phi = atan(vec[0] / vec[1]);
    } else {
      *phi = M_PI + atan(vec[0] / vec[1]);
    }
  }
  return;
}
