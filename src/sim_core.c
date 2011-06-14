using namespace std;
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <iostream>
#include "paramanage.h"
#include "el_pa_tool.h"
#include "math_tool.h"
#include "vector_tool.h"
#include "mag_pa_tool.h"
#include "sim_core.h"
#include "sim_help.h"
#include "sim_scatter.h"

/******************** calc_scatter_event ***********************/
void calc_scatter_event(struct particle_data *particle)
     // this function calculated energy and direction change
     // if particle gets scattered
{
  double vec1[3];
  double axis0[3];
  double rot0[3], rot1[3];
  double trans1[3];
  double final_vec[3];
  double alpha, beta, theta_new;
  double sin_square_theta_new;
  double e_kin_new;
  double e_para_new;
  double e_perp_new;

  // vec1 points in the particle's original direction
  // and has the length of kinetic energy
  vec1[0] = sqrt((*particle).e_para);
  vec1[1] = sqrt((*particle).e_perp);
  vec1[2] = 0.;
  //printf("vec1: %f %f %f\n",vec1[0],vec1[1],vec1[2]);

  // angle for phase rotation of particle vector 
  alpha = get_std_rand_1() * 2. * M_PI;
  //printf("alpha: %f\n",alpha*180./M_PI);

  // axis0 points in the particle's longitudinal direction
  // length is just the rotation angle
  axis0[0] = alpha;
  axis0[1] = 0.;
  axis0[2] = 0.;
  //printf("axis0: %f %f %f\n",axis0[0],axis0[1],axis0[2]);

  // rotates vector vec1 (output vec3) around rotation axis axis0 
  // by angle |axis0| = alpha
  // see Bronstein: section 2.6.5.2.3, page 215,216 
  rotate_vec(vec1, axis0, rot0);
  //printf("rot0: %f %f %f\n",rot0[0],rot0[1],rot0[2]);

  // now find an orthogonal vector (trans1) to rot0
  vector_times_scalar(rot0, 0., trans1);  // make NULL-Vector trans1
  if ((rot0[0] == 0.) || (rot0[1] == 0.) || (rot0[2] == 0.)) {
    if (rot0[0] == 0.) {
      trans1[0] = 1.;
    } else if (rot0[1] == 0.) {
      trans1[1] = 1.;
    } else if (rot0[2] == 0.) {
      trans1[2] = 1.;
    }
  } else {                       // make orthogonal just by (a,b,c) -> (-b,a,c)
    trans1[0] = -rot0[1];
    trans1[1] = rot0[0];
    trans1[2] = rot0[2];
  }

  //  same by use of cross product
  //
  //  cross_prod(rot0,axis0,trans1);
  //
  //  printf("trans1: %f %f %f\n",trans1[0],trans1[1],trans1[2]);

  // orthogonal vector trans1 has to have the right length to
  // reproduce the scattering angle
  // right length is squareroot of E_kin times tan(scattering angle)
  vector_times_scalar(trans1, sqrt((*particle).e_kin) * tan((M_PI / 180.) * scatter_get_angle()) / absvalue(trans1), trans1);
  //  printf("trans1': %f %f %f\n",trans1[0],trans1[1],trans1[2]);

  // to get a new vector with scatter angle included
  // one has to add trans1 to the rotated vector rot0
  vector_sum(trans1, rot0, rot1);
  //printf("rot1': %f %f %f\n",rot1[0],rot1[1],rot1[2]);

  // angle for phase rotation of scattering vector 
  beta = get_std_rand_1() * 2. * M_PI;
  //printf("beta: %f\n",beta*180./M_PI);

  // rot0 is now the vector to rotae around
  // length has to be the rotation angle
  vector_times_scalar(rot0, beta / absvalue(rot0), rot0);


  // rotates vector rot1 (output final_vec) around rotation axis rot0 
  // by angle |rot0| = beta
  // see Bronstein: section 2.6.5.2.3, page 215,216 
  rotate_vec(rot1, rot0, final_vec);
  //printf("rot0': %f %f %f\n",rot0[0],rot0[1],rot0[2]);

  // now get the new particle motion angle theta_new
  // it is the angle between final_vec and axis0
  theta_new = angle_rad(axis0, final_vec);
  //printf("theta_new :  %f\n",theta_new*180./M_PI);

  // now calc the altered kinetic energy
  e_kin_new = (*particle).e_kin - scatter_get_eloss();
  sin_square_theta_new = sin(theta_new) * sin(theta_new);
  e_perp_new = e_kin_new * sin_square_theta_new;
  e_para_new = e_kin_new - e_perp_new;

  //printf("kin energy', eloss:  %f %f\n",e_kin_new,scatter_get_eloss());


  // now put the energy change in particle structure for
  // further calculation in "include_energy_loss"
  // here, energy loss has to be positive to be substracted from
  // particle energy in "include_energy_loss"
  (*particle).e_scatter_para = (*particle).e_para - e_para_new;
  (*particle).e_scatter_perp = (*particle).e_perp - e_perp_new;

  //  printf("angle/energy change:  %f  %f\n",(*particle).e_scatter_para,
  //                                                                                                                                                                                                                                                            (*particle).e_scatter_perp);
  //exit(0);
}


/******************** get_gamma ***********************/
double get_gamma(struct particle_data *particle)
{                               // this function calculates the relativistic gamma factor
  return ((*particle).e_kin + (*particle).mass) / (*particle).mass;
}


/******************** get_gamma_cycl ***********************/
double get_gamma_cycl(struct particle_data *particle)
{                               // this function calculates the relativistic gamma factor
  return ((*particle).e_cycl + (*particle).mass) / (*particle).mass;
}


/******************** get_cyclrad ***********************/
double get_cyclrad(struct particle_data *particle)
{                               // this function calculates the cyclotron radius of particle in cm
  return (sqrt((*particle).e_cycl * (*particle).mass * (get_gamma(particle) + 1.))
          / (3000000. * (*particle).b_value));
}


/******************** get_syncrorad_para ***********************/
double get_syncrorad_para(struct particle_data *particle, double gain)
{                               // this function calculates the perp syncrotron radiation of particle in eV
  double p_syncrad;
  double e_syncrad;
  double electron_charge;
  double lightspeed;
  double delta_step;

  if ((*particle).delta_tof > 0.) {
    electron_charge = Echarge * SI2esE;  // Echarge in CGS system
    lightspeed = Clight * M2CM; // in cm per s for CGS system

    delta_step = absvalue((*particle).step_final);
    // final steplength of previous step

    // equation taken from Jackson Classical EDyn. page 791 (equ. 14.28)
    // important to transform from CGS to SI system !!!!!!!!!!!!!!!!!!
    p_syncrad = (2 * electron_charge * electron_charge * lightspeed * (*particle).delta_e_para * (*particle).delta_e_para)
        / (3 * (*particle).mass * (*particle).mass * delta_step * delta_step);

    e_syncrad = p_syncrad * (*particle).delta_tof * erg2eV;
    // improved numerics, result in eV

    return e_syncrad * gain;
  } else {
    return 0.;
  }
}

/******************** get_syncrorad_para ***********************/
double get_syncrorad_perp(struct particle_data *particle, double gain)
{                               // this function calculates the para syncrotron radiation of particle in eV
  double p_syncrad;
  double e_syncrad;
  double gamma;
  double beta;
  double cycrad;
  double electron_charge;
  double lightspeed;

  gamma = get_gamma_cycl(particle);
  beta = sqrt(1 - 1 / (gamma * gamma));
  if ((*particle).e_cycl > 0.) {
    cycrad = sqrt((*particle).e_cycl * (*particle).mass * (get_gamma(particle) + 1.))
        / (3000000. * (*particle).b_value);
    // in cm for CGS calculation

    electron_charge = Echarge * SI2esE;  // Echarge in CGS system
    lightspeed = Clight * M2CM; // in cm per s for CGS system

    // equation taken from Jackson Classical EDyn. page 791 (equ. 14.28)
    // important to transform from CGS to SI system !!!!!!!!!!!!!!!!!!
    p_syncrad = (2 * electron_charge * electron_charge * lightspeed * beta * beta * beta * beta * gamma * gamma * gamma * gamma) / (3 * cycrad * cycrad);

    e_syncrad = p_syncrad * (*particle).delta_tof * erg2eV;
    // improved numerics, result in eV
  } else {
    e_syncrad = 0.;
  }

  return e_syncrad * gain;
}


/******************** get_perp_energy_loss ***********************/
double get_perp_energy_loss(struct particle_data *particle)
{                               // this function adds up all perpendicular energy losses
  // just the syncrotron radiation energy loss for now 
  double energy_loss;
  double perp_gain = (double) PERP_ENERGY_LOSS;

  if (perp_gain > 0.) {
    energy_loss = get_syncrorad_perp(particle, perp_gain);
  } else {
    energy_loss = 0.;
  }

  (*particle).e_syncro = energy_loss;

  return energy_loss;
}


/******************** get_para_energy_trans ***********************/
double get_para_energy_loss(struct particle_data *particle)
{                               // this function calculates the energy loss due to long. acceleration
  double energy_loss;
  double para_gain = (double) PARA_ENERGY_LOSS;

  if (para_gain > 0.) {
    energy_loss = get_syncrorad_para(particle, para_gain);
  } else {
    energy_loss = 0.;
  }

  (*particle).e_accel = energy_loss;

  return energy_loss;
}


/******************** include_energy_loss ***********************/
void include_energy_loss(struct particle_data *particle, int *e_tag)
{                               // this function gathers and includes all calculated energy losses to current particle

  get_perp_energy_loss(particle);  // syncro loss now in particle structure 
  get_para_energy_loss(particle);  // longitudinal energy loss now in particle structure

  if (PERP_ENERGY_LOSS > 0.) {
    if ((*particle).e_syncro > (*particle).e_perp) {
      (*particle).e_syncro = (*particle).e_perp;
    }
  }

  if (PARA_ENERGY_LOSS > 0.) {
    if ((*particle).e_accel > (*particle).e_para) {
      (*particle).e_accel = (*particle).e_para;
    }
  }

  (*particle).e_perp -= (*particle).e_syncro;
  (*particle).e_cycl -= (*particle).e_syncro;
  (*particle).e_kin -= (*particle).e_syncro;

  (*particle).e_para -= (*particle).e_accel;
  (*particle).e_kin -= (*particle).e_accel;

  (*particle).e_perp -= (*particle).e_scatter_perp;
  (*particle).e_cycl -= (*particle).e_scatter_perp;
  (*particle).e_kin -= (*particle).e_scatter_perp;

  (*particle).e_para -= (*particle).e_scatter_para;
  (*particle).e_kin -= (*particle).e_scatter_para;

  if ((*particle).e_perp < 0.) {
    (*particle).e_perp = 0.;
  }
  if ((*particle).e_cycl < 0.) {
    (*particle).e_cycl = 0.;
  }
  if ((*particle).e_para < (*particle).e_para_min / 2.) {
    (*particle).e_para = (*particle).e_para_min / 2.;
  }

  (*particle).e_accel_sum += (*particle).e_accel;
  (*particle).e_syncro_sum += (*particle).e_syncro;
}


/******************** go_test_step ***********************/
void go_test_step(struct particle_data *particle, struct particle_data *probe)
     // this function adds to the position vector th estep vector
     // in other words it goes one step forward for testing
{
  (*probe) = (*particle);
  vector_sum((*probe).position, (*probe).step_final, (*probe).position);
}


/******************** do_real_step_now ***********************/
void make_step_real(struct particle_data *particle, struct particle_data *probe)
     // this function converts a confirmed test step to a real step
     // this is a real step forward
{
  (*particle) = (*probe);       // store data in probe back in particle
  (*particle).time_of_flight = (*particle).time_of_flight + (*particle).delta_tof * 1000000.;
  // add to the whole tof [in µs] value

  // keep track of cyclotron phase
  // OMEGA0 is in rad/sec and delta_tof is in s 
  (*particle).phase = (*particle).phase + OMEGA0 * (*particle).b_value / get_gamma(particle) * (*particle).delta_tof;  //in radians

  (*particle).omega = OMEGA0 * (*particle).b_value / get_gamma(particle) * 1e-6;  //radians/usec 
  (*particle).cyclrad = get_cyclrad(particle);


  // now stuff for scatter detection

  if (!scatter_check_event()) {
    scatter_add_step(absvalue((*particle).step_final), (*particle).e_kin);
    // add dx/path fraction of current step to total fraction
    (*particle).e_scatter_para = 0.;
    (*particle).e_scatter_perp = 0.;
  } else {
    printf("Particle scattered ");
    scatter_store_sigmas((*particle).e_kin);
    scatter_store_angle_eloss((*particle).e_kin);
    switch (scatter_get_process()) {
      case 0:
        printf("( elastic  ) : 0 ");
        break;
      case 1:
        printf("(ionization) : 1 ");
        break;
      case 2:
        printf("(excitation) : 2 ");
        break;
      default:
        printf("ERROR in scatter_get_eloss: process not set!\n");
    }
    printf("Eloss=%7.4feV tof=%eusec ang=%.2f\n", scatter_get_eloss(), (*particle).time_of_flight, scatter_get_angle());

    sim_help_store_scatter_data(particle, scatter_get_eloss(), scatter_get_angle(), scatter_get_process());

    if (scatter_get_eloss() > (*particle).e_kin) {
      printf("Energy violation! (sim_core.c:make_step_real)\n");
      //                                                                                                                                                                                                                                                             exit(0);
    }

    calc_scatter_event(particle);

    scatter_init_rel_cross_section();  // init scatter detection module
    // for next scattering process
  }
}

/******************** shrink_step ***********************/
void shrink_step(struct particle_data *particle)
{                               // this function halfs the shrinkfactor value each time it was called
  (*particle).shrinkfactor = (*particle).shrinkfactor * 0.5;
}

/******************** unshrink_step ***********************/
void unshrink_step(struct particle_data *particle)
{                               // this function doubles the shrinkfactor value each time it was called
  (*particle).shrinkfactor = (*particle).shrinkfactor * 2.;
}

/******************** turn_direction ***********************/
void turn_direction(struct particle_data *particle)
{                               // this function turns the direction of motion and counts the miror points
  (*particle).v_signum = -(*particle).v_signum;
  (*particle).mirrors++;
  if (((*particle).mirrors) % 100 == 0) {
    printf("reflections=%4.0d   tof=%10.6f   Etot=%10.4f \n", (*particle).mirrors, (*particle).time_of_flight, (*particle).e_kin - (*particle).e_pot);  // put out mirror count
  }
  //  scatter_debug();
}

/******************** test_electrode ***********************/
int test_electrode(struct particle_data *particle, struct particle_data *probe)
{
  double radius;
  double testvec[3];
  double testpos[3];
  double test_here[3];
  double origin[3];
  double gain;
  int e_tag;

  radius = sqrt(pow((*probe).position[1], 2) + pow((*probe).position[2], 2))
      + get_cyclrad(probe);

  vector_times_scalar((*probe).position, 1., testvec);
  vector_times_scalar((*probe).position, 1., origin);
  testvec[0] = 0.;
  origin[1] = 0.;
  origin[2] = 0.;

  if ((testvec[1] == 0.) && (testvec[2] == 0.)) {
    testvec[1] = 1.;
  }

  vector_times_scalar(testvec, MM_PER_UNIT / (2. * absvalue(testvec)), testvec);

  gain = 1.;
  e_tag = 0;

  do {
    vector_times_scalar(testvec, gain, testpos);
    vector_sum(testpos, origin, test_here);
    epot3d(test_here, &e_tag);
    gain++;
  } while ((e_tag == 0) && (absvalue(testpos) <= MAX_RADIUS));

  e_tag = 0;
  if (radius > sqrt(test_here[1] * test_here[1] + test_here[2] * test_here[2])) {
    printf("Particle at radius %5.2fcm splat on electrode at radius %5.2fcm\n", radius, sqrt(test_here[1] * test_here[1] + test_here[2] * test_here[2]));
    e_tag = 1;
  }
  if (radius >= MAX_RADIUS) {
    e_tag = 1;
  }
  return e_tag;
}



/******************** Runge Kutta Midpoint Methode ***********************/
void runge_kutta_midpoint(struct particle_data *particle)
  // this routine calculates the next point on bfield line by using
  // runge kutta midpoint method with second order.
{
  double pos1[3], b1[3], step1[3], step2[3];

  vector_times_scalar((*particle).step_final, 0.5, step1);  // half step value
  vector_sum((*particle).position, step1, pos1);  // go one step (of half value)

  get_bfield(pos1, b1);         // get bfield vector here

  vector_times_scalar(b1, (*particle).step_final[0] / b1[0], step2);  // go one step (whole value)
  // with bfield vector of middle point
  //MLL: vector_times_scalar(b1,absvalue((*particle).step_final)/absvalue(b1),step2); // go one step (whole value)

  vector_times_scalar(step2, 1., (*particle).step_final);  // new step vector in v_step
}


/******************** ExB Drift Calculation ***********************/
void ecrossb_drift(struct particle_data *particle)
{                               // this function calculates the EcrossB drift velocity and uses
  // the runge-kutta-method (2nd order) too
  double evec[3];
  int e_tag = 0;

  // starting point
  efield3d((*particle).position, evec, &e_tag);
  cross_prod(evec, (*particle).b_vec, (*particle).exb_vel);
  if ((*particle).b_value == 0.) {
    printf("Error: B-Value zero in ExB-drift calculation!\n");
  }
  vector_times_scalar((*particle).exb_vel, 1. / ((*particle).b_value * (*particle).b_value), (*particle).exb_vel);
}


/******************** Drift relativity Check ***********************/
int test_drift(struct particle_data *particle)
{ // this function calculates the EcrossB drift velocity and uses
  // the runge-kutta-method (2nd order) too
  int error = 0;

  if (absvalue((*particle).exb_vel) >= Clight / 2.) {
    printf("WARNING: ExB veclocity = %4.2f Clight\n", absvalue((*particle).exb_vel) / Clight);
    error = +1;
  }
  if (absvalue((*particle).dbxb_vel) >= Clight / 2.) {
    printf("WARNING: gradient drift velocity = %4.2f Clight\n", absvalue((*particle).exb_vel) / Clight);
    error = +2;
  }
  if (absvalue((*particle).rxb_vel) >= Clight / 2.) {
    printf("WARNING: curvature drift velocity = %4.2f Clight\n", absvalue((*particle).exb_vel) / Clight);
    error = +4;
  }
  return error;
}


/******************** magnetron drift ***********************/
void mag_drift(struct particle_data *particle)
{
  double pos1[3], pos3[3];      // positions on B field line 
  double bvec1[3], bvec3[3];    // B vector for gradient calculation
  double dbvec[3], db_perp_vec[3];  // vector of db and perpendicular db
  double curv_rad_vec[4];       // curvature radius vector + curvature in 4th component
  double delta = MAG_MM_PER_UNIT;  // delta value for gradient calculation
  double gdrift[3];             // final gradient drift vector
  double rdrift[3];             // curvature drift vector
  double b_unitary_old[3];
  double b_unitary[3];
  int komp;                     // component counter

  // now  computing the gradient and curvature drift

  delta = delta * MM2CM;        // delta in mm must be used in cm 
  if (delta == 0.) {
    printf("Error: MM_PER_UNIT constant not set!\n");
  }
  for (komp = 0; komp <= 2; komp++) {  // this is a simple gradient calculation for Bfield value
    vector_times_scalar((*particle).position, 1., pos1);
    vector_times_scalar((*particle).position, 1., pos3);
    pos1[komp] = pos1[komp] - delta;
    pos3[komp] = pos3[komp] + delta;
    get_bfield(pos1, bvec1);
    get_bfield(pos3, bvec3);
    dbvec[komp] = (absvalue(bvec3) - absvalue(bvec1)) * M2CM / (2. * delta);
  }

  vector_times_scalar((*particle).b_vec_old, 1. / absvalue((*particle).b_vec_old), b_unitary_old);
  vector_times_scalar((*particle).b_vec, 1. / absvalue((*particle).b_vec), b_unitary);
  vector_curvation(b_unitary_old, b_unitary, curv_rad_vec);
  (*particle).curv_rad = 1. / curv_rad_vec[3];
  // find b_vec to b_vec_old transfer in go step function

  gradB_perp((*particle).b_vec, dbvec, db_perp_vec);  // use only perpendicular part of gradB

  cross_prod((*particle).b_vec, db_perp_vec, gdrift);  // calc gradB_perp cross Bvector

  cross_prod(curv_rad_vec, (*particle).b_vec, rdrift);  // calc RxB

  if ((*particle).b_value == 0.) {
    printf("Error: B-Value zero in mag-drift calculation!\n");
  }
  vector_times_scalar(gdrift, -(*particle).charge * (*particle).e_perp / ((*particle).b_value * (*particle).b_value * (*particle).b_value), (*particle).dbxb_vel);  // finishing gradient drift velocity  
//  vector_times_scalar((*particle).dbxb_vel,0.,(*particle).dbxb_vel);
  vector_times_scalar(gdrift, -(*particle).charge * 2. * (*particle).e_para / ((*particle).b_value * (*particle).b_value * (*particle).b_value), (*particle).rxb_vel);  // finishing curvature drift velocity 
  //vector_times_scalar((*particle).rxb_vel,0.,(*particle).rxb_vel);

}


/******************** check_mirror  ***********************/
void check_mirror(struct particle_data *particle, int *e_tag)
{                               // this routine checks if particle is mirrored by potential
  if ((*particle).e_para < (*particle).e_para_min) {
    if ((*particle).shrinkfactor > MIN_SHRINK_FACTOR) {
      *e_tag = *e_tag + 2;
      // when parallel energy below zero
      // there has to be a mirror point, therefor e_tag is set to 2
    } else {
      //                                                                                                                                                                                                                                                             (*particle).e_para = (*particle).e_para_min/2.;
      // don't know why this was used before
      (*particle).shrinkfactor = 1.;
    }
  }
}


/******************** check_total_energy  ***********************/
void check_total_energy(struct particle_data *particle, int *e_tag)
{                               // this routine checks if particle has stopped
  if (parameter.e_min_cooling > 0.) {
    if ((*particle).e_kin - (*particle).e_pot < parameter.e_min_cooling) {
      *e_tag += 100;
    }
  }
}


/******************** calc_drifts ***********************/
void calc_drifts(struct particle_data *particle)
{
  vector_times_scalar((*particle).exb_vel, 0., (*particle).exb_vel);
  vector_times_scalar((*particle).rxb_vel, 0., (*particle).rxb_vel);
  vector_times_scalar((*particle).dbxb_vel, 0., (*particle).dbxb_vel);

  if ((*particle).calc_order > 0) {
    if ((*particle).calc_order == 1) {
      ecrossb_drift(particle);  // calculate ExB drift
    }
    if ((*particle).calc_order == 2) {
      mag_drift(particle);      // calculate curvature and gradient drift
    }
  } else {
    ecrossb_drift(particle);    // calculate ExB drift
    mag_drift(particle);        // calculate curvature and gradient drift
  }
}


/******************** init_particle  ***********************/
void init_particle(struct particle_data *particle)
{                               // this routine is used to init all importent data before simulation begins
  double v_value;
  double v_vec[3];
  int e_tag = 0;

  // calculate start values now
  (*particle).u_start = epot3d((*particle).start_pos, &e_tag);  // potential at start position
  (*particle).e_pot = (*particle).u_start;
  vector_times_scalar((*particle).start_pos, 1., (*particle).position);
  // here starting position equal to position
  (*particle).e_kin = (*particle).e_start;  // e_kin equal to e_start
  (*particle).gamma_start = get_gamma(particle);  // calc rel gamma factor at start
  if ((*particle).starting_theta < 90) {
    (*particle).v_signum = +1.; // start in positive direction
  } else {
    (*particle).v_signum = -1.; // start in negative direction
  }
  (*particle).time_of_flight = 0.;  // reset tof
  (*particle).delta_tof = 0.;   // reset delta_tof
  (*particle).mirrors = 0;      // reset mirror count
  (*particle).shrinkfactor = 1.;  // reset shrinkfactor
  (*particle).e_curv = 0.;
  (*particle).e_grad = 0.;
  (*particle).e_ExB = 0.;
  (*particle).e_syncro = 0.;    // reset syncrotron energy loss
  (*particle).e_syncro_sum = 0.;  // reset accumulated syncrotron energy loss
  (*particle).e_accel = 0.;     // reset accel. energy loss
  (*particle).e_accel_sum = 0.;
  vector_times_scalar((*particle).step_final, 0., (*particle).step_final);  // reset step_final

  v_value = sqrt((*particle).e_start * (*particle).e_start +  // velocity here
                 2. * (*particle).e_start * (*particle).mass) * Clight / ((*particle).e_start + (*particle).mass);

  if (parameter.enable_rel_start_angle == 1) {
    double beta;                // angle between x-direction and B-vector
    double x_unit[3] = { 1., 0., 0. };  // direction of x-axis

    cout << "Change Angle " << (*particle).starting_theta;

    get_bfield((*particle).start_pos, (*particle).b_start);
    // starting bvec vector at start position
    beta = angle((*particle).b_start, x_unit);  // get angle of B-field-vector
    if ((*particle).b_start[1] < 0.) {
      beta = -beta;             // test B-Vector pointing down
    }
    (*particle).starting_phi = 0.;  // no phi needed
    (*particle).starting_theta = beta + (*particle).starting_theta;  // adjust starting_theta
    //if (((*particle).starting_theta > 89.) && ((*particle).starting_theta < 91.)) 
//                                                                                                                                                                                                                                                             (*particle).starting_theta = 89.;

    cout << "deg relative to Bfield to " << (*particle).starting_theta;
    cout << "deg absolute!" << endl;
  }


  spher2kart(v_vec, v_value, (*particle).starting_theta, (*particle).starting_phi);
  // transform sherical vector in kartesian vector
  if (v_vec[0] < 0.) {
    (*particle).v_signum = -1.;
  }                             // turn direction if velo_vec pointing to neg.

  get_bfield((*particle).start_pos, (*particle).b_start);
  // starting bvec vector at start position
  (*particle).b_start_value = absvalue((*particle).b_start);
  // and its value 

  vector_times_scalar((*particle).b_start, 1., (*particle).b_vec_old);
  // store b_start in b_vec used for curv. calc.
  (*particle).b_value_old = (*particle).b_start_value;  // and absolute value

  vector_times_scalar((*particle).b_start, 1., (*particle).b_vec);
  // store b_start in b_vec used for curv. calc.
  (*particle).b_value = (*particle).b_start_value;  // and absolute value

  (*particle).sin2_alpha_start = sin(angle_rad((*particle).b_start, v_vec))
      * sin(angle_rad((*particle).b_start, v_vec));
  // angle between bfield vector and velocity vector

  vector_times_scalar((*particle).b_start, scalar_prod((*particle).b_start, v_vec) / ((*particle).b_start_value * (*particle).b_start_value), (*particle).v_para);
  vector_times_scalar((*particle).v_para, -1., (*particle).v_para);  // turn direction of v_para
  vector_sum(v_vec, (*particle).v_para, (*particle).v_perp_start);
  // to bfield perpendicular part of velocity vector
  vector_times_scalar((*particle).v_para, -1., (*particle).v_para);  // turn back direction of v_para

  (*particle).v_perp_value = absvalue((*particle).v_perp_start);
  (*particle).v_perp_start_value = (*particle).v_perp_value;  // value of perpendicular velocity
  (*particle).v_para_value = absvalue((*particle).v_para);  // value of parallel velocity

  (*particle).e_perp = (*particle).e_start * (*particle).sin2_alpha_start;
  (*particle).e_para = (*particle).e_start - (*particle).e_perp;
  (*particle).omega = OMEGA0 * (*particle).b_value / get_gamma(particle) * 1e-6;  //rad/usec 
}


/******************** calc_energies  ***********************/
void calc_energies(struct particle_data *particle, int *e_tag)
{                               // this routine calculates the energies, perpendicular part and the parallel part
  double gamma;
  double delta_e_perp;
  double e_pot_old;

  e_pot_old = (*particle).e_pot;
  (*particle).e_pot = epot3d((*particle).position, e_tag);

  vector_times_scalar((*particle).b_vec, 1., (*particle).b_vec_old);
  (*particle).b_value_old = (*particle).b_value;  // store old values

  get_bfield((*particle).position, (*particle).b_vec);  // get new values
  (*particle).b_value = absvalue((*particle).b_vec);

  gamma = get_gamma(particle);

  calc_drifts(particle);        // calculate all recommended drift velocities
  // but only of previous step !!!

  (*particle).e_curv = gamma * gamma * (*particle).mass * absvalue((*particle).rxb_vel) * absvalue((*particle).rxb_vel)
      / (Clight * Clight * (gamma + 1));

  (*particle).e_grad = gamma * gamma * (*particle).mass * absvalue((*particle).dbxb_vel) * absvalue((*particle).dbxb_vel)
      / (Clight * Clight * (gamma + 1));

  (*particle).e_ExB = gamma * gamma * (*particle).mass * absvalue((*particle).exb_vel) * absvalue((*particle).exb_vel)
      / (Clight * Clight * (gamma + 1));

  delta_e_perp = (*particle).e_perp * (((*particle).b_value / (*particle).b_value_old) - 1.);

  (*particle).e_perp += delta_e_perp;

  // (*particle).e_perp =
  //    (*particle).e_start*(*particle).sin2_alpha_start*(*particle).b_value*
  //   ((*particle).gamma_start + 1.)/((*particle).b_start_value*(gamma + 1.));
  // old eperp transformation calculation
  // not used since eloss is included

  (*particle).delta_e_para = (*particle).charge * ((*particle).e_pot - e_pot_old) - delta_e_perp;
  // store change of e_para for later calculation of e_accel
  // can also be used here for further calculation

  //  printf(" %20.12f",(*particle).e_para+(*particle).delta_e_para);

  (*particle).e_para += (*particle).delta_e_para;

  //printf("   %20.12f   %20.12f \n",(*particle).e_para,(*particle).delta_e_para);

  (*particle).e_cycl = (*particle).e_perp - (*particle).e_grad;

  (*particle).e_kin = (*particle).e_para + (*particle).e_perp;
  // kinetic energy
}



/******************** calc_velocities  ***********************/
void calc_velocities(struct particle_data *particle)
{
  double v_value;
  double v_trafo;

  if ((*particle).shrinkfactor < 1.) {  // expand shrinkfactor if lower than 1
    unshrink_step(particle);
  } else {
    (*particle).shrinkfactor = 1.;  // set default value for shrinkfactor
  }

  //  (*particle).v_perp_value = (*particle).v_perp_start_value* 
  //                           sqrt((*particle).b_value/(*particle).b_start_value)*
  //                           get_gamma(particle)/(*particle).gamma_start;
  // calculate perpendicular velocity value by squareroot of the fraction of B to Bstart

  (*particle).v_perp_value = sqrt((get_gamma(particle) + 1.) * (*particle).e_perp / ((*particle).mass * get_gamma(particle) * get_gamma(particle))) * Clight;  // new relativistic version
  // calculate perp velocity in each step with respect to current e_perp
  // this depends no more on the start values

  v_value = sqrt((*particle).e_kin * (*particle).e_kin +  // velocity here 
                 2. * (*particle).e_kin * (*particle).mass) * Clight / ((*particle).e_kin + (*particle).mass);
  // calculate value of velocity by using relativistic expression

  v_trafo = (get_gamma(particle) + 1.) * (*particle).e_para / ((*particle).mass * get_gamma(particle) * get_gamma(particle));  // new relativistic version

  vector_times_scalar((*particle).b_vec, (*particle).v_signum * sqrt(v_trafo) * Clight / (*particle).b_value, (*particle).v_para);  // make parallel velocity vector

  (*particle).v_para_value = absvalue((*particle).v_para);  // value of v_para
}

/******************** calc_step ***********************/
void calc_step(struct particle_data *particle)
{
  double all_drifts[3];

  // cut the steplength to value given by max_step_length
  vector_times_scalar((*particle).v_para, (*particle).max_step_length / (*particle).v_para_value, (*particle).step_final);

  do {
    vector_times_scalar((*particle).step_final, (*particle).shrinkfactor, (*particle).step_final);
    // this is the place where the shrinkfactor is finaly used to shorten the step
    runge_kutta_midpoint(particle);  // use now the runge kutta approximation

    if ((*particle).v_para_value <= 0.) {
      printf("ERROR: devide by zero in TOF calculation!\n");
    }
    (*particle).delta_tof = absvalue((*particle).step_final) * 10000. / (*particle).v_para_value;  //in µs
    // tof for actual step in µs
    (*particle).delta_tof = (*particle).delta_tof / (1000000.);  // in sec
    // and transform into sec

    vector_times_scalar(all_drifts, 0., all_drifts);  // earase all_drift vector
    vector_sum(all_drifts, (*particle).exb_vel, all_drifts);  // add to all_drift vector
    vector_sum(all_drifts, (*particle).rxb_vel, all_drifts);  // add to all_drift vector
    vector_sum(all_drifts, (*particle).dbxb_vel, all_drifts);  // add to all_drift vector

    vector_times_scalar(all_drifts, (*particle).delta_tof * M2CM, (*particle).v_drift);
    // transform all_drift velocity vector to v_drift distance vector with factor 100 for cm
    vector_sum((*particle).step_final, (*particle).v_drift, (*particle).step_final);
    // add v_drift to step_final vector

    if (absvalue((*particle).step_final) > (*particle).max_step_length) {
      shrink_step(particle);
    }
    // test max_step_length again
  } while (absvalue((*particle).step_final) > (*particle).max_step_length);
  // repeat this until step_final is shorter then max_step_length
}


/******************** single_track ***********************/
void single_track(FILE * f_fly, struct particle_data *particle)
{
  // Routine for Particle Tracking, runs one single track

  int e_tag = 0;                // error tag
  int loop_counter;             // counter to prevent infinite loops
  //  double tof_backup = 0.; // time backup for pulsing control
  int got_turn = 0;             // stores turning point information of last move
  double dipol_value;           // backup for dipol value, used by pulsing control
  struct particle_data probe;

  //  scatter_debug();
  //  exit(0);
  //  only used for test of scatter module

  printf("\n");

  init_particle(particle);      // init particle data block

  scatter_init_rel_cross_section();  // init scatter detection module

  printf("Calculating energies \n");
  calc_energies(particle, &e_tag);  // energies for the first time
  // include_energy_loss(particle,&e_tag); // gather and include energy losses
  check_mirror(particle, &e_tag);
  check_total_energy(particle, &e_tag);

  printf("Calculating velocities\n");
  calc_velocities(particle);    // velocities for the first time

  (*particle).cyclrad = get_cyclrad(particle);
  sim_help_force_save_data();   // force data to be saved
  sim_help_save_particle_data(f_fly, particle);  // output first point
  sim_help_reset_save_data();
  dipol_value = get_dipole();
  loop_counter = 0;
  do {
    loop_counter++;
    calc_step(particle);        // do the step calculation
    go_test_step(particle, &probe);  // go forward

    e_tag = 0;                  // reset e_tag

    calc_energies(&probe, &e_tag);  // recalculate energies
    include_energy_loss(&probe, &e_tag);  // gather and include energy losses
    check_mirror(&probe, &e_tag);
    check_total_energy(&probe, &e_tag);

    if ((got_turn > 0) && (probe.e_para > probe.e_para_min)) {
      got_turn = 0;
    }
    if (e_tag < 2) {             // in case of no mirror point = energy minimum
                                 // if (e_tag == 1) printf("etag = 1 in single Track\n");
      if (e_tag == 0) {
        if ((probe.e_para <= probe.e_para_min) && (got_turn == 0)) {
          // no collision and energy lower minimum energy
          turn_direction(&probe);  // turn direction of motion
          got_turn = 1;
          loop_counter = 0;
          sim_help_force_save_data();  // force data to be saved
        }

        e_tag = test_electrode(particle, &probe);  // test collision with electrode
        if (test_drift(&probe) != 0) {
          e_tag = e_tag + 10 * test_drift(&probe);
        }

        make_step_real(particle, &probe);  // finally the step

        // output if no error and not first point
        sim_help_save_particle_data(f_fly, particle);
      }
      calc_velocities(particle);  // recalculate velocity
    } else {                     // in case of approaching mirror point shrink step length
      shrink_step(particle);
      e_tag = e_tag - 2;        // remove error code from error detection tag
    }
  } while ((((*particle).position[0] < SPEC_OUT)  // repeat until particle out of range
            && ((*particle).position[0] > SPEC_IN))
           && ((*particle).mirrors <= (*particle).max_mirrors)  // or mirror count reached
           && ((*particle).time_of_flight <= (*particle).max_tof * 1000000.)
           && (e_tag == 0)
           && (loop_counter <= MAX_LOOPS));  // or error happend

  sim_help_force_save_data();   // force data to be saved

  printf("Energy loss by syncrotron radiation:\n");
  printf("              perpendicular = %.10e meV\n", (*particle).e_syncro_sum * 1e3);
  printf("              longitudinal  = %.10e meV\n", (*particle).e_accel_sum * 1e3);

  if (loop_counter > MAX_LOOPS) {
    printf("EXIT by loop count overflow!\n");
  }
  if (((*particle).position[0] > SPEC_OUT) || ((*particle).position[0] < SPEC_IN)) {
    printf("EXIT by Position! --> splat on detector\n");
  }
  if ((*particle).mirrors > (*particle).max_mirrors) {
    printf("EXIT by Mirrors! --> trapped?\n");
  }
  if ((*particle).time_of_flight > (*particle).max_tof * 1000000.) {
    printf("EXIT by TOF! --> trapped?\n");
  }
  if (e_tag >= 90) {
    printf("EXIT by cooling --> particle stopped!\n");
  } else if (e_tag >= 8) {
    printf("EXIT by relativistic Drift --> ADIPARK unable to calculate!\n");
  } else if (e_tag != 0) {
    printf("EXIT by Energy Tag! --> splat on Electrode\n");
  }
}                               // end of routine single_track


/******************** single_track_trap  ***********************/
int single_track_trap(struct particle_data *particle)
{
  // Routine for Particle Tracking, runs one single track
  // special version for trapping calculation

  int loop_counter;             // counter to catch infinite loops
  int e_tag = 0;                // error tag
  int exit_tag = 0;             // exit tag for trapping volume routine 
  int got_turn = 0;             // detects a turn of direction
  struct particle_data probe;

  init_particle(particle);      // init particle data block

  scatter_init_rel_cross_section();  // init scatter detection

  calc_energies(particle, &e_tag);  // energies for the first time

  check_mirror(particle, &e_tag);
  check_total_energy(particle, &e_tag);
  calc_velocities(particle);    // velocities for the first time

  (*particle).cyclrad = get_cyclrad(particle);
  // not mandatory but I don't know if any subroutines use (*particle).cyclrad

  printf("%f %f %f %f %f\n", (*particle).position[0], (*particle).position[1], (*particle).position[2], (*particle).e_start, 180 * asin(sqrt((*particle).sin2_alpha_start)) / M_PI);

  loop_counter = 0;

  do {
    loop_counter++;
    calc_step(particle);        // do the step calculation
    go_test_step(particle, &probe);  // go forward

    e_tag = 0;                  // reset e_tag

    calc_energies(&probe, &e_tag);  // recalculate energies
    include_energy_loss(&probe, &e_tag);
    check_mirror(&probe, &e_tag);
    check_total_energy(&probe, &e_tag);

    if ((got_turn > 0) && (probe.e_para > probe.e_para_min)) {
      got_turn = 0;
    }

    if (e_tag < 2) {             // in case of no mirror point = energy minimum
      if (e_tag == 0) {
        if ((probe.e_para <= probe.e_para_min) && (got_turn == 0)) {
          // no collision and energy lower minimum energy
          turn_direction(&probe);  // turn direction of motion
          got_turn = 1;
          loop_counter = 0;
        }

        e_tag = test_electrode(particle, &probe);
        // test collision with electrode

        if (test_drift(&probe) != 0) {
          e_tag = e_tag + 10 * test_drift(&probe);
        }

        make_step_real(particle, &probe);  // finally the real step
      }
      calc_velocities(particle);  // recalculate velocity
    } else {                     // in case of approaching mirror point shrink step length
      shrink_step(particle);
      loop_counter++;
      //                                                                                                                                                                                                                                                             if (loop_counter <= MAX_LOOPS) e_tag = e_tag - 2;
      e_tag = e_tag - 2;
      // remove error code from error detection tag
    }
  } while ((((*particle).position[0] < SPEC_OUT)
            && ((*particle).position[0] > SPEC_IN))
           // repeat until particle out of range
           && ((*particle).mirrors <= (*particle).max_mirrors)
           // or mirror count reached
           && ((*particle).time_of_flight <= (*particle).max_tof * 1000000.)
           // or max flight time is reached
           && (e_tag == 0)      // or bad energy tag
           && (loop_counter <= MAX_LOOPS));
  // or infinite loop detected

  exit_tag = 0;                 // not trapped
  if (e_tag >= 10) {
    printf("EXIT by relativistic Drift --> ADIPARK unable to calculate!\n");
    exit_tag = 3;
    e_tag = 0;
  }
  if (loop_counter > MAX_LOOPS) {
    printf("possible infinite loop -> force exit! \n");
    exit_tag = 2;
  }
  if (((*particle).position[0] > SPEC_OUT) || ((*particle).position[0] < SPEC_IN)) {
    printf("Out of spectrometer!\n");
  }
  if ((*particle).mirrors > (*particle).max_mirrors) {
    printf("Particle trapped! (mirrors)\n");
    exit_tag = 1;
  }
  if ((*particle).time_of_flight > (*particle).max_tof * 1000000.) {
    printf("Particle trapped! (tof)\n");
    exit_tag = 1;
  }
  if (e_tag != 0) {
    printf("Impact on electrode! (e-tag)\n");
    exit_tag = 0;
  }
  return exit_tag;
}                               // end of routine single_track
