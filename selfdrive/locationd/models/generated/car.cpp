#include "car.h"

namespace {
#define DIM 9
#define EDIM 9
#define MEDIM 9
typedef void (*Hfun)(double *, double *, double *);

double mass;

void set_mass(double x){ mass = x;}

double rotational_inertia;

void set_rotational_inertia(double x){ rotational_inertia = x;}

double center_to_front;

void set_center_to_front(double x){ center_to_front = x;}

double center_to_rear;

void set_center_to_rear(double x){ center_to_rear = x;}

double stiffness_front;

void set_stiffness_front(double x){ stiffness_front = x;}

double stiffness_rear;

void set_stiffness_rear(double x){ stiffness_rear = x;}
const static double MAHA_THRESH_25 = 3.8414588206941227;
const static double MAHA_THRESH_24 = 5.991464547107981;
const static double MAHA_THRESH_30 = 3.8414588206941227;
const static double MAHA_THRESH_26 = 3.8414588206941227;
const static double MAHA_THRESH_27 = 3.8414588206941227;
const static double MAHA_THRESH_29 = 3.8414588206941227;
const static double MAHA_THRESH_28 = 3.8414588206941227;
const static double MAHA_THRESH_31 = 3.8414588206941227;

/******************************************************************************
 *                       Code generated with sympy 1.9                        *
 *                                                                            *
 *              See http://www.sympy.org/ for more information.               *
 *                                                                            *
 *                         This file is part of 'ekf'                         *
 ******************************************************************************/
void err_fun(double *nom_x, double *delta_x, double *out_1914130395499138379) {
   out_1914130395499138379[0] = delta_x[0] + nom_x[0];
   out_1914130395499138379[1] = delta_x[1] + nom_x[1];
   out_1914130395499138379[2] = delta_x[2] + nom_x[2];
   out_1914130395499138379[3] = delta_x[3] + nom_x[3];
   out_1914130395499138379[4] = delta_x[4] + nom_x[4];
   out_1914130395499138379[5] = delta_x[5] + nom_x[5];
   out_1914130395499138379[6] = delta_x[6] + nom_x[6];
   out_1914130395499138379[7] = delta_x[7] + nom_x[7];
   out_1914130395499138379[8] = delta_x[8] + nom_x[8];
}
void inv_err_fun(double *nom_x, double *true_x, double *out_2456965719472217437) {
   out_2456965719472217437[0] = -nom_x[0] + true_x[0];
   out_2456965719472217437[1] = -nom_x[1] + true_x[1];
   out_2456965719472217437[2] = -nom_x[2] + true_x[2];
   out_2456965719472217437[3] = -nom_x[3] + true_x[3];
   out_2456965719472217437[4] = -nom_x[4] + true_x[4];
   out_2456965719472217437[5] = -nom_x[5] + true_x[5];
   out_2456965719472217437[6] = -nom_x[6] + true_x[6];
   out_2456965719472217437[7] = -nom_x[7] + true_x[7];
   out_2456965719472217437[8] = -nom_x[8] + true_x[8];
}
void H_mod_fun(double *state, double *out_4779992119913662016) {
   out_4779992119913662016[0] = 1.0;
   out_4779992119913662016[1] = 0;
   out_4779992119913662016[2] = 0;
   out_4779992119913662016[3] = 0;
   out_4779992119913662016[4] = 0;
   out_4779992119913662016[5] = 0;
   out_4779992119913662016[6] = 0;
   out_4779992119913662016[7] = 0;
   out_4779992119913662016[8] = 0;
   out_4779992119913662016[9] = 0;
   out_4779992119913662016[10] = 1.0;
   out_4779992119913662016[11] = 0;
   out_4779992119913662016[12] = 0;
   out_4779992119913662016[13] = 0;
   out_4779992119913662016[14] = 0;
   out_4779992119913662016[15] = 0;
   out_4779992119913662016[16] = 0;
   out_4779992119913662016[17] = 0;
   out_4779992119913662016[18] = 0;
   out_4779992119913662016[19] = 0;
   out_4779992119913662016[20] = 1.0;
   out_4779992119913662016[21] = 0;
   out_4779992119913662016[22] = 0;
   out_4779992119913662016[23] = 0;
   out_4779992119913662016[24] = 0;
   out_4779992119913662016[25] = 0;
   out_4779992119913662016[26] = 0;
   out_4779992119913662016[27] = 0;
   out_4779992119913662016[28] = 0;
   out_4779992119913662016[29] = 0;
   out_4779992119913662016[30] = 1.0;
   out_4779992119913662016[31] = 0;
   out_4779992119913662016[32] = 0;
   out_4779992119913662016[33] = 0;
   out_4779992119913662016[34] = 0;
   out_4779992119913662016[35] = 0;
   out_4779992119913662016[36] = 0;
   out_4779992119913662016[37] = 0;
   out_4779992119913662016[38] = 0;
   out_4779992119913662016[39] = 0;
   out_4779992119913662016[40] = 1.0;
   out_4779992119913662016[41] = 0;
   out_4779992119913662016[42] = 0;
   out_4779992119913662016[43] = 0;
   out_4779992119913662016[44] = 0;
   out_4779992119913662016[45] = 0;
   out_4779992119913662016[46] = 0;
   out_4779992119913662016[47] = 0;
   out_4779992119913662016[48] = 0;
   out_4779992119913662016[49] = 0;
   out_4779992119913662016[50] = 1.0;
   out_4779992119913662016[51] = 0;
   out_4779992119913662016[52] = 0;
   out_4779992119913662016[53] = 0;
   out_4779992119913662016[54] = 0;
   out_4779992119913662016[55] = 0;
   out_4779992119913662016[56] = 0;
   out_4779992119913662016[57] = 0;
   out_4779992119913662016[58] = 0;
   out_4779992119913662016[59] = 0;
   out_4779992119913662016[60] = 1.0;
   out_4779992119913662016[61] = 0;
   out_4779992119913662016[62] = 0;
   out_4779992119913662016[63] = 0;
   out_4779992119913662016[64] = 0;
   out_4779992119913662016[65] = 0;
   out_4779992119913662016[66] = 0;
   out_4779992119913662016[67] = 0;
   out_4779992119913662016[68] = 0;
   out_4779992119913662016[69] = 0;
   out_4779992119913662016[70] = 1.0;
   out_4779992119913662016[71] = 0;
   out_4779992119913662016[72] = 0;
   out_4779992119913662016[73] = 0;
   out_4779992119913662016[74] = 0;
   out_4779992119913662016[75] = 0;
   out_4779992119913662016[76] = 0;
   out_4779992119913662016[77] = 0;
   out_4779992119913662016[78] = 0;
   out_4779992119913662016[79] = 0;
   out_4779992119913662016[80] = 1.0;
}
void f_fun(double *state, double dt, double *out_2238700801222790874) {
   out_2238700801222790874[0] = state[0];
   out_2238700801222790874[1] = state[1];
   out_2238700801222790874[2] = state[2];
   out_2238700801222790874[3] = state[3];
   out_2238700801222790874[4] = state[4];
   out_2238700801222790874[5] = dt*((-state[4] + (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(mass*state[4]))*state[6] - 9.8000000000000007*state[8] + stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(mass*state[1]) + (-stiffness_front*state[0] - stiffness_rear*state[0])*state[5]/(mass*state[4])) + state[5];
   out_2238700801222790874[6] = dt*(center_to_front*stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(rotational_inertia*state[1]) + (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])*state[5]/(rotational_inertia*state[4]) + (-pow(center_to_front, 2)*stiffness_front*state[0] - pow(center_to_rear, 2)*stiffness_rear*state[0])*state[6]/(rotational_inertia*state[4])) + state[6];
   out_2238700801222790874[7] = state[7];
   out_2238700801222790874[8] = state[8];
}
void F_fun(double *state, double dt, double *out_3188951629714217883) {
   out_3188951629714217883[0] = 1;
   out_3188951629714217883[1] = 0;
   out_3188951629714217883[2] = 0;
   out_3188951629714217883[3] = 0;
   out_3188951629714217883[4] = 0;
   out_3188951629714217883[5] = 0;
   out_3188951629714217883[6] = 0;
   out_3188951629714217883[7] = 0;
   out_3188951629714217883[8] = 0;
   out_3188951629714217883[9] = 0;
   out_3188951629714217883[10] = 1;
   out_3188951629714217883[11] = 0;
   out_3188951629714217883[12] = 0;
   out_3188951629714217883[13] = 0;
   out_3188951629714217883[14] = 0;
   out_3188951629714217883[15] = 0;
   out_3188951629714217883[16] = 0;
   out_3188951629714217883[17] = 0;
   out_3188951629714217883[18] = 0;
   out_3188951629714217883[19] = 0;
   out_3188951629714217883[20] = 1;
   out_3188951629714217883[21] = 0;
   out_3188951629714217883[22] = 0;
   out_3188951629714217883[23] = 0;
   out_3188951629714217883[24] = 0;
   out_3188951629714217883[25] = 0;
   out_3188951629714217883[26] = 0;
   out_3188951629714217883[27] = 0;
   out_3188951629714217883[28] = 0;
   out_3188951629714217883[29] = 0;
   out_3188951629714217883[30] = 1;
   out_3188951629714217883[31] = 0;
   out_3188951629714217883[32] = 0;
   out_3188951629714217883[33] = 0;
   out_3188951629714217883[34] = 0;
   out_3188951629714217883[35] = 0;
   out_3188951629714217883[36] = 0;
   out_3188951629714217883[37] = 0;
   out_3188951629714217883[38] = 0;
   out_3188951629714217883[39] = 0;
   out_3188951629714217883[40] = 1;
   out_3188951629714217883[41] = 0;
   out_3188951629714217883[42] = 0;
   out_3188951629714217883[43] = 0;
   out_3188951629714217883[44] = 0;
   out_3188951629714217883[45] = dt*(stiffness_front*(-state[2] - state[3] + state[7])/(mass*state[1]) + (-stiffness_front - stiffness_rear)*state[5]/(mass*state[4]) + (-center_to_front*stiffness_front + center_to_rear*stiffness_rear)*state[6]/(mass*state[4]));
   out_3188951629714217883[46] = -dt*stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(mass*pow(state[1], 2));
   out_3188951629714217883[47] = -dt*stiffness_front*state[0]/(mass*state[1]);
   out_3188951629714217883[48] = -dt*stiffness_front*state[0]/(mass*state[1]);
   out_3188951629714217883[49] = dt*((-1 - (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(mass*pow(state[4], 2)))*state[6] - (-stiffness_front*state[0] - stiffness_rear*state[0])*state[5]/(mass*pow(state[4], 2)));
   out_3188951629714217883[50] = dt*(-stiffness_front*state[0] - stiffness_rear*state[0])/(mass*state[4]) + 1;
   out_3188951629714217883[51] = dt*(-state[4] + (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(mass*state[4]));
   out_3188951629714217883[52] = dt*stiffness_front*state[0]/(mass*state[1]);
   out_3188951629714217883[53] = -9.8000000000000007*dt;
   out_3188951629714217883[54] = dt*(center_to_front*stiffness_front*(-state[2] - state[3] + state[7])/(rotational_inertia*state[1]) + (-center_to_front*stiffness_front + center_to_rear*stiffness_rear)*state[5]/(rotational_inertia*state[4]) + (-pow(center_to_front, 2)*stiffness_front - pow(center_to_rear, 2)*stiffness_rear)*state[6]/(rotational_inertia*state[4]));
   out_3188951629714217883[55] = -center_to_front*dt*stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(rotational_inertia*pow(state[1], 2));
   out_3188951629714217883[56] = -center_to_front*dt*stiffness_front*state[0]/(rotational_inertia*state[1]);
   out_3188951629714217883[57] = -center_to_front*dt*stiffness_front*state[0]/(rotational_inertia*state[1]);
   out_3188951629714217883[58] = dt*(-(-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])*state[5]/(rotational_inertia*pow(state[4], 2)) - (-pow(center_to_front, 2)*stiffness_front*state[0] - pow(center_to_rear, 2)*stiffness_rear*state[0])*state[6]/(rotational_inertia*pow(state[4], 2)));
   out_3188951629714217883[59] = dt*(-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(rotational_inertia*state[4]);
   out_3188951629714217883[60] = dt*(-pow(center_to_front, 2)*stiffness_front*state[0] - pow(center_to_rear, 2)*stiffness_rear*state[0])/(rotational_inertia*state[4]) + 1;
   out_3188951629714217883[61] = center_to_front*dt*stiffness_front*state[0]/(rotational_inertia*state[1]);
   out_3188951629714217883[62] = 0;
   out_3188951629714217883[63] = 0;
   out_3188951629714217883[64] = 0;
   out_3188951629714217883[65] = 0;
   out_3188951629714217883[66] = 0;
   out_3188951629714217883[67] = 0;
   out_3188951629714217883[68] = 0;
   out_3188951629714217883[69] = 0;
   out_3188951629714217883[70] = 1;
   out_3188951629714217883[71] = 0;
   out_3188951629714217883[72] = 0;
   out_3188951629714217883[73] = 0;
   out_3188951629714217883[74] = 0;
   out_3188951629714217883[75] = 0;
   out_3188951629714217883[76] = 0;
   out_3188951629714217883[77] = 0;
   out_3188951629714217883[78] = 0;
   out_3188951629714217883[79] = 0;
   out_3188951629714217883[80] = 1;
}
void h_25(double *state, double *unused, double *out_2955682509261511320) {
   out_2955682509261511320[0] = state[6];
}
void H_25(double *state, double *unused, double *out_7551971679374724871) {
   out_7551971679374724871[0] = 0;
   out_7551971679374724871[1] = 0;
   out_7551971679374724871[2] = 0;
   out_7551971679374724871[3] = 0;
   out_7551971679374724871[4] = 0;
   out_7551971679374724871[5] = 0;
   out_7551971679374724871[6] = 1;
   out_7551971679374724871[7] = 0;
   out_7551971679374724871[8] = 0;
}
void h_24(double *state, double *unused, double *out_4912051971117535801) {
   out_4912051971117535801[0] = state[4];
   out_4912051971117535801[1] = state[5];
}
void H_24(double *state, double *unused, double *out_976399872783206770) {
   out_976399872783206770[0] = 0;
   out_976399872783206770[1] = 0;
   out_976399872783206770[2] = 0;
   out_976399872783206770[3] = 0;
   out_976399872783206770[4] = 1;
   out_976399872783206770[5] = 0;
   out_976399872783206770[6] = 0;
   out_976399872783206770[7] = 0;
   out_976399872783206770[8] = 0;
   out_976399872783206770[9] = 0;
   out_976399872783206770[10] = 0;
   out_976399872783206770[11] = 0;
   out_976399872783206770[12] = 0;
   out_976399872783206770[13] = 0;
   out_976399872783206770[14] = 1;
   out_976399872783206770[15] = 0;
   out_976399872783206770[16] = 0;
   out_976399872783206770[17] = 0;
}
void h_30(double *state, double *unused, double *out_5047370505310505677) {
   out_5047370505310505677[0] = state[4];
}
void H_30(double *state, double *unused, double *out_635281337883108116) {
   out_635281337883108116[0] = 0;
   out_635281337883108116[1] = 0;
   out_635281337883108116[2] = 0;
   out_635281337883108116[3] = 0;
   out_635281337883108116[4] = 1;
   out_635281337883108116[5] = 0;
   out_635281337883108116[6] = 0;
   out_635281337883108116[7] = 0;
   out_635281337883108116[8] = 0;
}
void h_26(double *state, double *unused, double *out_1369527850446061802) {
   out_1369527850446061802[0] = state[7];
}
void H_26(double *state, double *unused, double *out_7153269075460770521) {
   out_7153269075460770521[0] = 0;
   out_7153269075460770521[1] = 0;
   out_7153269075460770521[2] = 0;
   out_7153269075460770521[3] = 0;
   out_7153269075460770521[4] = 0;
   out_7153269075460770521[5] = 0;
   out_7153269075460770521[6] = 0;
   out_7153269075460770521[7] = 1;
   out_7153269075460770521[8] = 0;
}
void h_27(double *state, double *unused, double *out_3899192362163173755) {
   out_3899192362163173755[0] = state[3];
}
void H_27(double *state, double *unused, double *out_2810044649683533027) {
   out_2810044649683533027[0] = 0;
   out_2810044649683533027[1] = 0;
   out_2810044649683533027[2] = 0;
   out_2810044649683533027[3] = 1;
   out_2810044649683533027[4] = 0;
   out_2810044649683533027[5] = 0;
   out_2810044649683533027[6] = 0;
   out_2810044649683533027[7] = 0;
   out_2810044649683533027[8] = 0;
}
void h_29(double *state, double *unused, double *out_595359365770451158) {
   out_595359365770451158[0] = state[1];
}
void H_29(double *state, double *unused, double *out_4523407376553084060) {
   out_4523407376553084060[0] = 0;
   out_4523407376553084060[1] = 1;
   out_4523407376553084060[2] = 0;
   out_4523407376553084060[3] = 0;
   out_4523407376553084060[4] = 0;
   out_4523407376553084060[5] = 0;
   out_4523407376553084060[6] = 0;
   out_4523407376553084060[7] = 0;
   out_4523407376553084060[8] = 0;
}
void h_28(double *state, double *unused, double *out_5177906949419811836) {
   out_5177906949419811836[0] = state[0];
}
void H_28(double *state, double *unused, double *out_8840937680086936982) {
   out_8840937680086936982[0] = 1;
   out_8840937680086936982[1] = 0;
   out_8840937680086936982[2] = 0;
   out_8840937680086936982[3] = 0;
   out_8840937680086936982[4] = 0;
   out_8840937680086936982[5] = 0;
   out_8840937680086936982[6] = 0;
   out_8840937680086936982[7] = 0;
   out_8840937680086936982[8] = 0;
}
void h_31(double *state, double *unused, double *out_3815152717088839616) {
   out_3815152717088839616[0] = state[8];
}
void H_31(double *state, double *unused, double *out_7521325717497764443) {
   out_7521325717497764443[0] = 0;
   out_7521325717497764443[1] = 0;
   out_7521325717497764443[2] = 0;
   out_7521325717497764443[3] = 0;
   out_7521325717497764443[4] = 0;
   out_7521325717497764443[5] = 0;
   out_7521325717497764443[6] = 0;
   out_7521325717497764443[7] = 0;
   out_7521325717497764443[8] = 1;
}
#include <eigen3/Eigen/Dense>
#include <iostream>

typedef Eigen::Matrix<double, DIM, DIM, Eigen::RowMajor> DDM;
typedef Eigen::Matrix<double, EDIM, EDIM, Eigen::RowMajor> EEM;
typedef Eigen::Matrix<double, DIM, EDIM, Eigen::RowMajor> DEM;

void predict(double *in_x, double *in_P, double *in_Q, double dt) {
  typedef Eigen::Matrix<double, MEDIM, MEDIM, Eigen::RowMajor> RRM;

  double nx[DIM] = {0};
  double in_F[EDIM*EDIM] = {0};

  // functions from sympy
  f_fun(in_x, dt, nx);
  F_fun(in_x, dt, in_F);


  EEM F(in_F);
  EEM P(in_P);
  EEM Q(in_Q);

  RRM F_main = F.topLeftCorner(MEDIM, MEDIM);
  P.topLeftCorner(MEDIM, MEDIM) = (F_main * P.topLeftCorner(MEDIM, MEDIM)) * F_main.transpose();
  P.topRightCorner(MEDIM, EDIM - MEDIM) = F_main * P.topRightCorner(MEDIM, EDIM - MEDIM);
  P.bottomLeftCorner(EDIM - MEDIM, MEDIM) = P.bottomLeftCorner(EDIM - MEDIM, MEDIM) * F_main.transpose();

  P = P + dt*Q;

  // copy out state
  memcpy(in_x, nx, DIM * sizeof(double));
  memcpy(in_P, P.data(), EDIM * EDIM * sizeof(double));
}

// note: extra_args dim only correct when null space projecting
// otherwise 1
template <int ZDIM, int EADIM, bool MAHA_TEST>
void update(double *in_x, double *in_P, Hfun h_fun, Hfun H_fun, Hfun Hea_fun, double *in_z, double *in_R, double *in_ea, double MAHA_THRESHOLD) {
  typedef Eigen::Matrix<double, ZDIM, ZDIM, Eigen::RowMajor> ZZM;
  typedef Eigen::Matrix<double, ZDIM, DIM, Eigen::RowMajor> ZDM;
  typedef Eigen::Matrix<double, Eigen::Dynamic, EDIM, Eigen::RowMajor> XEM;
  //typedef Eigen::Matrix<double, EDIM, ZDIM, Eigen::RowMajor> EZM;
  typedef Eigen::Matrix<double, Eigen::Dynamic, 1> X1M;
  typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> XXM;

  double in_hx[ZDIM] = {0};
  double in_H[ZDIM * DIM] = {0};
  double in_H_mod[EDIM * DIM] = {0};
  double delta_x[EDIM] = {0};
  double x_new[DIM] = {0};


  // state x, P
  Eigen::Matrix<double, ZDIM, 1> z(in_z);
  EEM P(in_P);
  ZZM pre_R(in_R);

  // functions from sympy
  h_fun(in_x, in_ea, in_hx);
  H_fun(in_x, in_ea, in_H);
  ZDM pre_H(in_H);

  // get y (y = z - hx)
  Eigen::Matrix<double, ZDIM, 1> pre_y(in_hx); pre_y = z - pre_y;
  X1M y; XXM H; XXM R;
  if (Hea_fun){
    typedef Eigen::Matrix<double, ZDIM, EADIM, Eigen::RowMajor> ZAM;
    double in_Hea[ZDIM * EADIM] = {0};
    Hea_fun(in_x, in_ea, in_Hea);
    ZAM Hea(in_Hea);
    XXM A = Hea.transpose().fullPivLu().kernel();


    y = A.transpose() * pre_y;
    H = A.transpose() * pre_H;
    R = A.transpose() * pre_R * A;
  } else {
    y = pre_y;
    H = pre_H;
    R = pre_R;
  }
  // get modified H
  H_mod_fun(in_x, in_H_mod);
  DEM H_mod(in_H_mod);
  XEM H_err = H * H_mod;

  // Do mahalobis distance test
  if (MAHA_TEST){
    XXM a = (H_err * P * H_err.transpose() + R).inverse();
    double maha_dist = y.transpose() * a * y;
    if (maha_dist > MAHA_THRESHOLD){
      R = 1.0e16 * R;
    }
  }

  // Outlier resilient weighting
  double weight = 1;//(1.5)/(1 + y.squaredNorm()/R.sum());

  // kalman gains and I_KH
  XXM S = ((H_err * P) * H_err.transpose()) + R/weight;
  XEM KT = S.fullPivLu().solve(H_err * P.transpose());
  //EZM K = KT.transpose(); TODO: WHY DOES THIS NOT COMPILE?
  //EZM K = S.fullPivLu().solve(H_err * P.transpose()).transpose();
  //std::cout << "Here is the matrix rot:\n" << K << std::endl;
  EEM I_KH = Eigen::Matrix<double, EDIM, EDIM>::Identity() - (KT.transpose() * H_err);

  // update state by injecting dx
  Eigen::Matrix<double, EDIM, 1> dx(delta_x);
  dx  = (KT.transpose() * y);
  memcpy(delta_x, dx.data(), EDIM * sizeof(double));
  err_fun(in_x, delta_x, x_new);
  Eigen::Matrix<double, DIM, 1> x(x_new);

  // update cov
  P = ((I_KH * P) * I_KH.transpose()) + ((KT.transpose() * R) * KT);

  // copy out state
  memcpy(in_x, x.data(), DIM * sizeof(double));
  memcpy(in_P, P.data(), EDIM * EDIM * sizeof(double));
  memcpy(in_z, y.data(), y.rows() * sizeof(double));
}




}
extern "C" {

void car_update_25(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
  update<1, 3, 0>(in_x, in_P, h_25, H_25, NULL, in_z, in_R, in_ea, MAHA_THRESH_25);
}
void car_update_24(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
  update<2, 3, 0>(in_x, in_P, h_24, H_24, NULL, in_z, in_R, in_ea, MAHA_THRESH_24);
}
void car_update_30(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
  update<1, 3, 0>(in_x, in_P, h_30, H_30, NULL, in_z, in_R, in_ea, MAHA_THRESH_30);
}
void car_update_26(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
  update<1, 3, 0>(in_x, in_P, h_26, H_26, NULL, in_z, in_R, in_ea, MAHA_THRESH_26);
}
void car_update_27(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
  update<1, 3, 0>(in_x, in_P, h_27, H_27, NULL, in_z, in_R, in_ea, MAHA_THRESH_27);
}
void car_update_29(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
  update<1, 3, 0>(in_x, in_P, h_29, H_29, NULL, in_z, in_R, in_ea, MAHA_THRESH_29);
}
void car_update_28(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
  update<1, 3, 0>(in_x, in_P, h_28, H_28, NULL, in_z, in_R, in_ea, MAHA_THRESH_28);
}
void car_update_31(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
  update<1, 3, 0>(in_x, in_P, h_31, H_31, NULL, in_z, in_R, in_ea, MAHA_THRESH_31);
}
void car_err_fun(double *nom_x, double *delta_x, double *out_1914130395499138379) {
  err_fun(nom_x, delta_x, out_1914130395499138379);
}
void car_inv_err_fun(double *nom_x, double *true_x, double *out_2456965719472217437) {
  inv_err_fun(nom_x, true_x, out_2456965719472217437);
}
void car_H_mod_fun(double *state, double *out_4779992119913662016) {
  H_mod_fun(state, out_4779992119913662016);
}
void car_f_fun(double *state, double dt, double *out_2238700801222790874) {
  f_fun(state,  dt, out_2238700801222790874);
}
void car_F_fun(double *state, double dt, double *out_3188951629714217883) {
  F_fun(state,  dt, out_3188951629714217883);
}
void car_h_25(double *state, double *unused, double *out_2955682509261511320) {
  h_25(state, unused, out_2955682509261511320);
}
void car_H_25(double *state, double *unused, double *out_7551971679374724871) {
  H_25(state, unused, out_7551971679374724871);
}
void car_h_24(double *state, double *unused, double *out_4912051971117535801) {
  h_24(state, unused, out_4912051971117535801);
}
void car_H_24(double *state, double *unused, double *out_976399872783206770) {
  H_24(state, unused, out_976399872783206770);
}
void car_h_30(double *state, double *unused, double *out_5047370505310505677) {
  h_30(state, unused, out_5047370505310505677);
}
void car_H_30(double *state, double *unused, double *out_635281337883108116) {
  H_30(state, unused, out_635281337883108116);
}
void car_h_26(double *state, double *unused, double *out_1369527850446061802) {
  h_26(state, unused, out_1369527850446061802);
}
void car_H_26(double *state, double *unused, double *out_7153269075460770521) {
  H_26(state, unused, out_7153269075460770521);
}
void car_h_27(double *state, double *unused, double *out_3899192362163173755) {
  h_27(state, unused, out_3899192362163173755);
}
void car_H_27(double *state, double *unused, double *out_2810044649683533027) {
  H_27(state, unused, out_2810044649683533027);
}
void car_h_29(double *state, double *unused, double *out_595359365770451158) {
  h_29(state, unused, out_595359365770451158);
}
void car_H_29(double *state, double *unused, double *out_4523407376553084060) {
  H_29(state, unused, out_4523407376553084060);
}
void car_h_28(double *state, double *unused, double *out_5177906949419811836) {
  h_28(state, unused, out_5177906949419811836);
}
void car_H_28(double *state, double *unused, double *out_8840937680086936982) {
  H_28(state, unused, out_8840937680086936982);
}
void car_h_31(double *state, double *unused, double *out_3815152717088839616) {
  h_31(state, unused, out_3815152717088839616);
}
void car_H_31(double *state, double *unused, double *out_7521325717497764443) {
  H_31(state, unused, out_7521325717497764443);
}
void car_predict(double *in_x, double *in_P, double *in_Q, double dt) {
  predict(in_x, in_P, in_Q, dt);
}
void car_set_mass(double x) {
  set_mass(x);
}
void car_set_rotational_inertia(double x) {
  set_rotational_inertia(x);
}
void car_set_center_to_front(double x) {
  set_center_to_front(x);
}
void car_set_center_to_rear(double x) {
  set_center_to_rear(x);
}
void car_set_stiffness_front(double x) {
  set_stiffness_front(x);
}
void car_set_stiffness_rear(double x) {
  set_stiffness_rear(x);
}
}

const EKF car = {
  .name = "car",
  .kinds = { 25, 24, 30, 26, 27, 29, 28, 31 },
  .feature_kinds = {  },
  .f_fun = car_f_fun,
  .F_fun = car_F_fun,
  .err_fun = car_err_fun,
  .inv_err_fun = car_inv_err_fun,
  .H_mod_fun = car_H_mod_fun,
  .predict = car_predict,
  .hs = {
    { 25, car_h_25 },
    { 24, car_h_24 },
    { 30, car_h_30 },
    { 26, car_h_26 },
    { 27, car_h_27 },
    { 29, car_h_29 },
    { 28, car_h_28 },
    { 31, car_h_31 },
  },
  .Hs = {
    { 25, car_H_25 },
    { 24, car_H_24 },
    { 30, car_H_30 },
    { 26, car_H_26 },
    { 27, car_H_27 },
    { 29, car_H_29 },
    { 28, car_H_28 },
    { 31, car_H_31 },
  },
  .updates = {
    { 25, car_update_25 },
    { 24, car_update_24 },
    { 30, car_update_30 },
    { 26, car_update_26 },
    { 27, car_update_27 },
    { 29, car_update_29 },
    { 28, car_update_28 },
    { 31, car_update_31 },
  },
  .Hes = {
  },
  .sets = {
    { "mass", car_set_mass },
    { "rotational_inertia", car_set_rotational_inertia },
    { "center_to_front", car_set_center_to_front },
    { "center_to_rear", car_set_center_to_rear },
    { "stiffness_front", car_set_stiffness_front },
    { "stiffness_rear", car_set_stiffness_rear },
  },
  .extra_routines = {
  },
};

ekf_init(car);
