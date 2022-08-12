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
void err_fun(double *nom_x, double *delta_x, double *out_2012210987460608007) {
   out_2012210987460608007[0] = delta_x[0] + nom_x[0];
   out_2012210987460608007[1] = delta_x[1] + nom_x[1];
   out_2012210987460608007[2] = delta_x[2] + nom_x[2];
   out_2012210987460608007[3] = delta_x[3] + nom_x[3];
   out_2012210987460608007[4] = delta_x[4] + nom_x[4];
   out_2012210987460608007[5] = delta_x[5] + nom_x[5];
   out_2012210987460608007[6] = delta_x[6] + nom_x[6];
   out_2012210987460608007[7] = delta_x[7] + nom_x[7];
   out_2012210987460608007[8] = delta_x[8] + nom_x[8];
}
void inv_err_fun(double *nom_x, double *true_x, double *out_4452402443006222173) {
   out_4452402443006222173[0] = -nom_x[0] + true_x[0];
   out_4452402443006222173[1] = -nom_x[1] + true_x[1];
   out_4452402443006222173[2] = -nom_x[2] + true_x[2];
   out_4452402443006222173[3] = -nom_x[3] + true_x[3];
   out_4452402443006222173[4] = -nom_x[4] + true_x[4];
   out_4452402443006222173[5] = -nom_x[5] + true_x[5];
   out_4452402443006222173[6] = -nom_x[6] + true_x[6];
   out_4452402443006222173[7] = -nom_x[7] + true_x[7];
   out_4452402443006222173[8] = -nom_x[8] + true_x[8];
}
void H_mod_fun(double *state, double *out_8269490911248935477) {
   out_8269490911248935477[0] = 1.0;
   out_8269490911248935477[1] = 0;
   out_8269490911248935477[2] = 0;
   out_8269490911248935477[3] = 0;
   out_8269490911248935477[4] = 0;
   out_8269490911248935477[5] = 0;
   out_8269490911248935477[6] = 0;
   out_8269490911248935477[7] = 0;
   out_8269490911248935477[8] = 0;
   out_8269490911248935477[9] = 0;
   out_8269490911248935477[10] = 1.0;
   out_8269490911248935477[11] = 0;
   out_8269490911248935477[12] = 0;
   out_8269490911248935477[13] = 0;
   out_8269490911248935477[14] = 0;
   out_8269490911248935477[15] = 0;
   out_8269490911248935477[16] = 0;
   out_8269490911248935477[17] = 0;
   out_8269490911248935477[18] = 0;
   out_8269490911248935477[19] = 0;
   out_8269490911248935477[20] = 1.0;
   out_8269490911248935477[21] = 0;
   out_8269490911248935477[22] = 0;
   out_8269490911248935477[23] = 0;
   out_8269490911248935477[24] = 0;
   out_8269490911248935477[25] = 0;
   out_8269490911248935477[26] = 0;
   out_8269490911248935477[27] = 0;
   out_8269490911248935477[28] = 0;
   out_8269490911248935477[29] = 0;
   out_8269490911248935477[30] = 1.0;
   out_8269490911248935477[31] = 0;
   out_8269490911248935477[32] = 0;
   out_8269490911248935477[33] = 0;
   out_8269490911248935477[34] = 0;
   out_8269490911248935477[35] = 0;
   out_8269490911248935477[36] = 0;
   out_8269490911248935477[37] = 0;
   out_8269490911248935477[38] = 0;
   out_8269490911248935477[39] = 0;
   out_8269490911248935477[40] = 1.0;
   out_8269490911248935477[41] = 0;
   out_8269490911248935477[42] = 0;
   out_8269490911248935477[43] = 0;
   out_8269490911248935477[44] = 0;
   out_8269490911248935477[45] = 0;
   out_8269490911248935477[46] = 0;
   out_8269490911248935477[47] = 0;
   out_8269490911248935477[48] = 0;
   out_8269490911248935477[49] = 0;
   out_8269490911248935477[50] = 1.0;
   out_8269490911248935477[51] = 0;
   out_8269490911248935477[52] = 0;
   out_8269490911248935477[53] = 0;
   out_8269490911248935477[54] = 0;
   out_8269490911248935477[55] = 0;
   out_8269490911248935477[56] = 0;
   out_8269490911248935477[57] = 0;
   out_8269490911248935477[58] = 0;
   out_8269490911248935477[59] = 0;
   out_8269490911248935477[60] = 1.0;
   out_8269490911248935477[61] = 0;
   out_8269490911248935477[62] = 0;
   out_8269490911248935477[63] = 0;
   out_8269490911248935477[64] = 0;
   out_8269490911248935477[65] = 0;
   out_8269490911248935477[66] = 0;
   out_8269490911248935477[67] = 0;
   out_8269490911248935477[68] = 0;
   out_8269490911248935477[69] = 0;
   out_8269490911248935477[70] = 1.0;
   out_8269490911248935477[71] = 0;
   out_8269490911248935477[72] = 0;
   out_8269490911248935477[73] = 0;
   out_8269490911248935477[74] = 0;
   out_8269490911248935477[75] = 0;
   out_8269490911248935477[76] = 0;
   out_8269490911248935477[77] = 0;
   out_8269490911248935477[78] = 0;
   out_8269490911248935477[79] = 0;
   out_8269490911248935477[80] = 1.0;
}
void f_fun(double *state, double dt, double *out_4375872833818994997) {
   out_4375872833818994997[0] = state[0];
   out_4375872833818994997[1] = state[1];
   out_4375872833818994997[2] = state[2];
   out_4375872833818994997[3] = state[3];
   out_4375872833818994997[4] = state[4];
   out_4375872833818994997[5] = dt*((-state[4] + (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(mass*state[4]))*state[6] - 9.8000000000000007*state[8] + stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(mass*state[1]) + (-stiffness_front*state[0] - stiffness_rear*state[0])*state[5]/(mass*state[4])) + state[5];
   out_4375872833818994997[6] = dt*(center_to_front*stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(rotational_inertia*state[1]) + (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])*state[5]/(rotational_inertia*state[4]) + (-pow(center_to_front, 2)*stiffness_front*state[0] - pow(center_to_rear, 2)*stiffness_rear*state[0])*state[6]/(rotational_inertia*state[4])) + state[6];
   out_4375872833818994997[7] = state[7];
   out_4375872833818994997[8] = state[8];
}
void F_fun(double *state, double dt, double *out_786040650709582685) {
   out_786040650709582685[0] = 1;
   out_786040650709582685[1] = 0;
   out_786040650709582685[2] = 0;
   out_786040650709582685[3] = 0;
   out_786040650709582685[4] = 0;
   out_786040650709582685[5] = 0;
   out_786040650709582685[6] = 0;
   out_786040650709582685[7] = 0;
   out_786040650709582685[8] = 0;
   out_786040650709582685[9] = 0;
   out_786040650709582685[10] = 1;
   out_786040650709582685[11] = 0;
   out_786040650709582685[12] = 0;
   out_786040650709582685[13] = 0;
   out_786040650709582685[14] = 0;
   out_786040650709582685[15] = 0;
   out_786040650709582685[16] = 0;
   out_786040650709582685[17] = 0;
   out_786040650709582685[18] = 0;
   out_786040650709582685[19] = 0;
   out_786040650709582685[20] = 1;
   out_786040650709582685[21] = 0;
   out_786040650709582685[22] = 0;
   out_786040650709582685[23] = 0;
   out_786040650709582685[24] = 0;
   out_786040650709582685[25] = 0;
   out_786040650709582685[26] = 0;
   out_786040650709582685[27] = 0;
   out_786040650709582685[28] = 0;
   out_786040650709582685[29] = 0;
   out_786040650709582685[30] = 1;
   out_786040650709582685[31] = 0;
   out_786040650709582685[32] = 0;
   out_786040650709582685[33] = 0;
   out_786040650709582685[34] = 0;
   out_786040650709582685[35] = 0;
   out_786040650709582685[36] = 0;
   out_786040650709582685[37] = 0;
   out_786040650709582685[38] = 0;
   out_786040650709582685[39] = 0;
   out_786040650709582685[40] = 1;
   out_786040650709582685[41] = 0;
   out_786040650709582685[42] = 0;
   out_786040650709582685[43] = 0;
   out_786040650709582685[44] = 0;
   out_786040650709582685[45] = dt*(stiffness_front*(-state[2] - state[3] + state[7])/(mass*state[1]) + (-stiffness_front - stiffness_rear)*state[5]/(mass*state[4]) + (-center_to_front*stiffness_front + center_to_rear*stiffness_rear)*state[6]/(mass*state[4]));
   out_786040650709582685[46] = -dt*stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(mass*pow(state[1], 2));
   out_786040650709582685[47] = -dt*stiffness_front*state[0]/(mass*state[1]);
   out_786040650709582685[48] = -dt*stiffness_front*state[0]/(mass*state[1]);
   out_786040650709582685[49] = dt*((-1 - (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(mass*pow(state[4], 2)))*state[6] - (-stiffness_front*state[0] - stiffness_rear*state[0])*state[5]/(mass*pow(state[4], 2)));
   out_786040650709582685[50] = dt*(-stiffness_front*state[0] - stiffness_rear*state[0])/(mass*state[4]) + 1;
   out_786040650709582685[51] = dt*(-state[4] + (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(mass*state[4]));
   out_786040650709582685[52] = dt*stiffness_front*state[0]/(mass*state[1]);
   out_786040650709582685[53] = -9.8000000000000007*dt;
   out_786040650709582685[54] = dt*(center_to_front*stiffness_front*(-state[2] - state[3] + state[7])/(rotational_inertia*state[1]) + (-center_to_front*stiffness_front + center_to_rear*stiffness_rear)*state[5]/(rotational_inertia*state[4]) + (-pow(center_to_front, 2)*stiffness_front - pow(center_to_rear, 2)*stiffness_rear)*state[6]/(rotational_inertia*state[4]));
   out_786040650709582685[55] = -center_to_front*dt*stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(rotational_inertia*pow(state[1], 2));
   out_786040650709582685[56] = -center_to_front*dt*stiffness_front*state[0]/(rotational_inertia*state[1]);
   out_786040650709582685[57] = -center_to_front*dt*stiffness_front*state[0]/(rotational_inertia*state[1]);
   out_786040650709582685[58] = dt*(-(-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])*state[5]/(rotational_inertia*pow(state[4], 2)) - (-pow(center_to_front, 2)*stiffness_front*state[0] - pow(center_to_rear, 2)*stiffness_rear*state[0])*state[6]/(rotational_inertia*pow(state[4], 2)));
   out_786040650709582685[59] = dt*(-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(rotational_inertia*state[4]);
   out_786040650709582685[60] = dt*(-pow(center_to_front, 2)*stiffness_front*state[0] - pow(center_to_rear, 2)*stiffness_rear*state[0])/(rotational_inertia*state[4]) + 1;
   out_786040650709582685[61] = center_to_front*dt*stiffness_front*state[0]/(rotational_inertia*state[1]);
   out_786040650709582685[62] = 0;
   out_786040650709582685[63] = 0;
   out_786040650709582685[64] = 0;
   out_786040650709582685[65] = 0;
   out_786040650709582685[66] = 0;
   out_786040650709582685[67] = 0;
   out_786040650709582685[68] = 0;
   out_786040650709582685[69] = 0;
   out_786040650709582685[70] = 1;
   out_786040650709582685[71] = 0;
   out_786040650709582685[72] = 0;
   out_786040650709582685[73] = 0;
   out_786040650709582685[74] = 0;
   out_786040650709582685[75] = 0;
   out_786040650709582685[76] = 0;
   out_786040650709582685[77] = 0;
   out_786040650709582685[78] = 0;
   out_786040650709582685[79] = 0;
   out_786040650709582685[80] = 1;
}
void h_25(double *state, double *unused, double *out_6497350029804809175) {
   out_6497350029804809175[0] = state[6];
}
void H_25(double *state, double *unused, double *out_7436575763441244356) {
   out_7436575763441244356[0] = 0;
   out_7436575763441244356[1] = 0;
   out_7436575763441244356[2] = 0;
   out_7436575763441244356[3] = 0;
   out_7436575763441244356[4] = 0;
   out_7436575763441244356[5] = 0;
   out_7436575763441244356[6] = 1;
   out_7436575763441244356[7] = 0;
   out_7436575763441244356[8] = 0;
}
void h_24(double *state, double *unused, double *out_3081959780956341134) {
   out_3081959780956341134[0] = state[4];
   out_3081959780956341134[1] = state[5];
}
void H_24(double *state, double *unused, double *out_7452242711465984124) {
   out_7452242711465984124[0] = 0;
   out_7452242711465984124[1] = 0;
   out_7452242711465984124[2] = 0;
   out_7452242711465984124[3] = 0;
   out_7452242711465984124[4] = 1;
   out_7452242711465984124[5] = 0;
   out_7452242711465984124[6] = 0;
   out_7452242711465984124[7] = 0;
   out_7452242711465984124[8] = 0;
   out_7452242711465984124[9] = 0;
   out_7452242711465984124[10] = 0;
   out_7452242711465984124[11] = 0;
   out_7452242711465984124[12] = 0;
   out_7452242711465984124[13] = 0;
   out_7452242711465984124[14] = 1;
   out_7452242711465984124[15] = 0;
   out_7452242711465984124[16] = 0;
   out_7452242711465984124[17] = 0;
}
void h_30(double *state, double *unused, double *out_4068121393150693256) {
   out_4068121393150693256[0] = state[4];
}
void H_30(double *state, double *unused, double *out_8491835351761058633) {
   out_8491835351761058633[0] = 0;
   out_8491835351761058633[1] = 0;
   out_8491835351761058633[2] = 0;
   out_8491835351761058633[3] = 0;
   out_8491835351761058633[4] = 1;
   out_8491835351761058633[5] = 0;
   out_8491835351761058633[6] = 0;
   out_8491835351761058633[7] = 0;
   out_8491835351761058633[8] = 0;
}
void h_26(double *state, double *unused, double *out_9069794660811538901) {
   out_9069794660811538901[0] = state[7];
}
void H_26(double *state, double *unused, double *out_3695072444567188132) {
   out_3695072444567188132[0] = 0;
   out_3695072444567188132[1] = 0;
   out_3695072444567188132[2] = 0;
   out_3695072444567188132[3] = 0;
   out_3695072444567188132[4] = 0;
   out_3695072444567188132[5] = 0;
   out_3695072444567188132[6] = 0;
   out_3695072444567188132[7] = 1;
   out_3695072444567188132[8] = 0;
}
void h_27(double *state, double *unused, double *out_1012265998463976720) {
   out_1012265998463976720[0] = state[3];
}
void H_27(double *state, double *unused, double *out_734116121513211247) {
   out_734116121513211247[0] = 0;
   out_734116121513211247[1] = 0;
   out_734116121513211247[2] = 0;
   out_734116121513211247[3] = 1;
   out_734116121513211247[4] = 0;
   out_734116121513211247[5] = 0;
   out_734116121513211247[6] = 0;
   out_734116121513211247[7] = 0;
   out_734116121513211247[8] = 0;
}
void h_29(double *state, double *unused, double *out_7484394682799763272) {
   out_7484394682799763272[0] = state[1];
}
void H_29(double *state, double *unused, double *out_3419110777628028342) {
   out_3419110777628028342[0] = 0;
   out_3419110777628028342[1] = 1;
   out_3419110777628028342[2] = 0;
   out_3419110777628028342[3] = 0;
   out_3419110777628028342[4] = 0;
   out_3419110777628028342[5] = 0;
   out_3419110777628028342[6] = 0;
   out_3419110777628028342[7] = 0;
   out_3419110777628028342[8] = 0;
}
void h_28(double *state, double *unused, double *out_6157156061579624257) {
   out_6157156061579624257[0] = state[0];
}
void H_28(double *state, double *unused, double *out_5382741049193354593) {
   out_5382741049193354593[0] = 1;
   out_5382741049193354593[1] = 0;
   out_5382741049193354593[2] = 0;
   out_5382741049193354593[3] = 0;
   out_5382741049193354593[4] = 0;
   out_5382741049193354593[5] = 0;
   out_5382741049193354593[6] = 0;
   out_5382741049193354593[7] = 0;
   out_5382741049193354593[8] = 0;
}
void h_31(double *state, double *unused, double *out_2859900501550786825) {
   out_2859900501550786825[0] = state[8];
}
void H_31(double *state, double *unused, double *out_3068864342333836656) {
   out_3068864342333836656[0] = 0;
   out_3068864342333836656[1] = 0;
   out_3068864342333836656[2] = 0;
   out_3068864342333836656[3] = 0;
   out_3068864342333836656[4] = 0;
   out_3068864342333836656[5] = 0;
   out_3068864342333836656[6] = 0;
   out_3068864342333836656[7] = 0;
   out_3068864342333836656[8] = 1;
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
void car_err_fun(double *nom_x, double *delta_x, double *out_2012210987460608007) {
  err_fun(nom_x, delta_x, out_2012210987460608007);
}
void car_inv_err_fun(double *nom_x, double *true_x, double *out_4452402443006222173) {
  inv_err_fun(nom_x, true_x, out_4452402443006222173);
}
void car_H_mod_fun(double *state, double *out_8269490911248935477) {
  H_mod_fun(state, out_8269490911248935477);
}
void car_f_fun(double *state, double dt, double *out_4375872833818994997) {
  f_fun(state,  dt, out_4375872833818994997);
}
void car_F_fun(double *state, double dt, double *out_786040650709582685) {
  F_fun(state,  dt, out_786040650709582685);
}
void car_h_25(double *state, double *unused, double *out_6497350029804809175) {
  h_25(state, unused, out_6497350029804809175);
}
void car_H_25(double *state, double *unused, double *out_7436575763441244356) {
  H_25(state, unused, out_7436575763441244356);
}
void car_h_24(double *state, double *unused, double *out_3081959780956341134) {
  h_24(state, unused, out_3081959780956341134);
}
void car_H_24(double *state, double *unused, double *out_7452242711465984124) {
  H_24(state, unused, out_7452242711465984124);
}
void car_h_30(double *state, double *unused, double *out_4068121393150693256) {
  h_30(state, unused, out_4068121393150693256);
}
void car_H_30(double *state, double *unused, double *out_8491835351761058633) {
  H_30(state, unused, out_8491835351761058633);
}
void car_h_26(double *state, double *unused, double *out_9069794660811538901) {
  h_26(state, unused, out_9069794660811538901);
}
void car_H_26(double *state, double *unused, double *out_3695072444567188132) {
  H_26(state, unused, out_3695072444567188132);
}
void car_h_27(double *state, double *unused, double *out_1012265998463976720) {
  h_27(state, unused, out_1012265998463976720);
}
void car_H_27(double *state, double *unused, double *out_734116121513211247) {
  H_27(state, unused, out_734116121513211247);
}
void car_h_29(double *state, double *unused, double *out_7484394682799763272) {
  h_29(state, unused, out_7484394682799763272);
}
void car_H_29(double *state, double *unused, double *out_3419110777628028342) {
  H_29(state, unused, out_3419110777628028342);
}
void car_h_28(double *state, double *unused, double *out_6157156061579624257) {
  h_28(state, unused, out_6157156061579624257);
}
void car_H_28(double *state, double *unused, double *out_5382741049193354593) {
  H_28(state, unused, out_5382741049193354593);
}
void car_h_31(double *state, double *unused, double *out_2859900501550786825) {
  h_31(state, unused, out_2859900501550786825);
}
void car_H_31(double *state, double *unused, double *out_3068864342333836656) {
  H_31(state, unused, out_3068864342333836656);
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
