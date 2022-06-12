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
void err_fun(double *nom_x, double *delta_x, double *out_4379403857117222854) {
   out_4379403857117222854[0] = delta_x[0] + nom_x[0];
   out_4379403857117222854[1] = delta_x[1] + nom_x[1];
   out_4379403857117222854[2] = delta_x[2] + nom_x[2];
   out_4379403857117222854[3] = delta_x[3] + nom_x[3];
   out_4379403857117222854[4] = delta_x[4] + nom_x[4];
   out_4379403857117222854[5] = delta_x[5] + nom_x[5];
   out_4379403857117222854[6] = delta_x[6] + nom_x[6];
   out_4379403857117222854[7] = delta_x[7] + nom_x[7];
   out_4379403857117222854[8] = delta_x[8] + nom_x[8];
}
void inv_err_fun(double *nom_x, double *true_x, double *out_3911892996156877123) {
   out_3911892996156877123[0] = -nom_x[0] + true_x[0];
   out_3911892996156877123[1] = -nom_x[1] + true_x[1];
   out_3911892996156877123[2] = -nom_x[2] + true_x[2];
   out_3911892996156877123[3] = -nom_x[3] + true_x[3];
   out_3911892996156877123[4] = -nom_x[4] + true_x[4];
   out_3911892996156877123[5] = -nom_x[5] + true_x[5];
   out_3911892996156877123[6] = -nom_x[6] + true_x[6];
   out_3911892996156877123[7] = -nom_x[7] + true_x[7];
   out_3911892996156877123[8] = -nom_x[8] + true_x[8];
}
void H_mod_fun(double *state, double *out_5860823269405666308) {
   out_5860823269405666308[0] = 1.0;
   out_5860823269405666308[1] = 0;
   out_5860823269405666308[2] = 0;
   out_5860823269405666308[3] = 0;
   out_5860823269405666308[4] = 0;
   out_5860823269405666308[5] = 0;
   out_5860823269405666308[6] = 0;
   out_5860823269405666308[7] = 0;
   out_5860823269405666308[8] = 0;
   out_5860823269405666308[9] = 0;
   out_5860823269405666308[10] = 1.0;
   out_5860823269405666308[11] = 0;
   out_5860823269405666308[12] = 0;
   out_5860823269405666308[13] = 0;
   out_5860823269405666308[14] = 0;
   out_5860823269405666308[15] = 0;
   out_5860823269405666308[16] = 0;
   out_5860823269405666308[17] = 0;
   out_5860823269405666308[18] = 0;
   out_5860823269405666308[19] = 0;
   out_5860823269405666308[20] = 1.0;
   out_5860823269405666308[21] = 0;
   out_5860823269405666308[22] = 0;
   out_5860823269405666308[23] = 0;
   out_5860823269405666308[24] = 0;
   out_5860823269405666308[25] = 0;
   out_5860823269405666308[26] = 0;
   out_5860823269405666308[27] = 0;
   out_5860823269405666308[28] = 0;
   out_5860823269405666308[29] = 0;
   out_5860823269405666308[30] = 1.0;
   out_5860823269405666308[31] = 0;
   out_5860823269405666308[32] = 0;
   out_5860823269405666308[33] = 0;
   out_5860823269405666308[34] = 0;
   out_5860823269405666308[35] = 0;
   out_5860823269405666308[36] = 0;
   out_5860823269405666308[37] = 0;
   out_5860823269405666308[38] = 0;
   out_5860823269405666308[39] = 0;
   out_5860823269405666308[40] = 1.0;
   out_5860823269405666308[41] = 0;
   out_5860823269405666308[42] = 0;
   out_5860823269405666308[43] = 0;
   out_5860823269405666308[44] = 0;
   out_5860823269405666308[45] = 0;
   out_5860823269405666308[46] = 0;
   out_5860823269405666308[47] = 0;
   out_5860823269405666308[48] = 0;
   out_5860823269405666308[49] = 0;
   out_5860823269405666308[50] = 1.0;
   out_5860823269405666308[51] = 0;
   out_5860823269405666308[52] = 0;
   out_5860823269405666308[53] = 0;
   out_5860823269405666308[54] = 0;
   out_5860823269405666308[55] = 0;
   out_5860823269405666308[56] = 0;
   out_5860823269405666308[57] = 0;
   out_5860823269405666308[58] = 0;
   out_5860823269405666308[59] = 0;
   out_5860823269405666308[60] = 1.0;
   out_5860823269405666308[61] = 0;
   out_5860823269405666308[62] = 0;
   out_5860823269405666308[63] = 0;
   out_5860823269405666308[64] = 0;
   out_5860823269405666308[65] = 0;
   out_5860823269405666308[66] = 0;
   out_5860823269405666308[67] = 0;
   out_5860823269405666308[68] = 0;
   out_5860823269405666308[69] = 0;
   out_5860823269405666308[70] = 1.0;
   out_5860823269405666308[71] = 0;
   out_5860823269405666308[72] = 0;
   out_5860823269405666308[73] = 0;
   out_5860823269405666308[74] = 0;
   out_5860823269405666308[75] = 0;
   out_5860823269405666308[76] = 0;
   out_5860823269405666308[77] = 0;
   out_5860823269405666308[78] = 0;
   out_5860823269405666308[79] = 0;
   out_5860823269405666308[80] = 1.0;
}
void f_fun(double *state, double dt, double *out_2647164349113807411) {
   out_2647164349113807411[0] = state[0];
   out_2647164349113807411[1] = state[1];
   out_2647164349113807411[2] = state[2];
   out_2647164349113807411[3] = state[3];
   out_2647164349113807411[4] = state[4];
   out_2647164349113807411[5] = dt*((-state[4] + (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(mass*state[4]))*state[6] - 9.8000000000000007*state[8] + stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(mass*state[1]) + (-stiffness_front*state[0] - stiffness_rear*state[0])*state[5]/(mass*state[4])) + state[5];
   out_2647164349113807411[6] = dt*(center_to_front*stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(rotational_inertia*state[1]) + (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])*state[5]/(rotational_inertia*state[4]) + (-pow(center_to_front, 2)*stiffness_front*state[0] - pow(center_to_rear, 2)*stiffness_rear*state[0])*state[6]/(rotational_inertia*state[4])) + state[6];
   out_2647164349113807411[7] = state[7];
   out_2647164349113807411[8] = state[8];
}
void F_fun(double *state, double dt, double *out_7007253849962333144) {
   out_7007253849962333144[0] = 1;
   out_7007253849962333144[1] = 0;
   out_7007253849962333144[2] = 0;
   out_7007253849962333144[3] = 0;
   out_7007253849962333144[4] = 0;
   out_7007253849962333144[5] = 0;
   out_7007253849962333144[6] = 0;
   out_7007253849962333144[7] = 0;
   out_7007253849962333144[8] = 0;
   out_7007253849962333144[9] = 0;
   out_7007253849962333144[10] = 1;
   out_7007253849962333144[11] = 0;
   out_7007253849962333144[12] = 0;
   out_7007253849962333144[13] = 0;
   out_7007253849962333144[14] = 0;
   out_7007253849962333144[15] = 0;
   out_7007253849962333144[16] = 0;
   out_7007253849962333144[17] = 0;
   out_7007253849962333144[18] = 0;
   out_7007253849962333144[19] = 0;
   out_7007253849962333144[20] = 1;
   out_7007253849962333144[21] = 0;
   out_7007253849962333144[22] = 0;
   out_7007253849962333144[23] = 0;
   out_7007253849962333144[24] = 0;
   out_7007253849962333144[25] = 0;
   out_7007253849962333144[26] = 0;
   out_7007253849962333144[27] = 0;
   out_7007253849962333144[28] = 0;
   out_7007253849962333144[29] = 0;
   out_7007253849962333144[30] = 1;
   out_7007253849962333144[31] = 0;
   out_7007253849962333144[32] = 0;
   out_7007253849962333144[33] = 0;
   out_7007253849962333144[34] = 0;
   out_7007253849962333144[35] = 0;
   out_7007253849962333144[36] = 0;
   out_7007253849962333144[37] = 0;
   out_7007253849962333144[38] = 0;
   out_7007253849962333144[39] = 0;
   out_7007253849962333144[40] = 1;
   out_7007253849962333144[41] = 0;
   out_7007253849962333144[42] = 0;
   out_7007253849962333144[43] = 0;
   out_7007253849962333144[44] = 0;
   out_7007253849962333144[45] = dt*(stiffness_front*(-state[2] - state[3] + state[7])/(mass*state[1]) + (-stiffness_front - stiffness_rear)*state[5]/(mass*state[4]) + (-center_to_front*stiffness_front + center_to_rear*stiffness_rear)*state[6]/(mass*state[4]));
   out_7007253849962333144[46] = -dt*stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(mass*pow(state[1], 2));
   out_7007253849962333144[47] = -dt*stiffness_front*state[0]/(mass*state[1]);
   out_7007253849962333144[48] = -dt*stiffness_front*state[0]/(mass*state[1]);
   out_7007253849962333144[49] = dt*((-1 - (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(mass*pow(state[4], 2)))*state[6] - (-stiffness_front*state[0] - stiffness_rear*state[0])*state[5]/(mass*pow(state[4], 2)));
   out_7007253849962333144[50] = dt*(-stiffness_front*state[0] - stiffness_rear*state[0])/(mass*state[4]) + 1;
   out_7007253849962333144[51] = dt*(-state[4] + (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(mass*state[4]));
   out_7007253849962333144[52] = dt*stiffness_front*state[0]/(mass*state[1]);
   out_7007253849962333144[53] = -9.8000000000000007*dt;
   out_7007253849962333144[54] = dt*(center_to_front*stiffness_front*(-state[2] - state[3] + state[7])/(rotational_inertia*state[1]) + (-center_to_front*stiffness_front + center_to_rear*stiffness_rear)*state[5]/(rotational_inertia*state[4]) + (-pow(center_to_front, 2)*stiffness_front - pow(center_to_rear, 2)*stiffness_rear)*state[6]/(rotational_inertia*state[4]));
   out_7007253849962333144[55] = -center_to_front*dt*stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(rotational_inertia*pow(state[1], 2));
   out_7007253849962333144[56] = -center_to_front*dt*stiffness_front*state[0]/(rotational_inertia*state[1]);
   out_7007253849962333144[57] = -center_to_front*dt*stiffness_front*state[0]/(rotational_inertia*state[1]);
   out_7007253849962333144[58] = dt*(-(-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])*state[5]/(rotational_inertia*pow(state[4], 2)) - (-pow(center_to_front, 2)*stiffness_front*state[0] - pow(center_to_rear, 2)*stiffness_rear*state[0])*state[6]/(rotational_inertia*pow(state[4], 2)));
   out_7007253849962333144[59] = dt*(-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(rotational_inertia*state[4]);
   out_7007253849962333144[60] = dt*(-pow(center_to_front, 2)*stiffness_front*state[0] - pow(center_to_rear, 2)*stiffness_rear*state[0])/(rotational_inertia*state[4]) + 1;
   out_7007253849962333144[61] = center_to_front*dt*stiffness_front*state[0]/(rotational_inertia*state[1]);
   out_7007253849962333144[62] = 0;
   out_7007253849962333144[63] = 0;
   out_7007253849962333144[64] = 0;
   out_7007253849962333144[65] = 0;
   out_7007253849962333144[66] = 0;
   out_7007253849962333144[67] = 0;
   out_7007253849962333144[68] = 0;
   out_7007253849962333144[69] = 0;
   out_7007253849962333144[70] = 1;
   out_7007253849962333144[71] = 0;
   out_7007253849962333144[72] = 0;
   out_7007253849962333144[73] = 0;
   out_7007253849962333144[74] = 0;
   out_7007253849962333144[75] = 0;
   out_7007253849962333144[76] = 0;
   out_7007253849962333144[77] = 0;
   out_7007253849962333144[78] = 0;
   out_7007253849962333144[79] = 0;
   out_7007253849962333144[80] = 1;
}
void h_25(double *state, double *unused, double *out_9084064332289274164) {
   out_9084064332289274164[0] = state[6];
}
void H_25(double *state, double *unused, double *out_3352101258583786588) {
   out_3352101258583786588[0] = 0;
   out_3352101258583786588[1] = 0;
   out_3352101258583786588[2] = 0;
   out_3352101258583786588[3] = 0;
   out_3352101258583786588[4] = 0;
   out_3352101258583786588[5] = 0;
   out_3352101258583786588[6] = 1;
   out_3352101258583786588[7] = 0;
   out_3352101258583786588[8] = 0;
}
void h_24(double *state, double *unused, double *out_4753543947970647466) {
   out_4753543947970647466[0] = state[4];
   out_4753543947970647466[1] = state[5];
}
void H_24(double *state, double *unused, double *out_1921393024329226006) {
   out_1921393024329226006[0] = 0;
   out_1921393024329226006[1] = 0;
   out_1921393024329226006[2] = 0;
   out_1921393024329226006[3] = 0;
   out_1921393024329226006[4] = 1;
   out_1921393024329226006[5] = 0;
   out_1921393024329226006[6] = 0;
   out_1921393024329226006[7] = 0;
   out_1921393024329226006[8] = 0;
   out_1921393024329226006[9] = 0;
   out_1921393024329226006[10] = 0;
   out_1921393024329226006[11] = 0;
   out_1921393024329226006[12] = 0;
   out_1921393024329226006[13] = 0;
   out_1921393024329226006[14] = 1;
   out_1921393024329226006[15] = 0;
   out_1921393024329226006[16] = 0;
   out_1921393024329226006[17] = 0;
}
void h_30(double *state, double *unused, double *out_6414953751064792491) {
   out_6414953751064792491[0] = state[4];
}
void H_30(double *state, double *unused, double *out_833768300076537961) {
   out_833768300076537961[0] = 0;
   out_833768300076537961[1] = 0;
   out_833768300076537961[2] = 0;
   out_833768300076537961[3] = 0;
   out_833768300076537961[4] = 1;
   out_833768300076537961[5] = 0;
   out_833768300076537961[6] = 0;
   out_833768300076537961[7] = 0;
   out_833768300076537961[8] = 0;
}
void h_26(double *state, double *unused, double *out_3378167699620459806) {
   out_3378167699620459806[0] = state[7];
}
void H_26(double *state, double *unused, double *out_7093604577457842812) {
   out_7093604577457842812[0] = 0;
   out_7093604577457842812[1] = 0;
   out_7093604577457842812[2] = 0;
   out_7093604577457842812[3] = 0;
   out_7093604577457842812[4] = 0;
   out_7093604577457842812[5] = 0;
   out_7093604577457842812[6] = 0;
   out_7093604577457842812[7] = 1;
   out_7093604577457842812[8] = 0;
}
void h_27(double *state, double *unused, double *out_8140554479387611729) {
   out_8140554479387611729[0] = state[3];
}
void H_27(double *state, double *unused, double *out_8392183173197731919) {
   out_8392183173197731919[0] = 0;
   out_8392183173197731919[1] = 0;
   out_8392183173197731919[2] = 0;
   out_8392183173197731919[3] = 1;
   out_8392183173197731919[4] = 0;
   out_8392183173197731919[5] = 0;
   out_8392183173197731919[6] = 0;
   out_8392183173197731919[7] = 0;
   out_8392183173197731919[8] = 0;
}
void h_29(double *state, double *unused, double *out_5811637866388314974) {
   out_5811637866388314974[0] = state[1];
}
void H_29(double *state, double *unused, double *out_7369566244397002602) {
   out_7369566244397002602[0] = 0;
   out_7369566244397002602[1] = 1;
   out_7369566244397002602[2] = 0;
   out_7369566244397002602[3] = 0;
   out_7369566244397002602[4] = 0;
   out_7369566244397002602[5] = 0;
   out_7369566244397002602[6] = 0;
   out_7369566244397002602[7] = 0;
   out_7369566244397002602[8] = 0;
}
void h_28(double *state, double *unused, double *out_5103768427943279044) {
   out_5103768427943279044[0] = state[0];
}
void H_28(double *state, double *unused, double *out_5405935972831676351) {
   out_5405935972831676351[0] = 1;
   out_5405935972831676351[1] = 0;
   out_5405935972831676351[2] = 0;
   out_5405935972831676351[3] = 0;
   out_5405935972831676351[4] = 0;
   out_5405935972831676351[5] = 0;
   out_5405935972831676351[6] = 0;
   out_5405935972831676351[7] = 0;
   out_5405935972831676351[8] = 0;
}
void h_31(double *state, double *unused, double *out_7960715283269820484) {
   out_7960715283269820484[0] = state[8];
}
void H_31(double *state, double *unused, double *out_7719812679691194288) {
   out_7719812679691194288[0] = 0;
   out_7719812679691194288[1] = 0;
   out_7719812679691194288[2] = 0;
   out_7719812679691194288[3] = 0;
   out_7719812679691194288[4] = 0;
   out_7719812679691194288[5] = 0;
   out_7719812679691194288[6] = 0;
   out_7719812679691194288[7] = 0;
   out_7719812679691194288[8] = 1;
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
void car_err_fun(double *nom_x, double *delta_x, double *out_4379403857117222854) {
  err_fun(nom_x, delta_x, out_4379403857117222854);
}
void car_inv_err_fun(double *nom_x, double *true_x, double *out_3911892996156877123) {
  inv_err_fun(nom_x, true_x, out_3911892996156877123);
}
void car_H_mod_fun(double *state, double *out_5860823269405666308) {
  H_mod_fun(state, out_5860823269405666308);
}
void car_f_fun(double *state, double dt, double *out_2647164349113807411) {
  f_fun(state,  dt, out_2647164349113807411);
}
void car_F_fun(double *state, double dt, double *out_7007253849962333144) {
  F_fun(state,  dt, out_7007253849962333144);
}
void car_h_25(double *state, double *unused, double *out_9084064332289274164) {
  h_25(state, unused, out_9084064332289274164);
}
void car_H_25(double *state, double *unused, double *out_3352101258583786588) {
  H_25(state, unused, out_3352101258583786588);
}
void car_h_24(double *state, double *unused, double *out_4753543947970647466) {
  h_24(state, unused, out_4753543947970647466);
}
void car_H_24(double *state, double *unused, double *out_1921393024329226006) {
  H_24(state, unused, out_1921393024329226006);
}
void car_h_30(double *state, double *unused, double *out_6414953751064792491) {
  h_30(state, unused, out_6414953751064792491);
}
void car_H_30(double *state, double *unused, double *out_833768300076537961) {
  H_30(state, unused, out_833768300076537961);
}
void car_h_26(double *state, double *unused, double *out_3378167699620459806) {
  h_26(state, unused, out_3378167699620459806);
}
void car_H_26(double *state, double *unused, double *out_7093604577457842812) {
  H_26(state, unused, out_7093604577457842812);
}
void car_h_27(double *state, double *unused, double *out_8140554479387611729) {
  h_27(state, unused, out_8140554479387611729);
}
void car_H_27(double *state, double *unused, double *out_8392183173197731919) {
  H_27(state, unused, out_8392183173197731919);
}
void car_h_29(double *state, double *unused, double *out_5811637866388314974) {
  h_29(state, unused, out_5811637866388314974);
}
void car_H_29(double *state, double *unused, double *out_7369566244397002602) {
  H_29(state, unused, out_7369566244397002602);
}
void car_h_28(double *state, double *unused, double *out_5103768427943279044) {
  h_28(state, unused, out_5103768427943279044);
}
void car_H_28(double *state, double *unused, double *out_5405935972831676351) {
  H_28(state, unused, out_5405935972831676351);
}
void car_h_31(double *state, double *unused, double *out_7960715283269820484) {
  h_31(state, unused, out_7960715283269820484);
}
void car_H_31(double *state, double *unused, double *out_7719812679691194288) {
  H_31(state, unused, out_7719812679691194288);
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
