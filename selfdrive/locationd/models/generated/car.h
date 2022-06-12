#pragma once
#include "rednose/helpers/common_ekf.h"
extern "C" {
void car_update_25(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void car_update_24(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void car_update_30(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void car_update_26(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void car_update_27(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void car_update_29(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void car_update_28(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void car_update_31(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void car_err_fun(double *nom_x, double *delta_x, double *out_4379403857117222854);
void car_inv_err_fun(double *nom_x, double *true_x, double *out_3911892996156877123);
void car_H_mod_fun(double *state, double *out_5860823269405666308);
void car_f_fun(double *state, double dt, double *out_2647164349113807411);
void car_F_fun(double *state, double dt, double *out_7007253849962333144);
void car_h_25(double *state, double *unused, double *out_9084064332289274164);
void car_H_25(double *state, double *unused, double *out_3352101258583786588);
void car_h_24(double *state, double *unused, double *out_4753543947970647466);
void car_H_24(double *state, double *unused, double *out_1921393024329226006);
void car_h_30(double *state, double *unused, double *out_6414953751064792491);
void car_H_30(double *state, double *unused, double *out_833768300076537961);
void car_h_26(double *state, double *unused, double *out_3378167699620459806);
void car_H_26(double *state, double *unused, double *out_7093604577457842812);
void car_h_27(double *state, double *unused, double *out_8140554479387611729);
void car_H_27(double *state, double *unused, double *out_8392183173197731919);
void car_h_29(double *state, double *unused, double *out_5811637866388314974);
void car_H_29(double *state, double *unused, double *out_7369566244397002602);
void car_h_28(double *state, double *unused, double *out_5103768427943279044);
void car_H_28(double *state, double *unused, double *out_5405935972831676351);
void car_h_31(double *state, double *unused, double *out_7960715283269820484);
void car_H_31(double *state, double *unused, double *out_7719812679691194288);
void car_predict(double *in_x, double *in_P, double *in_Q, double dt);
void car_set_mass(double x);
void car_set_rotational_inertia(double x);
void car_set_center_to_front(double x);
void car_set_center_to_rear(double x);
void car_set_stiffness_front(double x);
void car_set_stiffness_rear(double x);
}