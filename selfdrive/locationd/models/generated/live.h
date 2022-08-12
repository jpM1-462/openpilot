#pragma once
#include "rednose/helpers/common_ekf.h"
extern "C" {
void live_update_4(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void live_update_9(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void live_update_10(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void live_update_12(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void live_update_31(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void live_update_32(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void live_update_13(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void live_update_14(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void live_update_33(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void live_H(double *in_vec, double *out_7296470732038239115);
void live_err_fun(double *nom_x, double *delta_x, double *out_195928949731378623);
void live_inv_err_fun(double *nom_x, double *true_x, double *out_8048527368191010135);
void live_H_mod_fun(double *state, double *out_7785589846533456051);
void live_f_fun(double *state, double dt, double *out_6180571892829147108);
void live_F_fun(double *state, double dt, double *out_7713953993130752993);
void live_h_4(double *state, double *unused, double *out_9147143317146732236);
void live_H_4(double *state, double *unused, double *out_5140815643046409326);
void live_h_9(double *state, double *unused, double *out_3325647758290317255);
void live_H_9(double *state, double *unused, double *out_9148760694308364807);
void live_h_10(double *state, double *unused, double *out_7331233491917205272);
void live_H_10(double *state, double *unused, double *out_451170121987628608);
void live_h_12(double *state, double *unused, double *out_3549419012567469284);
void live_H_12(double *state, double *unused, double *out_4519716617998815659);
void live_h_31(double *state, double *unused, double *out_7170478479227021472);
void live_H_31(double *state, double *unused, double *out_1774153585673801950);
void live_h_32(double *state, double *unused, double *out_2516646216923074189);
void live_H_32(double *state, double *unused, double *out_8689018089367918306);
void live_h_13(double *state, double *unused, double *out_6753228495098360460);
void live_H_13(double *state, double *unused, double *out_3939471430803101762);
void live_h_14(double *state, double *unused, double *out_3325647758290317255);
void live_H_14(double *state, double *unused, double *out_9148760694308364807);
void live_h_33(double *state, double *unused, double *out_6932657845551980812);
void live_H_33(double *state, double *unused, double *out_1376403418965055654);
void live_predict(double *in_x, double *in_P, double *in_Q, double dt);
}