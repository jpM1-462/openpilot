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
void live_H(double *in_vec, double *out_8373769282202649548);
void live_err_fun(double *nom_x, double *delta_x, double *out_7843686144394736817);
void live_inv_err_fun(double *nom_x, double *true_x, double *out_4225469103014996173);
void live_H_mod_fun(double *state, double *out_8358142402501933144);
void live_f_fun(double *state, double dt, double *out_8210762158604147154);
void live_F_fun(double *state, double dt, double *out_2457475553137628519);
void live_h_4(double *state, double *unused, double *out_7478004446258192637);
void live_H_4(double *state, double *unused, double *out_3268929401193690479);
void live_h_9(double *state, double *unused, double *out_4165312183178748366);
void live_H_9(double *state, double *unused, double *out_7890595737251413667);
void live_h_10(double *state, double *unused, double *out_1153847551310806058);
void live_H_10(double *state, double *unused, double *out_6310129250639465020);
void live_h_12(double *state, double *unused, double *out_5975007637633465198);
void live_H_12(double *state, double *unused, double *out_3112328975849042517);
void live_h_31(double *state, double *unused, double *out_3876391844886862655);
void live_H_31(double *state, double *unused, double *out_6635591458566297855);
void live_h_32(double *state, double *unused, double *out_8014753800987825621);
void live_H_32(double *state, double *unused, double *out_2416892220894388785);
void live_h_13(double *state, double *unused, double *out_5843516187364789028);
void live_H_13(double *state, double *unused, double *out_8637572793181547907);
void live_h_14(double *state, double *unused, double *out_4165312183178748366);
void live_H_14(double *state, double *unused, double *out_7890595737251413667);
void live_h_33(double *state, double *unused, double *out_6095174459977264046);
void live_H_33(double *state, double *unused, double *out_8660595610504396157);
void live_predict(double *in_x, double *in_P, double *in_Q, double dt);
}