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
void live_H(double *in_vec, double *out_61720393214848937);
void live_err_fun(double *nom_x, double *delta_x, double *out_4824962865715069528);
void live_inv_err_fun(double *nom_x, double *true_x, double *out_457388796125979699);
void live_H_mod_fun(double *state, double *out_8995137806556538300);
void live_f_fun(double *state, double dt, double *out_7757873562176131886);
void live_F_fun(double *state, double dt, double *out_6457543617841644224);
void live_h_4(double *state, double *unused, double *out_7512130010912635017);
void live_H_4(double *state, double *unused, double *out_5757020154937823738);
void live_h_9(double *state, double *unused, double *out_8527478217063711845);
void live_H_9(double *state, double *unused, double *out_1530198780326623732);
void live_h_10(double *state, double *unused, double *out_7202546563935427918);
void live_H_10(double *state, double *unused, double *out_3104246774833475480);
void live_h_12(double *state, double *unused, double *out_866366390415833741);
void live_H_12(double *state, double *unused, double *out_6308465541728994882);
void live_h_31(double *state, double *unused, double *out_2206318948287845659);
void live_H_31(double *state, double *unused, double *out_2390358097565216362);
void live_h_32(double *state, double *unused, double *out_7210263444636478243);
void live_H_32(double *state, double *unused, double *out_7989768510432298595);
void live_h_13(double *state, double *unused, double *out_5085196629364501709);
void live_H_13(double *state, double *unused, double *out_8041898469067582234);
void live_h_14(double *state, double *unused, double *out_8527478217063711845);
void live_H_14(double *state, double *unused, double *out_1530198780326623732);
void live_h_33(double *state, double *unused, double *out_5902300547481657388);
void live_H_33(double *state, double *unused, double *out_7806228195708498067);
void live_predict(double *in_x, double *in_P, double *in_Q, double dt);
}