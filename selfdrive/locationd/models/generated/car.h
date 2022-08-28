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
void car_err_fun(double *nom_x, double *delta_x, double *out_1914130395499138379);
void car_inv_err_fun(double *nom_x, double *true_x, double *out_2456965719472217437);
void car_H_mod_fun(double *state, double *out_4779992119913662016);
void car_f_fun(double *state, double dt, double *out_2238700801222790874);
void car_F_fun(double *state, double dt, double *out_3188951629714217883);
void car_h_25(double *state, double *unused, double *out_2955682509261511320);
void car_H_25(double *state, double *unused, double *out_7551971679374724871);
void car_h_24(double *state, double *unused, double *out_4912051971117535801);
void car_H_24(double *state, double *unused, double *out_976399872783206770);
void car_h_30(double *state, double *unused, double *out_5047370505310505677);
void car_H_30(double *state, double *unused, double *out_635281337883108116);
void car_h_26(double *state, double *unused, double *out_1369527850446061802);
void car_H_26(double *state, double *unused, double *out_7153269075460770521);
void car_h_27(double *state, double *unused, double *out_3899192362163173755);
void car_H_27(double *state, double *unused, double *out_2810044649683533027);
void car_h_29(double *state, double *unused, double *out_595359365770451158);
void car_H_29(double *state, double *unused, double *out_4523407376553084060);
void car_h_28(double *state, double *unused, double *out_5177906949419811836);
void car_H_28(double *state, double *unused, double *out_8840937680086936982);
void car_h_31(double *state, double *unused, double *out_3815152717088839616);
void car_H_31(double *state, double *unused, double *out_7521325717497764443);
void car_predict(double *in_x, double *in_P, double *in_Q, double dt);
void car_set_mass(double x);
void car_set_rotational_inertia(double x);
void car_set_center_to_front(double x);
void car_set_center_to_rear(double x);
void car_set_stiffness_front(double x);
void car_set_stiffness_rear(double x);
}