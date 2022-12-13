#ifndef BTRAPZ_PY_CPP_H
#define BTRAPZ_PY_CPP_H

#include "solve_3d.h"

struct Params {
    double s_acc_weight;
    double s_jerk_weight;
    double l_acc_weight;
    double l_jerk_weight;

	double weight_s_ref;
    double weight_ds_ref;
    double weight_l_ref;
    double weight_dl_ref;

    double weight_end_s;
    double weight_end_l;

    int iteration;
};

struct Metrics {
    double s_avg_acc;
    double l_avg_acc;

    double s_max_acc;
    double l_max_acc;

    double a_cost;
};

#endif