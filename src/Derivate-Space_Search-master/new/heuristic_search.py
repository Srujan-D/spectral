from math import floor
from math import pow
import math
import numpy as np
import matplotlib.pyplot as plt

dp_cost_max = 1.0e7
safe_distance = 15.0
collision_distance = 2.0


def CalculateCost(start_index, end_index, state_last, state_sample, s_min, s_max, vel_ref):
    dt = 0.1
    jerk = (state_sample[2] - state_last[2]) / ((end_index - start_index) * dt)
    s = state_last[0]
    vel = state_last[1]
    acc = state_last[2]
    acc_cost = 0.0
    vel_cost = 0.0
    obstacle_cost = 0.0
    obstacle_cost_max = 2.0 * 1.0e4
    for i in range(start_index + 1, end_index + 1):
        acc_last = acc
        vel_last = vel
        s_last = s
        acc = acc_last + jerk * dt
        vel = vel_last + 1/2 * (acc_last + acc) * dt
        s = s_last + vel_last * dt + 1/3 * acc_last * dt * dt + 1/6 * acc * dt * dt

        ego_velocity_lowerbound = 0.0
        ego_velocity_upperbound = 20.0
        if vel < ego_velocity_lowerbound or vel > ego_velocity_upperbound:
            return dp_cost_max
        if s > s_max[i] - collision_distance or s < s_min[i]:
            return dp_cost_max

        acc_cost += acc * acc * dt
        vel_cost += (vel - vel_ref) * (vel - vel_ref) * dt
        if s > s_max[i] - safe_distance:
            delta_collision = s_max[i] - collision_distance - s + 1.0
            obstacle_cost += obstacle_cost_max * dt / \
                (delta_collision * delta_collision * delta_collision)
        if s < s_min[i] + 0.15 * safe_distance:
            delta_collision = s - s_min[i] - collision_distance + 1.0
            obstacle_cost += obstacle_cost_max * dt / \
                (delta_collision * delta_collision * delta_collision)

    vel_weight = 30.0
    acc_weight = 100.0
    dp_cost = acc_weight * acc_cost + vel_weight * vel_cost + obstacle_cost

    # print("dp_cost is: ", dp_cost)
    return dp_cost


def HeuristicSearch(step, s_initial, s_min, s_max, num_of_knots, vel_ref, s_ref):

    a_min = -5.0
    a_max = 1.5
    a_step = 0.25
    sample_a_num = floor((a_max - a_min) / a_step) + 1
    sample_a = np.zeros(sample_a_num, dtype=float)

    for i in range(sample_a_num):
        sample_a[i] = a_min + i * a_step

    # Calculate cols_num time terminals
    cols_num = 7
    delta_col = np.zeros(cols_num, dtype=int)
    if (num_of_knots - 1) % cols_num == 0:
        delta_col[0] = (num_of_knots - 1) / cols_num
    else:
        delta_col[0] = floor((num_of_knots - 1) / (cols_num - 1))

    if delta_col[0] == 0:
        return False
    for i in range(1, cols_num):
        if i != cols_num - 1:
            delta_col[i] = delta_col[0]
        else:
            delta_col[i] = (num_of_knots - 1) - i * delta_col[0]

    dt = 0.1
    s_search = np.zeros(cols_num + 1, dtype=float)
    v_search = np.zeros(cols_num + 1, dtype=float)
    a_search = np.zeros(cols_num + 1, dtype=float)
    s_search[0] = s_initial[0]
    v_search[0] = s_initial[1]
    a_search[0] = s_initial[2]
    s_last = s_initial[0]
    v_last = s_initial[1]
    a_last = s_initial[2]
    state_last = [s_last, v_last, a_last]
    start_index = 0
    end_index = 0

    # 1. Dynamic programming in control space
    for i in range(1, cols_num + 1):
        s_optimal = s_last
        v_optimal = v_last
        a_optimal = a_last
        cost_optimal = dp_cost_max

        unit_t = delta_col[i - 1] * dt
        start_index = end_index
        end_index = start_index + delta_col[i - 1]

        for j in range(sample_a_num):
            a_sample = sample_a[j]
            v_sample = v_last + 1/2 * (a_last + a_sample) * unit_t
            s_sample = s_last + v_last * unit_t + 1/3 * a_last * unit_t * unit_t + \
                1/6 * a_sample * unit_t * unit_t
            state_sample = [s_sample, v_sample, a_sample]
            cost_sample = CalculateCost(
                start_index, end_index, state_last, state_sample, s_min, s_max, vel_ref)

            if cost_sample < cost_optimal:
                cost_optimal = cost_sample
                s_optimal = s_sample
                v_optimal = v_sample
                a_optimal = a_sample
        s_last = s_optimal
        v_last = v_optimal
        a_last = a_optimal
        s_search[i] = s_optimal
        v_search[i] = v_optimal
        a_search[i] = a_optimal
        state_last = [s_last, v_last, a_last]

    # 2. Calculate s_ref by s_search

    v_ref = np.zeros(num_of_knots, dtype=float)
    a_ref = np.zeros(num_of_knots, dtype=float)
    sub_shift = 0

    for i in range(0, cols_num):
        s_ref[sub_shift] = s_search[i]
        v_ref[sub_shift] = v_search[i]
        a_ref[sub_shift] = a_search[i]
        unit_t = delta_col[i] * dt
        jerk = (a_search[i+1]-a_search[i])/unit_t

        for j in range(1, delta_col[i] + 1):
            a_ref[sub_shift+j] = a_ref[sub_shift+j-1] + jerk * dt
            v_ref[sub_shift+j] = v_ref[sub_shift+j-1] + \
                (a_ref[sub_shift+j-1] + a_ref[sub_shift+j])/2 * dt
            s_ref[sub_shift+j] = s_ref[sub_shift+j-1] + v_ref[sub_shift+j-1] * dt + \
                1/3 * a_ref[sub_shift+j-1] * dt * dt + \
                1/6 * a_ref[sub_shift+j] * dt * dt

        sub_shift = sub_shift + delta_col[i]

    plt.clf()
    plot_t = np.arange(0, num_of_knots, 1)
    plt.plot(plot_t, a_ref, 'r')
    plt.savefig('visualize/at_search' + str(step) + '.jpg')

    plt.clf()
    plot_t = np.arange(0, num_of_knots, 1)
    plt.plot(plot_t, v_ref, 'r')
    plt.savefig('visualize/vt_search' + str(step) + '.jpg')

    # n_coef = 6
    # n_order = 5
    # factorial = np.zeros(n_order + 1, dtype = float)
    # factorial[0] = 1.0
    # for i in range(1, n_order + 1):
    #     factorial[i] = factorial[i-1] * i

    # b_coe = np.zeros(6, dtype = float)
    # for i in range(n_order + 1):
    #     b_coe[i] = factorial[n_order] / (factorial[i] * factorial[n_order-i])

    # s_ref[0] = s_initial[0]
    # for i in range(cols_num):
    #     c = np.zeros(n_coef, dtype = float)
    #     c[0] = s_search[i]
    #     c[1] = c[0] + v_search[i] / n_order
    #     c[2] = 2 * c[1] - c[0] + a_search[i] / (n_order * (n_order - 1))
    #     c[5] = s_search[i+1]
    #     c[4] = c[5] - v_search[i+1] / n_order
    #     c[3] = 2 * c[4] - c[5] + a_search[i+1] / (n_order * (n_order - 1))

    #     for j in range(1, delta_col[i] + 1):
    #         for k in range(n_order + 1):
    #             s_ref[sub_shift+j] += c[k] * b_coe[k] * pow(j /delta_col[i] , k) * pow(1 - (j /delta_col[i]), n_order - k)

    #     sub_shift = sub_shift + delta_col[i]

    return True
