#!/usr/bin/python3
import numpy as np
from math import *
from math import pi, fmod

import math
from typing import List
import copy
import numpy as np

from commonroad.scenario.scenario import Scenario
from commonroad.scenario.trajectory import State
from commonroad.planning.planning_problem import PlanningProblem
from commonroad.geometry.shape import Rectangle

from commonroad_route_planner.route_planner import Route, RoutePlanner
from commonroad_route_planner.utility.visualization import obtain_plot_limits_from_reference_path, visualize_route
# from SMP.route_planner.route_planner.route import Route
# from SMP.route_planner.route_planner.route_planner import RoutePlanner
# from SMP.route_planner.route_planner.utils_visualization import get_plot_limits_from_reference_path, draw_route

FOLLOW_LANE = 0
FOLLOW_VEHICLE = 1
SHIFT_TO_LEFT = 2
SHIFT_TO_RIGHT = 3
STOP = 4

class Params:
    def __init__(self):
        self.planning_horizon = 70
        self.safe_distance = 2.5
        self.dt = scenario.dt

class BehaviourPlanner:
    def __init__(self, scenario: Scenario, planning_problem: PlanningProblem):
        self.scenario = scenario
        self.planning_problem = copy.copy(planning_problem)
        self.current_scenario = scenario
        self.list_routes = None

    def set_current_scenario(self, current_scenario: Scenario):
        self.current_scenario = current_scenario

    def get_routes(self, current_state: State) -> List[Route]:
        if not hasattr(current_state, 'slip_angle'):
            current_state.slip_angle = 0.0
        if not hasattr(current_state, 'steering_angle'):
            current_state.steering_angle = 0.0
        if not hasattr(current_state, 'yaw_rate'):
            current_state.yaw_rate = 0.0
        self.planning_problem.initial_state = current_state
        state = self.planning_problem.goal.state_list[0]
        if hasattr(state, 'position'):
            goal_position = copy.copy(state.position)
            if hasattr(goal_position, 'center'):
                new_goal_position = Rectangle(
                    0.1, 0.1, copy.copy(goal_position.center), copy.copy(goal_position.orientation))
                self.planning_problem.goal.state_list[0].position = new_goal_position

        # instantiate a route planner
        route_planner = RoutePlanner(
            self.current_scenario, self.planning_problem, backend=RoutePlanner.Backend.NETWORKX)

        # plan routes, and save the found routes in a route candidate holder
        candidate_holder = route_planner.plan_routes()

        # retrieving routes
        # option 1: retrieve all routes
        list_routes, num_routes = candidate_holder.retrieve_all_routes()
        if num_routes == 0:
            if hasattr(state, 'position'):
                self.planning_problem.goal.state_list[0].position = goal_position
            return self.list_routes
        self.list_routes = list_routes

        # for route in list_routes:
        #     plot_limits = get_plot_limits_from_reference_path(route)
        #     # option 2: plot limits from lanelets in the route
        #     # plot_limits = get_plot_limits_from_routes(route)
        #
        #     # determine the figure size for better visualization
        #     size_x = 6
        #     ratio_x_y = (plot_limits[1] - plot_limits[0]) / (plot_limits[3] - plot_limits[2])
        #     fig = plt.figure(figsize=(size_x, size_x / ratio_x_y))
        #     fig.gca().axis('equal')
        #
        #     draw_route(route, draw_route_lanelets=True, draw_reference_path=True, plot_limits=plot_limits)
        #     plt.show()
        #
        #     pass

        if hasattr(state, 'position'):
            self.planning_problem.goal.state_list[0].position = goal_position

        return list_routes

    def compute_reference(self, reference_path, x0, n, last: bool):
        X_ref = np.zeros((n + 1, 4, 1))
        X_ref[:, 0, 0] = reference_path[:n + 1, 0]
        X_ref[:, 1, 0] = reference_path[:n + 1, 1]

        Q = np.tile(
            np.array(
                [
                    [0.0, 0.0, 0.0, 0.0],
                    [0.0, 0.0, 0.0, 0.0],
                    [0.0, 0.0, 0.0, 0.0],
                    [0.0, 0.0, 0.0, 0.0],
                ]
            ),
            n + 1
        ).T.reshape(n + 1, 4, 4)
        Q[-1] = np.diag([200, 200, 0, 0])
        for i in range(len(X_ref)):
            lanelets = self.scenario.lanelet_network.find_lanelet_by_position([X_ref[i, 0:2, 0]])
            if len(lanelets) > 0 and len(lanelets[0]) > 0:
                lanelet = self.scenario.lanelet_network.find_lanelet_by_id(lanelets[0][0])
                if lanelet.adj_left is None and lanelet.adj_right is None:
                    Q[i, 0, 0] = 200
                    Q[i, 1, 1] = 200
                    # if i > 0:
                    #     X_ref[i, 3, 0] = math.atan2(
                    #         X_ref[i, 1, 0] - X_ref[i - 1, 1, 0], X_ref[i, 0, 0] - X_ref[i - 1, 0, 0])
                    #     Q[i, 3, 3] = 200000

        if last:
            Q[-1] = np.diag([2000, 2000, 0, 0])
            if not hasattr(self.planning_problem.goal.state_list[0], 'position'):
                Q[-1] = np.diag([0, 0, 0, 0])
            if hasattr(self.planning_problem.goal.state_list[0], 'velocity'):
                X_ref[-1, 2, 0] = \
                    (self.planning_problem.goal.state_list[0].velocity.start +
                     self.planning_problem.goal.state_list[0].velocity.end) / 2.0
                Q[-1, 2, 2] = 2000
            if hasattr(self.planning_problem.goal.state_list[0], 'orientation'):
                X_ref[-1, 3, 0] =\
                    (self.planning_problem.goal.state_list[0].orientation.start +
                     self.planning_problem.goal.state_list[0].orientation.end) / 2.0
                while abs(X_ref[-1, 3, 0] - x0[3]) > np.pi:
                    if X_ref[-1, 3, 0] > x0[3]:
                        X_ref[-1, 3, 0] -= 2 * np.pi
                    else:
                        X_ref[-1, 3, 0] += 2 * np.pi
                Q[-1, 3, 3] = 200000

        return X_ref, Q

def NormalizeAngle(angle):
    # if theta > pi:
    #     return theta - pi
    # else if theta < -pi:
    #     return theta + pi
    # return theta
    angle = fmod(angle + pi, 2.0 * pi)
    if(angle < 0.0):
        angle = angle + 2.0 * pi
    return angle - pi
    
def cartesian_to_frenet1D(rs, rx, ry, rtheta, x, y):
    s_condition = np.zeros(1)
    d_condition = np.zeros(1)
    
    dx = x - rx
    dy = y - ry
    
    if abs(dx) < 0.01 or abs(dy) < 0.01:
        return rs, 0.0
    cos_theta_r = cos(rtheta)
    sin_theta_r = sin(rtheta)
    
    cross_rd_nd = cos_theta_r * dy - sin_theta_r * dx

    # if np.isnan(dx) or np.isnan(dy):
    #     # print("x, y: ", x, y)
    #     # print("rx, ry: ", rx, ry)
    #     return nan, nan
    # try:
    #     dis = sqrt(dx * dx + dy * dy)
    # except:
    #     return nan, nan
    dis = sqrt(dx * dx + dy * dy)
    if dis < 0.01:
        d_condition[0] = 0.0
    else:
        d_condition[0] = copysign(dis, cross_rd_nd)    
    
    s_condition[0] = rs
    
    return s_condition, d_condition

def frenet_to_cartesian1D(rs, rx, ry, rtheta, s_condition, d_condition):
    if fabs(rs - s_condition[0])>= 1.0e-6:
        print("The reference point s and s_condition[0] don't match")
        
    cos_theta_r = cos(rtheta)
    sin_theta_r = sin(rtheta)
    
    x = rx - sin_theta_r * d_condition[0]
    y = ry + cos_theta_r * d_condition[0]    

    return x, y

def cartesian_to_frenet2D(rs, rx, ry, rtheta, rkappa, x, y, v, theta):
    s_condition = np.zeros(2)
    d_condition = np.zeros(2)
    
    dx = x - rx
    dy = y - ry
    
    cos_theta_r = cos(rtheta)
    sin_theta_r = sin(rtheta)
    
    cross_rd_nd = cos_theta_r * dy - sin_theta_r * dx
    dis = sqrt(dx * dx + dy * dy)
    if dis < 0.01:
        d_condition[0] = 0.0
    else:
        d_condition[0] = copysign(dis, cross_rd_nd)
    
    delta_theta = theta - rtheta
    tan_delta_theta = tan(delta_theta)
    cos_delta_theta = cos(delta_theta)
    
    one_minus_kappa_r_d = 1 - rkappa * d_condition[0]
    d_condition[1] = one_minus_kappa_r_d * tan_delta_theta
    
    
    s_condition[0] = rs
    s_condition[1] = v * cos_delta_theta / one_minus_kappa_r_d

    return s_condition, d_condition

def frenet_to_cartesian2D(rs, rx, ry, rtheta, rkappa, s_condition, d_condition):
    if fabs(rs - s_condition[0])>= 1.0e-6:
        print("The reference point s and s_condition[0] don't match")
        
    cos_theta_r = cos(rtheta)
    sin_theta_r = sin(rtheta)
    
    x = rx - sin_theta_r * d_condition[0]
    y = ry + cos_theta_r * d_condition[0]

    one_minus_kappa_r_d = 1 - rkappa * d_condition[0]
    tan_delta_theta = d_condition[1] / one_minus_kappa_r_d
    delta_theta = atan2(d_condition[1], one_minus_kappa_r_d)
    cos_delta_theta = cos(delta_theta)
    
    theta = NormalizeAngle(delta_theta + rtheta)    
    
    d_dot = d_condition[1] * s_condition[1]
    
    v = sqrt(one_minus_kappa_r_d * one_minus_kappa_r_d * s_condition[1] * s_condition[1] + d_dot * d_dot)   

    return x, y, v, theta

def cartesian_to_frenet3D(rs, rx, ry, rtheta, rkappa, rdkappa, x, y, v, a, theta, kappa):
    s_condition = np.zeros(3)
    d_condition = np.zeros(3)
    
    dx = x - rx
    dy = y - ry
    
    cos_theta_r = cos(rtheta)
    sin_theta_r = sin(rtheta)
    
    cross_rd_nd = cos_theta_r * dy - sin_theta_r * dx
    dis = sqrt(dx * dx + dy * dy)
    if dis < 0.01:
        d_condition[0] = 0.0
    else:
        d_condition[0] = copysign(dis, cross_rd_nd)
    
    delta_theta = theta - rtheta
    tan_delta_theta = tan(delta_theta)
    cos_delta_theta = cos(delta_theta)
    
    one_minus_kappa_r_d = 1 - rkappa * d_condition[0]
    d_condition[1] = one_minus_kappa_r_d * tan_delta_theta
    
    kappa_r_d_prime = rdkappa * d_condition[0] + rkappa * d_condition[1]
    
    d_condition[2] = (-kappa_r_d_prime * tan_delta_theta + 
      one_minus_kappa_r_d / cos_delta_theta / cos_delta_theta *
          (kappa * one_minus_kappa_r_d / cos_delta_theta - rkappa))
    
    s_condition[0] = rs
    s_condition[1] = v * cos_delta_theta / one_minus_kappa_r_d
    
    delta_theta_prime = one_minus_kappa_r_d / cos_delta_theta * kappa - rkappa
    s_condition[2] = ((a * cos_delta_theta -
                       s_condition[1] * s_condition[1] *
                       (d_condition[1] * delta_theta_prime - kappa_r_d_prime)) /
                          one_minus_kappa_r_d)
    return s_condition, d_condition

def frenet_to_cartesian3D(rs, rx, ry, rtheta, rkappa, rdkappa, s_condition, d_condition):
    if fabs(rs - s_condition[0])>= 1.0e-6:
        print("The reference point s and s_condition[0] don't match")
        
    cos_theta_r = cos(rtheta)
    sin_theta_r = sin(rtheta)
    
    x = rx - sin_theta_r * d_condition[0]
    y = ry + cos_theta_r * d_condition[0]

    one_minus_kappa_r_d = 1 - rkappa * d_condition[0]
    tan_delta_theta = d_condition[1] / one_minus_kappa_r_d
    delta_theta = atan2(d_condition[1], one_minus_kappa_r_d)
    cos_delta_theta = cos(delta_theta)
    
    theta = NormalizeAngle(delta_theta + rtheta)
    kappa_r_d_prime = rdkappa * d_condition[0] + rkappa * d_condition[1]
        
    kappa = ((((d_condition[2] + kappa_r_d_prime * tan_delta_theta) *
                 cos_delta_theta * cos_delta_theta) /
                    (one_minus_kappa_r_d) +
                rkappa) *
               cos_delta_theta / (one_minus_kappa_r_d))
    
    
    d_dot = d_condition[1] * s_condition[1]
    
    v = sqrt(one_minus_kappa_r_d * one_minus_kappa_r_d * s_condition[1] * s_condition[1] + d_dot * d_dot)
    
    delta_theta_prime = one_minus_kappa_r_d / cos_delta_theta * (kappa) - rkappa     
    a = (s_condition[2] * one_minus_kappa_r_d / cos_delta_theta +
           s_condition[1] * s_condition[1] / cos_delta_theta *
               (d_condition[1] * delta_theta_prime - kappa_r_d_prime))
    return x, y, v, a, theta, kappa

import os
import matplotlib.pyplot as plt
import numpy as np
from IPython import display

plt.rcParams.update({'figure.max_open_warning': 0})


# import functions to read xml file and visualize commonroad objects
from commonroad.common.file_reader import CommonRoadFileReader
from commonroad.visualization.mp_renderer import MPRenderer

from commonroad.geometry.shape import Rectangle
from commonroad.scenario.obstacle import StaticObstacle, DynamicObstacle, ObstacleType
from commonroad.scenario.trajectory import State, Trajectory
from commonroad.prediction.prediction import TrajectoryPrediction



# file_path = "/home/srujan_d/RISS/code/btrapz/src/commonroad/xml/ZAM_Tutorial-1_1_T-1.xml"
file_path = "/home/srujan_d/RISS/code/btrapz/src/commonroad/xml/s3.xml"
scenario, planning_problem_set = CommonRoadFileReader(file_path).open()

planning_problem = list(planning_problem_set.planning_problem_dict.values())[0]

# plot the scenario and the planning problem set
renderer = MPRenderer(figsize=(12, 12))

scenario.draw(renderer)
planning_problem.draw(renderer)

renderer.render()
# plt.margins(0,0)
plt.show()


# static_obstacle_id = scenario.generate_object_id()
# static_obstacle_type = ObstacleType.CONSTRUCTION_ZONE
# static_obstacle_shape = Rectangle(width = 2.5, length = 10.0)
# static_obstacle_initial_state = State(position = np.array([45.0, 0.0]), orientation = 0.0, time_step = 0)

# static_obstacle = StaticObstacle(static_obstacle_id, static_obstacle_type, static_obstacle_shape, static_obstacle_initial_state)

# scenario.add_objects(static_obstacle)

# dynamic_obs_init_state = State(position=np.array([30.0, 0.5]),
#                                 velocity = 4.08,
#                                 orientation = 0.197,
#                                 time_step = 0)

input_txt = "/home/srujan_d/RISS/code/btrapz/src/s1_slt_3d_2.txt"
f = open(input_txt)
ego_y = []
ego_x = []
ego_dx = []
ego_dy = []

for line in f:
    line = line.strip("\n")
    line = line.split(" ")

    ego_y.append(float(line[1]))
    ego_x.append(float(line[2]))
    ego_dy.append(float(line[3]))
    ego_dx.append(float(line[4]))

f.close

dynamic_obs_init_state = State(position=np.array([ego_y[0], ego_x[0]]),
                                velocity = (ego_dx[0]**2 + ego_dy[0]**2)**0.5,
                                orientation = 0.0,
                                time_step = 0)

state_list = []
for i in range(1, 71):
    new_pos = np.array([ego_y[i], ego_x[i]])
    new_state = State(position = new_pos, velocity = (ego_dx[i]**2 + ego_dy[i]**2)**0.5, orientation = 0.0, time_step = i)
    state_list.append(new_state)

dynamic_obs_traj = Trajectory(1, state_list)
dynamic_obs_shape = Rectangle(width = 2.0, length = 4.5)
dynamic_obs_pred = TrajectoryPrediction(dynamic_obs_traj, dynamic_obs_shape)

dynamic_obs_id = scenario.generate_object_id()
dynamic_obs_type = ObstacleType.CAR
dynamic_obs = DynamicObstacle(dynamic_obs_id,
                              dynamic_obs_type,
                              dynamic_obs_shape,
                              dynamic_obs_init_state,
                              dynamic_obs_pred)

scenario.add_objects(dynamic_obs)


# dynamic_obs_init_state = State(position=np.array([-5.0, 3.5]),
#                                 velocity = 5.0,
#                                 orientation = 0.0,
#                                 time_step = 0)

# state_list = []
# for i in range(1, 41):
#     new_pos = np.array([dynamic_obs_init_state.position[0] + scenario.dt*i*5.0, dynamic_obs_init_state.position[1] + scenario.dt*i*0.0])
#     new_state = State(position = new_pos, velocity = 5.0, orientation = 0.0, time_step = i)
#     state_list.append(new_state)

# dynamic_obs_traj = Trajectory(1, state_list)
# dynamic_obs_shape = Rectangle(width = 2.0, length = 4.5)
# dynamic_obs_pred = TrajectoryPrediction(dynamic_obs_traj, dynamic_obs_shape)

# dynamic_obs_id = scenario.generate_object_id()
# dynamic_obs_type = ObstacleType.CAR
# dynamic_obs = DynamicObstacle(dynamic_obs_id,
#                               dynamic_obs_type,
#                               dynamic_obs_shape,
#                               dynamic_obs_init_state,
#                               dynamic_obs_pred)

# scenario.add_objects(dynamic_obs)



# import necessary classes from different modules
# from commonroad.common.file_writer import CommonRoadFileWriter
# from commonroad.common.file_writer import OverwriteExistingFile
# from commonroad.scenario.scenario import Location
# from commonroad.scenario.scenario import Tag

# author = 'Srujan Deolasee'
# affiliation = 'CMU'
# source = ''
# tags = {Tag.CRITICAL, Tag.INTERSTATE}

# # write new scenario
# fw = CommonRoadFileWriter(scenario, planning_problem_set, author, affiliation, source, tags)

# filename = "ZAM_Tutorial-1_2_T-1.xml"
# fw.write_to_file(filename, OverwriteExistingFile.ALWAYS)


# file_path = "ZAM_Tutorial-1_2_T-1.xml"

# scenario, planning_problem_set = CommonRoadFileReader(file_path).open()


# behaviour_planner = BehaviourPlanner(scenario, planning_problem)

# state_current_ego = planning_problem.initial_state
# state_current_ego.steering_angle = 0.0
# current_ego_acceleration = 0.0
# obstacles = []
# obstacle_ids = []

# params = Params()

# end_time_step = planning_problem.goal.state_list[0].time_step.end
# frenet_coordinates_ego = np.zeros([3, 2], dtype = float)
# planned_ego_s = np.zeros(params.planning_horizon + 1, dtype = float)
# planned_ego_ds = np.zeros(params.planning_horizon + 1, dtype = float)
# planned_ego_dds = np.zeros(params.planning_horizon + 1, dtype = float)
# path_kappa = 0.0
# # print("end_time_step: ", end_time_step)
# state_list = [state_current_ego]
# s_initial = np.zeros(3)
# d_initial = np.zeros(3)
# for step in range(end_time_step):
#     behaviour_planner.set_current_scenario(scenario)
#     # retrieve the current state of the ego vehicle
#     X0 = np.array(
#         [
#             [state_current_ego.position[0]],
#             [state_current_ego.position[1]],
#             [state_current_ego.velocity],
#             [state_current_ego.orientation],
#         ]
#     )
#     # retrieve the reference path from behaviour planner
#     reference_path = behaviour_planner.get_routes(state_current_ego)[0].reference_path



# # instantiate a route planner with the scenario and the planning problem
# route_planner = RoutePlanner(scenario, planning_problem, backend=RoutePlanner.Backend.NETWORKX)
# print(route_planner.goal_region._state_list[0].position)
# # plan routes, and save the routes in a route candidate holder
# candidate_holder = route_planner.plan_routes()

# # option 1: retrieve all routes
# list_routes, num_route_candidates = candidate_holder.retrieve_all_routes()
# print(f"Number of route candidates: {num_route_candidates}")
# # here we retrieve the first route in the list, this is equivalent to: route = list_routes[0]
# route = candidate_holder.retrieve_first_route()

# # option 2: retrieve the best route by orientation metric
# # route = candidate_holder.retrieve_best_route_by_orientation()

# # print coordinates of the vertices of the reference path
# print("\nCoordinates [x, y]:")
# print(route.reference_path)
# print(len(route.reference_path))

# visualize_route(route, draw_route_lanelets=True, draw_reference_path=True, size_x=16)

# plt.ion()

for i in range(0, 71):
    plt.figure(figsize=(25, 10))
    rnd = MPRenderer()
    scenario.draw(rnd, draw_params={'time_begin': i})
    planning_problem_set.draw(rnd)
    rnd.render()
    # plt.draw()
    # plt.show()
    # plt.pause(0.0001)
    # plt.margins(0, 0)
    plt.show()
    

# plt.draw()
# plt.show()

    