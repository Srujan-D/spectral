#!/usr/bin/python3
from scipy.spatial import ConvexHull, distance
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
from mpl_toolkits.mplot3d import axes3d
from cmath import inf
from calendar import c
import time
from commonroad.scenario.scenario import Tag
from commonroad.scenario.scenario import Location
from commonroad.common.file_writer import OverwriteExistingFile
from commonroad.common.file_writer import CommonRoadFileWriter
from commonroad.prediction.prediction import TrajectoryPrediction
from commonroad.scenario.trajectory import State, Trajectory
from commonroad.scenario.obstacle import StaticObstacle, DynamicObstacle, ObstacleType
from commonroad.visualization.mp_renderer import MPRenderer
from commonroad.common.file_reader import CommonRoadFileReader
from IPython import display
import matplotlib.pyplot as plt
import os
from turtle import color
import numpy as np
from math import *
from math import pi, fmod

import math
from typing import List
import copy
import numpy as np
from sqlalchemy import true

# from trp_wrapper import find_traj

from cub_wrapper import find_traj

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

s_u_l = 50.0
s_l_l = 0.0

d_u_l = 8.0
d_l_l = -2.0

time_ = 7.0

num_of_knots = int(time_)*10+1

homotopy = 'yield'
cur_s = []
cur_l = []


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
			lanelets = self.scenario.lanelet_network.find_lanelet_by_position([
																			  X_ref[i, 0:2, 0]])
			if len(lanelets) > 0 and len(lanelets[0]) > 0:
				lanelet = self.scenario.lanelet_network.find_lanelet_by_id(
					lanelets[0][0])
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
	if (angle < 0.0):
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
	if fabs(rs - s_condition[0]) >= 1.0e-6:
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
	if fabs(rs - s_condition[0]) >= 1.0e-6:
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

	v = sqrt(one_minus_kappa_r_d * one_minus_kappa_r_d *
			 s_condition[1] * s_condition[1] + d_dot * d_dot)

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
	if fabs(rs - s_condition[0]) >= 1.0e-6:
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

	v = sqrt(one_minus_kappa_r_d * one_minus_kappa_r_d *
			 s_condition[1] * s_condition[1] + d_dot * d_dot)

	delta_theta_prime = one_minus_kappa_r_d / \
		cos_delta_theta * (kappa) - rkappa
	a = (s_condition[2] * one_minus_kappa_r_d / cos_delta_theta +
		 s_condition[1] * s_condition[1] / cos_delta_theta *
		 (d_condition[1] * delta_theta_prime - kappa_r_d_prime))
	return x, y, v, a, theta, kappa


def write_bounds(bounds, ego):
	corridor_file = open(
		"/home/srujan_d/RISS/code/btrapz/src/c_road_s1_3.txt", 'w')
	corridor_file.write("71 0.1\n\n")
	corridor_file.write(str(0.0) + " " + str(round(ego.velocity *
						cos(ego.orientation), 2)) + " " + str(0.0) + "\n\n")
	corridor_file.write(str(round(ego.position[1], 2)) + " " + str(
		ego.velocity*sin(ego.orientation)) + " " + str(0.0) + "\n\n")
	corridor_file.write(str(len(bounds)) + "\n\n")

	corridor_file.write("7.0 0.0 \n\n-2.0 2.0\n-30 30\n\n-0.7 0.7\n-10 10\n\n")

	for i in range(len(bounds)):
		for j in range(num_of_knots):
			if j % 10 == 9:
				corridor_file.write("\n")
			corridor_file.write(
				str(bounds[i][0][j][0]) + " " + str(bounds[i][0][j][1]) + " ")
		corridor_file.write("\n\n")
		for j in range(num_of_knots):
			if j % 10 == 9:
				corridor_file.write("\n")
			corridor_file.write(
				str(bounds[i][1][j][0]) + " " + str(bounds[i][1][j][1]) + " ")
		corridor_file.write("\n\n")
		# print("\n\n")

	corridor_file.write("0 10.0 0 10.0 0 10.0 0 10.0 0 10.0 0 10.0 0 10.0 0 10.0 0 10.0 0 10.0\n\
0 10.0 0 10.0 0 10.0 0 10.0 0 10.0 0 10.0 0 10.0 0 10.0 0 10.0 0 10.0\n\
0 10.0 0 10.0 0 10.0 0 10.0 0 10.0 0 10.0 0 10.0 0 10.0 0 10.0 0 10.0\n\
0 10.0 0 10.0 0 10.0 0 10.0 0 10.0 0 10.0 0 10.0 0 10.0 0 10.0 0 10.0\n\
0 10.0 0 10.0 0 10.0 0 10.0 0 10.0 0 10.0 0 10.0 0 10.0 0 10.0 0 10.0\n\
0 10.0 0 10.0 0 10.0 0 10.0 0 10.0 0 10.0 0 10.0 0 10.0 0 10.0 0 10.0\n\
0 10.0 0 10.0 0 10.0 0 10.0 0 10.0 0 10.0 0 10.0 0 10.0 0 10.0 0 10.0 0 10.0\n\n")

	corridor_file.write("-2.0 2.0 -2.0 2.0 -2.0 2.0 -2.0 2.0 -2.0 2.0 -2.0 2.0 -2.0 2.0 -2.0 2.0 -2.0 2.0 -2.0 2.0\n\
-2.0 2.0 -2.0 2.0 -2.0 2.0 -2.0 2.0 -2.0 2.0 -2.0 2.0 -2.0 2.0 -2.0 2.0 -2.0 2.0 -2.0 2.0\n\
-2.0 2.0 -2.0 2.0 -2.0 2.0 -2.0 2.0 -2.0 2.0 -2.0 2.0 -2.0 2.0 -2.0 2.0 -2.0 2.0 -2.0 2.0\n\
-2.0 2.0 -2.0 2.0 -2.0 2.0 -2.0 2.0 -2.0 2.0 -2.0 2.0 -2.0 2.0 -2.0 2.0 -2.0 2.0 -2.0 2.0\n\
-2.0 2.0 -2.0 2.0 -2.0 2.0 -2.0 2.0 -2.0 2.0 -2.0 2.0 -2.0 2.0 -2.0 2.0 -2.0 2.0 -2.0 2.0\n\
-2.0 2.0 -2.0 2.0 -2.0 2.0 -2.0 2.0 -2.0 2.0 -2.0 2.0 -2.0 2.0 -2.0 2.0 -2.0 2.0 -2.0 2.0\n\
-2.0 2.0 -2.0 2.0 -2.0 2.0 -2.0 2.0 -2.0 2.0 -2.0 2.0 -2.0 2.0 -2.0 2.0 -2.0 2.0 -2.0 2.0 -2.0 2.0\n\n")
	# print(bounds)

	write_reference_traj(corridor_file)
	# corridor_file.write("0.0 0.6 1.2 1.8 2.4 3.0 3.6 4.2 4.8 5.4\n\
	# 					6.0 6.6 7.2 7.8 8.4 9.0 9.6 10.2 10.8 11.4\n\
	# 					12.0 12.6 13.2 13.8 14.4 15.0 15.6 16.2 16.8 17.4\n\
	# 					18.0 18.6 19.2 19.8 20.4 21.0 21.6 22.2 22.8 23.4\n\
	# 					24.0 24.6 25.2 25.8 26.4 27.0 27.6 28.2 28.8 29.4\n\
	# 					30.0 30.6 31.2 31.8 32.4 33.0 33.6 34.2 34.8 35.4\n\
	# 					36.0 36.6 37.2 37.8 38.4 39.0 39.6 40.2 40.8 41.4\n\
	# 					42.0\n")

	# corridor_file.write("\n\n")

	# corridor_file.write("3.5 3.5 3.5 3.5 3.5 3.5 3.5 3.5 3.5\n\
	# 					3.5 3.5 3.41 3.31 3.22 3.12 3.03 2.94 2.84 2.75\n\
	# 					2.66 2.56 2.47 2.38 2.28 2.19 2.09 2.0 1.91 1.81\n\
	# 					1.72 1.62 1.53 1.44 1.34 1.25 1.16 1.06 0.97 0.88\n\
	# 					0.78 0.69 0.59 0.5 0.41 0.31 0.22 0.12 0.03 -0.06\n\
	# 					-0.16 -0.25 -0.25 -0.25 -0.25 -0.25 -0.25 -0.25 -0.25 -0.25\n\
	# 					-0.25 -0.25 -0.25 -0.25 -0.25 -0.25 -0.25 -0.25 -0.25 -0.25\n\
	# 					-0.25 -0.25\n\n")

	corridor_file.write(
		"0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0\n")

	corridor_file.write(
		"0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0\n")


x = []
y = []
z = []

ref_x = [-0.25, -0.25, -0.25, -0.22, -0.19, -0.17, -0.15, -0.12, -0.09, -0.03, -0.0, 0.7, 1.2, 1.223, 1.239, 1.259, 1.284, 1.315, 1.35, 1.392, 1.438, 1.489, 1.544, 1.604, 1.668, 1.735, 1.805, 1.877, 1.952, 2.027, 2.104, 2.181, 2.257, 2.333, 2.407, 2.479, 2.549, 2.615, 2.679, 2.738, 2.794, 2.844, 2.891, 2.932, 2.968, 3.0, 3.027, 3.049, 3.067, 3.081, 3.092, 3.1, 3.106, 3.11, 3.112, 3.114, 3.114, 3.115, 3.115, 3.115, 3.115, 3.115, 3.115, 3.115, 3.115, 3.115, 3.115, 3.115, 3.115, 3.115, 3.115, 3.115, 3.115, 3.115, 3.115, 3.115,
		 3.115, 3.115, 3.115, 3.115, 3.115, 3.115, 3.5, 3.5, 3.5, 3.5, 3.5, 3.5, 3.5, 3.5, 3.5, 3.5, 3.5, 3.41, 3.31, 3.22, 3.12, 3.03, 2.94, 2.84, 2.75, 2.66, 2.56, 2.47, 2.38, 2.28, 2.19, 2.09, 1.91, 1.72, 1.53, 1.34, 1.25, 1.16, 1.06, 0.97, 0.78, 0.59, 0.41, 0.22, 0.03, -0.06, -0.16, -0.25, -0.25, -0.25, -0.25, -0.25, -0.25, -0.25, -0.25, -0.25, -0.25, -0.25, -0.25, -0.25, -0.25, -0.25, -0.25, -0.25, -0.25, -0.25, -0.25, -0.25, -0.25, -0.25, -0.25, -0.25, -0.25, -0.25, -0.25, -0.25, -0.25, -0.25, -0.25, -0.25, -0.25, -0.25, -0.25, -0.25, -0.25, -0.25, -0.25, -0.25, -0.25, -0.25, -0.25, -0.25, -0.25, -0.25, -0.25, -0.25, -0.25, -0.25, -0.25]
ref_y = [0.0, 0.93, 1.855, 2.771, 3.673, 4.559, 5.426, 6.274, 7.1, 7.905, 8.689, 9.451, 10.193, 10.913, 11.612, 12.29, 12.947, 13.582, 14.197, 14.79, 15.363, 15.915, 16.447, 16.957, 17.447, 17.917, 18.365, 18.793, 19.201, 19.588, 19.955, 20.302, 20.628, 20.935, 21.223, 21.494, 21.751, 22.0, 22.246, 22.495, 22.753, 23.023, 23.309, 23.612, 23.934, 24.277, 24.639, 25.021, 25.424, 25.846, 26.288, 26.751, 27.233, 27.735, 28.257, 28.8, 29.362, 29.944, 30.547, 31.169, 31.811, 32.473, 33.155, 33.857, 34.578,
		35.319, 36.079, 36.857, 37.654, 38.469, 39.301, 40.15, 41.016, 41.899, 42.0, 42.6, 43.2, 43.8, 44.4, 45.0, 45.6, 46.2, 46.8, 47.4, 48.0, 48.6, 49.2, 49.8, 50.4, 51.0, 51.6, 52.2, 52.8, 53.4, 54.0, 54.6, 55.2, 55.8, 56.4, 57.0, 57.6, 58.2, 58.8, 59.4, 60.0, 60.6, 61.2, 61.8, 62.4, 63.0, 63.6, 64.2, 64.8, 65.4, 66.0, 66.6, 67.2, 67.8, 68.4, 69.0, 69.6, 70.2, 70.8, 71.4, 72.0, 72.6, 73.2, 73.8, 74.4, 75.0, 75.6, 76.2, 76.8, 77.4, 78.0, 78.6, 79.2, 79.8, 80.4, 81.0, 81.6, 82.2, 82.8, 83.4, 84.0, 84.6, 85.2, 85.8, 86.4, 87.0, 87.5, 88.0, 88.5, 89.0, 89.5, 90.0, 90.3, 90.6, 91.0, 91.3, 91.6, 92.0, 92.5, 93.0, 93.5, 93.6, 93.7, 93.8, 94.0, 94.3, 94.5, 95.0, 95.5, 96.0, 97.0]

def write_reference_traj(corridor_file):

	'''
	
	#using last local path as reference trajectory 

	input_txt = "/home/srujan_d/RISS/code/btrapz/src/s1_slt_3d_1.txt"
	f = open(input_txt)
	ref_y = []
	ref_x = []

	for line in f:
		line = line.strip("\n")
		line = line.split(" ")

		ref_y.append(float(line[1]))
		ref_x.append(float(line[2]))

	f.close

	res = 1
	# print("resolution: ", res)
	# print("before:", len(ref_y))
	# print(res)
	min_x = ref_y[0]
	for i in range(num_of_knots):
		x.append(round(ref_y[i*res]-min_x, 2))
		corridor_file.write(str(round(ref_y[i*res]-min_x, 2)) + ' ')

	corridor_file.write("\n\n")

	for i in range(num_of_knots):
		y.append(round(ref_x[i*res], 2))
		z.append(i)
		corridor_file.write(str(round(ref_x[i*res], 2)) + ' ')

	corridor_file.write("\n\n")

	return
	'''

	global ref_x, ref_y

	# res = int(len(route.reference_path)/(2*route.reference_path[-1][0]))
	# res = int(len(ref_y)/(2*ref_y[-1]))
	# res = 1
	res = int(len(ref_y)/(2*ref_y[-1])) if int(len(ref_y)/(2*ref_y[-1])) > 1 else 1
	print("resolution: ", res)
	# print("before:", len(ref_y))
	# print(res)
	min_x = ref_y[0]
	for i in range(num_of_knots):
		x.append(round(ref_y[i*res]-min_x, 2))
		corridor_file.write(str(round(ref_y[i*res]-min_x, 2)) + ' ')

	corridor_file.write("\n\n")

	for i in range(num_of_knots):
		y.append(round(ref_x[i*res], 2))
		z.append(i)
		corridor_file.write(str(round(ref_x[i*res], 2)) + ' ')

	corridor_file.write("\n\n")

	# ref_x = np.delete(ref_x, [i for i in range(res)], axis=0)
	# ref_y = np.delete(ref_y, [i for i in range(res)], axis=0)
	if len(ref_x) > num_of_knots:
		ref_x = ref_x[res:]
		ref_y = ref_y[res:]
	# print("after:", len(ref_x))
	# del route.reference_path[:num_of_knots]

	return

	'''
	using route planner as reference trajectory
	'''

	# res = int(len(route.reference_path)/(2*route.reference_path[-1][0]))
	# print("resolution: ", res)
	# print("before:",len(route.reference_path))
	# # print(res)
	# min_x = route.reference_path[0][0]
	# for i in range(num_of_knots):
	# 	x.append(round(route.reference_path[i*res][0]-min_x, 2))
	# 	corridor_file.write(str(round(route.reference_path[i*res][0]-min_x, 2))+ ' ')

	# corridor_file.write("\n\n")

	# for i in range(num_of_knots):
	# 	y.append(round(route.reference_path[i*res][1], 2))
	# 	z.append(i)
	# 	corridor_file.write(str(round(route.reference_path[i*res][1], 2))+ ' ')

	# corridor_file.write("\n\n")

	# route.reference_path = np.delete(route.reference_path, [i for i in range(res)], axis=0)
	# print("after:",len(route.reference_path))
	# del route.reference_path[:num_of_knots]

	return


class Ego:
	def __init__(self, centre, vel_s=0.0, vel_l=0.0, time=3.0):
		self.centre = centre

		self.forw_t_safe = time
		self.vel_s = vel_s
		self.vel_l = vel_l
		self.beg_l = self.centre[1]

		self.car = self.getCar()

	def getCar(self):
		if self.centre[2] != -1:
			l_car = 5
			l_ego = 5
			w_car = 2
			w_ego = 2

			l_safe = l_car/2 + l_ego/2
			w_safe = w_car/2 + w_ego/2

			# forw_t_safe = 2 # time interval forward = 3 sec
			# vel_s = 3 # const vel longitudinal s dir = 3 m/s
			# vel_l = 0.0 # const vel lateral l dir = 0 m/s

			forw_state = (self.centre[0]+self.vel_s*self.forw_t_safe, self.centre[1] +
						  self.vel_l*self.forw_t_safe, self.centre[2]+self.forw_t_safe)

			if self.vel_l >= 0:

				min_l = self.centre[1]-w_safe
				max_l = forw_state[1]+w_safe

				self.car = np.array([[self.centre[0]-l_safe, self.centre[1]-w_safe, self.centre[2]],
									 [self.centre[0]-l_safe, max_l, self.centre[2]],
									 [self.centre[0]+l_safe, self.centre[1] -
										 w_safe, self.centre[2]],
									 [self.centre[0]+l_safe, max_l, self.centre[2]],
									 [forw_state[0] - l_safe,
										 min_l, forw_state[2]],
									 [forw_state[0] - l_safe, forw_state[1] +
										 w_safe, forw_state[2]],
									 [forw_state[0] + l_safe,
										 min_l, forw_state[2]],
									 [forw_state[0] + l_safe, forw_state[1]+w_safe, forw_state[2]]])
				return self.car

			else:

				min_l = forw_state[1]-w_safe
				max_l = self.centre[1]+w_safe

				self.car = np.array([[self.centre[0]-l_safe, min_l, self.centre[2]],
									 [self.centre[0]-l_safe, self.centre[1] +
										 w_safe, self.centre[2]],
									 [self.centre[0]+l_safe, min_l, self.centre[2]],
									 [self.centre[0]+l_safe, self.centre[1] +
										 w_safe, self.centre[2]],
									 [forw_state[0] - l_safe, forw_state[1] -
										 w_safe, forw_state[2]],
									 [forw_state[0] - l_safe,
										 max_l, forw_state[2]],
									 [forw_state[0] + l_safe, forw_state[1] -
										 w_safe, forw_state[2]],
									 [forw_state[0] + l_safe, max_l, forw_state[2]]])
				return self.car

		else:
			self.car = np.array([[0, self.centre[1], 0], [0, self.centre[1]+0.01, 0], [30, self.centre[1], 0], [30, self.centre[1]+0.01, 0], [
								0, self.centre[1], 8], [0, self.centre[1]+0.01, 8], [30, self.centre[1], 8], [30, self.centre[1]+0.01, 8]])
			return self.car

	def getLBounds(self):
		return (min(self.car, key=lambda x: x[1])[1], max(self.car, key=lambda y: y[1])[1])


class Car:
	_lateral = []

	def get_lateral(self):
		return type(self)._lateral

	def set_lateral(self, ref, val):
		type(self)._lateral.append((self, ref, val))
		type(self)._lateral.sort(key=lambda x: x[2])
		# type(self)._lateral =  list(set(type(self)._lateral))

	lateral = property(get_lateral, set_lateral)

	# @property
	# def lateral(self):
	# 	return type(self)._lateral

	# @lateral.setter
	# def lateral(self,val):
	# 	type(self)._lateral.append(val)
	# 	type(self)._lateral.sort()

	def __init__(self, centre, vel_s=0.0, vel_l=0.0, time=3.0, ref=-1):
		self.centre = centre

		self.forw_t_safe = time
		self.vel_s = vel_s
		self.vel_l = vel_l
		self.beg_l = self.centre[1]

		self.ref = ref

		self.car = self.getCar()

		self.name = 'C' if self.centre == (5, 5, 3) else 'A'

		self.breaks = []

		# self.set_lateral(self.beg_l)

	def __lt__(self, a):
		return hash(self) < hash(a)

	def getCar(self):
		if self.centre[2] != -1:
			l_car = 5
			l_ego = 5
			w_car = 2
			w_ego = 2

			l_safe = l_car/3 + l_ego/3
			w_safe = w_car/3 + w_ego/3

			# forw_t_safe = 2 # time interval forward = 3 sec
			# vel_s = 3 # const vel longitudinal s dir = 3 m/s
			# vel_l = 0.0 # const vel lateral l dir = 0 m/s

			forw_state = (self.centre[0]+self.vel_s*self.forw_t_safe, self.centre[1] +
						  self.vel_l*self.forw_t_safe, self.centre[2]+self.forw_t_safe)

			if self.vel_l >= 0:

				min_l = self.centre[1]-w_safe
				max_l = forw_state[1]+w_safe
				# print(self, min_l, max_l)

				t = copy.deepcopy(Car._lateral)
				# print(t)

				for car in t:
					# print(car)
					if car[0].car[0][1] == min_l and car[0].car[1][1] == max_l:
						# print(car[0].car[0][1], min_l, car[0].car[1][1], max_l)
						# print("s9 -- same car")
						continue
					if car[0].car[1][1] < min_l:
						# print("s1 - car completely on right")
						continue
					if car[0].car[0][1] < min_l and car[0].car[1][1] == min_l:
						# print("s7 -- car just on right")
						continue
					if car[0].car[0][1] < min_l and car[0].car[1][1] > min_l:
						if car[0].car[1][1] < max_l:
							self.set_lateral(self.ref, car[0].car[1][1])
							car[0].set_lateral(car[0].ref, min_l)
							# print("s2")
							continue
						elif car[0].car[1][1] == max_l:
							# print("s8")
							continue
						else:
							car[0].set_lateral(car[0].ref, min_l)
							car[0].set_lateral(car[0].ref, max_l)
							# print("s6")
							continue
					if car[0].car[0][1] == min_l and car[0].car[1][1] < max_l:
						self.set_lateral(self.ref, car[0].car[1][1])
						# print("s10")
						continue
					if car[0].car[0][1] > min_l and car[0].car[1][1] < max_l:
						self.set_lateral(self.ref, car[0].car[0][1])
						self.set_lateral(self.ref, car[0].car[1][1])
						# print("s3")
						continue
					if car[0].car[0][1] > min_l and car[0].car[0][1] < max_l and car[0].car[1][1] > max_l:
						self.set_lateral(self.ref, car[0].car[0][1])
						car[0].set_lateral(car[0].ref, max_l)
						# print("s4")
						continue
					if car[0].car[0][1] > min_l and car[0].car[1][1] == max_l:
						self.set_lateral(self.ref, car[0].car[0][1])
						# print("s11")
						continue
					if car[0].car[0][1] > max_l:
						# print("s5")
						continue

				self.set_lateral(self.ref, min_l)
				self.set_lateral(self.ref, max_l)

				self.car = np.array([[self.centre[0]-l_safe, self.centre[1]-w_safe, self.centre[2]],
									 [self.centre[0]-l_safe, max_l, self.centre[2]],
									 [self.centre[0]+l_safe, self.centre[1] -
										 w_safe, self.centre[2]],
									 [self.centre[0]+l_safe, max_l, self.centre[2]],
									 [forw_state[0] - l_safe,
										 min_l, forw_state[2]],
									 [forw_state[0] - l_safe, forw_state[1] +
										 w_safe, forw_state[2]],
									 [forw_state[0] + l_safe,
										 min_l, forw_state[2]],
									 [forw_state[0] + l_safe, forw_state[1]+w_safe, forw_state[2]]])
				return self.car

			else:

				min_l = forw_state[1]-w_safe
				max_l = self.centre[1]+w_safe

				self.set_lateral(self.ref, min_l)
				self.set_lateral(self.ref, max_l)

				self.car = np.array([[self.centre[0]-l_safe, min_l, self.centre[2]],
									 [self.centre[0]-l_safe, self.centre[1] +
										 w_safe, self.centre[2]],
									 [self.centre[0]+l_safe, min_l, self.centre[2]],
									 [self.centre[0]+l_safe, self.centre[1] +
										 w_safe, self.centre[2]],
									 [forw_state[0] - l_safe, forw_state[1] -
										 w_safe, forw_state[2]],
									 [forw_state[0] - l_safe,
										 max_l, forw_state[2]],
									 [forw_state[0] + l_safe, forw_state[1] -
										 w_safe, forw_state[2]],
									 [forw_state[0] + l_safe, max_l, forw_state[2]]])
				return self.car

		else:
			self.car = np.array([[0, self.centre[1], 0], [0, self.centre[1]+0.01, 0], [30, self.centre[1], 0], [30, self.centre[1]+0.01, 0], [
								0, self.centre[1], 8], [0, self.centre[1]+0.01, 8], [30, self.centre[1], 8], [30, self.centre[1]+0.01, 8]])
			return self.car

	def getLBounds(self):
		return (min(self.car, key=lambda x: x[1])[1], max(self.car, key=lambda y: y[1])[1])


def delete_multiple_element(list_object, indices):
	indices = sorted(indices, reverse=True)
	for idx in indices:
		if idx < len(list_object):
			list_object.pop(idx)

	return list_object


def lineFromPoints(x1, y1, x2, y2):
	# print(x1,y1,x2,y2)
	a = y2 - y1
	b = x2 - x1
	c = a/b

	temp = []
	for i in range(num_of_knots):
		y = round(c*i/10-c*x1+y1, 2)
		temp.append(y)

	return temp


def get_bounds(lateral):
	s_bounds = []
	l_bounds = []

	# print(lateral)
	# for i in lateral:
	# 	print(i[0].car)

	# lateral = list(set([i for i in lateral]))
	# print(lateral)

	lateral = list(set(lateral))

	# print(lateral)
	# for i in lateral:
	# 	print(i[0].car)

	# print("afterw")

	lateral.sort(key=lambda x: x[2])
	# print("after sorting")
	# print(lateral)

	t = copy.deepcopy(lateral)
	t.sort(key=lambda x: x[2])
	# print("t sorting 0")
	# print(t)

	from collections import defaultdict
	cars = defaultdict(list)

	for car in t:
		cars[car[1]].append([car[0], car[2]])

	# print(cars)
	# print(len(cars))

	# print("after sorting")
	# for car in Car._lateral:
	# 	print(car[0].name, car[1])

	keys = list(cars.keys())

	flag = False
	if homotopy == 'yield':

		for car in (keys):
			for edge in range(len(cars[car])-1):
				if (cars[car][edge][1] > d_l_l and flag == False):
					cur_s = [[s_l_l, s_u_l] for i in range(num_of_knots)]
					cur_l = [[d_l_l, cars[car][edge][1]]
							 for i in range(num_of_knots)]
					s_bounds.append(cur_s)
					l_bounds.append(cur_l)
					flag = True
					# continue

				# if flag == True:
				# 	car = 0
				# 	flag = False

				# if (cars[car][edge][1] == cars[car][-1][1] and cars[car][edge][1] < d_u_l):
				# 	print("last------")
				# 	cur_s = [[s_l_l, s_u_l] for i in range(num_of_knots)]
				# 	cur_l = [[t[car][2], d_u_l] for i in range(num_of_knots)]
				# 	s_bounds.append(cur_s)
				# 	l_bounds.append(cur_l)
				# 	continue

				# if (car == len(t)-1 and car != max(t, key = lambda x: x[2])):
				# 	continue

				if (cars[car][edge][0].car[0][2] == 0):
					s_u = lineFromPoints(cars[car][edge][0].car[0][2], cars[car][edge][0].car[0]
										 [0], cars[car][edge][0].car[4][2], cars[car][edge][0].car[4][0])

					for i in range(num_of_knots):
						if i < cars[car][edge][0].car[0][2]*10:
							s_u[i] = s_u_l
						elif i > cars[car][edge][0].car[4][2]*10:
							s_u[i] = s_u_l

					cur_s = [[s_l_l, s_u[i]] for i in range(num_of_knots)]
					s_bounds.append(cur_s)

					try:
						# print("try1")
						cur_l = [[cars[car][edge][1], cars[car][edge+1][1]]
								 for i in range(num_of_knots)]
					except:
						# print("except1")
						cur_l = [[cars[car][edge][1], d_u_l]
								 for i in range(num_of_knots)]
					try:
						l_bounds.append(cur_l)
						if (cur_l[0][0] > cur_l[0][1] or cur_l[0][0] == cur_l[0][1]):
							# print("popping")
							l_bounds.pop()
							s_bounds.pop()
					except:
						# print("pass1")
						pass

				else:
					s_l = lineFromPoints(cars[car][edge][0].car[2][2], cars[car][edge][0].car[2]
										 [0], cars[car][edge][0].car[6][2], cars[car][edge][0].car[6][0])
					for i in range(num_of_knots):
						if i < cars[car][edge][0].car[0][2]*10:
							s_l[i] = s_l_l
						elif i > cars[car][edge][0].car[4][2]*10:
							s_l[i] = s_l_l

					cur_s = [[s_l[i], s_u_l] for i in range(num_of_knots)]
					s_bounds.append(cur_s)

					try:
						# print("tr2")
						cur_l = [[cars[car][edge][1], cars[car][edge+1][1]]
								 for i in range(num_of_knots)]
					except:
						# print("except2")
						cur_l = [[cars[car][edge][1], d_u_l]
								 for i in range(num_of_knots)]
					try:
						l_bounds.append(cur_l)
						if (cur_l[0][0] > cur_l[0][1] or cur_l[0][0] == cur_l[0][1]):
							# print("popping2")
							l_bounds.pop()
							s_bounds.pop()
					except:
						# print("pass2")
						pass

				cur_s = []
				cur_l = []

		max_l = max([cars[car][-1][1] for car in keys])
		if max_l < d_u_l:
			# print("last------")
			cur_s = [[s_l_l, s_u_l] for i in range(num_of_knots)]
			cur_l = [[max_l, d_u_l] for i in range(num_of_knots)]
			s_bounds.append(cur_s)
			l_bounds.append(cur_l)

		# print(len(s_bounds), len(l_bounds))
		# for i in range(len(l_bounds)):
		# 	print("i", i)
		# 	print(s_bounds[i])
		# 	print(l_bounds[i])
		# print(ds_bounds)
		# print(kappa)

		check = set()
		for l1 in range(len(l_bounds)-1):
			for l2 in range(l1+1, len(l_bounds)):
				if l_bounds[l1] == l_bounds[l2]:
					s_bounds[l1] = [[max(s_bounds[l1][i][0], s_bounds[l2][i][0]), min(
						s_bounds[l1][i][1], s_bounds[l2][i][1])] for i in range(len(s_bounds[l1]))]
					s_bounds[l2] = [[max(s_bounds[l1][i][0], s_bounds[l2][i][0]), min(
						s_bounds[l1][i][1], s_bounds[l2][i][1])] for i in range(len(s_bounds[l1]))]
					check.add(l2)

		# print(check)

		l_bounds = delete_multiple_element(l_bounds, check)
		s_bounds = delete_multiple_element(s_bounds, check)

		for l in range(len(l_bounds)-1):
			if l_bounds[l][0][1] != l_bounds[l+1][0][0]:
				# print("disconnect")
				cur_s = [[s_l_l, s_u_l] for i in range(num_of_knots)]
				cur_l = [[l_bounds[l][0][1], l_bounds[l+1][0][0]]
						 for i in range(num_of_knots)]
				s_bounds.insert(l+1, cur_s)
				l_bounds.insert(l+1, cur_l)

		# for i in check:
		# 	l_bounds.pop(i)
		# 	s_bounds.pop(i)

		# 	check.remove(i)

		# for i in range(len(l_bounds)):
		# 	print("i", i)
		# 	print(s_bounds[i])
		# 	print(l_bounds[i])

		bounds = []
		for i in range(len(l_bounds)):
			bounds.append([s_bounds[i], l_bounds[i]])

		Car._lateral = []

		return bounds


plt.rcParams.update({'figure.max_open_warning': 0})


# import functions to read xml file and visualize commonroad objects


# file_path = "/home/srujan_d/RISS/code/btrapz/src/commonroad/xml/ZAM_Tutorial-1_1_T-1.xml"
file_path = "/home/srujan_d/RISS/code/btrapz/src/commonroad/xml/s4.xml"
scenario, planning_problem_set = CommonRoadFileReader(file_path).open()

planning_problem = list(planning_problem_set.planning_problem_dict.values())[0]

# print(planning_problem._initial_state)

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

dynamic_obs_init_state = State(position=np.array([15.0, 0.5]),
							   velocity=5.1,
							   orientation=0.197,
							   time_step=0)

state_list = []
for i in range(1, 141):
	if i < 31:
		new_pos = np.array([dynamic_obs_init_state.position[0] + scenario.dt *
						   i*5.0, dynamic_obs_init_state.position[1] + scenario.dt*i*1])
		new_state = State(position=new_pos, velocity=5.1,
						  orientation=0.197, time_step=i)
		state_list.append(new_state)
	else:
		new_pos = np.array([state_list[29].position[0] +
						   scenario.dt*(i-31)*6.0, state_list[i-2].position[1]])
		new_state = State(position=new_pos, velocity=6.0,
						  orientation=0.0, time_step=i)
		state_list.append(new_state)

dynamic_obs_traj = Trajectory(1, state_list)
dynamic_obs_shape = Rectangle(width=2.0, length=4.5)
dynamic_obs_pred = TrajectoryPrediction(dynamic_obs_traj, dynamic_obs_shape)

dynamic_obs_id = scenario.generate_object_id()
dynamic_obs_type = ObstacleType.CAR
dynamic_obs = DynamicObstacle(dynamic_obs_id,
							  dynamic_obs_type,
							  dynamic_obs_shape,
							  dynamic_obs_init_state,
							  dynamic_obs_pred)

scenario.add_objects(dynamic_obs)


dynamic_obs_init_state = State(position=np.array([0.0, 3.5]),
							   velocity=4.5,
							   orientation=0.0,
							   time_step=0)

state_list = []
for i in range(1, 141):
	new_pos = np.array([dynamic_obs_init_state.position[0] + scenario.dt *
					   i*4.5, dynamic_obs_init_state.position[1] + scenario.dt*i*0.0])
	new_state = State(position=new_pos, velocity=5.0,
					  orientation=0.0, time_step=i)
	state_list.append(new_state)

dynamic_obs_traj = Trajectory(1, state_list)
dynamic_obs_shape = Rectangle(width=2.0, length=4.5)
dynamic_obs_pred = TrajectoryPrediction(dynamic_obs_traj, dynamic_obs_shape)

dynamic_obs_id = scenario.generate_object_id()
dynamic_obs_type = ObstacleType.CAR
dynamic_obs = DynamicObstacle(dynamic_obs_id,
							  dynamic_obs_type,
							  dynamic_obs_shape,
							  dynamic_obs_init_state,
							  dynamic_obs_pred)

scenario.add_objects(dynamic_obs)


def run_ego(ego, global_t, curr, succ, curr_state):
	input_txt = "/home/srujan_d/RISS/code/btrapz/src/s1_cub_3d_3.txt"
	f = open(input_txt)
	ego_y = []
	ego_x = []
	ego_dx = []
	ego_dy = []
	angle = []

	for line in f:
		line = line.strip("\n")
		line = line.split(" ")

		ego_y.append(float(line[1]))
		ego_x.append(float(line[2]))
		if (float(line[3]) > 5.0):
			ego_dy.append(float(line[3]))
		else:
			ego_dy.append(5.0)
		ego_dx.append(float(line[4]))

	f.close

	for i in range(len(ego_x)):
		try:
			angle.append(
				round((np.arctan2([ego_x[i + 1] - ego_x[i]], [ego_y[i + 1] - ego_y[i]])
					   )[0], 2)
			)
		except:
			angle.append(
				round((np.arctan2([ego_x[i] - ego_x[i - 1]], [ego_y[i] - ego_y[i - 1]])
					   )[0], 2)
			)
		if isnan(angle[i]):
			angle[i] = 0.0

	def scale_range(input, min, max):
		input += -(np.min(input))
		input /= np.max(input) / (max - min)
		input += min
		return input

	# ego_x = scale_range(ego_x, -0.25, 3.5)
	# ego_x = scale_range(ego_x, 0.0, max(ego_x))

	# print("init state = ", ego.state_at_time(64).position[0], ego.state_at_time(64).position[1])
	# ego_vehicle_init_state = State(position=np.array([ego.state_at_time(curr-1).position[0]+ scenario.dt*ego.state_at_time(curr-1).velocity*cos(ego.state_at_time(curr-1).orientation), ego.state_at_time(curr-1).position[1]]),
	# 								velocity = ego.state_at_time(curr-1).velocity,
	# 								orientation = ego.state_at_time(curr-1).orientation,
	# 								time_step = curr)
	
	# print("ego_y init ", ego_y[0])
	
	
	ego_vehicle_init_state = State(position=np.array([ego_y[0], ego_x[0]]),
								   velocity=(ego_dx[0]**2 + ego_dy[0]**2)**0.5,
								   orientation=angle[0],
								   time_step=curr)

	# ego_vehicle_init_state = curr_state

	if succ == True:
		print("found new traj!!!")
		state_list = [ego_vehicle_init_state]
		for i in range(1, len(ego_x)):
			new_pos = np.array(
				[ego_y[i]+curr_state.position[0], ego_x[i]])
			# new_pos = np.array([ego_vehicle_init_state.position[0] + scenario.dt*i*ego_dy[i], ego_x[i]])
			new_state = State(position=new_pos, velocity=(
				ego_dx[i]**2 + ego_dy[i]**2)**0.5, orientation=angle[i], time_step=i+curr)
			state_list.append(new_state)
			# print("new pos: ", ego_y[i], ego_vehicle_init_state.position[0])
			# if i < 41:
			# 	new_pos = np.array([ego_y[i]+ego_vehicle_init_state.position[0], ego_x[i]])
			# 	# new_pos = np.array([ego_vehicle_init_state.position[0] + scenario.dt*i*ego_dy[i], ego_x[i]])
			# 	new_state = State(position = new_pos, velocity = (ego_dx[i]**2 + ego_dy[i]**2)**0.5, orientation = angle[i], time_step = i+curr)
			# 	state_list.append(new_state)
			# else:
			# 	# new_pos = np.array([ego_y[i]+ego_vehicle_init_state.position[0], ego_x[i]])
			# 	try:
			# 		new_pos = np.array([state_list[39].position[0] + scenario.dt*(i-39)*ego_dy[i], ego_x[i]])
			# 		new_state = State(position = new_pos, velocity = (ego_dx[i]**2 + ego_dy[i]**2)**0.5, orientation = angle[i], time_step = i+curr)
			# 		state_list.append(new_state)
			# 	except:
			# 		print("bugggg----------------", len(state_list), 39, len(ego_dy), len(ego_x), i)
		ego_vehicle_traj = Trajectory(curr, state_list)
		ego_vehicle_shape = Rectangle(width=2.0, length=4.5)
		# ego_prediction = TrajectoryPrediction(ego_vehicle_traj, ego_vehicle_shape)

		ego.prediction = TrajectoryPrediction(
			ego_vehicle_traj, ego_vehicle_shape)

	# for i in range(71):
	# 	print("next ego:\n", ego.state_at_time(i+71))

	# ego_vehicle_id = scenario.generate_object_id()
	# ego_vehicle_type = ObstacleType.CAR

	# ego_vehicle = DynamicObstacle(ego_vehicle_id,
	# 							ego_vehicle_type,
	# 							ego.obstacle_shape,
	# 							ego_vehicle_init_state,
	# 							ego_prediction)

	# scenario.add_objects(ego_vehicle)
	# print("run ego: ", ego.state_at_time(62))
	else:
		print("continuing old traj!!!", curr, len(ego_x))
		state_list = []
		for i in range(curr, len(ego_x)):
			new_pos = np.array([ego_y[i], ego_x[i]])
			# new_pos = np.array([ego_vehicle_init_state.position[0] + scenario.dt*i*ego_dy[i], ego_x[i]])
			new_state = State(position=new_pos, velocity=(
				ego_dx[i]**2 + ego_dy[i]**2)**0.5, orientation=angle[i], time_step=i)
			state_list.append(new_state)
			# print("new pos: ", ego_y[i], ego_vehicle_init_state.position[0])
			# if i < 41:
			# 	new_pos = np.array([ego_y[i]+ego_vehicle_init_state.position[0], ego_x[i]])
			# 	# new_pos = np.array([ego_vehicle_init_state.position[0] + scenario.dt*i*ego_dy[i], ego_x[i]])
			# 	new_state = State(position = new_pos, velocity = (ego_dx[i]**2 + ego_dy[i]**2)**0.5, orientation = angle[i], time_step = i+curr)
			# 	state_list.append(new_state)
			# else:
			# 	# new_pos = np.array([ego_y[i]+ego_vehicle_init_state.position[0], ego_x[i]])
			# 	try:
			# 		new_pos = np.array([state_list[39].position[0] + scenario.dt*(i-39)*ego_dy[i], ego_x[i]])
			# 		new_state = State(position = new_pos, velocity = (ego_dx[i]**2 + ego_dy[i]**2)**0.5, orientation = angle[i], time_step = i+curr)
			# 		state_list.append(new_state)
			# 	except:
			# 		print("bugggg----------------", len(state_list), 39, len(ego_dy), len(ego_x), i)
		ego_vehicle_traj = Trajectory(curr, state_list)
		ego_vehicle_shape = Rectangle(width=2.0, length=4.5)
		# ego_prediction = TrajectoryPrediction(ego_vehicle_traj, ego_vehicle_shape)

	# 	ego.prediction = TrajectoryPrediction(ego_vehicle_traj, ego_vehicle_shape)
	return ego


input_txt = "/home/srujan_d/RISS/code/btrapz/src/s1_cub_3d_3.txt"
f = open(input_txt)
ego_y = []
ego_x = []
ego_dx = []
ego_dy = []
angle = []

for line in f:
	line = line.strip("\n")
	line = line.split(" ")

	ego_y.append(float(line[1]))
	ego_x.append(float(line[2]))
	ego_dy.append(float(line[3]))
	ego_dx.append(float(line[4]))

f.close

# print("start x", ego_x)
# print("start y", ego_y)

for i in range(len(ego_x)):
	try:
		angle.append(
			round((np.arctan2([ego_x[i + 1] - ego_x[i]], [ego_y[i + 1] - ego_y[i]])
				   )[0], 2)
		)
	except:
		angle.append(
			round((np.arctan2([ego_x[i] - ego_x[i - 1]], [ego_y[i] - ego_y[i - 1]])
				   )[0], 2)
		)


def scale_range(input, min, max):
	input += -(np.min(input))
	input /= np.max(input) / (max - min)
	input += min
	return input


ego_x = scale_range(ego_x, 0.0, 3.5)
# ego_x = scale_range(ego_x, 0.0, max(ego_x))

ego_vehicle_init_state = State(position=np.array([ego_y[0]+2.5, ego_x[0]]),
							   velocity=(ego_dx[0]**2 + ego_dy[0]**2)**0.5,
							   orientation=0.0,
							   time_step=0)

state_list = []
for i in range(1, 71):
	new_pos = np.array([ego_y[i]+2.5, ego_x[i]])
	new_state = State(position=new_pos, velocity=(
		ego_dx[i]**2 + ego_dy[i]**2)**0.5, orientation=angle[i], time_step=i)
	state_list.append(new_state)

ego_vehicle_traj = Trajectory(1, state_list)
ego_vehicle_shape = Rectangle(width=2.0, length=4.5)
ego_vehicle_pred = TrajectoryPrediction(ego_vehicle_traj, ego_vehicle_shape)

ego_vehicle_id = scenario.generate_object_id()
ego_vehicle_type = ObstacleType.CAR
ego_vehicle = DynamicObstacle(ego_vehicle_id,
							  ego_vehicle_type,
							  ego_vehicle_shape,
							  ego_vehicle_init_state,
							  ego_vehicle_pred)

# scenario.add_objects(ego_vehicle)


author = 'Srujan Deolasee'
affiliation = 'CMU'
source = ''

# write new scenario
# fw = CommonRoadFileWriter(scenario, planning_problem_set, author, affiliation, source, tags)

# filename = "/home/srujan_d/RISS/code/btrapz/src/commonroad/xml/s5.xml"
# fw.write_to_file(filename, OverwriteExistingFile.ALWAYS)


# file_path = "/home/srujan_d/RISS/code/btrapz/src/commonroad/xml/s4.xml"

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


# instantiate a route planner with the scenario and the planning problem

route_planner = RoutePlanner(
	scenario, planning_problem, backend=RoutePlanner.Backend.NETWORKX_REVERSED)


# print(route_planner.goal_region._state_list[0].position)
# plan routes, and save the routes in a route candidate holder
candidate_holder = route_planner.plan_routes()

# option 1: retrieve all routes
list_routes, num_route_candidates = candidate_holder.retrieve_all_routes()
print(f"Number of route candidates: {num_route_candidates}")
# here we retrieve the first route in the list, this is equivalent to: route = list_routes[0]
# route = candidate_holder.retrieve_first_route()

# option 2: retrieve the best route by orientation metric
route = candidate_holder.retrieve_best_route_by_orientation()


# print coordinates of the vertices of the reference path
print("\nCoordinates [x, y]:")
# for waypoint in route.reference_path:
# 	print(waypoint)
print(len(route.reference_path))

visualize_route(route, draw_route_lanelets=True,
				draw_reference_path=True, size_x=16)

# plt.ion()

k = 0
j = 0


def plot_car(centre, obj, fig, ax):
	# ax.scatter(centre[0], centre[1], centre[2] if centre[2] != -1 else 0, color = 'g')

	for cube, color in zip([obj], ['b']):
		# print(cube)
		hull = ConvexHull(cube)
		for s in hull.simplices:
			temp = cube[s].tolist()
			tri = Poly3DCollection([temp])
			color = 'b' if centre[2] != -1 else 'r'
			alpha = 0.5 if centre[2] != -1 else 0.2
			tri.set_color(color)
			tri.set_alpha(alpha)
			tri.set_edgecolor('none')
			ax.add_collection3d(tri)
			edges = []
			if distance.euclidean(cube[s[0]], cube[s[1]]) < distance.euclidean(cube[s[1]], cube[s[2]]):
				edges.append((s[0], s[1]))
				if distance.euclidean(cube[s[1]], cube[s[2]]) < distance.euclidean(cube[s[2]], cube[s[0]]):
					edges.append((s[1], s[2]))
				else:
					edges.append((s[2], s[0]))
			else:
				edges.append((s[1], s[2]))
				if distance.euclidean(cube[s[0]], cube[s[1]]) < distance.euclidean(cube[s[2]], cube[s[0]]):
					edges.append((s[0], s[1]))
				else:
					edges.append((s[2], s[0]))
			for v0, v1 in edges:
				ax.plot(xs=cube[[v0, v1], 0], ys=cube[[v0, v1],
						1], zs=cube[[v0, v1], 2], color='black')


def call_plot_car(centres, cars, ref):
	fig1 = plt.figure()
	ax = fig1.add_subplot(111, projection="3d")
	ax.set_xlabel("s")
	ax.set_ylabel("l")
	ax.set_zlabel("t")
	plt.margins(0)

	# print(x,y,z)

	for centre, car in zip(centres, cars):

		plot_car(centre, car, fig1, ax)

	for point in range(len(x)):
		ax.scatter(x[point], y[point], z[point]/10,
				   color='r', label='reference trajectory')
	x.clear()
	y.clear()
	z.clear()
	plt.savefig(
		"/home/srujan_d/RISS/code/btrapz/src/commonroad/trajectory/with_orient/replan_cont/slt/traj_s1_"+str(j)+".png")
	# plt.show()
	plt.close()

# fig = plt.figure()
# ax = fig.add_subplot(111, projection="3d")
# ax.set_xlabel("s")
# ax.set_ylabel("l")
# ax.set_zlabel("t")
# plt.margins(0)


replan = False
kplan = 0
# print(scenario.dynamic_obstacles)
for i in range(0, 141):
	# if kplan > 50:
	# 	print("enough")
	# 	print("breaking")
	plt.figure(figsize=(25, 10))
	rnd = MPRenderer()
	scenario.draw(rnd, draw_params={'time_begin': i})

	car1 = scenario.dynamic_obstacles[0]
	car2 = scenario.dynamic_obstacles[1]

	# if k == 0 or i < 71:
	# 	# print("old ego", ego_vehicle.state_at_time(i))
	# 	ego_vehicle.draw(rnd, draw_params={'time_begin': i, 'dynamic_obstacle': {
	# 		'vehicle_shape': {'occupancy': {'shape': {'rectangle': {
	# 			'facecolor': 'g'}}}}}})
	# else:
	# 	if replan == True:
	# 		new_ego = run_ego(ego_vehicle)
	# 		replan = False
	# 	# print("new ego", new_ego.state_at_time(i))
	# 	new_ego.draw(rnd, draw_params={'time_begin': i, 'dynamic_obstacle': {
	# 		'vehicle_shape': {'occupancy': {'shape': {'rectangle': {
	# 			'facecolor': 'g'}}}}}})

	if replan == True and j > 30:
		kplan += 1
		if succ == True:
			k = 0
		print("replanning...")
		ego_vehicle = run_ego(ego_vehicle, i, k, succ, curr_state)

	ego_vehicle.draw(rnd, draw_params={'time_begin': k, 'dynamic_obstacle': {
		'vehicle_shape': {'occupancy': {'shape': {'rectangle': {
			'facecolor': 'g'}}}}}})

	print("i", i)
	
	ego = ego_vehicle.state_at_time(k)
	if k == 0:
		print("ego state ", ego)
	# if i%65 == 64 and k < 1:
	# print("replanning...")
	seconds1 = time.time()
	state2 = car2.state_at_time(k)

	centre2 = (abs(state2.position[0]-ego.position[0]), state2.position[1],
			   0 if state2.position[0]-ego.position[0] > 0 else abs(state2.position[0]-ego.position[0])/5.0)
	c2 = Car(centre2, vel_l=0, vel_s=4.5, time=4.0, ref=11)
	obj2 = c2.getCar()

	state1 = car1.state_at_time(k)
	centre1 = (state1.position[0]-ego.position[0], state1.position[1], 0 if state1.position[0] -
			   ego.position[0] > 0 else abs(state1.position[0]-ego.position[0])/5.0)
	c1 = Car(centre1, vel_s=6.0, time=4.0, ref=10)
	obj1 = c1.getCar()

	# call_plot_car([centre1, centre2], [obj1, obj2], route.reference_path)

	# plot_car(centre1, obj1, fig, ax)
	# plot_car(centre2, obj2, fig, ax)
	# plt.savefig("/home/srujan_d/RISS/code/btrapz/src/commonroad/trajectory/with_orient/slt/traj_s1_"+str(j)+".png")
	# plt.show()

	# print(centre1, centre2)
	# print("getting bounds...")
	bounds = get_bounds(Car._lateral)

	write_bounds(bounds, ego)

	# call_plot_car([centre1,centre2], [obj1, obj2], route.reference_path)

	
	curr_state = ego_vehicle.state_at_time(k)
	if k > 2:
		print("finding trajectory...")
		succ = find_traj()
		replan = True
	else:
		print("k:", k)
		replan = False

	seconds2 = time.time()
	print("time passed: ", seconds2-seconds1)
	
	if len(ref_x) == num_of_knots:
		replan = False
	# new_ego = run_ego(ego_vehicle)

	# i = 0
	k += 1

	# if k == 2:
	# 	break

	j += 1

	planning_problem_set.draw(rnd)
	rnd.render()
	# plt.draw()
	# plt.show()
	# plt.pause(0.0001)
	# plt.margins(0, 0)
	plt.savefig(
		"/home/srujan_d/RISS/code/btrapz/src/commonroad/trajectory/with_orient/diff_replan_cont/traj_s1_"+str(j)+".png")


# plt.draw()
# plt.show()
