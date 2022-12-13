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



# generate path of the file to be opened
file_path = "/home/srujan_d/RISS/code/btrapz/src/commonroad/xml/ZAM_Tutorial-1_1_T-1.xml"

# read in the scenario and planning problem set
scenario, planning_problem_set = CommonRoadFileReader(file_path).open()

static_obstacle_id = scenario.generate_object_id()
static_obstacle_type = ObstacleType.CONSTRUCTION_ZONE
static_obstacle_shape = Rectangle(width = 2.5, length = 10.0)
static_obstacle_initial_state = State(position = np.array([45.0, 0.0]), orientation = 0.0, time_step = 0)

static_obstacle = StaticObstacle(static_obstacle_id, static_obstacle_type, static_obstacle_shape, static_obstacle_initial_state)

scenario.add_objects(static_obstacle)

dynamic_obs_init_state = State(position=np.array([30.0, 0.5]),
                                velocity = 3.1,
                                orientation = 0.26,
                                time_step = 0)

state_list = []
for i in range(1, 41):
    new_pos = np.array([dynamic_obs_init_state.position[0] + scenario.dt*i*3.0, dynamic_obs_init_state.position[1] + scenario.dt*i*0.8])
    new_state = State(position = new_pos, velocity = 3.1, orientation = 0.26, time_step = i)
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


dynamic_obs_init_state = State(position=np.array([10.0, 3.5]),
                                velocity = 5.0,
                                orientation = 0.0,
                                time_step = 0)

state_list = []
for i in range(1, 41):
    new_pos = np.array([dynamic_obs_init_state.position[0] + scenario.dt*i*5.0, dynamic_obs_init_state.position[1] + scenario.dt*i*0.0])
    new_state = State(position = new_pos, velocity = 5.0, orientation = 0.0, time_step = i)
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

dynamic_obs_init_state = State(position=np.array([10.0, 3.5]),
                                velocity = 5.0,
                                orientation = 0.0,
                                time_step = 0)

state_list = []
for i in range(1, 41):
    new_pos = np.array([dynamic_obs_init_state.position[0] + scenario.dt*i*5.0, dynamic_obs_init_state.position[1] + scenario.dt*i*0.0])
    new_state = State(position = new_pos, velocity = 5.0, orientation = 0.0, time_step = i)
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


for i in range(0, 41):
    plt.figure(figsize=(25, 10))
    rnd = MPRenderer()
    scenario.draw(rnd, draw_params={'time_begin': i})
    planning_problem_set.draw(rnd)
    rnd.render()
    plt.show()