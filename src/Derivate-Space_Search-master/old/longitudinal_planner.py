import numpy as np
import heapq
import math
# import flowstar
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.animation import FuncAnimation
from matplotlib.axis import Axis

from vehicle_state import LongitudinalState
from search_node import TranverseState, Node
from derivate_space_class import DerivateSpace
from predict_obstacles import PredictObstacles

vehicle_length = 5.0
vehicle_width = 2.2
vehicle_dis = 3.0
lane_width = 3.5

lane_y_lower_bound = -2
lane_y_middle = lane_y_lower_bound + lane_width
lane_y_upper_bound = lane_y_lower_bound + 2 * lane_width

initial_state_ego = [0, 0, 0, 10]
initial_state_obstacle = [10.0, 0, 0, 10]
current_state_ego = initial_state_ego
current_state_obstacle = initial_state_obstacle
planning_dt = 0.1
planning_horizon = 5.0
simulate_states_ego = []
simulate_states_obstacle = []

longitudinal_yield = 1
def GetSearchResult(final_node : Node):
    pass

def HeuristicSearch(derivate_space, current_state_ego, dds_initial, st_graph, behavior_state):
    
    nums_of_knots = len(st_graph)
    if behavior_state == longitudinal_yield:
        safe_distance = 2.0
        sf_ref = st_graph[nums_of_knots-1][1] - safe_distance
    
    # ds_final_ref = (st_graph[nums_of_knots-1] - st_graph[nums_of_knots-2]) / planning_dt
    ds_ref = 12.0
    dds_ref = 0.0
    s_ref = ds_ref * planning_horizon * 2

    # Heuristic Search on ddS-T Graph
    start_node = Node(LongitudinalState(current_state_ego[0], current_state_ego[3],\
        dds_initial, 0.0), TranverseState(0, 0, 0, 0))
    goal_state = LongitudinalState(s_ref, ds_ref, dds_ref, planning_horizon)
    goal_node = Node(goal_state, TranverseState(0, 0, 0, 0))
    final_node = goal_node
    # start_node.CalculateNodeCost(goal_state, ds_ref)
    open_set, closed_set = {}, {}
    priority_queue = []
    heapq.heappush(priority_queue, (start_node.GetNodeCost(), start_node))

    while(open_set):
        current_node_cost, current_node_index = heapq.heappop(priority_queue)
        if current_node_index in open_set:
            current_node = open_set[current_node_index]
            open_set.pop(current_node_index)
        else:
            continue
        if current_node_index in closed_set:
            continue
        else:
            closed_set[current_node_index] = current_node

        if current_node.GetTime() == planning_horizon:
            final_node = current_node
            break
        
        for dds_next_node in derivate_space.a_candidate_set:
            tranverse_ind = derivate_space.GetTranverseInd(current_node.longitudinal_state.dds, dds_next_node)
            tranverse_dds = derivate_space.a_tranverse_pairs[tranverse_ind[0]][tranverse_ind[1]]
            tranverse_ds = derivate_space.ds_tranverse_pairs[tranverse_ind[0]][tranverse_ind[1]]
            tranverse_s = derivate_space.s_tranverse_pairs[tranverse_ind[0]][tranverse_ind[1]]
            tranverse_cost = derivate_space.transition_cost_pairs[tranverse_ind[0]][tranverse_ind[1]]
            tranverse_state = TranverseState(tranverse_s, tranverse_ds, tranverse_dds, tranverse_cost)

            next_node = current_node.GenerateNextNode(tranverse_state)
            if next_node.IsCollisionFree(st_graph):
                next_node_cost = next_node.GetNodeCost()
                next_node_index = next_node.GetNodeIndex()
                if next_node_index in closed_set:
                    continue
                else:
                    heapq.heappush(priority_queue, (next_node.GetNodeCost(), next_node))
    for i in priority_queue:
        print(i[1].longitudinal_state.s, i[1].longitudinal_state.ds, i[1].longitudinal_state.dds)

    s_ref, ds_ref, dds_ref = [s_ref, s_ref], [ds_ref, ds_ref], [dds_ref, dds_ref]
    # s_ref, ds_ref, dds_ref = GetSearchResult(final_node)

    return s_ref, ds_ref, dds_ref

def PlanningSimulate(simulate_steps, simulate_delta_t):

    obstacle_acc = -2.0
    derivate_space = DerivateSpace()
    derivate_space.GenerateDdsCandidateSet()
    dds_ego = 0.0
    last_state_obstacle = current_state_obstacle
    for i in range(0, simulate_steps):
        # simulate ego vehicle
        last_state_ego = current_state_ego
        st_graph = PredictObstacles(last_state_obstacle, obstacle_acc)
        behavior_state = longitudinal_yield
        s_ref, ds_ref, dds_ref = HeuristicSearch(derivate_space, current_state_ego, dds_ego, st_graph, behavior_state)
        current_state_ego[0] = s_ref[1]
        current_state_ego[1] = last_state_ego[1]
        current_state_ego[2] = 0
        current_state_ego[3] = ds_ref[1]
        dds_geo = dds_ref[1]
        # current_state_ego[0] = last_state_ego[0] + last_state_ego[3] * math.cos(last_state_ego[2]) * simulate_delta_t
        # current_state_ego[1] = last_state_ego[1] + last_state_ego[3] * math.sin(last_state_ego[2]) * simulate_delta_t
        # current_state_ego[2] = 0
        # current_state_ego[3] = last_state_ego[3]
        
        # simulate obstacles
        
        # last_state_obstacle = current_state_obstacle
        current_state_obstacle[0] = last_state_obstacle[0] + last_state_obstacle[3] * math.cos(last_state_obstacle[2]) * simulate_delta_t
        current_state_obstacle[1] = last_state_obstacle[1] + last_state_obstacle[3] * math.sin(last_state_obstacle[2]) * simulate_delta_t
        current_state_obstacle[2] = 0
        current_state_obstacle[3] = last_state_obstacle[3] - obstacle_acc * simulate_delta_t
        
        last_state_obstacle = current_state_obstacle
        # to avoid yinyong
        simulate_states_ego.append([current_state_ego[0], current_state_ego[1], current_state_ego[2], current_state_ego[3]])
        simulate_states_obstacle.append([current_state_obstacle[0], current_state_obstacle[1], current_state_obstacle[2], current_state_obstacle[3]])

fig = plt.figure()
ax = fig.add_subplot()
# fig, ax = plt.subplots()

def ComputeArea(current_state_ego, range = 70):
    range_num_x = math.floor(current_state_ego[0] / range)
    axis_x_min = range_num_x * range - 5.0
    axis_x_max = (range_num_x + 1) * range + 5.0
    axis_y_min = lane_y_lower_bound - 1.0
    axis_y_max = 6.0
    return [axis_x_min, axis_x_max, axis_y_min, axis_y_max]

def DrawLanes(area):
    axis_x_min, axis_x_max, axis_y_min, axis_y_max = area[0], area[1], area[2], area[3]

    lane_x_lower = np.linspace(axis_x_min, axis_x_max, 1000)
    lane_y_lower = np.linspace(lane_y_lower_bound, lane_y_lower_bound, 1000)
    lane_middle_x = np.linspace(axis_x_min, axis_x_max, 1000)
    lane_middle_y = np.linspace(lane_y_middle, lane_y_middle, 1000)
    lane_x_upper = np.linspace(axis_x_min, axis_x_max, 1000)
    lane_y_upper = np.linspace(lane_y_upper_bound, lane_y_upper_bound, 1000)

    # plot three lanes
    plt.plot(lane_x_lower, lane_y_lower, color = 'k', linewidth = 2.0)
    plt.plot(lane_middle_x, lane_middle_y, linestyle = '--', color = 'k', linewidth = 2.0)
    plt.plot(lane_x_upper, lane_y_upper, color = 'k', linewidth = 2.0)

    # set scales of x-axis and y-axis
    unit_x = 10
    unit_y = 2
    x_tick_min = math.ceil(axis_x_min / unit_x) * unit_x
    x_tick_max = math.floor(axis_x_max / unit_x) * unit_x + unit_x
    x_ticks = np.arange(x_tick_min, x_tick_max, unit_x)
    y_ticks = np.arange(lane_y_lower_bound, lane_y_upper_bound + unit_y, unit_y)
    plt.xticks(x_ticks)
    plt.yticks(y_ticks)

    plt.axis(area)
    plt.gca().set_aspect(1)
    plt.title('Simulations of the ego vehicle and surrounding vehicles')
    plt.xlabel('x(m)')
    plt.ylabel('y(m)')


def Update(frame_number):
    current_state_ego = simulate_states_ego[frame_number]
    current_state_obstacle = simulate_states_obstacle[frame_number]
    area = ComputeArea(current_state_ego)
    DrawLanes(area)

    # ax.patches = []
    # ax.cla()
    # ego_patch = []
    # obstacle_patch = []
    ax.add_patch(patches.Rectangle((current_state_ego[0] - vehicle_length / 2, current_state_ego[1] - vehicle_width / 2), vehicle_length, vehicle_width, 
                    angle = current_state_ego[2], edgecolor = 'green', linewidth = 1.0, facecolor = 'green', fill = True))
    ax.add_patch(patches.Rectangle((current_state_obstacle[0] - vehicle_length / 2, current_state_obstacle[1] - vehicle_width / 2), vehicle_length, vehicle_width, 
                    angle = current_state_obstacle[2], edgecolor = 'blue', facecolor = 'blue', fill = True))
    
animation_delta_t = 0.2
animation_frames = 60
PlanningSimulate(animation_frames + 5, animation_delta_t)
animation = FuncAnimation(fig, Update, frames = animation_frames, interval = animation_delta_t * 1000)
animation.save("planning_video.mp4", writer = 'ffmpeg', fps = 20)

plt.show()
