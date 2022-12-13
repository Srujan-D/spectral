import math

def PredictObstacles(state_obstacle, obstacle_acc):
    prediction_horizon = 5.0
    prediction_dt = 0.1
    prediction_step_nums = math.ceil(prediction_horizon / prediction_dt) + 1
    st_graph = []
    current_state_obstacle = state_obstacle
    vehicle_length = 5.0
    for i in range(prediction_step_nums):
        # projection into S-T graph
        st_graph.append([0, current_state_obstacle[0] - vehicle_length/2])
        last_state_obstacle = current_state_obstacle
        # prediction
        current_state_obstacle[0] = last_state_obstacle[0] + last_state_obstacle[3] * math.cos(last_state_obstacle[2]) * prediction_dt
        current_state_obstacle[1] = last_state_obstacle[1] + last_state_obstacle[3] * math.sin(last_state_obstacle[2]) * prediction_dt
        current_state_obstacle[2] = 0
        current_state_obstacle[3] = last_state_obstacle[3] - obstacle_acc * prediction_dt
    
    return st_graph
