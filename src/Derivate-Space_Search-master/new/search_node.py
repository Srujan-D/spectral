import math
from re import S
from this import s
from vehicle_state import VehicleState
from derivate_space_class import DerivateSpace

class TranverseState:
    def __init__(self, tranverse_s, tranverse_ds, tranverse_dds, tranverse_l, tranverse_dl, tranverse_ddl, tranverse_cost):
        self.tranverse_s = tranverse_s
        self.tranverse_ds = tranverse_ds
        self.tranverse_dds = tranverse_dds

        self.tranverse_l = tranverse_l
        self.tranverse_dl = tranverse_dl
        self.tranverse_ddl = tranverse_ddl

        self.tranverse_cost = tranverse_cost

class Node:
    def __init__(self, vehicle_state: VehicleState, tranverse_state: TranverseState) -> None:
        self.vehicle_state = vehicle_state
        self.tranverse_state = tranverse_state
        self.s_resolution = 6.0
        self.ds_resolution = 4.0
        self.dds_resolution = 0.25

        self.s_ind = math.ceil(self.vehicle_state.s / self.s_resolution)
        self.ds_ind = math.ceil(self.vehicle_state.ds / self.ds_resolution)
        self.dds_ind = math.ceil(self.vehicle_state.dds / self.dds_resolution)

        self.l_resolution = 6.0
        self.dl_resolution = 4.0
        self.ddl_resolution = 0.25

        self.l_ind = math.ceil(self.vehicle_state.l / self.l_resolution)
        self.dl_ind = math.ceil(self.vehicle_state.dl / self.dl_resolution)
        self.ddl_ind = math.ceil(self.vehicle_state.ddl / self.ddl_resolution)

        self.id = 0
        self.total_cost = 0.0
        self.father_cost = 0.0
        self.comfort_cost = 0.0
        self.obstacle_cost = 0.0
        self.velocity_cost = 0.0
        self.heuristic_cost = 0.0

    def SetFatherNode(self, father_node)-> None:
        self.father = father_node.id
        self.father_cost = father_node.GetNodeCost()

    # def CalculateNodeCost(self, goal_state: LongitudinalState, ds_ref):
    #     self.total_cost = self.father_cost + self.comfort_cost + 

    def GetNodeCost(self):
        return self.total_cost

    def GetTime(self):
        return self.t
        
    def GetNodeIndex(self):
        return str(self.s_ind) + '_' + str(self.ds_ind) + '_' + str(self.dds_ind)

    def GenerateNextNode(self, tranverse_state):
        self.derivate_space = DerivateSpace()


