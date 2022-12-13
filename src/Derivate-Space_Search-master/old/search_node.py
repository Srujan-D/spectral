from __future__ import annotations # https://stackoverflow.com/questions/42845972/typed-python-using-the-classes-own-type-inside-class-definition
import math
from re import S
from this import s
from vehicle_state import LongitudinalState
from derivate_space_class import DerivateSpace

class TranverseState:
    def __init__(self, tranverse_s, tranverse_ds, tranverse_dds, tranverse_cost):
        self.tranverse_s = tranverse_s
        self.tranverse_ds = tranverse_ds
        self.tranverse_dds = tranverse_dds
        self.tranverse_cost = tranverse_cost

class Node:
    def __init__(self, longitudinal_state: LongitudinalState, tranverse_state: TranverseState) -> None:
        self.longitudinal_state = longitudinal_state
        self.tranverse_state = tranverse_state
        self.s_resolution = 6.0
        self.ds_resolution = 4.0
        self.dds_resolution = 0.25

        self.s_ind = math.ceil(self.longitudinal_state.s / self.s_resolution)
        self.ds_ind = math.ceil(self.longitudinal_state.ds / self.ds_resolution)
        self.dds_ind = math.ceil(self.longitudinal_state.dds / self.dds_resolution)
        self.id = 0
        self.total_cost = 0.0
        self.father_cost = 0.0
        self.comfort_cost = 0.0
        self.obstacle_cost = 0.0
        self.velocity_cost = 0.0
        self.heuristic_cost = 0.0

    def SetFatherNode(self, father_node: Node)-> None:
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


