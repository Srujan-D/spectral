import numpy as np
import math


class FittingPoints:
    def __init__(self):
        self.fitting_points = []


class DerivateSpace:
    def __init__(self):
        self.a_min = -4.5
        self.a_max = 1.5
        self.a_resolution = 0.25
        self.discretized_a_num = math.floor(
            (self.a_max - self.a_min) / self.a_resolution) + 1
        self.a_candidate_set = np.zeros(self.discretized_a_num, dtype=float)
        self.a_tranverse_pairs = []
        self.s_tranverse_pairs = []
        self.ds_tranverse_pairs = []
        self.transition_cost_pairs = []

        self.a_candidate_pairs = []

    def TraversePairsInit(self):
        for i in range(self.discretized_a_num):
            self.a_candidate_pairs.append([])
            self.s_tranverse_pairs.append([])
            self.ds_tranverse_pairs.append([])
            for j in range(self.discretized_a_num):
                self.a_candidate_pairs[i].append([])
                self.s_tranverse_pairs.append([])
                self.ds_tranverse_pairs.append([])

    def CalculateFittingPoints(self, start_point, end_point):
        w_2 = 0.2
        w_3 = 2.0
        lambda_opt = w_3 / (w_2 + 2 * w_3)
        # the number of fitting points
        total_fitting_points = 11
        nums_of_insert_points = total_fitting_points - 2
        A_determinant = np.zeros(nums_of_insert_points+1, dtype=float)
        A_determinant[1] = 1
        A_determinant[2] = 1 - lambda_opt * lambda_opt
        for i in range(3, nums_of_insert_points+1):
            A_determinant[i] = A_determinant[i-1] - \
                lambda_opt * lambda_opt * A_determinant[i-2]

        a_tranverse = FittingPoints()
        a_tranverse.fitting_points = list(range(total_fitting_points))
        a_tranverse.fitting_points[0] = start_point
        # dds[1] = lambda_opt
        if (nums_of_insert_points-1) % 2 == 0:
            a_tranverse.fitting_points[1] = (lambda_opt/A_determinant[nums_of_insert_points]) * (
                A_determinant[nums_of_insert_points-1] * start_point 
                + math.pow(lambda_opt, nums_of_insert_points-1) * end_point
                )
        else:
            a_tranverse.fitting_points[1] = (lambda_opt/A_determinant[nums_of_insert_points]) * (
                A_determinant[nums_of_insert_points-1] * start_point
                - math.pow(lambda_opt, nums_of_insert_points-1) * end_point
                )

        a_tranverse.fitting_points[total_fitting_points-1] = end_point

        for i in range(1, nums_of_insert_points):
            a_tranverse.fitting_points[i+1] = 1/lambda_opt * \
                a_tranverse.fitting_points[i] - \
                a_tranverse.fitting_points[i-1]

        return a_tranverse

    def GenerateDdsCandidateSet(self):
        
        for i in range(self.discretized_a_num):
            self.a_candidate_set[i] = self.a_min + i * self.a_resolution

        for i in range(self.discretized_a_num):
            for j in range(self.discretized_a_num):
                a_fitting_points = self.CalculateFittingPoints(
                    self.a_candidate_set[i], self.a_candidate_set[j])
                num_fitting_points = len(a_fitting_points.fitting_points)
                self.a_candidate_pairs[i][j].fitting_points =  list(range(num_fitting_points)) #[i for i in range(num_fitting_points)] 
                # list(range(num_fitting_points))
                for k in range(num_fitting_points):
                    self.a_candidate_pairs[i][j].fitting_points[k] = a_fitting_points.fitting_points[k]

    # def CalculateDdsInd(self, current_dds):
    #     return (current_dds - self.a_min) / self.a_resolution

    # def GetTranverseInd(self, current_dds, next_dds):
    #     current_ind = self.CalculateDdsInd(current_dds)
    #     next_ind = self.CalculateDdsInd(next_dds)
    #     return current_ind, next_ind

    # def GenerateNextNode():
