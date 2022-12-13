import numpy as np
import heapq
import math
# import flowstar
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.animation import FuncAnimation
from matplotlib.axis import Axis

class FittingPoints:
    def __init__(self):
        self.fitting_points = []
        self.fitting_points_cost = 0.0

class DerivateSpace:
    def __init__(self):
        self.a_min = -4.5
        self.a_max = 1.5
        self.a_resolution = 0.25
        self.discretized_a_num = math.floor((self.a_max - self.a_min) / self.a_resolution) + 1
        self.a_candidate_set = np.zeros(self.discretized_a_num, dtype = float)
        self.a_candidate_pairs = [[]]

        for i in range(self.discretized_a_num):
            for j in range(self.discretized_a_num):
                self.a_candidate_pairs[i].append(FittingPoints())
            self.a_candidate_pairs.append([])

    def CalculateFittingPoints(self, start_point, end_point):
        w_2 = 0.2
        w_3 = 2.0
        lambda_opt = w_3 / (w_2 + 2 * w_3)
        # the number of fitting points
        nums_of_fitting_points = 11
        n_f = nums_of_fitting_points - 2
        A_determinant = np.zeros(n_f+1, dtype = float)
        A_determinant[1] = 1
        A_determinant[2] = 1 - lambda_opt * lambda_opt
        for i in range(3, n_f+1):
            A_determinant[i] = A_determinant[i-1] - lambda_opt * lambda_opt * A_determinant[i-2]

        a_fitting_points = FittingPoints()
        a_fitting_points.fitting_points = list(range(nums_of_fitting_points))
        a_fitting_points.fitting_points[0] = start_point
        # dds[1] = lambda_opt
        if (n_f-1) % 2 == 0:
            a_fitting_points.fitting_points[1] = (lambda_opt/A_determinant[n_f]) * (A_determinant[n_f-1] * start_point +
                    math.pow(lambda_opt, n_f-1) * end_point)
        else:
            a_fitting_points.fitting_points[1] = (lambda_opt/A_determinant[n_f]) * (A_determinant[n_f-1] * start_point -
                    math.pow(lambda_opt, n_f-1) * end_point)
        a_fitting_points.fitting_points[nums_of_fitting_points-1] = end_point

        for i in range(1, n_f):
            a_fitting_points.fitting_points[i+1] = 1/lambda_opt * a_fitting_points.fitting_points[i] - a_fitting_points.fitting_points[i-1]
        
        return a_fitting_points

    def GenerateDdsCandidateSet(self):
        for i in range(self.discretized_a_num):
            self.a_candidate_set[i] = self.a_min + i * self.a_resolution

        for i in range(self.discretized_a_num):
            for j in range(self.discretized_a_num):
                a_fitting_points = self.CalculateFittingPoints(self.a_candidate_set[i], self.a_candidate_set[j])
                num_fitting_points = len(a_fitting_points.fitting_points)
                self.a_candidate_pairs[i][j].fitting_points = list(range(num_fitting_points))
                for k in range(num_fitting_points):
                    self.a_candidate_pairs[i][j].fitting_points[k] = a_fitting_points.fitting_points[k]

derivate_space = DerivateSpace()
derivate_space.GenerateDdsCandidateSet()
print(derivate_space.a_candidate_pairs[18][10].fitting_points)
nums_of_fitting_points = len(derivate_space.a_candidate_pairs[18][10].fitting_points)
plot_a = []
plot_t = []
delta_t = 0.1
for i in range(nums_of_fitting_points):
    plot_t.append(i * delta_t)
    plot_a.append(derivate_space.a_candidate_pairs[18][10].fitting_points[i])
plt.plot(plot_t, plot_a)
plt.show()
