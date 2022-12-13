#pragma once

#include "cube_type.h"
#include "piecewise_jerk_problem.h"

#include <iostream>
#include <fstream>
#include <utility>
#include <vector>
#include <eigen3/Eigen/Dense>

namespace navigation
{

	class PiecewiseJerkSpeedProblem : public PiecewiseJerkProblem
	{
	public:
		PiecewiseJerkSpeedProblem(const size_t num_of_knots, const double delta_s,
								  const std::array<double, 3> &x_init);

		virtual ~PiecewiseJerkSpeedProblem() = default;

		void set_dx_ref(const double weight_dx_ref, const double dx_ref);

		void set_penalty_dx(std::vector<double> penalty_dx);

		bool Optimize(const int max_iter = 4000) override;

		void CorridorGeneration();
		void PrintCorridor();
		void CorridorVisualize();

	protected:
		// naming convention follows osqp solver.
		void CalculateKernel(std::vector<c_float> *P_data,
							 std::vector<c_int> *P_indices,
							 std::vector<c_int> *P_indptr) override;

		void CalculateOffset(std::vector<c_float> *q) override;

		void CalculateAffineConstraint(std::vector<c_float> *A_data,
									   std::vector<c_int> *A_indices,
									   std::vector<c_int> *A_indptr,
									   std::vector<c_float> *lower_bounds,
									   std::vector<c_float> *upper_bounds) override;

		void S_CalculateKernel(std::vector<c_float> *P_data,
							   std::vector<c_int> *P_indices,
							   std::vector<c_int> *P_indptr) override;

		void S_CalculateOffset(std::vector<c_float> *q) override;

		void S_CalculateAffineConstraint(std::vector<c_float> *A_data,
										 std::vector<c_int> *A_indices,
										 std::vector<c_int> *A_indptr,
										 std::vector<c_float> *lower_bounds,
										 std::vector<c_float> *upper_bounds) override;

		void L_CalculateKernel(std::vector<c_float> *P_data,
							   std::vector<c_int> *P_indices,
							   std::vector<c_int> *P_indptr) override;

		void L_CalculateOffset(std::vector<c_float> *q) override;

		void L_CalculateAffineConstraint(std::vector<c_float> *A_data,
										 std::vector<c_int> *A_indices,
										 std::vector<c_int> *A_indptr,
										 std::vector<c_float> *lower_bounds,
										 std::vector<c_float> *upper_bounds) override;

		// void CorridorGeneration();
		void CorridorMerge();
		void CorridorSplit();
		// void CorridorVisualize();
		double CalculateCost(const double c[], double scale_t, int q_beg);

		OSQPSettings *SolverDefaultSettings() override;

		OSQPData *FormulateProblem() override;

		int n_of_knots;
		int segment_num;
		int traj_order;
		int n_poly;
		int n;

		std::vector<double> x_skew_;
		std::vector<double> x_bias_;

		bool has_dx_ref_ = false;
		double weight_dx_ref_ = 0.0;
		double dx_ref_ = 0.0;
		double weight_end_ref_ = 0.0;

		Eigen::MatrixXd M;
		Eigen::MatrixXd MQM[4];
		std::vector<c_float> q;

		std::vector<double> penalty_dx_;
		std::vector<Cube> corridor;
	};

}
