#pragma once

#include "cube_type.h"
#include "piecewise_jerk_problem.h"

#include <iostream>
#include <fstream>
#include <utility>
#include <vector>
#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Sparse>

namespace navigation
{

	class PiecewiseJerkSpeedProblem : public PiecewiseJerkProblem
	{
	public:
		PiecewiseJerkSpeedProblem(const size_t num_of_knots, const double delta_s,
								  const std::array<double, 3> &x_init, const std::array<double, 3> &y_init, const int num_of_obs);

		virtual ~PiecewiseJerkSpeedProblem() = default;

		void set_dx_ref(const double weight_dx_ref, const double dx_ref);

		void set_penalty_dx(std::vector<double> penalty_dx);

		void set_dy_ref(const double weight_dy_ref, const double dy_ref);

		void set_penalty_dy(std::vector<double> penalty_dy);

		void set_weight_end(const double weight_end_s, const double weight_end_l);

		bool Optimize(const int max_iter = 4000) override;
		// bool S_Optimize(const int max_iter = 4000);
		// bool L_Optimize(const int max_iter = 4000);

		void CorridorGeneration();
		void PrintCorridor();
		void CorridorVisualize();

		struct waypoint
		{
			double s;
			double l;
			double t;
		};

		void CollisionCheck();

		std::vector<waypoint> traj;

		std::vector<std::vector<Cube>> corridors;
		void FormNewCorridors();
		std::vector<Cube> new_corridor;

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

		// void S_CalculateKernel(std::vector<c_float>* P_data,
		// 				   std::vector<c_int>* P_indices,
		// 				   std::vector<c_int>* P_indptr) override;

		// void S_CalculateOffset(std::vector<c_float>* q) override;

		// void S_CalculateAffineConstraint(std::vector<c_float>* A_data,
		// 							std::vector<c_int>* A_indices,
		// 							std::vector<c_int>* A_indptr,
		// 							std::vector<c_float>* lower_bounds,
		// 							std::vector<c_float>* upper_bounds) override;

		// void L_CalculateKernel(std::vector<c_float>* P_data,
		// 				   std::vector<c_int>* P_indices,
		// 				   std::vector<c_int>* P_indptr) override;

		// void L_CalculateOffset(std::vector<c_float>* q) override;

		// void L_CalculateAffineConstraint(std::vector<c_float>* A_data,
		// 							std::vector<c_int>* A_indices,
		// 							std::vector<c_int>* A_indptr,
		// 							std::vector<c_float>* lower_bounds,
		// 							std::vector<c_float>* upper_bounds) override;

		// void CorridorGeneration();
		void CorridorMerge();
		void CorridorSplit();
		// void CorridorVisualize();
		double CalculateCost(const double c[], double scale_t, int q_beg);
		double CalculateCost(const double c[], double scale_t, int q_beg, std::vector<c_float> q);

		OSQPSettings *SolverDefaultSettings() override;

		OSQPData *FormulateProblem() override;

		// OSQPData* S_FormulateProblem();
		// OSQPData* L_FormulateProblem();

		int n_of_knots;
		int segment_num;
		int traj_order;
		int n_poly;
		int n;

		int num_of_points_ = 1;

		std::vector<double> x_skew_;
		std::vector<double> x_bias_;

		std::vector<double> y_skew_;
		std::vector<double> y_bias_;

		// std::vector<double> y_;

		bool has_dx_ref_ = false;
		double weight_dx_ref_ = 0.0;
		double dx_ref_ = 0.0;

		bool has_dy_ref_ = false;
		double weight_dy_ref_ = 0.0;
		double dy_ref_ = 0.0;

		double weight_end_ref_ = 0.0;
		double weight_end_s_ = 0.0;
		double weight_end_l_ = 0.0;

		typedef Eigen::Triplet<int> Trip;
		std::vector<Trip> trp;

		Eigen::MatrixXd M;

		Eigen::MatrixXd MQM_x[4];
		Eigen::MatrixXd MQM_y[4];

		std::vector<c_float> q;
		std::vector<c_float> S_q;
		std::vector<c_float> L_q;

		std::vector<double> penalty_dx_;
		std::vector<double> penalty_dy_;

		std::vector<Cube> corridor;

		// std::vector<std::vector <Cube>> corridors;
	};

}
