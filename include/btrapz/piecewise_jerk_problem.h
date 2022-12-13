#pragma once

#include <tuple>
#include <utility>
#include <vector>

// #include <osqp/include/osqp.h>
#include <osqp/osqp.h>

namespace navigation {

class PiecewiseJerkProblem {
 public:
	PiecewiseJerkProblem(const size_t num_of_knots, const double delta_s,
											 const std::array<double, 3>& x_init, const std::array<double, 3>& y_init, const int num_of_obs);

	virtual ~PiecewiseJerkProblem() = default;

	std::vector<std::pair<double, double>> get_x_bounds() const {
		return x_bounds_;
	}

	std::vector<std::pair<double, double>> get_dx_bounds() const {
		return dx_bounds_;
	}

	void set_x_bounds(std::vector<std::pair<double, double>> x_bounds);

	void set_x_bounds(const double x_lower_bound, const double x_upper_bound);

	void set_dx_bounds(std::vector<std::pair<double, double>> dx_bounds);

	void set_dx_bounds(const double dx_lower_bound, const double dx_upper_bound);

	void set_ddx_bounds(std::vector<std::pair<double, double>> ddx_bounds);

	void set_ddx_bounds(const double ddx_lower_bound,
						const double ddx_upper_bound);

	void set_dddx_bound(const double dddx_bound) {
		set_dddx_bound(-dddx_bound, dddx_bound);
	}

	void set_dddx_bound(const double dddx_lower_bound,
						const double dddx_upper_bound) {
		dddx_bound_.first = dddx_lower_bound;
		dddx_bound_.second = dddx_upper_bound;
	}

	void set_weight_x(const double weight_x) { weight_x_ = weight_x; }

	void set_weight_dx(const double weight_dx) { weight_dx_ = weight_dx; }

	void set_weight_ddx(const double weight_ddx) { weight_ddx_ = weight_ddx; }

	void set_weight_dddx(const double weight_dddx) { weight_dddx_ = weight_dddx; }

	void set_scale_factor(const std::array<double, 3>& scale_factor) {
		scale_factor_ = scale_factor;
	}

	/**
	 * @brief Set the x ref object and the uniform x_ref weighting
	 *
	 * @param weight_x_ref: uniform weighting for x_ref
	 * @param x_ref: objective value of x
	 */
	void set_x_ref(const double weight_x_ref, std::vector<double> x_ref);

	/**
	 * @brief Set the x ref object and piecewised x_ref weightings
	 *
	 * @param weight_x_ref_vec: piecewised x_ref weightings
	 * @param x_ref: objective value of x
	 */
	void set_x_ref(std::vector<double> weight_x_ref_vec,
					std::vector<double> x_ref);

	void set_end_state_ref(const std::array<double, 3>& weight_end_state,
							const std::array<double, 3>& end_state_ref);

	virtual bool Optimize(const int max_iter) = 0;

	const std::vector<double>& opt_x() const { return x_; }

	const std::vector<double>& opt_dx() const { return dx_; }

	const std::vector<double>& opt_ddx() const { return ddx_; }

// 

	std::vector<std::pair<double, double>> get_y_bounds() const {
		return y_bounds_;
	}

	std::vector<std::pair<double, double>> get_dy_bounds() const {
		return dy_bounds_;
	}

	void set_y_bounds(std::vector<std::pair<double, double>> y_bounds);

	void set_y_bounds(const double y_lower_bound, const double y_upper_bound);

	void set_dy_bounds(std::vector<std::pair<double, double>> dy_bounds);

	void set_dy_bounds(const double dy_lower_bound, const double dy_upper_bound);

	void set_ddy_bounds(std::vector<std::pair<double, double>> ddy_bounds);

	void set_ddy_bounds(const double ddy_lower_bound,
						const double ddy_upper_bound);

	void set_dddy_bound(const double dddy_bound) {
		set_dddy_bound(-dddy_bound, dddy_bound);
	}

	void set_dddy_bound(const double dddy_lower_bound,
						const double dddy_upper_bound) {
		dddy_bound_.first = dddy_lower_bound;
		dddy_bound_.second = dddy_upper_bound;
	}

	void set_weight_y(const double weight_y) { weight_y_ = weight_y; }

	void set_weight_dy(const double weight_dy) { weight_dy_ = weight_dy; }

	void set_weight_ddy(const double weight_ddy) { weight_ddy_ = weight_ddy; }

	void set_weight_dddy(const double weight_dddy) { weight_dddy_ = weight_dddy; }

	// void set_scale_factor(const std::array<double, 3>& scale_factor) {
	// 	scale_factor_ = scale_factor;
	// }

	/**
	 * @brief Set the y ref object and the uniform y_ref weighting
	 *
	 * @param weight_y_ref: uniform weighting for y_ref
	 * @param y_ref: objective value of y
	 */
	void set_y_ref(const double weight_y_ref, const std::vector<double> y_ref);

	/**
	 * @brief Set the y ref object and piecewised y_ref weightings
	 *
	 * @param weight_y_ref_vec: piecewised y_ref weightings
	 * @param y_ref: objective value of y
	 */
	void set_y_ref(std::vector<double> weight_y_ref_vec,
					std::vector<double> y_ref);

	// void set_end_state_ref(const std::array<double, 3>& weight_end_state,
	// 						const std::array<double, 3>& end_state_ref);

	const std::vector<double>& opt_y() const { return y_; }

	const std::vector<double>& opt_dy() const { return dy_; }

	const std::vector<double>& opt_ddy() const { return ddy_; }

 protected:
	// naming convention follows osqp solver.
	virtual void CalculateKernel(std::vector<c_float>* P_data,
								std::vector<c_int>* P_indices,
								std::vector<c_int>* P_indptr) = 0;

	virtual void CalculateOffset(std::vector<c_float>* q) = 0;

	virtual void CalculateAffineConstraint(std::vector<c_float>* A_data,
											std::vector<c_int>* A_indices,
											std::vector<c_int>* A_indptr,
											std::vector<c_float>* lower_bounds,
											std::vector<c_float>* upper_bounds) = 0;

	// virtual void S_CalculateKernel(std::vector<c_float>* P_data,
	// 							std::vector<c_int>* P_indices,
	// 							std::vector<c_int>* P_indptr) = 0;

	// virtual void S_CalculateOffset(std::vector<c_float>* q) = 0;

	// virtual void S_CalculateAffineConstraint(std::vector<c_float>* A_data,
	// 										std::vector<c_int>* A_indices,
	// 										std::vector<c_int>* A_indptr,
	// 										std::vector<c_float>* lower_bounds,
	// 										std::vector<c_float>* upper_bounds) = 0;

	// virtual void L_CalculateKernel(std::vector<c_float>* P_data,
	// 							std::vector<c_int>* P_indices,
	// 							std::vector<c_int>* P_indptr) = 0;

	// virtual void L_CalculateOffset(std::vector<c_float>* q) = 0;

	// virtual void L_CalculateAffineConstraint(std::vector<c_float>* A_data,
	// 										std::vector<c_int>* A_indices,
	// 										std::vector<c_int>* A_indptr,
	// 										std::vector<c_float>* lower_bounds,
	// 										std::vector<c_float>* upper_bounds) = 0;											

	virtual OSQPSettings* SolverDefaultSettings() = 0;

	virtual OSQPData* FormulateProblem() = 0;

	void FreeData(OSQPData* data);

	template <typename T>
	T* CopyData(const std::vector<T>& vec) {
		T* data = new T[vec.size()];
		memcpy(data, vec.data(), sizeof(T) * vec.size());
		return data;
	}

 protected:
	size_t num_of_knots_ = 0;

	// output
	std::vector<double> x_;
	std::vector<double> dx_;
	std::vector<double> ddx_;
	
	std::vector<double> y_;
	std::vector<double> dy_;
	std::vector<double> ddy_;

	std::array<double, 3> x_init_;
	std::array<double, 3> y_init_ = {{2.0, 0.0, 0.0}};
	std::array<double, 3> scale_factor_ = {{1.0, 1.0, 1.0}};

	std::vector<std::pair<double, double>> x_bounds_;
	std::vector<std::pair<double, double>> dx_bounds_;
	std::vector<std::pair<double, double>> ddx_bounds_;
	std::pair<double, double> dddx_bound_;

	std::vector<std::pair<double, double>> y_bounds_;
	std::vector<std::pair<double, double>> dy_bounds_;
	std::vector<std::pair<double, double>> ddy_bounds_;
	std::pair<double, double> dddy_bound_;

	double weight_x_ = 0.0;
	double weight_dx_ = 0.0;
	double weight_ddx_ = 0.0;
	double weight_dddx_ = 0.0;

	double weight_y_ = 0.0;
	double weight_dy_ = 0.0;
	double weight_ddy_ = 0.0;
	double weight_dddy_ = 0.0;

	int num_of_obs_;

	double delta_s_ = 1.0;

	bool has_x_ref_ = false;
	double weight_x_ref_ = 0.0;
	std::vector<double> x_ref_;
	// un-uniformed weighting
	std::vector<double> weight_x_ref_vec_;

	bool has_y_ref_ = false;
	double weight_y_ref_ = 0.0;
	std::vector<double> y_ref_;
	// un-uniformed weighting
	std::vector<double> weight_y_ref_vec_;

	bool has_end_state_ref_ = false;
	std::array<double, 3> weight_end_state_ = {{0.0, 0.0, 0.0}};
	std::array<double, 3> end_state_ref_;
};

}
