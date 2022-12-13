#include "../include/btrapz/piecewise_jerk_problem.h"
#include "../include/btrapz/logging.h"

#include <iostream>

namespace navigation {

namespace {
constexpr double kMaxVariableRange = 1.0e10;
}  // namespace

PiecewiseJerkProblem::PiecewiseJerkProblem(const size_t num_of_knots, const double delta_s,
											 const std::array<double, 3>& x_init, const std::array<double, 3>& y_init, const int num_of_obs) {
  CHECK_GE(num_of_knots, 2);
  num_of_knots_ = num_of_knots;

  x_init_ = x_init;

  y_init_ = y_init;

  num_of_obs_ = num_of_obs;

  delta_s_ = delta_s;

  x_bounds_.resize(num_of_knots_,std::make_pair(-kMaxVariableRange, kMaxVariableRange));
 
  dx_bounds_.resize(num_of_knots_, std::make_pair(-kMaxVariableRange, kMaxVariableRange));

  ddx_bounds_.resize(num_of_knots_, std::make_pair(-kMaxVariableRange, kMaxVariableRange));

  y_bounds_.resize(num_of_knots_,std::make_pair(-kMaxVariableRange, kMaxVariableRange));
 
  dy_bounds_.resize(num_of_knots_, std::make_pair(-kMaxVariableRange, kMaxVariableRange));

  ddy_bounds_.resize(num_of_knots_, std::make_pair(-kMaxVariableRange, kMaxVariableRange));

  weight_x_ref_vec_ = std::vector<double>(num_of_knots_, 0.0);
}


void PiecewiseJerkProblem::set_x_bounds(
    std::vector<std::pair<double, double>> x_bounds) {
  CHECK_EQ(x_bounds.size(), num_of_knots_);
  x_bounds_ = std::move(x_bounds);
}

void PiecewiseJerkProblem::set_dx_bounds(
    std::vector<std::pair<double, double>> dx_bounds) {
  CHECK_EQ(dx_bounds.size(), num_of_knots_);
  dx_bounds_ = std::move(dx_bounds);
}

void PiecewiseJerkProblem::set_ddx_bounds(
    std::vector<std::pair<double, double>> ddx_bounds) {
  CHECK_EQ(ddx_bounds.size(), num_of_knots_);
  ddx_bounds_ = std::move(ddx_bounds);
}

void PiecewiseJerkProblem::set_x_bounds(const double x_lower_bound,
                                        const double x_upper_bound) {
  for (auto& x : x_bounds_) {
    x.first = x_lower_bound;
    x.second = x_upper_bound;
  }
}

void PiecewiseJerkProblem::set_dx_bounds(const double dx_lower_bound,
                                         const double dx_upper_bound) {
  for (auto& x : dx_bounds_) {
    x.first = dx_lower_bound;
    x.second = dx_upper_bound;
  }
}

void PiecewiseJerkProblem::set_ddx_bounds(const double ddx_lower_bound,
                                          const double ddx_upper_bound) {
  for (auto& x : ddx_bounds_) {
    x.first = ddx_lower_bound;
    x.second = ddx_upper_bound;
  }
}

void PiecewiseJerkProblem::set_x_ref(const double weight_x_ref,
                                     std::vector<double> x_ref) {
  CHECK_EQ(x_ref.size(), num_of_knots_);
  weight_x_ref_ = weight_x_ref;
  // set uniform weighting
  weight_x_ref_vec_ = std::vector<double>(num_of_knots_, weight_x_ref);
  x_ref_ = std::move(x_ref);
  has_x_ref_ = true;
}

void PiecewiseJerkProblem::set_x_ref(std::vector<double> weight_x_ref_vec,
                                     std::vector<double> x_ref) {
  CHECK_EQ(x_ref.size(), num_of_knots_);
  CHECK_EQ(weight_x_ref_vec.size(), num_of_knots_);
  // set piecewise weighting
  weight_x_ref_vec_ = std::move(weight_x_ref_vec);
  x_ref_ = std::move(x_ref);
  has_x_ref_ = true;
}

void PiecewiseJerkProblem::set_y_bounds(
    std::vector<std::pair<double, double>> y_bounds) {
  CHECK_EQ(y_bounds.size(), num_of_knots_);
  y_bounds_ = std::move(y_bounds);
}

void PiecewiseJerkProblem::set_dy_bounds(
    std::vector<std::pair<double, double>> dy_bounds) {
  CHECK_EQ(dy_bounds.size(), num_of_knots_);
  dy_bounds_ = std::move(dy_bounds);
}

void PiecewiseJerkProblem::set_ddy_bounds(
    std::vector<std::pair<double, double>> ddy_bounds) {
  CHECK_EQ(ddy_bounds.size(), num_of_knots_);
  ddy_bounds_ = std::move(ddy_bounds);
}

void PiecewiseJerkProblem::set_y_bounds(const double y_lower_bound,
                                        const double y_upper_bound) {
  for (auto& y : y_bounds_) {
    y.first = y_lower_bound;
    y.second = y_upper_bound;
  }
}

void PiecewiseJerkProblem::set_dy_bounds(const double dy_lower_bound,
                                         const double dy_upper_bound) {
  for (auto& y : dy_bounds_) {
    y.first = dy_lower_bound;
    y.second = dy_upper_bound;
  }
}

void PiecewiseJerkProblem::set_ddy_bounds(const double ddy_lower_bound,
                                          const double ddy_upper_bound) {
  for (auto& y : ddy_bounds_) {
    y.first = ddy_lower_bound;
    y.second = ddy_upper_bound;
  }
}

void PiecewiseJerkProblem::set_y_ref(const double weight_y_ref,
                                     std::vector<double> y_ref) {
  CHECK_EQ(y_ref.size(), num_of_knots_);
  weight_y_ref_ = weight_y_ref;
  // set uniform weighting
  weight_y_ref_vec_ = std::vector<double>(num_of_knots_, weight_y_ref);
  y_ref_ = std::move(y_ref);
  has_y_ref_ = true;
}

void PiecewiseJerkProblem::set_y_ref(std::vector<double> weight_y_ref_vec,
                                     std::vector<double> y_ref) {
  CHECK_EQ(y_ref.size(), num_of_knots_);
  CHECK_EQ(weight_y_ref_vec.size(), num_of_knots_);
  // set piecewise weighting
  weight_y_ref_vec_ = std::move(weight_y_ref_vec);
  y_ref_ = std::move(y_ref);
  has_y_ref_ = true;
}

void PiecewiseJerkProblem::set_end_state_ref(
    const std::array<double, 3>& weight_end_state,
    const std::array<double, 3>& end_state_ref) {
  weight_end_state_ = weight_end_state;
  end_state_ref_ = end_state_ref;
  has_end_state_ref_ = true;
}

void PiecewiseJerkProblem::FreeData(OSQPData* data) {
  delete[] data->q;
  delete[] data->l;
  delete[] data->u;

  delete[] data->P->i;
  delete[] data->P->p;
  delete[] data->P->x;

  delete[] data->A->i;
  delete[] data->A->p;
  delete[] data->A->x;
}

}
