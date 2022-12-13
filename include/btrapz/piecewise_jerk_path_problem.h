#pragma once

#include <utility>
#include <vector>

#include "piecewise_jerk_problem.h"

namespace navigation {

class PiecewiseJerkPathProblem : public PiecewiseJerkProblem {
 public:
  PiecewiseJerkPathProblem(const size_t num_of_knots, const double delta_s, const std::array<double, 3>& x_init);

  bool Optimize(const int max_iter) override;
  virtual ~PiecewiseJerkPathProblem() = default;

 protected:
  OSQPData* FormulateProblem() override;
  void CalculateKernel(std::vector<c_float>* P_data,
					   std::vector<c_int>* P_indices,
					   std::vector<c_int>* P_indptr) override;

  void CalculateOffset(std::vector<c_float>* q) override;

  void CalculateAffineConstraint(std::vector<c_float>* A_data,
								std::vector<c_int>* A_indices,
								std::vector<c_int>* A_indptr,
								std::vector<c_float>* lower_bounds,
								std::vector<c_float>* upper_bounds) override;
  OSQPSettings* SolverDefaultSettings() override;
};
}
