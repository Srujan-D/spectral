#include "../include/btrapz/logging.h"
// #include "../include/btrapz/solve_3d.h"
#include "../include/btrapz/solve_3d.h"


// #include "piecewise_jerk_path_problem.h"

#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>

using namespace std;
using namespace navigation;
int main(int argc, char const *argv[])
{
	ifstream ifs("/home/srujan_d/RISS/code/riss/src/btrapz/src/c1.txt");
	// ifstream c2 ("/home/srujan_d/RISS/code/riss/src/btrapz/src/c2.txt");

	double mscale = 1.0;
	int num_of_knots; // = 71;
	double delta_t, acc_weight, jerk_weight;

	double s_acc_weight, s_jerk_weight, l_acc_weight, l_jerk_weight;

	std::array<double,2> dds_bounds, ddds_bounds, ddl_bounds, dddl_bounds;
	double weight_s_ref, weight_ds_ref, weight_l_ref, weight_dl_ref, ds_ref, dl_ref;
	std::array<double, 3> init_s; // {0, 1.0, 0.0};
	int num_of_obs;

	ifs >> num_of_knots >> delta_t;
	ifs >> init_s[0] >> init_s[1] >> init_s[2];
	ifs >> s_acc_weight >> s_jerk_weight;
	ifs >> l_acc_weight >> l_jerk_weight;
	ifs >> num_of_obs;

	ifs >> weight_s_ref >> weight_ds_ref;
	ifs >> weight_l_ref >> weight_dl_ref;
	ifs >> ds_ref >> dl_ref;

	ifs >> dds_bounds[0] >> dds_bounds[1];
	ifs >> ddds_bounds[0] >> ddds_bounds[1];
	ifs >> ddl_bounds[0] >> ddl_bounds[1];
	ifs >> dddl_bounds[0] >> dddl_bounds[1];

	std::cout<<weight_s_ref<<" "<<weight_ds_ref<<" "<<weight_l_ref<<" "<<weight_dl_ref<<" "<<ds_ref<<" "<<dl_ref<<"\n";
	std::cout << "starting " << num_of_knots << " " << delta_t << " " << init_s[0] << " " << init_s[1] << " " << init_s[2] << " " << s_acc_weight << " " << s_jerk_weight << " " << num_of_obs << " \n";

	std::vector<std::pair<double, double>> x_bounds, dx_bounds, dy_bounds;
	std::vector<double> x_ref, y_ref, x_kappa, y_kappa;
	std::vector<std::vector<std::pair<double, double>>> all_x_bounds;
	init_s[1] = init_s[1] * mscale;

	double val[2];
	// PiecewiseJerkSpeedProblem *piecewise_jerk_problem_0 = (PiecewiseJerkSpeedProblem *) malloc(sizeof(PiecewiseJerkSpeedProblem)*num_of_obs);

	// for (int i=0; i<num_of_obs; i++)    piecewise_jerk_problem[i] = PiecewiseJerkSpeedProblem(num_of_knots, delta_t, init_s);
	std::cout << "lessgoo\n";
	PiecewiseJerkSpeedProblem piecewise_jerk_problem_0 = PiecewiseJerkSpeedProblem(num_of_knots, delta_t, init_s);
	PiecewiseJerkSpeedProblem piecewise_jerk_problem_1 = PiecewiseJerkSpeedProblem(num_of_knots, delta_t, init_s);
	// PiecewiseJerkSpeedProblem piecewise_jerk_problem(num_of_knots, delta_t, init_s);

	// for (int obs=0; obs<num_of_obs; obs++){
	for (int i = 0; i < num_of_knots; i++)
	{
		// std::cout<<"hey\n";
		ifs >> val[0] >> val[1];
		x_bounds.emplace_back(mscale * val[0], mscale * val[1]);
	}

	piecewise_jerk_problem_0.set_x_bounds(x_bounds);
	all_x_bounds.push_back(x_bounds);

	// }

	std::cout << "corridor generate 0\n";
	piecewise_jerk_problem_0.CorridorGeneration();
	piecewise_jerk_problem_0.PrintCorridor();
	x_bounds.clear();

	for (int i = 0; i < num_of_knots; i++)
	{
		// std::cout<<"hey\n";
		ifs >> val[0] >> val[1];
		x_bounds.emplace_back(mscale * val[0], mscale * val[1]);
	}
	piecewise_jerk_problem_0.set_x_bounds(x_bounds);
	all_x_bounds.push_back(x_bounds);
	x_bounds.clear();

	std::cout << "corridor generate 1\n";

	piecewise_jerk_problem_0.CorridorGeneration();
	piecewise_jerk_problem_0.PrintCorridor();

	for (int i = 0; i < num_of_knots; i++)
	{
		ifs >> val[0] >> val[1];
		dx_bounds.emplace_back(mscale * val[0], mscale * val[1]);
	}

	for (int i = 0; i < num_of_knots; i++)
	{
		ifs >> val[0] >> val[1];
		dy_bounds.emplace_back(mscale * val[0], mscale * val[1]);
	}

	for (int i = 0; i < num_of_knots; i++) {
		ifs >> val[0];
		x_ref.emplace_back(mscale * val[0]);

		// if (x_ref[i] < x_bounds[i].second && x_ref[i] >= x_bounds[i].first)
		// 	cout << "\nx_ref valid\n";
		// else
		// 	cout << "\nx_ref invalid!!--!!\n";
	}

	for (int i = 0; i < num_of_knots; i++)
	{
		ifs >> val[0];
		y_ref.emplace_back(mscale * val[0]);
	}

	for (int i = 0; i < num_of_knots; i++)
	{
		ifs >> val[0];
		x_kappa.emplace_back(val[0]);
	}

	for (int i = 0; i < num_of_knots; i++)
	{
		ifs >> val[0];
		y_kappa.emplace_back(val[0]);
	}
	// for (int i = 0; i < num_of_knots; i++) {
	// 	ifs >> val[0] >> val[1];
	// 	x_bounds.emplace_back(mscale * val[0], mscale * val[1]);
	// }

	// piecewise_jerk_problem[i].set_x_bounds(x_bounds);
	// all_x_bounds.push_back(x_bounds);

	piecewise_jerk_problem_0.FormNewCorridors();

	// new_corridor.push_back(piecewise_jerk_problem_0.corridors[0][0]);
	// new_corridor[0].beg_l=1.0;
	// new_corridor[0].end_l=3.0;

	// new_corridor.push_back(piecewise_jerk_problem_0.corridors[0][1]);
	// new_corridor[1].beg_l=1.0;
	// new_corridor[1].end_l=3.0;

	// new_corridor.push_back(piecewise_jerk_problem_0.corridors[0][2]);
	// new_corridor[2].beg_l=1.0;
	// new_corridor[2].end_l=3.0;

	// Cube mcube = piecewise_jerk_problem_0.corridors[0][3];
	// mcube.end_t = 36;

	// new_corridor.push_back(mcube);
	// new_corridor[3].beg_l=1.0;
	// new_corridor[3].end_l=3.0;

	// mcube = piecewise_jerk_problem_0.corridors[1][3];
	// mcube.beg_t= 34;

	// new_corridor.push_back(mcube);
	// new_corridor[4].beg_l=3.0;
	// new_corridor[4].end_l=4.5;

	// new_corridor.push_back(piecewise_jerk_problem_0.corridors[1][4]);
	// new_corridor[4].beg_l=3.0;
	// new_corridor[4].end_l=4.5;

	// new_corridor.push_back(piecewise_jerk_problem_0.corridors[1][5]);
	// new_corridor[4].beg_l=3.0;
	// new_corridor[4].end_l=4.5;

	// new_corridor.push_back(piecewise_jerk_problem_0.corridors[1][6]);
	// new_corridor[4].beg_l=3.0;
	// new_corridor[4].end_l=4.5;

	// new_corridor.push_back(piecewise_jerk_problem_0.corridors[1][7]);
	// new_corridor[4].beg_l=3.0;
	// new_corridor[4].end_l=4.5;

	// call optimize method with these corridors

	// for (int i=0; i<piecewise_jerk_problem_0.corridors.size(); i++){
	//     for (int j=0; j<num_of_knots; j++){

	//     }
	// }

	piecewise_jerk_problem_0.set_scale_factor({1.0, 1.0, 1.0});

	piecewise_jerk_problem_0.set_weight_ddx(s_acc_weight);
	piecewise_jerk_problem_0.set_weight_dddx(s_jerk_weight);
	piecewise_jerk_problem_0.set_ddx_bounds(dds_bounds[0], dds_bounds[1]);
	piecewise_jerk_problem_0.set_dddx_bound(ddds_bounds[0], ddds_bounds[1]);

	// cout<<"hey x_bounds.size() is : "<<x_bounds.size()<<"\n";
	// piecewise_jerk_problem_0.set_x_bounds(x_bounds);
	piecewise_jerk_problem_0.set_dx_bounds(dx_bounds);

	piecewise_jerk_problem_0.set_x_ref(weight_s_ref, x_ref);
	piecewise_jerk_problem_0.set_dx_ref(weight_ds_ref, ds_ref * mscale);

	piecewise_jerk_problem_0.set_weight_ddy(l_acc_weight);
	piecewise_jerk_problem_0.set_weight_dddy(l_jerk_weight);
	piecewise_jerk_problem_0.set_ddy_bounds(ddl_bounds[0], ddl_bounds[1]);
	piecewise_jerk_problem_0.set_dddy_bound(dddl_bounds[0], dddl_bounds[1]);

	// piecewise_jerk_problem_0.set_y_bounds(y_bounds);
	piecewise_jerk_problem_0.set_dy_bounds(dy_bounds);

	piecewise_jerk_problem_0.set_y_ref(weight_l_ref, y_ref);
	piecewise_jerk_problem_0.set_dy_ref(weight_dl_ref, dl_ref* mscale);	// try 0.0 * mscale

	auto start_time = std::chrono::system_clock::now();
	bool success = piecewise_jerk_problem_0.Optimize(5000);
	auto end_time = std::chrono::system_clock::now();
	std::chrono::duration<double> diff = end_time - start_time;
	ADEBUG << "Speed Optimizer used time: " << diff.count() * 1000 << " ms.";
	if (!success) {
		std::string msg("Piecewise jerk speed optimizer failed!");
		AERROR << msg;
	}
	// Extract output
	const std::vector<double> &s = piecewise_jerk_problem_0.opt_x();
	const std::vector<double> &ds = piecewise_jerk_problem_0.opt_dx();
	const std::vector<double> &dds = piecewise_jerk_problem_0.opt_ddx();
	bool is_verbose = true;
	double a_cost = 0.0;
	double mmax_a = 0.0;
	for (int i = 0; i < num_of_knots && is_verbose; ++i) {
		double ddds = (i == 0) ? (dds[1] - dds[0]) / delta_t : (dds[i] - dds[i - 1]) / delta_t;

		a_cost += dds[i] * dds[i] * delta_t;
		mmax_a = max(mmax_a, abs(dds[i]));

		cout << std::fixed << std::setprecision(3) << "For t[" << i * delta_t << "], opt = " << s[i] << ", " << ds[i] << ", " << dds[i] << ", " << ddds << endl;
	}

	a_cost = acc_weight * sqrt(a_cost / (num_of_knots * delta_t));
	cout << "mmax_a " << mmax_a << endl;
	cout << "a_cost " << a_cost << endl;

	cout<<"\n\n\nL part of trajectory: \n\n";
	const std::vector<double> &l = piecewise_jerk_problem_0.opt_y();
	const std::vector<double> &dl = piecewise_jerk_problem_0.opt_dy();
	const std::vector<double> &ddl = piecewise_jerk_problem_0.opt_ddy();
	is_verbose = true;
	a_cost = 0.0;
	mmax_a = 0.0;
	cout<<"\n\n\nprinting L part: \n\nl size is\t"<<l.size()<<"\n\n";
	for (int i = 0; i < num_of_knots && is_verbose; ++i) {
		double dddl = (i == 0) ? (ddl[1] - ddl[0]) / delta_t : (ddl[i] - ddl[i - 1]) / delta_t;

		a_cost += ddl[i] * ddl[i] * delta_t;
		mmax_a = max(mmax_a, abs(ddl[i]));

		cout << std::fixed << std::setprecision(3) << "For t[" << i * delta_t << "], opt = " << l[i] << ", " << dl[i] << ", " << ddl[i] << ", " << dddl << endl;
	}

	ofstream ofs("/home/srujan_d/RISS/code/riss/src/btrapz/src/cuboid_slt.txt");
	if (!ofs.is_open()) {
		cout << "Not Found" << endl;
		return 1;
	}

	for (int i = 0; i < num_of_knots && is_verbose; ++i) {
		ofs << std::fixed << std::setprecision(3) << i * delta_t << " " << s[i]	<< " " << l[i] << " "<< ds[i] << " " << dl[i] << " " << dds[i] << " " << ddl[i] << endl;
	}
	// ADEBUG << "Speed Optimizer used time: " << diff.count() * 1000 << " ms.";
	return 0;
}
