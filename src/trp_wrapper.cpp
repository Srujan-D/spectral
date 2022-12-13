#include "../include/btrapz/logging.h"
#include "../include/btrapz/solve_3d.h"
#include "../include/btrapz/py_cpp_.h"

// #include "piecewise_jerk_path_problem.h"

#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>

using namespace std;
using namespace navigation;

extern "C"
{

	// Metrics * find_traj(Metrics *m)
	double find_traj(Params *p)
	{
		double acc;
		ifstream ifs("/home/srujan_d/RISS/code/btrapz/src/c_road_s1_2.txt");
		// ifstream ifs(path);

		double mscale = 1.0;
		int num_of_knots; // = 71;
		double delta_t;

		double s_acc_weight, s_jerk_weight, l_acc_weight, l_jerk_weight;

		std::array<double, 2> dds_bounds, ddds_bounds, ddl_bounds, dddl_bounds;
		double weight_s_ref, weight_ds_ref, weight_l_ref, weight_dl_ref, ds_ref, dl_ref, weight_end_s, weight_end_l;
		std::array<double, 3> init_s, init_l; // {0, 1.0, 0.0};
		int num_of_obs;

		int iteration = p->iteration;

		ifs >> num_of_knots >> delta_t;
		ifs >> init_s[0] >> init_s[1] >> init_s[2];
		ifs >> init_l[0] >> init_l[1] >> init_l[2];

		s_acc_weight = p->s_acc_weight;
		s_jerk_weight = p->s_jerk_weight;

		l_acc_weight = p->l_acc_weight;
		l_jerk_weight = p->l_jerk_weight;

		ifs >> num_of_obs;

		weight_dl_ref = p->weight_dl_ref;
		weight_ds_ref = p->weight_ds_ref;
		weight_l_ref = p->weight_l_ref;
		weight_s_ref = p->weight_s_ref;

		weight_end_s = p->weight_end_s;
		weight_end_l = p->weight_end_l;

		ifs >> ds_ref >> dl_ref;

		ifs >> dds_bounds[0] >> dds_bounds[1];
		ifs >> ddds_bounds[0] >> ddds_bounds[1];
		ifs >> ddl_bounds[0] >> ddl_bounds[1];
		ifs >> dddl_bounds[0] >> dddl_bounds[1];

		// std::cout << "starting " << num_of_knots << " " << delta_t << " " << init_s[0] << " " << init_s[1] << " " << init_s[2] << " " << s_acc_weight << " " << s_jerk_weight << " " << num_of_obs << " \n";

		std::vector<std::pair<double, double>> x_bounds, y_bounds, dx_bounds, dy_bounds;
		std::vector<double> x_ref, y_ref, x_kappa, y_kappa;
		std::vector<std::vector<std::pair<double, double>>> all_x_bounds, all_y_bounds;
		init_s[1] = init_s[1] * mscale;

		double val[2];

		// std::cout << "lessgoo\n";
		PiecewiseJerkSpeedProblem piecewise_jerk_problem_0 = PiecewiseJerkSpeedProblem(num_of_knots, delta_t, init_s, init_l, num_of_obs);
		for (int obs = 0; obs < num_of_obs; obs++)
		{
			for (int i = 0; i < num_of_knots; i++)
			{
				// std::cout<<"\nx_bounds\n";
				ifs >> val[0] >> val[1];
				// cout<<val[1]<<" ";
				x_bounds.emplace_back(mscale * val[0], mscale * val[1]);
			}

			all_x_bounds.push_back(x_bounds);

			for (int i = 0; i < num_of_knots; i++)
			{
				// std::cout<<"\ny_bounds\n";
				ifs >> val[0] >> val[1];
				// cout<<val[1]<<" ";
				y_bounds.emplace_back(mscale * val[0], mscale * val[1]);
			}
			// cout<<"\n----------==========----------\n";
			all_y_bounds.push_back(y_bounds);

			x_bounds.clear();
			y_bounds.clear();
		}

		for (int i = 0; i < num_of_knots; i++)
		{
			ifs >> val[0] >> val[1];
			dx_bounds.emplace_back(mscale * val[0], mscale * val[1]);
		}

		for (int i = 0; i < num_of_knots; i++)
		{
			ifs >> val[0] >> val[1];
			// cout<<val[1];
			dy_bounds.emplace_back(mscale * val[0], mscale * val[1]);
		}

		cout<<"\nx ref:\n";
		for (int i = 0; i < num_of_knots; i++)
		{
			ifs >> val[0];
			cout<<val[0]<<" ";
			x_ref.emplace_back(mscale * val[0]);
		}
	

		cout<<"\ny ref:\n";
		for (int i = 0; i < num_of_knots; i++)
		{
			ifs >> val[0];
			cout<<val[0]<<" ";
			y_ref.emplace_back(mscale * val[0]);
		}
		cout<<"\n\n";

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
		piecewise_jerk_problem_0.set_dy_ref(weight_dl_ref, dl_ref * mscale); // try 0.0 * mscale

		piecewise_jerk_problem_0.set_weight_end(weight_end_s, weight_end_l);

		auto start_time = std::chrono::system_clock::now();

		// FORM NEW CORRIDORS
		for (int i = 0; i < num_of_obs; i++)
		{
			piecewise_jerk_problem_0.set_x_bounds(all_x_bounds[i]);
			piecewise_jerk_problem_0.set_y_bounds(all_y_bounds[i]);

			// std::cout << "generate corridor number: " << i << "\n";
			piecewise_jerk_problem_0.CorridorGeneration();
			piecewise_jerk_problem_0.PrintCorridor();
		}

		// piecewise_jerk_problem_0.FormNewCorridors();

		piecewise_jerk_problem_0.CollisionCheck();

		// auto start_time = std::chrono::system_clock::now();
		bool success = piecewise_jerk_problem_0.Optimize(5000);
		auto end_time = std::chrono::system_clock::now();
		std::chrono::duration<double> diff = end_time - start_time;
		ADEBUG << "Speed Optimizer used time: " << diff.count() * 1000 << " ms.";
		if (!success)
		{
			std::string msg("Piecewise jerk speed optimizer failed!");
			AERROR << msg;
			return 100000000000;
		}

		// Extract output
		const std::vector<double> &s = piecewise_jerk_problem_0.opt_x();
		const std::vector<double> &ds = piecewise_jerk_problem_0.opt_dx();
		const std::vector<double> &dds = piecewise_jerk_problem_0.opt_ddx();

		double s_avg_acc, l_avg_acc, s_max_acc, l_max_acc;
		s_avg_acc = l_avg_acc = s_max_acc = l_max_acc = 0.0;

		bool is_verbose = true;
		double s_cost = 0.0;
		double a_cost = 0.0;
		double mmax_a = 0.0;

		int num_of_points = s.size();

		for (int i = 0; i < num_of_points && is_verbose; ++i)
		{
			double ddds = (i == 0) ? (dds[1] - dds[0]) / delta_t : (dds[i] - dds[i - 1]) / delta_t;

			s_cost += weight_s_ref*(s[i] - x_ref[i]) * (s[i] - x_ref[i]) * delta_t;
			s_cost += weight_ds_ref * ds[i] * ds[i] * delta_t;
			s_cost += s_acc_weight* dds[i] * dds[i] * delta_t;
			s_cost += s_jerk_weight* ddds * ddds * delta_t;
			// s_cost += dds[i] * dds[i] * dds[i] * dds[i] * delta_t;
			// s_cost += ddds * ddds * ddds * ddds * delta_t;
			mmax_a = max(mmax_a, abs(dds[i]));

			s_avg_acc += dds[i];
			cout << std::fixed << std::setprecision(3) << "For t[" << i * delta_t << "], opt = " << s[i] << ", " << ds[i] << ", " << dds[i] << ", " << ddds << endl;
		}

		s_avg_acc /= num_of_points;
		// m->s_max_acc = mmax_a;

		// s_cost = s_acc_weight * sqrt(s_cost / (num_of_points * delta_t));
		// s_cost += mmax_a * mmax_a * mmax_a * mmax_a;

		cout<<"\ns_cost is \t"<<s_cost<<"\n";

		cout << "mmax_a " << mmax_a << endl;
		// cout << "a_cost " << a_cost << endl;

		cout << "\n\n\nL part of trajectory: \n\n";
		const std::vector<double> &l = piecewise_jerk_problem_0.opt_y();
		const std::vector<double> &dl = piecewise_jerk_problem_0.opt_dy();
		const std::vector<double> &ddl = piecewise_jerk_problem_0.opt_ddy();
		is_verbose = true;
		double l_cost=0.0;
		// a_cost = 0.0;
		mmax_a = 0.0;
		cout << "\n\n\nprinting L part: \n\nl size is\t" << l.size() << "\n\n";
		for (int i = 0; i < num_of_points && is_verbose; ++i)
		{
			double dddl = (i == 0) ? (ddl[1] - ddl[0]) / delta_t : (ddl[i] - ddl[i - 1]) / delta_t;

			l_cost += weight_l_ref*(l[i] - y_ref[i]) * (l[i] - y_ref[i]) * delta_t;
			l_cost += weight_dl_ref* dl[i] * dl[i] * delta_t;
			l_cost += l_acc_weight* ddl[i] * ddl[i] * delta_t;
			l_cost += l_jerk_weight* dddl * dddl * delta_t;
			mmax_a = max(mmax_a, abs(ddl[i]));

			l_avg_acc += ddl[i];
			cout << std::fixed << std::setprecision(3) << "For t[" << i * delta_t << "], opt = " << l[i] << ", " << dl[i] << ", " << ddl[i] << ", " << dddl << endl;
		}

		// cout<<"weight_end_l = "<<weight_end_l<<"\tdelta_t = "<<delta_t<<"\t"<<(l[num_of_knots-1] - y_ref[num_of_knots-1])*(l[num_of_knots-1] - y_ref[num_of_knots-1])<<"\n"<<"l_cost = "<<l_cost<<" + weight_end_l*(l[num_of_knots-1] - y_ref[num_of_knots-1]) * (l[num_of_knots-1] - y_ref[num_of_knots-1]) * delta_t = ";

		l_cost += weight_end_l*(l[num_of_knots-1] - y_ref[num_of_knots-1]) * (l[num_of_knots-1] - y_ref[num_of_knots-1]) * delta_t;
		// cout<<l_cost<<"\n\n";
		l_avg_acc /= num_of_points;
		// m->l_max_acc = mmax_a;
		// l_cost = l_acc_weight * sqrt(l_cost / (num_of_points * delta_t));
		// l_cost += mmax_a * mmax_a; // * mmax_a * mmax_a;

		cout<<"\nl_cost is \t"<<l_cost<<"\n";

		a_cost = s_cost + l_cost;
		cout << "mmax_a " << mmax_a << endl;
		cout << "a_cost " << a_cost << endl;

		// m->a_cost = a_cost;
		// m->l_avg_acc = l_avg_acc;
		// m->s_avg_acc = s_avg_acc;

		acc = a_cost;

		std::string file = "/home/srujan_d/RISS/code/btrapz/src/s1_slt_3d_";
		file += std::to_string(iteration) + ".txt";
		ofstream ofs(file);
		if (!ofs.is_open())
		{
			// cout << "Not Found" << endl;
			ofs.open(file.c_str(), ios_base::in | ios_base::out | ios_base::trunc);
			// return acc;
		}

		for (int i = 0; i < num_of_points && is_verbose; ++i)
		{
			ofs << std::fixed << std::setprecision(3) << i * delta_t << " " << s[i] << " " << l[i] << " " << ds[i] << " " << dl[i] << " " << dds[i] << " " << ddl[i] << endl;
		}
		// ADEBUG << "Speed Optimizer used time: " << diff.count() * 1000 << " ms.";

		return acc;
	}
}