#include <cmath>
#include <algorithm>
#include "../include/btrapz/solve_3d.h"
#include "../include/btrapz/logging.h"
#include <cmath>
#include <iostream>



using namespace Eigen;
using namespace std;

namespace navigation
{

	PiecewiseJerkSpeedProblem::PiecewiseJerkSpeedProblem(
		const size_t num_of_knots, const double delta_s,
		const std::array<double, 3> &x_init,
		const std::array<double, 3> &y_init,
		const int num_of_obs)
		: PiecewiseJerkProblem(num_of_knots, delta_s, x_init, y_init, num_of_obs)
	{

		penalty_dx_.resize(num_of_knots_, 0.0);

		n = static_cast<int>(num_of_knots_);

		traj_order = 5;
		n_poly = traj_order + 1;
		weight_end_ref_ = 10.0;
	}


	void PiecewiseJerkSpeedProblem::set_weight_end(const double weight_end_s, const double weight_end_l){
		weight_end_s_ = weight_end_s;
		weight_end_l_ = weight_end_l;
	}

	void PiecewiseJerkSpeedProblem::set_dx_ref(const double weight_dx_ref,
											   const double dx_ref)
	{
		weight_dx_ref_ = weight_dx_ref;
		dx_ref_ = dx_ref;
		has_dx_ref_ = true;
	}

	void PiecewiseJerkSpeedProblem::set_penalty_dx(std::vector<double> penalty_dx)
	{
		CHECK_EQ(penalty_dx.size(), num_of_knots_);
		penalty_dx_ = std::move(penalty_dx);
	}

	void PiecewiseJerkSpeedProblem::set_dy_ref(const double weight_dy_ref,
											   const double dy_ref)
	{
		weight_dy_ref_ = weight_dy_ref;
		dy_ref_ = dy_ref;
		has_dy_ref_ = true;
	}

	void PiecewiseJerkSpeedProblem::set_penalty_dy(std::vector<double> penalty_dy)
	{
		CHECK_EQ(penalty_dy.size(), num_of_knots_);
		penalty_dy_ = std::move(penalty_dy);
	}

	void PiecewiseJerkSpeedProblem::CalculateKernel(std::vector<c_float> *P_data,
													std::vector<c_int> *P_indices,
													std::vector<c_int> *P_indptr)
	{
		MatrixXd pQp[4];

		MatrixXd pQp_x[4];
		MatrixXd pQp_y[4];

		for (int k = 0; k < 4; k++)
		{
			pQp_x[k] = MatrixXd::Zero(n_poly, n_poly);
			pQp_y[k] = MatrixXd::Zero(n_poly, n_poly);
		}

		// calculating integral of < < differentiation (df/dt)^2 to l_th order > from 0 to 1 > before hand:

		for (int i = 0; i < n_poly; i++)
		{
			for (int j = 0; j < n_poly; j++)
			{
				pQp_x[0](i, j) = double(weight_x_ref_) / (i + j + 1);
				pQp_y[0](i, j) = double(weight_y_ref_) / (i + j + 1);

				if (i >= 1 && j >= 1)
				{
					//pQp_x[1](i,j) = double((weight_dx_ref_ + penalty_dx_[i]) * i * j) / (i + j - 1);
					pQp_x[1](i, j) = double(weight_dx_ref_ * i * j) / (i + j - 1);
					pQp_y[1](i, j) = double(weight_dy_ref_ * i * j) / (i + j - 1);
				}

				if (i >= 2 && j >= 2)
				{
					pQp_x[2](i, j) = double(weight_ddx_ * i * j * (i - 1) * (j - 1)) / (i + j - 3);
					pQp_y[2](i, j) = double(weight_ddy_ * i * j * (i - 1) * (j - 1)) / (i + j - 3);
				}

				if (i >= 3 && j >= 3)
				{
					pQp_x[3](i, j) = double(weight_dddx_ * i * j * (i - 1) * (j - 1) * (i - 2) * (j - 2)) / (i + j - 5);
					pQp_y[3](i, j) = double(weight_dddy_ * i * j * (i - 1) * (j - 1) * (i - 2) * (j - 2)) / (i + j - 5);
				}
			}
		}

		// cout<<"\nweight_x_ref_:\t"<<weight_x_ref_<<"\tweight_dx_ref_:\t"<<weight_dx_ref_<<"\tweight_ddx_:\t"<<weight_ddx_<<"\tweight_dddx_:\t"<<weight_dddx_<<"\n";
		// cout<<"\nweight_y_ref_:\t"<<weight_y_ref_<<"\tweight_dy_ref_:\t"<<weight_dy_ref_<<"\tweight_ddy_:\t"<<weight_ddy_<<"\tweight_dddy_:\t"<<weight_dddy_<<"\n";

		//MatrixXd MQM[4];

		M.resize(n_poly, n_poly);

		M << 1, 0, 0, 0, 0, 0,
			-5, 5, 0, 0, 0, 0,
			10, -20, 10, 0, 0, 0,
			-10, 30, -30, 10, 0, 0,
			5, -20, 30, -20, 5, 0,
			-1, 5, -10, 10, -5, 1;

		// for (int k = 0; k < 4; k++)
		// {
		// 	MQM_x[k] = MatrixXd::Zero(n_poly, n_poly);
		// 	MQM_y[k] = MatrixXd::Zero(n_poly, n_poly);
		// }

		// cout << "\ncalc kernel \n";

		for (int k = 0; k < 4; k++)
		{
			MQM_x[k] = M.transpose() * pQp_x[k] * M;
			MQM_y[k] = M.transpose() * pQp_y[k] * M;
			// cout<<"\nMQM X is\n\n"<<pQp_x[k]<<"\n";
			// cout<<"\nMQM Y is\n\n"<<pQp_y[k]<<"\n";
		}


		int idx = 0;
		int sub_shift = 0;

		for (int k = 0; k < segment_num; k++)
		{
			for (int j = 0; j < n_poly; j++)
			{
				P_indptr->push_back(idx);

				for (int i = 0; i < n_poly; i++)
				{
					if (j >= i)
					{

						double mval = pow(new_corridor[k].t, 3) * MQM_x[0](i, j) + new_corridor[k].t * MQM_x[1](i, j) +
									  MQM_x[2](i, j) / new_corridor[k].t + MQM_x[3](i, j) / pow(new_corridor[k].t, 3);

						P_indices->push_back(sub_shift + i);

						if ((k == segment_num - 1) && (i == n_poly - 1) && (j == n_poly - 1))
						{	
							// std::cout<<"\nFor last control point mval before is \t"<<mval<<"\t";
							mval = mval + weight_end_s_ * new_corridor[k].t * new_corridor[k].t;
							// std::cout<<"and mval after is \t"<<mval<<"\n";
						}
						P_data->push_back(2.0 * mval);
						idx++;
						// trp.push_back(Trip(j, i, 2.0 * mval));
					}
				}
			}
			sub_shift += n_poly;
		}
		int mid = idx;
		// P_indptr->push_back(idx);

		// cout << "\nmid ind ptr is " << idx << "\n";

		// idx = 0;
		// sub_shift = 0;

		// cout << "seg num is\t" << segment_num << "\n";
		for (int k = 0; k < segment_num; k++)
		{
			for (int j = 0; j < n_poly; j++)
			{
				P_indptr->push_back(idx);

				for (int i = 0; i < n_poly; i++)
				{
					if (j >= i)
					{
						double mval = pow(new_corridor[k].t, 3) * MQM_y[0](i, j) + new_corridor[k].t * MQM_y[1](i, j) +
									 MQM_y[2](i, j) / new_corridor[k].t + MQM_y[3](i, j) / pow(new_corridor[k].t, 3); // new_corridor[k].t is h_k

						P_indices->push_back(sub_shift + i);

						if ((k == segment_num - 1) && (i == n_poly - 1) && (j == n_poly - 1))
						{
							// std::cout<<"\nFor last control point mval before is \t"<<mval<<"\t";
							mval = mval + weight_end_l_ * new_corridor[k].t * new_corridor[k].t;
							// std::cout<<"and mval after is \t"<<mval<<"\n";							
						}
						P_data->push_back(2.0 * mval);
						idx++;

						// trp.push_back(Trip(j+54, i+54, 2.0 * mval));
						// mid++;
					}
				}
			}
			sub_shift += n_poly;
		}

		P_indptr->push_back(idx);
		// cout << "P data size is : \t" << P_data->size() << "\tmid is :\t" << mid << "\t and idx is\t" << idx << "\n";
	}

	void PiecewiseJerkSpeedProblem::CalculateOffset(std::vector<c_float> *q)
	{
		CHECK_NOTNULL(q);
		q->resize(segment_num * n_poly * 2, 0.0);

		MatrixXd q_p, q_c;

		q_p = MatrixXd::Zero(1, n_poly);

		int p_shift = 0;

		for (int k = 0; k < segment_num; k++)
		{
			// std::cout<<"seg num \t"<<k<<"\n";
			for (int i = 0; i < n_poly; i++)
			{
				// std::cout<<"index i is \t"<<i<<"\n";
				q_p(0, i) = 0.0;

				if (has_x_ref_)
				{

					q_p(0, i) += -2.0 * pow(new_corridor[k].t, 3) * weight_x_ref_ * x_skew_[k] / (i + 2);
					q_p(0, i) += -2.0 * pow(new_corridor[k].t, 2) * weight_x_ref_ * x_bias_[k] / (i + 1);
				}

				if (has_dx_ref_ && i > 0)
				{
					q_p(0, i) += -2.0 * weight_dx_ref_ * dx_ref_ * new_corridor[k].t;
				}
			}
			q_c = q_p * M;

			for (int i = 0; i < n_poly; i++)
			{
				q->at(p_shift + i) = double(q_c(0, i));
			}
			p_shift += n_poly;
		}

		q->at(segment_num * n_poly - 1) -= dx_ref_ * 2.0 * x_ref_[num_of_knots_ - 1] * new_corridor[segment_num - 1].t;

		p_shift = 0;
		// std::cout << "\n LLLL \n";
		for (int k = 0; k < segment_num; k++)
		{
			// std::cout<<"seg num \t"<<k<<"\n";
			for (int i = 0; i < n_poly; i++)
			{
				// std::cout<<"index i is \t"<< i - n_poly<<"\n";
				q_p(0, i) = 0.0;

				if (has_y_ref_)
				{
					q_p(0, i) += -2.0 * pow(new_corridor[k].t, 3) * weight_y_ref_ * y_skew_[k] / (i + 2);
					q_p(0, i) += -2.0 * pow(new_corridor[k].t, 2) * weight_y_ref_ * y_bias_[k] / (i + 1);

				}

				if (has_dy_ref_ && i > 0)
				{
					q_p(0, i) += -2.0 * weight_dy_ref_ * dy_ref_ * new_corridor[k].t;
				}
			}
			q_c = q_p * M;

			for (int i = 0; i < n_poly; i++)
			{
				// std::cout<<"q size is "<<q->size(); std::cout<<"\t and index is "<<p_shift<<"+"<<i<<"\n";
				q->at(segment_num * n_poly + p_shift + i) = double(q_c(0, i));
			}
			p_shift += n_poly;
		}

		q->at(segment_num * 2 * n_poly - 1) -= dy_ref_ * 2.0 * y_ref_[num_of_knots_ - 1] * new_corridor[segment_num - 1].t ;
		// std::cout << "q is done finally with size \t" << q->size() << "\n";
		// for (int i=0; i<q->size(); i++)	std::cout<<i<<"\t"<<q->at(i)<<"\n";
		// std::cout << "q is done finally with size \t" << q->size() << "\n";
	}

	void PiecewiseJerkSpeedProblem::CorridorGeneration()
		{	

			int j = 0;
			{
				Cube mcube;

				mcube.beg_t = 0;

				mcube.down_skew = (x_bounds_[1].first - x_bounds_[0].first) / delta_s_;
				mcube.down_bias = x_bounds_[0].first;
				mcube.upp_skew = (x_bounds_[1].second - x_bounds_[0].second) / delta_s_;
				mcube.upp_bias = x_bounds_[0].second;
				// mcube.beg_s =  
				mcube.beg_l = y_bounds_[0].first;
				mcube.end_l = y_bounds_[0].second;
				corridor.push_back(mcube);
				j++;
			}

			for (int i = 2; i < n - 1; i++)
			{
				Cube mcube;

				double uskew, dskew;
				dskew = (x_bounds_[i].first - x_bounds_[i - 1].first) / delta_s_;
				uskew = (x_bounds_[i].second - x_bounds_[i - 1].second) / delta_s_;

				// cout<<"dskew : "<<dskew<<"\t uskew : "<<uskew<<"\n";

				//threshold
				double mthre = 0.2;

				if ((abs(dskew - corridor[j - 1].down_skew) > mthre) || (abs(uskew - corridor[j - 1].upp_skew) > mthre))
				{
					corridor[j - 1].end_t = i;

					mcube.beg_t = i;

					mcube.down_skew = (x_bounds_[i + 1].first - x_bounds_[i].first) / delta_s_;
					//mcube.down_bias = 2 * x_bounds_[i].first - x_bounds_[i + 1].first;
					mcube.down_bias = x_bounds_[i].first;

					mcube.upp_skew = (x_bounds_[i + 1].second - x_bounds_[i].second) / delta_s_;
					//mcube.upp_bias = 2 * x_bounds_[i].second - x_bounds_[i + 1].second;
					mcube.upp_bias = x_bounds_[i].second;

					mcube.beg_l = y_bounds_[i].first;
					mcube.end_l = y_bounds_[i].second;

					// cout<<"y bounds are "<<mcube.beg_l<<" and "<<mcube.end_l<<"\n";

					corridor.push_back(mcube);

					j++;
				}
			}


			// auto max_bias = std::max_element( corridor.begin(), corridor.end(),
			// 							[]( const Cube &a, const Cube &b )
			// 							{
			// 								return a.down_bias < b.down_bias;
			// 							} ); 
			// auto min_bias = std::min_element( corridor.begin(), corridor.end(),
			// 							[]( const Cube &a, const Cube &b )
			// 							{
			// 								return a.upp_bias > b.upp_bias;
			// 							} );
			// auto max_skew = std::max_element( corridor.begin(), corridor.end(),
			// 							[]( const Cube &a, const Cube &b )
			// 							{
			// 								return a.down_skew < b.down_skew;
			// 							} ); 
			// auto min_skew = std::min_element( corridor.begin(), corridor.end(),
			// 							[]( const Cube &a, const Cube &b )
			// 							{
			// 								return a.upp_skew > b.upp_skew;
			// 							} );
										
			corridor[j - 1].end_t = n - 1;

			for (int i = 0; i < int(corridor.size()); i++)
			{
				corridor[i].t = (corridor[i].end_t - corridor[i].beg_t) * delta_s_;
				
				// for cuboidal corridors, COMMENT FOR TRAPEZOIDAL CORRIDORS
				// corridor[i].down_bias = max_bias(corridor.begin(), corridor.end());
				// corridor[i].upp_bias = min_bias(corridor.begin(), corridor.end());
				// corridor[i].down_skew = max_skew(corridor.begin(), corridor.end());
				// corridor[i].upp_skew = min_skew(corridor.begin(), corridor.end());
			}
			CorridorSplit();
			CorridorMerge();

			
			corridors.push_back(corridor);
			// std::cout << "\n size of corridors is " << corridor.size() << "\n";
			// for (int i = 0; i < corridor.size(); i++)
			// {
			// 	std::cout << corridor[i].beg_t << " " << corridor[i].end_t << " " << corridor[i].down_bias << " " << corridor[i].upp_bias << " " << corridor[i].down_skew << " " << corridor[i].upp_skew << "\n";
			// }

			// // for(int i=0; i<corridor.size(); i++)	std::cout << "\t"<<corridor[i].t<<"\n";
			// // corridors.push_back(corridor);
			// corridor.clear();
		}

	void PiecewiseJerkSpeedProblem::CollisionCheck()
	{
		waypoint p1;
		for (int i = 0; i < x_ref_.size(); i++)
		{
			p1.s = x_ref_[i];
			// cout<<x_ref_[i]<<"\t";
			p1.l = y_ref_[i];
			// cout<<y_ref_[i]<<"\n";
			p1.t = i;
			traj.push_back(p1);
		}

		int total_num_corridors = 0;
		for (int i = 0; i < corridors.size(); i++)
		{
			total_num_corridors += corridors[i].size();
		}
		// waypoint p1;
		// p1.s = 1.0/delta_s_;	p1.l = 2.0; p1.t = 1.5/0.1;
		// traj.push_back(p1);
		std::vector<Cube> temp;

		int pos, neg;
		double d;
		pos = neg = 0;
		int i, j, k;
		i = j = k = 0;

		// cout<<"\n\npriting all \n\n";
		// for(j=0; j<corridors.size(); j++){
		// 	for(k=0; k<corridors[j].size(); k++){
		// 		std::cout << corridors[j][k].beg_t << " " << corridors[j][k].end_t << " " << corridors[j][k].down_bias << " " << corridors[j][k].upp_bias << " " << corridors[j][k].down_bias + corridors[j][k].down_skew*delta_s_ << " " << corridors[j][k].upp_skew*delta_s_+corridors[j][k].upp_bias <<" "<<  corridors[j][k].beg_l <<" " << corridors[j][k].end_l<< "\n";
		// 	}
		// }
		cout << "\n\nNEXT\n\n";
		int count = 0;
		for (j = 0; j < corridors.size(); j++)
		{
			for (k = 0; k < corridors[j].size(); k++)
			{
				for (i = 0; i < traj.size(); i++)
				{
					pos = neg = 0;
					// cout<<"i\t"<<i<<"\tj\t"<<j<<"\tk\t"<<k<<"\n";
					// cout<<"\ntraj[i].l = "<<traj[i].l<<"\tcorridors[j][k].beg_l = "<<corridors[j][k].beg_l<<"\tcorridors[j][k].end_l = "<<corridors[j][k].end_l<<"\n";
					if (traj[i].l <= corridors[j][k].end_l && traj[i].l >= corridors[j][k].beg_l)
					{
						d = (traj[i].s - corridors[j][k].down_bias) * (corridors[j][k].beg_t - corridors[j][k].beg_t) - (traj[i].t - corridors[j][k].beg_t) * (corridors[j][k].upp_bias - corridors[j][k].down_bias);
						if (d > 0)
							pos++;
						if (d < 0)
							neg++;
						if (pos > 0 && neg > 0)
						{
							// cout<<"\n"<<d<<"\t\tcontinuing\n";
							continue;
						}

						d = (traj[i].s - corridors[j][k].upp_bias) * (corridors[j][k].end_t - corridors[j][k].beg_t) - (traj[i].t - corridors[j][k].beg_t) * (corridors[j][k].upp_skew * delta_s_ + corridors[j][k].upp_bias - corridors[j][k].upp_bias);
						// cout<<"\npoint2:\ns1:\t"<<corridors[j][k].upp_bias<<"\tt1:\t"<<corridors[j][k].beg_t<<"\ts2:\t"<<corridors[j][k].upp_skew*delta_s_+corridors[j][k].upp_bias<<"\tt2:\t"<<corridors[j][k].end_t;
						if (d > 0)
							pos++;
						if (d < 0)
							neg++;
						if (pos > 0 && neg > 0)
						{
							// cout<<"\ncontinuing\n";
							continue;
						}

						d = (traj[i].s - corridors[j][k].upp_bias - corridors[j][k].upp_skew * delta_s_) * (corridors[j][k].end_t - corridors[j][k].end_t) - (traj[i].t - corridors[j][k].end_t) * (corridors[j][k].down_skew * delta_s_ + corridors[j][k].down_bias - corridors[j][k].upp_skew * delta_s_ - corridors[j][k].upp_bias);
						// cout<<"\npoint3:\ns1:\t"<<corridors[j][k].upp_bias<<"\tt1:\t"<<corridors[j][k].end_t<<"\ts2:\t"<<corridors[j][k].down_skew*delta_s_+corridors[j][k].down_bias<<"\tt2:\t"<<corridors[j][k].end_t;
						if (d > 0)
							pos++;
						if (d < 0)
							neg++;
						if (pos > 0 && neg > 0)
						{
							// cout<<"\ncontinuing\n";
							continue;
						}

						d = (traj[i].s - corridors[j][k].down_bias - corridors[j][k].down_skew * delta_s_) * (corridors[j][k].beg_t - corridors[j][k].end_t) - (traj[i].t - corridors[j][k].end_t) * (corridors[j][k].down_bias - corridors[j][k].down_skew * delta_s_ - corridors[j][k].down_bias);
						// cout<<"\npoint4:\ns1:\t"<<corridors[j][k].down_bias + corridors[j][k].down_skew*delta_s_<<"\tt1:\t"<<corridors[j][k].end_t<<"\ts2:\t"<<corridors[j][k].down_bias<<"\tt2:\t"<<corridors[j][k].beg_t;
						if (d > 0)
							pos++;
						if (d < 0)
							neg++;
						if (pos > 0 && neg > 0)
						{
							// cout<<"\ncontinuing\n";
							continue;
						}

						// cout<<"\npushing corridor\t"<<temp.size()<<"\n";
						count++;
						if (count > 2){
							corridors[j][k].count = count;
							temp.push_back(corridors[j][k]);
							count = 0;
							// k++;
						}
						else{
							corridors[j][k].count = count;
							// std::cout << corridors[j][k].beg_t << " " << corridors[j][k].end_t << " " << corridors[j][k].down_bias << " " << corridors[j][k].upp_bias << " " << corridors[j][k].down_bias + corridors[j][k].down_skew * delta_s_ << " " << corridors[j][k].upp_skew * delta_s_ + corridors[j][k].upp_bias << " " << corridors[j][k].beg_l << " " << corridors[j][k].end_l << "\t count is "<< corridors[j][k].count << "\n";
						}
						if (count > 2) 
							k++;
						// cout<<"\ndone, next\n";
						// continue;
						// }
					}
					if (j == corridors.size() || k == corridors[j].size()){
						count = 0;
						break;
					}
						
				}
				if (j == corridors.size()){
					count = 0;
					break;
				}
			}
		}
		// cout << "HEYEYYEYEYE\n\n"<< temp.size() << "\n\n";
		for (int i = 0; i < temp.size(); i++)
			std::cout << temp[i].beg_t << " " << temp[i].end_t << " " << temp[i].down_bias << " " << temp[i].upp_bias << " " << temp[i].down_bias + temp[i].down_skew * delta_s_ << " " << temp[i].upp_skew * delta_s_ + temp[i].upp_bias << " " << temp[i].beg_l << " " << temp[i].end_l << "\t count is "<< temp[i].count << "\n";

		for (int i = 0; i < temp.size() - 1; i++)
		{
			for (int j = i + 1; j < temp.size(); j++)
			{
				if (temp[i].beg_t == temp[j].beg_t && temp[i].end_t == temp[j].end_t && temp[i].down_bias == temp[j].down_bias && temp[i].down_skew == temp[j].down_skew && temp[i].upp_bias == temp[j].upp_bias && temp[i].upp_skew == temp[j].upp_skew && temp[i].beg_l == temp[j].beg_l && temp[i].end_l == temp[j].end_l)
				{	
					// std::cout << "erasing:\t"<<temp[j].beg_t << " " << temp[j].end_t << " " << temp[j].down_bias << " " << temp[j].upp_bias << " " << temp[j].down_bias + temp[j].down_skew * delta_s_ << " " << temp[j].upp_skew * delta_s_ + temp[j].upp_bias << " " << temp[j].beg_l << " " << temp[j].end_l << "\t count is "<< temp[j].count << "\n";
					temp.erase(temp.begin()+j);
					j--;
				}
			}
		}

		cout << "\n\n new corridors are \n\n";

		for (int i = 0; i < temp.size() - 1; i++)
		{
			for (int j = i + 1; j < temp.size(); j++)
			{	
				if (temp[i].beg_t == temp[j].beg_t && temp[i].end_t == temp[j].end_t)
				{
					int diff = (temp[i].end_t - temp[i].beg_t) / 3;
					temp[i].end_t = temp[i].end_t - (diff);
					temp[i].t = (temp[i].end_t - temp[i].beg_t) * delta_s_;

					temp[j].beg_t = temp[j].beg_t + (diff);
					temp[j].t = (temp[j].end_t - temp[j].beg_t) * delta_s_;
				}
			}
		}
		for (int i = 0; i < temp.size(); i++)
			std::cout << temp[i].beg_t << " " << temp[i].end_t << " " << temp[i].down_bias << " " << temp[i].upp_bias << " " << temp[i].down_bias + temp[i].down_skew * delta_s_ << " " << temp[i].upp_skew * delta_s_ + temp[i].upp_bias << " " << temp[i].beg_l << " " << temp[i].end_l << "\n";

		// cout << "HEYEYYEYEYE\n\n";
		new_corridor = std::move(temp);
	}

	void PiecewiseJerkSpeedProblem::PrintCorridor()
	{	
		// std::cout << "\n size of corridors is " << corridor.size() << "\n";
		// for (int i = 0; i < corridor.size(); i++)
		// {
		// 	std::cout << corridor[i].beg_t << " " << corridor[i].end_t << " " << corridor[i].down_bias << " " << corridor[i].upp_bias << " " << corridor[i].down_bias + corridor[i].down_skew*delta_s_ << " " << corridor[i].upp_skew*delta_s_+corridor[i].upp_bias <<" "<<  corridor[i].beg_l <<" " << corridor[i].end_l<< "\n";
		// }

		// for(int i=0; i<corridor.size(); i++)	std::cout << "\t"<<corridor[i].t<<"\n";
		// corridors.push_back(corridor);
		corridor.clear();
	}

	void PiecewiseJerkSpeedProblem::CorridorSplit()
	{
		int temp_num = int(corridor.size());

		for (int k = 0; k < temp_num; k++)
		{
			while (corridor[k].t > 1)
			{
				corridor[k].t = corridor[k].t - 1;

				Cube mcube;
				// construct
				mcube.beg_t = corridor[k].beg_t;

				corridor[k].beg_t = corridor[k].beg_t + 10;

				mcube.end_t = mcube.beg_t + 10;
				mcube.t = 1.0;
				mcube.down_skew = corridor[k].down_skew;
				mcube.down_bias = corridor[k].down_bias;

				corridor[k].down_bias = mcube.down_bias + 1.0 * mcube.down_skew;

				mcube.upp_skew = corridor[k].upp_skew;
				mcube.upp_bias = corridor[k].upp_bias;

				mcube.beg_l = corridor[k].beg_l;
				mcube.end_l = corridor[k].end_l;

				corridor[k].upp_bias = mcube.upp_bias + 1.0 * mcube.upp_skew;

				// insert(mcube)
				corridor.insert(corridor.begin() + k, mcube);
				temp_num++;
				k++;
			}
		}
	}

	void PiecewiseJerkSpeedProblem::CorridorMerge()
	{
		return;
	}

	void PiecewiseJerkSpeedProblem::CalculateAffineConstraint(
		std::vector<c_float> *A_data, std::vector<c_int> *A_indices,
		std::vector<c_int> *A_indptr, std::vector<c_float> *lower_bounds,
		std::vector<c_float> *upper_bounds)
	{

		int num_of_constraints = 2 * (segment_num * (3 * n_poly - 3 + n_poly - 3) + 3 + 3 * (segment_num - 1));

		int num_of_variables = segment_num * n_poly;

		vector<vector<pair<c_int, c_float>>> variables(segment_num * n_poly * 2);

		double aval[3];
		aval[0] = 1.0 * traj_order * (traj_order - 1);
		aval[1] = -2.0 * traj_order * (traj_order - 1);
		aval[2] = 1.0 * traj_order * (traj_order - 1);

		double jval[4];
		jval[0] = -1.0 * traj_order * (traj_order - 1) * (traj_order - 2);
		jval[1] = 3.0 * traj_order * (traj_order - 1) * (traj_order - 2);
		jval[2] = -3.0 * traj_order * (traj_order - 1) * (traj_order - 2);
		jval[3] = 1.0 * traj_order * (traj_order - 1) * (traj_order - 2);

		lower_bounds->resize(num_of_constraints);
		upper_bounds->resize(num_of_constraints);

		int var_shift = 0;

		int constraint_index = 0;

		MatrixXd inv_M;

		inv_M = MatrixXd::Zero(n_poly, n_poly);

		inv_M = M.inverse();

		double l_bound, u_bound;

		// cout << "\n S portion of constraints \n";

		for (int k = 0; k < segment_num; k++) // for each segment
		{
			// cout << "segment number \t" << k << "\n";
			// SAFETY CONSTRAINTS for trajectory:
			// cout << "\n\n-----------------safety constraints for corridor number "<<k<<"\t----------------------\n\n";
			l_bound=0;
			u_bound=100;
			for (int i = 0; i < n_poly; i++){
				l_bound = max(l_bound, new_corridor[k].down_bias + new_corridor[k].down_skew * inv_M(i, 1) * new_corridor[k].t);
				u_bound = min(u_bound, new_corridor[k].upp_bias + new_corridor[k].upp_skew * inv_M(i, 1) * new_corridor[k].t);
			}
				
			for (int i = 0; i < n_poly; i++) // bounds for each control point
			{
				variables[var_shift + i].emplace_back(constraint_index, 1.0 * new_corridor[k].t);
				
				lower_bounds->at(constraint_index) = l_bound; 
				upper_bounds->at(constraint_index) = u_bound; 

				// lower_bounds->at(constraint_index) = new_corridor[k].down_bias + new_corridor[k].down_skew * inv_M(i, 1) * new_corridor[k].t;
				// upper_bounds->at(constraint_index) = new_corridor[k].upp_bias + new_corridor[k].upp_skew * inv_M(i, 1) * new_corridor[k].t;
				// cout<<"for control point C_"<<i<<" \tvarshift = "<<var_shift<<" \t segment number = "<<k<<" \tconstraint index = "<<constraint_index<<"\t lower bounds = "<<lower_bounds->at(constraint_index)<<"\tupper bounds = "<<upper_bounds->at(constraint_index)<<"\n";
				// cout<<"start time:\t"<<new_corridor[k].beg_t<<"\tend time:\t"<<new_corridor[k].end_t<<"\t lower bound is : "<< lower_bounds->at(constraint_index)<<"\t and upper bound is: "<<upper_bounds->at(constraint_index)<<"\n";
				++constraint_index;

			}

			// PHYSICAL CONSTRAINTS :
			double dx_lower_bound = 0.0, dx_upper_bound = 1000.0;
			double ddx_lower_bound = -1000.0, ddx_upper_bound = 1000.0;

			for (int i = new_corridor[k].beg_t; i <= new_corridor[k].end_t; i++)
			{
				dx_lower_bound = max(dx_bounds_[i].first, dx_lower_bound);
				dx_upper_bound = min(dx_bounds_[i].second, dx_upper_bound);

				ddx_lower_bound = max(ddx_bounds_[i].first, ddx_lower_bound);
				ddx_upper_bound = min(ddx_bounds_[i].second, ddx_upper_bound);
			}

			// cout<<"\n\n----------------dx--------------\n\n";
			for (int i = 0; i < n_poly - 1; i++) // for 1st order trajectory (velocity)
			{
				// what are these two lines for? which constraints do they represent?
				variables[var_shift + i].emplace_back(constraint_index, -1.0 * traj_order);
				variables[var_shift + i + 1].emplace_back(constraint_index, 1.0 * traj_order);

				lower_bounds->at(constraint_index) = dx_lower_bound;
				upper_bounds->at(constraint_index) = dx_upper_bound;
		
				++constraint_index;
			}

			// cout<<"\n\n----------------ddx--------------\n\n";
			for (int i = 0; i < n_poly - 2; i++) // for 2nd order trajectory (acc)
			{
				variables[var_shift + i].emplace_back(constraint_index, aval[0]);
				variables[var_shift + i + 1].emplace_back(constraint_index, aval[1]);
				variables[var_shift + i + 2].emplace_back(constraint_index, aval[2]);

				lower_bounds->at(constraint_index) = ddx_lower_bound * new_corridor[k].t;
				upper_bounds->at(constraint_index) = ddx_upper_bound * new_corridor[k].t;
		
				++constraint_index;
			}

			// cout<<"\n\n----------------dddx--------------\n\n";
			for (int i = 0; i < n_poly - 3; i++) // for 3rd order trajectory (jerk)
			{
				variables[var_shift + i].emplace_back(constraint_index, jval[0]);
				variables[var_shift + i + 1].emplace_back(constraint_index, jval[1]);
				variables[var_shift + i + 2].emplace_back(constraint_index, jval[2]);
				variables[var_shift + i + 3].emplace_back(constraint_index, jval[3]);

				lower_bounds->at(constraint_index) = dddx_bound_.first * new_corridor[k].t * new_corridor[k].t;
				upper_bounds->at(constraint_index) = dddx_bound_.second * new_corridor[k].t * new_corridor[k].t;
		
				++constraint_index;
			}

			var_shift += n_poly;
		}

		//init_data << " initial_x = " << x_init_[0] << "initial_v = " << x_init_[1] << "initial_a = " << x_init_[2] << endl;
		// We are here.
		/*   Start position  */
		variables[0].emplace_back(constraint_index, 1.0 * new_corridor[0].t);
		lower_bounds->at(constraint_index) = x_init_[0];
		upper_bounds->at(constraint_index) = x_init_[0];
		++constraint_index;

		variables[0].emplace_back(constraint_index, -1.0 * traj_order);
		variables[1].emplace_back(constraint_index, 1.0 * traj_order);
		lower_bounds->at(constraint_index) = x_init_[1];
		upper_bounds->at(constraint_index) = x_init_[1];
		++constraint_index;

		variables[0].emplace_back(constraint_index, aval[0]);
		variables[1].emplace_back(constraint_index, aval[1]);
		variables[2].emplace_back(constraint_index, aval[2]);
		lower_bounds->at(constraint_index) = x_init_[2] * new_corridor[0].t;
		upper_bounds->at(constraint_index) = x_init_[2] * new_corridor[0].t;
		++constraint_index;

		/*   joint points  */
		// BOUNDARY CONSTRAINTS: (?)
		int sub_shift = 0;

		for (int k = 0; k < (segment_num - 1); k++)
		{

			sub_shift = (k + 1) * n_poly;
			// cout<<"subshift value = "<<sub_shift<<"\t\n";
			variables[sub_shift - 1].emplace_back(constraint_index, -1.0 * new_corridor[k].t);
			variables[sub_shift].emplace_back(constraint_index, 1.0 * new_corridor[k + 1].t);
			// cout<<"constraint index"<<constraint_index<<"\tsubshift -1 = "<<(sub_shift-1)<<"\tconstraint val = "<<(- 1.0 * new_corridor[k].t)<<"\tsubshift"<<sub_shift<<"\tconstraint val = "<<(1.0 * new_corridor[k + 1].t)<<"\n";
			lower_bounds->at(constraint_index) = 0.0;
			upper_bounds->at(constraint_index) = 0.0;
			++constraint_index;

			variables[sub_shift - 2].emplace_back(constraint_index, -1.0);
			variables[sub_shift - 1].emplace_back(constraint_index, 1.0);
			variables[sub_shift].emplace_back(constraint_index, 1.0);
			variables[sub_shift + 1].emplace_back(constraint_index, -1.0);

			lower_bounds->at(constraint_index) = 0.0;
			upper_bounds->at(constraint_index) = 0.0;
			constraint_index++;

			variables[sub_shift - 3].emplace_back(constraint_index, 1.0 * new_corridor[k + 1].t);
			variables[sub_shift - 2].emplace_back(constraint_index, -2.0 * new_corridor[k + 1].t);
			variables[sub_shift - 1].emplace_back(constraint_index, 1.0 * new_corridor[k + 1].t);
			variables[sub_shift].emplace_back(constraint_index, -1.0 * new_corridor[k].t);
			variables[sub_shift + 1].emplace_back(constraint_index, 2.0 * new_corridor[k].t);
			variables[sub_shift + 2].emplace_back(constraint_index, -1.0 * new_corridor[k].t);

			lower_bounds->at(constraint_index) = 0.0;
			upper_bounds->at(constraint_index) = 0.0;
			constraint_index++;
		}

		// L segment

		var_shift = 0;

		// cout << "\n L portion of constraints \n";
		// cout << "size of var vector; \t" << variables.size() << "\t and constraint index is: \t" << constraint_index << "\t " << num_of_variables << "\n";
		for (int k = 0; k < segment_num; k++) // for each segment
		{
			// SAFETY CONSTRAINTS for trajectory:
			// cout << "\n\n-----------------safety constraints for corridor number "<<k<<"\t----------------------\n\n";
			for (int i = 0; i < n_poly; i++) // bounds for each control point
			{
				variables[var_shift + i + num_of_variables].emplace_back(constraint_index, 1.0 * new_corridor[k].t);
				lower_bounds->at(constraint_index) = new_corridor[k].beg_l;
				upper_bounds->at(constraint_index) = new_corridor[k].end_l;
				// cout<<"y lower bound is : "<< lower_bounds->at(constraint_index)<<"\t and upper bound is: "<<upper_bounds->at(constraint_index)<<"\t and variable at index \t"<<var_shift + i + num_of_variables<<"\tis \t"<<1.0 * new_corridor[k].t<<"\n";
				// cout<<"for control point C_"<<i<<" \tvarshift = "<<var_shift<<" \t segment number = "<<k<<" \tconstraint index = "<<constraint_index<<"\t lower bounds = "<<lower_bounds->at(constraint_index)<<"\tupper bounds = "<<upper_bounds->at(constraint_index)<<"\n";
				++constraint_index;
			}

			// PHYSICAL CONSTRAINTS :
			// cout << "\n\n-----------------physical constraints----------------------\n\n";
			// double dy_lower_bound = 0.0, dy_upper_bound = 1000.0;
			// double ddy_lower_bound = -1000.0, ddy_upper_bound = 1000.0;

			double dy_lower_bound = 0.0, dy_upper_bound = 0.5;
			double ddy_lower_bound = -2.5, ddy_upper_bound = 1.5;

			// for (int i = new_corridor[k].beg_t; i <= new_corridor[k].end_t; i++)
			// {
			// 	cout << "size of dy bounds \t" << ddy_bounds_.size() << " and i is \t" << new_corridor[k].end_t << " " << new_corridor[k].beg_t << "\n";
			// 	dy_lower_bound = max(dy_bounds_[i].first, dy_lower_bound);
			// 	dy_upper_bound = min(dy_bounds_[i].second, dy_upper_bound);

			// 	ddy_lower_bound = max(ddy_bounds_[i].first, ddy_lower_bound);
			// 	ddy_upper_bound = min(ddy_bounds_[i].second, ddy_upper_bound);
			// }

			// cout << "\n\n----------------dy--------------\n\n";
			for (int i = 0; i < n_poly - 1; i++) // for 1st order trajectory (velocity)
			{
				// what are these two lines for? which constraints do they represent?
				variables[var_shift + i + num_of_variables].emplace_back(constraint_index, -1.0 * traj_order);
				variables[var_shift + i + 1 + num_of_variables].emplace_back(constraint_index, 1.0 * traj_order);

				// lower_bounds->at(constraint_index) = dy_lower_bound;
				// upper_bounds->at(constraint_index) = dy_upper_bound;

				lower_bounds->at(constraint_index) = dy_bounds_[i].first;
				upper_bounds->at(constraint_index) = dy_bounds_[i].second;
				// cout<<"lower bound is : "<< lower_bounds->at(constraint_index)<<"\t and upper bound is: "<<upper_bounds->at(constraint_index)<<"\n";
				++constraint_index;
			}

			// cout << "\n\n----------------ddy--------------\n\n";
			for (int i = 0; i < n_poly - 2; i++) // for 2nd order trajectory (acc)
			{
				variables[var_shift + i + num_of_variables].emplace_back(constraint_index, aval[0]);
				variables[var_shift + i + 1 + num_of_variables].emplace_back(constraint_index, aval[1]);
				variables[var_shift + i + 2 + num_of_variables].emplace_back(constraint_index, aval[2]);

				// lower_bounds->at(constraint_index) = ddy_lower_bound * new_corridor[k].t;
				// upper_bounds->at(constraint_index) = ddy_upper_bound * new_corridor[k].t;
				
				lower_bounds->at(constraint_index) = ddy_bounds_[i].first * new_corridor[k].t;
				upper_bounds->at(constraint_index) = ddy_bounds_[i].second * new_corridor[k].t;
				
				// cout<<"lower bound is : "<< lower_bounds->at(constraint_index)<<"\t and upper bound is: "<<upper_bounds->at(constraint_index)<<"\n";
				++constraint_index;
			}

			// cout << "\n\n----------------dddy--------------\n\n";
			for (int i = 0; i < n_poly - 3; i++) // for 3rd order trajectory (jerk)
			{
				variables[var_shift + i + num_of_variables].emplace_back(constraint_index, jval[0]);
				variables[var_shift + i + 1 + num_of_variables].emplace_back(constraint_index, jval[1]);
				variables[var_shift + i + 2 + num_of_variables].emplace_back(constraint_index, jval[2]);
				variables[var_shift + i + 3 + num_of_variables].emplace_back(constraint_index, jval[3]);

				lower_bounds->at(constraint_index) = dddy_bound_.first * new_corridor[k].t * new_corridor[k].t;
				upper_bounds->at(constraint_index) = dddy_bound_.second * new_corridor[k].t * new_corridor[k].t;
				// cout<<"lower bound is : "<< lower_bounds->at(constraint_index)<<"\t and upper bound is: "<<upper_bounds->at(constraint_index)<<"\n";
				++constraint_index;
			}

			var_shift += n_poly;
		}

		//init_data << " initial_x = " << x_init_[0] << "initial_v = " << x_init_[1] << "initial_a = " << x_init_[2] << endl;
		// We are here.
		/*   Start position  */
		// cout << "\n\n-----------------starting state constraints----------------------\n\n";
		variables[0 + num_of_variables].emplace_back(constraint_index, 1.0 * new_corridor[0].t);
		lower_bounds->at(constraint_index) = y_init_[0];
		upper_bounds->at(constraint_index) = y_init_[0];
		++constraint_index;

		variables[0 + num_of_variables].emplace_back(constraint_index, -1.0 * traj_order);
		variables[1 + num_of_variables].emplace_back(constraint_index, 1.0 * traj_order);
		lower_bounds->at(constraint_index) = y_init_[1];
		upper_bounds->at(constraint_index) = y_init_[1];
		++constraint_index;

		variables[0 + num_of_variables].emplace_back(constraint_index, aval[0]);
		variables[1 + num_of_variables].emplace_back(constraint_index, aval[1]);
		variables[2 + num_of_variables].emplace_back(constraint_index, aval[2]);
		lower_bounds->at(constraint_index) = y_init_[2] * new_corridor[0].t;
		upper_bounds->at(constraint_index) = y_init_[2] * new_corridor[0].t;
		++constraint_index;

		/*   joint points  */
		// BOUNDARY CONSTRAINTS: (?)
		sub_shift = 0;
		// int bruh = constraint_index;
		// cout << "\n\n-----------------boundary constraints----------------------\n\n";
		for (int k = 0; k < (segment_num - 1); k++)
		{

			sub_shift = (k + 1) * n_poly;
			// cout<<"subshift value = "<<sub_shift<<"\t\n";
			variables[sub_shift - 1 + num_of_variables].emplace_back(constraint_index, -1.0 * new_corridor[k].t);
			variables[sub_shift + num_of_variables].emplace_back(constraint_index, 1.0 * new_corridor[k + 1].t);
			// cout<<"constraint index"<<constraint_index<<"\tsubshift -1 = "<<(sub_shift-1)<<"\tconstraint val = "<<(- 1.0 * new_corridor[k].t)<<"\tsubshift"<<sub_shift<<"\tconstraint val = "<<(1.0 * new_corridor[k + 1].t)<<"\n";
			lower_bounds->at(constraint_index) = 0.0;
			upper_bounds->at(constraint_index) = 0.0;
			++constraint_index;

			variables[sub_shift - 2 + num_of_variables].emplace_back(constraint_index, -1.0);
			variables[sub_shift - 1 + num_of_variables].emplace_back(constraint_index, 1.0);
			variables[sub_shift + num_of_variables].emplace_back(constraint_index, 1.0);
			variables[sub_shift + 1 + num_of_variables].emplace_back(constraint_index, -1.0);

			lower_bounds->at(constraint_index) = 0.0;
			upper_bounds->at(constraint_index) = 0.0;
			constraint_index++;

			variables[sub_shift - 3 + num_of_variables].emplace_back(constraint_index, 1.0 * new_corridor[k + 1].t);
			variables[sub_shift - 2 + num_of_variables].emplace_back(constraint_index, -2.0 * new_corridor[k + 1].t);
			variables[sub_shift - 1 + num_of_variables].emplace_back(constraint_index, 1.0 * new_corridor[k + 1].t);
			variables[sub_shift + num_of_variables].emplace_back(constraint_index, -1.0 * new_corridor[k].t);
			variables[sub_shift + 1 + num_of_variables].emplace_back(constraint_index, 2.0 * new_corridor[k].t);
			variables[sub_shift + 2 + num_of_variables].emplace_back(constraint_index, -1.0 * new_corridor[k].t);

			lower_bounds->at(constraint_index) = 0.0;
			upper_bounds->at(constraint_index) = 0.0;
			constraint_index++;
		}

		// cout<<"diff is \t"<<constraint_index-bruh<<"\n";

		// cout << "size of var vector; \t" << variables.size() << "\t and index is: \t" << var_shift << " " << num_of_variables << "\n";

		// cout << "constraint index is \t" << constraint_index << " and num of constraints are \t" << num_of_constraints << "\n";
		CHECK_EQ(constraint_index, num_of_constraints);

		int ind_p = 0;
		for (int i = 0; i < 2*num_of_variables; ++i)
		{
			A_indptr->push_back(ind_p);
			for (const auto &variable_nz : variables[i])
			{
				// coefficient
				A_data->push_back(variable_nz.second);

				// constraint index
				A_indices->push_back(variable_nz.first);

				++ind_p;
				// trp.push_back(Trip(variable_nz.first, i, variable_nz.second));
			}
		}

		// We indeed need this line because of
		// https://github.com/oxfordcontrol/osqp/blob/master/src/cs.c#L255
		A_indptr->push_back(ind_p);
	}

	void PiecewiseJerkSpeedProblem::CorridorVisualize()
	{
		cout << "segment_num: " << segment_num << endl;
		for (int k = 0; k < segment_num; k++)
		{
			cout << k << "-th " << new_corridor[k].t;
			cout << " upp " << new_corridor[k].upp_skew << " " << new_corridor[k].upp_bias << " " << new_corridor[k].end_l;
			cout << " down " << new_corridor[k].down_skew << " " << new_corridor[k].down_bias << " " << new_corridor[k].beg_l;
			cout << endl;
		}
	}

	OSQPData *PiecewiseJerkSpeedProblem::FormulateProblem()
	{
		// calculate kernel
		// CorridorGeneration();

		segment_num = int(new_corridor.size());

		CorridorVisualize();

		x_skew_.resize(segment_num);
		x_bias_.resize(segment_num);

		y_skew_.resize(segment_num);
		y_bias_.resize(segment_num);


		// x_skew and x_bias are defined here
		for (int k = 0; k < segment_num; k++)
		{
			x_skew_[k] = (x_ref_[k * 10 + 1] - x_ref_[k * 10]) / delta_s_;
			x_bias_[k] = x_ref_[k * 10];

			y_skew_[k] = (y_ref_[k * 10 + 1] - y_ref_[k * 10]) / delta_s_;
			y_bias_[k] = y_ref_[k * 10];
		}

		std::vector<c_float> P_data;
		std::vector<c_int> P_indices;
		std::vector<c_int> P_indptr;
		//load data
		CalculateKernel(&P_data, &P_indices, &P_indptr);

		// int rows, cols;
		// rows = cols = segment_num*n_poly*2;
		// Eigen::SparseMatrix<double> temp(rows,cols);

		// temp.setFromTriplets(trp.begin(), trp.end());

		// std::cout <<"sparse to dense\n"<< temp << std::endl;


		// calculate affine constraints
		std::vector<c_float> A_data;
		std::vector<c_int> A_indices;
		std::vector<c_int> A_indptr;
		std::vector<c_float> lower_bounds;
		std::vector<c_float> upper_bounds;
		CalculateAffineConstraint(&A_data, &A_indices, &A_indptr, &lower_bounds, &upper_bounds);

		// int rows, cols;
		// rows = 2 * (segment_num * (3 * n_poly - 3 + n_poly - 3) + 3 + 3 * (segment_num - 1));
		// cols = segment_num*n_poly*2;
		// Eigen::SparseMatrix<double> temp(rows,cols);

		// temp.setFromTriplets(trp.begin(), trp.end());

		// std::cout <<"sparse to dense\n"<< temp << std::endl;

		// for (int i =0; i< lower_bounds.size(); i++){
		// 	std::cout<<"\nlower bounds\t"<<lower_bounds[i]<<"\tupper bounds\t"<<upper_bounds[i]<<"\n";
		// }

		// std::cout<<"\nlower bounds size\t"<<lower_bounds.size()<<"\tupper bounds size\t"<<upper_bounds.size()<<"\n";
		
		// calculate offset
		//std::vector<c_float> q;
		CalculateOffset(&q);

		OSQPData *data = reinterpret_cast<OSQPData *>(c_malloc(sizeof(OSQPData)));
		CHECK_EQ(lower_bounds.size(), upper_bounds.size());

		size_t kernel_dim;
		kernel_dim = segment_num * n_poly * 2;
		size_t num_affine_constraint = lower_bounds.size();

		data->n = kernel_dim;
		data->m = num_affine_constraint;
		data->P = csc_matrix(kernel_dim, kernel_dim, P_data.size(), CopyData(P_data), CopyData(P_indices), CopyData(P_indptr));
		// cout<<"\nP is: \n"<<data->P<<"\n";
		data->q = CopyData(q);
		// cout<<"\nq is: \n"<<data->q<<"\n";
		data->A = csc_matrix(num_affine_constraint, kernel_dim, A_data.size(), CopyData(A_data), CopyData(A_indices), CopyData(A_indptr));
		// cout<<"\nA is: \n"<<data->A<<"\n";
		data->l = CopyData(lower_bounds);
		data->u = CopyData(upper_bounds);
		return data;
	}

	bool PiecewiseJerkSpeedProblem::Optimize(const int max_iter)
	{
		OSQPData *data = FormulateProblem();

		OSQPSettings *settings = SolverDefaultSettings();
		settings->max_iter = max_iter;
		cout<<"\neps rel is \t"<<settings->eps_rel<<"\t and eps vel is \t"<<settings->scaling<<"\n";
		settings->eps_rel = 1e-5;
		settings->eps_abs = 1e-5;
		// settings->adaptive_rho_tolerance = 2;
		// settings->adaptive_rho = 0;
		settings->scaling = 4;
		settings->polish = 0;

		OSQPWorkspace *osqp_work = nullptr;
		osqp_work = osqp_setup(data, settings); // osqp-0.4.1, 0.5.0
		// osqp_setup(&osqp_work, data, settings); // osqp-0.6.0

		osqp_solve(osqp_work);

		auto status = osqp_work->info->status_val;

		if (status < 0 || (status != 1 && status != 2))
		{
			//AERROR << "failed optimization status:\t" << osqp_work->info->status;
			osqp_cleanup(osqp_work);
			FreeData(data);
			c_free(settings);
			return false;
		}
		else if (osqp_work->solution == nullptr)
		{
			//AERROR << "The solution from OSQP is nullptr";
			osqp_cleanup(osqp_work);
			FreeData(data);
			c_free(settings);
			return false;
		}

		
		for(int i=0; i<new_corridor.size(); i++){
			num_of_points_ += new_corridor[i].t/delta_s_ ;
		}
		cout<<"\n number of points are\t"<<num_of_points_<<"\n";
		x_.resize(num_of_points_);
		dx_.resize(num_of_points_);
		ddx_.resize(num_of_points_);

		y_.resize(num_of_points_);
		dy_.resize(num_of_points_);
		ddy_.resize(num_of_points_);

		auto mcost = osqp_work->info->obj_val;
		//auto mcost = osqp_work->info->pri_res;
		cout << " mcost: " << mcost << std::endl;

		double factorial[6];
		factorial[0] = 1.0;

		for (int i = 1; i < n_poly; i++)
		{
			factorial[i] = factorial[i - 1] * i;
		}
		double b_coe[6][3];

		for (int i = 0; i < n_poly; i++)
		{
			b_coe[i][0] = double(factorial[traj_order]) / (factorial[i] * factorial[traj_order - i]);
		}

		for (int i = 0; i < n_poly - 1; i++)
		{
			b_coe[i][1] = double(factorial[traj_order - 1]) / (factorial[i] * factorial[traj_order - 1 - i]);
		}

		for (int i = 0; i < n_poly - 2; i++)
		{
			b_coe[i][2] = double(factorial[traj_order - 2]) / (factorial[i] * factorial[traj_order - 2 - i]);
		}

		int sub_shift = 0;
		int var_index = 0;

		x_.at(var_index) = x_init_[0];
		dx_.at(var_index) = x_init_[1];
		ddx_.at(var_index) = x_init_[2];

		y_.at(var_index) = y_init_[0];
		dy_.at(var_index) = y_init_[1];
		ddy_.at(var_index) = y_init_[2];

		var_index++;

		double bcost = 0.0, Mcost = double(mcost);

		for (int k = 0; k < segment_num; ++k)
		{

			double c[12];
			for (int i = 0; i < n_poly; i++)
			{
				c[i] = osqp_work->solution->x[sub_shift + i];
				c[i+n_poly] = osqp_work->solution->x[sub_shift + i + segment_num*n_poly];
				//cout << c[coe_shift+i] << " ";
			}
			//cout << endl;
			bcost += CalculateCost(c, new_corridor[k].t, sub_shift);
			sub_shift += n_poly;

			int linter = new_corridor[k].t/delta_s_ ; //new_corridor[k].end_t - new_corridor[k].beg_t; // why not use new_corridor[k].t directly here, we used that above in place of h_k ?

			for (int l = 1; l <= linter; l++)
			{	
				// cout<<"l:\t"<<l<<"\tlinter:\t"<<linter<<"\tk:\t"<<k<<"\tvar_index:"<<var_index<<"\n";
				x_.at(var_index) = 0.0;
				dx_.at(var_index) = 0.0;
				ddx_.at(var_index) = 0.0;

				y_.at(var_index) = 0.0;
				dy_.at(var_index) = 0.0;
				ddy_.at(var_index) = 0.0;

				for (int i = 0; i < n_poly; i++)
				{
					x_.at(var_index) += c[i] * b_coe[i][0] * pow(double(l) / linter, i) * pow(1 - double(l) / linter, traj_order - i);
					y_.at(var_index) += c[i+n_poly] * b_coe[i][0] * pow(double(l) / linter, i) * pow(1 - double(l) / linter, traj_order - i);
				}

				x_.at(var_index) = x_.at(var_index) * new_corridor[k].t;
				y_.at(var_index) = y_.at(var_index) * new_corridor[k].t;
				// if(x_.at(var_index) > 50)	cout<<"\nx is \t"<<x_.at(var_index)<<" at new_corridor \t"<<k<<"\n";

				for (int i = 0; i < n_poly - 1; i++)
				{
					dx_.at(var_index) += traj_order * (c[i + 1] - c[i]) * b_coe[i][1] * pow(double(l) / linter, i) * pow(1 - double(l) / linter, traj_order - 1 - i);
					dy_.at(var_index) += traj_order * (c[i + 1+n_poly] - c[i+n_poly]) * b_coe[i][1] * pow(double(l) / linter, i) * pow(1 - double(l) / linter, traj_order - 1 - i);
				}

				//dx_.at(var_index) = dx_.at(var_index) * new_corridor[k].t;

				for (int i = 0; i < n_poly - 2; i++)
				{
					ddx_.at(var_index) += traj_order * (traj_order - 1) * (c[i + 2] - 2.0 * c[i + 1] + c[i]) * b_coe[i][2] * pow(double(l) / linter, i) * pow(1 - double(l) / linter, traj_order - 2 - i);
					ddy_.at(var_index) += traj_order * (traj_order - 1) * (c[i + 2+n_poly] - 2.0 * c[i + 1+n_poly] + c[i+n_poly]) * b_coe[i][2] * pow(double(l) / linter, i) * pow(1 - double(l) / linter, traj_order - 2 - i);
				}

				ddx_.at(var_index) = ddx_.at(var_index) / new_corridor[k].t;
				ddy_.at(var_index) = ddy_.at(var_index) / new_corridor[k].t;

				var_index++;
			}
			// cout<<"\nx size is \t"<<x_.size()<<"and y size is\t"<<y_.size()<<"\n";

			Mcost += weight_x_ref_ / (3.0 * x_skew_[k]) * (pow(x_skew_[k] * new_corridor[k].end_t * 0.1 + x_bias_[k], 3) - pow(x_skew_[k] * new_corridor[k].beg_t * 0.1 + x_bias_[k], 3));
		}

		double c_last = osqp_work->solution->x[segment_num * n_poly - 1];

		bcost += weight_end_ref_ * pow((new_corridor[segment_num - 1].t * c_last), 2); // bcost is same as mcost, which is the minimum cost of our objective function

		Mcost += weight_dx_ref_ * dx_ref_ * dx_ref_ * num_of_knots_ * delta_s_ + weight_end_ref_ * x_ref_[num_of_knots_ - 1] * x_ref_[num_of_knots_ - 1]; // what is M_cost and what does it denote semantically?

		cout << " bcost: " << bcost << std::endl;
		cout << " Mcost: " << Mcost << std::endl;

		CHECK_EQ(var_index, num_of_points_);

		// Cleanup
		osqp_cleanup(osqp_work);
		FreeData(data);
		c_free(settings);
		return true;
	}

	// return cost or value of our objective funciton
	double PiecewiseJerkSpeedProblem::CalculateCost(const double c[], double scale_t, int q_beg)
	{
		double bcost = 0.0;

		MatrixXd c_m, Q_m, temp_cost;

		temp_cost.resize(1, 1);
		c_m.resize(n_poly, 1);
		Q_m.resize(n_poly, n_poly);

		for (int i = 0; i < n_poly; i++)
		{
			c_m(i, 0) = c[i];
		}

		for (int i = 0; i < n_poly; i++)
		{
			for (int j = 0; j < n_poly; j++)
			{
				Q_m(i, j) = pow(scale_t, 3) * MQM_x[0](i, j) + scale_t * MQM_x[1](i, j) + MQM_x[2](i, j) / scale_t + MQM_x[3](i, j) / pow(scale_t, 3);
			}
			bcost += c_m(i, 0) * q[q_beg + i];
		}

		temp_cost = c_m.transpose() * Q_m * c_m;
		bcost = bcost + temp_cost(0, 0);
		return bcost;
	}

	OSQPSettings *PiecewiseJerkSpeedProblem::SolverDefaultSettings()
	{
		// Define Solver default settings
		OSQPSettings *settings =
			reinterpret_cast<OSQPSettings *>(c_malloc(sizeof(OSQPSettings)));
		osqp_set_default_settings(settings);
		settings->eps_abs = 0.00025;
		settings->eps_rel = 0.00025;
		settings->eps_prim_inf = 0.000025;
		settings->eps_dual_inf = 0.000025;
		settings->polish = true;
		//settings->verbose = FLAGS_enable_osqp_debug;
		settings->verbose = false;
		settings->scaled_termination = true;

		return settings;
	}

}