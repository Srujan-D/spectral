#include "../include/btrapz/logging.h"
#include "../include/btrapz/piecewise_jerk_speed_problem.h"
#include "FrenetOptimalTrajectory/FrenetOptimalTrajectory.h"
#include "FrenetOptimalTrajectory/FrenetPath.h"
#include "FrenetOptimalTrajectory/py_cpp_struct.h"

// #include "logging.h"
// #include "piecewise_jerk_speed_problem.h"

// #include "piecewise_jerk_path_problem.h"

#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>

#include <ros/ros.h>

// #include <dynamic_reconfigure/server.h>
// #include <dynamic_reconfigure/DoubleParameter.h>
// #include <dynamic_reconfigure/Reconfigure.h>
// #include <dynamic_reconfigure/Config.h>

#include <sensor_msgs/LaserScan.h>
#include <geometry_msgs/Twist.h>
#include <nav_msgs/Path.h>
#include <geometry_msgs/PoseStamped.h>
#include <nav_msgs/Odometry.h>

#include <tf/tf.h>

std::vector<std::pair<double, double>> x_bounds;
std::vector<std::pair<double, double>> dx_bounds;
std::vector<double> x_ref;
std::vector<double> x_kappa;
std::vector<double> penalty_dx;
nav_msgs::Path new_path;

std::vector<std::pair<double, double>> dx_y;

std::vector<std::tuple<double, double, double>> cart_path = {{1, 1, 0}, {2, 1, 0}, {2, 2, 0}};
std::vector<std::vector<double>> frenet_path;
std::vector<std::pair<double, double>> cart_waypoints = {{1, 1}, {2, 1}, {2, 2}};
std::vector<double> wx, wy;

nav_msgs::Path global_path;
geometry_msgs::PoseStamped g_waypoints;
// std::tuple<double, double,double> cur_odom;
// std::vector <double> cur_odom;
std::vector<std::pair<double, double>> btrapz_path;

using namespace std;
using namespace navigation;




// void getOdomCallback(const nav_msgs::Odometry::ConstPtr &msg){
// 	tf::Quaternion q;
// 	tf::Matrix3x3 m;
// 	double r, p, theta;
// 	q = tf::Quaternion(msg->pose.pose.orientation.x, msg->pose.pose.orientation.y, msg->pose.pose.orientation.z, msg->pose.pose.orientation.w);
// 	m = tf::Matrix3x3(q);
// 	m.getRPY(r, p, theta);
// 	cur_odom = {msg->pose.pose.position.x, msg->pose.pose.position.y, theta} ;
// }

void getPathCallback(const nav_msgs::Path::ConstPtr &msg)
{
	global_path = *msg;
	auto n_waypoints = global_path.poses.size();
	std::tuple<double, double, double> x_y_theta;
	tf::Quaternion q;
	tf::Matrix3x3 m;
	double r, p, theta;
	for (int i = 0; i < n_waypoints; i++)
	{
		q = tf::Quaternion(global_path.poses[i].pose.orientation.x, global_path.poses[i].pose.orientation.y, global_path.poses[i].pose.orientation.z, global_path.poses[i].pose.orientation.w);
		m = tf::Matrix3x3(q);
		m.getRPY(r, p, theta);
		cart_path.push_back(make_tuple(global_path.poses[i].pose.position.x, global_path.poses[i].pose.position.y, theta));
		// cart_waypoints.push_back(make_pair(global_path.poses[i].pose.position.x, global_path.poses[i].pose.position.y));
		wx.push_back(global_path.poses[i].pose.position.x);
		wy.push_back(global_path.poses[i].pose.position.y);
	}
}

inline constexpr double pi() { return M_PI; }
inline double deg2rad(double x) { return x * pi() / 180; }
inline double rad2deg(double x) { return x * 180 / pi(); }

double distance_(double x1, double y1, double x2, double y2)
{
	return sqrt((x2 - x1) * (x2 - x1) + (y2 - y1) * (y2 - y1));
}

// int ClosestWaypoint(double x, double y, std::vector<std::pair<double, double>> path)
// {
// 	cout<<"\n in ClosestWaypoint \n";
// 	double closestLen = 100000; // large number
// 	int closestWaypoint = 0;

// 	for (int i = 0; i < path.size(); i++)
// 	{
// 		cout<<"finding closest\n";
// 		double map_x = path[i].first;
// 		double map_y = path[i].second;
// 		double dist = distance_(x, y, map_x, map_y);
// 		if (dist < closestLen)
// 		{
// 			closestLen = dist;
// 			closestWaypoint = i;
// 			cout<<"curr dist - "<<closestLen<<"\t curr wayp - "<<closestWaypoint<<"\n";
// 		}
// 	}
// 	cout<<"found closest! - "<< closestWaypoint<<"\n";
// 	return closestWaypoint;
// }

// int NextWaypoint(double x, double y, double theta, std::vector<std::pair<double, double>> path)
// {
// 	cout<<"\n in NextWaypoint\n";
// 	int closestWaypoint = ClosestWaypoint(x, y, path);

// 	double map_x = path[closestWaypoint].first;
// 	double map_y = path[closestWaypoint].second;

// 	double heading = atan2((map_y - y), (map_x - x));

// 	double angle = abs(theta - heading);

// 	if (angle > pi() / 4)
// 	{
// 		closestWaypoint++;
// 	}
// 	cout<<"found next waypoint! - "<<closestWaypoint<<"\n";
// 	return closestWaypoint;
// }

// // Transform from Cartesian x,y coordinates to Frenet s,d coordinates
// vector<double> getFrenet(double x, double y, double theta, const std::vector<std::pair<double, double>> path)
// {
// 	cout<<"\n in get frenet\n";
// 	int next_wp = NextWaypoint(x, y, theta, path);
// 	int prev_wp;
// 	prev_wp = next_wp - 1;
// 	if (next_wp == 0)
// 	{
// 		prev_wp = path.size() - 1;
// 	}
// 	double n_x = path[next_wp].first - path[prev_wp].first;
// 	double n_y = path[next_wp].second - path[prev_wp].second;
// 	double x_x = x - path[prev_wp].first;
// 	double x_y = y - path[prev_wp].second;
// 	// find the projection of x onto n
// 	double proj_norm = (x_x * n_x + x_y * n_y) / (n_x * n_x + n_y * n_y);
// 	double proj_x = proj_norm * n_x;
// 	double proj_y = proj_norm * n_y;
// 	double frenet_d = distance_(x_x, x_y, proj_x, proj_y);
// 	// see if d value is positive or negative by comparing it to a center point
// 	double center_x = 1000 - path[prev_wp].first;
// 	double center_y = 2000 - path[prev_wp].second;
// 	double centerToPos = distance_(center_x, center_y, x_x, x_y);
// 	double centerToRef = distance_(center_x, center_y, proj_x, proj_y);
// 	if (centerToPos <= centerToRef)
// 	{
// 		frenet_d *= -1;
// 	}
// 	// calculate s value
// 	double frenet_s = 0;
// 	for (int i = 0; i < prev_wp; ++i)
// 	{
// 		frenet_s += distance_(path[i].first, path[i].second, path[i + 1].first, path[i + 1].second);
// 	}
// 	frenet_s += distance_(0, 0, proj_x, proj_y);
// 	return {frenet_s, frenet_d};
// }

// Transform from Frenet s,d coordinates to Cartesian x,y
vector<double> getXY(double s, double d, vector<double> maps_s, vector<double> maps_x, vector<double> maps_y)
{
	int prev_wp = -1;

	while (s > maps_s[prev_wp + 1] && (prev_wp < (int)(maps_s.size() - 1)))
	{
		prev_wp++;
	}

	int wp2 = (prev_wp + 1) % maps_x.size();
	
	double heading = atan2((maps_y[wp2] - maps_y[prev_wp]), (maps_x[wp2] - maps_x[prev_wp]));
	// the x,y,s along the segment
	double seg_s = (s - maps_s[prev_wp]);

	double seg_x = maps_x[prev_wp] + seg_s * cos(heading);
	double seg_y = maps_y[prev_wp] + seg_s * sin(heading);

	double perp_heading = heading - pi() / 2;

	double x = seg_x + d * cos(perp_heading);
	double y = seg_y + d * sin(perp_heading);

	return {x, y};
}

// void convertPath2Frenet(const std::vector<std::tuple<double, double, double>> cart_path)
// {
// 	// TODO
// 	std::vector<double> s_d;
// 	for(auto& point : cart_path){
// 		s_d = getFrenet(get<0>(point), get<1>(point), get<2>(point), cart_waypoints);
// 		frenet_path[0].push_back(s_d[0]);
// 		frenet_path[1].push_back(s_d[1]);
// 		cout<<"converting\n";
// 	}
// }

// void control(ros::Publisher pid, const std::vector<std::pair<double, double>> btrapz_path){

// 	for(int i=0; i<btrapz_path.size(); i++){
// 		while(distance_(cur_odom[0], cur_odom[1], btrapz_path[i].first, btrapz_path[i].second) < 0.1){

// 		}
// 	}

// 	return ;
// }

geometry_msgs::PoseStamped temp;

geometry_msgs::PoseStamped makePose(double x, double y)
{
	temp.pose.position.x = x;
	temp.pose.position.y = y;
	temp.header.frame_id = "map";
	return temp;
}
void convertXY2Path(const string path_id, const std::vector<std::pair<double, double>> btrapz_path)
{
	new_path.header.frame_id = path_id;
	for (int i = 0; i < btrapz_path.size(); i++)
	{
		new_path.poses.push_back(makePose(btrapz_path[i].first, btrapz_path[i].second));
	}
}



int main(int argc, char **argv)
{
	ros::init(argc, argv, "btrapz_node");
	ros::NodeHandle nh;


	string path_id = "odom";
	// cout<<"---"<<path_id<<"--\n";

	if(nh.hasParam("path")){
		cout<<"check\n";
		nh.getParam("path", path_id);
	}
	// else{
	// 	cout<<"set\n";
	// 	nh.param<std::string>("path", path_id, "odom");
	// }
	// cout<<"---"<<path_id<<"--\n";

	ros::Subscriber get_gp = nh.subscribe("/plan", 1000, getPathCallback);
	// ros::Subscriber get_odom = nh.subscribe("/odom", 1000, getOdomCallback);

	ros::Publisher pure_p = nh.advertise<nav_msgs::Path>("/path_segment", 1000);

	// convertPath2Frenet(cart_path);

	double o_llx[1] = {92.89};
	double o_lly[1] = {191.75};
	double o_urx[1] = {92.89};
	double o_ury[1] = {191.75};

	// wy = {0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 5.5, 6.0, 6.5, 7.0,
	// 	  7.5, 8.0, 8.5, 9.0, 9.5, 10.0, 10.5, 11.0, 11.5, 12.0, 12.5, 13.0, 13.5, 14.0, 14.5};
	// wx = {0.0, 2.12, 3.0, 3.67, 4.24, 4.74, 5.2, 5.61, 6.0, 6.36, 6.71, 7.04, 7.35,
	// 	  7.65, 7.94, 8.22, 8.49, 8.75, 9.0, 9.25, 9.49, 9.72, 9.95, 10.17, 10.39,
	// 	  10.61, 10.82, 11.02, 11.22, 11.42};
	// set up experiment
	FrenetInitialConditions fot_ic = {
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        10,
        wx,
        wy,
        1,
        o_llx,
        o_lly,
        o_urx,
        o_ury,
        1
    };
    FrenetHyperparameters fot_hp = {
        25.0,
        5.0,
        10.0,
        1.0,
        1.0,
        0.1,
        0.1,
        6.0,
        1.0,
        0.1,
        2.0,
        3.0,
        1.0,
        0.1,
        0.1,
        0.1,
        0.1,
        1.0,
        1.0,
        1.0
    };
	// run experiment
	FrenetOptimalTrajectory fot = FrenetOptimalTrajectory(&fot_ic, &fot_hp);
	FrenetPath *best_frenet_path = fot.getBestPath();
	cout << "\n best frenet path achieved ------ " << best_frenet_path->cf << "\n";

	ofstream ofs2("/home/srujan/RISS/code/riss/src/btrapz/src/bezier_path_ref.txt");
	if (!ofs2.is_open())
	{
		cout << "Not Found" << endl;
		return 1;
	}

	for (int i = 0; i < best_frenet_path->s.size(); ++i)
	{
		ofs2 << std::fixed << std::setprecision(3) << best_frenet_path->s[i] << " " << best_frenet_path->d[i] << endl;
	}
	ofs2.close();

	// double mscale = 10.0;
	// int num_of_knots = 71;
	// double delta_t = 0.1, acc_weight = 10, jerk_weight = 5;
	// std::array<double, 3> init_s{0, 1.0, 0};

	ifstream ifs("/home/srujan/RISS/code/riss/src/btrapz/src/jer_speed_info.txt");
	// ifstream ifs("piecewise/jerk_speed.txt");
	if (!ifs.is_open())
	{
		cout << "Not Found" << endl;
		return 1;
	}
	double mscale = 10.0;
	int num_of_knots;
	double delta_t, acc_weight, jerk_weight;
	std::array<double, 3> init_s;
	ifs >> num_of_knots >> delta_t;
	ifs >> init_s[0] >> init_s[1] >> init_s[2];
	ifs >> acc_weight >> jerk_weight;

	init_s[1] = init_s[1] * mscale;
	num_of_knots = best_frenet_path->s.size();
	AWARN << "knots:" << num_of_knots << " delta t:" << delta_t;
	AWARN << "s:" << init_s[0] << " v:" << init_s[1] << " a:" << init_s[2];
	AWARN << "w_acc: " << acc_weight << " w_jerk:" << jerk_weight;

	double val[5];
	for (int i = 0; i < num_of_knots; i++)
	{
		ifs >> val[0] >> val[1];
		x_bounds.emplace_back(mscale * val[0], mscale * val[1]);
	}
	for (int i = 0; i < num_of_knots; i++)
	{
		ifs >> val[0] >> val[1];
		dx_bounds.emplace_back(mscale * val[0], mscale * val[1]);
	}

	cout << "\n starting solving btrapz \n --------------------------- \n";

	PiecewiseJerkSpeedProblem piecewise_jerk_problem(num_of_knots, delta_t, init_s);
	cout << "\n setting vars \n --------------------------- \n";
	piecewise_jerk_problem.set_scale_factor({1.0, 1.0, 1.0});
	cout << "\n set_weight_ddx \n --------------------------- \n";
	piecewise_jerk_problem.set_weight_ddx(acc_weight);
	cout << "\n set_weight_dddx \n --------------------------- \n";
	piecewise_jerk_problem.set_weight_dddx(jerk_weight);
	cout << "\n set_ddx_bounds \n --------------------------- \n";
	piecewise_jerk_problem.set_ddx_bounds(-2.5, 1.5);
	cout << "\n set_dddx_bound \n --------------------------- \n";
	piecewise_jerk_problem.set_dddx_bound(-40, 40);
	cout << "\n set_x_bounds \n --------------------------- \n";
	piecewise_jerk_problem.set_x_bounds(x_bounds);
	cout << "\n set_dx_bounds \n --------------------------- \n";
	piecewise_jerk_problem.set_dx_bounds(dx_bounds);
	cout << "\n set_x_ref \n " << best_frenet_path->s[0] << "--------------------------- \n";

	// piecewise_jerk_problem.set_penalty_dx(penalty_dx);
	piecewise_jerk_problem.set_x_ref(0.10, best_frenet_path->s);
	cout << "\n set_dx_ref \n --------------------------- \n";
	piecewise_jerk_problem.set_dx_ref(0.10, 1.2 * mscale);
	cout << "\n starting time \n --------------------------- \n";
	auto start_time = std::chrono::system_clock::now();
	cout << "\n optimizing **** \n --------------------------- \n";
	bool success = piecewise_jerk_problem.Optimize(4000);
	auto end_time = std::chrono::system_clock::now();
	std::chrono::duration<double> diff = end_time - start_time;
	// ADEBUG << "Speed Optimizer used time: " << diff.count() * 1000 << " ms.";
	// if (!success)
	// {
	// 	std::string msg("Piecewise jerk speed optimizer failed!");
	// 	AERROR << msg;
	// }
	// Extract output
	const std::vector<double> &s = piecewise_jerk_problem.opt_x();
	const std::vector<double> &ds = piecewise_jerk_problem.opt_dx();
	const std::vector<double> &dds = piecewise_jerk_problem.opt_ddx();
	bool is_verbose = true;
	double a_cost = 0.0;
	double mmax_a = 0.0;
	for (int i = 0; i < num_of_knots && is_verbose; ++i)
	{
		double ddds = (i == 0) ? (dds[1] - dds[0]) / delta_t : (dds[i] - dds[i - 1]) / delta_t;

		a_cost += dds[i] * dds[i] * delta_t;
		mmax_a = max(mmax_a, abs(dds[i]));

		cout << std::fixed << std::setprecision(3) << "For t[" << i * delta_t << "], opt = " << s[i] << ", " << ds[i] << ", " << dds[i] << ", " << ddds << endl;
	}
	a_cost = acc_weight * sqrt(a_cost / (num_of_knots * delta_t));
	cout << "mmax_a " << mmax_a << endl;
	cout << "a_cost " << a_cost << endl;

	ADEBUG << "Speed Optimizer used time: " << diff.count() * 1000 << " ms.";

	// wx = {0.0, 1.5, 2.12, 2.6, 3.0, 3.35, 3.67, 3.97, 4.24, 4.5, 4.74, 4.97, 5.2, 5.41, 5.61, 5.81, 6.0, 6.18, 6.36, 6.54, 6.71, 6.87, 7.04, 7.19, 7.35, 7.5, 7.65, 7.79, 7.94, 8.08, 8.22, 8.35, 8.49, 8.62, 8.75, 8.87, 9.0, 9.12, 9.25, 9.37, 9.49, 9.6, 9.72, 9.84, 9.95, 10.06, 10.17, 10.28, 10.39, 10.5, 10.61, 10.71, 10.82, 10.92, 11.02, 11.12, 11.22, 11.32, 11.42, 11.52};
	// wy = {0.0, 0.25, 0.5, 0.75, 1.0, 1.25, 1.5, 1.75, 2.0, 2.25, 2.5, 2.75, 3.0, 3.25, 3.5, 3.75, 4.0, 4.25, 4.5, 4.75, 5.0, 5.25, 5.5, 5.75, 6.0, 6.25, 6.5, 6.75, 7.0, 7.25, 7.5, 7.75, 8.0, 8.25, 8.5, 8.75, 9.0, 9.25, 9.5, 9.75, 10.0, 10.25, 10.5, 10.75, 11.0, 11.25, 11.5, 11.75, 12.0, 12.25, 12.5, 12.75, 13.0, 13.25, 13.5, 13.75, 14.0, 14.25, 14.5, 14.75};
	std::vector<double> x_y_temp;
	for (int i = 0; i < s.size(); i++)
	{
		x_y_temp = getXY(s[i], best_frenet_path->d[i], s, wx, wy);
		btrapz_path.push_back(make_pair(x_y_temp[0], x_y_temp[1]));
	}

	ofstream ofs3("/home/srujan/RISS/code/riss/src/btrapz/src/btrapz_path_computed.txt");
	if (!ofs3.is_open())
	{
		cout << "Not Found" << endl;
		return 1;
	}

	for (int i = 0; i < s.size() && is_verbose; ++i)
	{
		ofs3 << std::fixed << std::setprecision(3) << i*0.1 << " " << s[i] << endl;
	}

	ofstream ofs("/home/srujan/RISS/code/riss/src/btrapz/src/x_y_path.txt");
	if (!ofs.is_open())
	{
		cout << "Not Found" << endl;
		return 1;
	}

	for (int i = 0; i < num_of_knots && is_verbose; ++i)
	{
		ofs << std::fixed << std::setprecision(3) << btrapz_path[i].first << " " << btrapz_path[i].second << endl;
	}
	
	// control(pid, btrapz_path);

	convertXY2Path(path_id, btrapz_path);
	pure_p.publish(new_path);

	ros::Rate rate(100);
	while (ros::ok())
	{
		pure_p.publish(new_path);
		ros::spinOnce();
		rate.sleep();
	}

	return 0;
}
