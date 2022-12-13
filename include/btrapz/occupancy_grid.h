#pragma once

#include <iostream>
#include <vector>
#include <set>
#include <algorithm>
#include <array>
#include <cmath>

using namespace std;

#include <CGAL/Simple_cartesian.h>
#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_traits.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/boost/graph/graph_traits_Polyhedron_3.h>
#include <CGAL/AABB_face_graph_triangle_primitive.h>
#include <CGAL/algorithm.h>
#include <CGAL/Side_of_triangle_mesh.h>
#include <CGAL/convex_hull_3.h>

typedef CGAL::Simple_cartesian<double> K;
typedef K::Point_3 Point;
typedef CGAL::Polyhedron_3<K> Polyhedron;
typedef CGAL::AABB_face_graph_triangle_primitive<Polyhedron> Primitive;
typedef CGAL::AABB_traits<K, Primitive> Traits;
typedef CGAL::AABB_tree<Traits> Tree;
typedef CGAL::Side_of_triangle_mesh<Polyhedron, K> Point_inside;



#define EPSILON  0.0000001
#define MODULUS(p) (sqrt(p.s*p.s + p.l*p.l + p.t*p.t))
#define TWOPI 6.283185307179586476925287
#define RTOD 57.2957795

class Grid{
    public:

		// Grid() {}
		~Grid() {}

        std::array<int, 3> map_size = {{300, 100, 81}};               // s, d, t
		
		std::array<double, 3> map_resolution = {{0.25, 0.2, 0.1}};  // m, m, s
		std::array<std::string, 3> axis_name = {{"s", "d", "t"}};

		typedef struct{
			double s;
			double l;
			double t;
		} SLT;

		typedef	std::vector<Point> corridor;
		corridor c;

		// int n = map_size[0]*map_size[1]*map_size[2];
		std::array<int, 300*100*81> flat_map{};	// (X,Y,Z) -> (X + Y * DX + Z * DY * DX), where DX and DY are the dimensions on X and Y, respectively
		std::array<int, 300*100*81> is_counted{};	// (X,Y,Z) -> (X + Y * DX + Z * DY * DX), where DX and DY are the dimensions on X and Y, respectively

		double l_car = 2.0;
        double l_ego = 2.0;

        double s_car = 5.0;
        double s_ego = 5.0;

        double l_safe = l_car/2 + l_ego/2;
        double s_safe = s_car/2 + s_ego/2;

		double time = 2.0;
		SLT forward_state;
		Point coord;

		double s_back_len = 0.0;
		double kMaxLongitudinalVel = 50.0;
		double kMinLongitudinalVel = 0.0;
		double kMaxLongitudinalAcc = 3.0;
		double kMaxLongitudinalDecel = -8.0;  // Avg. driver max
		double kMaxLateralVel = 3.0;
		double kMaxLateralAcc = 2.5;

		int kMaxNumOfGridAlongTime = 2;

		// std::array<int, 6> inflate_steps = {{20, 5, 10, 10, 1, 1}};

		std::array<int, 3> max_corridor_size = {30, 10, 10};
		std::array<int, 3> cur_corridor_size = {30, 10, 10};

		std::vector<Polyhedron> obstacles;
		Point temp;

		std::vector<corridor> convex_regions;

		void ConstructCorridors();

        void UpdateMapWithObstacle(const SLT obs, double vel_s=0.0, double vel_l=0.0);

		// double CalcAngleSum(SLT q, SLT *p, int n);

		SLT p[8] ;
		// void CalcObstacle(SLT centre, double vel_s, double vel_l);
		std::array<Point, 8> points;

		bool IsCubeFree(const SLT coord);

		Polyhedron poly;
		void CalcObstacle(SLT centre, double vel_s, double vel_l);
		bool pointInside(Point &query);
};