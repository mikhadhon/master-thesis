#include <iostream>

#include <CGAL/Epick_d.h>
#include <CGAL/K_neighbor_search.h>
#include <CGAL/Search_traits_d.h>
#include <CGAL/Search_traits_adapter.h>
#include "igl/readOBJ.h"
#include "igl/random_points_on_mesh.h"

typedef CGAL::Epick_d<CGAL::Dimension_tag<4>> K;
typedef K::Point_d Point;
typedef K::Vector_d Vector;
typedef K::Squared_distance_d squared_distance;
typedef K::FT FT;

typedef CGAL::Search_traits_d<K> TreeTraits;
typedef CGAL::K_neighbor_search<TreeTraits> search;
typedef search::iterator Search_iterator;
typedef search::Tree Tree;

void read_obj(const std::string &file_path, Eigen::MatrixXd &V, Eigen::MatrixXi &F) {
    igl::readOBJ(file_path, V, F);
}

FT hausdorff(Tree &tree, Eigen::MatrixXd test_points) {
    FT max_dist = 0;
    for (int i = 0; i < test_points.rows(); i++) {
        Point p(test_points.row(i)[0], test_points.row(i)[1], test_points.row(i)[2], test_points.row(i)[3]);
        search nearest(tree, p);

        if (FT current_dist = squared_distance()(p, nearest.begin()->first); current_dist > max_dist) {
            max_dist = squared_distance()(p, nearest.begin()->first);
        }
    }
    return max_dist;
}

FT hausdorff_two_sided(Tree& tree_0, Tree& tree_1, Eigen::MatrixXd V_0, Eigen::MatrixXd V_1) {
    FT max_dist_0 = hausdorff(tree_0, V_1);
    FT max_dist_1 = hausdorff(tree_1, V_0);

    return std::max(max_dist_0, max_dist_1);
}

int main(int argc, char *argv[]) {
    Eigen::MatrixXd V;
    Eigen::MatrixXi F;
    Eigen::MatrixXd V_flow;
    Eigen::MatrixXi F_flow;
    Eigen::MatrixXd B_flow;
    Eigen::MatrixXi FI_flow;
    Eigen::MatrixXd X_flow;
    Eigen::MatrixXd V_con;
    Eigen::MatrixXi F_con;
    Eigen::MatrixXd B_con;
    Eigen::MatrixXi FI_con;
    Eigen::MatrixXd X_con;

    read_obj("../hopf_torus_res_1024.obj", V, F);
    read_obj("../../hopf_torus_res_80.obj_6000_20260323_140829.obj", V_flow, F_flow);
    read_obj("../continuation_hopf_res_80.obj", V_con, F_con);

    std::cout << V.rows() << std::endl;
    std::cout << V_flow.rows() << std::endl;
    std::cout << V_con.rows() << std::endl;

    igl::random_points_on_mesh(10000, V_flow, F_flow, B_flow, FI_flow, X_flow);
    igl::random_points_on_mesh(10000, V_con, F_con, B_con, FI_con, X_con);

    std::vector<Point> points_truth;
    for (int i = 0; i < V.rows(); i++) {
        points_truth.push_back(Point(V.row(i)[0], V.row(i)[1], V.row(i)[2], V.row(i)[3]));
    }
    std::vector<Point> points_flow;
    for (int i = 0; i < V_flow.rows(); i++) {
        points_flow.push_back(Point(V_flow.row(i)[0], V_flow.row(i)[1], V_flow.row(i)[2], V_flow.row(i)[3]));
    }
    std::vector<Point> points_con;
    for (int i = 0; i < V_con.rows(); i++) {
        points_con.push_back(Point(V_con.row(i)[0], V_con.row(i)[1], V_con.row(i)[2], V_con.row(i)[3]));
    }

    Tree tree_truth(points_truth.begin(), points_truth.end());
    Tree tree_flow(points_flow.begin(), points_flow.end());
    Tree tree_con(points_con.begin(), points_con.end());

    FT truth_flow = hausdorff_two_sided(tree_truth, tree_flow, V, V_flow);
    FT truth_con = hausdorff_two_sided(tree_truth, tree_con, V, V_con);


    std::cout << truth_con << std::endl;
    std::cout << truth_flow << std::endl;

    return 0;
}