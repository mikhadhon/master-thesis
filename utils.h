#ifndef UTILS_H
#define UTILS_H

#define _USE_MATH_DEFINES

#include "Delaunay.h"

std::string get_timestamp();

void write_to_obj(Eigen::MatrixXd V, Eigen::MatrixXi F);

void read_ply(const std::string &file_path, Eigen::Matrix<double, Eigen::Dynamic, 3> &V, Eigen::Matrix<double, Eigen::Dynamic, 3> &F);

void gen_sphere_sample(int count, double radius, std::vector<Delaunay::Point> &points);

void gen_rectangle(int n, std::vector<Delaunay::Point> &points);

Point make_point(Eigen::VectorXd &eigen_point);

Eigen::VectorXd make_point_eigen(Point point);

std::vector<Delaunay::Point> get_points_from_handles(const std::vector<Delaunay::Vertex_handle> &handles);

void map_vertices_to_vector(
    Delaunay &delaunay,
    std::vector<Eigen::VectorXd> &vertices,
    std::map<Delaunay::Vertex_handle, size_t> &vertex_to_index
);

void write_faces_to_vector(
    Delaunay &delaunay,
    std::vector<std::array<size_t, 3> > &faces,
    std::map<Delaunay::Vertex_handle, size_t> &vertex_to_index
);

void generate_circle(int nsamples, int normal_windings, Eigen::Matrix<double, Eigen::Dynamic, 3> &V, Eigen::Matrix<double, Eigen::Dynamic, 3> &N1, Eigen::Matrix<double, Eigen::Dynamic, 3>  &N2, Eigen::Matrix<double, Eigen::Dynamic, 3> &T);

void generate_torus(int count, double radius, double rot_radius, std::vector<Delaunay::Point> &torus_samples);

bool is_point_in_triangle(const Point& p0, const Point& p1, const Point& p2, const Point& query);

#endif
