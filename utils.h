#ifndef UTILS_H
#define UTILS_H

#define _USE_MATH_DEFINES

#include "Delaunay.h"

void gen_sphere_sample(int count, double radius, std::vector<Delaunay::Point> &points);

void gen_rectangle(int n, std::vector<Delaunay::Point> &points);

Point make_point(Eigen::Vector3d &eigen_point);

Eigen::VectorXd make_point_eigen(Point point);

void map_vertices_to_vector(
    Delaunay &delaunay,
    std::vector<std::array<double, 3> > &vertices,
    std::map<Delaunay::Vertex_handle, size_t> &vertex_to_index
);

void write_faces_to_vector(
    Delaunay &delaunay,
    std::vector<std::array<size_t, 3> > &faces,
    std::map<Delaunay::Vertex_handle, size_t> &vertex_to_index
);

void generate_circle(int nsamples, int normal_windings, Eigen::Matrix<double, Eigen::Dynamic, 3> &V, Eigen::Matrix<double, Eigen::Dynamic, 3> &N1, Eigen::Matrix<double, Eigen::Dynamic, 3>  &N2, Eigen::Matrix<double, Eigen::Dynamic, 3> &T);

void generate_torus(int count, double radius, double rot_radius, std::vector<Delaunay::Point> &torus_samples);

void triangle_circumcircle(Eigen::Vector3d &i, Eigen::Vector3d &j, Eigen::Vector3d &l, Eigen::VectorXd &center, double &radius);

void simplex_circumsphere(Delaunay::Full_cell_handle simplex, double &radius, Eigen::VectorXd &center);

bool intersect_ray_segment(
    Eigen::VectorXd &ray_origin,
    Eigen::VectorXd &ray_pass_through,
    Eigen::VectorXd &segment_start,
    Eigen::VectorXd &segment_end,
    Eigen::VectorXd &intersection
    );

bool is_inside_triangle(
    const Eigen::VectorXd& p,
    const Eigen::VectorXd& a,
    const Eigen::VectorXd& b,
    const Eigen::VectorXd& c
);

#endif
