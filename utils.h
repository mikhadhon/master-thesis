#ifndef UTILS_H
#define UTILS_H

#define _USE_MATH_DEFINES

#include "Delaunay.h"

std::string get_timestamp();

void write_to_obj(Eigen::MatrixXd V, Eigen::MatrixXi F, std::string identifier);

void read_ply(const std::string &file_path, Eigen::MatrixXd &V, Eigen::MatrixXi &F);

void read_obj(const std::string &file_path, Eigen::MatrixXd &V, Eigen::MatrixXi &F);

void gen_sphere_sample(int count, double radius, std::vector<Delaunay::Point> &points);

Eigen::VectorXd make_point_eigen(Point point);

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

void generate_torus(int count, double radius, double rot_radius, std::vector<Delaunay::Point> &torus_samples);

void generate_trefoil(int nsamples, int normal_windings, Eigen::Matrix<double, Eigen::Dynamic, 3> &V, Eigen::Matrix<double, Eigen::Dynamic, 3>  &N1, Eigen::Matrix<double, Eigen::Dynamic, 3> &N2, Eigen::Matrix<double, Eigen::Dynamic, 3> &T);

bool intersect_segment_ray(std::pair<Eigen::VectorXd, Eigen::VectorXd> ray, std::pair<Eigen::VectorXd, Eigen::VectorXd> segment, Eigen::VectorXd &intersection);

#endif
