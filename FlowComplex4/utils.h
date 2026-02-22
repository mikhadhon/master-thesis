#ifndef UTILS_H
#define UTILS_H

#define _USE_MATH_DEFINES

#include "Delaunay.h"

Eigen::MatrixXd cliffordgen(int nsamples);

Eigen::MatrixXd stereo_projection(Eigen::MatrixXd object);

std::vector<Eigen::VectorXd> stereo_projection(std::vector<Eigen::VectorXd> object);

std::string get_timestamp();

void write_to_obj(Eigen::MatrixXd V, Eigen::MatrixXi F, std::string identifier);

void read_ply(const std::string &file_path, Eigen::MatrixXd &V, Eigen::MatrixXi &F);

void read_obj(const std::string &file_path, Eigen::MatrixXd &V, Eigen::MatrixXi &F);


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

struct Vertex_to_point {
    Point operator()(Delaunay::Vertex_handle vh) const { return vh->point(); }
};

#endif
