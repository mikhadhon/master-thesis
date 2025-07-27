#ifndef UTILS_H
#define UTILS_H

#define _USE_MATH_DEFINES

#include "Delaunay.h"

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

void extract_edges(Delaunay &delaunay, std::map<Delaunay::Vertex_handle, size_t> vertex_to_index, std::vector<std::array<size_t, 2>> &edges);

void generate_circle(int nsamples, int normal_windings, Eigen::Matrix<double, Eigen::Dynamic, 3> &V, Eigen::Matrix<double, Eigen::Dynamic, 3> &N1, Eigen::Matrix<double, Eigen::Dynamic, 3>  &N2, Eigen::Matrix<double, Eigen::Dynamic, 3> &T);

void generate_torus(int count, double radius, double rot_radius, std::vector<Delaunay::Point> &torus_samples);

#endif
