#ifndef FLOWCOMPLEX_H
#define FLOWCOMPLEX_H
#include "Delaunay.h"

struct Saddle_2 {
    Face df;
    Voronoi_face dual;
};

struct stable_manifold_2 {
    std::vector<Edge> gabriel_edges;
    Saddle_2 saddle_2;
};

struct helper_data {
    int vertex_count;
    std::map<Delaunay::Vertex_handle, size_t> vertex_to_index;
    std::map<std::array<double, 4>, size_t> fc_vertex_to_index;
    std::vector<std::array<size_t, 2>> gabriel_edges;
    std::vector<stable_manifold_2> sm;
};

void flow_complex(Delaunay &delaunay, std::vector<Eigen::VectorXd> &vertices, std::vector<std::array<size_t, 3>> &faces, helper_data &data);
#endif
