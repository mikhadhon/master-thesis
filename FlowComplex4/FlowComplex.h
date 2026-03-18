#ifndef FLOWCOMPLEX_H
#define FLOWCOMPLEX_H
#include "Delaunay.h"

struct Saddle_2 {
    Face df;
    Voronoi_face dual;

    bool operator==(const Saddle_2 & other) const {
        return df == other.df;
    }

    Saddle_2(Face df, Voronoi_face dual) : df(df), dual(dual) {}
};

struct stable_manifold_2 {
    std::vector<Edge> gabriel_edges;
    Saddle_2 saddle_2;
    std::vector<std::array<size_t, 3>> faces;

    bool operator==(const stable_manifold_2 & other) const {
        return saddle_2 == other.saddle_2;
    }

    stable_manifold_2(std::vector<Edge> gabriel_edges, Saddle_2 saddle_2, std::vector<std::array<size_t, 3>> faces) : gabriel_edges(gabriel_edges), saddle_2(saddle_2), faces(faces) {}
};

struct valid_pair {
    stable_manifold_2 sm;
    Voronoi_vertex v;
    FT distance;

    bool operator==(const valid_pair & other) const {
        return sm == other.sm && v == other.v;
    }

    bool operator<(const valid_pair & other) const {
        return distance < other.distance;
    }

    valid_pair(stable_manifold_2 sm, Voronoi_vertex v, FT distance) : sm(sm), v(v), distance(distance) {}
};

struct helper_data {
    int vertex_count;
    std::map<Delaunay::Vertex_handle, size_t> vertex_to_index;
    std::map<std::array<double, 4>, size_t> fc_vertex_to_index;
    std::vector<std::array<size_t, 2>> gabriel_edges;
    std::vector<stable_manifold_2> stable_manifolds;
    std::vector<valid_pair> valid_pairs;
    std::map<Edge, int> gabriel_topology;
};

void flow_complex(Delaunay &delaunay, std::vector<Eigen::VectorXd> &vertices, std::vector<std::array<size_t, 3>> &faces, helper_data &data);

void reduce_flow_complex(Delaunay &dt, helper_data &data);
#endif
