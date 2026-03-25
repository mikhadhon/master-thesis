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
    FT max_distance;

    bool operator==(const stable_manifold_2 & other) const {
        return saddle_2 == other.saddle_2;
    }
    bool operator<(const stable_manifold_2 & other) const {
        return max_distance < other.max_distance;
    }

    stable_manifold_2(std::vector<Edge> gabriel_edges, Saddle_2 saddle_2, std::vector<std::array<size_t, 3>> faces, FT max_distance) : gabriel_edges(gabriel_edges), saddle_2(saddle_2), faces(faces), max_distance(max_distance) {}
};

struct fc_edge {
    Point vertex1;
    Point vertex2;

    fc_edge(const Point &vertex1, const Point &vertex2) : vertex1(vertex1), vertex2(vertex2) {}

    bool operator==(const fc_edge & other) const {
        return (vertex1 == other.vertex1 && vertex2 == other.vertex2) ||
            (vertex1 == other.vertex2 && vertex2 == other.vertex1);
    }
};

struct fc_edge_hash {
    std::size_t operator()(const fc_edge& e) const {
        std::size_t seed = 0;
        boost::hash_combine(seed, &e.vertex1);
        boost::hash_combine(seed, &e.vertex2);
        return seed;
    }
};

struct helper_data {
    int vertex_count;
    std::map<Delaunay::Vertex_handle, size_t> vertex_to_index;
    std::map<std::array<double, 4>, size_t> fc_vertex_to_index;
    std::vector<std::array<size_t, 2>> gabriel_edges;
    std::vector<stable_manifold_2> stable_manifolds;
    std::unordered_map<fc_edge, int, fc_edge_hash> edge_incidences;
    FT threshold;
    std::vector<Point> index_2;
};

void flow_complex(Delaunay &delaunay, std::vector<Eigen::VectorXd> &vertices, std::vector<std::array<size_t, 3>> &faces, helper_data &data);

#endif
