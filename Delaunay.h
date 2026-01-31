#ifndef DELAUNAY_H
#define DELAUNAY_H

#include <CGAL/Epick_d.h>
#include <CGAL/Delaunay_triangulation.h>
#include <boost/functional/hash.hpp>

typedef CGAL::Dimension_tag<3> Dimension;
typedef CGAL::Epick_d<Dimension> K;
typedef K::Point_d Point;
typedef K::FT FT;

typedef CGAL::Delaunay_triangulation<K> Delaunay;
typedef Delaunay::Geom_traits::Midpoint_d midpoint;
typedef Delaunay::Geom_traits::Construct_circumcenter_d circumcenter;
typedef Delaunay::Geom_traits::Squared_distance_d squared_distance;
typedef Delaunay::Geom_traits::Contained_in_simplex_d contained_in_simplex;

struct Edge {
    Delaunay::Vertex_handle vertex1;
    Delaunay::Vertex_handle vertex2;

    Delaunay::Vertex_handle co_vertex;

    Edge(Delaunay::Vertex_handle vertex1, Delaunay::Vertex_handle vertex2, Delaunay::Vertex_handle co_vertex) : vertex1(vertex1), vertex2(vertex2), co_vertex(co_vertex) {}

    bool operator==(const Edge & other) const {
        return (vertex1 == other.vertex1 && vertex2 == other.vertex2) ||
            (vertex1 == other.vertex2 && vertex2 == other.vertex1);
    }

    bool operator<(const Edge& other) const {
        std::pair<Delaunay::Vertex_handle, Delaunay::Vertex_handle> this_canonical = (vertex1 < vertex2) ? std::make_pair(vertex1, vertex2) : std::make_pair(vertex2, vertex1);
        std::pair<Delaunay::Vertex_handle, Delaunay::Vertex_handle> other_canonical = (other.vertex1 < other.vertex2) ? std::make_pair(other.vertex1, other.vertex2) : std::make_pair(other.vertex2, other.vertex1);
        return this_canonical < other_canonical;
    }
};

struct EdgeHash {
    std::size_t operator()(const Edge& e) const {
        std::size_t seed = 0;
        boost::hash_combine(seed, &*e.vertex1);
        boost::hash_combine(seed, &*e.vertex2);
        return seed;
    }
};

void find_gabriel_edges(Delaunay& dt, std::vector<Edge>& gabriel_edges);

struct Face {
    Delaunay::Face face;
    int index_of_covertex;

    Face(Delaunay::Face face, int index_of_covertex) : face(face), index_of_covertex(index_of_covertex) {}
};

inline bool lexicographic_less_vector(const Eigen::VectorXd &a, const Eigen::VectorXd &b) {
    int dimension = a.size();
    for (int i = 0; i < dimension; i++) {
        if (a[i] < b[i]) { return true; }
        if (a[i] > b[i]) { return false; }
    }
}

struct Voronoi_vertex {
    bool is_infinite;
    Eigen::VectorXd point;

    bool operator==(const Voronoi_vertex & other) const {
        return (is_infinite == other.is_infinite && point == other.point);
    }
};

struct Voronoi_edge {
    Voronoi_vertex vertex1;
    Voronoi_vertex vertex2;

    Delaunay::Full_cell_handle cell1;
    Delaunay::Full_cell_handle cell2;

    bool operator==(const Voronoi_edge & other) const {
        return (cell1 == other.cell1 && cell2 == other.cell2) ||
            (cell1 == other.cell2 && cell2 == other.cell1);
    }

    Voronoi_edge(Voronoi_vertex vertex1, Voronoi_vertex vertex2, Delaunay::Full_cell_handle cell1, Delaunay::Full_cell_handle cell2) : vertex1(vertex1), vertex2(vertex2), cell1(cell1), cell2(cell2) {}
};

struct Voronoi_face {
    std::vector<Voronoi_vertex> voronoi_vertices;
    std::vector<Voronoi_edge> voronoi_edges;
    std::vector<Eigen::VectorXd> ps_vertices;
    std::vector<std::array<size_t, 2>> ps_edges;
};

Voronoi_face delaunay_edge_dual(Edge &edge, Face &df, Delaunay &dt);

Voronoi_edge delaunay_face_dual(Face &face, Delaunay &dt);

Face voronoi_edge_dual(Voronoi_edge &voronoi_edge);

Eigen::VectorXd simplex_circumsphere(Delaunay::Full_cell_handle simplex);

void insert_points(std::vector<Point> &points, Delaunay &delaunay);

bool is_index_two_critical_point(Face &face, Delaunay &dt);

bool is_gabriel(Edge &edge, Delaunay &dt);

void get_incident_cells_to_vertices(Edge &edge, Delaunay &dt, std::vector<Delaunay::Full_cell_handle> &cell_neighbors);

void get_facet_vertices(
    const Delaunay &delaunay,
    const Delaunay::Facet_iterator &facet,
    std::vector<Delaunay::Vertex_handle> &facet_vertices
);

void get_facet_normal(std::vector<Eigen::VectorXd> &facet_points, Eigen::VectorXd &normal);

void orient_voronoi_edge(std::vector<Eigen::VectorXd> shared_facet_points, Eigen::VectorXd finite_voronoi_vertex, Eigen::VectorXd &voronoi_edge_direction);

void extract_edges(Delaunay &delaunay, std::map<Delaunay::Vertex_handle, size_t> vertex_to_index, std::vector<std::array<size_t, 2>> &edges);

void calculate_driver(const Eigen::VectorXd &voronoi_vertex, const Edge &delaunay_edge, Eigen::VectorXd &driver);


#endif
