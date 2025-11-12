#ifndef DELAUNAY_H
#define DELAUNAY_H

#include <CGAL/Epick_d.h>
#include <CGAL/Delaunay_triangulation.h>

typedef CGAL::Dimension_tag<3> Dimension;
typedef CGAL::Epick_d<Dimension> K;
typedef K::Point_d Point;
typedef K::FT FT;

typedef CGAL::Delaunay_triangulation<K> Delaunay;

struct Edge {
    Delaunay::Vertex_handle vertex1;
    Delaunay::Vertex_handle vertex2;

    Edge(Delaunay::Vertex_handle vertex1, Delaunay::Vertex_handle vertex2) : vertex1(vertex1), vertex2(vertex2) {}

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

inline bool lexicographic_less_vector(const Eigen::VectorXd &a, const Eigen::VectorXd &b) {
    int dimension = a.size();
    for (int i = 0; i < dimension; i++) {
        if (a[i] < b[i]) { return true; }
        if (a[i] > b[i]) { return false; }
    }
}

struct Voronoi_edge {
    Eigen::VectorXd vertex1;
    Eigen::VectorXd vertex2;

    Delaunay::Full_cell_handle cell1;
    Delaunay::Full_cell_handle cell2;


    Voronoi_edge(Eigen::VectorXd vertex1, Eigen::VectorXd vertex2, Delaunay::Full_cell_handle cell1, Delaunay::Full_cell_handle cell2) : vertex1(vertex1), vertex2(vertex2), cell1(cell1), cell2(cell2) {}

    bool operator==(const Voronoi_edge & other) const {
        return (vertex1 == other.vertex1 && vertex2 == other.vertex2) ||
            (vertex1 == other.vertex2 && vertex2 == other.vertex1);
    }

    bool operator<(const Voronoi_edge& other) const {
        auto this_canonical = lexicographic_less_vector(vertex1, vertex2) ? std::make_pair(vertex1, vertex2) : std::make_pair(vertex2, vertex1);
        auto other_canonical = lexicographic_less_vector(other.vertex1, other.vertex2) ? std::make_pair(other.vertex1, other.vertex2) : std::make_pair(other.vertex2, other.vertex1);
        if (this_canonical.first == other_canonical.first) {
            return lexicographic_less_vector(this_canonical.second, other_canonical.second);
        }
        return lexicographic_less_vector(this_canonical.first, other_canonical.first);
    }
};

struct Voronoi_face {
    std::vector<Voronoi_edge> voronoi_edges;
};

void insert_points(std::vector<Point> &points, Delaunay &delaunay);

bool is_index_two_critical_point(const std::vector<Delaunay::Vertex_handle> &facet_vertices);

bool is_gabriel(Edge &edge);

void get_incident_cells_to_vertices(Edge &edge, std::vector<Delaunay::Full_cell_handle> &cell_neighbors);

void voronoi_facet_from_edge(Edge &edge, Voronoi_face &voronoi_face, Delaunay &delaunay);

void get_voronoi_facet_vertices(const Edge& edge, const Delaunay& delaunay, std::vector<Eigen::VectorXd>& facet_vertices);

void get_facet_vertices(
    const Delaunay &delaunay,
    const Delaunay::Facet_iterator &facet,
    std::vector<Delaunay::Vertex_handle> &facet_vertices
);

void extract_edges(Delaunay &delaunay, std::map<Delaunay::Vertex_handle, size_t> vertex_to_index, std::vector<std::array<size_t, 2>> &edges);

void calculate_driver(const Eigen::VectorXd &voronoi_vertex, const Edge &delaunay_edge, Eigen::VectorXd &driver);

bool is_intersection_in_facet(
    const Eigen::VectorXd& intersection_point,
    const std::vector<Eigen::VectorXd>& facet_vertices);

#endif
