#ifndef DELAUNAY_H
#define DELAUNAY_H

#include <CGAL/Epeck_d.h>
#include <CGAL/Delaunay_triangulation.h>
#include <boost/functional/hash.hpp>

typedef CGAL::Dimension_tag<4> Dimension;
typedef CGAL::Epeck_d<Dimension> K;
typedef K::Point_d Point;
typedef K::Vector_d Vector;
typedef K::FT FT;

typedef CGAL::Delaunay_triangulation<K> Delaunay;
typedef K::Midpoint_d midpoint;
typedef K::Construct_circumcenter_d circumcenter;
typedef K::Squared_distance_d squared_distance;
typedef K::Contained_in_simplex_d contained_in_simplex;
typedef K::Scalar_product_d dot;
typedef K::Equal_d equal;

struct Edge {
    Delaunay::Vertex_handle vertex1;
    Delaunay::Vertex_handle vertex2;

    Delaunay::Vertex_handle co_vertex;

    Edge(const Delaunay::Vertex_handle vertex1, const Delaunay::Vertex_handle vertex2, const Delaunay::Vertex_handle co_vertex) : vertex1(vertex1), vertex2(vertex2), co_vertex(co_vertex) {}

    bool operator==(const Edge & other) const {
        return (vertex1 == other.vertex1 && vertex2 == other.vertex2) ||
            (vertex1 == other.vertex2 && vertex2 == other.vertex1);
    }

    bool operator<(const Edge& other) const {
        const std::pair<Delaunay::Vertex_handle, Delaunay::Vertex_handle> this_canonical = (vertex1 < vertex2) ? std::make_pair(vertex1, vertex2) : std::make_pair(vertex2, vertex1);
        const std::pair<Delaunay::Vertex_handle, Delaunay::Vertex_handle> other_canonical = (other.vertex1 < other.vertex2) ? std::make_pair(other.vertex1, other.vertex2) : std::make_pair(other.vertex2, other.vertex1);
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
    int index_of_covertex_0;
    int index_of_covertex_1;

    bool operator==(const Face & other) const {
        return (face.vertex(0) == other.face.vertex(0) && face.vertex(1) == other.face.vertex(1)) && face.vertex(2) == other.face.vertex(2);
    }

    Face(const Delaunay::Face &face, const int index_of_covertex_0, const int index_of_covertex_1) : face(face), index_of_covertex_0(index_of_covertex_0), index_of_covertex_1(index_of_covertex_1) {}
};

struct FaceHash {
    std::size_t operator()(const Face& f) const {
        std::size_t seed = 0;
        boost::hash_combine(seed, &*f.face.vertex(0));
        boost::hash_combine(seed, &*f.face.vertex(1));
        boost::hash_combine(seed, &*f.face.vertex(2));
        return seed;
    }
};

struct Voronoi_vertex {
    bool is_infinite;
    Point point;
    Vector infinite_direction;

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

    Voronoi_edge(const Voronoi_vertex &vertex1, const Voronoi_vertex &vertex2, const Delaunay::Full_cell_handle cell1, const Delaunay::Full_cell_handle cell2) : vertex1(vertex1), vertex2(vertex2), cell1(cell1), cell2(cell2) {}
};

struct Voronoi_face {
    std::vector<Voronoi_vertex> voronoi_vertices;
    std::vector<Voronoi_edge> voronoi_edges;
    std::vector<Eigen::VectorXd> ps_vertices;
    std::vector<std::array<size_t, 2>> ps_edges;
    Face dual;
};

std::vector<Voronoi_face> delaunay_edge_dual(const Edge &edge, const Face &df, const Delaunay &dt);

Voronoi_face delaunay_face_dual(const Face &face, const Delaunay &dt);

Eigen::VectorXd simplex_circumsphere(const Delaunay::Full_cell_handle &simplex);

void insert_points(const std::vector<Point> &points, Delaunay &delaunay);

bool is_index_two_critical_point(const Face &face, const Delaunay &dt);

bool is_gabriel(const Edge &edge, const Delaunay &dt);

void get_incident_cells_to_edge(const Edge &edge, const Delaunay &dt, std::vector<Delaunay::Full_cell_handle> &cell_neighbors);

void get_incident_cells_to_face(const Face &face, const Delaunay &dt, std::vector<Delaunay::Full_cell_handle> &cell_neighbors);

void get_facet_vertices(
    const Delaunay &delaunay,
    const Delaunay::Facet_iterator &facet,
    std::vector<Delaunay::Vertex_handle> &facet_vertices
);

void get_tetrahedron_normal(const std::vector<Eigen::VectorXd> &tetrahedron_points, Eigen::VectorXd &normal);

std::unordered_set<Face, FaceHash> get_delaunay_faces(Delaunay &dt);

#endif
