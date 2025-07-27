#ifndef DELAUNAY_H
#define DELAUNAY_H

#include <CGAL/Epick_d.h>
#include <CGAL/Delaunay_triangulation.h>

typedef CGAL::Dimension_tag<3> Dimension;
typedef CGAL::Epick_d<Dimension> K;
typedef K::Point_d Point;
typedef K::FT FT;

typedef CGAL::Delaunay_triangulation<CGAL::Epick_d<Dimension>> Delaunay;

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

void insert_points(std::vector<Point> &points, Delaunay &delaunay);

bool is_index_two_critical_point(const std::vector<Delaunay::Vertex_handle> &facet_vertices);

bool is_gabriel(Edge &edge);

void get_incident_cells(Edge &edge, std::vector<Delaunay::Full_cell_handle> &cell_neighbors);

void get_facet_vertices(
    const Delaunay &delaunay,
    const Delaunay::Facet_iterator &facet,
    std::vector<Delaunay::Vertex_handle> &facet_vertices
);

void extract_edges(Delaunay &delaunay, std::map<Delaunay::Vertex_handle, size_t> vertex_to_index, std::vector<std::array<size_t, 2>> &edges);

#endif
