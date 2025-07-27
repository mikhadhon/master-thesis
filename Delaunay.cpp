#include <CGAL/Linear_algebraCd.h>

#include "Delaunay.h"

typedef CGAL::Linear_algebraCd<FT> Linear_algebra;
typedef Linear_algebra::Matrix Matrix;
typedef Linear_algebra::Vector Vector;
typedef Linear_algebra::RT RT;

void insert_points(std::vector<Point> &points, Delaunay &delaunay) {
    Delaunay::Vertex_handle hint;
    int i = 0;
    for (auto it = points.begin(); it != points.end(); ++it) {
        if (Delaunay::Vertex_handle() != hint) {
            hint = delaunay.insert(*it, hint);
        }
        else {
            hint = delaunay.insert(*it);
        }
        printf("Processing: %d/%d\n", ++i, static_cast<int>(points.size()));
    }
    if (!delaunay.is_valid()) {
        std::cerr << "Triangulation is invalid!" << std::endl;
    }
}

bool is_index_two_critical_point(const std::vector<Delaunay::Vertex_handle> &facet_vertices) {

    Eigen::Vector3d i(facet_vertices[0]->point()[0], facet_vertices[0]->point()[1], facet_vertices[0]->point()[2]);
    Eigen::Vector3d j(facet_vertices[1]->point()[0], facet_vertices[1]->point()[1], facet_vertices[1]->point()[2]);
    Eigen::Vector3d l(facet_vertices[2]->point()[0], facet_vertices[2]->point()[1], facet_vertices[2]->point()[2]);

    return i.dot(j) > 0 && i.dot(l) > 0 && j.dot(l) > 0;
}

bool is_gabriel(Edge) {

}

void get_facet_vertices(const Delaunay &delaunay, const Delaunay::Facet_iterator &facet,
    std::vector<Delaunay::Vertex_handle> &facet_vertices
    ) {
    const int dimension = facet->full_cell()->maximal_dimension();
    const int co_vertex_index = facet->index_of_covertex();

    for (int i = 0; i <= dimension; ++i) {
        if (i != co_vertex_index) {
            if (const auto vertex_handle = facet->full_cell()->vertex(i);
                vertex_handle != nullptr && !delaunay.is_infinite(vertex_handle)
            ) {
                facet_vertices.push_back(vertex_handle);
            }
        }
    }
}

void extract_edges(Delaunay &delaunay, std::map<Delaunay::Vertex_handle, size_t> vertex_to_index, std::vector<std::array<size_t, 2>> &edges) {
    std::set<Edge> unique_edges;
    for (auto facet = delaunay.facets_begin(); facet != delaunay.facets_end(); ++facet) {
        if (!delaunay.is_infinite(*facet)) {
            std::vector<Delaunay::Vertex_handle> face_vertices;
            get_facet_vertices(delaunay, facet, face_vertices);
            if (face_vertices.size() == 2) {
                unique_edges.insert(Edge(face_vertices[0], face_vertices[1]));
            }
        }
    }
    for (const auto &edge : unique_edges) {
        edges.push_back(std::array{vertex_to_index[edge.vertex1], vertex_to_index[edge.vertex2]});
    }
}
