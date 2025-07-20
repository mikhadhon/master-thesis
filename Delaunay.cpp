#include <CGAL/Linear_algebraCd.h>

#include "Delaunay.h"

typedef CGAL::Linear_algebraCd<FT> Linear_algebra;
typedef Linear_algebra::Matrix Matrix;
typedef Linear_algebra::Vector Vector;
typedef Linear_algebra::RT RT;

void insertPoints(std::vector<Point> &points, Delaunay &delaunay) {
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

void indexTwoCriticalPoint(const Delaunay &delaunay, const Delaunay::Facet_iterator &facet, Point &criticalPoint) {
    std::vector<Delaunay::Vertex_handle> facet_vertices;
    getFacetVertices(delaunay, facet, facet_vertices);

    //calc critical point
    std::cout << "calc critical point" << std::endl;
}

void getFacetVertices(const Delaunay &delaunay, const Delaunay::Facet_iterator &facet,
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