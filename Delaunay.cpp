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

void indexTwoCriticalPoint(const Delaunay * delaunay, const Delaunay::Facet_iterator &facet, Point * criticalPoint) {
    std::cout << "----------------------" << std::endl;
    std::vector<Point> facet_vertices;
    getFacetVertices(delaunay, facet, &facet_vertices);

    std::cout << "facet_vertices: " << facet_vertices.size() << std::endl;
    if (facet_vertices.size() != facet->full_cell()->maximal_dimension()) {
        return;
    }
    //calc critical point
    std::cout << "calc critical point" << std::endl;
}

void getFacetVertices(const Delaunay * delaunay, const Delaunay::Facet_iterator &facet, std::vector<Point> * facet_vertices) {
    const int dimension = facet->full_cell()->maximal_dimension();
    const int co_vertex_index = facet->index_of_covertex();
    const Delaunay::Full_cell_const_handle full_cell = facet->full_cell();

    /*
    if (delaunay->is_infinite(facet->full_cell()->vertex(co_vertex_index))) {
        std::cout << "co_vertex is infinite" << std::endl;
    }

    int vertex_count = 0;
    int total = 0;
    for (auto v = facet->full_cell()->vertices_begin(); v != facet->full_cell()->vertices_end(); ++v) {
        if (*v != nullptr) vertex_count++;
        total++;
    }
    std::cout << "Number of non-null vertices in full_cell: " << vertex_count << std::endl;
    std::cout << "Total number of vertices in full_cell: " << total << std::endl;
    */

    for (int i = 0; i <= dimension; ++i) {
        if (i != co_vertex_index) {
            const auto vertex_handle = facet->full_cell()->vertex(i);

            if (vertex_handle != nullptr && !delaunay->is_infinite(vertex_handle)) {
                facet_vertices->push_back(vertex_handle->point());
            }
        }
    }
}