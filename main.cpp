#include <iostream>
#include <polyscope/polyscope.h>
#include <polyscope/surface_mesh.h>

#include "Delaunay.h"
#include "utils.h"

int main() {
    double samples[][3] = {
        {0, 0, 0},
        {2, 0, 0},
        {1, 2, 0},
        {1, 0.5, 1}
    };


    Delaunay delaunay3D(3);

    std::vector<Delaunay::Point> points;
    points.reserve(13);

    for (auto & sample : samples) {
        Delaunay::Point p(&sample[0], &sample[3]);
        points.push_back(p);
    }

    insertPoints(points, delaunay3D);

    std::cout << "Triangulation dimension: " << delaunay3D.current_dimension() << std::endl;

    polyscope::init();

    std::vector<std::array<double, 3> > vertices;
    std::vector<std::array<size_t, 3> > faces;
    std::map<Delaunay::Vertex_handle, size_t> vertex_to_index;

    map_vertices_to_vector(delaunay3D, vertices, vertex_to_index);
    write_faces_to_vector(delaunay3D, faces, vertex_to_index);

    std::cout << "Number of vertices: " << vertices.size() << std::endl;
    std::cout << "Number of faces: " << faces.size() << std::endl;

    auto *psMesh = polyscope::registerSurfaceMesh("delaunay mesh", vertices, faces);

    polyscope::show();

    Point criticalPoint;
    for (auto facet = delaunay3D.facets_begin(); facet != delaunay3D.facets_end(); ++facet) {
        if (!delaunay3D.is_infinite(*facet)) {
            indexTwoCriticalPoint(&delaunay3D, facet, &criticalPoint);
        }
    }

    return 0;
}
