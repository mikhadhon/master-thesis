#include <iostream>
#include <polyscope/polyscope.h>
#include <polyscope/surface_mesh.h>

#include "Delaunay.h"
#include "FlowComplex.h"
#include "utils.h"
#include "polyscope/curve_network.h"
#include "polyscope/point_cloud.h"


int main() {
    double samples[][3] = {
        {0, 0, 0},
        {2, 0, 0},
        // {1, 2, 0},
        // {1, 0.5, 1}
    };

    Eigen::Matrix<double, Eigen::Dynamic, 3> V;
    Eigen::Matrix<double, Eigen::Dynamic, 3> N1;
    Eigen::Matrix<double, Eigen::Dynamic, 3> N2;
    Eigen::Matrix<double, Eigen::Dynamic, 3> T;

    // generate_circle(100, 0, V, N1, N2, T);

    Delaunay delaunay3D(3);

    std::vector<Delaunay::Point> points;
    //gen_rectangle(10, points);
    // generate_torus(1000, 0.5, 1, points);
    gen_sphere_sample(5000, 1, points);


    // for (int i = 0; i < V.rows(); i++) {
    //     Delaunay::Point p(V(i, 0), V(i, 1), V(i, 2));
    //     points.push_back(p);
    // }

    // for (auto & sample : samples) {
    //     Delaunay::Point p(&sample[0], &sample[3]);
    //     points.push_back(p);
    // }

    insert_points(points, delaunay3D);

    std::cout << "Triangulation dimension: " << delaunay3D.current_dimension() << std::endl;

    polyscope::init();

    std::vector<std::array<double, 3> > vertices;
    std::vector<std::array<size_t, 3> > faces;
    std::vector<std::array<size_t, 3> > reconstructed_faces;
    std::vector<std::array<size_t, 2> > edges;
    std::map<Delaunay::Vertex_handle, size_t> vertex_to_index;

    map_vertices_to_vector(delaunay3D, vertices, vertex_to_index);
    write_faces_to_vector(delaunay3D, faces, vertex_to_index);
    extract_edges(delaunay3D, vertex_to_index, edges);
    // auto *psMesh = polyscope::registerSurfaceMesh("delaunay mesh", vertices, faces);
    auto *pcCloud = polyscope::registerPointCloud("vertices", vertices);

    std::vector<Eigen::Vector3d> centers;
    std::vector<Eigen::Vector3d> centers2;
    flow_complex(delaunay3D, vertices, reconstructed_faces, vertex_to_index);

    std::cout << "Number of vertices: " << vertices.size() << std::endl;
    std::cout << "Number of faces: " << faces.size() << std::endl;
    std::cout << "Number of edges: " << edges.size() << std::endl;

    if (!faces.empty()) {
        auto *psRecMesh = polyscope::registerSurfaceMesh("rec mesh", vertices, reconstructed_faces);
        // auto *pcCloud = polyscope::registerPointCloud("origins", centers);
        // auto *pcCloud2 = polyscope::registerPointCloud("s", centers2);
    }
    if (!edges.empty()) {
        auto *psCurve = polyscope::registerCurveNetwork("delaunay vertices", vertices, edges);
    }

    polyscope::show();

    return 0;
}
