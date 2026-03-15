#include <iostream>
#include <polyscope/polyscope.h>
#include <polyscope/surface_mesh.h>

#include "FlowComplex.h"
#include "utils.h"
#include "polyscope/point_cloud.h"
#include "igl/random_points_on_mesh.h"


int main() {
    const int dim = 4;

    Eigen::MatrixXd V;
    Eigen::MatrixXi F;
    // V = cliffordgen(1000);
    //V = sphere2gen(1000);
    // V = cliffordgengrid(gen_rectangle(1000));
    read_obj("../hopf_torus.obj", V, F);

    std::vector<Point> points;
    for (int i = 0; i < V.rows(); i++) {
        Point p(V(i, 0), V(i, 1), V(i, 2), V(i, 3));
        points.push_back(p);
    }

    Delaunay dt(dim);
    insert_points(points, dt);

    std::vector<Eigen::VectorXd> vertices;
    std::vector<std::array<size_t, 3> > faces;
    std::vector<std::array<size_t, 3> > reconstructed_faces;
    std::map<Delaunay::Vertex_handle, size_t> vertex_to_index;

    map_vertices_to_vector(dt, vertices, vertex_to_index);
    write_faces_to_vector(dt, faces, vertex_to_index);
    Eigen::MatrixXd V_projected = stereo_projection(V);

    //write_to_obj(V, F, "cliffordgen.obj");
    polyscope::init();
    // auto *psMesh = polyscope::registerSurfaceMesh("delaunay mesh", V_projected, faces);
    auto *pcCloud = polyscope::registerPointCloud("delaunay vertices", V_projected);

    std::vector<Eigen::VectorXd> reconstructed_vertices = vertices;
    int vertex_count = static_cast<int>(vertices.size());
    std::map<std::array<double, 4>, size_t> fc_vertex_to_index;
    std::vector<std::array<size_t, 2>> gabriel_edges;
    std::vector<stable_manifold_2> sm;

    helper_data data = {vertex_count, vertex_to_index, fc_vertex_to_index, gabriel_edges, sm};

    const auto time_start = std::chrono::high_resolution_clock::now();
    flow_complex(dt, reconstructed_vertices, reconstructed_faces, data);
    const std::vector<Eigen::VectorXd> V_reconstructed = stereo_projection(reconstructed_vertices);
    polyscope::registerPointCloud("index 2", V_reconstructed);
    const auto time_end = std::chrono::high_resolution_clock::now();

    const std::chrono::duration<double, std::milli> float_ms = time_end - time_start;

    std::cout << "Elapsed time: " << float_ms.count() << " ms" << std::endl;

    std::cout << "Number of vertices: " << reconstructed_vertices.size() << std::endl;
    std::cout << "Number of faces: " << faces.size() << std::endl;

    if (!reconstructed_faces.empty()) {
        auto *psRecMesh = polyscope::registerSurfaceMesh("reconstructed mesh", V_reconstructed, reconstructed_faces);
    }

    polyscope::show();

    return 0;
}
