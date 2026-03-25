#include <iostream>
#include <polyscope/polyscope.h>
#include <polyscope/surface_mesh.h>

#include "FlowComplex.h"
#include "utils.h"
#include "polyscope/point_cloud.h"
#include "igl/random_points_on_mesh.h"


int main(int argc, char *argv[]) {
    const int dim = 4;

    Eigen::MatrixXd V;
    Eigen::MatrixXi F;
    std::string filename = "output";
    if (argc > 1) {
        filename = argv[1];
        read_obj(argv[1], V, F);
        std::cout << "Loaded " << V.rows() << " points" << std::endl;
    } else {
        V = cliffordgen(1000);
        // write_to_obj(V, F, "../million_ground_truth.obj");
        // read_obj("../hopf_torus.obj", V, F);
    }

    std::vector<Point> points;
    for (int i = 0; i < V.rows(); i++) {
        Point p(V(i, 0), V(i, 1), V(i, 2), V(i, 3));
        points.push_back(p);
    }

    Delaunay dt(dim);
    insert_points(points, dt);

    std::cout << dt.current_dimension() << std::endl;

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
    pcCloud->setEnabled(false);

    std::vector<Eigen::VectorXd> reconstructed_vertices = vertices;
    int vertex_count = static_cast<int>(vertices.size());
    std::map<std::array<double, 4>, size_t> fc_vertex_to_index;
    std::vector<std::array<size_t, 2>> gabriel_edges;
    std::vector<stable_manifold_2> stable_manifolds;
    std::unordered_map<fc_edge, int, fc_edge_hash> edge_incidences;
    FT threshold = 0.06;
    std::vector<Point> index_2;

    helper_data data = {vertex_count, vertex_to_index, fc_vertex_to_index, gabriel_edges, stable_manifolds, edge_incidences, threshold, index_2};

    const auto time_start = std::chrono::high_resolution_clock::now();
    flow_complex(dt, reconstructed_vertices, reconstructed_faces, data);
    const auto time_end = std::chrono::high_resolution_clock::now();

    const std::chrono::duration<double, std::milli> float_ms = time_end - time_start;

    std::cout << "Elapsed time: " << float_ms.count() << " ms" << std::endl;

    std::cout << "Number of vertices: " << reconstructed_vertices.size() << std::endl;
    std::cout << "Number of faces: " << faces.size() << std::endl;

    std::vector<Eigen::VectorXd> flagged_vertices;
    std::unordered_map<fc_edge, std::vector<stable_manifold_2>, fc_edge_hash> more_postprocessing;

    for (auto [e, count] : data.edge_incidences) {
        if (count != 2) {
            flagged_vertices.push_back(Eigen::Vector4d(CGAL::to_double(e.vertex1.cartesian(0)), CGAL::to_double(e.vertex1.cartesian(1)), CGAL::to_double(e.vertex1.cartesian(2)), CGAL::to_double(e.vertex1.cartesian(3))));
            flagged_vertices.push_back(Eigen::Vector4d(CGAL::to_double(e.vertex2.cartesian(0)), CGAL::to_double(e.vertex2.cartesian(1)), CGAL::to_double(e.vertex2.cartesian(2)), CGAL::to_double(e.vertex2.cartesian(3))));
        }
    }

    std::vector<std::array<size_t, 3> > reduced_faces;
    for (auto sm : data.stable_manifolds) {
        for (auto fc_face : sm.faces) {
            reduced_faces.push_back(fc_face);
        }
    }

    const std::vector<Eigen::VectorXd> V_reconstructed = stereo_projection(reconstructed_vertices);
    polyscope::registerPointCloud("reconstructed vertices", V_reconstructed)->setEnabled(false);

    std::vector<Eigen::VectorXd> index_2_m;
    for (int i = 0; i < data.index_2.size(); i++) {
        index_2_m.push_back(make_point_eigen(data.index_2[i]));
    }
    std::vector<Eigen::VectorXd> projected_index_2 = stereo_projection(index_2_m);
    polyscope::registerPointCloud("index_2", projected_index_2)->setEnabled(false);

    std::vector<Eigen::VectorXd> flagged_vertices_projected = stereo_projection(flagged_vertices);
    polyscope::registerPointCloud("flagged vertices", flagged_vertices_projected);

    if (!reconstructed_faces.empty()) {
        auto *psRecMesh = polyscope::registerSurfaceMesh("reconstructed mesh", V_reconstructed, reconstructed_faces);
        psRecMesh->setEnabled(false);
        polyscope::registerSurfaceMesh("reduced mesh", V_reconstructed, reduced_faces);

        Eigen::MatrixXd V_out(reconstructed_vertices.size(), 4);
        Eigen::MatrixXi F_out(reduced_faces.size(), 3);

        for (int i = 0; i < reconstructed_vertices.size(); i++) {
            V_out.row(i) = reconstructed_vertices[i];
        }
        for (int i = 0; i < reduced_faces.size(); i++) {
            Eigen::VectorXi face(3);
            face << reduced_faces[i][0], reduced_faces[i][1], reduced_faces[i][2];
            F_out.row(i) = face;
        }
        write_to_obj(V_out, F_out, filename + "_" + std::to_string(points.size()));
    }

    polyscope::show();

    return 0;
}
