#include <iostream>
#include <polyscope/polyscope.h>
#include <polyscope/surface_mesh.h>

#include "Delaunay.h"
#include "FlowComplex.h"
#include "utils.h"
#include "data_gen.h"
#include "polyscope/curve_network.h"
#include "polyscope/point_cloud.h"
#include "igl/random_points_on_mesh.h"

void load_obj_points(const std::string& filename, std::vector<Delaunay::Point>& points) {
    std::ifstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Error: Could not open file " << filename << std::endl;
        return;
    }
    std::string line;
    while (std::getline(file, line)) {
        if (line.substr(0, 2) == "v ") {
            std::istringstream s(line.substr(2));
            double x, y, z;
            s >> x >> y >> z;
            points.emplace_back(x, y, z);
        }
    }
    file.close();
}

void load_xyz_points(const std::string& filename, std::vector<Delaunay::Point>& points) {
    std::ifstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Error: Could not open file " << filename << std::endl;
        return;
    }
    std::string line;
    while (std::getline(file, line)) {
        std::istringstream s(line);
        double x, y, z;
        s >> x >> y >> z;
        points.emplace_back(x, y, z);
    }
    file.close();
}

int main() {
    std::string testfile = "happy_vrip_res";
    std::ofstream timestamps("../output/" + testfile + "_" + get_timestamp() + ".csv");

    timestamps << "Vertices,Duration\n";

    std::vector nvertices = {1000000};
    // for (int i = 1; i >= 1; i--) {
    for (auto n : nvertices) {
        Eigen::MatrixXi F;
        Eigen::Matrix<double, Eigen::Dynamic, 3> N1;
        Eigen::Matrix<double, Eigen::Dynamic, 3> N2;
        Eigen::Matrix<double, Eigen::Dynamic, 3> T;

        // Eigen::Matrix<double, Eigen::Dynamic, 3> V;
        // generate_trefoil(10000, 3, V, N1, N2, T);
        std::string filename = testfile + std::to_string(n);
        // std::string filepath = "../" + filename + ".ply";
        //std::string filepath = "../output/happy_vrip_res1000000_1000000_20260204_213045.obj";
        std::string filepath = "../happy_vrip_res1.ply";
        Eigen::MatrixXd V;
        read_ply(filepath, V, F);
        //read_obj(filepath, V, F);
        std::vector<Delaunay::Point> points;
        //gen_sphere_sample(10000, 1, points);
        //generate_torus(1000, 0.5, 1, points);


        Eigen::MatrixXd B;
        Eigen::MatrixXi FI;
        Eigen::MatrixXd X;
        igl::random_points_on_mesh(n, V, F, B, FI, X);

        Delaunay delaunay3D(3);

        //load_obj_points("stanford-bunny.obj", points);
        //load_xyz_points("../bunny_raw.xyz", points);

        for (int i = 0; i < X.rows(); i++) {
            Delaunay::Point p(X(i, 0), X(i, 1), X(i, 2));
            points.push_back(p);
        }
        // for (int i = 0; i < V.rows(); i++) {
        //     Delaunay::Point p(V(i, 0), V(i, 1), V(i, 2));
        //     points.push_back(p);
        // }

        insert_points(points, delaunay3D);

        std::cout << "Triangulation dimension: " << delaunay3D.current_dimension() << std::endl;
        std::cout << "Vertices: " << points.size() << std::endl;

        polyscope::init();
        polyscope::options::groundPlaneMode = polyscope::GroundPlaneMode::None;

        std::vector<Eigen::VectorXd> vertices;
        std::vector<std::array<size_t, 3> > faces;
        std::vector<std::array<size_t, 3> > reconstructed_faces;
        std::map<Delaunay::Vertex_handle, size_t> vertex_to_index;

        map_vertices_to_vector(delaunay3D, vertices, vertex_to_index);
        write_faces_to_vector(delaunay3D, faces, vertex_to_index);
        //auto *psMesh = polyscope::registerSurfaceMesh("delaunay mesh", vertices, faces);
        auto *pcCloud = polyscope::registerPointCloud("delaunay vertices", vertices);


        std::vector<Eigen::VectorXd> reconstructed_vertices = vertices;
        auto time_start = std::chrono::high_resolution_clock::now();
        std::vector<std::vector<std::array<size_t, 3>>> index_2_stable_manifolds;
        flow_complex(delaunay3D, reconstructed_vertices, reconstructed_faces, vertex_to_index, index_2_stable_manifolds);
        auto time_end = std::chrono::high_resolution_clock::now();

        std::chrono::duration<double, std::milli> float_ms = time_end - time_start;

        std::cout << "Elapsed time: " << float_ms.count() << " ms" << std::endl;
        timestamps << points.size() << "," << float_ms.count() << std::endl;

        std::vector<Edge> gabriel_edges;
        //find_gabriel_edges(delaunay3D, gabriel_edges);

        std::vector<std::array<size_t, 2>> ps_edges;
        for (auto edge: gabriel_edges) {
            ps_edges.push_back({vertex_to_index[edge.vertex1], vertex_to_index[edge.vertex2]});
        }

        //create reduced flow complex
        //find pairs of gabriel edge + index 2 stable manifolds where stable manifold has gabriel edge incident to 2 stable manifolds


        std::cout << "Number of vertices: " << reconstructed_vertices.size() << std::endl;
        std::cout << "Number of faces: " << faces.size() << std::endl;

        if (!reconstructed_faces.empty()) {
            auto *psRecMesh = polyscope::registerSurfaceMesh("reconstructed mesh", reconstructed_vertices, reconstructed_faces);

            Eigen::MatrixXd V_out(reconstructed_vertices.size(), 3);
            Eigen::MatrixXi F_out(reconstructed_faces.size(), 3);

            for (int i = 0; i < reconstructed_vertices.size(); i++) {
                V_out.row(i) = reconstructed_vertices[i];
            }

            for (int i = 0; i < reconstructed_faces.size(); i++) {
                Eigen::VectorXi face(3);
                face << reconstructed_faces[i][0], reconstructed_faces[i][1], reconstructed_faces[i][2];
                F_out.row(i) = face;
            }

            write_to_obj(V_out, F_out, filename + "_" + std::to_string(points.size()));
        }

        if (!ps_edges.empty()) {
            polyscope::registerCurveNetwork("Gabriel edges", vertices, ps_edges);

            Eigen::MatrixXd V_out(vertices.size(), 3);
            Eigen::MatrixXi F_out(ps_edges.size(), 2);

            for (int i = 0; i < vertices.size(); i++) {
                V_out.row(i) = vertices[i];
            }

            for (int i = 0; i < ps_edges.size(); i++) {
                Eigen::VectorXi edge(2);
                edge << ps_edges[i][0], ps_edges[i][1];
                F_out.row(i) = edge;
            }

            write_to_obj(V_out, F_out, filename + "_" + std::to_string(points.size()));
        }
    }

    timestamps.close();
    polyscope::show();

    return 0;
}
