#include <iostream>
#include <polyscope/polyscope.h>
#include <polyscope/surface_mesh.h>

#include "Delaunay.h"
#include "FlowComplex.h"
#include "utils.h"
#include "polyscope/curve_network.h"
#include "polyscope/point_cloud.h"

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
    // gen_sphere_sample(1000, 1, points);

    load_obj_points("bunny.obj", points);

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
    std::map<Delaunay::Vertex_handle, size_t> vertex_to_index;

    map_vertices_to_vector(delaunay3D, vertices, vertex_to_index);
    write_faces_to_vector(delaunay3D, faces, vertex_to_index);
    auto *psMesh = polyscope::registerSurfaceMesh("delaunay mesh", vertices, faces);
    //auto *pcCloud = polyscope::registerPointCloud("delaunay vertices", vertices);

    std::vector<std::array<double, 3> > reconstructed_vertices = vertices;
    flow_complex(delaunay3D, reconstructed_vertices, reconstructed_faces, vertex_to_index);

    std::cout << "Number of vertices: " << reconstructed_vertices.size() << std::endl;
    std::cout << "Number of faces: " << faces.size() << std::endl;

    if (!faces.empty()) {
        auto *psRecMesh = polyscope::registerSurfaceMesh("reconstructed mesh", reconstructed_vertices, reconstructed_faces);
    }

    polyscope::show();

    return 0;
}
