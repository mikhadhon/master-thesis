#include <iostream>
#include <polyscope/polyscope.h>
#include <polyscope/surface_mesh.h>

#include "Delaunay.h"
#include "FlowComplex.h"
#include "utils.h"
#include "polyscope/curve_network.h"
#include "polyscope/point_cloud.h"
#include "igl/random_points_on_mesh.h"

Eigen::MatrixXd cliffordgengrid(Eigen::MatrixXd samples) {
    int size = int(sqrt(samples.rows())) - 1;
    double radius = 0.7;

    Eigen::MatrixXd torus_samples(0, 4);

    for (int i = 0; i < samples.rows() - (2 * size + 1); i++) {
        double q = samples.row(i)(0) * 2 * M_PI;
        double p = samples.row(i)(1) * 2 * M_PI;
        torus_samples.conservativeResize(torus_samples.rows() + 1, Eigen::NoChange);
        torus_samples.row(torus_samples.rows() - 1) = Eigen::RowVector4d(radius * cos(q), radius * sin(q), radius * cos(p), radius * sin(p));
    }

    return torus_samples;
}

Eigen::MatrixXd gen_rectangle(int n) {
    std::random_device rd{};
    std::mt19937 gen{ rd() };

    std::uniform_real_distribution<> noise{ -0.00000001, 0.000000001 };

    int size = n;

    Eigen::MatrixXd samples((size + 1) * (size + 1), 4);

    for (int i = 0; i < size; i++) {
        for (int j = 0; j < size; j++) {
            double x = -1, y = -1;
            if (i == 0 || j == 0) {
                x = i * (1 / (double)size);
                y = j * (1 / (double)size);
            }
            else {
                while (x > 1 || x < 0) {
                    x = i * (1 / (double)size) + noise(gen);
                }
                while (y > 1 || y < 0) {
                    y = j * (1 / (double)size) + noise(gen);
                }
            }
            samples.row((size) * i + j) = Eigen::RowVector4d(x, y, 0, 0);
        }
    }
    for (int j = 0; j < size; j++) {
        double x, y;
        x = samples.row(j)(0) + 1;
        y = samples.row(j)(1);
        samples.row(size * size + j) = Eigen::RowVector4d(x, y, 0, 0);
    }
    for (int i = 0; i <= size; i++) {
        double x, y;
        x = samples.row((size) * i)(0);
        y = samples.row((size) * i)(1) + 1;
        samples.row(size * (size + 1) + i) = Eigen::RowVector4d(x, y, 0, 0);
    }
    return samples;
}

Eigen::MatrixXd sphere3gen(int nsamples) {
    std::random_device rd{};
    std::mt19937 gen{ rd() };

    std::uniform_real_distribution<> angle{ 0, M_PI };
    std::uniform_real_distribution<> angle_r{ 0, 2 * M_PI };
    Eigen::MatrixXd sphere_samples(nsamples, 4);

    double radius = 1;

    for (int i = 0; i < nsamples; i++) {
        double p = angle(gen);
        double q = angle(gen);
        double r = angle_r(gen);

        double x = radius * cos(p);
        double y = radius * sin(p) * cos(q);
        double z = radius * sin(p) * sin(q) * cos(r);
        double w = radius * sin(p) * sin(q) * sin(r);

        sphere_samples.row(i) = Eigen::Vector4d(x, y, z, w).transpose();
    }
    return sphere_samples;
}



int main() {
    int dim = 4;

    Eigen::MatrixXd V;
    Eigen::MatrixXi F;
    // V = cliffordgen(1000);
    V = sphere3gen(1000);
    // V = cliffordgengrid(gen_rectangle(1000));
    //read_obj("../output/cliffordgen.obj_20260309_110635.obj", V, F);

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

    polyscope::show();

    std::vector<Eigen::VectorXd> reconstructed_vertices = vertices;
    auto time_start = std::chrono::high_resolution_clock::now();
    flow_complex(dt, reconstructed_vertices, reconstructed_faces, vertex_to_index);
    std::vector<Eigen::VectorXd> V_reconstructed = stereo_projection(reconstructed_vertices);
    polyscope::registerPointCloud("index 2", V_reconstructed);
    auto time_end = std::chrono::high_resolution_clock::now();

    std::chrono::duration<double, std::milli> float_ms = time_end - time_start;

    std::cout << "Elapsed time: " << float_ms.count() << " ms" << std::endl;

    std::cout << "Number of vertices: " << reconstructed_vertices.size() << std::endl;
    std::cout << "Number of faces: " << faces.size() << std::endl;

    if (!reconstructed_faces.empty()) {
        auto *psRecMesh = polyscope::registerSurfaceMesh("reconstructed mesh", V_reconstructed, reconstructed_faces);
    }

    polyscope::show();

    return 0;
}
