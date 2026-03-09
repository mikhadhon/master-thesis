#include <random>

#include "utils.h"
#include "Delaunay.h"
#include "polyscope/point_cloud.h"
#include "polyscope/surface_mesh.h"
#include "polyscope/render/engine.h"
#include "igl/writeOBJ.h"
#include "igl/readPLY.h"
#include "igl/readOBJ.h"

Eigen::MatrixXd cliffordgen(int nsamples) {
    std::random_device rd{};
    std::mt19937 gen{ rd() };
    const double radius = 1 / sqrt(2);

    std::uniform_real_distribution<> angle_t{ 0, 2 * M_PI }, angle_p{ 0, 2 * M_PI };
    std::uniform_real_distribution<> noise{ -0.0001, 0.0001};

    Eigen::MatrixXd samples(nsamples, 4);

    for (int i = 0; i < nsamples; i++) {
        double p = angle_t(gen);
        double q = angle_p(gen);

        double a = 1;
        double b = 1 /*+ noise(gen)*/;

        auto sample = Eigen::Vector4d(a * cos(p), a * sin(p), b * cos(q), b * sin(q));
        sample *= radius;

        samples.row(i) = sample.transpose();
    }
    return samples;
}

Eigen::MatrixXd stereo_projection(Eigen::MatrixXd object) {
    double l = 1;
    double w = 0;
    Eigen::MatrixXd projection_matrix(3, 4);
    projection_matrix << l, 0.0, 0.0, 0.0,
                         0.0, l, 0.0, 0.0,
                         0.0, 0.0, l, 0.0;

    Eigen::MatrixXd m(object.rows(), 3);
    for (int i = 0; i < object.rows(); i++) {
        m.row(i) = ((projection_matrix / (l - object.row(i)(3))) * object.row(i).transpose()).transpose();
    }
    return m;
}

std::vector<Eigen::VectorXd> stereo_projection(std::vector<Eigen::VectorXd> object) {
    double l = 1;
    double w = 0;
    Eigen::MatrixXd projection_matrix(3, 4);
    projection_matrix << l, 0.0, 0.0, 0.0,
                         0.0, l, 0.0, 0.0,
                         0.0, 0.0, l, 0.0;

    std::vector<Eigen::VectorXd> m;
    for (auto &o : object) {
        m.push_back(((projection_matrix / (l - o[3])) * o));
    }
    return m;
}

std::string get_timestamp() {
    auto now = std::chrono::system_clock::now();
    auto in_time_t = std::chrono::system_clock::to_time_t(now);
    std::stringstream ss;
    ss << std::put_time(std::localtime(&in_time_t), "%Y%m%d_%H%M%S");
    return ss.str();
}

void read_ply(const std::string &file_path, Eigen::MatrixXd &V, Eigen::MatrixXi &F) {
    igl::readPLY(file_path, V, F);
}

void read_obj(const std::string &file_path, Eigen::MatrixXd &V, Eigen::MatrixXi &F) {
    igl::readOBJ(file_path, V, F);
}

void write_to_obj(Eigen::MatrixXd V, Eigen::MatrixXi F, std::string identifier) {
    std::string filename = "../output/" + identifier + "_" + get_timestamp() + ".obj";
    if (igl::writeOBJ(filename, V, F)) {
        std::cout << "Curve network saved to: " << filename << std::endl;
    } else {
        std::cerr << "Failed to save curve to: " << filename << std::endl;
    }
}

Eigen::VectorXd make_point_eigen(Point point) {
    int dimension = point.dimension();
    Eigen::VectorXd eigen_point(dimension);
    for (int i = 0; i < dimension; i++) {
        eigen_point(i) = CGAL::to_double(point[i]);
    }
    return eigen_point;
}

void map_vertices_to_vector(
    Delaunay &delaunay,
    std::vector<Eigen::VectorXd> &vertices,
    std::map<Delaunay::Vertex_handle, size_t> &vertex_to_index
) {
    size_t vertex_index = 0;
    int dim = delaunay.maximal_dimension();
    for (auto vertex_it = delaunay.vertices_begin();
        vertex_it != delaunay.vertices_end(); ++vertex_it) {
        if (!delaunay.is_infinite(vertex_it)) {
            Delaunay::Point p = vertex_it->point();
            Eigen::VectorXd v(dim);
            v << CGAL::to_double(p[0]), CGAL::to_double(p[1]), CGAL::to_double(p[2]) , CGAL::to_double(p[3]);
            vertices.push_back(v);
            vertex_to_index[vertex_it] = vertex_index++;
        }
    }
}

void write_faces_to_vector(
    Delaunay &delaunay,
    std::vector<std::array<size_t, 3> > &faces,
    std::map<Delaunay::Vertex_handle, size_t> &vertex_to_index
) {
    auto dt_faces = get_delaunay_faces(delaunay);
    for (auto dt_face = dt_faces.begin(); dt_face != dt_faces.end(); ++dt_face) {

        std::array<size_t, 3> face{};

        for (int i = 0; i < 3; ++i) {
            face[i] = vertex_to_index[dt_face->face.vertex(i)];
        }
        faces.push_back(face);
    }
}



