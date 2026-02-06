#include <random>

#include "utils.h"
#include "Delaunay.h"
#include "polyscope/point_cloud.h"
#include "polyscope/surface_mesh.h"
#include "polyscope/render/engine.h"
#include "igl/writeOBJ.h"
#include "igl/readPLY.h"
#include "igl/readOBJ.h"

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

std::vector<Delaunay::Point> get_points_from_handles(const std::vector<Delaunay::Vertex_handle> &handles) {
    std::vector<Delaunay::Point> points;
    points.reserve(handles.size());
    for (const auto& handle : handles) {
        points.push_back(handle->point());
    }
    return points;
}

void gen_sphere_sample(int count, double radius, std::vector<Delaunay::Point> &points) {
    std::random_device rd{};
    std::mt19937 gen{ rd() };

    std::uniform_real_distribution<> angle_t{ 0, 2 * M_PI }, angle_p{ 0, 2 * M_PI };
    std::uniform_real_distribution<> v_dist{-1.0, 1.0};
    std::uniform_real_distribution<> noise{ 0, 0.01 };
    Eigen::MatrixXd sphere_samples(count, 3);

    for (int i = 0; i < count; i++) {
        double p = angle_t(gen);
        double q = v_dist(gen);

        double horizontal_radius = sqrt(1 - q * q);
        double r = radius + noise(gen);

        double x = r * horizontal_radius * cos(p);
        double y = r * horizontal_radius * sin(p);
        double z = r * q;

        sphere_samples.row(i) = Eigen::Vector3d(x, y, z).transpose();
    }
    for (int i = 0; i < sphere_samples.rows(); i++) {
        points.push_back(Delaunay::Point(sphere_samples.row(i)[0], sphere_samples.row(i)[1], sphere_samples.row(i)[2] ));
    }
}

void gen_rectangle(int n, std::vector<Delaunay::Point> &points) {
    std::random_device rd{};
    std::mt19937 gen{ rd() };

    std::uniform_real_distribution<> noise{ -0.0000001, 0.00000001 };

    int size = n;

    Eigen::MatrixXd samples((size + 1) * (size + 1), 3);

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
            samples.row((size)*i + j) = Eigen::RowVector3d(x, y, 0);
        }
    }
    for (int j = 0; j < size; j++) {
        double x, y;
        x = samples.row(j)(0) + 1;
        y = samples.row(j)(1);
        samples.row(size * size + j) = Eigen::RowVector3d(x, y, 0);
    }
    for (int i = 0; i <= size; i++) {
        double x, y;
        x = samples.row((size)*i)(0);
        y = samples.row((size)*i)(1) + 1;
        samples.row(size * (size + 1) + i) = Eigen::RowVector3d(x, y, 0);
    }

    for (int i = 0; i < samples.rows(); i++) {
        points.push_back(Delaunay::Point(samples.row(i)(0), samples.row(i)(1), samples.row(i)(2)));
    }
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
            v << CGAL::to_double(p[0]), CGAL::to_double(p[1]), CGAL::to_double(p[2]);
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
    for (auto facet = delaunay.facets_begin(); facet != delaunay.facets_end(); ++facet) {

            std::array<size_t, 3> face{};
            size_t face_vertex_idx = 0;
            std::vector<Delaunay::Vertex_handle> face_vertices;

            get_facet_vertices(delaunay, facet, face_vertices);
            if (face_vertices.size() == delaunay.maximal_dimension()) {
                for (int i = 0; i < delaunay.maximal_dimension(); ++i) {
                    face[face_vertex_idx++] = vertex_to_index[face_vertices[i]];
                }
                faces.push_back(face);
            }

    }
}

void generate_torus(int count, double radius, double rot_radius, std::vector<Delaunay::Point> &torus_samples) {
    std::random_device rd{};
    std::mt19937 gen{ rd() };

    //std::normal_distribution<> angle_t{M_PI, M_PI_2}, angle_p{M_PI, M_PI_2};
    std::uniform_real_distribution<> angle_t{ 0, 2 * M_PI }, angle_p{ 0, 2 * M_PI };
    std::uniform_real_distribution<> noise{ -0.0000001, 0.00000001 };

    for (int i = 0; i < count; i++) {
        double t = angle_t(gen);
        double p = angle_p(gen);

        double x = (rot_radius + radius * cos(p)) * cos(t) + noise(gen);
        double y = (rot_radius + radius * cos(p)) * sin(t) + noise(gen);
        double z = radius * sin(p) + noise(gen);

        torus_samples.push_back(Delaunay::Point(x, y, z));
    }
}

// bool is_point_in_triangle(const Point& p0, const Point& p1, const Point& p2, const Point& query) {
//     Eigen::Vector3d v0(p0[0], p0[1], p0[2]);
//     Eigen::Vector3d v1(p1[0], p1[1], p1[2]);
//     Eigen::Vector3d v2(p2[0], p2[1], p2[2]);
//     Eigen::Vector3d p(query[0], query[1], query[2]);
//
//     Eigen::Vector3d v0v1 = v1 - v0;
//     Eigen::Vector3d v0v2 = v2 - v0;
//     Eigen::Vector3d v0p = p - v0;
//
//     double dot00 = v0v1.dot(v0v1);
//     double dot01 = v0v1.dot(v0v2);
//     double dot02 = v0v1.dot(v0p);
//     double dot11 = v0v2.dot(v0v2);
//     double dot12 = v0v2.dot(v0p);
//
//     double invDenom = 1.0 / (dot00 * dot11 - dot01 * dot01);
//     double u = (dot11 * dot02 - dot01 * dot12) * invDenom;
//     double v = (dot00 * dot12 - dot01 * dot02) * invDenom;
//
//     const double eps = 1e-10;
//     return (u >= -eps) && (v >= -eps) && (u + v <= 1.0 + eps);
// }

void generate_trefoil(int nsamples, int normal_windings, Eigen::Matrix<double, Eigen::Dynamic, 3> &V, Eigen::Matrix<double, Eigen::Dynamic, 3>  &N1, Eigen::Matrix<double, Eigen::Dynamic, 3> &N2, Eigen::Matrix<double, Eigen::Dynamic, 3> &T){
    double r = 1.;

    V.resize(nsamples, 3);
    N1.resize(nsamples, 3);
    N2.resize(nsamples, 3);
    T.resize(nsamples, 3);

    // double stretchz = 2.7;
    double stretchz = 1.;

    for (int i=0; i<nsamples; i++) {
        double t = (1.*i) / nsamples;
        double dtheta = (2.*M_PI);
        double theta  = dtheta*t;

        V.row(i)  = Eigen::Vector3d(    cos(theta) + 2*cos(2*theta),
                                        sin(theta) - 2*sin(2*theta),
                                       -sin(3*theta) * stretchz);
        Eigen::Vector3d Tangent = Eigen::Vector3d( -sin(theta) - 4*sin(2*theta),
                                                    cos(theta) - 4*cos(2*theta),
                                                 -3*cos(3*theta) * stretchz ).normalized();
        T.row(i)  = Tangent;
        Eigen::Vector3d dTudt = Eigen::Vector3d(  -cos(theta) - 8*cos(2*theta),
                                                  -sin(theta) + 8*sin(2*theta),
                                                 9*sin(3*theta) * stretchz);
        N1.row(i) = (dTudt - dTudt.dot(Tangent) * Tangent).normalized();

        N2.row(i) = T.row(i).cross(N1.row(i));
        N2.row(i).normalize();

        if (normal_windings > 0) {
            theta = (2.*M_PI/nsamples)*i * normal_windings;
            Eigen::Matrix3d frame;
            frame.row(0) = N1.row(i);
            frame.row(1) = N2.row(i);
            frame.row(2) =  T.row(i);
            Eigen::Matrix3d Rz; Rz << cos(theta), -sin(theta), 0, sin(theta), cos(theta), 0, 0, 0, 1;
            Eigen::Matrix3d R =  frame.transpose() * Rz * frame;

            N1.row(i) = N1.row(i) * R.transpose();
            N2.row(i) = N2.row(i) * R.transpose();
        }
    }

    std::cout << "trefoil, min/max per dim: " << std::endl;;
    std::cout << V.col(0).minCoeff() << "-" << V.col(0).maxCoeff() << std::endl;
    std::cout << V.col(1).minCoeff() << "-" << V.col(1).maxCoeff() << std::endl;
    std::cout << V.col(2).minCoeff() << "-" << V.col(2).maxCoeff() << std::endl;

}

bool intersect_segment_ray(std::pair<Eigen::VectorXd, Eigen::VectorXd> ray, std::pair<Eigen::VectorXd, Eigen::VectorXd> segment, Eigen::VectorXd &intersection) {
    const Eigen::Vector3d segment_dir = segment.second - segment.first;

    const Eigen::VectorXd ray_direction = (ray.second - ray.first).normalized();
    const Eigen::VectorXd w0 = ray.first - segment.first;
    const double a = ray_direction.dot(ray_direction);
    const double b = ray_direction.dot(segment_dir);
    const double c = segment_dir.dot(segment_dir);
    const double d_coef = ray_direction.dot(w0);
    const double e = segment_dir.dot(w0);

    const double denom = a * c - b * b;

    const double t_ray = (b * e - c * d_coef) / denom;

    if (const double t_edge = (a * e - b * d_coef) / denom; t_ray >= 0.0 && t_edge >= 0.0 && t_edge <= 1.0) {
        intersection = ray.first + t_ray * ray_direction;
        return true;
    }
    return false;
}
