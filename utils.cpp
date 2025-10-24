#include <random>

#include "utils.h"
#include "Delaunay.h"

Point make_point(Eigen::Vector3d &eigen_point) {
    double coords[3] = {eigen_point[0], eigen_point[1], eigen_point[2]};
    return Point(3, coords, coords+3);
}

Eigen::VectorXd make_point_eigen(Point point) {
    int dimension = point.dimension();
    Eigen::VectorXd eigen_point(dimension);
    for (int i = 0; i < dimension; i++) {
        eigen_point(i) = point[i];
    }
    return eigen_point;
}

void gen_sphere_sample(int count, double radius, std::vector<Delaunay::Point> &points) {
    std::random_device rd{};
    std::mt19937 gen{ rd() };

    std::uniform_real_distribution<> angle_t{ 0, 2 * M_PI }, angle_p{ 0, 2 * M_PI };
    std::uniform_real_distribution<> noise{ -0.01, 0.01 };
    Eigen::MatrixXd sphere_samples(count, 3);

    for (int i = 0; i < count; i++) {
        double p = angle_t(gen);
        double q = angle_p(gen);

        double x = radius * sin(p) * cos(q) + noise(gen);
        double y = radius * sin(p) * sin(q) + noise(gen);
        double z = radius * cos(p) + noise(gen);

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
    std::vector<std::array<double, 3> > &vertices,
    std::map<Delaunay::Vertex_handle, size_t> &vertex_to_index
) {
    size_t vertex_index = 0;
    for (auto vertex_it = delaunay.vertices_begin();
        vertex_it != delaunay.vertices_end(); ++vertex_it) {
        if (!delaunay.is_infinite(vertex_it)) {
            Delaunay::Point p = vertex_it->point();
            vertices.push_back({p[0], p[1], p[2]});
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
        if (!delaunay.is_infinite(*facet)) {

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
}

void generate_circle(int nsamples, int normal_windings, Eigen::Matrix<double, Eigen::Dynamic, 3> &V, Eigen::Matrix<double, Eigen::Dynamic, 3> &N1, Eigen::Matrix<double, Eigen::Dynamic, 3>  &N2, Eigen::Matrix<double, Eigen::Dynamic, 3> &T){
    double r = 1.;
    V.resize(nsamples, 3);
    N1.resize(nsamples, 3);
    N2.resize(nsamples, 3);
    T.resize(nsamples, 3);
    for (int i=0; i<nsamples; i++) {
        double theta = (2.*M_PI/nsamples)*i;
        V.row(i)  = Eigen::Vector3d(cos(theta), sin(theta), 1.);
        N1.row(i) = Eigen::Vector3d(cos(theta), sin(theta), 0.).normalized();
        N2.row(i) = Eigen::Vector3d(        0.,         0., 1.).normalized();
        T.row(i)  = N1.row(i).cross(N2.row(i));
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

//circumcircle calculating without cross-product
void triangle_circumcircle(Eigen::Vector3d &i, Eigen::Vector3d &j, Eigen::Vector3d &l, Eigen::VectorXd &center, double &radius) {
    const Eigen::Vector3d a = i - l;
    const Eigen::Vector3d b = j - l;

    const double cos_theta = (a.dot(b)) / (a.norm() * b.norm());
    radius = ((i - j).norm()) / (2 * sin(acos(cos_theta)));

    const Eigen::Vector3d helper_term = (pow(a.norm(), 2) * b - pow(b.norm(), 2) * a);
    center = ((helper_term.dot(b) * a - helper_term.dot(a) * b) / (2 * (pow(a.norm(), 2) * pow(b.norm(), 2) - pow(a.dot(b), 2)))) + l;
}

void simplex_circumsphere(Delaunay::Full_cell_handle simplex, double &radius, Eigen::VectorXd &center) {
    int dimension = simplex->maximal_dimension();

    Eigen::MatrixXd A(dimension, dimension);
    Eigen::VectorXd b(dimension);

    Eigen::VectorXd p0(dimension);
    for (int i = 0; i < dimension; i++) p0(i) = simplex->vertex(0)->point()[i];
    double p0_squared_norm = p0.squaredNorm();

    for (int i = 1; i <= dimension; i++) {
        Eigen::VectorXd pi(dimension);
        for (int j = 0; j < dimension; j++) pi(j) = simplex->vertex(i)->point()[j];
        A.row(i - 1) = 2 * (pi - p0).transpose();
        b(i - 1) = pi.squaredNorm() - p0_squared_norm;
    }

    center = A.colPivHouseholderQr().solve(b);
    radius = (center - p0).norm();

    for (auto vertex = simplex->vertices_begin(); vertex != simplex->vertices_end(); vertex++) {
        double eps = 1e-10;
        Eigen::Vector3d vi((*vertex)->point()[0], (*vertex)->point()[1], (*vertex)->point()[2]);
        assert(std::abs((vi - center).norm() - radius ) < eps);
    }
}

bool intersect_ray_segment(
    Eigen::VectorXd &ray_origin,
    Eigen::VectorXd &ray_pass_through,
    Eigen::VectorXd &segment_start,
    Eigen::VectorXd &segment_end,
    Eigen::VectorXd &intersection
    ) {
    int dimension = ray_origin.size();
    Eigen::VectorXd ray_direction = ray_pass_through - ray_origin;
    Eigen::VectorXd segment_direction = segment_end - segment_start;
    Eigen::VectorXd right_hand_side = segment_start - ray_origin;

    Eigen::MatrixXd A(dimension, 2);
    A.col(0) = ray_direction;
    A.col(1) = -segment_direction;

    Eigen::Vector2d st = A.colPivHouseholderQr().solve(right_hand_side);

    double s = st(0);
    double t = st(1);

    Eigen::VectorXd point_on_ray = ray_origin + s * ray_direction;
    Eigen::VectorXd point_on_segment = segment_start + t * segment_direction;

    double eps = 1e-10;

    if ((point_on_ray - point_on_segment).norm() > eps) {
        return false;
    }

    if (s < -eps || t < -eps || t > 1 + eps ) {
        return false;
    }

    intersection = point_on_ray;
    return true;
}
