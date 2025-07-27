#include <random>

#include "utils.h"
#include "Delaunay.h"

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

    for (int i = 0; i < count; i++) {
        double t = angle_t(gen);
        double p = angle_p(gen);

        double x = (rot_radius + radius * cos(p)) * cos(t);
        double y = (rot_radius + radius * cos(p)) * sin(t);
        double z = radius * sin(p);

        torus_samples.push_back(Delaunay::Point(x, y, z));
    }
}
