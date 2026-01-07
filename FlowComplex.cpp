#include <queue>
#include <Eigen/StdVector>
#include <polyscope/polyscope.h>
#include "polyscope/curve_network.h"
#include <polyscope/surface_mesh.h>
#include "polyscope/point_cloud.h"


#include "FlowComplex.h"
#include "utils.h"

void flow_complex(Delaunay &delaunay, std::vector<std::array<double, 3> > &vertices, std::vector<std::array<size_t, 3>> &faces, std::map<Delaunay::Vertex_handle, size_t> &vertex_to_index) {
    std::map<std::array<double, 3>, size_t> fc_vertex_to_index;
    int vertex_count = vertices.size();
    std::vector<Point> s_primes;
    for (auto facet = delaunay.facets_begin(); facet != delaunay.facets_end(); ++facet) {
            Delaunay::Face face_no_co_vertex((*facet).full_cell());
            int face_vertex_count = 0;
            for (auto vertex = face_no_co_vertex.full_cell()->vertices_begin(); vertex != face_no_co_vertex.full_cell()->vertices_end(); ++vertex) {
                int index_in_cell = face_no_co_vertex.full_cell()->index(*vertex);
                if (index_in_cell != facet->index_of_covertex()) {
                    face_no_co_vertex.set_index(face_vertex_count++, index_in_cell);
                }
            }
            Face face(face_no_co_vertex, facet->index_of_covertex());

            if (is_index_two_critical_point(face, delaunay)) {
                std::queue<std::tuple<Eigen::VectorXd, Edge, Face>> edge_queue;

                std::vector<Point> cgal_points = {face.face.vertex(0)->point(), face.face.vertex(1)->point(), face.face.vertex(2)->point()};
                Point cgal_center = circumcenter()(cgal_points.begin(), cgal_points.end());

                Eigen::VectorXd center = make_point_eigen(cgal_center);

                edge_queue.emplace(center, Edge(face.face.vertex(0), face.face.vertex(1), face.face.vertex(2)), face);
                edge_queue.emplace(center, Edge(face.face.vertex(0), face.face.vertex(2), face.face.vertex(1)), face);
                edge_queue.emplace(center, Edge(face.face.vertex(1), face.face.vertex(2), face.face.vertex(0)), face);

                while (!edge_queue.empty()) {
                    std::tuple<Eigen::VectorXd, Edge, Face> current = edge_queue.front();
                    edge_queue.pop();
                    Eigen::VectorXd current_center = std::get<0>(current);
                    Edge current_edge = std::get<1>(current);
                    Face current_delaunay_face = std::get<2>(current);

                    if (is_gabriel(current_edge)) {
                        // std::vector<Point> a = {current_edge.vertex1->point(), current_edge.vertex2->point(), current_edge.co_vertex->point()};
                        // a.push_back(Point(current.first[0], current.first[1], current.first[2]));
                        // std::vector<std::array<size_t, 3>> b = {{3, 0, 1}};
                        // polyscope::registerSurfaceMesh("index 1 triangle", a, b);
                        // polyscope::registerPointCloud("triangle points", a);
                        // std::vector<Eigen::VectorXd> vv;
                        // vv.push_back(edge.vertex1.point);
                        // vv.push_back(edge.vertex2.point);
                        // std::vector<std::array<size_t, 2>> ve = {{0, 1}};
                        // polyscope::registerCurveNetwork("ve", vv, ve);
                        // polyscope::show();

                        std::array fc_center = {current_center[0], current_center[1], current_center[2]};
                        auto fc_index = fc_vertex_to_index.find(fc_center);
                        std::array<size_t, 3> fc_face{};
                        if (fc_index == fc_vertex_to_index.end()) {
                            fc_vertex_to_index[fc_center] = vertex_count++;
                        }
                        fc_face[0] = vertex_to_index[current_edge.vertex1];
                        fc_face[1] = vertex_to_index[current_edge.vertex2];
                        fc_face[2] = fc_vertex_to_index[fc_center];
                        faces.push_back(fc_face);
                    }
                    else {
                        Voronoi_edge current_voronoi_edge = delaunay_face_dual(current_delaunay_face, delaunay);
                        Voronoi_face voronoi_face = delaunay_edge_dual(current_edge, current_delaunay_face, delaunay);
                        Eigen::Vector3d face_normal = make_point_eigen(current_edge.vertex2->point()) - make_point_eigen(current_edge.vertex1->point());
                        Point driver = midpoint()(current_edge.vertex1->point(), current_edge.vertex2->point());
                        Point s_prime;

                        std::optional<Face> next_face;

                        for (auto v_edge : voronoi_face.voronoi_edges) {
                            Eigen::Vector3d edge_dir = v_edge.vertex2.point - v_edge.vertex1.point;
                            Eigen::VectorXd edge_normal = edge_dir.cross(face_normal);

                            double side_d = edge_normal.dot(make_point_eigen(driver) - v_edge.vertex1.point);
                            double side_v = edge_normal.dot(current_center - v_edge.vertex1.point);

                            if (side_d * side_v < 0) {
                                Eigen::VectorXd ray_direction = (current_center - make_point_eigen(driver)).normalized();
                                Eigen::VectorXd w0 = make_point_eigen(driver) - v_edge.vertex1.point;
                                double a = ray_direction.dot(ray_direction);      // always 1 if normalized
                                double b = ray_direction.dot(edge_dir);
                                double c = edge_dir.dot(edge_dir);
                                double d_coef = ray_direction.dot(w0);
                                double e = edge_dir.dot(w0);

                                double denom = a * c - b * b;
                                double t = (b * e - c * d_coef) / denom;

                                Eigen::VectorXd temp = make_point_eigen(driver) + t * ray_direction;
                                s_prime = Point(temp(0), temp(1), temp(2));
                                break;
                            }
                            // if (v_edge == current_voronoi_edge) continue;
                            // Face dual_delaunay_face = voronoi_edge_dual(v_edge);
                            // std::vector dual_delaunay_face_vertices = {dual_delaunay_face.face.vertex(0)->point(), dual_delaunay_face.face.vertex(1)->point(), dual_delaunay_face.face.vertex(2)->point()};
                            // std::vector<Point> current_delaunay_face_vertices = {current_delaunay_face.face.vertex(0)->point(), current_delaunay_face.face.vertex(1)->point(), current_delaunay_face.face.vertex(2)->point()};
                            // Point dual_circumcenter = circumcenter()(dual_delaunay_face_vertices.begin(), dual_delaunay_face_vertices.end());
                            // if (is_point_in_triangle(current_delaunay_face_vertices[0], current_delaunay_face_vertices[1], current_delaunay_face_vertices[2], dual_circumcenter)) {
                            //     s_prime = dual_circumcenter;
                            //     s_primes.push_back(s_prime);
                            //     next_face = dual_delaunay_face;
                            //     break;
                            // }
                        }
                        std::array<size_t, 3> fc_face1{}, fc_face2{};
                        std::array fc_center = {current_center[0], current_center[1], current_center[2]};
                        std::array fc_s_prime = {s_prime[0], s_prime[1], s_prime[2]};
                        auto fc_index = fc_vertex_to_index.find(fc_center);
                        if (fc_index == fc_vertex_to_index.end()) {
                            fc_vertex_to_index[fc_center] = vertex_count++;
                        }
                        auto fc_s_prime_index = fc_vertex_to_index.find(fc_s_prime);
                        if (fc_s_prime_index == fc_vertex_to_index.end()) {
                            fc_vertex_to_index[fc_s_prime] = vertex_count++;
                        }
                        fc_face1[0] = vertex_to_index[current_edge.vertex1];
                        fc_face1[1] = fc_vertex_to_index[fc_s_prime];
                        fc_face1[2] = fc_vertex_to_index[fc_center];
                        fc_face2[0] = vertex_to_index[current_edge.vertex2];
                        fc_face2[1] = fc_vertex_to_index[fc_s_prime];
                        fc_face2[2] = fc_vertex_to_index[fc_center];
                        faces.push_back(fc_face1);
                        faces.push_back(fc_face2);

                        if (next_face.has_value()) {
                            for (int i = 0; i <= next_face->face.face_dimension(); i++) {
                                Edge next_edge(next_face->face.vertex(i), next_face->face.vertex((i + 1)%3), next_face->face.vertex((i + 2)%3));
                                if (next_edge == current_edge) continue;
                                edge_queue.emplace(make_point_eigen(s_prime), next_edge, next_face.value());
                            }
                        }
                    }
                }
            }
        }
    vertices.resize(vertex_count);
    for (auto fc_vertex : fc_vertex_to_index) {
        vertices[fc_vertex.second] = fc_vertex.first;
    }
    polyscope::registerPointCloud("s primes", s_primes);
}