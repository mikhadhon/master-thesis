#include <queue>
#include <Eigen/StdVector>
#include <polyscope/polyscope.h>
#include "polyscope/curve_network.h"

#include "FlowComplex.h"
#include "utils.h"

void flow_complex(Delaunay &delaunay, std::vector<std::array<double, 3> > &vertices, std::vector<std::array<size_t, 3>> &faces, std::map<Delaunay::Vertex_handle, size_t> &vertex_to_index, std::vector<Eigen::Vector3d> &centers, std::vector<Eigen::Vector3d> &centers2) {
    std::map<std::array<double, 3>, size_t> fc_vertex_to_index;
    int vertex_count = vertices.size();
    std::set<Edge> visited;

    for (auto face = delaunay.facets_begin(); face != delaunay.facets_end(); ++face) {
        std::vector<Delaunay::Vertex_handle> facet_vertices;
        get_facet_vertices(delaunay, face, facet_vertices);
        if (facet_vertices.size() == delaunay.maximal_dimension() && is_index_two_critical_point(facet_vertices)) {
            std::queue<std::pair<Eigen::Vector3d, Edge>> edge_queue;

            Eigen::VectorXd center(delaunay.maximal_dimension());
            double radius;
            Eigen::Vector3d i(facet_vertices[0]->point()[0], facet_vertices[0]->point()[1], facet_vertices[0]->point()[2]);
            Eigen::Vector3d j(facet_vertices[1]->point()[0], facet_vertices[1]->point()[1], facet_vertices[1]->point()[2]);
            Eigen::Vector3d l(facet_vertices[2]->point()[0], facet_vertices[2]->point()[1], facet_vertices[2]->point()[2]);

            triangle_circumcircle(i, j, l, center, radius);

            edge_queue.emplace(center, Edge(facet_vertices[0], facet_vertices[1]));
            edge_queue.emplace(center, Edge(facet_vertices[0], facet_vertices[2]));
            edge_queue.emplace(center, Edge(facet_vertices[1], facet_vertices[2]));

            visited.insert(Edge(facet_vertices[0], facet_vertices[1]));
            visited.insert(Edge(facet_vertices[1], facet_vertices[2]));
            visited.insert(Edge(facet_vertices[0], facet_vertices[2]));

            while (!edge_queue.empty()) {
                auto current = edge_queue.front();
                edge_queue.pop();
                Edge current_edge = current.second;
                if (is_gabriel(current_edge)) {
                    std::array fc_center = {center[0], center[1], center[2]};
                    auto fc_index = fc_vertex_to_index.find(fc_center);
                    std::array<size_t, 3> fc_face{};
                    if (fc_index == fc_vertex_to_index.end()) {
                        fc_vertex_to_index[fc_center] = vertex_count++;
                    }
                    fc_face[0] = vertex_to_index[current_edge.vertex1];
                    fc_face[1] = vertex_to_index[current_edge.vertex2];
                    fc_face[2] = fc_vertex_to_index[fc_center];
                    faces.push_back(fc_face);
                    centers2.push_back(center);
                }
                else {
                    std::vector<std::pair<Delaunay::Full_cell_handle, Delaunay::Full_cell_handle>> neighboring_cells;
                    voronoi_facet_from_edge(current_edge, neighboring_cells, delaunay);

                    std::vector<Voronoi_edge> voronoi_edges;

                    for (auto cell_pair : neighboring_cells) {
                        if (!delaunay.is_infinite(cell_pair.first) && !delaunay.is_infinite(cell_pair.second)) {
                            double radius1, radius2;
                            Eigen::VectorXd center1, center2;
                            simplex_circumsphere(cell_pair.first, radius1, center1);
                            simplex_circumsphere(cell_pair.second, radius2, center2);
                            voronoi_edges.push_back(Voronoi_edge(center1, center2, cell_pair.first, cell_pair.second));
                        }
                    }

                    if (!voronoi_edges.empty()) {
                        Eigen::VectorXd driver;
                        Eigen::VectorXd point_on_plane = voronoi_edges[0].vertex1;
                        calculate_driver(point_on_plane, current_edge, driver);
                        Eigen::VectorXd ray_direction = (center - driver).normalized();

                        Eigen::VectorXd origin = (make_point_eigen(current_edge.vertex1->point()) + make_point_eigen(current_edge.vertex2->point())) / 2;

                        double debug_dist = (origin - driver).norm();

                        double closest_intersection_dist = INFINITY;
                        Voronoi_edge *closest_edge = nullptr;
                        Eigen::VectorXd s_prime(delaunay.maximal_dimension());
                        for (auto edge: voronoi_edges) {
                            Eigen::VectorXd intersection(delaunay.maximal_dimension());

                            if (intersect_ray_segment(driver, center, edge.vertex1, edge.vertex2, intersection)) {
                                double current_intersection_dist = (intersection - driver).norm();
                                if (current_intersection_dist < closest_intersection_dist) {
                                    s_prime = intersection;
                                    closest_intersection_dist = current_intersection_dist;
                                    closest_edge = &edge;
                                }
                            }
                        }

                        if (closest_intersection_dist < INFINITY) {
                            std::array<size_t, 3> fc_face1{}, fc_face2{};
                            std::array fc_center = {center[0], center[1], center[2]};
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
                            centers2.push_back(center);
                            centers.push_back(s_prime);

                            std::vector<Delaunay::Vertex_handle> next_face;
                            for (int i = 0; i <= delaunay.maximal_dimension(); i++) {
                                if (closest_edge->cell2->has_vertex(closest_edge->cell1->vertex(i))) {
                                    next_face.push_back(closest_edge->cell1->vertex(i));
                                }
                            }
                            Eigen::VectorXd next_center(delaunay.maximal_dimension());
                            double next_radius;
                            Eigen::Vector3d next_i(next_face[0]->point()[0], next_face[0]->point()[1], next_face[0]->point()[2]);
                            Eigen::Vector3d next_j(next_face[1]->point()[0], next_face[1]->point()[1], next_face[1]->point()[2]);
                            Eigen::Vector3d next_l(next_face[2]->point()[0], next_face[2]->point()[1], next_face[2]->point()[2]);

                            triangle_circumcircle(next_i, next_j, next_l, next_center, next_radius);

                            for (int i = 0; i < delaunay.maximal_dimension(); i++) {
                                for (int j = i + 1; j < delaunay.maximal_dimension(); j++) {
                                    Edge new_edge(next_face[i], next_face[j]);
                                    if (!(new_edge == current_edge) && !visited.contains(new_edge)) {
                                        edge_queue.emplace(std::make_pair(next_center, new_edge));
                                        //visited.insert(new_edge);
                                    }
                                }
                            }
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
}
