#include <queue>
#include <polyscope/polyscope.h>
#include <polyscope/surface_mesh.h>
#include <polyscope/curve_network.h>
#include <polyscope/point_cloud.h>

#include "FlowComplex.h"
#include "utils.h"

void flow_complex(Delaunay &delaunay, std::vector<Eigen::VectorXd> &vertices, std::vector<std::array<size_t, 3>> &faces, std::map<Delaunay::Vertex_handle, size_t> &vertex_to_index, std::vector<std::vector<std::array<size_t, 3>>> &index_2_stable_manifolds) {
    std::map<std::array<double, 3>, size_t> fc_vertex_to_index;
    int vertex_count = static_cast<int>(vertices.size());
    std::vector<Eigen::VectorXd> index_two_points;
    std::vector<std::array<size_t, 3>> delaunay_faces;
    std::vector<std::array<size_t, 3>> non_gabriel_faces;
    std::vector<std::array<size_t, 3>> index_two_triangles;
    std::vector<Eigen::VectorXd> circumcenters;

    for (auto facet = delaunay.facets_begin(); facet != delaunay.facets_end(); ++facet) {
        if (!delaunay.is_infinite(*facet)) {
            Delaunay::Face face_no_co_vertex(facet->full_cell());
            int face_vertex_count = 0;
            for (auto vertex = face_no_co_vertex.full_cell()->vertices_begin(); vertex != face_no_co_vertex.full_cell()->vertices_end(); ++vertex) {
                int index_in_cell = face_no_co_vertex.full_cell()->index(*vertex);
                if (index_in_cell != facet->index_of_covertex()) {
                    face_no_co_vertex.set_index(face_vertex_count++, index_in_cell);
                }
            }
            Face face(face_no_co_vertex, facet->index_of_covertex());
            delaunay_faces.push_back({vertex_to_index[face.face.vertex(0)], vertex_to_index[face.face.vertex((1))], vertex_to_index[face.face.vertex((2))]});
            int debug = 5991;
            if (vertex_to_index[face.face.vertex(0)] == debug || vertex_to_index[face.face.vertex(1)] == debug || vertex_to_index[face.face.vertex(2)] == debug) {
                std::vector face_vertices = {face.face.vertex(0)->point(), face.face.vertex(1)->point(), face.face.vertex(2)->point()};
                circumcenters.push_back(make_point_eigen(circumcenter()(face_vertices.begin(), face_vertices.end())));
                std::vector<Eigen::VectorXd> vv;
                auto current_voronoi_edge = delaunay_face_dual(face, delaunay);
                vv.push_back(current_voronoi_edge.vertex1.point);
                vv.push_back(current_voronoi_edge.vertex2.point);
                std::vector<std::array<size_t, 2>> ve = {{0, 1}};
                //polyscope::registerCurveNetwork(std::to_string(count++), vv, ve);
            }
            if (is_index_two_critical_point(face, delaunay)) {
                std::vector<std::array<size_t, 3>> stable_manifold;
                std::queue<std::tuple<Eigen::VectorXd, Edge, Face>> edge_queue;

                std::vector<Point> cgal_points = {face.face.vertex(0)->point(), face.face.vertex(1)->point(), face.face.vertex(2)->point()};
                Point cgal_center = circumcenter()(cgal_points.begin(), cgal_points.end());

                Eigen::VectorXd center = make_point_eigen(cgal_center);
                index_two_points.push_back(center);

                edge_queue.emplace(center, Edge(face.face.vertex(0), face.face.vertex(1), face.face.vertex(2)), face);
                edge_queue.emplace(center, Edge(face.face.vertex(0), face.face.vertex(2), face.face.vertex(1)), face);
                edge_queue.emplace(center, Edge(face.face.vertex(1), face.face.vertex(2), face.face.vertex(0)), face);

                while (!edge_queue.empty()) {
                    std::tuple<Eigen::VectorXd, Edge, Face> current = edge_queue.front();
                    edge_queue.pop();
                    Eigen::VectorXd current_center = std::get<0>(current);
                    Edge current_edge = std::get<1>(current);
                    Face current_delaunay_face = std::get<2>(current);

                    if (is_gabriel(current_edge, delaunay)) {
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
                        stable_manifold.push_back(fc_face);
                    }
                    else {
                        Voronoi_edge current_voronoi_edge = delaunay_face_dual(current_delaunay_face, delaunay);
                        index_two_triangles.push_back({vertex_to_index[current_edge.vertex1], vertex_to_index[current_edge.vertex2], vertex_to_index[current_edge.co_vertex]});

                        Voronoi_face voronoi_face = delaunay_edge_dual(current_edge, current_delaunay_face, delaunay);
                        Point driver = midpoint()(Point(current_edge.vertex1->point()), Point(current_edge.vertex2->point()));
                        std::optional<Eigen::VectorXd> s_prime;

                        std::optional<Face> next_face;
                        for (auto v_edge : voronoi_face.voronoi_edges) {
                            if (v_edge == current_voronoi_edge) {
                                continue;
                            }
                            Eigen::VectorXd intersection;
                            if (intersect_segment_ray(std::make_pair(make_point_eigen(driver), current_center), std::make_pair(v_edge.vertex1.point, v_edge.vertex2.point) , intersection)) {
                                s_prime = intersection;

                                Face dual_delaunay_face = voronoi_edge_dual(v_edge);
                                next_face = dual_delaunay_face;

                                std::array<size_t, 3> fc_face1{}, fc_face2{};
                                std::array fc_center = {current_center[0], current_center[1], current_center[2]};
                                std::array fc_s_prime = {s_prime.value()[0], s_prime.value()[1], s_prime.value()[2]};
                                if (auto fc_index = fc_vertex_to_index.find(fc_center); fc_index == fc_vertex_to_index.end()) {
                                    fc_vertex_to_index[fc_center] = vertex_count++;
                                }
                                if (auto fc_s_prime_index = fc_vertex_to_index.find(fc_s_prime); fc_s_prime_index == fc_vertex_to_index.end()) {
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
                                non_gabriel_faces.push_back(fc_face1);
                                non_gabriel_faces.push_back(fc_face2);
                            }
                        }

                        if (next_face.has_value()) {
                            for (int i = 0; i <= next_face->face.face_dimension(); i++) {
                                Edge next_edge(next_face->face.vertex(i), next_face->face.vertex((i + 1)%3), next_face->face.vertex((i + 2)%3));
                                if (next_edge == current_edge) continue;
                                edge_queue.emplace(s_prime.value(), next_edge, next_face.value());
                            }
                        }
                    }
                }
                index_2_stable_manifolds.push_back(stable_manifold);
            }
        }
    }
    auto indexTwoPoints = polyscope::registerPointCloud("index two points", index_two_points);
    //polyscope::registerSurfaceMesh("delaunay readout", vertices, delaunay_faces);
    vertices.resize(vertex_count);
    for (auto fc_vertex : fc_vertex_to_index) {
        Eigen::VectorXd new_vertex(delaunay.maximal_dimension());
        new_vertex << fc_vertex.first[0], fc_vertex.first[1], fc_vertex.first[2];
        vertices[fc_vertex.second] = new_vertex;
    }
    //polyscope::registerSurfaceMesh("non-gabriel", vertices, non_gabriel_faces);
    polyscope::registerPointCloud("circumcenters", circumcenters);
}