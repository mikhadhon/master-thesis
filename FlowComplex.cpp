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
    std::set<Edge> visited;

    for (auto face = delaunay.facets_begin(); face != delaunay.facets_end(); ++face) {
            bool test = delaunay.is_infinite(*face);
            std::vector<Delaunay::Vertex_handle> facet_vertices;
            get_facet_vertices(delaunay, face, facet_vertices);
            if (facet_vertices.size() == delaunay.maximal_dimension() && is_index_two_critical_point(face, delaunay)) {
                std::queue<std::pair<Eigen::VectorXd, Edge>> edge_queue;

                std::vector<Point> cgal_points = get_points_from_handles(facet_vertices);
                Point cgal_center = circumcenter()(cgal_points.begin(), cgal_points.end());

                Eigen::VectorXd center = make_point_eigen(cgal_center);

                edge_queue.emplace(center, Edge(facet_vertices[0], facet_vertices[1], facet_vertices[2]));
                edge_queue.emplace(center, Edge(facet_vertices[0], facet_vertices[2], facet_vertices[1]));
                edge_queue.emplace(center, Edge(facet_vertices[1], facet_vertices[2], facet_vertices[0]));

                visited.insert(Edge(facet_vertices[0], facet_vertices[1], facet_vertices[2]));
                visited.insert(Edge(facet_vertices[1], facet_vertices[2], facet_vertices[0]));
                visited.insert(Edge(facet_vertices[0], facet_vertices[2], facet_vertices[1]));

                while (!edge_queue.empty()) {
                    auto current = edge_queue.front();
                    edge_queue.pop();
                    Edge current_edge = current.second;
                    if (is_gabriel(current_edge)) {
                        // std::vector<Point> a = {current_edge.vertex1->point(), current_edge.vertex2->point(), current_edge.co_vertex->point()};
                        // a.push_back(Point(current.first[0], current.first[1], current.first[2]));
                        // std::vector<std::array<size_t, 3>> b = {{3, 0, 1}};
                        // polyscope::registerSurfaceMesh("index 1 triangle", a, b);
                        // polyscope::registerPointCloud("triangle points", a);
                        Voronoi_edge edge = delaunay_face_dual(face, delaunay);
                        // std::vector<Eigen::VectorXd> vv;
                        // vv.push_back(edge.vertex1.point);
                        // vv.push_back(edge.vertex2.point);
                        // std::vector<std::array<size_t, 2>> ve = {{0, 1}};
                        // polyscope::registerCurveNetwork("ve", vv, ve);
                        // polyscope::show();

                        std::array fc_center = {current.first[0], current.first[1], current.first[2]};
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
                        Voronoi_face voronoi_face = delaunay_edge_dual(current_edge, face, delaunay);
                        Point driver = midpoint()(current_edge.vertex1->point(), current_edge.vertex2->point());

                        for (auto v_edge : voronoi_face.voronoi_edges) {
                            Delaunay::Face dual_delaunay_face = voronoi_edge_dual(v_edge);
                            std::vector dual_delaunay_face_vertices = {dual_delaunay_face.vertex(0)->point(), dual_delaunay_face.vertex(1)->point(), dual_delaunay_face.vertex(2)->point()};
                            Point dual_circumcenter = circumcenter()(dual_delaunay_face_vertices.begin(), dual_delaunay_face_vertices.end());


                        }
                    }
                    // else {
                    //     Voronoi_face voronoi_face;
                    //     voronoi_facet_from_edge(current_edge, voronoi_face, delaunay);
                    //
                    //     if (!voronoi_face.voronoi_edges.empty()) {
                    //         Eigen::VectorXd origin = (make_point_eigen(current_edge.vertex1->point()) + make_point_eigen(current_edge.vertex2->point())) / 2;
                    //
                    //         Eigen::VectorXd s_prime(delaunay.maximal_dimension());
                    //         Voronoi_edge *nextEdge = nullptr;
                    //
                    //         for (auto voronoi_edge : voronoi_face.voronoi_edges) {
                    //             Eigen::VectorXd intersection;
                    //             Eigen::VectorXd curr_edge_1 = make_point_eigen(current_edge.vertex1->point());
                    //             Eigen::VectorXd curr_edge_2 = make_point_eigen(current_edge.vertex2->point());
                    //             bool intersection_exists = intersect_voronoi_edge_triangle(voronoi_edge, curr_edge_1, curr_edge_2, current.first, intersection);
                    //
                    //             if (intersection_exists) {
                    //                 if ((intersection - current.first).norm() > 10e-8 || (intersection - current.first).norm() < -10e-8) {
                    //                     s_prime = intersection;
                    //                     nextEdge = &voronoi_edge;
                    //                     break;
                    //                 }
                    //             }
                    //         }
                    //         std::vector<std::array<double, 3>> vv;
                    //         std::vector<std::array<size_t, 2>> ve;
                    //         size_t step = 0;
                    //         for (auto voronoi_edge : voronoi_face.voronoi_edges) {
                    //             if (!voronoi_edge.vertex1.is_infinite && !voronoi_edge.vertex2.is_infinite) {
                    //                 vv.push_back({voronoi_edge.vertex1.point(0), voronoi_edge.vertex1.point(1), voronoi_edge.vertex1.point(2)});
                    //                 vv.push_back({voronoi_edge.vertex2.point(0), voronoi_edge.vertex2.point(1), voronoi_edge.vertex2.point(2)});
                    //             }
                    //             else {
                    //                 Voronoi_vertex* infinite_vertex = voronoi_edge.vertex1.is_infinite ? &voronoi_edge.vertex1 : &voronoi_edge.vertex2;
                    //                 Voronoi_vertex* finite_vertex = voronoi_edge.vertex1.is_infinite ? &voronoi_edge.vertex2 : &voronoi_edge.vertex1;
                    //                 vv.push_back({finite_vertex->point(0), finite_vertex->point(1), finite_vertex->point(2)});
                    //                 Eigen::VectorXd helper_vertex(3);
                    //                 helper_vertex << infinite_vertex->point(0), infinite_vertex->point(1), infinite_vertex->point(2);
                    //                 helper_vertex = helper_vertex + voronoi_edge.voronoi_edge_direction;
                    //                 vv.push_back({helper_vertex(0), helper_vertex(1), helper_vertex(2)});
                    //             }
                    //             ve.push_back({step, step + 1});
                    //             step += 2;
                    //         }
                    //         Eigen::VectorXd curr_del_1 = make_point_eigen(facet_vertices[0]->point());
                    //         Eigen::VectorXd curr_del_2 = make_point_eigen(facet_vertices[1]->point());
                    //         Eigen::VectorXd curr_del_3 = make_point_eigen(facet_vertices[2]->point());
                    //
                    //         std::vector<std::array<double, 3>> current_delaunay_triangle = {{curr_del_1[0], curr_del_1[1], curr_del_1[2]}, {curr_del_2[0], curr_del_2[1], curr_del_2[2]}, {curr_del_3[0], curr_del_3[1], curr_del_3[2]}};
                    //         std::vector<std::array<size_t, 3>> current_delaunay_triangle_face = {{0, 1, 2}};
                    //         std::vector<std::array<double, 3>> current_s = {{current.first(0), current.first(1), current.first(2)}};
                    //         std::vector<std::array<double, 3>> current_s_prime = {{s_prime(0), s_prime(1), s_prime(2)}};
                    //         std::vector<std::array<double, 3>> current_driver = {{origin(0), origin(1), origin(2)}};
                    //         auto *psMesh = polyscope::registerSurfaceMesh("current delaunay triangle", current_delaunay_triangle, current_delaunay_triangle_face);
                    //         auto *pcCloud = polyscope::registerPointCloud("s", current_s);
                    //         // auto *pcCloud2 = polyscope::registerPointCloud("s'", current_s_prime);
                    //         // //auto *pcCloud3 = polyscope::registerPointCloud("driver", current_driver);
                    //         // auto *psCurve = polyscope::registerCurveNetwork("voronoi edges", vv, ve);
                    //         polyscope::show();
                    //
                    //         if (nextEdge != nullptr) {
                    //             std::array<size_t, 3> fc_face1{}, fc_face2{};
                    //             std::array fc_center = {current.first[0], current.first[1], current.first[2]};
                    //             std::array fc_s_prime = {s_prime[0], s_prime[1], s_prime[2]};
                    //             auto fc_index = fc_vertex_to_index.find(fc_center);
                    //             if (fc_index == fc_vertex_to_index.end()) {
                    //                 fc_vertex_to_index[fc_center] = vertex_count++;
                    //             }
                    //             auto fc_s_prime_index = fc_vertex_to_index.find(fc_s_prime);
                    //             if (fc_s_prime_index == fc_vertex_to_index.end()) {
                    //                 fc_vertex_to_index[fc_s_prime] = vertex_count++;
                    //             }
                    //             fc_face1[0] = vertex_to_index[current_edge.vertex1];
                    //             fc_face1[1] = fc_vertex_to_index[fc_s_prime];
                    //             fc_face1[2] = fc_vertex_to_index[fc_center];
                    //             fc_face2[0] = vertex_to_index[current_edge.vertex2];
                    //             fc_face2[1] = fc_vertex_to_index[fc_s_prime];
                    //             fc_face2[2] = fc_vertex_to_index[fc_center];
                    //             faces.push_back(fc_face1);
                    //             faces.push_back(fc_face2);
                    //
                    //             std::vector<Delaunay::Vertex_handle> next_face;
                    //             for (int i = 0; i <= delaunay.maximal_dimension(); i++) {
                    //                 if (nextEdge->cell2->has_vertex(nextEdge->cell1->vertex(i))) {
                    //                     next_face.push_back(nextEdge->cell1->vertex(i));
                    //                 }
                    //             }
                    //             Eigen::VectorXd next_center(delaunay.maximal_dimension());
                    //             double next_radius;
                    //             Eigen::VectorXd next_i(delaunay.maximal_dimension()), next_j(delaunay.maximal_dimension()), next_l(delaunay.maximal_dimension());
                    //             next_i << next_face[0]->point()[0], next_face[0]->point()[1], next_face[0]->point()[2];
                    //             next_j << next_face[1]->point()[0], next_face[1]->point()[1], next_face[1]->point()[2];
                    //             next_l << next_face[2]->point()[0], next_face[2]->point()[1], next_face[2]->point()[2];
                    //
                    //             triangle_circumcircle(next_i, next_j, next_l, next_center, next_radius);
                    //
                    //             for (int i = 0; i < delaunay.maximal_dimension(); i++) {
                    //                 for (int j = i + 1; j < delaunay.maximal_dimension(); j++) {
                    //                     Edge new_edge(next_face[i], next_face[j]);
                    //                     if (!(new_edge == current_edge) && !visited.contains(new_edge)) {
                    //                         //edge_queue.emplace(std::make_pair(next_center, new_edge));
                    //                         visited.insert(new_edge);
                    //                     }
                    //                 }
                    //             }
                    //         }
                    //
                    //     }
                    // }
                }
            }
            // else {
            //     std::vector<Voronoi_face> voronoi_faces;
            //     if (!delaunay.is_infinite(*face)) {
            //         polyscope::registerPointCloud("face vertices", get_points_from_handles(facet_vertices));
            //
            //         Voronoi_edge edge = delaunay_face_dual(face, delaunay);
            //         std::vector<Eigen::VectorXd> vv;
            //         vv.push_back(edge.vertex1.point);
            //         vv.push_back(edge.vertex2.point);
            //         std::vector<std::array<size_t, 2>> ve = {{0, 1}};
            //         polyscope::registerCurveNetwork("ve", vv, ve);
            //
            //         polyscope::show();
            //     }
            // }
        }

    vertices.resize(vertex_count);
    for (auto fc_vertex : fc_vertex_to_index) {
        vertices[fc_vertex.second] = fc_vertex.first;
    }

}
