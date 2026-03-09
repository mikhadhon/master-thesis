#include <queue>
#include <polyscope/polyscope.h>
#include <polyscope/surface_mesh.h>
#include <polyscope/curve_network.h>
#include <polyscope/point_cloud.h>

#include "FlowComplex.h"
#include "utils.h"

bool intersect_df_vp(const Face &df, Voronoi_face &vp, Point &out) {
    Point pi = df.face.vertex(0)->point();
    Point pj = df.face.vertex(1)->point();
    Point pl = df.face.vertex(2)->point();

    const int nv = static_cast<int>(vp.voronoi_vertices.size());

    int b0 = -1, b1 = -1, b2 = -1;
    for (int i = 0; i < nv; ++i) {
        if (b0 == -1) b0 = i;
        else if (b1 == -1) b1 = i;
        else { b2 = i; break; }
    }

    LA_Matrix M(4,4);
    for (int k = 0; k < 4; ++k) {
        auto row = M.row_begin(k);
        *row++ = pj[k] - pi[k];
        *row++ = pl[k] - pi[k];
        *row++ = -(vp.voronoi_vertices[b1].point[k] - vp.voronoi_vertices[b0].point[k]);
        *row = -(vp.voronoi_vertices[b2].point[k] - vp.voronoi_vertices[b0].point[k]);
    }
    LA_Vector b(4);
    auto coord = b.begin();
    for (int k = 0; k < 4; ++k) {
        *coord++ = vp.voronoi_vertices[b0].point[k] - pi[k];
    }

    LA_Vector x(4);
    FT D;
    LA::linear_solver(M, b, x, D);

    FT s = x[0], t = x[1], u = x[2], w = x[3];

    // Triangle containment: s >= 0, t >= 0, s + t <= 1
    if (!(s > FT(0) && t > FT(0) && s + t < FT(1))) {
        return false;
    }

    Vector pjpi = (Vector(pj) - pi);
    Vector spjpi = Vector(s * pjpi[0], s * pjpi[1], s * pjpi[2], s * pjpi[3]);
    Vector plpi = (Vector(pl) - pi);
    Vector tpjpi = Vector(t * plpi[0], t * plpi[1], t * plpi[2], t * plpi[3]);
    Vector p_int = pi + spjpi + tpjpi;
    out = p_int;

    Voronoi_vertex v_vertex_base = vp.voronoi_vertices.front();
    auto v_vertex_it = vp.voronoi_vertices.begin();
    while (v_vertex_base.is_infinite) {
        v_vertex_base = *(++v_vertex_it);
    }
    LA_Vector v_0(v_vertex_base.point.cartesian_begin(), v_vertex_base.point.cartesian_end());

    for (int i = 0; i < vp.voronoi_vertices.size(); ++i) {
        if (vp.voronoi_vertices[i] == v_vertex_base) continue;
        int next = (i + 1) % vp.voronoi_vertices.size();
        if (vp.voronoi_vertices[next] == v_vertex_base) next = (i + 2) % vp.voronoi_vertices.size();
        LA_Vector v_1(vp.voronoi_vertices[i].point.cartesian_begin(), vp.voronoi_vertices[i].point.cartesian_end());
        LA_Vector v_2(vp.voronoi_vertices[next].point.cartesian_begin(), vp.voronoi_vertices[next].point.cartesian_end());
        std::vector matrix_cols = {v_1 - v_0, v_2 - v_0};

        LA_Matrix M_P(matrix_cols.begin(), matrix_cols.end());
        LA_Vector b_P(4);
        auto c = b_P.begin();
        for (int k = 0; k < 4; ++k) {
            *c++ = out[k] - v_0[k];
        }
        LA_Vector x_p(4);
        FT D_p;

        LA::linear_solver(M_P, b_P, x_p, D_p);

        if (x_p[0] >= FT(0) && x_p[1] >= FT(0) && x_p[0] + x_p[1] <= FT(1)) {
            return true;
        }
    }

    return false;
}

void flow_complex(Delaunay &delaunay, std::vector<Eigen::VectorXd> &vertices, std::vector<std::array<size_t, 3>> &faces, std::map<Delaunay::Vertex_handle, size_t> &vertex_to_index) {
    std::map<std::array<double, 4>, size_t> fc_vertex_to_index;
    int vertex_count = static_cast<int>(vertices.size());
    std::vector<Eigen::VectorXd> index_two_points;
    std::vector<std::array<size_t, 3>> non_gabriel_faces;
    std::vector<std::array<size_t, 2>> gabriel_edges;
    std::vector<Eigen::VectorXd> V_projected = stereo_projection(vertices);
    int missed_intersections = 0;
    int mulitple = 0;

    for (auto delaunay_faces = get_delaunay_faces(delaunay); const auto & delaunay_face : delaunay_faces) {
        if (!delaunay.is_infinite(delaunay_face.face)) {
            if (is_index_two_critical_point(delaunay_face, delaunay)) {
                std::queue<std::tuple<Point, Edge, Face>> edge_queue;

                std::vector cgal_points = {delaunay_face.face.vertex(0)->point(), delaunay_face.face.vertex(1)->point(), delaunay_face.face.vertex(2)->point()};
                Point cgal_center = circumcenter()(cgal_points.begin(), cgal_points.end());

                edge_queue.emplace(cgal_center, Edge(delaunay_face.face.vertex(0), delaunay_face.face.vertex(1), delaunay_face.face.vertex(2)), delaunay_face);
                edge_queue.emplace(cgal_center, Edge(delaunay_face.face.vertex(0), delaunay_face.face.vertex(2), delaunay_face.face.vertex(1)), delaunay_face);
                edge_queue.emplace(cgal_center, Edge(delaunay_face.face.vertex(1), delaunay_face.face.vertex(2), delaunay_face.face.vertex(0)), delaunay_face);

                while (!edge_queue.empty()) {
                    std::tuple<Point, Edge, Face> current = edge_queue.front();
                    edge_queue.pop();
                    Point current_center = std::get<0>(current);
                    Edge current_edge = std::get<1>(current);
                    Face current_delaunay_face = std::get<2>(current);
                    if (is_gabriel(current_edge, delaunay)) {
                        std::array fc_center = {CGAL::to_double(current_center[0]), CGAL::to_double(current_center[1]), CGAL::to_double(current_center[2]), CGAL::to_double(current_center[3])};
                        auto fc_index = fc_vertex_to_index.find(fc_center);
                        std::array<size_t, 3> fc_face{};
                        if (fc_index == fc_vertex_to_index.end()) {
                            fc_vertex_to_index[fc_center] = vertex_count++;
                        }
                        fc_face[0] = vertex_to_index[current_edge.vertex1];
                        fc_face[1] = vertex_to_index[current_edge.vertex2];
                        fc_face[2] = fc_vertex_to_index[fc_center];
                        faces.push_back(fc_face);
                        gabriel_edges.push_back({vertex_to_index[current_edge.vertex1], vertex_to_index[current_edge.vertex2]});
                    }
                    else {
                        std::vector<Delaunay::Full_cell_handle> neighboring_cells;

                        Point edge_mid = midpoint()(current_edge.vertex1->point(), current_edge.vertex2->point());

                        std::vector<std::pair<Point, Face>> next;

                        Voronoi_face current_voronoi_face = delaunay_face_dual(current_delaunay_face, delaunay);

                        std::vector<Voronoi_face> voronoi_faces = delaunay_edge_dual(current_edge, current_delaunay_face, delaunay);

                        bool intersected = false;
                        for (auto v_face : voronoi_faces) {
                            Point out;
                            if (intersect_df_vp(current_delaunay_face, v_face, out)) {
                                intersected = true;
                                next.push_back({out, v_face.dual});
                            }
                        }

                        if (!intersected) {
                            missed_intersections++;
                        }
                        if (next.size() > 1) {
                            std::cout << "Multiple intersections: " << next.size() << std::endl;
                            mulitple++;
                        }

                        for (const auto& [s_prime, next_face]: next) {
                            std::array<size_t, 3> fc_face1{}, fc_face2{};
                            std::array fc_center = {CGAL::to_double(current_center[0]), CGAL::to_double(current_center[1]), CGAL::to_double(current_center[2]), CGAL::to_double(current_center[3])};
                            std::array fc_s_prime = {CGAL::to_double(s_prime[0]), CGAL::to_double(s_prime[1]), CGAL::to_double(s_prime[2]), CGAL::to_double(s_prime[3])};
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

                            for (int i = 0; i <= next_face.face.face_dimension(); i++) {
                                Edge next_edge(next_face.face.vertex(i), next_face.face.vertex((i + 1)%3), next_face.face.vertex((i + 2)%3));
                                if (next_edge == current_edge) {
                                    continue;
                                }
                                edge_queue.emplace(s_prime, next_edge, next_face);
                            }
                        }
                    }
                }
            }
        }
    }
    vertices.resize(vertex_count);
    for (auto [vertex, index] : fc_vertex_to_index) {
        Eigen::VectorXd new_vertex(delaunay.maximal_dimension());
        new_vertex << vertex[0], vertex[1], vertex[2], vertex[3];
        vertices[index] = new_vertex;
    }
    polyscope::registerCurveNetwork("Gabriel edges", V_projected, gabriel_edges);
    //polyscope::registerSurfaceMesh("non-gabriel", vertices, non_gabriel_faces);
    std::cout << "Missed intersections: " << missed_intersections << std::endl;
    std::cout << "Multiple intersections: " << mulitple << std::endl;
}