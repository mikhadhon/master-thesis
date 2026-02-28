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

    // Find three finite basis vertices for the polygon plane
    int b0 = -1, b1 = -1, b2 = -1;
    for (int i = 0; i < nv; ++i) {
        if (b0 == -1) b0 = i;
        else if (b1 == -1) b1 = i;
        else if (b2 == -1) { b2 = i; break; }
    }

    // Build 4x4 system in exact FT:
    // [e1 | e2 | -f1 | -f2] * [s; t; u; w] = V[b0] - pi
    // e1 = pj - pi, e2 = pl - pi
    // f1 = V[b1] - V[b0], f2 = V[b2] - V[b0]
    FT m[4][5]; // augmented matrix [M | rhs]
    for (int k = 0; k < 4; ++k) {
        m[k][0] = pj[k] - pi[k];
        m[k][1] = pl[k] - pi[k];
        m[k][2] = -(vp.voronoi_vertices[b1].point[k] - vp.voronoi_vertices[b0].point[k]);
        m[k][3] = -(vp.voronoi_vertices[b2].point[k] - vp.voronoi_vertices[b0].point[k]);
        m[k][4] = vp.voronoi_vertices[b0].point[k] - pi[k];
    }

    // Gaussian elimination
    for (int col = 0; col < 4; ++col) {
        int pivot = -1;
        for (int row = col; row < 4; ++row) {
            if (m[row][col] != FT(0)) { pivot = row; break; }
        }
        if (pivot == -1) {
            std::cout << "pivot 0" << std::endl;
            return false;
        }
        if (pivot != col)
            for (int j = col; j < 5; ++j) std::swap(m[col][j], m[pivot][j]);
        for (int row = col + 1; row < 4; ++row) {
            FT factor = m[row][col] / m[col][col];
            for (int j = col; j < 5; ++j)
                m[row][j] -= factor * m[col][j];
        }
    }

    // Back substitution
    FT sol[4];
    for (int i = 3; i >= 0; --i) {
        sol[i] = m[i][4];
        for (int j = i + 1; j < 4; ++j)
            sol[i] -= m[i][j] * sol[j];
        sol[i] /= m[i][i];
    }

    FT s = sol[0], t = sol[1], u = sol[2], w = sol[3];

    // Triangle containment: s >= 0, t >= 0, s + t <= 1
    if (!(s >= FT(0) && t >= FT(0) && s + t <= FT(1))) {
        return false;
    }

    Vector pjpi = (Vector(pj) - pi);
    Vector spjpi = Vector(s * pjpi[0], s * pjpi[1], s * pjpi[2], s * pjpi[3]);
    Vector plpi = (Vector(pl) - pi);
    Vector tpjpi = Vector(t * plpi[0], t * plpi[1], t * plpi[2], t * plpi[3]);
    Vector p_int = pi + spjpi + tpjpi;
    out = p_int;
    // Voronoi polygon containment via SVD projection to 2D
    // Basis vectors from finite vertices
    Eigen::Matrix<double, 2, 4> A;
    for (int k = 0; k < 4; ++k) {
        A(0, k) = CGAL::to_double(vp.voronoi_vertices[b1].point[k] - vp.voronoi_vertices[b0].point[k]);
        A(1, k) = CGAL::to_double(vp.voronoi_vertices[b2].point[k] - vp.voronoi_vertices[b0].point[k]);
    }

    Eigen::JacobiSVD<Eigen::Matrix<double, 2, 4>> svd(A, Eigen::ComputeFullV);
    // Last two columns of V are the orthonormal basis plane's null space
    Eigen::Matrix<double, 2, 4> VT = svd.matrixV().rightCols<2>().transpose();

    // Debug: centroid of finite vertices is guaranteed inside the polygon
    FT cx(0), cy(0), cz(0), cw(0);
    int finite_count = 0;
    for (int i = 0; i < nv; ++i) {
        if (!vp.voronoi_vertices[i].is_infinite) {
            cx += vp.voronoi_vertices[i].point[0];
            cy += vp.voronoi_vertices[i].point[1];
            cz += vp.voronoi_vertices[i].point[2];
            cw += vp.voronoi_vertices[i].point[3];
            finite_count++;
        }
    }
    Point centroid(cx / finite_count, cy / finite_count, cz / finite_count, cw / finite_count);
    Point n1(VT.row(0).begin(), VT.row(0).end());
    Point n2(VT.row(1).begin(), VT.row(1).end());

    // Check containment: all edge determinants must have the same sign
    int sign = 0;
    for (const auto &edge : vp.voronoi_edges) {
        Point ref;
        Vector vec;

        if (!edge.vertex1.is_infinite && !edge.vertex2.is_infinite) {
            ref = edge.vertex1.point;
            vec = Vector(edge.vertex2.point) - edge.vertex1.point;
        }
        else if (!edge.vertex1.is_infinite && edge.vertex2.is_infinite) {
            ref = edge.vertex1.point;
            vec = edge.vertex2.infinite_direction;
        }
        else if (edge.vertex1.is_infinite && !edge.vertex2.is_infinite) {
            ref = edge.vertex2.point;
            vec = -edge.vertex1.infinite_direction;
        }
        else continue;
        auto orientation_test = [&](const Point& p) {
            std::vector<Point> simplex = {
                p,
                ref,
                Point(Vector(ref) + vec),
                Point(Vector(ref) + pi),
                Point(Vector(ref) + pj)
            };
            return orientation()(simplex.begin(), simplex.end());
        };

        auto o_point = orientation_test(p_int);
        auto o_centroid = orientation_test(centroid);
        if (o_point == CGAL::ZERO) {
            continue;
        }

        // if (o_point != o_centroid && o_centroid != CGAL::ZERO) {
        //     return false;
        // }
        if (sign == 0) {
            sign = (o_point > 0) ? 1 : -1;
        }
        else if ((sign == 1 && o_point < 0) || (sign == -1 && o_point > 0)) {
            return false;
        }
    }

    return true;
}

bool intersect_alt(const Face &df, Voronoi_face &vp, Point &out) {
    Vector O1 = df.face.vertex(0)->point();
    Vector u1 = Vector(df.face.vertex(1)->point()) - df.face.vertex(0)->point();
    Vector v1 = Vector(df.face.vertex(2)->point()) - df.face.vertex(0)->point();
    Vector O2 = vp.voronoi_vertices[0].point;
    Vector u2 = Vector(vp.voronoi_vertices[1].point) - vp.voronoi_vertices[0].point;
    Vector v2 = Vector(vp.voronoi_vertices[2].point) - vp.voronoi_vertices[0].point;
    // Set up the right-hand side and the negated vectors for the matrix
    Vector rhs = O2 - O1;
    Vector neg_u2 = -u2;
    Vector neg_v2 = -v2;

    // Calculate the main determinant
    FT det = CGAL::determinant(
        u1[0], v1[0], neg_u2[0], neg_v2[0],
        u1[1], v1[1], neg_u2[1], neg_v2[1],
        u1[2], v1[2], neg_u2[2], neg_v2[2],
        u1[3], v1[3], neg_u2[3], neg_v2[3]
    );

    // If determinant is 0, the planes are parallel or share a line/plane
    if (det == FT(0)) {
        return false;
    }

    // Cramer's rule for alpha (a)
    FT det_a = CGAL::determinant(
        rhs[0], v1[0], neg_u2[0], neg_v2[0],
        rhs[1], v1[1], neg_u2[1], neg_v2[1],
        rhs[2], v1[2], neg_u2[2], neg_v2[2],
        rhs[3], v1[3], neg_u2[3], neg_v2[3]
    );
    FT a = det_a / det;

    // Cramer's rule for beta (b)
    FT det_b = CGAL::determinant(
        u1[0], rhs[0], neg_u2[0], neg_v2[0],
        u1[1], rhs[1], neg_u2[1], neg_v2[1],
        u1[2], rhs[2], neg_u2[2], neg_v2[2],
        u1[3], rhs[3], neg_u2[3], neg_v2[3]
    );
    FT b = det_b / det;

    // Calculate the exact intersection point using Plane 1's parameters
    std::vector<FT> coords(4);
    for (int i = 0; i < 4; ++i) {
        coords[i] = O1[i] + a * u1[i] + b * v1[i];
    }

    Point p_int = Point(coords.begin(), coords.end());
    out = p_int;

    Eigen::Matrix<double, 2, 4> A;
    for (int k = 0; k < 4; ++k) {
        A(0, k) = CGAL::to_double(vp.voronoi_vertices[1].point[k] - vp.voronoi_vertices[0].point[k]);
        A(1, k) = CGAL::to_double(vp.voronoi_vertices[2].point[k] - vp.voronoi_vertices[0].point[k]);
    }

    Eigen::JacobiSVD<Eigen::Matrix<double, 2, 4>> svd(A, Eigen::ComputeFullV);
    // Last two columns of V are the orthonormal basis plane's null space
    Eigen::Matrix<double, 2, 4> VT = svd.matrixV().rightCols<2>().transpose();

    // Debug: centroid of finite vertices is guaranteed inside the polygon
    FT cx(0), cy(0), cz(0), cw(0);
    int finite_count = 0;
    const int nv = static_cast<int>(vp.voronoi_vertices.size());
    for (int i = 0; i < nv; ++i) {
        if (!vp.voronoi_vertices[i].is_infinite) {
            cx += vp.voronoi_vertices[i].point[0];
            cy += vp.voronoi_vertices[i].point[1];
            cz += vp.voronoi_vertices[i].point[2];
            cw += vp.voronoi_vertices[i].point[3];
            finite_count++;
        }
    }
    Point centroid(cx / finite_count, cy / finite_count, cz / finite_count, cw / finite_count);
    Point n1(VT.row(0).begin(), VT.row(0).end());
    Point n2(VT.row(1).begin(), VT.row(1).end());

    // Check containment: all edge determinants must have the same sign
    int sign = 0;
    for (const auto &edge : vp.voronoi_edges) {
        Point ref;
        Vector vec;

        if (!edge.vertex1.is_infinite && !edge.vertex2.is_infinite) {
            ref = edge.vertex1.point;
            vec = Vector(edge.vertex2.point) - edge.vertex1.point;
        }
        else if (!edge.vertex1.is_infinite && edge.vertex2.is_infinite) {
            ref = edge.vertex1.point;
            vec = edge.vertex2.infinite_direction;
        }
        else if (edge.vertex1.is_infinite && !edge.vertex2.is_infinite) {
            ref = edge.vertex2.point;
            vec = -edge.vertex1.infinite_direction;
        }
        else continue;
        auto orientation_test = [&](const Point& p) {
            std::vector<Point> simplex = {
                p,
                ref,
                Point(Vector(ref) + vec),
                Point(Vector(ref) + u1),
                Point(Vector(ref) + v1)
            };
            return orientation()(simplex.begin(), simplex.end());
        };

        auto o_point = orientation_test(p_int);
        auto o_centroid = orientation_test(centroid);
        if (o_point == CGAL::ZERO) {
            continue;
        }

        // if (o_point != o_centroid && o_centroid != CGAL::ZERO) {
        //     return false;
        // }
        if (sign == 0) {
            sign = (o_point > 0) ? 1 : -1;
        }
        else if ((sign == 1 && o_point < 0) || (sign == -1 && o_point > 0)) {
            return false;
        }
    }

    return true;
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
                        Voronoi_face current_voronoi_face = delaunay_face_dual(current_delaunay_face, delaunay);

                        std::vector<Voronoi_face> voronoi_faces = delaunay_edge_dual(current_edge, current_delaunay_face, delaunay);

                        std::vector<std::pair<Point, Face>> next;

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
                        // if (next.size() > 0) {
                        //     mulitple++;
                        //     // std::cout << next[0].first << std::endl;
                        //     // std::cout << next[1].first << std::endl;
                        //     // std::cout << "-----------" << std::endl;
                        //
                        //     // --- Debug visualization for multiple intersections ---
                        //
                        //     // Current Delaunay triangle
                        //     std::vector<Eigen::VectorXd> dbg_tri_4d = {
                        //         make_point_eigen(current_delaunay_face.face.vertex(0)->point()),
                        //         make_point_eigen(current_delaunay_face.face.vertex(1)->point()),
                        //         make_point_eigen(current_delaunay_face.face.vertex(2)->point())
                        //     };
                        //     std::vector<Eigen::VectorXd> dbg_tri_3d = stereo_projection(dbg_tri_4d);
                        //     std::vector<std::array<size_t, 3>> dbg_tri_faces = {{0, 1, 2}};
                        //     auto *dbg_tri = polyscope::registerSurfaceMesh("debug: current triangle", dbg_tri_3d, dbg_tri_faces);
                        //     dbg_tri->setSurfaceColor({0.2, 0.5, 1.0});
                        //     dbg_tri->setTransparency(0.5);
                        //
                        //     // Circumcenter of current Delaunay triangle
                        //     std::vector<Eigen::VectorXd> dbg_cc_4d = {make_point_eigen(current_center)};
                        //     std::vector<Eigen::VectorXd> dbg_cc_3d = stereo_projection(dbg_cc_4d);
                        //     auto *dbg_cc = polyscope::registerPointCloud("debug: circumcenter", dbg_cc_3d);
                        //     dbg_cc->setPointRadius(0.008);
                        //     dbg_cc->setPointColor({1.0, 0.0, 0.5});
                        //
                        //     // Current edge (highlighted)
                        //     std::vector<Eigen::VectorXd> dbg_edge_4d = {
                        //         make_point_eigen(current_edge.vertex1->point()),
                        //         make_point_eigen(current_edge.vertex2->point())
                        //     };
                        //     std::vector<Eigen::VectorXd> dbg_edge_3d = stereo_projection(dbg_edge_4d);
                        //     std::vector<std::array<size_t, 2>> dbg_edge_idx = {{0, 1}};
                        //     auto *dbg_edge = polyscope::registerCurveNetwork("debug: current edge", dbg_edge_3d, dbg_edge_idx);
                        //     dbg_edge->setColor({1.0, 0.0, 0.0});
                        //     dbg_edge->setRadius(0.005);
                        //
                        //     // s_prime intersection points
                        //     std::vector<Eigen::VectorXd> dbg_sp_4d;
                        //     for (const auto& [pt, _] : next) {
                        //         dbg_sp_4d.push_back(make_point_eigen(pt));
                        //     }
                        //     std::vector<Eigen::VectorXd> dbg_sp_3d = stereo_projection(dbg_sp_4d);
                        //     auto *dbg_sp = polyscope::registerPointCloud("debug: s_prime points", dbg_sp_3d);
                        //     dbg_sp->setPointRadius(0.008);
                        //     dbg_sp->setPointColor({1.0, 1.0, 0.0});
                        //
                        //     // Neighboring Delaunay triangles from next
                        //     std::vector<Eigen::VectorXd> dbg_nb_4d;
                        //     std::vector<std::array<size_t, 3>> dbg_nb_faces;
                        //     for (size_t ni = 0; ni < next.size(); ni++) {
                        //         size_t base = dbg_nb_4d.size();
                        //         dbg_nb_4d.push_back(make_point_eigen(next[ni].second.face.vertex(0)->point()));
                        //         dbg_nb_4d.push_back(make_point_eigen(next[ni].second.face.vertex(1)->point()));
                        //         dbg_nb_4d.push_back(make_point_eigen(next[ni].second.face.vertex(2)->point()));
                        //         dbg_nb_faces.push_back({base, base + 1, base + 2});
                        //     }
                        //     std::vector<Eigen::VectorXd> dbg_nb_3d = stereo_projection(dbg_nb_4d);
                        //     auto *dbg_nb = polyscope::registerSurfaceMesh("debug: neighbor triangles", dbg_nb_3d, dbg_nb_faces);
                        //     dbg_nb->setSurfaceColor({0.0, 1.0, 0.3});
                        //     dbg_nb->setTransparency(0.5);
                        //
                        //     polyscope::show();
                        // }

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
                                if (next_edge == current_edge) continue;
                                // edge_queue.emplace(s_prime, next_edge, next_face);
                                std::vector next_v = {next_face.face.vertex(0)->point(), next_face.face.vertex(1)->point(), next_face.face.vertex(2)->point()};
                                // edge_queue.emplace(circumcenter()(next_v.begin(), next_v.end()), next_edge, next_face);
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