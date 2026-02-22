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
        if (!vp.voronoi_vertices[i].is_infinite) {
            if (b0 == -1) b0 = i;
            else if (b1 == -1) b1 = i;
            else if (b2 == -1) { b2 = i; break; }
        }
    }
    if (b2 == -1) return false;

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
    if (!(s > FT(0) && t > FT(0) && s + t < FT(1))) {
        return false;
    }

    // Voronoi polygon containment in exact FT
    // Basis vectors from finite vertices
    FT fv1[4], fv2[4];
    for (int k = 0; k < 4; ++k) {
        fv1[k] = vp.voronoi_vertices[b1].point[k] - vp.voronoi_vertices[b0].point[k];
        fv2[k] = vp.voronoi_vertices[b2].point[k] - vp.voronoi_vertices[b0].point[k];
    }

    // Find two rows where [f1|f2] has a non-zero 2x2 determinant
    int r1 = -1, r2 = -1;
    for (int a = 0; a < 4 && r1 == -1; ++a)
        for (int b = a + 1; b < 4 && r1 == -1; ++b)
            if (fv1[a] * fv2[b] - fv1[b] * fv2[a] != FT(0))
                { r1 = a; r2 = b; }
    if (r1 == -1) {
        return false;
    }

    FT det_f = fv1[r1] * fv2[r2] - fv1[r2] * fv2[r1];

    // Project finite polygon vertices to 2D (u_k, w_k) coordinates
    std::vector<std::pair<FT, FT>> P(nv);
    for (int i = 0; i < nv; ++i) {
        if (vp.voronoi_vertices[i].is_infinite) continue;
        FT dv_r1 = vp.voronoi_vertices[i].point[r1] - vp.voronoi_vertices[b0].point[r1];
        FT dv_r2 = vp.voronoi_vertices[i].point[r2] - vp.voronoi_vertices[b0].point[r2];
        P[i] = { (dv_r1 * fv2[r2] - dv_r2 * fv2[r1]) / det_f,
                  (fv1[r1] * dv_r2 - fv1[r2] * dv_r1) / det_f };
    }

    // Helper: project a 4D direction vector to 2D
    auto project_dir = [&](const Vector &dir) -> std::pair<FT, FT> {
        FT d_r1 = dir[r1];
        FT d_r2 = dir[r2];
        return { (d_r1 * fv2[r2] - d_r2 * fv2[r1]) / det_f,
                 (fv1[r1] * d_r2 - fv1[r2] * d_r1) / det_f };
    };

    // Winding is always CCW with basis (b0,b1,b2) → P[b0]=(0,0), P[b1]=(1,0), P[b2]=(0,1)
    // Check that (u, w) is inside every edge
    for (int i = 0; i < nv; ++i) {
        int j = (i + 1) % nv;
        bool inf_i = vp.voronoi_vertices[i].is_infinite;
        bool inf_j = vp.voronoi_vertices[j].is_infinite;

        if (inf_i && inf_j) continue;

        FT cross;
        if (!inf_i && !inf_j) {
            // Both finite: standard cross product
            FT ex = P[j].first  - P[i].first;
            FT ey = P[j].second - P[i].second;
            cross = ex * (w - P[i].second) - ey * (u - P[i].first);
        } else if (!inf_i && inf_j) {
            // Finite → infinite: edge direction is the infinite direction of V_j
            auto [dx, dy] = project_dir(vp.voronoi_vertices[j].infinite_direction);
            cross = dx * (w - P[i].second) - dy * (u - P[i].first);
        } else {
            // Infinite → finite: edge direction is -d_i (from infinity toward V_j)
            auto [dx, dy] = project_dir(vp.voronoi_vertices[i].infinite_direction);
            cross = (-dx) * (w - P[j].second) - (-dy) * (u - P[j].first);
        }

        if (cross < FT(0)) {
            return false;
        }
    }

    Vector pjpi = (Vector(pj) - pi);
    Vector spjpi = Vector(s * pjpi[0], s * pjpi[1], s * pjpi[2], s * pjpi[3]);
    Vector plpi = (Vector(pl) - pi);
    Vector tpjpi = Vector(t * plpi[0], t * plpi[1], t * plpi[2], t * plpi[3]);
    out = pi + spjpi + tpjpi;
    return true;
}

bool intersect_ray_vp(const Point &origin, const Point &through, Voronoi_face &vp, Point &out) {
    const int nv = static_cast<int>(vp.voronoi_vertices.size());

    // Find three finite basis vertices for the polygon plane
    int b0 = -1, b1 = -1, b2 = -1;
    for (int i = 0; i < nv; ++i) {
        if (!vp.voronoi_vertices[i].is_infinite) {
            if (b0 == -1) b0 = i;
            else if (b1 == -1) b1 = i;
            else if (b2 == -1) { b2 = i; break; }
        }
    }
    if (b2 == -1) return false;

    // Ray: origin + λ * (through - origin), λ ≥ 0
    // Voronoi plane: V[b0] + u * f1 + w * f2
    // where f1 = V[b1] - V[b0], f2 = V[b2] - V[b0]
    //
    // Intersection: λ*d - u*f1 - w*f2 = V[b0] - origin
    // where d = through - origin
    // This is 4 equations (one per dimension), 3 unknowns (λ, u, w)
    FT m[4][4]; // augmented matrix [M | rhs]
    for (int k = 0; k < 4; ++k) {
        m[k][0] = through[k] - origin[k];
        m[k][1] = -(vp.voronoi_vertices[b1].point[k] - vp.voronoi_vertices[b0].point[k]);
        m[k][2] = -(vp.voronoi_vertices[b2].point[k] - vp.voronoi_vertices[b0].point[k]);
        m[k][3] = vp.voronoi_vertices[b0].point[k] - origin[k];
    }

    // Gaussian elimination on 4x3 coefficient part (3 pivot columns, 4 rows)
    for (int col = 0; col < 3; ++col) {
        int pivot = -1;
        for (int row = col; row < 4; ++row) {
            if (m[row][col] != FT(0)) { pivot = row; break; }
        }
        if (pivot == -1) {
            return false;
        } // rank < 3, degenerate
        if (pivot != col)
            for (int j = col; j < 4; ++j) std::swap(m[col][j], m[pivot][j]);
        for (int row = col + 1; row < 4; ++row) {
            FT factor = m[row][col] / m[col][col];
            for (int j = col; j < 4; ++j)
                m[row][j] -= factor * m[col][j];
        }
    }

    // After elimination, row 3 has columns 0-2 zeroed.
    // For the system to be consistent, the rhs must also be zero.
    if (m[3][3] != FT(0)) {
        return false;
    }

    // Back substitution on rows 0-2
    FT sol[3];
    for (int i = 2; i >= 0; --i) {
        sol[i] = m[i][3];
        for (int j = i + 1; j < 3; ++j)
            sol[i] -= m[i][j] * sol[j];
        sol[i] /= m[i][i];
    }

    FT lambda = sol[0], u = sol[1], w = sol[2];

    // Ray: λ ≥ 0
    if (lambda < FT(0)) {
        return false;
    }

    // Voronoi polygon containment in exact FT
    // Basis vectors from finite vertices
    FT fv1[4], fv2[4];
    for (int k = 0; k < 4; ++k) {
        fv1[k] = vp.voronoi_vertices[b1].point[k] - vp.voronoi_vertices[b0].point[k];
        fv2[k] = vp.voronoi_vertices[b2].point[k] - vp.voronoi_vertices[b0].point[k];
    }

    // Find two rows where [f1|f2] has a non-zero 2x2 determinant
    int r1 = -1, r2 = -1;
    for (int a = 0; a < 4 && r1 == -1; ++a)
        for (int b = a + 1; b < 4 && r1 == -1; ++b)
            if (fv1[a] * fv2[b] - fv1[b] * fv2[a] != FT(0))
                { r1 = a; r2 = b; }
    if (r1 == -1) {
        return false;
    }

    FT det_f = fv1[r1] * fv2[r2] - fv1[r2] * fv2[r1];

    // Project finite polygon vertices to 2D (u_k, w_k) coordinates
    std::vector<std::pair<FT, FT>> P(nv);
    for (int i = 0; i < nv; ++i) {
        if (vp.voronoi_vertices[i].is_infinite) continue;
        FT dv_r1 = vp.voronoi_vertices[i].point[r1] - vp.voronoi_vertices[b0].point[r1];
        FT dv_r2 = vp.voronoi_vertices[i].point[r2] - vp.voronoi_vertices[b0].point[r2];
        P[i] = { (dv_r1 * fv2[r2] - dv_r2 * fv2[r1]) / det_f,
                  (fv1[r1] * dv_r2 - fv1[r2] * dv_r1) / det_f };
    }

    // Helper: project a 4D direction vector to 2D
    auto project_dir = [&](const Vector &dir) -> std::pair<FT, FT> {
        FT d_r1 = dir[r1];
        FT d_r2 = dir[r2];
        return { (d_r1 * fv2[r2] - d_r2 * fv2[r1]) / det_f,
                 (fv1[r1] * d_r2 - fv1[r2] * d_r1) / det_f };
    };

    // Winding is always CCW with basis (b0,b1,b2) → P[b0]=(0,0), P[b1]=(1,0), P[b2]=(0,1)
    // Check that (u, w) is inside every edge
    for (int i = 0; i < nv; ++i) {
        int j = (i + 1) % nv;
        bool inf_i = vp.voronoi_vertices[i].is_infinite;
        bool inf_j = vp.voronoi_vertices[j].is_infinite;

        if (inf_i && inf_j) continue;

        FT cross;
        if (!inf_i && !inf_j) {
            // Both finite: standard cross product
            FT ex = P[j].first  - P[i].first;
            FT ey = P[j].second - P[i].second;
            cross = ex * (w - P[i].second) - ey * (u - P[i].first);
        } else if (!inf_i && inf_j) {
            // Finite → infinite: edge direction is the infinite direction of V_j
            auto [dx, dy] = project_dir(vp.voronoi_vertices[j].infinite_direction);
            cross = dx * (w - P[i].second) - dy * (u - P[i].first);
        } else {
            // Infinite → finite: edge direction is -d_i (from infinity toward V_j)
            auto [dx, dy] = project_dir(vp.voronoi_vertices[i].infinite_direction);
            cross = (-dx) * (w - P[j].second) - (-dy) * (u - P[j].first);
        }

        if (cross <= FT(0)) {
            return false;
        }
    }

    Vector d = (Vector(through) - origin);
    Vector ld = Vector(lambda * d[0], lambda * d[1], lambda * d[2], lambda * d[3]);
    out = origin + ld;
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
                        Point driver = midpoint()(Point(current_edge.vertex1->point()), Point(current_edge.vertex2->point()));
                        std::vector df_points = {current_delaunay_face.face.vertex(0)->point(), current_delaunay_face.face.vertex(1)->point(), current_delaunay_face.face.vertex(2)->point()};
                        Point through = circumcenter()(df_points.begin(), df_points.end());
                        std::optional<Point> s_prime;

                        std::optional<Face> next_face;
                        bool intersected = false;
                        std::optional<Voronoi_face> debug;
                        for (auto v_face : voronoi_faces) {
                            if (Point out; intersect_df_vp(current_delaunay_face, v_face, out)) {
                                intersected = true;
                                debug = v_face;
                                next_face = v_face.dual;
                                s_prime = out;
                                break;
                            }
                        }
                        // for (auto v_face : voronoi_faces) {
                        //     if (Point out; intersect_ray_vp(driver, through, v_face, out)) {
                        //         intersected = true;
                        //         next_face = v_face.dual;
                        //         s_prime = out;
                        //         break;
                        //     }
                        // }
                        if (!intersected) {
                            missed_intersections++;
                        }

                        if (next_face.has_value() && s_prime.has_value()) {
                            std::array<size_t, 3> fc_face1{}, fc_face2{};
                            std::array fc_center = {CGAL::to_double(current_center[0]), CGAL::to_double(current_center[1]), CGAL::to_double(current_center[2]), CGAL::to_double(current_center[3])};
                            std::array fc_s_prime = {CGAL::to_double(s_prime.value()[0]), CGAL::to_double(s_prime.value()[1]), CGAL::to_double(s_prime.value()[2]), CGAL::to_double(s_prime.value()[3])};
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

                            for (int i = 0; i <= next_face->face.face_dimension(); i++) {
                                Edge next_edge(next_face->face.vertex(i), next_face->face.vertex((i + 1)%3), next_face->face.vertex((i + 2)%3));
                                if (next_edge == current_edge) continue;
                                edge_queue.emplace(s_prime.value(), next_edge, next_face.value());
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
}