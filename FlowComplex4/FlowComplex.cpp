#include <queue>
#include <polyscope/polyscope.h>
#include <polyscope/surface_mesh.h>
#include <polyscope/curve_network.h>

#include "FlowComplex.h"
#include "utils.h"

bool intersect_df_vp(const Face &df, Voronoi_face &vp, Point &out, const Edge &current_edge, const Point &current_center, const Delaunay &dt) {
    Point pi = current_edge.vertex1->point();
    Point pj = current_edge.vertex2->point();
    Point pl = current_center;

    const int nv = static_cast<int>(vp.voronoi_vertices.size());

    int b0 = -1, b1 = -1, b2 = -1;
    for (int i = 0; i < nv; ++i) {
        if (b0 == -1) b0 = i;
        else if (b1 == -1) b1 = i;
        else { b2 = i; break; }
    }

    LA_Matrix M(4,4);
    for (int k = 0; k < 4; ++k) {
        LA_Matrix::row_iterator row = M.row_begin(k);
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

    const FT s = x[0]/D;
    const FT t = x[1]/D;

    Vector pjpi = (Vector(pj) - pi);
    Vector spjpi = Vector(s * pjpi[0], s * pjpi[1], s * pjpi[2], s * pjpi[3]);
    Vector plpi = (Vector(pl) - pi);
    Vector tplpi = Vector(t * plpi[0], t * plpi[1], t * plpi[2], t * plpi[3]);
    Vector p_int = pi + spjpi + tplpi;
    out = p_int;

    // Triangle containment: s >= 0, t >= 0, s + t <= 1
    if (!(s >= FT(0) && t >= FT(0) && s + t <= FT(1))) {
        return false;
    }

    Delaunay::Vertex_handle dual_point;
    for (int i = 0; i < 3; i++) {
        if (vp.dual.face.vertex(i) != current_edge.vertex1 && vp.dual.face.vertex(i) != current_edge.vertex2) {
            dual_point = vp.dual.face.vertex(i);
        }
    }
    std::vector<Delaunay::Full_cell_handle> cell_neighbors;
    get_incident_cells_to_edge(current_edge, dt, cell_neighbors);
    std::set<Delaunay::Vertex_handle> neighboring_vertices;

    for (const auto cell : cell_neighbors) {
        for (auto v = cell->vertices_begin(); v != cell->vertices_end(); ++v) {
            neighboring_vertices.insert(*v);
        }
    }
    const FT sq_dist_out_dual = squared_distance()(out, current_edge.vertex1->point());
    for (Delaunay::Vertex_handle v :neighboring_vertices) {
        if (v == current_edge.vertex1 || v == current_edge.vertex2 || v == dual_point) continue;
        if (FT sq_dist_out_neighbor = squared_distance()(out, v->point()); sq_dist_out_neighbor < sq_dist_out_dual) {
            return false;
        }
    }

    return true;
}

void index_constructed_vertex(const std::array<double, 4> &constructed_vertex, std::map<std::array<double, 4>, size_t> &map, int &total_count) {
    if (const auto already_indexed = map.find(constructed_vertex); already_indexed == map.end()) {
        map[constructed_vertex] = total_count++;
    }
}

std::array<size_t, 3> make_fc_face_gabriel(
    const Delaunay::Vertex_handle e1,
    const Delaunay::Vertex_handle e2,
    const std::array<double, 4> &center,
    std::map<Delaunay::Vertex_handle, size_t> delaunay_map,
    std::map<std::array<double, 4>, size_t> constructed_map
) {
    std::array<size_t, 3> fc_face{};
    fc_face[0] = delaunay_map[e1];
    fc_face[1] = delaunay_map[e2];
    fc_face[2] = constructed_map[center];

    return fc_face;
}
bool intersect_ray_vp(const Face &df, Voronoi_face &vp, Point &out, const Edge &current_edge, const Point &current_center, const Delaunay &dt) {
    // Point pi = df.face.vertex(0)->point();
    // Point pj = df.face.vertex(1)->point();
    // Point pl = df.face.vertex(2)->point();
    Point pi = current_edge.vertex1->point();
    Point pj = current_edge.vertex2->point();
    Point pl = current_center;

    Point driver = midpoint()(current_edge.vertex1->point(), current_edge.vertex2->point());
    Vector ray_dir = Vector(current_center) - driver;

    std::vector v_debug = {pi, pj, pl};

    const int nv = static_cast<int>(vp.voronoi_vertices.size());
    int b0 = -1;
    int b1 = -1;
    int b2 = -1;
    for (int i = 0; i < nv; ++i) {
        if (b0 == -1) {
            b0 = i;
        }
        else if (b1 == -1) {b1 = i;}
        else { b2 = i; break; }
    }

    LA_Matrix M(3,4);
    for (int k = 0; k < 3; ++k) {
        auto row = M.row_begin(k);
        *row++ = ray_dir[k];
        *row++ = -(vp.voronoi_vertices[b1].point[k] - vp.voronoi_vertices[b0].point[k]);
        *row = -(vp.voronoi_vertices[b2].point[k] - vp.voronoi_vertices[b0].point[k]);
    }
    LA_Vector b(4);
    auto coord = b.begin();
    for (int k = 0; k < 4; ++k) {
        *coord++ = vp.voronoi_vertices[b0].point[k] - driver[k];
    }

    LA_Vector x(4);
    FT D;
    LA::linear_solver(M, b, x, D);

    FT s = x[0]/D, t = x[1]/D;

    Vector s_driver_dir = Vector(s * ray_dir[0], s * ray_dir[1], s * ray_dir[2], s * ray_dir[3]);

    Vector p_int = driver + s_driver_dir;
    out = p_int;

    // Triangle containment: s >= 0, t >= 0, s + t <= 1
    if (!(s >= FT(0) && s <= FT(1))) {
        return false;
    }

    Delaunay::Vertex_handle dual_point;
    for (int i = 0; i < 3; i++) {
        if (vp.dual.face.vertex(i) != current_edge.vertex1 && vp.dual.face.vertex(i) != current_edge.vertex2) {
            dual_point = vp.dual.face.vertex(i);
        }
    }
    std::vector<Delaunay::Full_cell_handle> cell_neighbors;
    get_incident_cells_to_edge(current_edge, dt, cell_neighbors);
    std::set<Delaunay::Vertex_handle> neighboring_vertices;

    for (const auto cell : cell_neighbors) {
        for (auto v = cell->vertices_begin(); v != cell->vertices_end(); ++v) {
            neighboring_vertices.insert(*v);
        }
    }
    FT sq_dist_out_dual = squared_distance()(out, current_edge.vertex1->point());
    for (Delaunay::Vertex_handle v :neighboring_vertices) {
        if (v == current_edge.vertex1 || v == current_edge.vertex2 || v == dual_point) continue;
        FT sq_dist_out_neighbor = squared_distance()(out, v->point());
        if (sq_dist_out_neighbor < sq_dist_out_dual) {
            return false;
        }
    }

    return true;
}

std::array<size_t, 3> make_fc_face_non_gabriel(
    const Delaunay::Vertex_handle e,
    const std::array<double, 4> &center,
    const std::array<double, 4> &s_prime,
    std::map<Delaunay::Vertex_handle, size_t> delaunay_map,
    std::map<std::array<double, 4>, size_t> constructed_map
) {
    std::array<size_t, 3> fc_face{};
    fc_face[0] = delaunay_map[e];
    fc_face[1] = constructed_map[s_prime];
    fc_face[2] = constructed_map[center];

    return fc_face;
}

void index_two_stable_manifold(const Delaunay &dt, const Face &starting_face, std::vector<std::array<size_t, 3>> &faces, helper_data &data) {
    std::queue<std::tuple<Point, Edge, Face>> edge_queue;
    std::vector<Edge> boundary;

    std::vector cgal_points = {starting_face.face.vertex(0)->point(), starting_face.face.vertex(1)->point(), starting_face.face.vertex(2)->point()};
    Point cgal_center = circumcenter()(cgal_points.begin(), cgal_points.end());

    edge_queue.emplace(cgal_center, Edge(starting_face.face.vertex(0), starting_face.face.vertex(1), starting_face.face.vertex(2)), starting_face);
    edge_queue.emplace(cgal_center, Edge(starting_face.face.vertex(0), starting_face.face.vertex(2), starting_face.face.vertex(1)), starting_face);
    edge_queue.emplace(cgal_center, Edge(starting_face.face.vertex(1), starting_face.face.vertex(2), starting_face.face.vertex(0)), starting_face);

    std::vector<std::array<size_t, 3>> sm_faces;

    while (!edge_queue.empty()) {
        std::tuple<Point, Edge, Face> current = edge_queue.front();
        edge_queue.pop();
        Point current_center = std::get<0>(current);
        Edge current_edge = std::get<1>(current);
        Face current_delaunay_face = std::get<2>(current);

        if (is_gabriel(current_edge, dt)) {
            std::array center_converted = {CGAL::to_double(current_center[0]), CGAL::to_double(current_center[1]), CGAL::to_double(current_center[2]), CGAL::to_double(current_center[3])};
            index_constructed_vertex(center_converted, data.fc_vertex_to_index, data.vertex_count);

            std::array<size_t, 3> fc_face = make_fc_face_gabriel(current_edge.vertex1, current_edge.vertex2, center_converted, data.vertex_to_index, data.fc_vertex_to_index);

            faces.push_back(fc_face);
            sm_faces.push_back(fc_face);
            data.gabriel_edges.push_back({data.vertex_to_index[current_edge.vertex1], data.vertex_to_index[current_edge.vertex2]});
            boundary.push_back(current_edge);
            if (!data.gabriel_topology.contains(current_edge)) {
                data.gabriel_topology[current_edge] = 1;
            }
            else {
                data.gabriel_topology[current_edge]++;
            }
        }
        else {
            Point edge_mid = midpoint()(current_edge.vertex1->point(), current_edge.vertex2->point());

            Voronoi_face current_voronoi_face = delaunay_face_dual(current_delaunay_face, dt);
            std::vector<Voronoi_face> voronoi_faces = delaunay_edge_dual(current_edge, current_delaunay_face, dt);

            Point s_prime;
            std::optional<Face> next_face;
            for (Voronoi_face v_face : voronoi_faces) {
                if (Point out; intersect_df_vp(current_delaunay_face, v_face, out, current_edge, current_center, dt)) {
                    s_prime = out;
                    next_face = v_face.dual;
                    break;
                }
            }

            std::array fc_center = {CGAL::to_double(current_center[0]), CGAL::to_double(current_center[1]), CGAL::to_double(current_center[2]), CGAL::to_double(current_center[3])};
            std::array fc_s_prime = {CGAL::to_double(s_prime[0]), CGAL::to_double(s_prime[1]), CGAL::to_double(s_prime[2]), CGAL::to_double(s_prime[3])};

            index_constructed_vertex(fc_center, data.fc_vertex_to_index, data.vertex_count);
            index_constructed_vertex(fc_s_prime, data.fc_vertex_to_index, data.vertex_count);

            std::array<size_t, 3> fc_face1 = make_fc_face_non_gabriel(current_edge.vertex1, fc_center, fc_s_prime, data.vertex_to_index, data.fc_vertex_to_index);
            std::array<size_t, 3> fc_face2 = make_fc_face_non_gabriel(current_edge.vertex2, fc_center, fc_s_prime, data.vertex_to_index, data.fc_vertex_to_index);

            faces.push_back(fc_face1);
            faces.push_back(fc_face2);
            sm_faces.push_back(fc_face1);
            sm_faces.push_back(fc_face2);

            for (int i = 0; i <= next_face->face.face_dimension(); i++) {
                Edge next_edge(next_face->face.vertex(i), next_face->face.vertex((i + 1)%3), next_face->face.vertex((i + 2)%3));
                if (next_edge == current_edge) {
                    continue;
                }
                edge_queue.emplace(s_prime, next_edge, *next_face);
            }
        }
    }
    Saddle_2 saddle_2(starting_face, delaunay_face_dual(starting_face, dt));
    stable_manifold_2 sm(boundary, saddle_2, sm_faces);
    data.stable_manifolds.push_back(sm);

    Voronoi_face triangle_dual = delaunay_face_dual(starting_face, dt);
    FT distance_saddle = squared_distance()(cgal_center, starting_face.face.vertex(0)->point());
    for (Voronoi_vertex v : triangle_dual.voronoi_vertices) {
        FT distance_voronoi = v.is_infinite ? FT(std::numeric_limits<double>::max()) : squared_distance()(v.point, starting_face.face.vertex(0)->point());
        FT distance_heuristic = distance_voronoi - distance_saddle;

        data.valid_pairs.push_back(valid_pair(sm, v, CGAL::abs(distance_heuristic)));
    }
}

void flow_complex(Delaunay &delaunay, std::vector<Eigen::VectorXd> &vertices, std::vector<std::array<size_t, 3>> &faces, helper_data &data) {
    const std::vector<Eigen::VectorXd> V_projected = stereo_projection(vertices);

    for (const std::unordered_set<Face, FaceHash> delaunay_faces = get_delaunay_faces(delaunay); const Face &delaunay_face : delaunay_faces) {
        if (is_index_two_critical_point(delaunay_face, delaunay)) {
            index_two_stable_manifold(delaunay, delaunay_face, faces, data);
        }
    }
    vertices.resize(data.vertex_count);
    for (auto [vertex, index] : data.fc_vertex_to_index) {
        Eigen::VectorXd new_vertex(delaunay.maximal_dimension());
        new_vertex << vertex[0], vertex[1], vertex[2], vertex[3];
        vertices[index] = new_vertex;
    }
    polyscope::registerCurveNetwork("Gabriel edges", V_projected, data.gabriel_edges);
}

void reduce_flow_complex(Delaunay &dt, helper_data &data) {
    std::sort(data.valid_pairs.begin(), data.valid_pairs.end());

    while (true) {
        // Step 5: Find the valid pair with minimum distance that fulfills the topology constraint
        auto best_it = std::find_if(data.valid_pairs.begin(), data.valid_pairs.end(),
            [&](const valid_pair &vp) {
                return std::any_of(vp.sm.gabriel_edges.begin(), vp.sm.gabriel_edges.end(),
                    [&](const Edge &e) { return data.gabriel_topology[e] > 2; });
            });

        if (best_it == data.valid_pairs.end())
            break;

        valid_pair selected = *best_it;
        const Voronoi_vertex &b = selected.v;

        // Step 7: Collect all other maxima paired with saddle a (all voronoi vertices except b)
        std::vector<Voronoi_vertex> other_maxima;
        for (const auto &v : selected.sm.saddle_2.dual.voronoi_vertices) {
            if (!(v == b))
                other_maxima.push_back(v);
        }

        // Step 6: Remove stable manifold of a from F
        auto sm_it = std::find(data.stable_manifolds.begin(), data.stable_manifolds.end(), selected.sm);
        if (sm_it != data.stable_manifolds.end()) {
            data.stable_manifolds.erase(sm_it);
        }

        // Update gabriel_topology: a's boundary is removed
        for (const Edge &e : selected.sm.gabriel_edges) {
            data.gabriel_topology[e]--;
        }

        // Remove all of saddle a's own pairs from V (a is removed from F)
        // and update pairs for other saddles paired with b
        std::vector<valid_pair> to_remove;
        std::vector<valid_pair> to_add;

        for (const auto &vp : data.valid_pairs) {
            // Remove all pairs involving saddle a
            if (vp.sm == selected.sm) {
                to_remove.push_back(vp);
                continue;
            }
            // For each saddle s ≠ a, if (s, b) is valid, replace b with each other maximum c
            if (vp.v == b) {
                to_remove.push_back(vp);
                for (const auto &c : other_maxima) {
                    auto sc_it = std::find(data.valid_pairs.begin(), data.valid_pairs.end(), valid_pair(vp.sm, c, FT(0)));
                    if (sc_it != data.valid_pairs.end()) {
                        to_remove.push_back(*sc_it);
                    } else {
                        const auto &face_s = vp.sm.saddle_2.df.face;
                        std::vector face_s_vertices = {face_s.vertex(0)->point(), face_s.vertex(1)->point(), face_s.vertex(2)->point()};
                        Point cc_s = circumcenter()(face_s_vertices.begin(), face_s_vertices.end());
                        FT distance_saddle = squared_distance()(cc_s, face_s.vertex(0)->point());
                        const auto &face_a = selected.sm.saddle_2.df.face;
                        FT distance_voronoi = c.is_infinite ? FT(std::numeric_limits<double>::max()) : squared_distance()(c.point, face_a.vertex(0)->point());
                        FT distance_heuristic = distance_voronoi - distance_saddle;
                        to_add.push_back(valid_pair(vp.sm, c, CGAL::abs(distance_heuristic)));
                    }
                }
            }
        }

        for (const auto &r : to_remove) {
            auto it = std::find(data.valid_pairs.begin(), data.valid_pairs.end(), r);
            if (it != data.valid_pairs.end())
                data.valid_pairs.erase(it);
        }

        for (const auto &a : to_add) {
            data.valid_pairs.push_back(a);
        }

        std::sort(data.valid_pairs.begin(), data.valid_pairs.end());
    }
}