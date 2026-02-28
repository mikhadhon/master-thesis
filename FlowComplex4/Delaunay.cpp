#include <CGAL/Linear_algebraCd.h>

#include "Delaunay.h"

#include "utils.h"
#include "polyscope/curve_network.h"
#include "polyscope/point_cloud.h"
#include "polyscope/surface_mesh.h"

Delaunay::Full_cell_handle get_next_cell(const Delaunay::Full_cell_handle &current_cell, const Face &face, const Delaunay::Vertex_handle current_opposite, const Delaunay &dt) {
    for (int i = 0; i <= dt.maximal_dimension(); ++i) {
        const auto neighbor = current_cell->neighbor(i);

        if (auto vertex_at_i = current_cell->vertex(i); vertex_at_i != face.face.vertex(0) && vertex_at_i != face.face.vertex(1) && vertex_at_i != face.face.vertex(2) && vertex_at_i != current_opposite) {
            if (neighbor->has_vertex(face.face.vertex(0)) && neighbor->has_vertex(face.face.vertex(1)) && neighbor->has_vertex(face.face.vertex(2))) {
                return neighbor;
            }
        }
    }
    return {};
}

void get_shared_delaunay_tetrahedron(const Delaunay::Full_cell_handle cell1, const Delaunay::Full_cell_handle cell2, std::vector<Point> &tetrahedron_vertices) {
    for (int i = 0; i <= cell1->maximal_dimension(); ++i) {
        if (auto current_vertex = cell1->vertex(i); cell2->has_vertex(current_vertex)) {
            tetrahedron_vertices.push_back(current_vertex->point());
        }
    }
}

void get_incident_faces_to_cell(const Delaunay::Full_cell_handle &cell, const Delaunay &dt, std::unordered_set<Face, FaceHash> &faces) {
    if (!dt.is_infinite(*cell)) {
        for (int i = 0; i < dt.maximal_dimension(); ++i) {
            for (int j = i + 1; j <= dt.maximal_dimension(); ++j) {
                for (int k = j + 1; k <= dt.maximal_dimension(); ++k) {
                    std::vector face_vertices = {cell->vertex(i), cell->vertex(j), cell->vertex(k)};
                    std::ranges::sort(face_vertices);
                    Delaunay::Face dt_face(cell);
                    for (int l = 0; l < 3; ++l){
                        dt_face.set_index(l, cell->index(face_vertices[l]));
                    }
                    Face face(dt_face, -1, -1);
                    for (auto vertex = cell->vertices_begin(); vertex != cell->vertices_end(); ++vertex) {
                        if (*vertex != face.face.vertex(0) && *vertex != face.face.vertex(1) && *vertex != face.face.vertex(2)) {
                            if (face.index_of_covertex_0 == -1) face.index_of_covertex_0 = cell->index(*vertex);
                            else if (face.index_of_covertex_1 == -1) face.index_of_covertex_1 = cell->index(*vertex);
                        }
                    }
                    faces.insert(face);
                }
            }
        }
    }
}

std::unordered_set<Face, FaceHash> get_incident_faces_to_edge(const Edge &edge, const Delaunay &dt) {
    std::vector<Delaunay::Full_cell_handle> incident_cells;
    get_incident_cells_to_edge(edge, dt, incident_cells);

    std::unordered_set<Face, FaceHash> faces;
    for (auto cell : incident_cells) {
        if (dt.is_infinite(cell)) continue;
        for (int i = 0; i <= dt.maximal_dimension(); ++i) {
            auto current_vertex = cell->vertex(i);
            if (current_vertex != edge.vertex1 && current_vertex != edge.vertex2) {
                std::vector face_vertices = {current_vertex, edge.vertex1, edge.vertex2};
                std::sort(face_vertices.begin(), face_vertices.end());
                Delaunay::Face dt_face(cell);
                for (int l = 0; l < 3; ++l){
                    dt_face.set_index(l, cell->index(face_vertices[l]));
                }
                Face face(dt_face, -1, -1);
                for (auto vertex = cell->vertices_begin(); vertex != cell->vertices_end(); ++vertex) {
                    if (*vertex != face.face.vertex(0) && *vertex != face.face.vertex(1) && *vertex != face.face.vertex(2)) {
                        if (face.index_of_covertex_0 == -1) face.index_of_covertex_0 = cell->index(*vertex);
                        else if (face.index_of_covertex_1 == -1) face.index_of_covertex_1 = cell->index(*vertex);
                    }
                }
                faces.insert(face);
            }
        }

    }
    return faces;
}

std::vector<Voronoi_face> delaunay_edge_dual(const Edge &edge, const Face &df, const Delaunay &dt) {
    const std::unordered_set<Face, FaceHash> dual_faces = get_incident_faces_to_edge(edge, dt);
    std::vector<Voronoi_face> dual_voronoi_faces;
    for (const auto& dual_face : dual_faces) {
        if (dual_face == df) continue;
        dual_voronoi_faces.push_back(delaunay_face_dual(dual_face, dt));
    }
    return dual_voronoi_faces;
}

Voronoi_face delaunay_face_dual(const Face &face, const Delaunay &dt) {
    std::vector<Delaunay::Full_cell_handle> incident_cells;
    get_incident_cells_to_face(face, dt, incident_cells);

    std::vector<Delaunay::Full_cell_handle> ordered_incident_cells;
    std::vector<bool> is_infinite_cell;

    Delaunay::Full_cell_handle start_cell = face.face.full_cell();
    ordered_incident_cells.push_back(start_cell);
    is_infinite_cell.push_back(dt.is_infinite(start_cell));

    Delaunay::Vertex_handle current_opposite_vertex = start_cell->vertex(face.index_of_covertex_0);
    Delaunay::Full_cell_handle current_cell = start_cell;

    do {
        Delaunay::Full_cell_handle next_cell = get_next_cell(current_cell, face, current_opposite_vertex, dt);

        if (next_cell == Delaunay::Full_cell_handle() || next_cell == start_cell) {
            break;
        }

        Delaunay::Vertex_handle new_opposite = nullptr;
        for (auto vertex = next_cell->vertices_begin(); vertex != next_cell->vertices_end(); ++vertex) {
            if (*vertex != face.face.vertex(0) && *vertex != face.face.vertex(1) && *vertex != face.face.vertex(2) && !current_cell->has_vertex(*vertex)) {
                new_opposite = *vertex;
                break;
            }
        }

        ordered_incident_cells.push_back(next_cell);
        is_infinite_cell.push_back(dt.is_infinite(next_cell));

        current_cell = next_cell;
        current_opposite_vertex = new_opposite;
    } while (current_cell != start_cell && ordered_incident_cells.size() < incident_cells.size());

    std::vector<Voronoi_vertex> voronoi_vertices;
    std::vector<Voronoi_edge> voronoi_edges;
    std::vector<Eigen::VectorXd> ps_vertices;
    std::vector<std::array<size_t, 2>> ps_edges;
    std::vector<Point> infinite_vertices;

    for (size_t i = 0; i < ordered_incident_cells.size(); ++i) {
        Voronoi_vertex voronoi_vertex;

        if (!is_infinite_cell[i]) {
            Point cell_circumcenter = circumcenter()(boost::make_transform_iterator(ordered_incident_cells[i]->vertices_begin(), Vertex_to_point()), boost::make_transform_iterator(ordered_incident_cells[i]->vertices_end(), Vertex_to_point()));
            voronoi_vertex = Voronoi_vertex{false, cell_circumcenter, Vector()};
        }
        else {
            Delaunay::Full_cell_handle finite_neighbor = is_infinite_cell[(i + 1) % is_infinite_cell.size()] ? ordered_incident_cells[(i - 1) % is_infinite_cell.size()] : ordered_incident_cells[(i + 1) % is_infinite_cell.size()];
            Point finite_circumcenter = circumcenter()(boost::make_transform_iterator(finite_neighbor->vertices_begin(), Vertex_to_point()), boost::make_transform_iterator(finite_neighbor->vertices_end(), Vertex_to_point()));

            std::vector<Point> shared_tetrahedron_points;
            get_shared_delaunay_tetrahedron(ordered_incident_cells[i], finite_neighbor, shared_tetrahedron_points);

            // Compute exact tetrahedron normal via cofactor expansion of the 3x4 edge matrix
            Point p0 = shared_tetrahedron_points[0];
            FT v1[4], v2[4], v3[4];
            for (int k = 0; k < 4; ++k) {
                v1[k] = shared_tetrahedron_points[1][k] - p0[k];
                v2[k] = shared_tetrahedron_points[2][k] - p0[k];
                v3[k] = shared_tetrahedron_points[3][k] - p0[k];
            }

            // 4D generalized cross product: n[k] = (-1)^k * det of 3x3 minor deleting column k
            auto det3 = [](FT a00, FT a01, FT a02,
                           FT a10, FT a11, FT a12,
                           FT a20, FT a21, FT a22) -> FT {
                return a00 * (a11 * a22 - a12 * a21)
                     - a01 * (a10 * a22 - a12 * a20)
                     + a02 * (a10 * a21 - a11 * a20);
            };

            FT n0 =  det3(v1[1],v1[2],v1[3], v2[1],v2[2],v2[3], v3[1],v3[2],v3[3]);
            FT n1 = -det3(v1[0],v1[2],v1[3], v2[0],v2[2],v2[3], v3[0],v3[2],v3[3]);
            FT n2 =  det3(v1[0],v1[1],v1[3], v2[0],v2[1],v2[3], v3[0],v3[1],v3[3]);
            FT n3 = -det3(v1[0],v1[1],v1[2], v2[0],v2[1],v2[2], v3[0],v3[1],v3[2]);

            // Orient normal to point away from the opposite vertex in the finite neighbor
            Point opposite_point;
            for (auto v = finite_neighbor->vertices_begin(); v != finite_neighbor->vertices_end(); ++v) {
                if (!dt.is_infinite(*v) && !ordered_incident_cells[i]->has_vertex(*v)) {
                    opposite_point = (*v)->point();
                    break;
                }
            }
            FT to_opp_dot = n0 * (opposite_point[0] - p0[0])
                          + n1 * (opposite_point[1] - p0[1])
                          + n2 * (opposite_point[2] - p0[2])
                          + n3 * (opposite_point[3] - p0[3]);
            if (to_opp_dot > FT(0)) {
                n0 = -n0; n1 = -n1; n2 = -n2; n3 = -n3;
            }

            Vector exact_direction(n0, n1, n2, n3);

            // Keep approximate point for Polyscope visualization
            double nd0 = CGAL::to_double(n0), nd1 = CGAL::to_double(n1),
                   nd2 = CGAL::to_double(n2), nd3 = CGAL::to_double(n3);
            double len = std::sqrt(nd0*nd0 + nd1*nd1 + nd2*nd2 + nd3*nd3);
            if (len > 0) { nd0 /= len; nd1 /= len; nd2 /= len; nd3 /= len; }
            double scale = 2.0;
            Point infinite_point = Vector(finite_circumcenter) + Point(nd0*scale, nd1*scale, nd2*scale, nd3*scale);

            voronoi_vertex = Voronoi_vertex{true, infinite_point, exact_direction};
            infinite_vertices.push_back(infinite_point);
        }

        voronoi_vertices.push_back(voronoi_vertex);
        ps_vertices.push_back(make_point_eigen(voronoi_vertex.point));
    }

    for (size_t i = 0; i < ordered_incident_cells.size(); ++i) {
        size_t next_i = (i + 1) % ordered_incident_cells.size();

        if (!(is_infinite_cell[i] && is_infinite_cell[next_i])) {
            voronoi_edges.push_back(Voronoi_edge(
               voronoi_vertices[i],
               voronoi_vertices[next_i],
               ordered_incident_cells[i],
               ordered_incident_cells[next_i]
           ));
            ps_edges.push_back({i, next_i});
        }
    }
    //polyscope::registerPointCloud("infinite voronoi vertices", infinite_vertices);
    return Voronoi_face{voronoi_vertices, voronoi_edges, ps_vertices, ps_edges, face};
}

// Face voronoi_edge_dual(Voronoi_edge &voronoi_edge) {
//
// }

Eigen::VectorXd simplex_circumsphere(const Delaunay::Full_cell_handle &simplex) {
    std::vector<Point> simplex_points;
    for (auto point = simplex->vertices_begin(); point != simplex->vertices_end(); ++point) {
        simplex_points.push_back((*point)->point());
    }
    const Point cgal_center = circumcenter()(simplex_points.begin(), simplex_points.end());
    return make_point_eigen(cgal_center);
}

void insert_points(const std::vector<Point> &points, Delaunay &delaunay) {
    Delaunay::Vertex_handle hint;
    for (auto & point : points) {
        if (Delaunay::Vertex_handle() != hint) {
            hint = delaunay.insert(point, hint);
        }
        else {
            hint = delaunay.insert(point);
        }
        // printf("Processing: %d/%d\n", ++i, static_cast<int>(points.size()));
    }
    if (!delaunay.is_valid()) {
        std::cerr << "Triangulation is invalid!" << std::endl;
    }
}

bool intersect_df_vp_d(const Face &df, Voronoi_face &vp, Point &out) {
    Point pi = df.face.vertex(0)->point();
    Point pj = df.face.vertex(1)->point();
    Point pl = df.face.vertex(2)->point();

    const int nv = static_cast<int>(vp.voronoi_vertices.size());
    //
    // // Find three finite basis vertices for the polygon plane
    int b0 = -1, b1 = -1, b2 = -1;
    for (int i = 0; i < nv; ++i) {
        if (b0 == -1) b0 = i;
        else if (b1 == -1) b1 = i;
        else if (b2 == -1) { b2 = i; break; }
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
    if (!(s >= FT(0) && t >= FT(0) && s + t <= FT(1))) {
        return false;
    }

    // std::vector df_vertices = {df.face.vertex(0)->point(), df.face.vertex(1)->point(), df.face.vertex(2)->point()};
    // Point p_int = circumcenter()(df_vertices.begin(), df_vertices.end());

    if (!(dot()(Vector(pj) - pi, Vector(pl) - pi) >= 0 && dot()(Vector(pi) - pj, Vector(pl) -pj ) >= 0 && dot()(Vector(pi) - pl, Vector(pj) - pl) >= 0)) {
        return false;
    };

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


    Point n1(VT.row(0).begin(), VT.row(0).end());
    Point n2(VT.row(1).begin(), VT.row(1).end());

    Vector pjpi = (Vector(pj) - pi);
    Vector spjpi = Vector(s * pjpi[0], s * pjpi[1], s * pjpi[2], s * pjpi[3]);
    Vector plpi = (Vector(pl) - pi);
    Vector tpjpi = Vector(t * plpi[0], t * plpi[1], t * plpi[2], t * plpi[3]);
    Vector p_int = pi + spjpi + tpjpi;

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
            vec = edge.vertex1.infinite_direction;
        }
        else continue;
        // std::vector simplex_vertices  = {Point(p_int), ref, Point(Vector(ref) + vec), Point(Vector(ref) + pjpi), Point(Vector(ref) + plpi)};
        //
        // auto o = orientation()(simplex_vertices.begin(), simplex_vertices.end());

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

        if (o_point != o_centroid && o_centroid != CGAL::ZERO) {
            return false;
        }
        // int current_sign = (o == CGAL::POSITIVE) ? 1 : -1;
        // if (sign == 0) sign = current_sign;
        // else if (sign != current_sign) {
        //     return false;
        // }
    }

    // Vector pjpi = (Vector(pj) - pi);
    // Vector spjpi = Vector(s * pjpi[0], s * pjpi[1], s * pjpi[2], s * pjpi[3]);
    // Vector plpi = (Vector(pl) - pi);
    // Vector tpjpi = Vector(t * plpi[0], t * plpi[1], t * plpi[2], t * plpi[3]);
    out = p_int;
    return true;
}


bool is_index_two_critical_point(const Face &face, const Delaunay &dt) {
    Voronoi_face f = delaunay_face_dual(face, dt);
    // for (auto &edge : f.voronoi_edges) {
    //     if (edge.vertex1.is_infinite || edge.vertex2.is_infinite) {
    //         return false;
    //     }
    // }
    Point out;
    bool test2 = intersect_df_vp_d(face, f, out);
    std::vector face_vertices = {face.face.vertex(0)->point(), face.face.vertex(1)->point(), face.face.vertex(2)->point()};
    Point cc = circumcenter()(face_vertices.begin(), face_vertices.end());
    const FT sq_radius = squared_distance()(cc, face.face.vertex(0)->point());

    std::vector<Delaunay::Full_cell_handle> cell_neighbors;
    get_incident_cells_to_face(face, dt, cell_neighbors);
    std::set<Delaunay::Vertex_handle> neighboring_vertices;

    for (const auto cell : cell_neighbors) {
        for (auto v = cell->vertices_begin(); v != cell->vertices_end(); ++v) {
            neighboring_vertices.insert(*v);
        }
    }

    bool test = true;
    for (auto v : neighboring_vertices) {
        if (v == face.face.vertex(0) || v == face.face.vertex(1) || v == face.face.vertex(2)) continue;
        if (dt.is_infinite(v)) continue;

        FT sq_dist = squared_distance()(cc, v->point());
        if (sq_dist < sq_radius) {
            test = false;
        }
    }

    Point pi = face.face.vertex(0)->point();
    Point pj = face.face.vertex(1)->point();
    Point pl = face.face.vertex(2)->point();

    bool result = test && dot()(Vector(pj) - pi, Vector(pl) - pi) >= 0 && dot()(Vector(pi) - pj, Vector(pl) -pj ) >= 0 && dot()(Vector(pi) - pl, Vector(pj) - pl) >= 0;
    if (result && !test2) {
        std::cerr << "Test failed!" << std::endl;
    }
    return result;
}

bool is_gabriel(const Edge &edge, const Delaunay &dt) {
    Point p1 = edge.vertex1->point();
    Point p2 = edge.vertex2->point();

    Point edge_midpoint = midpoint()(p1, p2);
    const FT sq_radius = squared_distance()(p1, edge_midpoint);

    std::vector<Delaunay::Full_cell_handle> cell_neighbors;
    get_incident_cells_to_edge(edge, dt, cell_neighbors);
    std::set<Delaunay::Vertex_handle> neighboring_vertices;

    for (const auto cell : cell_neighbors) {
        for (auto v = cell->vertices_begin(); v != cell->vertices_end(); ++v) {
            neighboring_vertices.insert(*v);
        }
    }

    for (auto v : neighboring_vertices) {
        if (v == edge.vertex1 || v == edge.vertex2) continue;
        if (dt.is_infinite(v)) continue;

        FT sq_dist = squared_distance()(edge_midpoint, v->point());
        if (sq_dist < sq_radius) {
            return false;
        }
    }

    return true;
}

void get_incident_cells_to_edge(const Edge &edge, const Delaunay &dt, std::vector<Delaunay::Full_cell_handle> &cell_neighbors) {
    cell_neighbors.clear();

    std::vector<Delaunay::Full_cell_handle> v1_cells;
    dt.incident_full_cells(edge.vertex1, std::back_inserter(v1_cells));

    for (auto cell: v1_cells) {
        if (cell->has_vertex(edge.vertex2)) {
            cell_neighbors.emplace_back(cell);
        }
    }
}

void get_incident_cells_to_face(const Face &face, const Delaunay &dt, std::vector<Delaunay::Full_cell_handle> &cell_neighbors) {
    cell_neighbors.clear();

    std::vector<Delaunay::Full_cell_handle> v0_cells;
    dt.incident_full_cells(face.face.vertex(0), std::back_inserter(v0_cells));

    for (auto cell: v0_cells) {
        if (cell->has_vertex(face.face.vertex(1)) && cell->has_vertex(face.face.vertex(2))) {
            cell_neighbors.emplace_back(cell);
        }
    }
}

void get_tetrahedron_normal(const std::vector<Eigen::VectorXd> &tetrahedron_points, Eigen::VectorXd &normal) {
    const Eigen::VectorXd v1 = tetrahedron_points[1] - tetrahedron_points[0];
    const Eigen::VectorXd v2 = tetrahedron_points[2] - tetrahedron_points[0];
    const Eigen::VectorXd v3 = tetrahedron_points[3] - tetrahedron_points[0];

    Eigen::MatrixXd A(4, 3);
    A.col(0) = v1;
    A.col(1) = v2;
    A.col(2) = v3;

    const Eigen::JacobiSVD<Eigen::MatrixXd> svd(A, Eigen::ComputeFullU);
    normal = svd.matrixU().col(3);
    normal = normal.normalized();
}

void get_facet_vertices(const Delaunay &delaunay, const Delaunay::Facet_iterator &facet,
    std::vector<Delaunay::Vertex_handle> &facet_vertices
    ) {
    const int dimension = facet->full_cell()->maximal_dimension();
    const int co_vertex_index = facet->index_of_covertex();

    for (int i = 0; i <= dimension; ++i) {
        if (i != co_vertex_index) {
            if (const auto vertex_handle = facet->full_cell()->vertex(i);
                vertex_handle != nullptr && !delaunay.is_infinite(vertex_handle)
            ) {
                facet_vertices.push_back(vertex_handle);
            }
        }
    }
}

void find_gabriel_edges(Delaunay &dt, std::vector<Edge> &gabriel_edges) {
    std::unordered_set<Edge, EdgeHash> visited_edges;
    const int dim = dt.maximal_dimension();

    for (auto cell = dt.finite_full_cells_begin(); cell != dt.finite_full_cells_end(); ++cell) {
        for (int i = 0; i < dim; ++i) {
            for (int j = i + 1; j <= dim; ++j) {
                Delaunay::Vertex_handle v1 = cell->vertex(i);
                Delaunay::Vertex_handle v2 = cell->vertex(j);

                if (v2 < v1) std::swap(v1, v2);
                Edge e(v1, v2, Delaunay::Vertex_handle());
                if (visited_edges.insert(e).second && is_gabriel(e, dt)) {
                    gabriel_edges.push_back(e);
                }
            }
        }
    }
}

std::unordered_set<Face, FaceHash> get_delaunay_faces(Delaunay &dt) {
    std::unordered_set<Face, FaceHash> faces;
    for (auto cell = dt.full_cells_begin(); cell != dt.full_cells_end(); ++cell) {
        get_incident_faces_to_cell(cell, dt, faces);
    }
    return faces;
}
