#include <CGAL/Linear_algebraCd.h>

#include "Delaunay.h"

#include "utils.h"
#include "polyscope/curve_network.h"
#include "polyscope/point_cloud.h"

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
            constexpr int skip_cols[4] = {0, 1, 2, 3};
            LA_Matrix minors[4] = {LA_Matrix(3, 3), LA_Matrix(3, 3), LA_Matrix(3, 3), LA_Matrix(3, 3)};

            for (int m = 0; m < 4; ++m) {
                int col_idx = 0;
                for (int k = 0; k < 4; ++k) {
                    if (k == skip_cols[m]) continue;
                    minors[m](0, col_idx) = v1[k];
                    minors[m](1, col_idx) = v2[k];
                    minors[m](2, col_idx) = v3[k];
                    col_idx++;
                }
            }

            FT n0 =  LA::determinant(minors[0]);
            FT n1 = -LA::determinant(minors[1]);
            FT n2 =  LA::determinant(minors[2]);
            FT n3 = -LA::determinant(minors[3]);

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
            Vector exact_direction(n0, n1, n2, n3);
            if (dot()(exact_direction, Vector(opposite_point) - p0) > FT(0)) {
                exact_direction = -exact_direction;
            }

            Point infinite_point = Vector(finite_circumcenter) + Point(99999999 * exact_direction[0], 99999999 * exact_direction[1], 99999999 * exact_direction[2], 99999999 * exact_direction[3]);

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
    return Voronoi_face{voronoi_vertices, voronoi_edges, ps_vertices, ps_edges, face};
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
    }
    if (!delaunay.is_valid()) {
        std::cerr << "Triangulation is invalid!" << std::endl;
    }
}

bool is_index_two_critical_point(const Face &face, const Delaunay &dt) {
    Voronoi_face f = delaunay_face_dual(face, dt);
    std::vector face_vertices = {face.face.vertex(0)->point(), face.face.vertex(1)->point(), face.face.vertex(2)->point()};
    Point cc = circumcenter()(face_vertices.begin(), face_vertices.end());
    const FT sq_radius = squared_distance()(cc, face.face.vertex(0)->point());

    std::set<Delaunay::Vertex_handle> neighboring_vertices;

    for (int i = 0; i < 3; ++i) {
        std::vector<Delaunay::Full_cell_handle> v_cells;
        dt.incident_full_cells(face.face.vertex(i), std::back_inserter(v_cells));
        for (const auto& cell : v_cells) {
            for (auto v = cell->vertices_begin(); v != cell->vertices_end(); ++v) {
                neighboring_vertices.insert(*v);
            }
        }
    }

    for (auto v : neighboring_vertices) {
        if (v == face.face.vertex(0) || v == face.face.vertex(1) || v == face.face.vertex(2)) continue;
        if (dt.is_infinite(v)) continue;

        FT sq_dist = squared_distance()(cc, v->point());
        if (sq_dist < sq_radius) {
            return false;
        }
    }

    Point pi = face.face.vertex(0)->point();
    Point pj = face.face.vertex(1)->point();
    Point pl = face.face.vertex(2)->point();

    return dot()(Vector(pj) - pi, Vector(pl) - pi) >= 0 && dot()(Vector(pi) - pj, Vector(pl) -pj ) >= 0 && dot()(Vector(pi) - pl, Vector(pj) - pl) >= 0;
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
