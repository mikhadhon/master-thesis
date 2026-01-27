#include <CGAL/Linear_algebraCd.h>

#include "Delaunay.h"

#include "utils.h"
#include "polyscope/curve_network.h"
#include "polyscope/point_cloud.h"
#include "polyscope/surface_mesh.h"

typedef CGAL::Linear_algebraCd<FT> Linear_algebra;
typedef Linear_algebra::Matrix Matrix;
typedef Linear_algebra::Vector Vector;
typedef Linear_algebra::RT RT;

// Voronoi_face delaunay_edge_dual(Edge &edge, Face &df, Delaunay &dt) {
//     std::vector<Delaunay::Full_cell_handle> incident_cells;
//     std::vector<Delaunay::Vertex_handle> current_face_vertices;
//
//     incident_cells.push_back(df.face.full_cell());
//     incident_cells.push_back(df.face.full_cell()->neighbor(df.index_of_covertex));
//
//     current_face_vertices.push_back(df.face.vertex(0));
//     current_face_vertices.push_back(df.face.vertex(1));
//     current_face_vertices.push_back(df.face.vertex(2));
//
//     Delaunay::Full_cell_handle start_cell;
//     Delaunay::Full_cell_handle current_cell;
//     Delaunay::Full_cell_handle next_cell;
//
//     if (dt.is_infinite(incident_cells[0]) && dt.is_infinite(incident_cells[1])) {
//         Eigen::VectorXd voronoi_vertex_start(dt.maximal_dimension());
//         Eigen::VectorXd voronoi_vertex_end(dt.maximal_dimension());
//         std::vector<Eigen::VectorXd> facet_points = {make_point_eigen(current_face_vertices[0]->point()), make_point_eigen(current_face_vertices[1]->point()), make_point_eigen(current_face_vertices[2]->point())};
//         Eigen::VectorXd normal;
//         get_facet_normal(facet_points, normal);
//         simplex_circumsphere(df.face.full_cell(), voronoi_vertex_start);
//         orient_voronoi_edge(facet_points, voronoi_vertex_start, normal);
//         voronoi_vertex_start = voronoi_vertex_start + 2*normal;
//         voronoi_vertex_end = voronoi_vertex_start - 4*normal;
//         Voronoi_vertex voronoi_vertex1(true, voronoi_vertex_start);
//         Voronoi_vertex voronoi_vertex2(true, voronoi_vertex_end);
//         std::vector<Voronoi_vertex> voronoi_vertices = {voronoi_vertex1, voronoi_vertex2};
//         std::vector<Eigen::VectorXd> ps_vertices = {voronoi_vertex1.point, voronoi_vertex2.point};
//         std::vector<std::array<size_t, 2>> ps_edges = {{0, 1}};
//         std::vector<Voronoi_edge> edges = {{Voronoi_edge(voronoi_vertex1, voronoi_vertex2, df.face.full_cell(), df.face.full_cell())}};
//         return Voronoi_face(voronoi_vertices, edges, ps_vertices, ps_edges);
//     }
//     else if (dt.is_infinite(incident_cells[0])) {
//         start_cell = incident_cells[0];
//         current_cell = incident_cells[0];
//         next_cell = incident_cells[1];
//     }
//     else {
//         start_cell = incident_cells[1];
//         current_cell = incident_cells[1];
//         next_cell = incident_cells[0];
//     }
//
//     Delaunay::Vertex_handle co_vertex_current_next;
//     for (auto vertex = start_cell->vertices_begin(); vertex != start_cell->vertices_end(); ++vertex) {
//         if (*vertex != current_face_vertices[0] && *vertex != current_face_vertices[1] && *vertex != current_face_vertices[2]) {
//             co_vertex_current_next = *vertex;
//         }
//     }
//
//     bool test = dt.is_infinite(co_vertex_current_next);
//
//     std::vector<Voronoi_vertex> voronoi_vertices;
//     std::vector<Voronoi_edge> voronoi_edges;
//     std::vector<Eigen::VectorXd> ps_vertices;
//     std::vector<std::array<size_t, 2>> ps_edges;
//     Eigen::VectorXd voronoi_vertex_start(dt.maximal_dimension());
//
//     if (!dt.is_infinite(co_vertex_current_next)) {
//         simplex_circumsphere(start_cell, voronoi_vertex_start);
//         Voronoi_vertex voronoi_vertex(false, voronoi_vertex_start);
//         voronoi_vertices.push_back(voronoi_vertex);
//         ps_vertices.push_back(voronoi_vertex.point);
//     }
//     else {
//         std::vector<Eigen::VectorXd> facet_points = {make_point_eigen(current_face_vertices[0]->point()), make_point_eigen(current_face_vertices[1]->point()), make_point_eigen(current_face_vertices[2]->point())};
//         Eigen::VectorXd normal;
//         get_facet_normal(facet_points, normal);
//         simplex_circumsphere(next_cell, voronoi_vertex_start);
//         orient_voronoi_edge(facet_points, voronoi_vertex_start, normal);
//         voronoi_vertex_start = voronoi_vertex_start + 2*normal;
//         Voronoi_vertex voronoi_vertex(true, voronoi_vertex_start);
//         voronoi_vertices.push_back(voronoi_vertex);
//         ps_vertices.push_back(voronoi_vertex.point);
//     }
//
//     while (next_cell != start_cell && !dt.is_infinite(next_cell)) {
//         Eigen::VectorXd voronoi_vertex_next(3);
//
//         simplex_circumsphere(next_cell, voronoi_vertex_next);
//         Voronoi_vertex voronoi_vertex(false, voronoi_vertex_next);
//         voronoi_vertices.push_back(voronoi_vertex);
//         ps_vertices.push_back(voronoi_vertex.point);
//         voronoi_edges.push_back(
//             Voronoi_edge(
//                 voronoi_vertices[voronoi_vertices.size() - 2],
//                 voronoi_vertices[voronoi_vertices.size() - 1],
//                 current_cell,
//                 next_cell
//             )
//         );
//         ps_edges.push_back({voronoi_vertices.size() - 2, voronoi_vertices.size() - 1});
//
//         for (int i = 0; i < 3; ++i) {
//             if (current_face_vertices[i] != edge.vertex1 && current_face_vertices[i] != edge.vertex2) {
//                 Delaunay::Vertex_handle temp = current_cell->mirror_vertex(current_cell->index(co_vertex_current_next), dt.maximal_dimension());
//                 co_vertex_current_next = current_face_vertices[i];
//                 current_face_vertices[i] = temp;
//                 current_cell = next_cell;
//                 next_cell = next_cell->neighbor(next_cell->index(co_vertex_current_next));
//                 break;
//             }
//         }
//     }
//
//     if (!dt.is_infinite(next_cell)) {
//         voronoi_edges.push_back(
//             Voronoi_edge(
//                 voronoi_vertices[voronoi_vertices.size() - 1],
//                 voronoi_vertices[0],
//                 current_cell,
//                 start_cell
//             )
//         );
//         ps_edges.push_back({voronoi_vertices.size() - 1, 0});
//     }
//     else {
//         Eigen::VectorXd voronoi_vertex_next(3);
//
//         simplex_circumsphere(current_cell, voronoi_vertex_next);
//         std::vector<Eigen::VectorXd> final_facet_vertices;
//         for (auto vertex = current_cell->vertices_begin(); vertex != current_cell->vertices_end(); ++vertex) {
//             if (*vertex != co_vertex_current_next) {
//                 final_facet_vertices.push_back(make_point_eigen((*vertex)->point()));
//             }
//         }
//         Eigen::VectorXd normal;
//         get_facet_normal(final_facet_vertices, normal);
//         orient_voronoi_edge(final_facet_vertices, voronoi_vertex_next, normal);
//         Voronoi_vertex voronoi_vertex(true, voronoi_vertex_next + 2*normal);
//         voronoi_vertices.push_back(voronoi_vertex);
//         ps_vertices.push_back(voronoi_vertex.point);
//         voronoi_edges.push_back(
//             Voronoi_edge(
//                 voronoi_vertices[voronoi_vertices.size() - 2],
//                 voronoi_vertices[voronoi_vertices.size() - 1],
//                 current_cell,
//                 current_cell->neighbor(current_cell->index(co_vertex_current_next))
//             )
//         );
//         ps_edges.push_back({voronoi_vertices.size() - 2, voronoi_vertices.size() - 1});
//     }
//     std::vector<Eigen::VectorXd> triangle = {make_point_eigen(current_face_vertices[0]->point()), make_point_eigen(current_face_vertices[1]->point()), make_point_eigen(current_face_vertices[2]->point())};
//     std::vector<std::array<size_t, 3>> triangle_face = {{0, 1, 2}};
//     //polyscope::registerSurfaceMesh("current Delaunay face", triangle, triangle_face);
//     std::vector<Eigen::VectorXd> voronoi_points;
//     for (auto & vertex : voronoi_vertices) {
//         voronoi_points.push_back(vertex.point);
//     }
//     std::vector<std::array<size_t, 2>> voronoi_edges_indices;
//     for (auto & edge : voronoi_edges) {
//         size_t index_first = std::find(voronoi_vertices.begin(), voronoi_vertices.end(), edge.vertex1) - voronoi_vertices.begin();
//         size_t index_second = std::find(voronoi_vertices.begin(), voronoi_vertices.end(), edge.vertex2) - voronoi_vertices.begin();
//         voronoi_edges_indices.push_back({index_first, index_second});
//     }
//     voronoi_edges_indices.push_back({voronoi_edges.size() - 1, 0});
//     //polyscope::registerCurveNetwork("current Voronoi face", voronoi_points, voronoi_edges_indices);
//     // polyscope::show();
//     return Voronoi_face(voronoi_vertices, voronoi_edges, ps_vertices, ps_edges);
// }

Delaunay::Full_cell_handle get_next_tetrahedron(Delaunay::Full_cell_handle &current_tet, Edge &edge, Delaunay::Vertex_handle current_opposite, Delaunay &dt) {
    for (int i = 0; i <= dt.maximal_dimension(); ++i) {
        auto neighbor = current_tet->neighbor(i);
        auto vertex_at_i = current_tet->vertex(i);

        if (vertex_at_i != edge.vertex1 && vertex_at_i != edge.vertex2 && vertex_at_i != current_opposite) {
            if (neighbor->has_vertex(edge.vertex1) && neighbor->has_vertex(edge.vertex2)) {
                return neighbor;
            }
        }
    }
    return Delaunay::Full_cell_handle();
}

Eigen::VectorXd calculate_infinite_voronoi_vertex_3d(
    Delaunay::Full_cell_handle finite_tetrahedron,
    const std::vector<Eigen::VectorXd> &shared_triangle_vertices,
    const Delaunay &dt) {

    // Berechne Circumcenter des finiten Tetraeders
    Eigen::VectorXd circumcenter;
    simplex_circumsphere(finite_tetrahedron, circumcenter);

    // Berechne Normale des gemeinsamen Dreiecks
    Eigen::VectorXd normal;
    get_facet_normal(const_cast<std::vector<Eigen::VectorXd>&>(shared_triangle_vertices), normal);

    // Berechne Schwerpunkt des Dreiecks
    Eigen::VectorXd triangle_center = Eigen::VectorXd::Zero(3);
    for (const auto &vertex : shared_triangle_vertices) {
        triangle_center += vertex;
    }
    triangle_center /= shared_triangle_vertices.size();

    // Orientiere Normale weg vom Circumcenter
    Eigen::VectorXd to_triangle = triangle_center - circumcenter;
    if (normal.dot(to_triangle) < 0) {
        normal = -normal;
    }

    // Berechne infinite vertex mit großer aber finiter Distanz
    return circumcenter + 2 * normal;
}

Voronoi_face delaunay_edge_dual(Edge &edge, Face &df, Delaunay &dt) {
    std::vector<Delaunay::Full_cell_handle> incident_cells_to_edge;
    get_incident_cells_to_vertices(edge, dt, incident_cells_to_edge);

    std::vector<Delaunay::Full_cell_handle> ordered_incident_cells;
    std::vector<bool> is_infinite_tet;

    Delaunay::Full_cell_handle start_tet = df.face.full_cell();
    ordered_incident_cells.push_back(start_tet);
    is_infinite_tet.push_back(dt.is_infinite(start_tet));

    Delaunay::Vertex_handle current_opposite_vertex = edge.co_vertex;
    Delaunay::Full_cell_handle current_tet = start_tet;

    do {
        Delaunay::Full_cell_handle next_tet = get_next_tetrahedron(current_tet, edge, current_opposite_vertex, dt);

        if (next_tet == Delaunay::Full_cell_handle() || next_tet == start_tet) {
            break;
        }

        Delaunay::Vertex_handle new_opposite = nullptr;
        for (auto vertex = next_tet->vertices_begin(); vertex != next_tet->vertices_end(); ++vertex) {
            if (*vertex != edge.vertex1 && *vertex != edge.vertex2 && !current_tet->has_vertex(*vertex)) {
                new_opposite = *vertex;
                break;
            }
        }

        ordered_incident_cells.push_back(next_tet);
        is_infinite_tet.push_back(dt.is_infinite(next_tet));

        current_tet = next_tet;
        current_opposite_vertex = new_opposite;
    } while (current_tet != start_tet && ordered_incident_cells.size() < incident_cells_to_edge.size());

    // if (std::count(is_infinite_tet.begin(), is_infinite_tet.end(), true) == 2) {
    //     std::cout << "[ ";
    //     for (bool i : is_infinite_tet) {
    //        std::cout << i << ",";
    //     }
    //     std::cout << "]" << std::endl;
    // }
    
    std::vector<Voronoi_vertex> voronoi_vertices;
    std::vector<Voronoi_edge> voronoi_edges;
    std::vector<Eigen::VectorXd> ps_vertices;
    std::vector<std::array<size_t, 2>> ps_edges;
    std::vector<Eigen::VectorXd> infinite_vertices;

    std::vector<Eigen::VectorXd> edge_triangle_vertices = {
        make_point_eigen(edge.vertex1->point()),
        make_point_eigen(edge.vertex2->point()),
        make_point_eigen(edge.co_vertex->point())
    };

    for (size_t i = 0; i < ordered_incident_cells.size(); ++i) {
        Voronoi_vertex v_vertex;

        if (!is_infinite_tet[i]) {
            Eigen::VectorXd circumcenter;
            simplex_circumsphere(ordered_incident_cells[i], circumcenter);
            v_vertex = Voronoi_vertex{false, circumcenter};
        } else {
            Delaunay::Full_cell_handle finite_neighbor = nullptr;

            for (size_t j = 0; j < ordered_incident_cells.size(); ++j) {
                if (j != i && !is_infinite_tet[j]) {
                    bool shares_facet = false;
                    for (int k = 0; k <= dt.maximal_dimension(); ++k) {
                        if (ordered_incident_cells[i]->neighbor(k) == ordered_incident_cells[j]) {
                            finite_neighbor = ordered_incident_cells[j];
                            shares_facet = true;
                            break;
                        }
                    }
                    if (shares_facet) break;
                }
            }

            if (finite_neighbor != nullptr) {
                Eigen::VectorXd infinite_point = calculate_infinite_voronoi_vertex_3d(
                    finite_neighbor, edge_triangle_vertices, dt);
                v_vertex = Voronoi_vertex{true, infinite_point};
                infinite_vertices.push_back(infinite_point);
            }
        }

        voronoi_vertices.push_back(v_vertex);
        ps_vertices.push_back(v_vertex.point);
    }

    // 4. Erstelle Voronoi-Edges zwischen aufeinanderfolgenden Vertices
    for (size_t i = 0; i < ordered_incident_cells.size(); ++i) {
        size_t next_i = (i + 1) % ordered_incident_cells.size();

        if (!(is_infinite_tet[i] && is_infinite_tet[next_i])) {
            voronoi_edges.push_back(Voronoi_edge(
               voronoi_vertices[i],
               voronoi_vertices[next_i],
               ordered_incident_cells[i],
               ordered_incident_cells[next_i]
           ));
            ps_edges.push_back({i, next_i});
        }
    }
    return Voronoi_face{voronoi_vertices, voronoi_edges, ps_vertices, ps_edges};

}

Voronoi_edge delaunay_face_dual(Face &face, Delaunay &dt) {
    Delaunay::Full_cell_handle incident_cell0, incident_cell1;
    incident_cell0 = face.face.full_cell();
    incident_cell1 = face.face.full_cell()->neighbor(face.index_of_covertex);

    bool infinite0 = dt.is_infinite(incident_cell0);
    bool infinite1 = dt.is_infinite(incident_cell1);

    Eigen::VectorXd voronoi_vertex0, voronoi_vertex1;
    if (!infinite0 && !infinite1) {
        simplex_circumsphere(incident_cell0, voronoi_vertex0);
        simplex_circumsphere(incident_cell1, voronoi_vertex1);
        return Voronoi_edge(Voronoi_vertex(infinite0, voronoi_vertex0), Voronoi_vertex(infinite1, voronoi_vertex1), incident_cell0, incident_cell1);
    }
    else {
        std::vector<Eigen::VectorXd> final_facet_vertices;

        final_facet_vertices.push_back(make_point_eigen(face.face.vertex(0)->point()));
        final_facet_vertices.push_back(make_point_eigen(face.face.vertex(1)->point()));
        final_facet_vertices.push_back(make_point_eigen(face.face.vertex(2)->point()));

        Delaunay::Full_cell_handle finite_cell = infinite0 ? incident_cell1 : incident_cell0;
        Delaunay::Full_cell_handle infinite_cell = infinite0 ? incident_cell0 : incident_cell1;

        Eigen::VectorXd normal;
        get_facet_normal(final_facet_vertices, normal);
        Eigen::VectorXd voronoi_vertex_next(3);
        simplex_circumsphere(finite_cell, voronoi_vertex_next);
        orient_voronoi_edge(final_facet_vertices, voronoi_vertex_next, normal);
        Voronoi_vertex infinite_voronoi_vertex(true, voronoi_vertex_next + 2*normal);
        Voronoi_vertex finite_voronoi_vertex(false, voronoi_vertex_next);

        return Voronoi_edge(finite_voronoi_vertex, infinite_voronoi_vertex, finite_cell, infinite_cell);
    }
}

Face voronoi_edge_dual(Voronoi_edge &voronoi_edge) {
    Delaunay::Face face_no_covertex(voronoi_edge.cell1);
    int index_of_covertex;
    int vertex_count = 0;

    for (auto vertex = voronoi_edge.cell1->vertices_begin(); vertex != voronoi_edge.cell1->vertices_end(); ++vertex) {
        int vertex_index_in_cell1 = voronoi_edge.cell1->index(*vertex);
        if (voronoi_edge.cell1->neighbor(vertex_index_in_cell1) == voronoi_edge.cell2) index_of_covertex = vertex_index_in_cell1;
        else {
            face_no_covertex.set_index(vertex_count++, vertex_index_in_cell1);
        }
    }

    Face face(face_no_covertex, index_of_covertex);
    return face;
}

void simplex_circumsphere(Delaunay::Full_cell_handle simplex, Eigen::VectorXd &center) {
    std::vector<Point> simplex_points;
    for (auto point = simplex->vertices_begin(); point != simplex->vertices_end(); ++point) {
        simplex_points.push_back((*point)->point());
    }
    Point cgal_center = circumcenter()(simplex_points.begin(), simplex_points.end());
    center = make_point_eigen(cgal_center);
}

void insert_points(std::vector<Point> &points, Delaunay &delaunay) {
    Delaunay::Vertex_handle hint;
    int i = 0;
    for (auto it = points.begin(); it != points.end(); ++it) {
        if (Delaunay::Vertex_handle() != hint) {
            hint = delaunay.insert(*it, hint);
        }
        else {
            hint = delaunay.insert(*it);
        }
        printf("Processing: %d/%d\n", ++i, static_cast<int>(points.size()));
    }
    if (!delaunay.is_valid()) {
        std::cerr << "Triangulation is invalid!" << std::endl;
    }
}

bool is_index_two_critical_point(Face &face, Delaunay &dt) {
    Voronoi_edge face_dual = delaunay_face_dual(face, dt);

    Eigen::VectorXd face_normal;

    std::vector<Eigen::VectorXd> face_vertices = {make_point_eigen(face.face.vertex(0)->point()), make_point_eigen(face.face.vertex(1)->point()), make_point_eigen(face.face.vertex(2)->point())};
    get_facet_normal(face_vertices, face_normal);

    double vertex0_direction = face_normal.dot(face_dual.vertex1.point - face_vertices[0]);
    double vertex1_direction = face_normal.dot(face_dual.vertex2.point - face_vertices[1]);

    Eigen::VectorXd i = make_point_eigen(face.face.vertex(0)->point());
    Eigen::VectorXd j = make_point_eigen(face.face.vertex(1)->point());
    Eigen::VectorXd l = make_point_eigen(face.face.vertex(2)->point());

    return (vertex0_direction * vertex1_direction < 0) && (j-i).dot(l-i) > 0 && (i-j).dot(l-j) > 0 && (i-l).dot(j-l) > 0;
}

bool is_gabriel(Edge &edge, Delaunay &dt) {
    std::vector<Delaunay::Point> triangle_points;
    triangle_points.push_back(edge.vertex1->point());
    triangle_points.push_back(edge.vertex2->point());
    triangle_points.push_back(edge.co_vertex->point());

    Delaunay::Point circumcenter_point = circumcenter()(triangle_points.begin(), triangle_points.end());
    double sq_radius = squared_distance()(circumcenter_point, triangle_points[0]);

    std::vector<Delaunay::Full_cell_handle> cell_neighbors;
    get_incident_cells_to_vertices(edge, dt, cell_neighbors);
    std::set<Delaunay::Vertex_handle> neighboring_vertices;

    for (auto cell : cell_neighbors) {
        for (auto v = cell->vertices_begin(); v != cell->vertices_end(); ++v) {
            neighboring_vertices.insert(*v);
        }
    }

    const double eps = 1e-8;
    for (auto v : neighboring_vertices) {
        if (v == edge.vertex1 || v == edge.vertex2 || v == edge.co_vertex) continue;

        double sq_dist = squared_distance()(v->point(), circumcenter_point);
        if (sq_dist < sq_radius - eps) {
            return false;
        }
    }

    return true;
}

void get_incident_cells_to_vertices(Edge &edge, Delaunay &dt, std::vector<Delaunay::Full_cell_handle> &cell_neighbors) {
    cell_neighbors.clear();

    std::vector<Delaunay::Full_cell_handle> v1_cells;
    dt.incident_full_cells(edge.vertex1, std::back_inserter(v1_cells));

    for (auto cell: v1_cells) {
        if (cell->has_vertex(edge.vertex2)) {
            cell_neighbors.emplace_back(cell);
        }
    }
    // std::set<Delaunay::Full_cell_handle> visited_cells;
    // std::queue<Delaunay::Full_cell_handle> cell_queue;
    //
    // Delaunay::Full_cell_handle start_cell = vertex1->full_cell();
    //
    // cell_queue.push(start_cell);
    // visited_cells.insert(start_cell);
    //
    // while (!cell_queue.empty()) {
    //     Delaunay::Full_cell_handle current_cell = cell_queue.front();
    //     cell_queue.pop();
    //
    //     if (current_cell->has_vertex(vertex2)) {
    //         cell_neighbors.push_back(current_cell);
    //     }
    //
    //     for (int i = 0; i < current_cell->maximal_dimension() + 1; ++i) {
    //         Delaunay::Full_cell_handle neighbor = current_cell->neighbor(i);
    //         if (neighbor->has_vertex(vertex1) && !visited_cells.contains(neighbor)) {
    //             visited_cells.insert(neighbor);
    //             cell_queue.push(neighbor);
    //         }
    //     }
    // }
}

void orient_voronoi_edge(std::vector<Eigen::VectorXd> shared_facet_points, Eigen::VectorXd finite_voronoi_vertex, Eigen::VectorXd &voronoi_edge_direction) {
    Eigen::VectorXd facet_center_estimate = shared_facet_points[0] + shared_facet_points[1] + shared_facet_points[2];
    facet_center_estimate /= 3;
    Eigen::VectorXd to_facet = facet_center_estimate - finite_voronoi_vertex;
    int a = voronoi_edge_direction.size();
    int b = to_facet.size();
    if (voronoi_edge_direction.dot(to_facet) < 0) {
        voronoi_edge_direction = -voronoi_edge_direction;
    }
}

void get_facet_normal(std::vector<Eigen::VectorXd> &facet_points, Eigen::VectorXd &normal) {
    Eigen::VectorXd v1 = facet_points[1] - facet_points[0];
    Eigen::VectorXd v2 = facet_points[2] - facet_points[0];

    Eigen::MatrixXd A(3, 2);
    A.col(0) = v1;
    A.col(1) = v2;

    Eigen::JacobiSVD<Eigen::MatrixXd> svd(A, Eigen::ComputeFullU);
    normal = svd.matrixU().col(2);
    normal = normal.normalized();
}

void get_shared_delaunay_facet(Delaunay::Full_cell_handle cell1, Delaunay::Full_cell_handle cell2, Delaunay &delaunay, std::array<Eigen::VectorXd, 3> &face_vertices) {
    int count = 0;
    for (int i = 0; i <= cell1->maximal_dimension(); ++i) {
        auto current_vertex = cell1->vertex(i);
        if (cell2->has_vertex(current_vertex)) {
            Eigen::VectorXd vertex_point(delaunay.maximal_dimension());
            vertex_point << current_vertex->point()[0], current_vertex->point()[1], current_vertex->point()[2];
            face_vertices[count++] = vertex_point;
        }
    }
}

void get_voronoi_facet_vertices(const Edge& edge, const Delaunay& delaunay, std::vector<Eigen::VectorXd>& facet_vertices) {
    Delaunay::Vertex_handle vertex1 = edge.vertex1;
    Delaunay::Vertex_handle vertex2 = edge.vertex2;

    std::vector<Delaunay::Full_cell_handle> v1_incident;
    delaunay.incident_full_cells(vertex1, std::back_inserter(v1_incident));

    std::set<Delaunay::Full_cell_handle> both_incident;

    for (auto incident_cell : v1_incident) {
        if (incident_cell->has_vertex(vertex2)) {
            both_incident.insert(incident_cell);
        }
    }

    for (auto cell : both_incident) {
        Eigen::VectorXd center;
        simplex_circumsphere(cell, center);
        facet_vertices.push_back(center);
    }
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

void extract_edges(Delaunay &delaunay, std::map<Delaunay::Vertex_handle, size_t> vertex_to_index, std::vector<std::array<size_t, 2>> &edges) {
    std::set<Edge> unique_edges;
    for (auto facet = delaunay.facets_begin(); facet != delaunay.facets_end(); ++facet) {
        if (!delaunay.is_infinite(*facet)) {
            std::vector<Delaunay::Vertex_handle> face_vertices;
            get_facet_vertices(delaunay, facet, face_vertices);
            if (face_vertices.size() == 2) {
                unique_edges.insert(Edge(face_vertices[0], face_vertices[1], face_vertices[2]));
            }
        }
    }
    for (const auto &edge : unique_edges) {
        edges.push_back(std::array{vertex_to_index[edge.vertex1], vertex_to_index[edge.vertex2]});
    }
}


void calculate_driver(const Eigen::VectorXd &voronoi_vertex, const Edge &delaunay_edge, Eigen::VectorXd &driver) {
    Eigen::VectorXd p1 = make_point_eigen(delaunay_edge.vertex1->point());
    Eigen::VectorXd p2 = make_point_eigen(delaunay_edge.vertex2->point());

    Eigen::VectorXd line_vec = p2 - p1;
    Eigen::VectorXd point_vec = voronoi_vertex - p1;

    double t = point_vec.dot(line_vec) / line_vec.squaredNorm();
    driver = p1 + t * line_vec;
}
