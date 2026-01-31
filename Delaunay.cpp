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

void get_shared_delaunay_facet(Delaunay::Full_cell_handle cell1, Delaunay::Full_cell_handle cell2, Delaunay &delaunay, std::vector<Point> &face_vertices) {
    int count = 0;
    for (int i = 0; i <= cell1->maximal_dimension(); ++i) {
        auto current_vertex = cell1->vertex(i);
        if (cell2->has_vertex(current_vertex)) {
            face_vertices.push_back(current_vertex->point());
        }
    }
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

    for (size_t i = 0; i < ordered_incident_cells.size(); ++i) {
        Voronoi_vertex v_vertex;

        if (!is_infinite_tet[i]) {
            Eigen::VectorXd circumcenter = simplex_circumsphere(ordered_incident_cells[i]);
            v_vertex = Voronoi_vertex{false, circumcenter};
        } else {
            Delaunay::Full_cell_handle finite_neighbor = is_infinite_tet[(i + 1) % is_infinite_tet.size()] ? ordered_incident_cells[(i - 1) % is_infinite_tet.size()] : ordered_incident_cells[(i + 1) % is_infinite_tet.size()];
            Eigen::VectorXd finite_circumcenter = simplex_circumsphere(finite_neighbor);

            std::vector<Point> shared_facet_points;
            get_shared_delaunay_facet(ordered_incident_cells[i], finite_neighbor, dt, shared_facet_points);

            std::vector<Eigen::VectorXd> facet_points_eigen;
            for (const auto& p : shared_facet_points) {
                facet_points_eigen.push_back(make_point_eigen(p));
            }

            Eigen::VectorXd normal;
            get_facet_normal(facet_points_eigen, normal);

            Eigen::VectorXd opposite_vertex;
            for (auto v = finite_neighbor->vertices_begin(); v != finite_neighbor->vertices_end(); ++v) {
                if (!dt.is_infinite(*v) && !ordered_incident_cells[i]->has_vertex(*v)) {
                    opposite_vertex = make_point_eigen((*v)->point());
                    break;
                }
            }

            Eigen::VectorXd to_opposite = opposite_vertex - facet_points_eigen[0];
            if (normal.dot(to_opposite) > 0) {
                normal = -normal;
            }

            Eigen::VectorXd infinite_point = finite_circumcenter + 2 * normal;
            v_vertex = Voronoi_vertex{true, infinite_point};
            infinite_vertices.push_back(infinite_point);
        }

        voronoi_vertices.push_back(v_vertex);
        ps_vertices.push_back(v_vertex.point);
    }

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
    polyscope::registerPointCloud("infinite voronoi vertices", infinite_vertices);
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
        voronoi_vertex0 = simplex_circumsphere(incident_cell0);
        voronoi_vertex1 = simplex_circumsphere(incident_cell1);
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
        Eigen::VectorXd finite_circumcenter = simplex_circumsphere(finite_cell);

        Eigen::VectorXd opposite_vertex;
        for (auto v = finite_cell->vertices_begin(); v != finite_cell->vertices_end(); ++v) {
            if (!dt.is_infinite(*v) && !infinite_cell->has_vertex(*v)) {
                opposite_vertex = make_point_eigen((*v)->point());
                break;
            }
        }

        Eigen::VectorXd to_opposite = opposite_vertex - final_facet_vertices[0];
        if (normal.dot(to_opposite) > 0) {
            normal = -normal;
        }

        Voronoi_vertex infinite_voronoi_vertex(true, finite_circumcenter + 1000 * normal);
        Voronoi_vertex finite_voronoi_vertex(false, finite_circumcenter);

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

Eigen::VectorXd simplex_circumsphere(Delaunay::Full_cell_handle simplex) {
    std::vector<Point> simplex_points;
    for (auto point = simplex->vertices_begin(); point != simplex->vertices_end(); ++point) {
        simplex_points.push_back((*point)->point());
    }
    Point cgal_center = circumcenter()(simplex_points.begin(), simplex_points.end());
    return make_point_eigen(cgal_center);
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
    Eigen::VectorXd p1 = make_point_eigen(edge.vertex1->point());
    Eigen::VectorXd p2 = make_point_eigen(edge.vertex2->point());

    Eigen::VectorXd edge_midpoint = (p1 + p2) / 2.0;
    double sq_radius = (p1 - edge_midpoint).squaredNorm();

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
        if (v == edge.vertex1 || v == edge.vertex2) continue;
        if (dt.is_infinite(v)) continue;

        double sq_dist = (make_point_eigen(v->point()) - edge_midpoint).squaredNorm();
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
        Eigen::VectorXd center = simplex_circumsphere(cell);
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

void find_gabriel_edges(Delaunay &dt, std::vector<Edge> &gabriel_edges) {
     std::unordered_set<Edge, EdgeHash> visited_edges;
    int dim = dt.maximal_dimension();

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
