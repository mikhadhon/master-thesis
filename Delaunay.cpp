#include <CGAL/Linear_algebraCd.h>

#include "Delaunay.h"

#include "utils.h"

typedef CGAL::Linear_algebraCd<FT> Linear_algebra;
typedef Linear_algebra::Matrix Matrix;
typedef Linear_algebra::Vector Vector;
typedef Linear_algebra::RT RT;

Voronoi_face delaunay_edge_dual(Edge &edge, Delaunay::Facet_iterator &df, Delaunay &dt) {
    std::vector<Delaunay::Vertex_handle> face_vertices;
    get_facet_vertices(dt, df, face_vertices);

    std::vector<Delaunay::Full_cell_handle> incident_cells;
    incident_cells.push_back(df->full_cell());
    incident_cells.push_back(incident_cells[0]->neighbor(dt.index_of_covertex(*df)));

    Delaunay::Full_cell_handle start_cell;
    Delaunay::Full_cell_handle current_cell;
    Delaunay::Full_cell_handle next_cell;

    if (dt.is_infinite(incident_cells[0])) {
        start_cell = incident_cells[0];
        current_cell = incident_cells[0];
        next_cell = incident_cells[1];
    }
    else {
        start_cell = incident_cells[1];
        current_cell = incident_cells[1];
        next_cell = incident_cells[0];
    }

    Delaunay::Vertex_handle co_vertex_current_next;
    for (auto vertex = start_cell->vertices_begin(); vertex != start_cell->vertices_end(); ++vertex) {
        if (*vertex != face_vertices[0] && *vertex != face_vertices[1] && *vertex != face_vertices[2]) {
            co_vertex_current_next = *vertex;
        }
    }

    bool test = dt.is_infinite(co_vertex_current_next);

    std::vector<Voronoi_vertex> voronoi_vertices;
    std::vector<Voronoi_edge> voronoi_edges;
    std::vector<Eigen::VectorXd> ps_vertices;
    std::vector<std::array<size_t, 2>> ps_edges;
    Eigen::VectorXd voronoi_vertex_start(dt.maximal_dimension());

    if (!dt.is_infinite(co_vertex_current_next)) {
        simplex_circumsphere(start_cell, voronoi_vertex_start);
        Voronoi_vertex voronoi_vertex(false, voronoi_vertex_start);
        voronoi_vertices.push_back(voronoi_vertex);
        ps_vertices.push_back(voronoi_vertex.point);
    }
    else {
        std::vector<Eigen::VectorXd> facet_points = {make_point_eigen(face_vertices[0]->point()), make_point_eigen(face_vertices[1]->point()), make_point_eigen(face_vertices[2]->point())};
        Eigen::VectorXd normal;
        get_facet_normal(facet_points, normal);
        simplex_circumsphere(next_cell, voronoi_vertex_start);
        orient_voronoi_edge(facet_points, voronoi_vertex_start, normal);
        voronoi_vertex_start = voronoi_vertex_start + 2*normal;
        Voronoi_vertex voronoi_vertex(true, voronoi_vertex_start);
        voronoi_vertices.push_back(voronoi_vertex);
        ps_vertices.push_back(voronoi_vertex.point);
    }

    while (next_cell != start_cell && !dt.is_infinite(next_cell)) {
        Eigen::VectorXd voronoi_vertex_next(3);

        simplex_circumsphere(next_cell, voronoi_vertex_next);
        Voronoi_vertex voronoi_vertex(false, voronoi_vertex_next);
        voronoi_vertices.push_back(voronoi_vertex);
        ps_vertices.push_back(voronoi_vertex.point);
        voronoi_edges.push_back(
            Voronoi_edge(
                voronoi_vertices[voronoi_vertices.size() - 2],
                voronoi_vertices[voronoi_vertices.size() - 1],
                current_cell,
                next_cell
            )
        );
        ps_edges.push_back({voronoi_vertices.size() - 2, voronoi_vertices.size() - 1});

        for (int i = 0; i < 3; ++i) {
            if (face_vertices[i] != edge.vertex1 && face_vertices[i] != edge.vertex2) {
                Delaunay::Vertex_handle temp = current_cell->mirror_vertex(current_cell->index(co_vertex_current_next), dt.maximal_dimension());
                co_vertex_current_next = face_vertices[i];
                face_vertices[i] = temp;
                current_cell = next_cell;
                next_cell = next_cell->neighbor(next_cell->index(co_vertex_current_next));
                break;
            }
        }
    }

    if (!dt.is_infinite(next_cell)) {
        voronoi_edges.push_back(
            Voronoi_edge(
                voronoi_vertices[voronoi_vertices.size() - 1],
                voronoi_vertices[0],
                next_cell,
                start_cell
            )
        );
        ps_edges.push_back({voronoi_vertices.size() - 1, 0});
    }
    else {
        Eigen::VectorXd voronoi_vertex_next(3);

        simplex_circumsphere(current_cell, voronoi_vertex_next);
        std::vector<Eigen::VectorXd> final_facet_vertices;
        for (auto vertex = current_cell->vertices_begin(); vertex != current_cell->vertices_end(); ++vertex) {
            if (*vertex != co_vertex_current_next) {
                final_facet_vertices.push_back(make_point_eigen((*vertex)->point()));
            }
        }
        Eigen::VectorXd normal;
        get_facet_normal(final_facet_vertices, normal);
        orient_voronoi_edge(final_facet_vertices, voronoi_vertex_next, normal);
        Voronoi_vertex voronoi_vertex(true, voronoi_vertex_next + 2*normal);
        voronoi_vertices.push_back(voronoi_vertex);
        ps_vertices.push_back(voronoi_vertex.point);
        voronoi_edges.push_back(
            Voronoi_edge(
                voronoi_vertices[voronoi_vertices.size() - 2],
                voronoi_vertices[voronoi_vertices.size() - 1],
                current_cell,
                current_cell->neighbor(current_cell->index(co_vertex_current_next))
            )
        );
        ps_edges.push_back({voronoi_vertices.size() - 2, voronoi_vertices.size() - 1});
    }
    return Voronoi_face(voronoi_vertices, voronoi_edges, ps_vertices, ps_edges);
}

Voronoi_edge delaunay_face_dual(Delaunay::Facet_iterator &face, Delaunay &dt) {
    Delaunay::Full_cell_handle incident_cell0, incident_cell1;
    incident_cell0 = face->full_cell();
    incident_cell1 = incident_cell0->neighbor(face->index_of_covertex());

    bool infinite0 = dt.is_infinite(incident_cell0);
    bool infinite1 = dt.is_infinite(incident_cell1);

    Eigen::VectorXd voronoi_vertex0, voronoi_vertex1;
    if (!infinite0 && !infinite1) {
        simplex_circumsphere(incident_cell0, voronoi_vertex0);
        simplex_circumsphere(incident_cell1, voronoi_vertex1);
        return Voronoi_edge(Voronoi_vertex(infinite0, voronoi_vertex0), Voronoi_vertex(infinite1, voronoi_vertex1), incident_cell0, incident_cell1);
    }
    else {
        std::vector<Delaunay::Vertex_handle> face_vertex_handles;
        get_facet_vertices(dt, face, face_vertex_handles);
        std::vector<Eigen::VectorXd> final_facet_vertices;

        final_facet_vertices.push_back(make_point_eigen(face_vertex_handles[0]->point()));
        final_facet_vertices.push_back(make_point_eigen(face_vertex_handles[1]->point()));
        final_facet_vertices.push_back(make_point_eigen(face_vertex_handles[2]->point()));

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

Delaunay::Face voronoi_edge_dual(Voronoi_edge &voronoi_edge) {
    Delaunay::Face face(voronoi_edge.cell1);
    int vertex_count = 0;

    for (auto vertex = voronoi_edge.cell1->vertices_begin(); vertex != voronoi_edge.cell1->vertices_end(); ++vertex) {
        int vertex_index_in_cell1 = voronoi_edge.cell1->index(*vertex);
        if (voronoi_edge.cell1->neighbor(vertex_index_in_cell1) == voronoi_edge.cell2) continue;
        face.set_index(vertex_count, vertex_index_in_cell1);
    }
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

bool is_index_two_critical_point(Delaunay::Facet_iterator &face, Delaunay &dt) {
    Voronoi_edge face_dual = delaunay_face_dual(face, dt);

    std::vector<Delaunay::Vertex_handle> face_vertex_handles;
    Eigen::VectorXd face_normal;
    get_facet_vertices(dt, face, face_vertex_handles);

    std::vector<Eigen::VectorXd> face_vertices = {make_point_eigen(face_vertex_handles[0]->point()), make_point_eigen(face_vertex_handles[1]->point()), make_point_eigen(face_vertex_handles[2]->point())};
    get_facet_normal(face_vertices, face_normal);

    double vertex0_direction = face_normal.dot(face_dual.vertex1.point - face_vertices[0]);
    double vertex1_direction = face_normal.dot(face_dual.vertex2.point - face_vertices[1]);

    Eigen::VectorXd i = make_point_eigen(face_vertex_handles[0]->point());
    Eigen::VectorXd j = make_point_eigen(face_vertex_handles[1]->point());
    Eigen::VectorXd l = make_point_eigen(face_vertex_handles[2]->point());

    return (vertex0_direction * vertex1_direction < 0) && (j-i).dot(l-i) > 0 && (i-j).dot(l-j) > 0 && (i-l).dot(j-l) > 0;
}

bool is_gabriel(Edge &edge) {
    std::vector<Delaunay::Point> triangle_points;
    triangle_points.push_back(edge.vertex1->point());
    triangle_points.push_back(edge.vertex2->point());
    triangle_points.push_back(edge.co_vertex->point());

    Delaunay::Point circumcenter_point = circumcenter()(triangle_points.begin(), triangle_points.end());
    double sq_radius = squared_distance()(circumcenter_point, triangle_points[0]);

    std::vector<Delaunay::Full_cell_handle> cell_neighbors;
    get_incident_cells_to_vertices(edge, cell_neighbors);
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

void get_incident_cells_to_vertices(Edge &edge, std::vector<Delaunay::Full_cell_handle> &cell_neighbors) {
    Delaunay::Vertex_handle vertex1 = edge.vertex1;
    Delaunay::Vertex_handle vertex2 = edge.vertex2;

    std::set<Delaunay::Full_cell_handle> visited_cells;
    std::queue<Delaunay::Full_cell_handle> cell_queue;

    Delaunay::Full_cell_handle start_cell = vertex1->full_cell();

    cell_queue.push(start_cell);
    visited_cells.insert(start_cell);

    while (!cell_queue.empty()) {
        Delaunay::Full_cell_handle current_cell = cell_queue.front();
        cell_queue.pop();

        if (current_cell->has_vertex(vertex2)) {
            cell_neighbors.push_back(current_cell);
        }

        for (int i = 0; i < current_cell->maximal_dimension() + 1; ++i) {
            Delaunay::Full_cell_handle neighbor = current_cell->neighbor(i);
            if (neighbor->has_vertex(vertex1) && !visited_cells.contains(neighbor)) {
                visited_cells.insert(neighbor);
                cell_queue.push(neighbor);
            }
        }
    }
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
