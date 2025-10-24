#include <CGAL/Linear_algebraCd.h>

#include "Delaunay.h"

typedef CGAL::Linear_algebraCd<FT> Linear_algebra;
typedef Linear_algebra::Matrix Matrix;
typedef Linear_algebra::Vector Vector;
typedef Linear_algebra::RT RT;

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

bool is_index_two_critical_point(const std::vector<Delaunay::Vertex_handle> &facet_vertices) {

    Eigen::Vector3d i(facet_vertices[0]->point()[0], facet_vertices[0]->point()[1], facet_vertices[0]->point()[2]);
    Eigen::Vector3d j(facet_vertices[1]->point()[0], facet_vertices[1]->point()[1], facet_vertices[1]->point()[2]);
    Eigen::Vector3d l(facet_vertices[2]->point()[0], facet_vertices[2]->point()[1], facet_vertices[2]->point()[2]);

    return (j-i).dot(l-i) > 0 && (i-j).dot(l-j) > 0 && (i-l).dot(j-l) > 0;
}

bool is_gabriel(Edge &edge) {
    std::vector<Delaunay::Full_cell_handle> cell_neighbors;
    get_incident_cells_to_vertices(edge, cell_neighbors);
    std::set<Delaunay::Vertex_handle> neighboring_vertices;

    for (auto cell : cell_neighbors) {
        for (auto v = cell->vertices_begin(); v != cell->vertices_end(); ++v) {
            neighboring_vertices.insert(*v);
        }
    }

    Eigen::Vector3d i(edge.vertex1->point()[0], edge.vertex1->point()[1], edge.vertex1->point()[2]);
    Eigen::Vector3d j(edge.vertex2->point()[0], edge.vertex2->point()[1], edge.vertex2->point()[2]);

    Eigen::Vector3d center = (i + j) / 2;
    double radius = (i - j).norm() / 2;

    for (auto vertex : neighboring_vertices) {
        Eigen::Vector3d v(vertex->point()[0], vertex->point()[1], vertex->point()[2]);
        if ((v - center).squaredNorm() < radius * radius) {
            return false;
        }
    }

    return true;
}

void get_incident_cells_to_vertices(Edge &edge, std::vector<Delaunay::Full_cell_handle> &cell_neighbors) {
    Delaunay::Vertex_handle vertex1 = edge.vertex1;
    Delaunay::Vertex_handle vertex2 = edge.vertex2;

    std::set<Delaunay::Full_cell_handle> visited_cells;
    std::queue<std::pair<Delaunay::Vertex_handle, Delaunay::Full_cell_handle>> cell_queue;

    Delaunay::Full_cell_handle v1_full_cell = vertex1->full_cell();
    Delaunay::Full_cell_handle v2_full_cell = vertex2->full_cell();

    cell_queue.push(std::make_pair(vertex1, v1_full_cell));
    cell_queue.push(std::make_pair(vertex2, v2_full_cell));

    visited_cells.insert(v1_full_cell);
    visited_cells.insert(v2_full_cell);

    while (!cell_queue.empty()) {
        Delaunay::Vertex_handle current_vertex = cell_queue.front().first;
        Delaunay::Full_cell_handle current_cell = cell_queue.front().second;
        cell_queue.pop();

        for (int i = 0; i < current_cell->maximal_dimension() + 1; ++i) {
            Delaunay::Full_cell_handle neighbor = current_cell->neighbor(i);
            if (neighbor->has_vertex(current_vertex) && !visited_cells.contains(neighbor)) {
                visited_cells.insert(neighbor);
                cell_neighbors.push_back(neighbor);
                cell_queue.push(std::make_pair(current_vertex, neighbor));
            }
        }
    }
}

void voronoi_facet_from_edge(Edge &edge, std::vector<std::pair<Delaunay::Full_cell_handle, Delaunay::Full_cell_handle>> &facet_edges, Delaunay &delaunay) {
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

    int d = delaunay.maximal_dimension();
    for (auto cell : both_incident) {
        for (int i = 0; i <= d; ++i) {
            auto neighbor = cell->neighbor(i);
            if (both_incident.contains(neighbor) && cell < neighbor) {
                facet_edges.push_back(std::make_pair(neighbor, cell));
            }
        }
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
                unique_edges.insert(Edge(face_vertices[0], face_vertices[1]));
            }
        }
    }
    for (const auto &edge : unique_edges) {
        edges.push_back(std::array{vertex_to_index[edge.vertex1], vertex_to_index[edge.vertex2]});
    }
}
