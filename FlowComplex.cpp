#include <queue>
#include <Eigen/StdVector>

#include "FlowComplex.h"
#include "utils.h"

void flow_complex(Delaunay &delaunay, std::vector<std::array<double, 3> > &vertices, std::vector<std::array<size_t, 3>> &faces, std::map<Delaunay::Vertex_handle, size_t> &vertex_to_index) {
    std::map<std::array<double, 3>, size_t> fc_vertex_to_index;
    int vertex_count = vertices.size();

    for (auto face = delaunay.facets_begin(); face != delaunay.facets_end(); ++face) {
        std::vector<Delaunay::Vertex_handle> facet_vertices;
        get_facet_vertices(delaunay, face, facet_vertices);
        if (facet_vertices.size() == delaunay.maximal_dimension() && is_index_two_critical_point(facet_vertices)) {
            std::queue<std::pair<Eigen::Vector3d, Edge>> edge_queue;

            Eigen::Vector3d center;
            double radius;
            Eigen::Vector3d i(facet_vertices[0]->point()[0], facet_vertices[0]->point()[1], facet_vertices[0]->point()[2]);
            Eigen::Vector3d j(facet_vertices[1]->point()[0], facet_vertices[1]->point()[1], facet_vertices[1]->point()[2]);
            Eigen::Vector3d l(facet_vertices[2]->point()[0], facet_vertices[2]->point()[1], facet_vertices[2]->point()[2]);

            circumcircle(i, j, l, center, radius);

            edge_queue.push(std::make_pair(center, Edge(facet_vertices[0], facet_vertices[1])));
            edge_queue.push(std::make_pair(center, Edge(facet_vertices[0], facet_vertices[2])));
            edge_queue.push(std::make_pair(center, Edge(facet_vertices[1], facet_vertices[2])));

            while (!edge_queue.empty()) {
                auto current = edge_queue.front();
                edge_queue.pop();
                Edge current_edge = current.second;
                if (is_gabriel(current_edge)) {
                    std::array fc_center = {center[0], center[1], center[2]};
                    auto fc_index = fc_vertex_to_index.find(fc_center);
                    std::array<size_t, 3> fc_face{};
                    if (fc_index == fc_vertex_to_index.end()) {
                        fc_vertex_to_index[fc_center] = vertex_count++;
                    }
                    fc_face[0] = vertex_to_index[current_edge.vertex1];
                    fc_face[1] = vertex_to_index[current_edge.vertex2];
                    fc_face[2] = fc_vertex_to_index[fc_center];
                    faces.push_back(fc_face);
                }
            }
        }
    }

    vertices.resize(vertex_count);
    for (auto fc_vertex : fc_vertex_to_index) {
        vertices[fc_vertex.second] = fc_vertex.first;
    }
}
