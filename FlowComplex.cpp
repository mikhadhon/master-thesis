#include <queue>

#include "FlowComplex.h"
#include "utils.h"

void flow_complex(Delaunay &delaunay, std::vector<std::array<size_t, 3>> &faces) {
    for (auto face = delaunay.facets_begin(); face != delaunay.facets_end(); ++face) {
        std::vector<Delaunay::Vertex_handle> facet_vertices;
        get_facet_vertices(delaunay, face, facet_vertices);
        if (is_index_two_critical_point(facet_vertices)) {
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

            }
        }
    }

}
