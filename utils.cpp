#include "utils.h"

void map_vertices_to_vector(
    Delaunay &delaunay,
    std::vector<std::array<double, 3> > &vertices,
    std::map<Delaunay::Vertex_handle, size_t> &vertex_to_index
) {
    size_t vertex_index = 0;
    for (auto vertex_it = delaunay.vertices_begin();
        vertex_it != delaunay.vertices_end(); ++vertex_it) {
        if (!delaunay.is_infinite(vertex_it)) {
            Delaunay::Point p = vertex_it->point();
            vertices.push_back({p[0], p[1], p[2]});
            vertex_to_index[vertex_it] = vertex_index++;
        }
    }
}

void write_faces_to_vector(
    Delaunay &delaunay,
    std::vector<std::array<size_t, 3> > &faces,
    std::map<Delaunay::Vertex_handle, size_t> &vertex_to_index
) {
    for (auto facet = delaunay.facets_begin();
        facet != delaunay.facets_end(); ++facet) {
        if (!delaunay.is_infinite(*facet)){
            const Delaunay::Full_cell_handle cell = facet->full_cell();
            const int opposite_vertex = facet->index_of_covertex();

            std::array<size_t, 3> face;
            size_t face_vertex_idx = 0;

            for (int i = 0; i < 4; ++i) {
                if (i != opposite_vertex) {
                    Delaunay::Vertex_handle vh = cell->vertex(i);
                    face[face_vertex_idx++] = vertex_to_index[vh];
                    }
                }

                faces.push_back(face);
            }
    }
}
