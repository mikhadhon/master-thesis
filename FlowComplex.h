#ifndef FLOWCOMPLEX_H
#define FLOWCOMPLEX_H
#include "Delaunay.h"

void flow_complex(Delaunay &delaunay, std::vector<Eigen::VectorXd> &vertices, std::vector<std::array<size_t, 3>> &faces, std::map<Delaunay::Vertex_handle, size_t> &vertex_to_index, std::vector<std::vector<std::array<size_t, 3>>> &index_2_stable_manifolds);
#endif
