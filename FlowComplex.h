#ifndef FLOWCOMPLEX_H
#define FLOWCOMPLEX_H
#include "Delaunay.h"

void flow_complex(Delaunay &delaunay, std::vector<std::array<double, 3> > &vertices, std::vector<std::array<size_t, 3>> &faces, std::map<Delaunay::Vertex_handle, size_t> &vertex_to_index, std::vector<Eigen::Vector3d> &centers, std::vector<Eigen::Vector3d> &centers2);
#endif
