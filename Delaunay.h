#ifndef DELAUNAY_H
#define DELAUNAY_H

#include <CGAL/Epick_d.h>
#include <CGAL/Delaunay_triangulation.h>

typedef CGAL::Dimension_tag<3> Dimension;
typedef CGAL::Epick_d<Dimension> K;
typedef K::Point_d Point;
typedef K::FT FT;

typedef CGAL::Delaunay_triangulation<CGAL::Epick_d<Dimension>> Delaunay;

void insertPoints(std::vector<Point> &points, Delaunay &delaunay);

void indexTwoCriticalPoint(
    const Delaunay &delaunay,
    const Delaunay::Facet_iterator &facet,
    Point &criticalPoint
);

void getFacetVertices(
    const Delaunay &delaunay,
    const Delaunay::Facet_iterator &facet,
    std::vector<Delaunay::Vertex_handle> &facet_vertices
);

#endif
