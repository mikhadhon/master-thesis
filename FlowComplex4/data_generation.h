#ifndef FLOWCOMPLEX4_DATA_GENERATION_H
#define FLOWCOMPLEX4_DATA_GENERATION_H

#define _USE_MATH_DEFINES

#include <random>

#include "Delaunay.h"

Eigen::MatrixXd sphere3gen(int nsamples);

Eigen::MatrixXd sphere2gen(int nsamples);

Eigen::MatrixXd cliffordgen(int nsamples);

#endif //FLOWCOMPLEX4_DATA_GENERATION_H