#ifndef FLOWCOMPLEX_DATA_GEN_H
#define FLOWCOMPLEX_DATA_GEN_H

#include <Eigen/Core>

void generate_trefoil(int nsamples, int normal_windings, Eigen::Matrix<double, Eigen::Dynamic, 3> &V, Eigen::Matrix<double, Eigen::Dynamic, 3>  &N1, Eigen::Matrix<double, Eigen::Dynamic, 3> &N2, Eigen::Matrix<double, Eigen::Dynamic, 3> &T);

void generate_spiral(int nsamples, int normal_windings, Eigen::Matrix<double, Eigen::Dynamic, 3> &V, Eigen::Matrix<double, Eigen::Dynamic, 3>  &N1, Eigen::Matrix<double, Eigen::Dynamic, 3> &N2,  Eigen::Matrix<double, Eigen::Dynamic, 3> &T);

void generate_circle(int nsamples, int normal_windings, Eigen::Matrix<double, Eigen::Dynamic, 3> &V, Eigen::Matrix<double, Eigen::Dynamic, 3> &N1, Eigen::Matrix<double, Eigen::Dynamic, 3>  &N2, Eigen::Matrix<double, Eigen::Dynamic, 3> &T);

void read_point_data(const std::string filename, Eigen::MatrixXd &V, Eigen::MatrixXd &N1, Eigen::MatrixXd &N2 );

void write_point_data(const std::string filename, const Eigen::MatrixXd &V, const Eigen::MatrixXd &N1, const Eigen::MatrixXd &N2 );

#endif //FLOWCOMPLEX_DATA_GEN_H