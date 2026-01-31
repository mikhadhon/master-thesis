#include "data_gen.h"
#include <fstream>
#include <iostream>


void generate_trefoil(int nsamples, int normal_windings, Eigen::Matrix<double, Eigen::Dynamic, 3> &V, Eigen::Matrix<double, Eigen::Dynamic, 3>  &N1, Eigen::Matrix<double, Eigen::Dynamic, 3> &N2, Eigen::Matrix<double, Eigen::Dynamic, 3> &T){
    double r = 1.;

    V.resize(nsamples, 3);
    N1.resize(nsamples, 3);
    N2.resize(nsamples, 3);
    T.resize(nsamples, 3);

    // double stretchz = 2.7;
    double stretchz = 1.;

    for (int i=0; i<nsamples; i++) {
        double t = (1.*i) / nsamples;
        double dtheta = (2.*M_PI);
        double theta  = dtheta*t;

        V.row(i)  = Eigen::Vector3d(    cos(theta) + 2*cos(2*theta),
                                        sin(theta) - 2*sin(2*theta),
                                       -sin(3*theta) * stretchz);
        Eigen::Vector3d Tangent = Eigen::Vector3d( -sin(theta) - 4*sin(2*theta),
                                                    cos(theta) - 4*cos(2*theta),
                                                 -3*cos(3*theta) * stretchz ).normalized(); 
        T.row(i)  = Tangent;
        Eigen::Vector3d dTudt = Eigen::Vector3d(  -cos(theta) - 8*cos(2*theta),
                                                  -sin(theta) + 8*sin(2*theta),
                                                 9*sin(3*theta) * stretchz);
        N1.row(i) = (dTudt - dTudt.dot(Tangent) * Tangent).normalized();

        Eigen::Vector3d v1 = T.row(i);
        Eigen::Vector3d v2 = N1.row(i);
        N2.row(i) = v1.cross(v2);
        N2.row(i).normalize();

        if (normal_windings > 0) {
            theta = (2.*M_PI/nsamples)*i * normal_windings;
            Eigen::Matrix3d frame;
            frame.row(0) = N1.row(i);
            frame.row(1) = N2.row(i);
            frame.row(2) =  T.row(i);
            Eigen::Matrix3d Rz; Rz << cos(theta), -sin(theta), 0, sin(theta), cos(theta), 0, 0, 0, 1;
            Eigen::Matrix3d R =  frame.transpose() * Rz * frame;

            N1.row(i) = N1.row(i) * R.transpose();
            N2.row(i) = N2.row(i) * R.transpose();
        }
    }

    std::cout << "trefoil, min/max per dim: " << std::endl;;
    std::cout << V.col(0).minCoeff() << "-" << V.col(0).maxCoeff() << std::endl;
    std::cout << V.col(1).minCoeff() << "-" << V.col(1).maxCoeff() << std::endl;
    std::cout << V.col(2).minCoeff() << "-" << V.col(2).maxCoeff() << std::endl;

}
