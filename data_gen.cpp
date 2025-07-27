#include "data_gen.hpp"
#include <fstream>
#include <iostream>

void write_point_data(const std::string filename, const Eigen::MatrixXd &V, const Eigen::MatrixXd &N1, const Eigen::MatrixXd &N2 ){

    int n = V.rows();

    std::ofstream ofile(filename);

    ofile << n << std::endl;

    for (int i=0; i<n; i++) {
        ofile << V(i,0) << std::endl; 
        ofile << V(i,1) << std::endl; 
        ofile << V(i,2) << std::endl; 

        ofile << N1(i,0) << std::endl; 
        ofile << N1(i,1) << std::endl; 
        ofile << N1(i,2) << std::endl; 

        ofile << N2(i,0) << std::endl; 
        ofile << N2(i,1) << std::endl; 
        ofile << N2(i,2) << std::endl; 
    }
    ofile.close();
}

void read_point_data(const std::string filename, Eigen::MatrixXd &V, Eigen::MatrixXd &N1, Eigen::MatrixXd &N2 ){

    std::ifstream ifile(filename);
    int n;
    ifile >> n;

    V.resize(n,3);
    N1.resize(n,3);
    N2.resize(n,3);

    for (int i=0; i<n; i++) {
        ifile >> V(i,0); 
        ifile >> V(i,1); 
        ifile >> V(i,2); 

        ifile >> N1(i,0); 
        ifile >> N1(i,1); 
        ifile >> N1(i,2); 

        ifile >> N2(i,0); 
        ifile >> N2(i,1); 
        ifile >> N2(i,2); 
    }
    ifile.close();
}

void generate_circle(int nsamples, int normal_windings, Eigen::Matrix<double, Eigen::Dynamic, 3> &V, Eigen::Matrix<double, Eigen::Dynamic, 3> &N1, Eigen::Matrix<double, Eigen::Dynamic, 3>  &N2, Eigen::Matrix<double, Eigen::Dynamic, 3> &T){
    double r = 1.;
    V.resize(nsamples, 3);
    N1.resize(nsamples, 3);
    N2.resize(nsamples, 3);
    T.resize(nsamples, 3);
    for (int i=0; i<nsamples; i++) {
        double theta = (2.*M_PI/nsamples)*i;
        V.row(i)  = Eigen::Vector3d(cos(theta), sin(theta), 1.);
        N1.row(i) = Eigen::Vector3d(cos(theta), sin(theta), 0.).normalized();
        N2.row(i) = Eigen::Vector3d(        0.,         0., 1.).normalized();
        T.row(i)  = N1.row(i).cross(N2.row(i));
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
}

void generate_spiral(int nsamples, int normal_windings, Eigen::Matrix<double, Eigen::Dynamic, 3> &V, Eigen::Matrix<double, Eigen::Dynamic, 3>  &N1, Eigen::Matrix<double, Eigen::Dynamic, 3> &N2,  Eigen::Matrix<double, Eigen::Dynamic, 3> &T){
    double r = 1.;
    int turns = 2;
    double height = 2.;

    V.resize(nsamples, 3);
    N1.resize(nsamples, 3);
    N2.resize(nsamples, 3);
    T.resize(nsamples, 3);
    for (int i=0; i<nsamples; i++) {
        double dtheta = (turns*2.*M_PI/nsamples);
        double theta  = dtheta*i;
        double dz = (height / nsamples);
        double z  = dz*i;
        V.row(i)  = Eigen::Vector3d(cos(theta),                sin(theta), z);
        T.row(i)  = Eigen::Vector3d(-sin(theta)*dtheta, cos(theta)*dtheta,dz).normalized();
        N1.row(i) = Eigen::Vector3d(-cos(theta)*dtheta, -sin(theta)*dtheta, 0.);
        // Eigen::Vector3d N = (N1.row(i) - N1.row(i).dot(T.row(i))*T.row(i)).normalized();
        // N1.row(i) = N;
        Eigen::Vector3d N = (N1.row(i) - N1.row(i).dot(T.row(i))*T.row(i)).normalized();
        N1.row(i) = N;
        N2.row(i) = T.row(i).cross(N1.row(i));

        if (normal_windings > 0) {
            theta = (2.*M_PI/nsamples)*i * normal_windings;
            Eigen::Matrix3d frame;
            frame.row(0) = N1.row(i);
            frame.row(1) = N2.row(i);
            frame.row(2) =  T.row(i); // tangent vector
            Eigen::Matrix3d Rz; Rz << cos(theta), -sin(theta), 0, sin(theta), cos(theta), 0, 0, 0, 1;
            Eigen::Matrix3d R =  frame.transpose() * Rz * frame;

            N1.row(i) = N1.row(i) * R.transpose();
            N2.row(i) = N2.row(i) * R.transpose();
        }
    }
}

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

        N2.row(i) = T.row(i).cross(N1.row(i));
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
