#include "data_generation.h"

Eigen::MatrixXd cliffordgen(int nsamples) {
    std::random_device rd{};
    std::mt19937 gen{ rd() };
    const double radius = 1 / sqrt(2);

    std::uniform_real_distribution<> angle_t{ 0, 2 * M_PI }, angle_p{ 0, 2 * M_PI };
    std::uniform_real_distribution<> noise{ -0.0001, 0.0001};

    Eigen::MatrixXd samples(nsamples, 4);

    for (int i = 0; i < nsamples; i++) {
        double p = angle_t(gen);
        double q = angle_p(gen);

        double a = 1;
        double b = 1 /*+ noise(gen)*/;

        auto sample = Eigen::Vector4d(a * cos(p), a * sin(p), b * cos(q), b * sin(q));
        sample *= radius;

        samples.row(i) = sample.transpose();
    }
    return samples;
}

Eigen::MatrixXd sphere2gen(int nsamples = 500) {
    //generate 3 dimensional orthonormal system
    Eigen::Vector4d w0;
    Eigen::Vector4d w1;
    Eigen::Vector4d w2;
    std::vector<Eigen::Vector4d> system(3);

    std::random_device rd{};
    std::mt19937 gen{ rd() };
    std::uniform_real_distribution<float> dist(-1, 1);

    bool dependent = true;
    while (dependent) {
        w0 = Eigen::Vector4d(dist(gen), dist(gen), dist(gen), dist(gen));
        w1 = Eigen::Vector4d(dist(gen), dist(gen), dist(gen), dist(gen));
        w2 = Eigen::Vector4d(dist(gen), dist(gen), dist(gen), dist(gen));

        dependent = false;
        for (Eigen::Vector4d i : system) {
            for (Eigen::Vector4d j : system) {
                if (i != j && i.normalized() == j.normalized()) dependent = true;
            }
        }
    }

    //Gram-Schmidt
    system[0] = w0.normalized();
    system[1] = w1 - (system[0].dot(w1) / system[0].dot(system[0])) * system[0];
    system[2] = w2 - (system[0].dot(w2) / system[0].dot(system[0])) * system[0] - (system[1].dot(w2) / system[1].dot(system[1])) * system[1];

    system[1].normalize();
    system[2].normalize();

    double radius = 1.0;

    //std::normal_distribution<> angle_t{ M_PI, M_PI_2 }, angle_p{ M_PI, M_PI_2 };
    std::uniform_real_distribution<> angle_t{ 0, 2 * M_PI }, angle_p{ 0, 2 * M_PI };
    Eigen::MatrixXd sphere_samples(nsamples, 4);

    for (int i = 0; i < nsamples; i++) {
        double p = angle_t(gen);
        double q = angle_p(gen);

        Eigen::Vector4d sample = cos(q) * (sin(p) * system[1] + cos(p) * system[0]) + sin(q) * system[2];
        sample *= radius;

        sphere_samples.row(i) = sample.transpose();
    }
    return sphere_samples;
}

Eigen::MatrixXd sphere3gen(int nsamples) {
    std::random_device rd{};
    std::mt19937 gen{ rd() };

    std::uniform_real_distribution<> angle{ 0, M_PI };
    std::uniform_real_distribution<> angle_r{ 0, 2 * M_PI };
    Eigen::MatrixXd sphere_samples(nsamples, 4);

    double radius = 1;

    for (int i = 0; i < nsamples; i++) {
        double p = angle(gen);
        double q = angle(gen);
        double r = angle_r(gen);

        double x = radius * cos(p);
        double y = radius * sin(p) * cos(q);
        double z = radius * sin(p) * sin(q) * cos(r);
        double w = radius * sin(p) * sin(q) * sin(r);

        sphere_samples.row(i) = Eigen::Vector4d(x, y, z, w).transpose();
    }
    return sphere_samples;
}
