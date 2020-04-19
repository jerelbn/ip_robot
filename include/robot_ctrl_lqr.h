#pragma once

#include <fstream>

#include "common_cpp/common.h"
#include "common_cpp/logger.h"
#include "robot_common.h"
#include "lin_alg_tools/care.h"

namespace robot
{

class LQRController
{

public:
    LQRController();
    LQRController(const std::string &filename);
    ~LQRController();

    void load(const std::string &filename);
    void computeControl(const double &t, const xVector &x, const double &dx_d, const double &dpsi_d);

    const Eigen::Vector2d &u() const { return u_; }

private:
    CareSolver<6, NUM_INPUTS> solver_;

    bool initialized_;
    double t_prev_;

    uVector u_;
    Eigen::Matrix<double,6,6> A_, P_, Q_;
    Eigen::Matrix<double,6,2> B_;
    Eigen::Matrix2d R_, R_inv_;
    Eigen::Matrix<double, NUM_INPUTS, 6> K_;

    int update_rate_;
    double mc_;
    double mp_;
    double Jy_;
    double Jz_;
    double Jp_;
    double Jm_;
    double Km_;
    double Lm_;
    double Rm_;
    double L_;
    double l_;
    double r_;
    double bm_;
    double bp_;

    common::Logger command_log_;

    void f(const Eigen::Matrix<double,6,1> &x, const uVector &u, Eigen::Matrix<double,6,1> &dx);
    void numericalAB(const Eigen::Matrix<double,6,1> &x, const uVector &u, Eigen::Matrix<double,6,6> &A, Eigen::Matrix<double,6,2> &B);
    void log(const double &t);
};

} // namespace robot