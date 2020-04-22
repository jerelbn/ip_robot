#pragma once

#include "common_cpp/common.h"
#include "common_cpp/logger.h"
#include "robot_common.h"
#include "lin_alg_tools/care.h"

#define NS 8 // number of states

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
    CareSolver<NS, NUM_INPUTS> solver_;

    bool initialized_;
    double t_prev_;

    uVector u_;
    Eigen::Matrix<double,NS,NS> A_, P_, Q_;
    Eigen::Matrix<double,NS,2> B_;
    Eigen::Matrix2d R_, R_inv_;
    Eigen::Matrix<double, NUM_INPUTS, NS> K_;
    double dx_d_, dx_int_;
    double dpsi_d_, dpsi_int_;

    int update_rate_;
    double max_voltage_;
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

    void f(const Eigen::Matrix<double,NS,1> &x, const uVector &u, Eigen::Matrix<double,NS,1> &dx) const;
    void ftilde(const double &dt, const Eigen::Matrix<double,NS,1> &xref, const Eigen::Matrix<double,NS,1> &xtilde,
                const uVector &u, Eigen::Matrix<double,NS,1> &dxtilde) const;
    void numericalAB(const Eigen::Matrix<double, NS, 1> &x, const Eigen::Matrix<double, NS, 1> &xref, const uVector &u);
    void log(const double &t);
};

} // namespace robot