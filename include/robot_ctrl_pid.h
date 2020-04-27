#pragma once

#include "common_cpp/common.h"
#include "common_cpp/logger.h"
#include "robot_common.h"

namespace robot
{

class PIDController
{

public:
    PIDController();
    PIDController(const std::string &filename);
    ~PIDController();

    void load(const std::string &filename);
    void computeControl(const double &t, const xVector &x, const double &dx_d, const double &dpsi_d);

    const Eigen::Vector2d &u() const { return u_; }

private:
    bool initialized_;
    double t_prev_;

    uVector u_;

    int update_rate_;
    double L_;
    double r_;
    double kp_theta_;
    double ki_theta_;
    double kd_theta_;
    double tau_theta_;
    double kp_V_;
    double ki_V_;
    double kd_V_;
    double tau_V_;
    
    common::PID<double> pid_theta_;
    common::PID<double> pid_v_;

    common::Logger command_log_;

    void log(const double &t);
};

} // namespace robot