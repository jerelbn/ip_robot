#pragma once

#include "common_cpp/common.h"
#include "common_cpp/logger.h"
#include "robot_common.h"

namespace robot
{

class Robot
{

public:
    Robot();
    Robot(const std::string &filename);
    ~Robot();

    void load(const std::string &filename);
    void propagate(const double &t, const uVector &u);

    const xVector &x() const { return x_; }

private:
    void f(const xVector &x, const uVector &u, xVector &dx);
    void log(const double &t);

    double t_prev_;
    xVector x_;
    xVector dx_;

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

    common::Logger state_log_;
};

} // namespace robot