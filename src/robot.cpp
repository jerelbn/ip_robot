#include "robot.h"

namespace robot
{

Robot::Robot() : t_prev_(0) {}

Robot::Robot(const std::string &filename) : t_prev_(0)
{
    load(filename);
}

Robot::~Robot() {}

void Robot::load(const std::string &filename)
{
    // Load all Robot parameters
    common::getYamlEigen("x0", filename, x_);
    common::getYamlNode("mc", filename, mc_);
    common::getYamlNode("mp", filename, mp_);
    common::getYamlNode("Jy", filename, Jy_);
    common::getYamlNode("Jz", filename, Jz_);
    common::getYamlNode("Jp", filename, Jp_);
    common::getYamlNode("Jm", filename, Jm_);
    common::getYamlNode("Km", filename, Km_);
    common::getYamlNode("Lm", filename, Lm_);
    common::getYamlNode("Rm", filename, Rm_);
    common::getYamlNode("L", filename, L_);
    common::getYamlNode("l", filename, l_);
    common::getYamlNode("r", filename, r_);
    common::getYamlNode("b", filename, b_);

    // Initialize loggers and log initial data
    std::string logname_true_state;
    common::getYamlNode("logname_true_state", filename, logname_true_state);
    state_log_.open(logname_true_state);
}

void Robot::propagate(const double &t, const uVector &u)
{
    double dt = t - t_prev_;
    t_prev_ = t;

    if (t > 0)
    {
        // 4th order Runge-Kutta integration
        common::rk4<double, NUM_STATES, NUM_INPUTS>(std::bind(&Robot::f, this,
                                                              std::placeholders::_1,
                                                              std::placeholders::_2,
                                                              std::placeholders::_3),
                                                    dt, x_, u, dx_);
        x_ += dx_;
    }

    log(t);
}

void Robot::f(const xVector &x, const uVector &u, xVector &dx)
{
    // Constants
    static double mp2 = mp_*mp_;
    static double l2 = l_*l_;
    static double M = mc_ + mp_;
    static double Jmp = Jy_ + mp_*l2;
    static double g = common::gravity;

    // Unpack states/inputs for readability
    double theta = x_(THETA);
    double dtheta = x_(DTHETA);
    double omegal = x_(OMEGAL);
    double omegar = x_(OMEGAR);
    double ql = x_(QL);
    double qr = x_(QR);
    double Vl = u(VL);
    double Vr = u(VR);

    // Common calcs
    double sin_theta = sin(theta);
    double cos_theta = cos(theta);
    double denom = r_*(M*Jmp - mp2*l2*cos_theta*cos_theta);
    double qpq = ql + qr;
    double dtheta2 = dtheta*dtheta;
    
    // Equations of motion
    dx(DX) = ((mp_*g*l_*cos_theta - dtheta2*Jmp)*r_*mp_*l_*sin_theta - Km_*qpq*Jmp)/denom;
    dx(DPSI) = L_*Km_*(qr - ql)/(r_*(Jz_ + Jp_*theta));
    dx(THETA) = dtheta;
    dx(DTHETA) = (mp_*l_*(r_*M*g*sin_theta - (Km_*qpq + r_*mp_*l_*dtheta2*sin_theta)*cos_theta))/denom;
    dx(OMEGAL) = (Km_*ql - b_*omegal)/Jm_;
    dx(OMEGAR) = (Km_*qr - b_*omegar)/Jm_;
    dx(QL) = (Vl - Rm_*ql - Km_*omegal)/Lm_;
    dx(QR) = (Vr - Rm_*qr - Km_*omegar)/Lm_;
}

void Robot::log(const double &t)
{
    state_log_.log(t);
    state_log_.logMatrix(x_);
}

} // namespace robot