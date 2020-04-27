#include "robot_ctrl_pid.h"

namespace robot
{

PIDController::PIDController() : t_prev_(0), initialized_(false)
{
    u_.setZero();
}

PIDController::PIDController(const std::string &filename) : t_prev_(0), initialized_(false)
{
    u_.setZero();
    load(filename);
}

PIDController::~PIDController() {}

void PIDController::load(const std::string &filename)
{
    // Load parameters
    double max_voltage;
    common::getYamlNode("ctrl_update_rate", filename, update_rate_);
    common::getYamlNode("max_voltage", filename, max_voltage);
    common::getYamlNode("L", filename, L_);
    common::getYamlNode("r", filename, r_);
    common::getYamlNode("kp_theta", filename, kp_theta_);
    common::getYamlNode("ki_theta", filename, ki_theta_);
    common::getYamlNode("kd_theta", filename, kd_theta_);
    common::getYamlNode("tau_theta", filename, tau_theta_);
    common::getYamlNode("kp_V", filename, kp_V_);
    common::getYamlNode("ki_V", filename, ki_V_);
    common::getYamlNode("kd_V", filename, kd_V_);
    common::getYamlNode("tau_V", filename, tau_V_);

    pid_theta_.init(kp_theta_, ki_theta_, kd_theta_, M_PI/4, -M_PI/4, tau_theta_);
    pid_v_.init(kp_V_, ki_V_, kd_V_, max_voltage, -max_voltage, tau_V_);

    // Initialize logger
    command_log_.open("/tmp/robot_pid_command.log");
}

void PIDController::computeControl(const double &t, const xVector &x, const double &dx_d, const double &dpsi_d)
{
    double dt = common::roundDecimal(t - t_prev_, 6);
    if (t == 0 || dt >= 1.0 / update_rate_)
    {
        // Store time for next iteration
        t_prev_ = t;

        // Unpack relevant states
        double theta = x(THETA);
        double dtheta = x(DTHETA);
        double omegal = x(OMEGAL);
        double omegar = x(OMEGAR);

        // Reference wheel velocities from desired forward/heading velocities
        double omegal_d = -(dx_d + L_ / 2.0 * dpsi_d) / r_;
        double omegar_d = -(dx_d - L_ / 2.0 * dpsi_d) / r_;
        double omega_nom = (omegal_d + omegar_d) / 2.0;

        // Outer loop velocity control to produce desired pendulum angle
        // Inner loop pendulum angle control to produce nominal voltage
        double dx = -r_ / 2.0 * (omegal + omegar);
        double theta_d = pid_theta_.run(dt, dx, dx_d, false);
        double v_nom = pid_v_.run(dt, -theta, theta_d, false, -dtheta);

        // Compute control by weighting voltage difference to each wheel
        u_(VL) = omegal_d / omega_nom * v_nom;
        u_(VR) = omegar_d / omega_nom * v_nom;
    }

    // Log all data
    log(t);
}

void PIDController::log(const double &t)
{
    // Write data to binary files and plot in another program
    command_log_.log(t);
    command_log_.logMatrix(u_);
}

} // namespace robot