#include "robot_ctrl_lqr.h"

namespace robot
{

LQRController::LQRController() : t_prev_(0), initialized_(false)
{
    u_.setZero();
    dx_int_ = 0;
}

LQRController::LQRController(const std::string &filename) : t_prev_(0), initialized_(false)
{
    u_.setZero();
    dx_int_ = 0;
    load(filename);
}

LQRController::~LQRController() {}

void LQRController::load(const std::string &filename)
{
    // Load parameters
    common::getYamlNode("ctrl_update_rate", filename, update_rate_);
    common::getYamlNode("max_voltage", filename, max_voltage_);
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
    common::getYamlNode("bm", filename, bm_);
    common::getYamlNode("bp", filename, bp_);
    common::getYamlEigenDiag("ctrl_Q", filename, Q_);
    common::getYamlEigenDiag("ctrl_R", filename, R_);
    R_inv_ = R_.inverse();

    // Initialize logger
    command_log_.open("/tmp/robot_lqr_command.log");
}

void LQRController::computeControl(const double &t, const xVector &x, const double &dx_d, const double &dpsi_d)
{
    double dt = common::roundDecimal(t - t_prev_, 6);
    if (t == 0 || dt >= 1.0 / update_rate_)
    {
        // Declarations
        Eigen::Matrix<double, NS, 1> xref, xtilde;

        // Copy the current state and time
        t_prev_ = t;

        // Extract desired wheel rates
        double omegal_d = -(dx_d + L_ / 2.0 * dpsi_d) / r_;
        double omegar_d = -(dx_d - L_ / 2.0 * dpsi_d) / r_;

        // Create LQR state
        Eigen::Matrix<double, NS, 1> x2;
        x2(0) = x(THETA);
        x2(1) = x(DTHETA);
        x2(2) = x(OMEGAL);
        x2(3) = x(OMEGAR);
        x2(4) = x(QL);
        x2(5) = x(QR);
        x2(6) = dx_int_;

        // Reference state
        xref.setZero();
        xref(2) = omegal_d;
        xref(3) = omegar_d;

        // Error state
        xtilde = xref - x2;
        // xtilde(2) = common::saturate(xtilde(2), 0.01, -0.01);
        // xtilde(3) = common::saturate(xtilde(3), 0.01, -0.01);

        // Jacobians w.r.t. state and input
        numericalAB(x2, xref, u_);

        // Solve CARE
        solver_.solve(P_, A_, B_, Q_, R_);
        K_ = R_inv_ * B_.transpose() * P_;

        // Compute control gain
        u_ = -K_ * xtilde;
        u_(VL) = common::saturate(u_(VL), max_voltage_, -max_voltage_);
        u_(VR) = common::saturate(u_(VR), max_voltage_, -max_voltage_);

        // Integrate velocity error
        if (u_(VL) != max_voltage_ && u_(VR) != max_voltage_)
            dx_int_ += 10.0 * (-r_ / 2.0 * (x(OMEGAL) + x(OMEGAR)) - dx_d_) * dt;
    }

    // Log all data
    log(t);
}

void LQRController::f(const Eigen::Matrix<double, NS, 1> &x, const uVector &u, Eigen::Matrix<double, NS, 1> &dx) const
{
    // Constants
    static double g = common::gravity;
    static double c0 = mp_ * l_;
    static double c1 = Jy_ + c0 * l_;
    static double c2 = c0 * r_;
    static double c22 = c2 * c2;
    static double c3 = 2.0 * Jm_ * c1;
    static double c4 = 2.0 * Jm_ * c2;
    static double c5 = c0 * g;

    // Unpack states/inputs for readability
    double theta = x(0);
    double dtheta = x(1);
    double omegal = x(2);
    double omegar = x(3);
    double ql = x(4);
    double qr = x(5);
    double Vl = u(VL);
    double Vr = u(VR);

    // Common calcs
    double sin_theta = sin(theta);
    double cos_theta = cos(theta);
    double cos2_theta = cos_theta * cos_theta;
    double dtheta2 = dtheta * dtheta;
    double v0 = c22 * cos2_theta;
    double v1 = (c1 * dtheta2 - c5) * sin_theta - bp_ * dtheta;
    double v2 = c4 * cos_theta * v1;
    double v3 = Km_ * (ql + qr);
    double v4 = bm_ * (omegal + omegar);
    double v5 = Jm_ * c2 * cos_theta;
    double v6 = c2 * cos_theta;
    double denom = Jm_ * (Jm_ * c1 - v0);

    // Equations of motion
    dx(0) = dtheta;
    dx(1) = (c5 * sin_theta - bp_ * dtheta) / c1 - v6 / (2.0 * c1) * ((c3 + v0 * (Jm_ - 1.0)) * (v3 - v4) + v2) / denom;
    dx(2) = ((c3 - v0) * (Km_ * ql - bm_ * omegal) + v5 * (v6 * (Km_ * qr - bm_ * omegar) + v1)) / denom;
    dx(3) = ((c3 - v0) * (Km_ * qr - bm_ * omegar) + v5 * (v6 * (Km_ * ql - bm_ * omegal) + v1)) / denom;
    dx(4) = (Vl - Rm_ * ql - Km_ * omegal) / Lm_;
    dx(5) = (Vr - Rm_ * qr - Km_ * omegar) / Lm_;
    dx(6) = -r_ / 2.0 * (omegal + omegar) - dx_d_;
}

void LQRController::numericalAB(const Eigen::Matrix<double, NS, 1> &x, const Eigen::Matrix<double, NS, 1> &xref, const uVector &u)
{
    static const double eps = 1e-5;
    static const Eigen::Matrix<double, NS, NS> Ix = Eigen::Matrix<double, NS, NS>::Identity();
    static const Eigen::Matrix2d Iu = Eigen::Matrix2d::Identity();
    static Eigen::Matrix<double, NS, 1> xtildep, xtildem, dxtildep, dxtildem;
    static uVector up, um;

    // Error state
    Eigen::Matrix<double, NS, 1> xtilde = xref - x;

    for (int i = 0; i < A_.cols(); ++i)
    {
        // Poke the error state
        xtildep = xtilde + eps * Ix.col(i);
        xtildem = xtilde + -eps * Ix.col(i);

        // Error state derivatives
        ftilde(eps, xref, xtildep, u, dxtildep);
        ftilde(eps, xref, xtildem, u, dxtildem);

        // Derivative of dxtilde w.r.t. xtilde
        A_.col(i) = (dxtildep - dxtildem) / (2.0 * eps);
    }

    for (int i = 0; i < B_.cols(); ++i)
    {
        // Poke the command vector
        uVector up = u + eps * Iu.col(i);
        uVector um = u + -eps * Iu.col(i);

        // Error state derivatives
        ftilde(eps, xref, xtilde, up, dxtildep);
        ftilde(eps, xref, xtilde, um, dxtildem);

        // Derivative of dxtilde w.r.t. u_ref
        B_.col(i) = (dxtildep - dxtildem) / (2.0 * eps);
    }
}

void LQRController::ftilde(const double &dt, const Eigen::Matrix<double,NS,1> &xref, const Eigen::Matrix<double,NS,1> &xtilde,
                           const uVector &u, Eigen::Matrix<double,NS,1> &dxtilde) const
{
    // Declarations
    static Eigen::Matrix<double, NS, 1> x, dx, xp, xm, xrefp, xrefm, xtildep, xtildem;

    // 'True state'
    x = xref + -xtilde;

    // Derivative at current time
    f(x, u, dx);

    // Future and previous states
    xp = x + dx * dt;
    xm = x + -dx * dt;

    // Future and previous error states
    xtildep = xref - xp;
    xtildem = xref - xm;

    // Error state derivative
    dxtilde = (xtildep - xtildem) / (2.0 * dt);
}

void LQRController::log(const double &t)
{
    // Write data to binary files and plot in another program
    command_log_.log(t);
    command_log_.logMatrix(u_);
}

} // namespace robot