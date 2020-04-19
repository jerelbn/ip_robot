#include "robot_ctrl_lqr.h"

namespace robot
{

LQRController::LQRController() : t_prev_(0), initialized_(false)
{
    u_.setZero();
}

LQRController::LQRController(const std::string &filename) : t_prev_(0), initialized_(false)
{
    u_.setZero();
    load(filename);
}

LQRController::~LQRController() {}

void LQRController::load(const std::string &filename)
{
    // Load parameters
    common::getYamlNode("ctrl_update_rate", filename, update_rate_);
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
    command_log_.open("robot_lqr_command.log");
}

void LQRController::computeControl(const double &t, const xVector &x, const double &dx_d, const double &dpsi_d)
{
    double dt = common::roundDecimal(t - t_prev_, 6);
    if (t == 0 || dt >= 1.0 / update_rate_)
    {
        Eigen::Matrix<double,6,1> x2 = x.segment<6>(THETA);

        // Copy the current state and time
        t_prev_ = t;

        // Jacobians w.r.t. state and input
        numericalAB(x2, u_, A_, B_);

        // Solve CARE
        solver_.solve(P_, A_, B_, Q_, R_);
        K_ = R_inv_ * B_.transpose() * P_;

        // Define reference state and compute control gain
        Eigen::Matrix<double,6,1> x_ref = x.segment<6>(THETA);
        u_ = -K_ * x_ref;
    }

    // Log all data
    log(t);
}

// void LQRController::f(const xVector &x, const uVector &u, xVector &dx)
// {
//     // Constants
//     static double g = common::gravity;
//     static double c0 = mp_ * l_;
//     static double c1 = Jy_ + c0 * l_;
//     static double c2 = c0 * r_;
//     static double c22 = c2 * c2;
//     static double c3 = 2.0 * Jm_ * c1;
//     static double c4 = 2.0 * Jm_ * c2;
//     static double c5 = c0 * g;

//     // Unpack states/inputs for readability
//     double theta = x(THETA);
//     double dtheta = x(DTHETA);
//     double omegal = x(OMEGAL);
//     double omegar = x(OMEGAR);
//     double ql = x(QL);
//     double qr = x(QR);
//     double Vl = u(VL);
//     double Vr = u(VR);

//     // Common calcs
//     double sin_theta = sin(theta);
//     double cos_theta = cos(theta);
//     double cos2_theta = cos_theta * cos_theta;
//     double dtheta2 = dtheta * dtheta;
//     double v0 = c22 * cos2_theta;
//     double v1 = (c1 * dtheta2 - c5) * sin_theta - bp_ * dtheta;
//     double v2 = c4 * cos_theta * v1;
//     double v3 = Km_ * (ql + qr);
//     double v4 = bm_ * (omegal + omegar);
//     double v5 = Jm_ * c2 * cos_theta;
//     double v6 = c2 * cos_theta;
//     double denom = Jm_ * (Jm_ * c1 - v0);

//     // Equations of motion
//     dx(DX) = -(r_ / 2.0) * ((c3 + v0 * (Jm_ - 1.0)) * (v3 - v4) + v2) / denom;
//     dx(DPSI) = (r_ / L_) * ((c3 - v0 * (Jm_ + 1.0)) * (Km_ * (qr - ql) + bm_ * (omegal - omegar))) / denom;
//     dx(THETA) = dtheta;
//     dx(DTHETA) = (c5 * sin_theta - bp_ * dtheta) / c1 - v6 / (2.0 * c1) * ((c3 + v0 * (Jm_ - 1.0)) * (v3 - v4) + v2) / denom;
//     dx(OMEGAL) = ((c3 - v0) * (Km_ * ql - bm_ * omegal) + v5 * (v6 * (Km_ * qr - bm_ * omegar) + v1)) / denom;
//     dx(OMEGAR) = ((c3 - v0) * (Km_ * qr - bm_ * omegar) + v5 * (v6 * (Km_ * ql - bm_ * omegal) + v1)) / denom;
//     dx(QL) = (Vl - Rm_ * ql - Km_ * omegal) / Lm_;
//     dx(QR) = (Vr - Rm_ * qr - Km_ * omegar) / Lm_;
// }

void LQRController::f(const Eigen::Matrix<double,6,1> &x, const uVector &u, Eigen::Matrix<double,6,1> &dx)
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
}

void LQRController::numericalAB(const Eigen::Matrix<double,6,1> &x, const uVector &u, Eigen::Matrix<double,6,6> &A, Eigen::Matrix<double,6,2> &B)
{
    static const double eps = 1e-5;
    static const Eigen::Matrix<double,6,6> Ix = Eigen::Matrix<double,6,6>::Identity();
    static const Eigen::Matrix2d Iu = Eigen::Matrix2d::Identity();
    static Eigen::Matrix<double,6,1> xp, xm, dxp, dxm;
    static uVector up, um;
    for (int i = 0; i < A.cols(); ++i)
    {
        xp = x + eps * Ix.col(i);
        xm = x + -eps * Ix.col(i);

        f(xp, u, dxp);
        f(xm, u, dxm);

        A.col(i) = (dxp - dxm) / (2.0 * eps);
    }
    for (int i = 0; i < B.cols(); ++i)
    {
        up = u + eps * Iu.col(i);
        um = u + -eps * Iu.col(i);

        f(x, up, dxp);
        f(x, um, dxm);

        B.col(i) = (dxp - dxm) / (2.0 * eps);
    }
}

void LQRController::log(const double &t)
{
    // Write data to binary files and plot in another program
    command_log_.log(t);
    command_log_.logMatrix(u_);
}

} // namespace robot