#include "common_cpp/common.h"
#include "common_cpp/progress_bar.h"
#include "robot.h"
#include "robot_ctrl_lqr.h"
#include "robot_ctrl_pid.h"

int main()
{
    int seed;
    double t(0), tf, dt;
    common::getYamlNode("seed", std::string(PARAM_DIR)+"/robot.yaml", seed);
    common::getYamlNode("tf", std::string(PARAM_DIR)+"/robot.yaml", tf);
    common::getYamlNode("dt", std::string(PARAM_DIR)+"/robot.yaml", dt);
    if (seed < 0)
        seed = std::chrono::system_clock::now().time_since_epoch().count();
    std::default_random_engine rng(seed);
    std::srand(seed);

    // Create progress bar
    common::ProgressBar prog_bar;
    prog_bar.init(tf / dt, 40);

    // Create vehicles, controllers, estimators, sensor packages
    robot::Robot robot(std::string(PARAM_DIR)+"/robot.yaml");
    // robot::LQRController robot_ctrl(std::string(PARAM_DIR)+"/robot.yaml");
    robot::PIDController robot_ctrl(std::string(PARAM_DIR)+"/robot.yaml");

    // Main simulation loop
    while (t <= tf+dt)
    {
        // Update vehicles, controllers, sensors, estimators
        robot_ctrl.computeControl(t, robot.x(), 1.0, 0.5);
        robot.propagate(t, robot_ctrl.u());

        // Update time step
        t += dt;
        prog_bar.print(t / dt);
    }
    prog_bar.finished();

    return 0;
}