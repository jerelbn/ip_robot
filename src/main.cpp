#include "common_cpp/common.h"
#include "common_cpp/progress_bar.h"
#include "robot.h"

int main()
{
    int seed;
    double t(0), tf, dt;
    common::getYamlNode("seed", "../param/robot.yaml", seed);
    common::getYamlNode("tf", "../param/robot.yaml", tf);
    common::getYamlNode("dt", "../param/robot.yaml", dt);
    if (seed < 0)
        seed = std::chrono::system_clock::now().time_since_epoch().count();
    std::default_random_engine rng(seed);
    std::srand(seed);

    // Create progress bar
    common::ProgressBar prog_bar;
    prog_bar.init(tf / dt, 40);

    // Create vehicles, controllers, estimators, sensor packages
    robot::Robot robot("../param/robot.yaml");
    //   robot::Controller robot_ctrl("../params/robot.yaml", rng);

    // Main simulation loop
    while (t <= tf+dt)
    {
        // Update vehicles, controllers, sensors, estimators
        Eigen::Vector2d u;
        u << 0, 0;
        robot.propagate(t, u);
        // robot_ctrl.computeControl(t, robot.x());
        // robot.updateAccelerations(robot_ctrl.u_);

        // Update time step
        t += dt;
        prog_bar.print(t / dt);
    }
    prog_bar.finished();

    return 0;
}