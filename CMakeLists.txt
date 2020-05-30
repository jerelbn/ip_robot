cmake_minimum_required(VERSION 3.10)
project(ip_robot)

# default to release build
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif(NOT CMAKE_BUILD_TYPE)
set(CMAKE_CXX_FLAGS "-std=c++14")

find_package(Eigen3 REQUIRED)
find_package(yaml-cpp REQUIRED)
find_package(Ceres REQUIRED)

add_subdirectory(lib/common_cpp)
add_subdirectory(lib/lin_alg_tools)

add_definitions(-DPARAM_DIR="${CMAKE_CURRENT_SOURCE_DIR}/param")

include_directories(include
    lib/common_cpp/include
    lib/lin_alg_tools/include
    ${EIGEN3_INCLUDE_DIRS}
    ${YAML_CPP_INCLUDE_DIR}
    ${CERES_INCLUDE_DIRS}
)

add_executable(main
    src/main.cpp
    src/robot.cpp
    src/robot_ctrl_lqr.cpp
    src/robot_ctrl_pid.cpp
    )
target_link_libraries(main
    ${CERES_LIBRARIES}
    ${YAML_CPP_LIBRARIES}
    )

