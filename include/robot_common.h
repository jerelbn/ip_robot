#pragma once

#include "common_cpp/common.h"
#include "common_cpp/logger.h"

namespace robot
{

enum
{
    DX,
    DPSI,
    THETA,
    DTHETA,
    OMEGAL,
    OMEGAR,
    QL,
    QR,
    NUM_STATES
};

enum
{
    VL,
    VR,
    NUM_INPUTS
};

typedef Eigen::Matrix<double, NUM_STATES, 1> xVector;
typedef Eigen::Matrix<double, NUM_INPUTS, 1> uVector;

typedef Eigen::Matrix<double, NUM_STATES, NUM_STATES> xMatrix;
typedef Eigen::Matrix<double, NUM_STATES, NUM_INPUTS> uMatrix;

} // namespace robot