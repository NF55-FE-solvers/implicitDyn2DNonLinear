#ifndef REARRANGE_H
#define REARRANGE_H

#include <vector>
#include <algorithm>
#include <cmath>
#include "Eigen/Core"
#include "settings.h"

void rearrange(Eigen::Matrix<Settings::precision, Eigen::Dynamic, 1uz>& disp
    , Eigen::Matrix<Settings::precision, Eigen::Dynamic, 1uz>& force
    , Eigen::Matrix<int, Eigen::Dynamic, 1uz>& mapDof);

#endif