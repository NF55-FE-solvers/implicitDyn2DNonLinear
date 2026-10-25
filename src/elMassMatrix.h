#ifndef EL_MASS_MATRIX_H
#define EL_MASS_MATRIX_H

#include <array>
#include "Eigen/Core"
#include "settings.h"
#include "material.h"

namespace mMat{
    using MelEigen = Eigen::Matrix<Settings::precision,Settings::nDofEl
        ,Settings::nDofEl>;
    using NodesCoordEigen = Eigen::Matrix<Settings::precision
        ,Settings::nNodesEl,Settings::dofNode>;
    using NEigen = Eigen::Matrix<Settings::precision,Settings::dofNode
        ,Settings::nDofEl>;

    constexpr const std::array<Settings::precision,2uz> gaussianPts 
        {{-1.0/Settings::sqrt3, 1.0/Settings::sqrt3}};

    void computeElM(const NodesCoordEigen& coords, NEigen& N, MelEigen& mEl);
}

#endif