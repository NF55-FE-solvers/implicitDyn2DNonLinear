#ifndef EL_MATRICES_H
#define EL_MATRICES_H

#include <array>
#include "Eigen/Core"
#include "settings.h"
#include "material.h"

namespace kmMat{
    using KMelEigen = Eigen::Matrix<Settings::precision,Settings::nDofEl
        ,Settings::nDofEl>;
    using NodesCoordEigen = Eigen::Matrix<Settings::precision
        ,Settings::nNodesEl,Settings::dofNode>;
    using NEigen = Eigen::Matrix<Settings::precision,Settings::dofNode
        ,Settings::nDofEl>;
    using dNEigen = Eigen::Matrix<Settings::precision,Settings::nNodesEl
        ,Settings::nDofEl>;
    using invJEigen = Eigen::Matrix<Settings::precision,Settings::nNodesEl
        ,Settings::nNodesEl>;
    using BEigen = Eigen::Matrix<Settings::precision,3uz,Settings::nDofEl>;

    constexpr const std::array<Settings::precision,2uz> gaussianPts 
        {{-1.0/Settings::sqrt3, 1.0/Settings::sqrt3}};

    constexpr const Eigen::Matrix<Settings::precision,3uz,4uz> arrangeMat 
    {
        {1.0, 0.0, 0.0, 0.0},
        {0.0, 0.0, 0.0, 1.0},
        {0.0, 1.0, 1.0, 0.0}
    };

    void computeElKM(const NodesCoordEigen& coords, NEigen& N, dNEigen& dN
        , invJEigen& invJ, BEigen& bEl, KMelEigen& kEl, KMelEigen& mEl);
}


#endif