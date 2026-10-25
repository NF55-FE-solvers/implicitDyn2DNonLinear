#ifndef INTERNAL_FORCE_AND_KT_H
#define INTERNAL_FORCE_AND_KT_H

#include "Eigen/Core"
#include "settings.h"
#include "material.h"
#include "Mesh.h"
#include "elMassMatrix.h"

namespace FintKtInt{
    using SizT = std::size_t;
    using nDofVectf = Eigen::Vector<Settings::precision,Eigen::Dynamic>;
    using nDofElVectf = Eigen::Vector<Settings::precision,Settings::nDofEl>;
    using sparseMat = Eigen::SparseMatrix<Settings::precision>;
    using elMat = Eigen::Matrix<Settings::precision,Settings::nDofEl
        ,Settings::nDofEl>;
    using dNEigen = Eigen::Matrix<Settings::precision,Settings::nNodesEl
        ,Settings::nDofEl>;
    using invJEigen = Eigen::Matrix<Settings::precision,Settings::nNodesEl
        ,Settings::nNodesEl>;
    using B0Eigen = Eigen::Matrix<Settings::precision,Settings::nNodesEl
        ,Settings::nDofEl>;
    using voigtTensorFull = Eigen::Vector<Settings::precision
        ,Settings::dofNode*Settings::dofNode>;
    using voigtTensor = Eigen::Vector<Settings::precision,3>;
    using tensor = Eigen::Matrix<Settings::precision,Settings::dofNode
        ,Settings::dofNode>;
    using matMultipliedByVel = Eigen::Matrix<Settings::precision,3
        ,Settings::nDofEl>;
    using matMultipliedByVelFull = Eigen::Matrix<Settings::precision,4
        ,Settings::nDofEl>;

    void compute(const Mesh& mesh, nDofVectf& dispNew, nDofVectf& fInt
        , sparseMat& ktInt);
}


#endif