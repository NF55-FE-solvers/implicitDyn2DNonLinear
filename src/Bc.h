#ifndef BC_H
#define BC_H

#include <vector>
#include <array>
#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>
#include "Eigen/Core"
#include "settings.h"

class Bc
{
private:
    using sizT = std::size_t;
    using ndofVectf = Eigen::Vector<Settings::precision,Eigen::Dynamic>;
    using ndofVecti = Eigen::Vector<int, Eigen::Dynamic>;
    using sparseMat = Eigen::SparseMatrix<Settings::precision>;
    using eigenIterator=Eigen::SparseMatrix<Settings::precision>::InnerIterator;

    enum BCnames : sizT{
        BCx,
        BCy,
        Fx,
        Fy,
        countBCnames
    };

    std::array<std::vector<int>,countBCnames> bcs {};

public:
    Bc() = default;

    void extractBCSalome();
    void extractBCAbaqus();

    void setBcInit(ndofVectf& disp, ndofVectf& force) const;
    
    void setBcTimeDependent(ndofVectf& dispNewOrForce
        , const Settings::precision time) const;
    
    void embedBCintoKt(sparseMat& ktEff) const;
    void embedBCintoResidual(ndofVectf& residual) const;
};

#endif