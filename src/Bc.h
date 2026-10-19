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
    using ndofVectf = Eigen::Vector<Settings::precision,Eigen::Dynamic>;
    using ndofVecti = Eigen::Vector<int, Eigen::Dynamic>;

    enum BCnames : std::size_t{
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
    
    void setBcTimeDependent(ndofVectf& dispNewOrForce, const ndofVecti& mapDof
        , const Settings::precision time, const int mapDofBlockSize);
};

#endif