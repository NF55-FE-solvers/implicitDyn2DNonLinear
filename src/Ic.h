#ifndef IC_H
#define IC_H

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include "Eigen/Eigen"
#include "settings.h"

namespace Ic{
    using ndofVectf = Eigen::Vector<Settings::precision, Eigen::Dynamic>;
    using ndofVecti = Eigen::Vector<int,Eigen::Dynamic>;
    using sparseMat = Eigen::SparseMatrix<Settings::precision>;
    using namespace std::string_literals;

    constexpr void usedPrecision(double val, std::string& str){ str="double"s;}
    constexpr void usedPrecision(float val, std::string& str){ str="float"s; }

    void setIC(ndofVectf& dispOld, ndofVectf& velOld, ndofVectf& accOld
        , ndofVectf& force, sparseMat& mG, sparseMat& kG, int nUnkDisp
        , ndofVecti& mapDof)
    {
        std::cout << "\nEmbedding initial conditions @t=" << Settings::t0 
            << "s ...\n\n";
        if (Settings::importICFromStatic){
            std::ifstream inFile {Settings::dispFromStaticPath
                , std::ios::binary};

            std::string tempStr {""s};
            std::stringstream str2int;
            std::getline(inFile, tempStr);
            str2int << tempStr;
            int dispVectSize {};
            str2int >> dispVectSize;
            if (dispVectSize != dispOld.size()) throw "The imported nodal displacement vector has an incompatible size - Ic.h";
            tempStr = ""s;
            std::getline(inFile, tempStr);
            std::string currentPrecision {""s};
            usedPrecision(Settings::sqrt3, currentPrecision);
            if (tempStr != currentPrecision) throw "Incompatible precision between static and dynamic solver - Ic.h";

            inFile.read(reinterpret_cast<char *>(dispOld.data())
                , dispOld.size() * sizeof(Settings::precision));
            ndofVectf tempVect = dispOld(mapDof); // due to aliasing
            dispOld = tempVect;
            velOld = ndofVectf::Zero(dispOld.size());
            accOld = ndofVectf::Zero(dispOld.size()); // just initialized
            Settings::solverType solver;
            solver.compute(mG.topLeftCorner(nUnkDisp,nUnkDisp));
            if(solver.info()!=Eigen::Success){
                throw "Linear solver failed decomposition.";
            }
            ndofVectf rhsEq = force.head(nUnkDisp);
            rhsEq -= kG.topLeftCorner(nUnkDisp,nUnkDisp)*dispOld.head(nUnkDisp);
            accOld.head(nUnkDisp) = solver.solve(rhsEq);
            if(solver.info()!=Eigen::Success){
                throw "Linear solver failed solving(1).";
            }
        }else{
            dispOld.head(nUnkDisp) = ndofVectf::Zero(nUnkDisp);
            velOld = ndofVectf::Zero(dispOld.size());
            accOld = ndofVectf::Zero(dispOld.size());
        }
            
    }
}

#endif