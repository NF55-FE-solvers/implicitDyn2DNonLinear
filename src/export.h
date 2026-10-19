#ifndef EXPORT_H
#define EXPORT_H

#include <iostream>
#include <fstream>
#include <sstream>
#include "Eigen/Core"
#include "Mesh.h"
#include "settings.h"
#include "material.h"

namespace Export{
    constexpr const int stressComp {3}; // number of stress components
    using DispEigen = Eigen::Vector<Settings::precision, Eigen::Dynamic>;
    using DispMatEigen = Eigen::Matrix<Settings::precision,Eigen::Dynamic,2>;
    using SigmaMat = Eigen::Matrix<Settings::precision,Eigen::Dynamic
        ,stressComp>;
    using NodesCoordEigen = Eigen::Matrix<Settings::precision
        ,Settings::nNodesEl,Settings::dofNode>;
    using dNEigen = Eigen::Matrix<Settings::precision,Settings::nNodesEl
        ,Settings::nDofEl>;
    using InvJEigen = Eigen::Matrix<Settings::precision,Settings::nNodesEl
        ,Settings::nNodesEl>;
    using BEigen = Eigen::Matrix<Settings::precision,3uz,Settings::nDofEl>;

    constexpr const std::array<Settings::precision,Settings::nNodesEl> xiV
        {{-1.0/Settings::sqrt3, 1.0/Settings::sqrt3, 1.0/Settings::sqrt3
        , -1.0/Settings::sqrt3}};
    constexpr const std::array<Settings::precision,Settings::nNodesEl> etaV
        {{-1.0/Settings::sqrt3, -1.0/Settings::sqrt3, 1.0/Settings::sqrt3
        , 1.0/Settings::sqrt3}};
    constexpr const std::array<Settings::precision,Settings::nNodesEl> rV
        {{-Settings::sqrt3, Settings::sqrt3, Settings::sqrt3,-Settings::sqrt3}};
    constexpr const std::array<Settings::precision,Settings::nNodesEl> sV
        {{-Settings::sqrt3, -Settings::sqrt3, Settings::sqrt3,Settings::sqrt3}};
    constexpr const Eigen::Matrix<Settings::precision,3uz,4uz> arrangeMat 
    {
        {1.0, 0.0, 0.0, 0.0},
        {0.0, 0.0, 0.0, 1.0},
        {0.0, 1.0, 1.0, 0.0}
    };

    void exportResultsBinary(const Mesh& mesh, const DispEigen& disp
        , const int counterVTK, const Settings::precision time);
    void computeStressField(const Mesh& mesh, DispMatEigen& dispMat
            , SigmaMat& sigmaM);
    void generateFileName(std::string& fileName, const int counterVTK);

    
    // Little-endian to Big-endian (vtk legacy) conversion, template function
    template <typename T>
    T swapEndian(T val){
        union{ // data structure that share same memory address
            T val; // value stored in function argument val
            uint8_t bytes[sizeof(T)]; // empty C-vector with length equal to 
            // the number of bytes used to describe val 
        } src, dest; // defining 2 variables of type union just defined
        src.val = val;
        // reversing the order of bytes used to describe function argument val
        for (size_t i = 0; i < sizeof(T); ++i) {
            dest.bytes[i] = src.bytes[sizeof(T) - 1 - i];
        }
        return dest.val;
    }

    // Overloaded function to find which floating point type is used
    constexpr bool isDouble(double floatingP){ return true; }
    constexpr bool isDouble(float floatingP){ return false; }
}

#endif