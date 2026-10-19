#ifndef MESH_H
#define MESH_H

#include "Eigen/Core"
#include "settings.h"

struct Mesh
{
    int nDof {};
    int nElements {};
    
    // nodal coordinates
    Eigen::Matrix<Settings::precision,Eigen::Dynamic,2,Eigen::RowMajor> 
        nodes;
    // connectivity table
    Eigen::Matrix<int,Eigen::Dynamic,4,Eigen::RowMajor> elements;
};


#endif