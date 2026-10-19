#ifndef IMPORT_H
#define IMPORT_H

#include <iostream>
#include <fstream>
#include <sstream>
#include <array>
#include <vector>
#include <filesystem>
#include "Eigen/Core"
#include "Mesh.h"
#include "settings.h"
#include "Bc.h"

void importMesh(Mesh& mesh, Bc& bc);

void extractMeshFromAbaqusInp(std::ifstream& file, Mesh& mesh);
void extractMeshFromSalomeDat(std::ifstream& file, Mesh& mesh);

#endif