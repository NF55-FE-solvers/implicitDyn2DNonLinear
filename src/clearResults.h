#ifndef CLEAR_RESULTS_H
#define CLEAR_RESULTS_H

#include <filesystem>
#include "settings.h"

void clearPriorResults(){
    // current working directory has been already checked against the expected 
    // one in import.cpp
    Settings::fPath resultsPathName {Settings::resultsDirName};
    if (!resultsPathName.empty()){
        for(auto& file : std::filesystem::directory_iterator {resultsPathName}){
            if (file.path().extension() == ".vtk"){
                std::filesystem::remove(file.path());
            }
        }
    }
}

#endif