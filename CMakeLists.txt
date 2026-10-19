cmake_minimum_required(VERSION 3.10.0)
project(implicitDyn2DNonLinear VERSION 0.1.0 LANGUAGES C CXX)

set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

set(MKL_INTERFACE "lp64")
set(MKL_THREADING "intel_thread")
set(MKL_LINK "dynamic")

find_package(MKL CONFIG REQUIRED)

file(GLOB SOURCES
    "src/*.cpp" 
    "src/*.h"     
)

add_executable(${PROJECT_NAME} ${SOURCES})

add_compile_options( 
    -Wall 
    -Wsign-conversion 
)

add_subdirectory(EigenLib/eigen-5.0.0/)

target_link_libraries(${PROJECT_NAME} 
    PRIVATE
        Eigen3::Eigen 
        MKL::MKL  
)

# adding necessary preprocessor macros for Eigen-MKL coupling
target_compile_definitions(${PROJECT_NAME} 
    PRIVATE
        EIGEN_USE_MKL_ALL)
