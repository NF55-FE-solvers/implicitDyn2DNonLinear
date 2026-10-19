#ifndef MATERIAL_H
#define MATERIAL_H

#include <Eigen/Core>
#include "settings.h"

namespace Material{

    constexpr const Settings::precision rho {7.8e-9}; // Elements/part density
    constexpr const Settings::precision t {2.0}; // Elements/part thickness
    constexpr const Settings::precision E {200000.0}; // Young's modulus
    constexpr const Settings::precision nu {0.33}; // Poisson's ratio
    constexpr const Settings::precision G {E/(2*(1+nu))}; // Shear modulus

    #define PLANE_STRESS 1

    #if PLANE_STRESS
        constexpr const Settings::precision nu2 {nu*nu};

        constexpr const Eigen::Matrix<Settings::precision,3,3> C
        {
            {E/(1-nu2),         (E*nu)/(1-nu2),     0.0},
            {(E*nu)/(1-nu2),    E/(1-nu2),          0.0},
            {0.0,               0.0,                G}
        };
    #else // PLANE STRAIN CASE
        constexpr const Settings::precision const4Ccalc {E/((1+nu)*(1-2*nu))};

        constexpr const Eigen::Matrix<Settings::precision,3,3> C
        {
            {const4Ccalc*(1-nu),  const4Ccalc*nu,         0.0},
            {const4Ccalc*nu,      const4Ccalc*(1-nu),     0.0},
            {0.0,                 0.0,                      G}
        };
    #endif

}

#endif