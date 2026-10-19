#include "rearrangeDof.h"

void rearrange(Eigen::Matrix<Settings::precision, Eigen::Dynamic, 1uz>& disp
    , Eigen::Matrix<Settings::precision, Eigen::Dynamic, 1uz>& force
    , Eigen::Matrix<int, Eigen::Dynamic, 1uz>& mapDof)
{
    // Creating the mapDof vector based on NaN and !NaN distinction
    std::vector<int> mapDofTemp (mapDof.data(), mapDof.data()+mapDof.rows());
    std::stable_partition(mapDofTemp.begin(), mapDofTemp.end(), [&](int ind)
        ->bool { return std::isnan(disp(ind)); });
    mapDof = Eigen::Map<Eigen::Matrix<int, Eigen::Dynamic, 1uz>>(
        mapDofTemp.data(), mapDofTemp.size());

    // Re-arranging disp vector
    disp = disp(mapDof);

    // Re-arranging force vector
    force = force(mapDof);
}