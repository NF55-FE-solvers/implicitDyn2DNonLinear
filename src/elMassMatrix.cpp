#include "elMassMatrix.h"

namespace mMat{
    void computeElM(const NodesCoordEigen& coords, NEigen& N, MelEigen& mEl)
    {
        mEl.setZero();
        Settings::precision dxdxi   {};
        Settings::precision dydxi   {};
        Settings::precision dxdeta  {};
        Settings::precision dydeta  {};
        Settings::precision detJ    {};
        Settings::precision scalarValueM {};
        for (const Settings::precision xi : gaussianPts){
            for (const Settings::precision eta : gaussianPts){
                N << (1.0-xi)*(1.0-eta),0.0,(1.0+xi)*(1.0-eta),0.0
                    ,(1.0+xi)*(1.0+eta),0.0,(1.0-xi)*(1+eta),0.0
                    ,0.0,(1.0-xi)*(1.0-eta),0.0,(1.0+xi)*(1.0-eta)
                    ,0.0,(1.0+xi)*(1.0+eta),0.0,(1.0-xi)*(1+eta);
                N *= 0.25;
                dxdxi = 0.25*(-coords(0,0)+coords(1,0)+coords(2,0)-coords(3,0)
                    +eta*(coords(0,0)-coords(1,0)+coords(2,0)-coords(3,0)));
                dydxi = 0.25*(-coords(0,1)+coords(1,1)+coords(2,1)-coords(3,1)
                    +eta*(coords(0,1)-coords(1,1)+coords(2,1)-coords(3,1)));
                dxdeta = 0.25*(-coords(0,0)-coords(1,0)+coords(2,0)+coords(3,0)
                    +xi*(coords(0,0)-coords(1,0)+coords(2,0)-coords(3,0)));
                dydeta = 0.25*(-coords(0,1)-coords(1,1)+coords(2,1)+coords(3,1)
                    +xi*(coords(0,1)-coords(1,1)+coords(2,1)-coords(3,1)));
                detJ = dxdxi*dydeta - dydxi*dxdeta;
                scalarValueM = Material::t*detJ*Material::rho; // gauss w = 1
                mEl.noalias() += scalarValueM * (N.transpose() * N);
            }
        }
    }
} 