#include "internalForceAndKt.h"

namespace FintKtInt{
    void compute(const Mesh& mesh, nDofVectf& dispNew, nDofVectf& fInt
        , sparseMat& ktInt)
    {
        Settings::precision dXdxi   {};
        Settings::precision dYdxi   {};
        Settings::precision dXdeta  {};
        Settings::precision dYdeta  {};
        Settings::precision detJ    {};
        Settings::precision invDetJ {};
        dNEigen dN;
        invJEigen invJ;
        B0Eigen b0mat;
        voigtTensorFull fVect;
        tensor fTensor; // [F]
        tensor eTensor; // [E]
        voigtTensor eVect;
        tensor pk2Tensor; // 2nd Piola-Kirchhoff stress tensor
        voigtTensor pk2Vect;
        tensor pTensor; // Nominal stress tensor
        voigtTensorFull pVect;
        nDofElVectf b0p;
        Settings::precision scalar4GaussInt {};
        matMultipliedByVel eDotStar;
        matMultipliedByVel pk2DotStar;
        matMultipliedByVelFull pDotStar;
        matMultipliedByVelFull fDotStar;
        std::vector<Eigen::Triplet<Settings::precision>> tripletListKtInt;
        tripletListKtInt.reserve(mesh.nElements*Settings::nDofKel);
        elMat ktMat, ktGeo;
        mMat::NodesCoordEigen coordsL;
        nDofElVectf coordsE;
        Eigen::Matrix<int,Settings::dofNode,Settings::nNodesEl> 
            dofEliMat;
        int elMatContiguousIndx {};

        for (SizT e=0uz; e<static_cast<SizT>(mesh.nElements); ++e){
            ktMat.setZero();
            ktGeo.setZero();
            // Lagrangian coordinates (upper case X and Y in the notes)
            coordsL = mesh.nodes(mesh.elements.row(e),Eigen::placeholders::all);
            dofEliMat.row(0) = mesh.elements.row(e)*Settings::dofNode;
            dofEliMat.row(1) = mesh.elements.row(e)*Settings::dofNode 
                + Eigen::RowVector<int,Settings::nNodesEl>::Ones();
            const auto dofEli = dofEliMat.reshaped(); // view
            // Eulerian coordinates (lower case x and y in the notes)
            coordsE = coordsL.transpose().reshaped() + dispNew(dofEli);
            
            for (const auto xi : mMat::gaussianPts){
                for (const auto eta : mMat::gaussianPts){
                    dN <<(-1.0+eta),0.0,(1.0-eta),0.0,(1.0+eta),0.0,(-1.0-eta)
                        ,0.0,(-1.0+xi),0.0,(-1.0-xi),0.0,(1.0+xi),0.0,(1.0-xi)
                        ,0.0,0.0,(-1.0+eta),0.0,(1.0-eta),0.0,(1.0+eta),0.0
                        ,(-1.0-eta),0.0,(-1.0+xi),0.0,(-1.0-xi),0.0,(1.0+xi)
                        ,0.0,(1.0-xi);
                    dN *= 0.25;
                    dXdxi = 0.25*(-coordsL(0,0)+coordsL(1,0)+coordsL(2,0)
                        -coordsL(3,0)+eta*(coordsL(0,0)-coordsL(1,0)
                        +coordsL(2,0)-coordsL(3,0)));
                    dYdxi = 0.25*(-coordsL(0,1)+coordsL(1,1)+coordsL(2,1)
                        -coordsL(3,1)+eta*(coordsL(0,1)-coordsL(1,1)
                        +coordsL(2,1)-coordsL(3,1)));
                    dXdeta = 0.25*(-coordsL(0,0)-coordsL(1,0)+coordsL(2,0)
                        +coordsL(3,0)+xi*(coordsL(0,0)-coordsL(1,0)+coordsL(2,0)
                        -coordsL(3,0)));
                    dYdeta = 0.25*(-coordsL(0,1)-coordsL(1,1)+coordsL(2,1)
                        +coordsL(3,1)+xi*(coordsL(0,1)-coordsL(1,1)+coordsL(2,1)
                        -coordsL(3,1)));
                    detJ = dXdxi*dYdeta - dYdxi*dXdeta;
                    invDetJ = 1.0/detJ;
                    invJ << dYdeta, -dYdxi, 0, 0,
                            -dXdeta, dXdxi, 0, 0,
                            0, 0, dYdeta, -dYdxi,
                            0, 0, -dXdeta, dXdxi;
                    invJ *= invDetJ;
                    b0mat.noalias() = invJ * dN;
                    fVect.noalias() = b0mat * coordsE;
                    fTensor << fVect(0),fVect(1),fVect(2),fVect(3);
                    eTensor.noalias() = 0.5 * fTensor.transpose() * fTensor;
                    eTensor.noalias() -= 0.5 
                        * tensor::Identity(Settings::dofNode,Settings::dofNode);
                    eVect << eTensor(0,0), eTensor(1,1), 2.0*eTensor(0,1);
                    pk2Vect.noalias() = Material::C * eVect;
                    pk2Tensor << pk2Vect(0), pk2Vect(2), pk2Vect(2), pk2Vect(1);
                    pTensor.noalias() = fTensor * pk2Tensor;
                    pVect <<pTensor(0,0),pTensor(0,1),pTensor(1,0),pTensor(1,1);
                    scalar4GaussInt = Material::t * detJ; // gauss weights = 1
                    b0p.noalias() = b0mat.transpose() * pVect;
                    // computing F_int for this element
                    for (int i=0; i<Settings::nDofEl; ++i){
                        fInt(dofEli(i)) += scalar4GaussInt * b0p(i);
                    }
                    // computing material tangent stiffness matrix
                    eDotStar.row(0) = b0mat.row(0)*fTensor(0,0) + b0mat.row(2)
                        *fTensor(1,0);
                    eDotStar.row(1) = b0mat.row(1)*fTensor(0,1) + b0mat.row(3)
                        *fTensor(1,1);
                    eDotStar.row(2) = b0mat.row(0)*fTensor(0,1) + b0mat.row(2)
                        *fTensor(1,1) + b0mat.row(1)*fTensor(0,0) +b0mat.row(3)
                        *fTensor(1,0);
                    pk2DotStar = Material::C * eDotStar;
                    pDotStar.row(0) = pk2DotStar.row(0)*fTensor(0,0) 
                        + pk2DotStar.row(2)*fTensor(0,1);
                    pDotStar.row(1) = pk2DotStar.row(1)*fTensor(0,0)
                        + pk2DotStar.row(3)*fTensor(0,1);
                    pDotStar.row(2) = pk2DotStar.row(0)*fTensor(1,0)
                        + pk2DotStar.row(2)*fTensor(1,1);
                    pDotStar.row(3) = pk2DotStar.row(1)*fTensor(1,0)
                        + pk2DotStar.row(3)*fTensor(1,1);
                    ktMat.noalias() += scalar4GaussInt * b0mat.transpose() 
                        * pDotStar;
                    // computing geometric tangent stiffness matrix
                    fDotStar.row(0) = b0mat.row(0)*pk2Tensor(0,0) 
                        + b0mat.row(1)*pk2Tensor(1,0);
                    fDotStar.row(1) = b0mat.row(0)*pk2Tensor(0,1)
                        + b0mat.row(1)*pk2Tensor(1,1);
                    fDotStar.row(2) = b0mat.row(2)*pk2Tensor(0,0)
                        + b0mat.row(3)*pk2Tensor(1,0);
                    fDotStar.row(3) = b0mat.row(2)*pk2Tensor(0,1)
                        + b0mat.row(3)*pk2Tensor(1,1);
                    ktGeo.noalias() += scalar4GaussInt * b0mat.transpose()
                        * fDotStar;
                }
            }
            
            elMatContiguousIndx = 0;
            for (const int j : dofEli){
                for (const int i : dofEli){
                    tripletListKtInt.emplace_back(i,j
                        ,ktMat(elMatContiguousIndx)+ktGeo(elMatContiguousIndx));
                    ++elMatContiguousIndx;
                }
            }
        }
        ktInt.setFromTriplets(tripletListKtInt.begin(), tripletListKtInt.end());
    }
}