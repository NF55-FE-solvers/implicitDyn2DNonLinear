#include <iostream>
#include <fstream>
#include <algorithm>
#include "Eigen/Eigen"
#include "Eigen/PardisoSupport"
#include "Mesh.h"
#include "import.h"
#include "settings.h"
#include "TimerHighRes.h"
#include "material.h"
#include "Bc.h"
#include "rearrangeDof.h"
#include "elMatrices.h"
#include "export.h"
#include "Ic.h"
#include "clearResults.h"

int main(){
    //Timer timer {};
    //timer.start();
    Mesh mesh {};
    Bc bc {}; // boundary conditions class
    using ndofVectf = Eigen::Vector<Settings::precision, Eigen::Dynamic>;
    using ndofVecti = Eigen::Vector<int, Eigen::Dynamic>;
    using SizT = std::size_t;
    ndofVectf dispOld;
    ndofVectf force;
    
    // To handle exceptions in any function/method
    try{
        //-----Importing mesh and boundary conditions
        importMesh(mesh, bc);
        std::cout << "Problem size: " << mesh.nDof << " d.o.f.\n";

        //-----(global mass matrix + global stiffness matrix) computation
        Eigen::SparseMatrix<Settings::precision> kG (mesh.nDof, mesh.nDof);
        Eigen::SparseMatrix<Settings::precision> mG (mesh.nDof, mesh.nDof);
        {
            std::vector<Eigen::Triplet<Settings::precision>> tripletListK;
            std::vector<Eigen::Triplet<Settings::precision>> tripletListM;
            tripletListK.reserve(mesh.nElements*Settings::nDofKel);
            tripletListM.reserve(mesh.nElements*Settings::nDofKel);
            kmMat::KMelEigen kEl;
            kmMat::KMelEigen mEl;
            kmMat::NodesCoordEigen nodesCoords;
            kmMat::NEigen N;
            kmMat::dNEigen dN;
            kmMat::invJEigen invJ;
            kmMat::BEigen bEl;
            Eigen::Matrix<int,Settings::dofNode,Settings::nNodesEl> 
                dofEliMat;
            int kmElContiguousIndx {};
            const float numEl1perc {mesh.nElements/20.0f};
            float numEl2Message {numEl1perc};
            for (SizT e=0uz; e<static_cast<SizT>(mesh.nElements); ++e){
                nodesCoords = mesh.nodes(mesh.elements.row(e)
                    ,Eigen::placeholders::all);
                kmMat::computeElKM(nodesCoords,N,dN,invJ,bEl,kEl,mEl);
                dofEliMat.row(0) = mesh.elements.row(e)*2;
                dofEliMat.row(1) = mesh.elements.row(e)*2 
                    + Eigen::RowVector<int,Settings::nNodesEl>::Ones();
                const auto dofEli = dofEliMat.reshaped(); // view
                // as Eigen Matrices are column-major, stored contiguously 
                // column by column, it's more efficient to read them also 
                // in this way
                kmElContiguousIndx = 0;
                for (const int j : dofEli){
                    for (const int i : dofEli){
                        tripletListM.emplace_back(i,j
                            ,mEl(kmElContiguousIndx));
                        tripletListK.emplace_back(i,j
                            ,kEl(kmElContiguousIndx));
                        ++kmElContiguousIndx;
                    }
                }
                if (e > numEl2Message){
                    std::cout << "Sparse global mass & stiffness matrix " 
                        << "computation... " << numEl2Message/numEl1perc
                        *5.0f << "%\n";
                    numEl2Message += numEl1perc;
                }
            }
            kG.setFromTriplets(tripletListK.begin(), tripletListK.end());
            mG.setFromTriplets(tripletListM.begin(), tripletListM.end());
            std::cout << "\nGlobal Mass [M] and Stiffness [K] matrices computed"
                << " successfully!\n\n";
        } // Triplet list and all other helper matrices deallocated

        //-----Solving system of eqs. (splitting method) 
        //     & time integration (Newmark-beta method)

        dispOld = ndofVectf::Constant(mesh.nDof,Settings::NaN);
        force = ndofVectf::Zero(mesh.nDof);

        bc.setBcInit(dispOld, force);

        // Re-arranging vectors (splitting method) + saving mapping vector
        ndofVecti mapDof;
        mapDof = ndofVecti::LinSpaced(mesh.nDof, 0, mesh.nDof-1);
        rearrange(dispOld, force, mapDof);
        Eigen::PermutationMatrix<Eigen::Dynamic, Eigen::Dynamic, int> 
            mapDofMat {mapDof};
        kG = kG.twistedBy(mapDofMat.transpose()); // Re-arranging
        mG = mG.twistedBy(mapDofMat.transpose()); // Re-arranging

        int lastNanIndx {mesh.nDof-1};
        for (; lastNanIndx>0; --lastNanIndx){
            if (std::isnan(dispOld(lastNanIndx))) break;
        }
        int nUnknownDisp {lastNanIndx+1};
        int nKnownDisp {mesh.nDof - nUnknownDisp};

        ndofVectf dispNew, velOld, velNew, accOld, accNew, dispReord;
        
        Ic::setIC(dispOld, velOld, accOld, force, mG, kG, nUnknownDisp, mapDof);
        dispNew = dispOld;
        velNew  = velOld;
        accNew  = accOld;

        Settings::precision time {Settings::t0 + Settings::dt};
        Settings::precision tSaveVTK {Settings::t0 + Settings::dtSaveVTK};
        int countVTK {0};

        clearPriorResults();

        constexpr const Settings::precision dt2 {Settings::dt*Settings::dt};
        constexpr const Settings::precision coeffM {1.0/(Settings::beta*dt2)};
        constexpr const Settings::precision coeffAcc {0.5-Settings::beta};
        constexpr const Settings::precision coeffVel0 {Settings::gamma
            /(Settings::beta*Settings::dt)};
        constexpr const Settings::precision coeffVel1 {1.0-Settings::gamma
            /Settings::beta};
        constexpr const Settings::precision coeffVel2 {1.0-0.5*Settings::gamma
            /Settings::beta};
        
        Eigen::SparseMatrix<Settings::precision> effStiffMat 
            (nUnknownDisp,nUnknownDisp);
        effStiffMat = coeffM*mG.block(0,0,nUnknownDisp,nUnknownDisp) 
            + kG.block(0,0,nUnknownDisp,nUnknownDisp);
        Settings::solverType solver;
        solver.compute(effStiffMat);
        if(solver.info()!=Eigen::Success) {
            throw "Linear solver failed decomposition.";
        }

        Eigen::Vector<Settings::precision,Eigen::Dynamic> rhsEquation;

        while (time < (Settings::tFinal + Settings::dt/2.0)){
            if (!Settings::importICFromStatic){
                bc.setBcTimeDependent(force,mapDof,time,nUnknownDisp);
            }

            rhsEquation.noalias() = force.head(nUnknownDisp);
            rhsEquation.noalias() += coeffM*mG.topLeftCorner(nUnknownDisp
                ,nUnknownDisp)*(dispOld.head(nUnknownDisp) 
                + Settings::dt*velOld.head(nUnknownDisp) 
                + coeffAcc*dt2*accOld.head(nUnknownDisp));
            rhsEquation.noalias() += coeffM*mG.topRightCorner(nUnknownDisp
                ,nKnownDisp)*(-dispNew.tail(nKnownDisp) 
                + dispOld.tail(nKnownDisp) 
                + Settings::dt*velOld.tail(nKnownDisp) 
                + coeffAcc*dt2*accOld.tail(nKnownDisp));
            rhsEquation.noalias() -= kG.topRightCorner(nUnknownDisp,nKnownDisp)
                *dispNew.tail(nKnownDisp);
            
            dispNew.head(nUnknownDisp) = solver.solve(rhsEquation);
            if(solver.info()!=Eigen::Success){
                throw "Linear solver failed solving(1).";
            }

            velNew.head(nUnknownDisp).noalias() = coeffVel0
                *(dispNew.head(nUnknownDisp)-dispOld.head(nUnknownDisp))
                + coeffVel1*velOld.head(nUnknownDisp) + coeffVel2*Settings::dt
                *accOld.head(nUnknownDisp);
            
            accNew.head(nUnknownDisp).noalias() = coeffM
                *(dispNew.head(nUnknownDisp) - dispOld.head(nUnknownDisp) 
                - velOld.head(nUnknownDisp)*Settings::dt - coeffAcc
                *dt2*accOld.head(nUnknownDisp));

            if (time > (tSaveVTK-Settings::dt/2.0)){
                dispReord = mapDofMat * dispNew;
                Export::exportResultsBinary(mesh, dispReord, countVTK, time);
                tSaveVTK += Settings::dtSaveVTK;
                ++countVTK;
                std::cout << "Time instant " << time << " / " 
                    << Settings::tFinal << " s saved!\n";
            }

            std::swap(dispOld, dispNew);
            std::swap(velOld, velNew);
            std::swap(accOld, accNew);

            time += Settings::dt;
        }

        std::cout << "\nResults saved in " << Settings::resultsDirName << '\n';
        std::cout << "To display time in paraview: filter -> python annotation"
            << " -> 'Time: %.3e' % (input.FieldData['TIME'][0])\n";
        //timer.stop();
        //std::cout << "Elapsed: " << std::setprecision(10) 
        //    << timer.getDurationMicro() << " micro sec.\n";
    }
    catch(const char* exceptionMsg){
        std::cerr << "Error: " << exceptionMsg << '\n';
    }
    catch(...){
        std::cerr << "Unrecognized error occured!\n";
    }

    return 0;
}
