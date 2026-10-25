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
#include "elMassMatrix.h"
#include "export.h"
#include "Ic.h"
#include "clearResults.h"
#include "internalForceAndKt.h"

int main(){
    //Timer timer {};
    //timer.start();
    Mesh mesh {};
    Bc bc {}; // boundary conditions class
    using ndofVectf = Eigen::Vector<Settings::precision, Eigen::Dynamic>;
    using ndofVecti = Eigen::Vector<int, Eigen::Dynamic>;
    using SizT = std::size_t;
    
    // To handle exceptions in any function/method
    try{
        //-----Importing mesh and boundary conditions
        importMesh(mesh, bc);
        std::cout << "Problem size: " << mesh.nDof << " d.o.f.\n";

        //-----global mass matrix computation
        Eigen::SparseMatrix<Settings::precision> mG (mesh.nDof, mesh.nDof);
        {
            std::vector<Eigen::Triplet<Settings::precision>> tripletListM;
            tripletListM.reserve(mesh.nElements*Settings::nDofKel);
            mMat::MelEigen mEl;
            mMat::NodesCoordEigen nodesCoords;
            mMat::NEigen N;
            Eigen::Matrix<int,Settings::dofNode,Settings::nNodesEl> 
                dofEliMat;
            int mElContiguousIndx {};
            const float numEl1perc {mesh.nElements/20.0f};
            float numEl2Message {numEl1perc};
            for (SizT e=0uz; e<static_cast<SizT>(mesh.nElements); ++e){
                nodesCoords = mesh.nodes(mesh.elements.row(e)
                    ,Eigen::placeholders::all);
                mMat::computeElM(nodesCoords,N,mEl);
                dofEliMat.row(0) = mesh.elements.row(e)*2;
                dofEliMat.row(1) = mesh.elements.row(e)*2 
                    + Eigen::RowVector<int,Settings::nNodesEl>::Ones();
                const auto dofEli = dofEliMat.reshaped(); // view
                // as Eigen Matrices are column-major, stored contiguously 
                // column by column, it's more efficient to read them also 
                // in this way
                mElContiguousIndx = 0;
                for (const int j : dofEli){
                    for (const int i : dofEli){
                        tripletListM.emplace_back(i,j
                            ,mEl(mElContiguousIndx++));
                    }
                }
                if (e > numEl2Message){
                    std::cout << "Sparse global mass & stiffness matrix " 
                        << "computation... " << numEl2Message/numEl1perc
                        *5.0f << "%\n";
                    numEl2Message += numEl1perc;
                }
            }
            mG.setFromTriplets(tripletListM.begin(), tripletListM.end());
            std::cout << "\nGlobal Mass [M] matrix computed"
                << " successfully!\n\n";
        } // Triplet list and all other helper matrices deallocated

        //-----Solving system of eqs. (full matrix method, no partitioning) 
        //     + Newton-Raphson + time integration (Newmark-beta method)

        ndofVectf force, dispOld, dispNew, velOld, velNew, accOld, accNew;
        
        Ic::setIC(dispOld, velOld, accOld, force, mG);
        bc.setBcInit(dispOld, force);
        dispNew = dispOld;
        velNew  = velOld;
        accNew  = accOld;

        Settings::precision time {Settings::t0 + Settings::dt};
        Settings::precision tSaveVTK {Settings::t0 + Settings::dtSaveVTK};
        int countVTK {0};
        Settings::precision normR {100.0};
        int newtonRaphsonCounter {0};

        clearPriorResults();

        constexpr const Settings::precision dt2 {Settings::dt*Settings::dt};
        constexpr const Settings::precision coeffM {1.0/(Settings::beta*dt2)};
        constexpr const Settings::precision coeffAcc {0.5-Settings::beta};

        ndofVectf forceInt;
        ndofVectf residual;
        Eigen::SparseMatrix<Settings::precision> ktInt (mesh.nDof, mesh.nDof);
        Eigen::SparseMatrix<Settings::precision> ktEff (mesh.nDof, mesh.nDof);

        while (time < (Settings::tFinal + Settings::dt/2.0))
        {
            if (!Settings::importICFromStatic && (newtonRaphsonCounter == 0)){
                bc.setBcTimeDependent(force,time);
            }

            accNew = coeffM*(dispNew - dispOld - velOld*Settings::dt 
                - coeffAcc*dt2*accOld);
            velNew = velOld + (1.0-Settings::gamma)*Settings::dt*accOld
                + Settings::gamma*Settings::dt*accNew;
            
            FintKtInt::compute(mesh, dispNew, forceInt, ktInt);
            
            residual.noalias() = mG*accNew;
            residual.noalias() += (forceInt - force);
            
            normR = residual.cwiseAbs().maxCoeff()/force.cwiseAbs().maxCoeff();
            
            if (normR < Settings::residualTargetTol){
                if (time > (tSaveVTK-Settings::dt/2.0)){
                    Export::exportResultsBinary(mesh, dispNew, countVTK, time);
                    tSaveVTK += Settings::dtSaveVTK;
                    ++countVTK;
                    std::cout << "Time instant " << time << " / " 
                        << Settings::tFinal << " s saved!\n";
                }

                std::swap(dispOld, dispNew);
                std::swap(velOld, velNew);
                std::swap(accOld, accNew);

                time += Settings::dt;
                newtonRaphsonCounter = 0;
                continue;
            }
            
            ktEff = coeffM * mG + ktInt;

            
            // modify KtEff as well as residual to avoid using partitioning
            // solve for dispNew
            if (++newtonRaphsonCounter > Settings::maxNewtonIters){
                throw "Maximum number of Newton-Raphson iterations reached!";
            }
        }

        std::cout << "\nResults saved in " << Settings::resultsDirName << '\n';
        std::cout << "To display time in paraview: filter -> python annotation"
            << " -> 'Time: %.3e' % (input.FieldData['TIME'][0])\n";
        //timer.stop();
        //std::cout<<"Elapsed: "<<timer.getDurationMicro()<<" micro sec.\n";
    }
    catch(const char* exceptionMsg){
        std::cerr << "Error: " << exceptionMsg << '\n';
    }
    catch(...){
        std::cerr << "Unrecognized error occured!\n";
    }

    return 0;
}