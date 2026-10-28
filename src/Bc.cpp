#include "Bc.h"

void Bc::setBcInit(ndofVectf& disp, ndofVectf& force) const
{
    std::cout << "\nEmbedding boundary conditions @t=" << Settings::t0 
        << "s ...\n";

    // PAY ATTENTION NOT TO DO INTEGER DIVISION (e.g. 1000/3 = 333 not 333.333)
    // as the .size() of the BC/load vectors is an integral type

    // setting BCx:
    disp(bcs[BCx]).setConstant(0.0);
    // setting BCy:
    disp(bcs[BCy]).setConstant(0.0);
    // setting Fx:
    force(bcs[Fx]).setConstant(0.0);
    // setting Fy:
    force(bcs[Fy]).setConstant(0.0);

    // same can be done for any other bc...
}

void Bc::setBcTimeDependent(ndofVectf& dispNewOrForce
    , const Settings::precision time) const
{
    // use a bc or alternatively create one only for the time dependent bcs
    const auto chosenBc {Fy};

    // creation of custom amplitude
    constexpr const Settings::precision amp {-300.0};
    Settings::precision dispOrFcurrent 
        {amp*time/(Settings::tFinal-Settings::t0)};

    if (bcs[chosenBc].empty()) throw "The time dependent boundary condition selected is empty - Bc.cpp";

    dispNewOrForce(bcs[chosenBc]).setConstant(dispOrFcurrent);
}

void Bc::embedBCintoKt(sparseMat& ktEff) const
{
    std::vector<bool> isDofBC (ktEff.rows(), false);

    for (const auto dof : bcs[BCx]){
        ktEff.coeffRef(dof, dof) = 1.0;
        isDofBC[dof] = true;
    }
    for (const auto dof : bcs[BCy]){
        ktEff.coeffRef(dof, dof) = 1.0;
        isDofBC[dof] = true;
    }

    // clearing rows (and columns to keep symmetry or improve factorization)
    int row {};
    for (int col=0; col<ktEff.rows(); ++col){
        for (eigenIterator it(ktEff,col); it; ++it){
            row = it.row();

            if (row == col) continue;

            if (isDofBC[col] || isDofBC[row]) it.valueRef() = 0.0;

            // following lines commented out, assuming no time-dependent BC disp
            // residual(it.row()) -= it.value()*incrementDisp(dofBC);

        }
    }
}

void Bc::embedBCintoResidual(ndofVectf& residual) const{
    // Assuming for now, no time-dependent BC applied to disp. vector, this 
    // leads to optimization: no need to subtract from residual vector the 
    // product between the ktEff element and the known displacement increment. 
    // As the latter would be 0. 
    // We set the constrained dof of residual to 0.0
    
    // make sure to include all constrained dof (for displacements)
    for (const auto dof : bcs[BCx]){
        residual(dof) = 0.0;
    }
    for (const auto dof : bcs[BCy]){
        residual(dof) = 0.0;
    }
}

void Bc::extractBCAbaqus(){
    constexpr std::array<const char*, Bc::countBCnames> nameBCfiles 
        {"./MESH/SETS/BCx.txt", "./MESH/SETS/BCy.txt"
        , "./MESH/SETS/Fx.txt", "./MESH/SETS/Fy.txt"};
    
    using namespace std::string_literals;
    std::ifstream bcFile;
    std::string line {};
    std::string tempStr;
    tempStr.reserve(10uz);
    int tempVal {};
    std::stringstream tempSStream {};
    
    for (std::size_t nFile{0uz}; nFile<Bc::countBCnames; ++nFile){
        bcFile.open(nameBCfiles[nFile]);
        if (!bcFile)
            throw "Cannot open a BC (boundary condition) file in ./MESH/SETS/";
        while (bcFile >> line){
            for (std::size_t ind{0}; ind<line.length(); ++ind){
                if (line[ind] == ','){
                PushVal:
                    tempSStream << tempStr;
                    tempSStream >> tempVal;
                    // the -1 is to use 0-based indexing
                    if ((nFile == Bc::BCx) || (nFile == Bc::Fx))
                        bcs[nFile].push_back((tempVal-1)*2);
                    else
                        bcs[nFile].push_back((tempVal-1)*2 + 1);
                    tempSStream.clear();
                    tempSStream.str("");
                    tempStr = ""s;
                    continue;
                }
                tempStr.push_back(line[ind]);
                if (ind == (line.length()-1)) goto PushVal;
            }
        }
        bcFile.close();
    }
}

void Bc::extractBCSalome(){
    constexpr std::array<const char*, Bc::countBCnames> nameBCfiles 
        {"./MESH/SETS/BCx.dat", "./MESH/SETS/BCy.dat"
        , "./MESH/SETS/Fx.dat", "./MESH/SETS/Fy.dat"};
    
    using namespace std::string_literals;
    std::ifstream bcFile;
    std::string line {};
    std::string tempStr;
    tempStr.reserve(10uz);
    int tempVal {};
    std::stringstream tempSStream {};
    
    for (std::size_t nFile{0uz}; nFile<Bc::countBCnames; ++nFile){
        bcFile.open(nameBCfiles[nFile]);
        if (!bcFile)
            throw "Cannot open a BC (boundary condition) file in ./MESH/SETS/";
        // Increasing capacity of std::vector according to number of nodes given
        // avoiding expensive dynamic memory re-allocations
        if (!(bcFile >> line)){
            bcFile.close();
            continue;
        }
        tempSStream << line;
        tempSStream >> tempVal;
        bcs[nFile].reserve(static_cast<std::size_t>(tempVal));
        tempSStream.clear();
        tempSStream.str("");

        std::getline(bcFile,line); // skipping the rest of first line
        
        while (std::getline(bcFile,line)){
            for (std::size_t ind{0}; ind<line.length(); ++ind){
                if (line[ind] == ' '){
                    tempSStream << tempStr;
                    tempSStream >> tempVal;
                    // the -1 is to use 0-based indexing
                    if ((nFile == Bc::BCx) || (nFile == Bc::Fx))
                        bcs[nFile].push_back((tempVal-1)*2);
                    else
                        bcs[nFile].push_back((tempVal-1)*2 + 1);
                    tempSStream.clear();
                    tempSStream.str("");
                    tempStr = ""s;
                    break;
                }
                tempStr.push_back(line[ind]);
            }
        }
        bcFile.close();
    }
}