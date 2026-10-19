#include "Bc.h"

void Bc::setBcInit(Bc::ndofVectf& disp, Bc::ndofVectf& force) const
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

void Bc::setBcTimeDependent(Bc::ndofVectf& dispNewOrForce
    , const Bc::ndofVecti& mapDof, const Settings::precision time
    , const int mapDofBlockSize)
{
    // creation of custom amplitude
    constexpr const Settings::precision FAmp {-20.0};
    Settings::precision Fcurrent {FAmp*time/(Settings::tFinal-Settings::t0)};
    
    // use a bc or alternatively create one only for the time dependent bcs
    const auto chosenBc {Fy};

    if (bcs[chosenBc].empty()) throw "The time dependent boundary condition selected is empty - Bc.cpp";

    // making sure the elements of bcs[chosenBc] are sorted (see comments below)
    using sizT = std::size_t;
    if (bcs[chosenBc].size() > 1){
        bool isSorted {true};
        for (sizT i=0uz; i<(bcs[chosenBc].size()-1uz); ++i){
            if (bcs[chosenBc][i] > bcs[chosenBc][i+1uz]){
                isSorted = false;
                break;
            }
        }
        if (!isSorted){
            std::sort(bcs[chosenBc].begin(), bcs[chosenBc].end());
        }
    }
    // here using a more "complex" algorithm, but should be *faster than: 
    // force(mapDof) = force; 
    // force(bcs[chosenBc]) = Fcurrent; 
    // force = force(mapDof);
    // *if the number of dof interested by the time dep. bc are a small part of
    // total dof.
    int indx2findBc {0};
    // exploiting the fact that mapDof is sorted in blocks
    int startIndxMap {bcs[chosenBc][indx2findBc] 
        - mapDof.head(mapDofBlockSize)(0)};
    for (int indx {startIndxMap}; indx>=0; --indx){
        if (mapDof(indx) == bcs[chosenBc][indx2findBc]){
            dispNewOrForce(indx) = Fcurrent;
            startIndxMap = indx + 1;
            ++indx2findBc;
            break;
        }else if (indx==0) throw "ID/index of time dependent bc couldn't be found in the specified block of mapDof - Bc.cpp";
    }
    // now continuing with the rest of indices
    for (int indx {startIndxMap}; indx<mapDofBlockSize; ++indx){
        if (indx2findBc >= static_cast<int>(bcs[chosenBc].size())) break;
        if (mapDof(indx) == bcs[chosenBc][indx2findBc]){
            dispNewOrForce(indx) = Fcurrent;
            ++indx2findBc;
        }
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