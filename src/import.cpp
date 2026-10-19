#include "import.h"


// main function for mesh import
void importMesh(Mesh& mesh, Bc& bc){
    std::cout << "\nNOTICE: for now, Abaqus (.inp) and SALOME (.dat) input/mesh"
        << " files are supported, \nrename these with the name of the program "
        << "(Abaqus.inp | SALOME.dat)\n\n";
    
    // checking current working directory against expected one
    Settings::fPath cwd {std::filesystem::current_path().filename()};
    if (cwd != Settings::expectedRootProj){
        throw "Current Working Directory different from expected one (see settings), change it accordingly.";
    }
    
    constexpr std::array<const char*,2> allowedMeshFiles {
        {"./MESH/Abaqus.inp", "./MESH/SALOME.dat"}
    };

    enum MeshFileNames : std::size_t{
        Abaqus,
        SALOME,
        maxNames,
    };

    static_assert(MeshFileNames::maxNames == allowedMeshFiles.size());

    std::ifstream meshFile;

    for (std::size_t ind {0}; ind<allowedMeshFiles.size(); ++ind){
        meshFile.open(allowedMeshFiles[ind]);
        if (meshFile){
            switch (ind)
            {
            case Abaqus:
                extractMeshFromAbaqusInp(meshFile, mesh);
                bc.extractBCAbaqus();
                std::cout << "\nAbaqus mesh imported successfully!\n\n";
                meshFile.close();
                return;
            case SALOME:
                extractMeshFromSalomeDat(meshFile, mesh);
                bc.extractBCSalome();
                std::cout << "\nSalome mesh imported successfully!\n\n";
                meshFile.close();
                return;
            }
        }
    }

    meshFile.close();

    throw "In /MESH/ there's no file with a compatible name - importMesh()";
    
}


// Functions used to extract mesh for each different mesh file:

void extractMeshFromAbaqusInp(std::ifstream& file, Mesh& mesh){
    using Mat2d = std::vector<std::vector<Settings::precision>>;
    using Mat2dI = std::vector<std::vector<int>>;
    using namespace std::string_literals;
    Mat2d mat4Nodes;
    mat4Nodes.reserve(100uz); // reserving an initial temptative capacity
    Mat2dI mat4Elements;
    mat4Elements.reserve(100uz); // reserving an initial temptative capacity
    std::string word {};
    std::string tempStr;
    tempStr.reserve(10uz); // reserving some capacity to avoid too frequent
    std::stringstream tempSStream {};
    Settings::precision tempValx {};
    Settings::precision tempValy {}; 
    std::array<int,4> tempValEl {};
    // dynamic memory reallocation
    // Key positions: header of node list -> header of element list -> end of 
    // element list
    bool nodeListHeaderFound {false};
    bool elementListHeaderFound {false};

    // Iterating through the meshFile
    while (file >> word){
        if (word == "*Node"s){
            nodeListHeaderFound = true;
            file >> word; // next word
        }
        if (word == "*Element,"s){
            elementListHeaderFound = true;
            std::getline(file,word); // getting the entire line so that the next
            // word is surely in the next line
            file >> word; // next word
        }
        if ((word[0uz] == '*') && (elementListHeaderFound)){
            break;
        }

        
        // Node list
        if (nodeListHeaderFound && !elementListHeaderFound){
            // Iterating a single line (assuming 3 entries, 2D case)
            file >> word; // skipping the node number
            tempStr = ""s;
            // removing ','
            for (std::size_t i{0uz}; i<word.length(); ++i){
                if (word[i] == ',') break;
                tempStr.push_back(word[i]);
            }
            // Entering the coordinates of current node
            tempSStream << tempStr;
            tempSStream >> tempValx;
            tempSStream.clear();
            tempSStream.str("");    
            file >> word;
            tempSStream << word;
            tempSStream >> tempValy; //possible runtime error if there's a ','
            tempSStream.clear();
            tempSStream.str("");
            mat4Nodes.push_back({{tempValx, tempValy}});
        }
        
        // Element list
        if (elementListHeaderFound){
            // Iterating a single line (assuming 4 entries, quad element)
            file >> word;
            for (std::size_t nodei{0uz}; nodei<3uz; ++nodei){
                tempStr = ""s;
                for (std::size_t i{0uz}; i<word.length(); ++i){
                    if (word[i] == ',') break;
                    tempStr.push_back(word[i]);
                }
                // Entering the node numbers of current element
                tempSStream << tempStr;
                tempSStream >> tempValEl[nodei];
                tempSStream.clear();
                tempSStream.str("");
                file >> word;
            }
            tempSStream << word;
            //possible runtime error if there's a ','
            tempSStream >>tempValEl[3uz];
            tempSStream.clear();
            tempSStream.str("");
            mat4Elements.push_back({tempValEl[0uz],tempValEl[1uz],tempValEl[2uz]
                , tempValEl[3uz]});
        }
    }
    // populating the mesh struct with the acquired data
    mesh.nodes.resize(mat4Nodes.size(), mat4Nodes[0uz].size());
    for (std::size_t indRow{0uz}; indRow<mat4Nodes.size(); ++indRow){
        mesh.nodes.row(indRow) = Eigen::Matrix<Settings::precision,1,2
            ,Eigen::RowMajor>::Map(&mat4Nodes[indRow][0uz]);
    }
    mesh.elements.resize(mat4Elements.size(), mat4Elements[0].size());
    for (std::size_t indRow{0uz}; indRow<mat4Elements.size(); ++indRow){
        mesh.elements.row(indRow) = Eigen::Matrix<int,1,4,Eigen::RowMajor>::Map
            (&mat4Elements[indRow][0uz]);
    }
    mesh.elements.array() -= 1; // so to refer to nodes with 0-based indexing
    mesh.nDof = mesh.nodes.rows() * Settings::dofNode;
    mesh.nElements = mesh.elements.rows();
}


void extractMeshFromSalomeDat(std::ifstream& file, Mesh& mesh){
    // Exploiting the fact that SALOME .dat mesh file gives the node count on 
    // top of the file:
    int nodeCount {};
    int maxElementCount {};
    std::string word;
    word.resize(25uz);
    file >> word;
    std::stringstream tempSStream {};
    tempSStream << word;
    tempSStream >> nodeCount;
    tempSStream.clear();
    tempSStream.str("");
    file >> word; // this is NOT the element count, but used to define the max
    // size of the element matrix
    tempSStream << word;
    tempSStream >> maxElementCount;
    tempSStream.clear();
    tempSStream.str("");
    // Resizing Eigen dynamic matrices accordingly (so to avoid many dynamic 
    // reallocations)
    mesh.nodes.resize(nodeCount, 2);
    Eigen::Matrix<int,Eigen::Dynamic,4,Eigen::RowMajor> elements;
    elements.resize(maxElementCount,4);

    // Extracting node list
    for (int ind{0}; ind<nodeCount; ++ind){
        file >> word; // skipping the node number
        for (int i{0}; i<2; ++i){
            file >> word;
            tempSStream << word;
            tempSStream >> mesh.nodes(ind,i);
            tempSStream.clear();
            tempSStream.str("");
        }
        file >> word; // skipping the third dimension
    }
    
    using namespace std::string_literals;

    // Extracting element list
    // Iterating forward until the 2nd value in line is 204 (quad element)
    bool isElementID204 {false};
    while (!isElementID204){
        std::getline(file, word);
        for (std::size_t i{0}; i<word.length(); ++i){
            if (word[i] == ' '){ // assuming we're not at then end
                if ((word[i+1uz]=='2') && (word[i+2uz]=='0') 
                    && (word[i+3uz]=='4'))
                {
                    std::string stringBuffer {};
                    int col {0};
                    isElementID204 = true;
                    for (std::size_t j{i+5uz}; j<word.length(); ++j){
                        if (word[j] == ' '){
                            tempSStream << stringBuffer;
                            tempSStream >> elements(0,col++);
                            tempSStream.clear();
                            tempSStream.str("");
                            stringBuffer = ""s;
                            continue;
                        }
                        stringBuffer.push_back(word[j]);
                    }
                }
                break;
            }
        }
    }
    // continue element list extraction in sequence
    int row {1};
    while (file >> word){
        file >> word; // skipping elementID
        for (int col{0}; col<4; ++col){
            file >> word;
            tempSStream << word;
            tempSStream >> elements(row,col);
            tempSStream.clear();
            tempSStream.str("");
        }
        ++row;
    }
    mesh.elements.resize(row,4);
    mesh.elements = elements.block(0,0,row,4);
    mesh.elements.array() -= 1; // so to refer to nodes with 0-based indexing
    mesh.nDof = mesh.nodes.rows() * Settings::dofNode;
    mesh.nElements = mesh.elements.rows();
}
