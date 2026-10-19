# Structural-FE-codes
(For info about this project go to https://nf55-55.github.io/projects/Efficient%20structural%20FE%20code/)

## Requirements:
- Eigen library: create 'EigenLib' folder in the root (already ignored with .gitignored) and paste the Eigen library inside of it (e.g. EigenLib/eigen-5.0.0/)
- (Not essential but suggested) Intel oneMKL to use a state of the art linear
solver. If not, you can still use the (slower) Eigen built-in linear solvers.

## Current workflow:
- create the mesh in either: Abaqus or SALOME.
- create the sets/groups of nodes on which the boundary conditions (BC) and loads will be applied (fundamental for large/complex meshes to do this
	step in the mesher itself, as the alternative is to smartly select the nodes by using their x,y coordinates info).
- for Abaqus write out the .inp file (which contains both the mesh and sets of nodes). 
	While for SALOME export the mesh as .dat file as well as the mesh groups containing the BC and load nodes as .dat files.
- Add the mesh .dat or .inp file to the ./MESH folder, and rename the file with the name of the software used, e.g.: "Abaqus.inp", "SALOME.dat". The only valid names that should be included in the name of the file are: "Abaqus" or "SALOME". This allows the program to recognize automatically the type of mesh file.
- For Abaqus the node sets should be created manually, by copy-pasting the list of node IDs of the .inp file in new .txt files making sure all IDs
	are on the same line, and added to "./MESH/SETS" folder. Sometimes Abaqus recognizes patterns in the node IDs
	so it's possible that only the first, last and the increment between consecutive ones is given, in that case use the python script in the ./MESH/SETS folder to generate the list of node IDs.
	For SALOME just add the .dat files containing the node groups in the ./MESH/SETS folder.
- SET material properties, solver settings, IC and BC + Load amplitudes by modifying the following files: material.h, settings.h, Bc.cpp and Bc.h.
- Build (CMake) the program.
- execute the program (in ./build/ folder) from the root folder.
- the set of .vtk files containing the results will be saved in "./VTK RESULTS" folder.

[To set the IC from a static solution activate the corresponding option in settings.h, same should be done in the static solver before simulation. A binary file will be generated in ./static2dynamic after static simulation, copy this file and paste it in the same folder of the dynamic solver.]
