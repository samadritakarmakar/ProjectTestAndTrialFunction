# ProjectTestFunction
## The goal of the project is to build an open source general finite element solver.  
The basic inputs shall be a .msh file generated by Gmsh (http://gmsh.info/), a free 3D finite element mesh generator.  
The weak form of can be declared in the main.cpp itself. All boundary conditions, Neumann as well as Dirichlet can be delared within the main file.  
The Nodes for Nuemann and Dirichlet boundary conditions are detected on the nodes declared as 'Physical Entities' in gmsh file.  
All the mesh data including physical entities must be declared within one file. For this, after meshing the domain, use 'Export..' option of Gmsh and save as '.msh' file.  
The Project is using Armadillo (http://arma.sourceforge.net) and Gmsh (gmsh.info/) APIs.   
  
The Project is a work in progress.  
  
To compile the code, you need to download the latest Gmsh source code from gmsh.info/ and compile  by typing the following in your Gmsh build directory. This will build a static library of Gmsh.  

    cmake -DDEFAULT=0 -DENABLE_BUILD_LIB=1 -DENABLE_POST=1 -DENABLE_PARSER=1 <path gmsh source directory>
    make lib
    
Copy the libgmsh.a file generated to the folder ProjectTestFunction/Cpp/libGmshReader/GmshApi/  
Run the following command in your ProjectTestFunction Build directory.

    cmake <path ProjectTestFunction/Cpp/ source directory>
    make
    
A new file ProjectTestFunction will be created

Run by typing

    Usage: ./ProjectTestFunction <.msh Filename> <Dimension>
    
where <.msh Filename> is the '.msh' filename and <Dimension> is the Dimension of the domain, 



