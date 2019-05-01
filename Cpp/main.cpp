#include <iostream>
#include "LagrangeShapeFunctionAllElementTypes.hpp"
#include "TrialFunction.hpp"
using namespace arma;
int main(int argc, char *argv[])
{
    if(argc==1 || argc>3)
    {
        std::cout<<"Usage: ./ProjectShapeFunction <.msh Filename> <Dimension>\n";
        return 0;
    }
    std::string FileName(argv[1]);
    int Dimension=*argv[2]-'0';
    //cout<<"File Name= "<<FileName<<"\n";
    //cout<<"Dimension= "<<Dimension<<"\n";
    libGmshReader::MeshReader Mesh(FileName, Dimension);
    /*LagrangeShapeFunctionAllElementTypes trial(Mesh);
    mat GaussPointx={0};
    mat GaussPointy={0};
    mat GaussPointz={0};
    //cout<<"Element Name= "<<Mesh.GmshElementNameOnly[0]<<"\n";
    std::vector<mat> N=trial.GetShapeFunction(GaussPointx, GaussPointy);
    for (int i=0; i<N.size(); i++)
    {
        cout<<N[i];
    }*/
    TrialFunction u(Mesh, 2);
    return 0;
}
