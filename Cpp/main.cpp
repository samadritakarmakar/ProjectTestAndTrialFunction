#include <iostream>
#include "LagrangeShapeFunctionAllElementTypes.hpp"
#include "TrialFunction.hpp"
#include "TestFunctionGalerkin.hpp"
#include "Form.hpp"
#include "CrossProduct.hpp"
#include "TrialFunctionNeumannSurface.hpp"
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
    int vectorLevel=3;
    TrialFunction u(Mesh,vectorLevel);
    //TrialFunction u(FileName, Dimension, vectorLevel);
    TestFunctionGalerkin v(u);
    Form a;

    cout<<"u= \n"<<mat(a.u(u));
    cout<<"Grad(u) =\n"<<mat(a.grad(u));
    cout<<"v:u =\n"<<mat(a.inner(v,u));
    cout<<"grad(v):grad(u) =\n"<<mat (a.inner(a.grad(v),a.grad(u)));
    vec vector={1, 2, 3};
    a.set_u_Internal(u);
    cout<<"v.Vectr.grad u= \n"<<mat(a.inner(v,a.dot(vector,a.grad(u))));
    //cout<<"curl_v.curl_u=\n"<<mat(a.inner(a.curl(v), a.curl(u)));

    mat Matrix1={{5,6,7}};
    mat Matrix2={{8,9,10}};
    cout<<"CrossProduct =\n"<<CrossProduct(Matrix1, Matrix2);

    TrialFunctionNeumannSurface u_surf(u,0);
    return 0;
}
