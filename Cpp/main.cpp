#include <iostream>
#include "LagrangeShapeFunctionAllElementTypes.hpp"
#include "TrialFunction.hpp"
#include "TestFunctionGalerkin.hpp"
#include "Form.hpp"
#include "TrialFunctionNeumannSurface.hpp"
#include "TrialFunctionNeumannLine.hpp"
#include "LocalIntegration.hpp"
#include "SystemAssembly.hpp"
using namespace arma;
class new_LocalIntegrator: public LocalIntegrator<TrialFunction>
{
public:
    new_LocalIntegrator(Form<TrialFunction>& a, TrialFunction& u, TestFunctionGalerkin<TrialFunction>& v):
        LocalIntegrator (a,u,v)
    {
    }
    sp_mat weak_form(Form<TrialFunction>& a, TrialFunction& u, TestFunctionGalerkin<TrialFunction>& v)
    {
        //return a.inner(v,u)*a.dX(u);
        //vec vector={1, 2, 3};
        //return a=a.inner(v,a.dot(vector,a.grad(u)))*a.dX(u);
        return a.inner(a.grad(v),a.grad(u))*a.dX(u);
        //return a.dot(a.curl(v),a.curl(u));
    }
};

class new_Neu_Surf_LclIntgrtr: public LocalIntegrator<TrialFunctionNeumannLine>
{
public:
    new_Neu_Surf_LclIntgrtr(Form<TrialFunctionNeumannLine>& a, TrialFunctionNeumannLine& u,
                            TestFunctionGalerkin<TrialFunctionNeumannLine>& v):
        LocalIntegrator (a,u,v)
    {
    }
    sp_mat weak_form(Form<TrialFunctionNeumannLine>& a, TrialFunctionNeumannLine& u,
                     TestFunctionGalerkin<TrialFunctionNeumannLine>& v)
    {
        //return a.inner(v,u)*a.dX(u);

        //vec vector={1, 2, 3};
        //return a=a.inner(v,a.dot(vector,a.grad(u)))*a.dX(u);
        return a.dot(v,sp_mat(ones(a.v(v).n_cols ,1)))*a.dL(u);
    }
};

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

    int vectorLevel=3;
    TrialFunction u(Mesh,vectorLevel);
    TestFunctionGalerkin<TrialFunction> v(u);
    Form<TrialFunction> a;
    //std::shared_ptr<LocalIntegrator<TrialFunction>> intgrt (new LocalIntegrator<TrialFunction>(a,u,v));
    std::shared_ptr<new_LocalIntegrator> lcl_intgrt(new new_LocalIntegrator(a,u,v));
    SystemAssembler<new_LocalIntegrator> systmAssmbly(a,u,v);
    systmAssmbly.SetLocalIntegrator(lcl_intgrt);
    sp_mat A;
    systmAssmbly.RunSystemAssembly(A);
    //cout<<(A);
    //intgrt=lcl_intgrt;
    //intgrt->local_intergrator();
    //cout<<"grad(v):grad(u)*dx =\n"<<mat(a.ResultingMat);

    //cout<<"curl_v.curl_u=\n"<<mat(a.inner(a.curl(v), a.curl(u)));

    /*TrialFunctionNeumannLine u_surf(u,0);
    Form<TrialFunctionNeumannLine> a2;
    TestFunctionGalerkin<TrialFunctionNeumannLine> v2(u_surf);
    std::shared_ptr<LocalIntegrator<TrialFunctionNeumannLine>> intgrt2
            (new LocalIntegrator<TrialFunctionNeumannLine>(a2,u_surf,v2));
    std::shared_ptr<new_Neu_Surf_LclIntgrtr> lcl_intgrt2(new new_Neu_Surf_LclIntgrtr(a2,u_surf,v2));
    intgrt2=lcl_intgrt2;
    intgrt2->local_intergrator();
    cout<<"grad(v):grad(u)*dx Surface =\n"<<mat(a2.ResultingMat);
   // u.NoOfElementTypes;
   // cout<<"Element Number = "<<a.ElementNumber
   //    <<"\n where the Element Nodes are\n" <<u.Msh->ElementNodes[0].row(a.ElementNumber);
    //Form<TrialFunctionNeumannLine> a3;
    /*TrialFunctionNeumannLine u_line(u,0);
    TestFunctionGalerkin<TrialFunctionNeumannLine> v_line(u_line);
    cout<<"u_line=  \n"<<mat(a2.u(u_line));
    cout<<"Grad(u_line) =\n"<<mat(a2.grad(u_line));
    a2=a2.inner(v_line,u_line)*a2.dL(u_line);
    cout<<"v_line:u_line*dl =\n"<<mat(a2.ResultingMat);
    a2=a2.inner(a2.grad(v_line),a2.grad(u_line))*a2.dL(u_line);
    cout<<"grad(v_line):grad(u_line)*dl =\n"<<mat (a2.ResultingMat);
    */
    return 0;
}
