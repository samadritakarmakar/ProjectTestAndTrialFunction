#ifndef LOCALINTEGRATION_HPP
#define LOCALINTEGRATION_HPP
#include "Form.hpp"
#include "TrialFunction.hpp"
#include "TrialFunctionNeumannSurface.hpp"
#include "TrialFunctionNeumannLine.hpp"
#include "TestFunctionGalerkin.hpp"
template <class GenericTrialFunction>
class LocalIntegrator
{
public:
    LocalIntegrator(Form<GenericTrialFunction>& a, GenericTrialFunction& u, TestFunctionGalerkin<GenericTrialFunction>& v):
        a(a), u(u),v(v)
    {

    }

    virtual ~LocalIntegrator()
    {
        //----Do-Nothing
    }

    virtual sp_mat weak_form(Form<GenericTrialFunction>& a, GenericTrialFunction& u, TestFunctionGalerkin<GenericTrialFunction>& v);

    void local_intergrator()
    {
        a.set_u_Internal(u);
        a.set_v_Internal(v);
        //cout<<"Element Type ="<<a.ElementType<<"Gauss pt is at "<<a.GaussPntr<<"\n";
        a=weak_form(a,u,v);
        int NoOfGaussPnts=u.GetNumberOfGaussPoints(a.ElementType);
        for (int i=0; i<NoOfGaussPnts-1; i++)
        {
            a.NextGaussPntr();
            //cout<<"Element Type ="<<a.ElementType<<"Gauss pt is at "<<a.GaussPntr<<"i= "<<i<<"\n";
            a=a.ResultingMat+weak_form(a,u,v);
        }
        a.GaussPntr=0;
    }


    Form<GenericTrialFunction>& a;
    GenericTrialFunction& u;
    TestFunctionGalerkin<GenericTrialFunction>& v;

};

template <class GenericTrialFunction>
sp_mat LocalIntegrator<GenericTrialFunction>::weak_form(Form<GenericTrialFunction>& a,
                                                        GenericTrialFunction& u, TestFunctionGalerkin<GenericTrialFunction>& v)
{
    return a.inner(a.grad(v),a.grad(u))*a.dX(u);;
}


#endif // LOCALINTEGRATION_HPP
