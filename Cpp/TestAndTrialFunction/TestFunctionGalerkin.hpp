#ifndef TESTFUNCTIONGALERKIN_HPP
#define TESTFUNCTIONGALERKIN_HPP
#include "TrialFunction.hpp"
class TestFunctionGalerkin
{
public:
    TestFunctionGalerkin(TrialFunction u1)
    {
        u=u1;
        //cout<<"v =\n"<<mat(Get_v(0,0))<<"\n";
        //cout<<"Grad_v\n"<<mat(Get_grad_v(0,0,0))<<"\n";
    }
    /// Return Shape Function of v in Matrix form
    inline sp_mat Get_v(int ElementType, int GaussPntr)
    {
        return u.Get_u(ElementType, GaussPntr).t();
    }
    ///
    inline sp_mat Get_grad_v(int ElementType, int ElementNumber, int GaussPntr)
    {
        return u.Get_grad_u(ElementType, ElementNumber, GaussPntr).t();
    }

    void Get_F(int ElementType, int ElementNumber, int GaussPntr, mat& F)
    {
        u.Get_F(ElementType, ElementNumber, GaussPntr, F);
    }

private:
    TrialFunction u;
};

#endif // TESTFUNCTION_HPP
