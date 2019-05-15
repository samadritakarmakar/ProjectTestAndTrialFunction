#ifndef FORM_HPP
#define FORM_HPP
#include "TrialFunction.hpp"
#include "TestFunctionGalerkin.hpp"
#include <armadillo>


class Form : public TrialFunction
{
public:
    Form()
    {
        ElementType=0;
        ElementNumber=0;
        GaussPntr=0;
    }
    inline sp_mat u(TrialFunction u)
    {
        return u.Get_u(ElementType,GaussPntr);
    }

    inline sp_mat v(TestFunctionGalerkin v)
    {
        return v.Get_v(ElementType, GaussPntr);
    }

    inline sp_mat grad(TrialFunction u)
    {
        return u.Get_grad_u(ElementType, ElementNumber, GaussPntr);
    }

    inline sp_mat grad(TestFunctionGalerkin v)
    {
        return v.Get_grad_v(ElementType, ElementNumber, GaussPntr);
    }

    inline sp_mat curl(TrialFunction u)
    {
        return u.Get_curl_u(ElementType, ElementNumber, GaussPntr);
    }

    inline sp_mat curl(TestFunctionGalerkin v)
    {
        return v.Get_curl_v(ElementType, ElementNumber, GaussPntr);
    }

    inline sp_mat inner(TestFunctionGalerkin v, TrialFunction u)
    {
        return v.Get_v(ElementType,GaussPntr)*u.Get_u(ElementType,GaussPntr);
    }

    inline sp_mat inner(TestFunctionGalerkin v, sp_mat u)
    {
        return v.Get_v(ElementType,GaussPntr)*u;
    }

    inline sp_mat inner(sp_mat grad_v, sp_mat grad_u)
    {
        return grad_v*grad_u;
    }

    inline sp_mat dot(TestFunctionGalerkin v, sp_mat grad_u)
    {
        return v.Get_v(ElementType,GaussPntr)*grad_u;
    }

    inline sp_mat dot(sp_mat grad_v, sp_mat grad_u)
    {
        return grad_v*grad_u;
    }

    inline sp_mat dot(vec a, sp_mat grad_u)
    {
        return dot_vectrLvl_grad_u(a,grad_u);
    }

    int ElementType;
    int ElementNumber;
    int GaussPntr;
};



#endif // FORM_HPP
