#ifndef FORM_HPP
#define FORM_HPP
#include "TrialFunction.hpp"
#include "TestFunctionGalerkin.hpp"
#include <armadillo>


class Form//: public TrialFunction
{
public:

    Form()
    {
        ElementType=0;
        ElementNumber=0;
        GaussPntr=0;

       /* u_down=std::vector<TrialFunction> (u.MeshDimension-1);
        v_down=std::vector<TestFunctionGalerkin> (u.MeshDimension-1);
        libGmshReader::MeshReader Mesh(u.Msh->ElementData::fileName,u.MeshDimension-1);
        u_down[0]=TrialFunction(Mesh, u.vectorLvl-1);
        v_down[0]=TestFunctionGalerkin(u_down[0]);

        cout<<mat(inner(grad(v_down[0]),grad(u_down[0]))); */
    }

    Form(const TrialFunction& u): u_Internal(u)
    {
        ElementType=0;
        ElementNumber=0;
        GaussPntr=0;

       /* u_down=std::vector<TrialFunction> (u.MeshDimension-1);
        v_down=std::vector<TestFunctionGalerkin> (u.MeshDimension-1);
        libGmshReader::MeshReader Mesh(u.Msh->ElementData::fileName,u.MeshDimension-1);
        u_down[0]=TrialFunction(Mesh, u.vectorLvl-1);
        v_down[0]=TestFunctionGalerkin(u_down[0]);

        cout<<mat(inner(grad(v_down[0]),grad(u_down[0]))); */
    }

    void set_u_Internal(TrialFunction& u)
    {
        u_Internal=u;
    }

    inline sp_mat u(TrialFunction& u)
    {
        return u.Get_u(ElementType,GaussPntr);
    }

    inline sp_mat v(TestFunctionGalerkin& v)
    {
        return v.Get_v(ElementType, GaussPntr);
    }

    inline sp_mat grad(TrialFunction& u)
    {
        return u.Get_grad_u(ElementType, ElementNumber, GaussPntr);
    }

    inline sp_mat grad(TestFunctionGalerkin& v)
    {
        return v.Get_grad_v(ElementType, ElementNumber, GaussPntr);
    }

    inline sp_mat curl(TrialFunction& u)
    {
        return u.Get_curl_u(ElementType, ElementNumber, GaussPntr);
    }

    inline sp_mat curl(TestFunctionGalerkin& v)
    {
        return v.Get_curl_v(ElementType, ElementNumber, GaussPntr);
    }

    inline sp_mat inner(TestFunctionGalerkin& v, TrialFunction& u)
    {
        return v.Get_v(ElementType,GaussPntr)*u.Get_u(ElementType,GaussPntr);
    }

    inline sp_mat inner(TestFunctionGalerkin& v, sp_mat u)
    {
        return v.Get_v(ElementType,GaussPntr)*u;
    }

    inline sp_mat inner(sp_mat grad_v, sp_mat grad_u)
    {
        return grad_v*grad_u;
    }

    inline sp_mat dot(TestFunctionGalerkin& v, sp_mat grad_u)
    {
        return v.Get_v(ElementType,GaussPntr)*grad_u;
    }

    inline sp_mat dot(sp_mat grad_v, sp_mat grad_u)
    {
        return grad_v*grad_u;
    }

    inline sp_mat dot(vec a, sp_mat grad_u)
    {
        //return dot_vectrLvl_grad_u(a,grad_u);
        int& vectorLvl = u_Internal.vectorLvl;
        if (vectorLvl==1)
        {
            //mat aMatrx=repmat(a.t(),grad_u.n_rows,1);
            //return sp_mat(aMatrx%grad_u);
            mat aMatrx=a.rows(1,u_Internal.MeshDimension).t();
            //cout<<grad_u;
            return sp_mat(aMatrx*grad_u);
        }
        else
        {
            mat vctr=a(span(0,u_Internal.vectorLvl-1)).t();
            sp_mat vecMatrx(vctr.n_cols,vectorLvl*vctr.n_cols);
            for (int i=0; i<vectorLvl; i++)
            {
                //cout<<"col ="<<i*vectorLvl<<":"<<i*vectorLvl+vectorLvl-1;
                vecMatrx.cols(i*vectorLvl,i*vectorLvl+vectorLvl-1)=sp_mat(diagmat(vctr.t()));
            }
            return vecMatrx*grad_u;
        }

    }

    int ElementType;
    int ElementNumber;
    int GaussPntr;

protected:
    TrialFunction u_Internal;
    TestFunctionGalerkin v_Internal;
};



#endif // FORM_HPP
