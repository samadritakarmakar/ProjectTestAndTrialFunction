#ifndef TRIALFUNCTION_HPP
#define TRIALFUNCTION_HPP
#include "libGmshReader.h"
#include "FEMtools.h"
#include "LagrangeShapeFunction.hpp"
#include <vector>
#include <armadillo>
#include <assert.h>
using namespace arma;
class TrialFunction
{
public:
    TrialFunction(libGmshReader::MeshReader& Mesh, int vectorLevel)
    {
        vectorLvl=vectorLevel;
        N=std::vector<mat> (Mesh.NumOfElementTypes);
        Phi=std::vector<LagrangeShapeFunction> (Mesh.NumOfElementTypes);
        GaussData=std::vector<mat> (Mesh.NumOfElementTypes);
        Weight=std::vector<mat> (Mesh.NumOfElementTypes);
        GaussPointx=std::vector<mat> (Mesh.NumOfElementTypes);
        GaussPointy=std::vector<mat> (Mesh.NumOfElementTypes);
        GaussPointz=std::vector<mat> (Mesh.NumOfElementTypes);
        u=std::vector<std::vector<sp_mat>>(Mesh.NumOfElementTypes);
        for (int ElementType = 0; ElementType<Mesh.NumOfElementTypes; ++ElementType)
        {
            Generate_GaussPoints_Weights_ShapeFunctions(Mesh,ElementType);
            int NoOfGaussPts=GaussPointx[ElementType].n_rows;
            u[ElementType]=std::vector<sp_mat>(NoOfGaussPts);
            for (int GaussPt=0;GaussPt<NoOfGaussPts;GaussPt++)
            {
                mat Ncol=N[ElementType].col(GaussPt);
                u[ElementType][GaussPt]=vectorizeQuantity(Ncol,vectorLvl);
                //cout<<mat(u[ElementType][GaussPt])<<"\n";
            }
        }
    }
protected:
    /// 'Variable' must be in form of a vector
    /// This function builds the representation similar to
    /// how shape functions are represented generally in text books
    inline sp_mat vectorizeQuantity(mat& Variable, int VectorLevel)
    {
        vec VariableVector=vectorise(Variable);
        int NoOfRows=VariableVector.n_rows;
        sp_mat I=speye(VectorLevel, VectorLevel);
        //cout<<"I=\n"<<mat(I)<<"\n";
        //cout<<"N=\n"<<Variable<<"\n";
        sp_mat Matrx(VectorLevel, NoOfRows*VectorLevel);
        for (int i=0;i<NoOfRows;i++)
        {
            Matrx.cols(VectorLevel*i,VectorLevel*i+(VectorLevel-1))=VariableVector(i)*I;
        }
        return Matrx;
    }
    int vectorLvl;
    std::vector<mat> GaussData;
    std::vector<mat> Weight;
    std::vector<mat> GaussPointx;
    std::vector<mat> GaussPointy;
    std::vector<mat> GaussPointz;
    std::vector<std::vector<sp_mat>> u;
    std::vector<LagrangeShapeFunction> Phi;
    std::vector<mat> N;

    void CalculateGaussPointsAndWeights(mat &wt, mat& GaussPntx, mat& GaussPnty, int ElementType)
    {
        wt.set_size(pow(GaussData[ElementType].n_rows,2),1);
        //cout<<"Wt size =("<<wt.n_rows<<","<<wt.n_cols<<")\n";
        GaussPntx.set_size(wt.n_rows,1);
        GaussPnty.set_size(wt.n_rows,1);
        int pos=0;
        for(int i=0; i<GaussData[ElementType].n_rows; i++)
        {
            for (int j=0; j<GaussData[ElementType].n_rows; j++)
            {
                wt(pos,0)=GaussData[ElementType](i,1)*GaussData[ElementType](j,1);
                GaussPntx(pos,0)=GaussData[ElementType](i,2);
                GaussPnty(pos,0)=GaussData[ElementType](j,2);
                //cout<<"i= "<<i<<" j= "<<j<<" pos="<<pos<<"\n";
                pos++;
            }
        }
    }

    void CalculateGaussPointsAndWeights(mat &wt, mat& GaussPntx, mat& GaussPnty, mat& GaussPntz, int ElementType)
    {
        wt.set_size(pow(GaussData[ElementType].n_rows,3),1);
        GaussPntx.set_size(wt.n_rows,1);
        GaussPnty.set_size(wt.n_rows,1);
        GaussPntz.set_size(wt.n_rows,1);
        int pos=0;
        for(int i=0; i<GaussData[ElementType].n_rows; i++)
        {
            for (int j=0; j<GaussData[ElementType].n_rows; j++)
            {
                for (int k=0; k<GaussData[ElementType].n_rows; k++)
                {
                    wt(pos,0)=GaussData[ElementType](i,1)*GaussData[ElementType](j,1)*GaussData[ElementType](k,1);
                    GaussPntx(pos,0)=GaussData[ElementType](i,2);
                    GaussPnty(pos,0)=GaussData[ElementType](j,2);
                    GaussPntz(pos,0)=GaussData[ElementType](k,2);
                    pos++;
                }
            }
        }
    }

    void Generate_GaussPoints_Weights_ShapeFunctions(libGmshReader::MeshReader &Mesh, int i)
    {
        Phi[i]=LagrangeShapeFunction(Mesh,i);
        GaussData[i].load(FEMtools::LoadGaussFile(Mesh,i));
        if(Mesh.GmshElementNameOnly[i].compare("Line")==0)
        {
            Weight[i]=GaussData[i].col(1);
            GaussPointx[i]=GaussData[i].col(2);
            N[i]=Phi[i].GetShapeFunction(GaussPointx[i]);
        }
        else if(Mesh.GmshElementNameOnly[i].compare("Triangle")==0)
        {
            Weight[i]=GaussData[i].col(1);
            GaussPointx[i]=GaussData[i].col(2);
            GaussPointy[i]=GaussData[i].col(3);
            N[i]=Phi[i].GetShapeFunction(GaussPointx[i],GaussPointy[i]);
        }
        else if (Mesh.GmshElementNameOnly[i].compare("Quadrilateral")==0)
        {
            CalculateGaussPointsAndWeights(Weight[i],GaussPointx[i],GaussPointy[i],i);
            N[i]=Phi[i].GetShapeFunction(GaussPointx[i],GaussPointy[i]);
        }
        else if (Mesh.GmshElementNameOnly[i].compare("Tetrahedron")==0)
        {
            Weight[i]=GaussData[i].col(1);
            GaussPointx[i]=GaussData[i].col(2);
            GaussPointy[i]=GaussData[i].col(3);
            GaussPointz[i]=GaussData[i].col(4);
            N[i]=Phi[i].GetShapeFunction(GaussPointx[i],GaussPointy[i],GaussPointz[i]);
        }
        else if(Mesh.GmshElementNameOnly[i].compare("Hexahedron")==0
                ||Mesh.GmshElementNameOnly[i].compare("Prism")==0
                ||Mesh.GmshElementNameOnly[i].compare("Pyramid")==0 )
        {
            CalculateGaussPointsAndWeights(Weight[i],GaussPointx[i],GaussPointy[i],GaussPointz[i],i);
            N[i]=Phi[i].GetShapeFunction(GaussPointx[i],GaussPointy[i],GaussPointz[i]);
        }
    }
};

#endif // TRIALFUNCTION_HPP
