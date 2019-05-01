#ifndef TRIALFUNCTION_HPP
#define TRIALFUNCTION_HPP
#include "libGmshReader.h"
#include "FEMtools.h"
#include "LagrangeShapeFunction.hpp"
#include "Jacobian.hpp"
#include <vector>
#include <armadillo>
#include <assert.h>
using namespace arma;
class OtherData
{
public:
    int ElementType=0, GuassPntCounter=0;
    //libGmshReader::MeshReader Mesh;
};
class TrialFunction
{
public:
    mat (TrialFunction::*Pntr_Calc_N)(mat , OtherData);
    std::vector<int> NoOfGaussPts;
    TrialFunction(libGmshReader::MeshReader& Mesh, int vectorLevel)
    {
        Msh=&Mesh;
        vectorLvl=vectorLevel;
        N=std::vector<mat> (Msh->NumOfElementTypes);
        Phi=std::vector<LagrangeShapeFunction> (Msh->NumOfElementTypes);
        GaussData=std::vector<mat> (Msh->NumOfElementTypes);
        Weight=std::vector<mat> (Msh->NumOfElementTypes);
        GaussPointx=std::vector<mat> (Msh->NumOfElementTypes);
        GaussPointy=std::vector<mat> (Msh->NumOfElementTypes);
        GaussPointz=std::vector<mat> (Msh->NumOfElementTypes);
        NoOfGaussPts=std::vector<int> (Msh->NumOfElementTypes);
        u=std::vector<std::vector<sp_mat>>(Msh->NumOfElementTypes);
        dN_by_dEps=std::vector<std::vector<mat>>(Msh->NumOfElementTypes);
        for (int ElementType = 0; ElementType<Msh->NumOfElementTypes; ++ElementType)
        {
            Generate_GaussPoints_Weights_ShapeFunctions(ElementType);
            NoOfGaussPts[ElementType]=GaussPointx[ElementType].n_rows;
            u[ElementType]=std::vector<sp_mat>(NoOfGaussPts[ElementType]);
            dN_by_dEps[ElementType]=std::vector<mat>(NoOfGaussPts[ElementType]);
            /*N[ElementType]=Phi[ElementType].GetShapeFunction(GaussPointx[ElementType],GaussPointy[ElementType],
                                                             GaussPointz[ElementType]);*/
            for (int GaussPt=0;GaussPt<NoOfGaussPts[ElementType];GaussPt++)
            {
                mat Ncol=N[ElementType].col(GaussPt);
                u[ElementType][GaussPt]=vectorizeQuantity(Ncol,vectorLvl);
                //cout<<mat(u[ElementType][GaussPt])<<"\n";
            }
        }
        Generate_dN_by_dEps ();
    }

    void Get_F(int ElementType, int ElementNumber, int GaussPntr, mat& F)
    {
        umat NodesAtElmntNmbr=Msh->ElementNodes[ElementType].row(ElementNumber);
        //cout<<"Gmsh Node Tags are "<<Msh->GmshNodeTag[ElementType].row(ElementNumber);
        //cout<<"ElementNodes =\n"<<NodesAtElmntNmbr;
        //cout<<"Coodinates of "<<ElementNumber<<" are \n"<<Msh->NodalCoordinates.rows(NodesAtElmntNmbr);
        mat Coordinates=Msh->NodalCoordinates.rows(NodesAtElmntNmbr);
        F=Coordinates.t()*dN_by_dEps[ElementType][GaussPntr];
    }

    void Generate_dN_by_dEps ()
    {
        OtherData Data;
        int dim=Msh->ElementData::dim;
        mat x;
        x.set_size(1,dim);
        for (int ElmntTypCountr=0; ElmntTypCountr<Msh->NumOfElementTypes; ElmntTypCountr++)
        {
            int &i=ElmntTypCountr;
            Data.ElementType=ElmntTypCountr;
            for (int GaussPntr=0; GaussPntr<GaussPointx[i].n_rows; GaussPntr++)
            {
                if(Msh->ElementData::dim==1)
                {
                    x=GaussPointx[i].row(GaussPntr);
                }
                else if (Msh->ElementData::dim==2)
                {
                    x=join_vert(GaussPointx[i].row(GaussPntr), GaussPointy[i].row(GaussPntr));
                }
                else if (Msh->ElementData::dim==3)
                {
                    mat x1=join_vert(GaussPointx[i].row(GaussPntr), GaussPointy[i].row(GaussPntr));
                    x=join_vert(x1,GaussPointz[i].row(GaussPntr));
                }
                Pntr_Calc_N=&TrialFunction::Calculate_N_GaussPointWise;
                dN_by_dEps[i][GaussPntr]= Jacobian(this,(this->Pntr_Calc_N), x, Data);
                //cout<<"dN_by_dEps at Gauss Pnt "<<GaussPntr<<" =\n"<<dN_by_dEps[i][GaussPntr];
            }
        }
    }

    inline mat Calculate_N_GaussPointWise(mat x1, OtherData Data)
    {
        vec x=vectorise(x1);
        //int &ElemntTyp=Data.ElementType;
        int &GaussPt=Data.GuassPntCounter;
        mat GaussPtx,GaussPty,GaussPtz, N_AtGaussPt;
        if (Msh->ElementData::dim==1)
        {
            //cout<<"x= "<<x<<"\n";
            GaussPtx=x.row(0);
            N_AtGaussPt=Phi[Data.ElementType].GetShapeFunction(GaussPtx);
        }
        else if (Msh->ElementData::dim==2)
        {
            //cout<<"x= "<<x<<"\n";
            GaussPtx=x.row(0);
            GaussPty=x.row(1);
            N_AtGaussPt=Phi[Data.ElementType].GetShapeFunction(GaussPtx,GaussPty);
        }
        else if(Msh->ElementData::dim=3)
        {
            //cout<<"x= "<<x<<"\n";
            GaussPtx=x.row(0);
            GaussPty=x.row(1);
            GaussPtz=x.row(2);
            N_AtGaussPt=Phi[Data.ElementType].GetShapeFunction(GaussPtx,GaussPty,GaussPtz);
        }
        //sp_mat u_AtGaussPt=vectorizeQuantity(N_AtGaussPt, vectorLvl);
        return N_AtGaussPt;
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
    std::vector<std::vector<mat>> dN_by_dEps;
    std::vector<LagrangeShapeFunction> Phi;
    std::vector<mat> N;
    libGmshReader::MeshReader *Msh;

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

    void Generate_GaussPoints_Weights_ShapeFunctions(int ElementType)
    {
        Phi[ElementType]=LagrangeShapeFunction(*Msh,ElementType);
        GaussData[ElementType].load(FEMtools::LoadGaussFile(*Msh,ElementType));
        if(Msh->GmshElementNameOnly[ElementType].compare("Line")==0)
        {
            Weight[ElementType]=GaussData[ElementType].col(1);
            GaussPointx[ElementType]=GaussData[ElementType].col(2);
            N[ElementType]=Phi[ElementType].GetShapeFunction(GaussPointx[ElementType]);
        }
        else if(Msh->GmshElementNameOnly[ElementType].compare("Triangle")==0)
        {
            Weight[ElementType]=GaussData[ElementType].col(1);
            GaussPointx[ElementType]=GaussData[ElementType].col(2);
            GaussPointy[ElementType]=GaussData[ElementType].col(3);
            N[ElementType]=Phi[ElementType].GetShapeFunction(GaussPointx[ElementType],GaussPointy[ElementType]);
        }
        else if (Msh->GmshElementNameOnly[ElementType].compare("Quadrilateral")==0)
        {
            CalculateGaussPointsAndWeights(Weight[ElementType],GaussPointx[ElementType],GaussPointy[ElementType],ElementType);
            N[ElementType]=Phi[ElementType].GetShapeFunction(GaussPointx[ElementType],GaussPointy[ElementType]);
        }
        else if (Msh->GmshElementNameOnly[ElementType].compare("Tetrahedron")==0)
        {
            Weight[ElementType]=GaussData[ElementType].col(1);
            GaussPointx[ElementType]=GaussData[ElementType].col(2);
            GaussPointy[ElementType]=GaussData[ElementType].col(3);
            GaussPointz[ElementType]=GaussData[ElementType].col(4);
            N[ElementType]=Phi[ElementType].GetShapeFunction(GaussPointx[ElementType],GaussPointy[ElementType],
                                                             GaussPointz[ElementType]);
        }
        else if(Msh->GmshElementNameOnly[ElementType].compare("Hexahedron")==0
                ||Msh->GmshElementNameOnly[ElementType].compare("Prism")==0
                ||Msh->GmshElementNameOnly[ElementType].compare("Pyramid")==0 )
        {
            CalculateGaussPointsAndWeights(Weight[ElementType],GaussPointx[ElementType],GaussPointy[ElementType],
                                           GaussPointz[ElementType],ElementType);
            N[ElementType]=Phi[ElementType].GetShapeFunction(GaussPointx[ElementType],GaussPointy[ElementType],
                                                             GaussPointz[ElementType]);
        }
    }
};

#endif // TRIALFUNCTION_HPP
