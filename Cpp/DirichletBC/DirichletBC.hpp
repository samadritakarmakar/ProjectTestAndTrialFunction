#ifndef DIRICHLETBC_HPP
#define DIRICHLETBC_HPP
#include "TrialFunction.hpp"
#include "libGmshReader.h"
class DirichletBC
{
public:
    DirichletBC(TrialFunction &u, int PhysicalGroupNumber, umat& boolDiricletNodes):
        u_Internal(u), PhysclGrpNum (PhysicalGroupNumber), Msh(u.Msh),
        vctrLvlInternal(u.originalVctrLvl), currentNodeNumber(0)
    {
        DirichletBCNodes=std::vector<umat>(Msh->NumOfElementTypes);
        DirichletBCNodesFill=std::vector<umat>(Msh->NumOfElementTypes);
        NoOfNodes=std::vector<int>(Msh->NumOfElementTypes);
        NodePositions=std::vector<umat>(Msh->NumOfElementTypes);
        for (int ElementType=0; ElementType<Msh->NumOfElementTypes; ElementType++)
        {
            DirichletBCNodes[ElementType]=unique(Msh->ElmntPhysclGrpNodes[ElementType][PhysclGrpNum]);
            NoOfNodes[ElementType]=DirichletBCNodes[ElementType].n_rows;
            GetNodePostions(NodePositions[ElementType],DirichletBCNodes[ElementType],vctrLvlInternal);
            DirichletBCNodesFill[ElementType].set_size(NodePositions[ElementType].n_rows,1);
            int DiricletNodeSize=1;
            int boolDiriclet=0;
            for (int colNo=0;colNo<NodePositions[ElementType].n_cols;colNo++)
            {
                if(bool(boolDiricletNodes(0,boolDiriclet)))
                {
                    DirichletBCNodesFill[ElementType].resize(NodePositions[ElementType].n_rows, DiricletNodeSize);
                    DirichletBCNodesFill[ElementType].col(DiricletNodeSize-1)=NodePositions[ElementType].col(colNo);
                    DiricletNodeSize++;
                }
                boolDiriclet++;
                if(boolDiriclet>=vctrLvlInternal)
                {
                    boolDiriclet=0;
                }
            }
        }
        //cout<<NodePositions[0];
        //cout<<DirichletBCNodesFill[0];
    }

    virtual mat Expression(mat& x)
    {
        mat Exprssn=zeros(vctrLvlInternal,1);
        return Exprssn;
    }

    ///Returns a matrix of coordinates [x, y, z] of all current node number
    /// within the Physical Enitity.
     mat x()
    {
         for (int ElementType=0; ElementType<Msh->NumOfElementTypes; ElementType++)
         {
             mat Coordinates=Msh->NodalCoordinates.rows(DirichletBCNodes[ElementType]);
             return Coordinates.cols(currentNodeNumber,vctrLvlInternal-1);
         }
     }

     void ApplyBC(sp_mat &A, mat& b)
     {
         mat x={0};
         for (int ElementType=0; ElementType<Msh->NumOfElementTypes; ElementType++)
         {
             umat NodePosition;
             GetNodePostions(NodePosition, DirichletBCNodes[ElementType], vctrLvlInternal);
             for (int NodeNumber=0; NodeNumber<NoOfNodes[ElementType]; NodeNumber++)
             {
                 currentNodeNumber=NodeNumber;
                 umat Positions1=NodePosition.row(NodeNumber);
                 b.rows(Positions1)=Expression(x);

             }
         }
     }
private:
    TrialFunction& u_Internal;
    int& PhysclGrpNum;
    int& vctrLvlInternal;
    int currentNodeNumber;
    libGmshReader::MeshReader *Msh;
    std::vector<umat> DirichletBCNodes, DirichletBCNodesFill;
    std::vector<int> NoOfNodes;
    std::vector<umat> NodePositions;
};

#endif // DIRICHLETBC_HPP
