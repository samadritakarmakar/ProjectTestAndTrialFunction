#ifndef SYSTEMASSEMBLY_HPP
#define SYSTEMASSEMBLY_HPP
#include "LocalIntegration.hpp"
#include "FEMtools.h"
template<class GenericLocalIntegrator>
class SystemAssembler
{
public:
    SystemAssembler(Form<TrialFunction>& a, TrialFunction& u, TestFunctionGalerkin<TrialFunction>& v):
        a_Internal(a), u_Internal(u), v_Internal(v)
    {
        integrate=std::
                shared_ptr<LocalIntegrator<TrialFunction>>(new LocalIntegrator<TrialFunction>(a_Internal,u_Internal,v_Internal));
    }

    void SetLocalIntegrator(std::shared_ptr<GenericLocalIntegrator>& Integrate)
    {
        integrate=Integrate;
    }

    void RunLocalIntegration()
    {
        integrate->local_intergrator();
    }

    void SetMatrixSize(sp_mat& A)
    {
        A_Internal=&A;
        A_Internal->set_size
                (v_Internal.vectorLvl*u_Internal.Msh->NodeTag.n_rows,u_Internal.vectorLvl*u_Internal.Msh->NodeTag.n_rows);
        //cout<<"No of Rows of NodeTags= "<<u_Internal.Msh->NodalCoordinates;
    }

    void RunSystemAssembly(sp_mat& A)
    {
        SetMatrixSize(A);
        cout<<"Size of A= "<<A.n_rows<<","<<A.n_cols<<"\n";
        NodePositions=std::vector<umat>(u_Internal.NoOfElementTypes);

        //Configuration for Batch addition of matrix to global matrix
        bool add_values=false;
        bool sort_locations=true;
        bool check_for_zeros=true;
        //------------------------------------------------------------
        for (int ElmntTyp=0; ElmntTyp<u_Internal.NoOfElementTypes; ElmntTyp++)
        {
            GetNodePostions(NodePositions[ElmntTyp], u_Internal.Msh->ElementNodes[ElmntTyp], u_Internal.vectorLvl);
            //u_Internal.Msh->ElementNodes[ElmntTyp].n_rows;
            for (int ElmntNmbr=0; ElmntNmbr<u_Internal.Msh->ElementNodes[ElmntTyp].n_rows; ElmntNmbr++)
            {
                RunLocalIntegration();
                umat positions=NodePositions[ElmntTyp].row(a_Internal.ElementNumber);
                //cout<<"positions= "<<positions<<"\n";
                umat locations(2,positions.n_cols*positions.n_cols);
                int locationPtr=0;
                for (int row=0;row<positions.n_cols;row++)
                {
                    for (int col=0;col<positions.n_cols;col++)
                    {
                        umat locations1;
                        locations1<<positions(row)<<endr<<positions(col)<<endr;
                        locations.col(locationPtr)=locations1;
                        locationPtr++;
                    }
                }
                vec values=vectorise(mat(a_Internal.ResultingMat.t()));
                //cout<<"Size of locations ="<<locations.n_rows<<","<<locations.n_cols<<"\n";
                //cout<<"Size of values ="<<values.n_rows<<","<<values.n_cols<<"\n";
                sp_mat A_temp=sp_mat(add_values, locations, values, A.n_rows, A.n_cols, sort_locations, check_for_zeros);
                A=A+A_temp;
                /*cout<<mat(a_Internal.ResultingMat);
                for(int lcl=0; lcl<locations.n_cols; lcl++)
                {
                    cout<<A(locations(0,lcl),locations(1,lcl))<<" ";
                    if ((lcl+1)%positions.n_cols==0)
                        cout<<"\n";
                }*/
                a_Internal.NextElementNumber();
            }
            a_Internal.NextElementType();
        }

    }


private:
    Form<TrialFunction>& a_Internal;
    TrialFunction& u_Internal;
    TestFunctionGalerkin<TrialFunction>& v_Internal;
    std::shared_ptr<LocalIntegrator<TrialFunction>> integrate;
    std::vector<umat> NodePositions;
    sp_mat* A_Internal;
};

#endif // SYSTEMASSEMBLY_HPP
