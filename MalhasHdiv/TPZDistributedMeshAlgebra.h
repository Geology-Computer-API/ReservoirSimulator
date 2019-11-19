//
//  TPZDistributedMeshAlgebra.h
//  ReservoirSimulation
//
//  Created by Jorge Paúl Ordóñez Andrade on 05/11/19.
//

#ifndef TPZDistributedMeshAlgebra_h
#define TPZDistributedMeshAlgebra_h

#include "TPZCondensedMeshAlgebra.h"
#include "tpzverysparsematrix.h"
#include "pzblockdiag.h"
#include "pzmatred.h"
#include "pzbdstrmatrix.h"


class TPZDistributedMeshAlgebra
{
protected:
    //** @brief represents the condensable stiffness matrix
    //this data structure will be computed by a structmatrix type object at the end level
    //(when the are no children)
    //otherwise fStiffness = 0
    TPZMatrix<STATE> *fStiffness;
    //represents the right hand side
    TPZFMatrix<STATE> fRhs;
    //represents the block structure of the matrix
    TPZBlock<STATE> fBlock;
    //represents the residual
    TPZFMatrix<STATE> fResidual;
    //solution vector
    TPZFMatrix<STATE> fSolution;
    //block diagonal matrix used for block diagonal preconditioner
    TPZBlockDiagonal<STATE> fBlockDiagonal;
    //index vector of the degrees of freedom in the father mesh
    //the size of the vector determines the number of external equations
    TPZVec<int64_t> fFatherBlockIndexes;
    //indexes in the global vector
    TPZVec<int64_t> fGlobalIndexes;
    //pointer to the preconditioning matrix datastructure
//    TPZCondensedMeshAlgebra *fCondensed;
    //matrix to transfer the solution/rhs from the fine mesh to the coarse mesh
    TPZVerySparseMatrix<STATE> fTransferCoarseToFine;
    //pointer to the father datastructure
    TPZDistributedMeshAlgebra *fFather;
    //vector of leafs under this object
    TPZVec<TPZDistributedMeshAlgebra *> fChildren;
    
    
public:
    
    TPZDistributedMeshAlgebra();
    ~TPZDistributedMeshAlgebra();
    
    
    //compute the values of the blockdiagonal
    //fill in the block diagonal contribution to the father mesh
    void AssembleDiagonal(TPZBlockDiagonalStructMatrix &father_diagonal);
    //method for computing the residual Res=F-KU
    //returns the contribution to the father residual
    void ComputeResidual(TPZFMatrix<STATE> &father_rhs);
    
    //updates de solution based on the residual and block diagonal
    void BlockDiagonalSolve();
    
    //assemble the square of the residual
    double ResidualNormSquared();
    
    //copy the solution of the global vector
    void LoadSolution(TPZFMatrix<STATE> &global_sol);
    
    //export the residual vector
    void ExportResidual(TPZFMatrix<STATE> &global_residual);
    
};
#endif /* TPZDistributedMeshAlgebra_h */
