//
//  TPZDistributedMeshAlgebra.cpp
//
//  Created by Jorge Paúl Ordóñez Andrade on 05/11/19.
//

#include "TPZDistributedMeshAlgebra.h"

TPZDistributedMeshAlgebra::TPZDistributedMeshAlgebra(){
    
}
TPZDistributedMeshAlgebra::~TPZDistributedMeshAlgebra(){
    
}

void TPZDistributedMeshAlgebra::AssembleDiagonal(TPZBlockDiagonalStructMatrix &father_diagonal){
    
    
}

void TPZDistributedMeshAlgebra::ComputeResidual(TPZFMatrix<STATE> &father_rhs){
     DebugStop();
}
void TPZDistributedMeshAlgebra::BlockDiagonalSolve(){
     DebugStop();
}
double TPZDistributedMeshAlgebra::ResidualNormSquared(){
     DebugStop();
    return 0;
}
void TPZDistributedMeshAlgebra::LoadSolution(TPZFMatrix<STATE> &global_sol){
     DebugStop();
}
void TPZDistributedMeshAlgebra::ExportResidual(TPZFMatrix<STATE> &global_residual){
     DebugStop();
}
