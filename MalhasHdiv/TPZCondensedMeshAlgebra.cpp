//
//  TPZCondensedMeshAlgebra.cpp
//
//  Created by Jorge Paúl Ordóñez Andrade on 05/11/19.
//

#include "TPZCondensedMeshAlgebra.h"

TPZCondensedMeshAlgebra::TPZCondensedMeshAlgebra(){
        
}
TPZCondensedMeshAlgebra::~TPZCondensedMeshAlgebra(){
    
}

void TPZCondensedMeshAlgebra::CalcStiff(TPZElementMatrix &ek, TPZElementMatrix &ef){
     DebugStop();
}
void TPZCondensedMeshAlgebra::Assemble(){
    DebugStop();
    
}
void TPZCondensedMeshAlgebra::AssembleRhs(TPZFMatrix<STATE> &father_rhs){
     DebugStop();
}
//void TPZCondensedMeshAlgebra::LoadSolution(TPZMatrix<STATE> $father_solution){
//     DebugStop();
//}
void TPZCondensedMeshAlgebra::TransferStiffness(){
     DebugStop();
}
