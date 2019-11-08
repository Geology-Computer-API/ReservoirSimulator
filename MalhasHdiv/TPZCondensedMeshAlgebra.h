//
//  TPZCondensedMeshAlgebra.h
//  ReservoirSimulation
//
//  ** Data structure of the preconditioner **
//
//  Created by Jorge Paúl Ordóñez Andrade on 05/11/19.
//

#ifndef TPZCondensedMeshAlgebra_h
#define TPZCondensedMeshAlgebra_h

#include "pzmatrix.h"
#include "pzblock.h"
#include "pzmatred.h"
#include "pzelmat.h"
#include "pzvec.h"
#include "TPZDistributedMeshAlgebra.h"



class TPZCondensedMeshAlgebra
{
protected:
    //represents the condensable stiffness matrix
    //this data structure will be computed by a structmatrix type object
    TPZMatRed<REAL, TPZFMatrix<REAL>> fStiffness;
    //block structure of the stiffness
    TPZBlock<STATE> fBlock;
    //index vector of the degrees of freedom in the father mesh
    //the size of the vector determines the number of external equations
    TPZVec<int64_t> fFatherBlockIndexes;
    //pointer to the father datastructure
    TPZCondensedMeshAlgebra *fFather;
    //vector of leafs under this object
    TPZVec<TPZCondensedMeshAlgebra *> fChildren;
    //the fine mesh object
//    TPZDistributedMeshAlgebra *fReference;
    
public:
    /** @brief default constructor  */
    TPZCondensedMeshAlgebra();
    
    /** @brief default desconstructor  */
    ~TPZCondensedMeshAlgebra();
    
    /** @brief Copy constructor */
    TPZCondensedMeshAlgebra(const TPZCondensedMeshAlgebra &cp);
    
    //compute the condensed stiffness matrix using the fStiffness data structure
    //this method will be called by the fFather object to assemble its stiffness matrix
    void CalcStiff(TPZElementMatrix &ek, TPZElementMatrix &ef);
    
    //method for assembling the fStiffness object if the object has a children
    //the method calls assemble CalcStiff for all children and assembles the stiffness and rhs
    void Assemble();
    
    //method for condense/assembling the rhs
    //the method will reduce the rhs using the known stiffness matrix
    void AssembleRhs(TPZFMatrix<STATE> &father_rhs);
    
    //method for computing the solution of the internal degrees of freedom
    void LoadSolution(TPZMatrix<STATE> &father_solution);
    
    //compute the stiffness matrix by selecting elements from the fine mesh matrix
    //build the coarse mesh stiffness matrix from the fine mesh stiffness matrix
    void TransferStiffness();
    
};

#endif /* TPZCondensedMeshAlgebra_h */
