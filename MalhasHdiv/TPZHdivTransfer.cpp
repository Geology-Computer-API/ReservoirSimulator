//
//  TPZHdivTransfer.cpp
//  reservoirlib
//
//  Created by Jorge Paúl Ordóñez Andrade on 15/07/19.
//

#include "TPZHdivTransfer.h"


template<class TVar>
TPZHdivTransfer<TVar>::TPZHdivTransfer() : TPZMatrix<TVar>() {
    DebugStop();
}


template<class TVar>
TPZHdivTransfer<TVar>::TPZHdivTransfer(const TPZHdivTransfer<TVar> &cp) : TPZMatrix<TVar>(cp){
    DebugStop();
}

template<class TVar>
TPZHdivTransfer<TVar>::TPZHdivTransfer(TPZVec<int64_t> &Indexes) : TPZMatrix<TVar>() {
    fIndexes = Indexes;
}

~TPZHdivTransfer(){
    
}

template<class TVar>
void TPZHdivTransfer<TVar>::Print(const char *name = NULL, std::ostream &out = std::cout , const MatrixOutputFormat form = EFormatted) const {
    DebugStop();
}

/**
 * @brief Gather a vector
 * @param y: Compressed vector which stores the information of the full vector
 * @param x: Vector which is going to provide the values for x
 */
template<class TVar>
void TPZHdivTransfer<TVar>::Gather(TPZFMatrix<TVar> &y, TPZFMatrix<TVar> &x){
    DebugStop();
}

/**
 * @brief Scatter a vector
 * @param y: Full vector which stores the information of the compressed vector
 * @param x: Vector which is going to provide the values for y
 */
template<class TVar>
void TPZHdivTransfer<TVar>::Scatter(TPZFMatrix<TVar> &y, TPZFMatrix<TVar> &x){
    DebugStop();
}


template<class TVar>
void TPZHdivTransfer<TVar>::SetIndexes(TPZVec<int64_t> &Indexes){
    fIndexes = Indexes;
}


template<class TVar>
TPZVec<int64_t> & TPZHdivTransfer<TVar>::GetIndexes(){
    return fIndexes;
}
