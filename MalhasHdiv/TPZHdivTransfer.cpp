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

template<class TVar>
void TPZHdivTransfer<TVar>::Gather(TPZFMatrix<TVar> &y, TPZFMatrix<TVar> &x){
    DebugStop();
}

template<class TVar>
void TPZHdivTransfer<TVar>::Scatter(TPZFMatrix<TVar> &y, TPZFMatrix<TVar> &x){
    DebugStop();
}
