//
//  TPZHdivTransfer.h
//  reservoirlib
//
//  Created by Jorge Paúl Ordóñez Andrade on 15/07/19.
//

#ifndef TPZHdivTransfer_h
#define TPZHdivTransfer_h

#include <stdio.h>
#include "pzvec.h"
#include "pzmatrix.h"

template<class TVar>
class TPZHdivTransfer : public TPZMatrix<TVar> {
    
    TPZVec<int64_t> fIndexes;
    
public:
    
    /** @brief Default constructor */
    TPZHdivTransfer();
    
    /** @brief Copy constructor */
    TPZHdivTransfer(const TPZHdivTransfer<TVar> &cp);
    
    /** @brief Constructor based on indexes */
    TPZHdivTransfer(TPZVec<int64_t> &Indexes);
    
    /** @brief Default constructor */
    ~TPZHdivTransfer();
    
    /** @brief Copy constructor */
    virtual void Print(const char *name = NULL, std::ostream &out = std::cout , const MatrixOutputFormat form = EFormatted) const override;
    
    /** @brief Gather y elements into x vector */
    void Gather(TPZFMatrix<TVar> &y, TPZFMatrix<TVar> &x);
    
    /** @brief Scatter y elements into x vector */
    void Scatter(TPZFMatrix<TVar> &y, TPZFMatrix<TVar> &x);
    
};

#endif /* TPZHdivTransfer_h */
