//
//  ConfigurateCase.h
//  ReservoirSimulation
//
//  Created by Jorge Paúl Ordóñez Andrade on 10/11/19.

#include "pzcmesh.h"
#include "pzintel.h"
#include "pzcompel.h"
#include "pzfmatrix.h"
#include "pzvec.h"
#include "pzmanvector.h"
#include "pzstack.h"

#include "pzanalysis.h"
#include "pzfstrmatrix.h"
#include "pzskylstrmatrix.h"
#include "pzstepsolver.h"
#include "pzgeopyramid.h"
#include "TPZGeoLinear.h"
#include "TPZMHMixedMeshControl.h"
#include <TPZRefPattern.h>

#include "TPZMaterial.h"
//#include "pzelasmat.h"
//#include "pzlog.h"
#include "pzgengrid.h"

//#include <time.h>
//#include <stdio.h>

#include <math.h>

#include <iostream>
#include <fstream>
#include <string.h>
#include <sstream>

#include <set>
#include <map>
#include <vector>
#include "TPZRefPatternDataBase.h"
#include "TPZRefPatternTools.h"
#include "TPZVTKGeoMesh.h"
#include "TPZExtendGridDimension.h"

//#include "TPZMatLaplacian.h"
#include "pzpoisson3d.h"
#include "TPZNullMaterial.h"
#include "pzbndcond.h"
#include "TPZSSpStructMatrix.h"
#include "pzskylstrmatrix.h"
#include "TPZSkylineNSymStructMatrix.h"

#include "pzanalysis.h"
#include "pzfstrmatrix.h"
#include "pzskylstrmatrix.h"
#include "pzsolve.h"
#include "TPZPersistenceManager.h"
#include "pzbuildmultiphysicsmesh.h"
#include "TPZMultiphysicsCompMesh.h"
#include "TPZCompMeshTools.h"
#include "TPZMixedDarcyWithFourSpaces.h"
#include "pzcmesh.h"
#include "pzbdstrmatrix.h"

#include "TPZHdivTransfer.h"
#include "pzseqsolver.h"
#include "pzmgsolver.h"
#include "TPZTimer.h"

#include "TPZCondensedMeshAlgebra.h"
#include "TPZDistributedMeshAlgebra.h"
#include "RSimulationCase.h"
#include "pzbuildmultiphysicsmesh.h"
#include "TPZMHMixedMesh4SpacesControl.h"

#ifndef ConfigurateCase_h
#define ConfigurateCase_h

class ConfigurateCase {
    TPZGeoMesh *m_gmesh;
    SimulationCase fsim_case;
    int m_fineorder = 2;
    int m_coarseorder = 1;
    public:
    
    ConfigurateCase();
    
    TPZGeoMesh * CreateUniformMesh(int nx, REAL L, int ny=0, REAL h=0, int nz=0, REAL w=0);
    
    TPZCompMesh * HDivMesh(TPZGeoMesh * gmesh, int orderfine, int ordercoarse);
    TPZCompMesh * DiscontinuousMesh(TPZGeoMesh * gmesh,int order, int lagrangemult);
    TPZMultiphysicsCompMesh *CreateMultCompMesh();
    TPZAnalysis *CreateAnalysis(TPZMultiphysicsCompMesh *mcmesh);
    
    void SetSimulationCase(SimulationCase simcase){
        fsim_case =simcase;
    }
    
    void SetFineOrder(int fineorder){
        m_fineorder= fineorder;
    }
    
    int GetFineOrder(int fineorder){
        return m_fineorder;
    }
    
    void SetCoarseOrder(int coarseorder){
        m_coarseorder= coarseorder;
    }
    
    int GetCoarseOrder(int coarseorder){
        return m_coarseorder;
    }
    void SetGmesh(TPZGeoMesh *gmesh){
        m_gmesh= gmesh;
    }
    
    TPZGeoMesh* GetGmesh(){
        return m_gmesh;
    }
    SimulationCase GetSimulationCase(){
        if(fsim_case.omega_ids.size() == 0){
            std::cout<<"Please set simulation case!"<<std::endl;
        }
        return fsim_case;
    }
   TPZAutoPointer<TPZMHMixedMesh4SpacesControl> CreateMHMMixedMesh();
//    void ComputeCoarseIndices(TPZGeoMesh *gmesh, TPZVec<int64_t> &coarseindices);
    void InsertMaterialObjects(TPZMHMeshControl &control);
    TPZGeoMesh *CreateMHMGeoMesh(int ncoarse_elx, int ncooars_eley, int nfine_elx, int nfine_ely);
    void CreateRefPattern();
    TPZGeoMesh *CreateGeowithRefPattern();
    
};


#endif /* ConfigurateCase_h */
