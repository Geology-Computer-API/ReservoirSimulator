//$Id: main.cc,v 1.21 2010-07-22 17:43:43 caju Exp $
#ifdef HAVE_CONFIG_H
#include <pz_config.h>
#endif

#include "pzgeoel.h"
#include "pzgnode.h"
#include "pzgmesh.h"
#include "pzbiharmonic.h"
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

#include <TPZRefPattern.h>

#include "TPZMaterial.h"
#include "pzelasmat.h"
#include "pzlog.h"

#include "pzgengrid.h"

#include <time.h>
#include <stdio.h>

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

#include <opencv2/opencv.hpp>

#include "TRSLinearInterpolator.h"
#include "TPZMatLaplacian.h"
#include "pzpoisson3d.h"
#include "mixedpoisson.h"
#include "pzbndcond.h"
#include "TPZSSpStructMatrix.h"
#include "pzskylstrmatrix.h"
#include "TPZSkylineNSymStructMatrix.h"
#include "TPZHybridizeHDiv.h"
#include "pzbuildmultiphysicsmesh.h"
#include "TPZMultiphasicFlowAnalysis.h"
#include "TPZMixedDarcyFlow.h"
#include "TPZPrimalPoisson.h"
#include "TPZGmshReader.h"
#include "pzgeopoint.h"
#include "tpzgeoelrefpattern.h"


//#include "tpzgeoelrefpattern.h.h"
using namespace std;
using namespace cv;

TPZGeoMesh *CreateGeoMesh(int order, int nx, int ny);
TPZCompMesh *CreateCompMesh(TPZGeoMesh *gmesh, int order);
TPZCompMesh *CreateCompMeshNewMat(TPZGeoMesh *gmesh, int order);
TPZCompMesh *CMeshMultphysics(TPZGeoMesh * gmesh, TPZVec<TPZCompMesh *> meshvec);
TPZCompMesh *CMeshMultphysicsNewMat(TPZGeoMesh * gmesh, TPZVec<TPZCompMesh *> meshvec);

TPZCompMesh *CreatePressureCmesh(TPZGeoMesh *gmesh, int order);
TPZCompMesh *CreateFluxCmesh(TPZGeoMesh *gmesh, int order);
TPZGeoMesh *gmesh1D(TPZVec<double> x1,TPZVec<double> x2, int n);
void H1TestMat();
void MixedTest();
void MixedTestNewMat();
int main(){
    
    MixedTest();
    return 0;
}

void H1TestMat(){
    TPZGeoMesh *gmesh = CreateGeoMesh(1, 10,10);
    std::ofstream out("testprint.vtk");
    TPZVTKGeoMesh::PrintGMeshVTK(gmesh, out);
    
    //   auto cmesh = CreateCompMesh(gmesh,1);
    auto cmesh = CreateCompMeshNewMat(gmesh, 2);
    std::ofstream outcmesh("testprintcmesh.vtk");
    
    TPZVTKGeoMesh::PrintCMeshVTK(cmesh, outcmesh);
    
    TPZAnalysis *an = new TPZAnalysis(cmesh,true);
    
#ifdef USING_MKL2
    TPZSymetricSpStructMatrix struct_mat(cmesh);
    struct_mat.SetNumThreads(4);
    an->SetStructuralMatrix(struct_mat);
#else
    TPZSkylineStructMatrix struct_mat(cmesh);
    struct_mat.SetNumThreads(4);
    an->SetStructuralMatrix(struct_mat);
#endif
    TPZStepSolver<STATE> step;
    step.SetDirect(ECholesky);
    an->SetSolver(step);
    
    
    // Solving the LS
    an->Assemble();
    an->Solve();
    
    // post-processing step
    {
        
        const int dim = an->Mesh()->Dimension();
        int div = 0;
        std::string plotfile = "h1_approximation.vtk";
        TPZStack<std::string> scalar_names, vec_names;
        scalar_names.push_back("Solution");
        vec_names.push_back("Derivative");
        an->DefineGraphMesh(dim,scalar_names,vec_names,plotfile);
        an->PostProcess(div,dim);
        std::cout << "Standard post-processing finished." << std::endl;
        
    }
}
void MixedTest(){
    
 //   TPZGeoMesh *gmesh = CreateGeoMesh(1,2,1);
    TPZVec<REAL> x1(3,0.0);
    TPZVec<REAL> x2(3,0.0);
    x2[0]=1.0;
  
    TPZGeoMesh *gmesh = gmesh1D(x1, x2, 3);
    int flux_order = 1;
    int p_order = 1;
    
    {
#ifdef PZDEBUG
        std::ofstream file("maze.txt");
        gmesh->Print(file);
        
        std::ofstream out("maze.vtk");
        TPZVTKGeoMesh::PrintGMeshVTK(gmesh, out, true);
#endif
    }
    
    TPZCompMesh *cmesh_flux =CreateFluxCmesh(gmesh,flux_order);
    TPZCompMesh *cmesh_presure =CreatePressureCmesh(gmesh,p_order);
    
    std::ofstream outflux_antes("Flux_antes.txt");
    cmesh_flux->Print(outflux_antes);
    std::ofstream outpressure_antes("Pressure_antes.txt");
    cmesh_presure->Print(outpressure_antes);
    
    TPZVec<TPZCompMesh *> fmeshvec(2);
    fmeshvec[0]=cmesh_flux;
    fmeshvec[1]=cmesh_presure;
    gmesh->ResetReference();
    
    TPZCompMesh *MixedMesh = CMeshMultphysics(gmesh,fmeshvec);
    std::ofstream file("MixedCMesh.vtk");
    TPZVTKGeoMesh::PrintCMeshVTK(MixedMesh, file);
    
    std::ofstream out("MixedCMesh.txt");
    MixedMesh->Print(out);
    
    
    
    std::cout << "number of equations = " << MixedMesh->NEquations() << std::endl;
    //Solving the system:
    bool optimizeBandwidth = true;
    MixedMesh->InitializeBlock();
    
    TPZCompMesh * cmesh_m_Hybrid;
//    TPZManVector<TPZCompMesh*, 3> meshvector_Hybrid(3);
//    TPZHybridizeHDiv hybridizer;
//    cmesh_m_Hybrid = hybridizer.Hybridize(MixedMesh, false, 1.);
//    cmesh_m_Hybrid->InitializeBlock();
   
    cmesh_m_Hybrid =MixedMesh;
    meshvector_Hybrid = fmeshvec;
    
    std::ofstream outhy("MixedCMesh_Hy.txt");
    cmesh_m_Hybrid->Print(outhy);
    
  
    
    bool must_opt_band_width_Q = true;
    int number_threads = 4;
    TPZAnalysis *an = new TPZAnalysis(cmesh_m_Hybrid,must_opt_band_width_Q);
    TPZElementMatrix mat;
    TPZElementMatrix ef;
    
    MixedMesh->Element(0)->CalcStiff(mat, ef);
    mat.Print(std::cout);
    
    //
    TPZSymetricSpStructMatrix sparse_matrix(cmesh_m_Hybrid);
    TPZStepSolver<STATE> step;
    sparse_matrix.SetNumThreads(number_threads);
    step.SetDirect(ELDLt);
    an->SetStructuralMatrix(sparse_matrix);
    an->SetSolver(step);
    an->Assemble();
    an->Solve();
    std::ofstream out2("global.txt");
    
    std::cout<<"order: "<<an->Mesh()->GetDefaultOrder()<<std::endl;
//    virtual void Print(const char *name, std::ostream &out = std::cout ,const MatrixOutputFormat form = EFormatted) const;
    
    an->Solver().Matrix()->Print("EK= ",std::cout,EMathematicaInput);
    an->Rhs().Print("EF= ",std::cout,EMathematicaInput);
    an->Mesh()->Solution().Print("alphas= ",std::cout,EMathematicaInput);
    std::ofstream filemalla("malla_sol.txt");
    an->Mesh()->Print(filemalla);
    std::ofstream outflux_despues("Flux_despues.txt");
    meshvector_Hybrid[0]->Print(outflux_despues);
    std::ofstream outpressure_despues("Pressure_despues.txt");
    meshvector_Hybrid[1]->Print(outpressure_despues);
    
    
    TPZBuildMultiphysicsMesh::TransferFromMeshes(meshvector_Hybrid, MixedMesh);
    //POS
    TPZManVector<std::string,10> scalnames(3), vecnames(1);
    vecnames[0]  = "Flux";
    
    scalnames[0] = "Pressure";
    scalnames[1] = "Permeability";
    scalnames[2] = "Divergence";
    
    
    const int dim = an->Mesh()->Dimension();
    int div = 0;
      std::string plotfile = "hdiv_approximation.vtk";
    an->DefineGraphMesh(dim,scalnames,vecnames,plotfile);
    an->PostProcess(div,dim);
    std::cout << "Standard post-processing finished." << std::endl;
    
    return 0;
    
}
TPZGeoMesh *CreateGeoMesh(int order, int nx, int ny){
    TPZGeoMesh *gmesh = new TPZGeoMesh;
    // Creating the Geo mesh
    TPZManVector<REAL,3> x0(3,0.),x1(3,0.0);
    x1[0] = 1;
    x1[1] = 0;
    TPZManVector<int,2> nel(2,2);
    nel[0] = nx;
    nel[1] = ny;
    TPZGenGrid gengrid(nel,x0,x1);
    gengrid.SetElementType(EOned);
   
    gmesh->SetDimension(1);
    gengrid.Read(gmesh);
    //gengrid.Read(gmesh,2);
    
  //  gengrid.SetBC(TPZGeoMesh *gr, int side, int bc)
    gengrid.SetBC(gmesh, 0, -1);
    gengrid.SetBC(gmesh, 1, -2);
//    gengrid.SetBC(gmesh, 6, -3);
//    gengrid.SetBC(gmesh, 7, -4);
    
    gmesh->BuildConnectivity();
    return gmesh;
   
}
TPZCompMesh *CreateCompMesh(TPZGeoMesh *gmesh, int order){
    
    int impervious_mat = 1;
    int permeable_mat = 2;
    int dim = gmesh->Dimension();
    
    REAL perm_0 = 1.0e-3;
    REAL perm_1 = 1000.0;
    
    REAL conv=0;
    TPZVec<REAL> convdir(dim,0.0);
    
    TPZCompMesh *cmesh = new TPZCompMesh(gmesh);
    cmesh->SetDimModel(dim);
    
    TPZMatPoisson3d *mat_0 = new TPZMatPoisson3d(impervious_mat,dim);
    TPZMatPoisson3d *mat_1 = new TPZMatPoisson3d(permeable_mat,dim);
    
    
//    TPZDarcyFlow *mat_0 = new TPZDarcyFlow(impervious_mat,dim);
//    TPZDarcyFlow *mat_1 = new TPZDarcyFlow(permeable_mat,dim);
//
    
    mat_0->SetParameters(perm_0, conv, convdir);
    mat_1->SetParameters(perm_1, conv, convdir);
    
    
    
    //  inserting volumetric materials objects
    cmesh->InsertMaterialObject(mat_0);
  //  cmesh->InsertMaterialObject(mat_1);
    
    
    int type_D = 0;
    int type_N = 1;
    TPZFMatrix<STATE> val1(1, 1, 0.), val2(1, 1, 0.);
    
    // Insert boundary conditions
    //Neumann boundary conditions (flux = 0)
    int right_bc_id = -2;
    val2(0,0) = 0.0;
    TPZMaterial * right_bc = mat_0->CreateBC(mat_0, right_bc_id, type_D, val1, val2);
    cmesh->InsertMaterialObject(right_bc);
    
    int left_bc_id = -4;
    val2(0,0) = 100.0;
    TPZMaterial * left_bc = mat_0->CreateBC(mat_0, left_bc_id, type_D, val1, val2);
    cmesh->InsertMaterialObject(left_bc);
    
    int bottom_bc_1id = -1;
    val2(0,0) = 0.0;
    TPZMaterial * bottom_bc_1 = mat_0->CreateBC(mat_0, bottom_bc_1id, type_N, val1, val2);
    cmesh->InsertMaterialObject(bottom_bc_1);
    
    int top_bc_1id = -3;
    val2(0,0) = 0.0;
    TPZMaterial * top_bc_1 = mat_0->CreateBC(mat_0, top_bc_1id, type_N, val1, val2);
    cmesh->InsertMaterialObject(top_bc_1);
    
    
    cmesh->SetName("LaberintoTest");
    cmesh->SetAllCreateFunctionsContinuous();
    cmesh->SetDefaultOrder(order);
    cmesh->AutoBuild();
    cmesh->Print();
    
    std::ofstream file("cmesh_h.txt");
    cmesh->Print(file);
    
    
    return cmesh;
    
    
}
TPZCompMesh *CreateCompMeshNewMat(TPZGeoMesh *gmesh, int order){
    
    int mat_id = 1;
    
    int dim = gmesh->Dimension();
    
    TPZCompMesh *cmesh = new TPZCompMesh(gmesh);
    cmesh->SetDimModel(dim);
    
    
    
    TPZMixedDarcyFlow *mat_0 = new TPZMixedDarcyFlow(mat_id,2);
    
    cmesh->InsertMaterialObject(mat_0);
    
    int type_D = 0;
    int type_N = 1;
    TPZFMatrix<STATE> val1(1, 1, 0.), val2(1, 1, 0.);
    
    int right_bc_id = -1;
    val2(0,0) = 0.0;
    TPZMaterial * right_bc = mat_0->CreateBC(mat_0, right_bc_id, type_D, val1, val2);
    cmesh->InsertMaterialObject(right_bc);
    
     right_bc_id = -2;
    val2(0,0) = 0.0;
    TPZMaterial * right_bc1 = mat_0->CreateBC(mat_0, right_bc_id, type_D, val1, val2);
    cmesh->InsertMaterialObject(right_bc1);
    
     right_bc_id = -3;
    val2(0,0) = 0.0;
    TPZMaterial * right_bc2 = mat_0->CreateBC(mat_0, right_bc_id, type_D, val1, val2);
    cmesh->InsertMaterialObject(right_bc2);
     right_bc_id = -4;
    val2(0,0) = 0.0;
    TPZMaterial * right_bc3 = mat_0->CreateBC(mat_0, right_bc_id, type_D, val1, val2);
    cmesh->InsertMaterialObject(right_bc3);
    
    
    cmesh->SetName("LaberintoTest");
    cmesh->SetDefaultOrder(order);
    cmesh->AutoBuild();
    cmesh->SetAllCreateFunctionsContinuous();
    
    
    std::cout<<"N elements: "<<cmesh->NElements()<<std::endl;
    std::cout<<"N elements: "<<gmesh->NElements()<<std::endl;
    
    //  cmesh->InsertMaterialObject(mat_1);
    return cmesh;
}
void sourceterm( const TPZVec<REAL> &pt, TPZVec<STATE> &disp){
    
    STATE x = pt[0];
    STATE y = pt[1];
    STATE z = pt[2];
    double fx= 0.0*2*M_PI*M_PI*sin(M_PI*x)*sin(M_PI*y);
    
    disp[0]=fx;
}

TPZCompMesh *CMeshMultphysics(TPZGeoMesh * gmesh, TPZVec<TPZCompMesh *> meshvec){
    
    //Creating computational mesh for multiphysic elements
    TPZCompMesh *mphysics = new TPZCompMesh(gmesh);
    
    int impervious_mat = 1;
  
    int dim = gmesh->Dimension();
    
    REAL perm_0 = 1.0;
 
    
    REAL conv=0;
    TPZVec<REAL> convdir(dim,0.0);
    
    std::cout<<mphysics->NMaterials();
    
    TPZMixedPoisson *mat_0 = new TPZMixedPoisson(impervious_mat,dim);
    mat_0->SetPermeability(perm_0);
    
   // mat->SetPermeability(1.e-3);
    TPZFNMatrix<9,REAL> K(3,3,0.),KInv(3,3,0.);
    K(0,0) = 1.;
    K(1,1) = 1.;
    K(2,2) = 1.;
    KInv(0,0) = 1.;
    KInv(1,1) = 1;
    KInv(2,2) = 1;
    mat_0->SetPermeabilityTensor(K, KInv);
   
    

    
    mat_0->SetParameters(perm_0, conv, convdir);
   
    mphysics->InsertMaterialObject(mat_0);
    
    //Inserir condicoes de contorno
    
    
    int type_D = 0;
    int type_N = 1;
    TPZFMatrix<STATE> val1(1, 1, 0.), val2(1, 1, 0.);
    
    // Insert boundary conditions
    //Neumann boundary conditions (flux = 0)
    
    int bottom_bc_1id = -1;
    val2(0,0) = 10;
    TPZMaterial * bottom_bc_1 = mat_0->CreateBC(mat_0, bottom_bc_1id, type_D, val1, val2);
    mphysics->InsertMaterialObject(bottom_bc_1);
    
    int right_bc_id = -2;
    val2(0,0) = 100.0;
    TPZMaterial * right_bc = mat_0->CreateBC(mat_0, right_bc_id, type_D, val1, val2);
    mphysics->InsertMaterialObject(right_bc);
    

//    int top_bc_1id = -3;
//    val2(0,0) = 0.0;
//    TPZMaterial * top_bc_1 = mat_0->CreateBC(mat_0, top_bc_1id, type_N, val1, val2);
//    mphysics->InsertMaterialObject(top_bc_1);
//
//
//    int left_bc_id = -4;
//    val2(0,0) = 10.0;
//    TPZMaterial * left_bc = mat_0->CreateBC(mat_0, left_bc_id, type_D, val1, val2);
//    mphysics->InsertMaterialObject(left_bc);
    

    mphysics->SetAllCreateFunctionsMultiphysicElem();
    mphysics->SetDimModel(gmesh->Dimension());
    mphysics->AutoBuild();
    
    TPZBuildMultiphysicsMesh::AddElements(meshvec, mphysics);
    TPZBuildMultiphysicsMesh::AddConnects(meshvec,mphysics);
    TPZBuildMultiphysicsMesh::TransferFromMeshes(meshvec, mphysics);
    
#ifdef PZDEBUG
    std::ofstream file("cmesh_mphysics.txt");
    mphysics->Print(file);
#endif
    
    return mphysics;
}
TPZCompMesh *CreatePressureCmesh(TPZGeoMesh *gmesh, int pOrder){
    int MatId = 1;
    int dim = gmesh->Dimension();
    REAL perm_0 = 1.0;
   
    REAL conv=0;
    TPZVec<REAL> convdir(dim,0.0);
    
    TPZCompMesh *cmesh = new TPZCompMesh(gmesh);
    cmesh->SetDimModel(dim);
    
    cmesh->SetDefaultOrder(pOrder); //Default polynomial order of the approximation
    cmesh->SetAllCreateFunctionsDiscontinuous();
    cmesh->ApproxSpace().CreateDisconnectedElements(true);
    
        TPZMatPoisson3d *mat_0 = new TPZMatPoisson3d(MatId,dim);
        mat_0->SetParameters(perm_0, conv, convdir);
        
        
        //  inserting volumetric materials objects
        cmesh->InsertMaterialObject(mat_0);
        
    
    
  
    
    cmesh->SetName("Pressure");
    cmesh->AutoBuild();
    
    int ncon = cmesh->NConnects();
    for(int i=0; i<ncon; i++)
    {
        TPZConnect &newnod = cmesh->ConnectVec()[i];
        newnod.SetLagrangeMultiplier(1);
    }
    
#ifdef PZDEBUG
    std::ofstream out("cmeshPress.txt");
    cmesh->Print(out);
#endif
    
    return cmesh;
}
TPZCompMesh *CreateFluxCmesh(TPZGeoMesh *gmesh, int pOrder){
    int Mat_Id = 1;
 
    int dim = gmesh->Dimension();
    
    REAL perm_0 = 1.0;
 
    
    REAL conv=0;
    TPZVec<REAL> convdir(dim,0.0);
    
    TPZCompMesh *cmesh = new TPZCompMesh(gmesh);
    cmesh->SetDimModel(dim);
    cmesh->SetDefaultOrder(pOrder); //Default polynomial order of the approximation
    
    //Definition of the approximation space:
    
    
        TPZMatPoisson3d *mat_0 = new TPZMatPoisson3d(Mat_Id,dim);
        mat_0->SetParameters(perm_0, conv, convdir);
        
        
        //  inserting volumetric materials objects
        cmesh->InsertMaterialObject(mat_0);
        
        
        cmesh->SetAllCreateFunctionsHDiv(); //Creating H(div) functions
        
        
        int type_D = 0;
        int type_N = 1;
        TPZFMatrix<STATE> val1(1, 1, 0.), val2(1, 1, 0.);
        
        // Insert boundary conditions
        //Neumann boundary conditions (flux = 0)
        int right_bc_id = -2;
    val2(0,0)=100.0;
        TPZMaterial * right_bc = mat_0->CreateBC(mat_0, right_bc_id, type_D, val1, val2);
        cmesh->InsertMaterialObject(right_bc);
        
//        int left_bc_id = -4;
//        TPZMaterial * left_bc = mat_0->CreateBC(mat_0, left_bc_id, type_D, val1, val2);
//        cmesh->InsertMaterialObject(left_bc);
    
        int bottom_bc_1id = -1;
        val2(0,0)=10.0;
        TPZMaterial * bottom_bc_1 = mat_0->CreateBC(mat_0, bottom_bc_1id, type_D, val1, val2);
        cmesh->InsertMaterialObject(bottom_bc_1);
        
//        int top_bc_1id = -3;
//        TPZMaterial * top_bc_1 = mat_0->CreateBC(mat_0, top_bc_1id, type_D, val1, val2);
//        cmesh->InsertMaterialObject(top_bc_1);
    

    cmesh->SetName("LaberintoTest");
    cmesh->AutoBuild();
    
#ifdef PZDEBUG
    std::ofstream file("cmesh_flux.txt");
    cmesh->Print(file);
#endif
    
    return cmesh;
}

TPZCompMesh *CMeshMultphysicsNewMat(TPZGeoMesh * gmesh, TPZVec<TPZCompMesh *> meshvec){
    
    //Creating computational mesh for multiphysic elements
    TPZCompMesh *mphysics = new TPZCompMesh(gmesh);
    
    int impervious_mat = 1;
    
    int dim = gmesh->Dimension();
    
    REAL perm_0 = 1.0;
    
    
    REAL conv=0;
    TPZVec<REAL> convdir(dim,0.0);
    
    std::cout<<mphysics->NMaterials();
    
    TPZMixedDarcyFlow *mat_0 = new TPZMixedDarcyFlow(impervious_mat,dim);
    mat_0->SetPermeability(perm_0);
    
    
    TPZAutoPointer<TPZFunction<STATE> > sourcetermm = new TPZDummyFunction<STATE>(sourceterm, 5);
    
    mat_0->SetForcingFunction(sourcetermm);
    
    
 //   mat_0->SetParameters(perm_0, conv, convdir);
    
    mphysics->InsertMaterialObject(mat_0);
    
    //Inserir condicoes de contorno
    
    
    int type_D = 0;
    int type_N = 1;
    TPZFMatrix<STATE> val1(1, 1, 0.), val2(1, 1, 0.);
    
    // Insert boundary conditions
    //Neumann boundary conditions (flux = 0)
    
    int bottom_bc_1id = -1;
    val2(0,0) = 0;
    TPZMaterial * bottom_bc_1 = mat_0->CreateBC(mat_0, bottom_bc_1id, type_D, val1, val2);
    mphysics->InsertMaterialObject(bottom_bc_1);
    
    int right_bc_id = -2;
    val2(0,0) = 0.0;
    TPZMaterial * right_bc = mat_0->CreateBC(mat_0, right_bc_id, type_D, val1, val2);
    mphysics->InsertMaterialObject(right_bc);
    
    
    int top_bc_1id = -3;
    val2(0,0) = 0.0;
    TPZMaterial * top_bc_1 = mat_0->CreateBC(mat_0, top_bc_1id, type_D, val1, val2);
    mphysics->InsertMaterialObject(top_bc_1);
    
    
    int left_bc_id = -4;
    val2(0,0) = 0.0;
    TPZMaterial * left_bc = mat_0->CreateBC(mat_0, left_bc_id, type_D, val1, val2);
    mphysics->InsertMaterialObject(left_bc);
    
    
    mphysics->SetAllCreateFunctionsMultiphysicElem();
    mphysics->SetDimModel(gmesh->Dimension());
    mphysics->AutoBuild();
    
    TPZBuildMultiphysicsMesh::AddElements(meshvec, mphysics);
    TPZBuildMultiphysicsMesh::AddConnects(meshvec,mphysics);
    TPZBuildMultiphysicsMesh::TransferFromMeshes(meshvec, mphysics);
    
#ifdef PZDEBUG
    std::ofstream file("cmesh_mphysics.txt");
    //mphysics->Print(file);
#endif
    
    return mphysics;
}
TPZGeoMesh *gmesh1D(TPZVec<REAL> x1,TPZVec<REAL> x2, int n){
    TPZGeoMesh *gmesh = new TPZGeoMesh;
//    gmesh->NodeVec().Resize(nNods+4);
//    gmesh->NodeVec()[nNods]=Node[0];
    
    TPZGeoNode Node;
    Node.SetCoord(x1);
    gmesh->NodeVec().Resize(n+1);
    gmesh->NodeVec()[0]=Node;
    for (int i=0; i<=n ; i++) {
        TPZGeoNode Node;
        TPZVec<REAL> xi(3,0.0);
        xi[0] = x1[0] + ((x2[0]-x1[0])/n)*i;
        Node.SetCoord(xi);
        gmesh->NodeVec()[i]=Node;
    }
    
    for (int i=0; i<n; i++) {
        int64_t elnum=i;
        TPZVec<int64_t> cords(2);
        cords[0]=i;
        cords[1]=i+1;
        gmesh->CreateGeoElement(EOned, cords, 1, elnum);
    }
    
    gmesh->BuildConnectivity();
    TPZGeoElSide gelside(gmesh->Element(0),0);
    TPZGeoElBC gbc(gelside, -1);
    int nel = gmesh->NElements();
    TPZGeoElSide gelside2(gmesh->Element(nel-2),1);
    TPZGeoElBC gbc2(gelside2, -2);
    gmesh->BuildConnectivity();
    std::ofstream file("meshoned.txt");
    gmesh->Print(file);
    gmesh->SetDimension(1);
    return gmesh;
}
TPZCompMesh *CmeshPressure(TPZGeoMesh *gmesh){
    
}
