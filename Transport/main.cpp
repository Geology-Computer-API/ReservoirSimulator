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
#include "pzpoisson3d.h"
#include "mixedpoisson.h"
#include "pzl2projection.h"
#include "TPZPrimalPoisson.h"
#include "pzmatmixedpoisson3d.h"

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
#include "pzbuildmultiphysicsmesh.h"
#include <opencv2/opencv.hpp>

#include "TRSLinearInterpolator.h"
#include "TPZMatLaplacian.h"
#include "pzpoisson3d.h"
#include "pzbndcond.h"
#include "TPZSSpStructMatrix.h"
#include "pzskylstrmatrix.h"
#include "TPZSkylineNSymStructMatrix.h"
#include "TRSRibFrac.h"
#include "TRSRibs.h"
#include "TRSFace.h"
#include "TPZMixedDarcyFlow.h"
#include "TPZMultiphysicsCompMesh.h"
TPZGeoMesh *GenerateGmesh(int nx,int ny,double l,double h);
TPZCompMesh * FluxMesh(TPZGeoMesh * geometry);
TPZCompMesh * PressureMesh(TPZGeoMesh * geometry);
TPZCompMesh * SaturationWMesh(TPZGeoMesh * geometry);
TPZMultiphysicsCompMesh * MPCMeshMixed(TPZGeoMesh * geometry, TPZVec<TPZCompMesh *> &meshvec);
void forcing(const TPZVec<REAL> &p, TPZVec<STATE> &f);
TPZAnalysis * CreateAnalysis(TPZCompMesh * cmesh);
void CreateInterfaces(TPZMultiphysicsCompMesh *MPCompMesh);

using namespace std;
int main(){
    
    TPZGeoMesh *gmesh = GenerateGmesh(2, 1, 2, 2);
    TPZVec<TPZCompMesh *> meshvec;
    TPZMultiphysicsCompMesh *cmixedmesh = NULL;
    cmixedmesh = MPCMeshMixed(gmesh, meshvec);
    
    CreateInterfaces(cmixedmesh);
    std::ofstream filemixed("mixedMesh.txt");
    cmixedmesh->Print(filemixed);
    
    TPZCompMesh *cmeshm =NULL;

    cmeshm=cmixedmesh;
  
    TPZAnalysis *an = CreateAnalysis(cmeshm);
    
    std::cout << "Assembly neq = " << cmeshm->NEquations() << std::endl;
    an->Assemble();
    
    std::cout << "Solution of the system" << std::endl;
    an->Solve();
     TPZBuildMultiphysicsMesh::TransferFromMultiPhysics(meshvec, cmeshm);
    TPZStack<std::string,10> scalnames, vecnames;
    vecnames.Push("q");
    scalnames.Push("p");
//    scalnames.Push("Permeability");
    
    int div = 0;
    std::string fileresult("case_1_1k.vtk");
    an->DefineGraphMesh(2,scalnames,vecnames,fileresult);
    an->PostProcess(div,2);
    
    
}
TPZGeoMesh *GenerateGmesh(int nx,int ny,double l,double h){

    // Creating the Geo mesh
    TPZManVector<REAL,3> x0(3,0.),x1(3,0.0);
    x1[0] = l;
    x1[1] = h;
    TPZManVector<int,2> nelx(2,0);
    nelx[0] = nx;
    nelx[1] = ny;
    TPZGenGrid gengrid(nelx,x0,x1);
    gengrid.SetElementType(EQuadrilateral);
    TPZGeoMesh *gmesh = new TPZGeoMesh;
    gmesh->SetDimension(2);
    gengrid.Read(gmesh);
    //gengrid.Read(gmesh,2);
    
    //MatsID
    int nels = gmesh->NElements();
    
    
    //gengrid.SetBC(TPZGeoMesh *gr, int side, int bc)
    gengrid.SetBC(gmesh, 4, -1);
    gengrid.SetBC(gmesh, 5, -2);
    gengrid.SetBC(gmesh, 6, -3);
    gengrid.SetBC(gmesh, 7, -4);
    
    gmesh->BuildConnectivity();
    return gmesh;
}
TPZCompMesh *GenerateFluxMesh(TPZGeoMesh *gmesh){
    TPZCompMesh *cmesh = new TPZCompMesh(gmesh);
  
    
    int impervious_mat=1;
    int dim = 2;
    //Definition of the approximation space:
    
    
    
        TPZMatMixedPoisson3D * mat_0 = new TPZMatMixedPoisson3D(1,dim);
        mat_0->SetPermeability(1.0);
        cmesh->InsertMaterialObject(mat_0);
    
   


   
   
    int type_D = 0;
    int type_N = 1;
    TPZFMatrix<STATE> val1(1, 1, 0.), val2(1, 1, 0.);
    
    // Insert boundary conditions
    //Neumann boundary conditions (flux = 0)
    int right_bc_id = -2;
    val2(0,0)=2.0;
    TPZMaterial * right_bc = mat_0->CreateBC(mat_0, right_bc_id, type_D, val1, val2);
    cmesh->InsertMaterialObject(right_bc);
    
    int left_bc_id = -4;
    val2(0,0)=1.0;
    TPZMaterial * left_bc = mat_0->CreateBC(mat_0, left_bc_id, type_D, val1, val2);
    cmesh->InsertMaterialObject(left_bc);
    
    int bottom_bc_1id = -1;
     val2(0,0)=0.0;
    TPZMaterial * bottom_bc_1 = mat_0->CreateBC(mat_0, bottom_bc_1id, type_N, val1, val2);
    cmesh->InsertMaterialObject(bottom_bc_1);
    
    int top_bc_1id = -3;
    val2(0,0)=0.0;
    TPZMaterial * top_bc_1 = mat_0->CreateBC(mat_0, top_bc_1id, type_N, val1, val2);
    cmesh->InsertMaterialObject(top_bc_1);
    
    
    cmesh->SetDimModel(2);
    cmesh->SetDefaultOrder(1);
    cmesh->SetAllCreateFunctionsHDiv();
    cmesh->AutoBuild();
    cmesh->InitializeBlock();
    
    std::ofstream file("test_flux.vtk");
    TPZVTKGeoMesh::PrintCMeshVTK(cmesh, file);
    return cmesh;
}
TPZCompMesh * PressureMesh(TPZGeoMesh * geometry){
    
    int dimension = geometry->Dimension();
    TPZCompMesh *cmesh = new TPZCompMesh(geometry);
    
    TPZFMatrix<STATE> val1(dimension,dimension,0.0),val2(dimension,1,0.0);
    

        TPZMatMixedPoisson3D * volume = new TPZMatMixedPoisson3D(1,dimension);
        volume->SetPermeability(1.0);
        cmesh->InsertMaterialObject(volume);
  
    
    cmesh->SetDimModel(dimension);
    cmesh->SetDefaultOrder(1);
    cmesh->SetAllCreateFunctionsContinuous();
    cmesh->ApproxSpace().CreateDisconnectedElements(true);
    
    cmesh->AutoBuild();
    cmesh->InitializeBlock();
    
    int ncon = cmesh->NConnects();
    for(int i=0; i<ncon; i++)
    {
        TPZConnect &newnod = cmesh->ConnectVec()[i];
        newnod.SetLagrangeMultiplier(1);
    }

    
    return cmesh;
    
}

TPZCompMesh * FluxMesh(TPZGeoMesh * geometry){
    
    int dimension = geometry->Dimension();

    
    
    TPZCompMesh *cmesh = new TPZCompMesh(geometry);
    
    TPZFMatrix<STATE> val1(dimension,dimension,0.0),val2(dimension,1,0.0);
    

        TPZMatMixedPoisson3D * volume = new TPZMatMixedPoisson3D(1,dimension);
        volume->SetPermeability(1.0);
        
        TPZDummyFunction<STATE> * rhs_exact = new TPZDummyFunction<STATE>(forcing, 5);
        cmesh->InsertMaterialObject(volume);
        
    
    TPZMaterial * face = volume->CreateBC(volume,-1,1,val1,val2);
    cmesh->InsertMaterialObject(face);
    
    TPZMaterial * face_2 = volume->CreateBC(volume,-2,0,val1,val2);
    cmesh->InsertMaterialObject(face_2);
    
    TPZMaterial * face_3 = volume->CreateBC(volume,-3,1,val1,val2);
    cmesh->InsertMaterialObject(face_3);
    
    TPZMaterial * face_4 = volume->CreateBC(volume,-4,0,val1,val2);
    cmesh->InsertMaterialObject(face_4);
    
    cmesh->SetDimModel(dimension);
    cmesh->SetDefaultOrder(1);
    cmesh->SetAllCreateFunctionsHDiv();
    cmesh->AutoBuild();
    cmesh->InitializeBlock();
    

    
    return cmesh;
    
}
void forcing(const TPZVec<REAL> &p, TPZVec<STATE> &f){
    REAL x = p[0];
    REAL y = p[0];
    REAL z = p[0];
    f[0]=0.0*x*y*z;
}
TPZMultiphysicsCompMesh * MPCMeshMixed(TPZGeoMesh * geometry, TPZVec<TPZCompMesh *> &meshvec){
    
    int dimension = geometry->Dimension();
  
    TPZMultiphysicsCompMesh *cmesh = new TPZMultiphysicsCompMesh(geometry);
    
    TPZFNMatrix<9,STATE> val1(dimension,dimension,0.0),val2(dimension,1,0.0);
    
  
        
    TPZMixedDarcyFlow * volume = new TPZMixedDarcyFlow(1,dimension);
    volume->SetPermeability(1.0);

    cmesh->InsertMaterialObject(volume);
    
    val2(0,0)=0.0;
    TPZMaterial * face = volume->CreateBC(volume,-1,1,val1,val2);
    cmesh->InsertMaterialObject(face);
    
    val2(0,0)=2.0;
    TPZMaterial * face_2 = volume->CreateBC(volume,-2,0,val1,val2);
    cmesh->InsertMaterialObject(face_2);
    
    val2(0,0)=0.0;
    TPZMaterial * face_3 = volume->CreateBC(volume,-3,1,val1,val2);
    cmesh->InsertMaterialObject(face_3);
    
    val2(0,0)=1.0;
    TPZMaterial * face_4 = volume->CreateBC(volume,-4,0,val1,val2);
    cmesh->InsertMaterialObject(face_4);
    
    
    cmesh->SetDimModel(dimension);
    
    TPZManVector<TPZCompMesh * ,3> mesh_vec(3);
    mesh_vec[0] = FluxMesh(geometry);
    mesh_vec[1] = PressureMesh(geometry);
    mesh_vec[2] = SaturationWMesh(geometry);
    TPZManVector<int,5> active_approx_spaces(3); /// 1 stands for an active approximation spaces
    active_approx_spaces[0] = 1;
    active_approx_spaces[1] = 1;
    active_approx_spaces[2] = 1;
    cmesh->BuildMultiphysicsSpace(active_approx_spaces,mesh_vec);
    
    
    meshvec = mesh_vec;
    return cmesh;
}

TPZAnalysis * CreateAnalysis(TPZCompMesh * cmesh){
    
    TPZAnalysis * analysis = new TPZAnalysis(cmesh, true);
    
 
        
        TPZSymetricSpStructMatrix matrix(cmesh);
        //        TPZSkylineStructMatrix matrix(cmesh);
        matrix.SetNumThreads(1);
        analysis->SetStructuralMatrix(matrix);
        TPZStepSolver<STATE> step;
        step.SetDirect(ELDLt);
        analysis->SetSolver(step);
        
        return analysis;
   
  
    
}
TPZCompMesh * SaturationWMesh(TPZGeoMesh * geometry){
    
    int dim = geometry->Dimension();
    TPZFMatrix<STATE> val1(3,2,0.), val2(3,1,0.);
    
    // Malha computacional
    TPZCompMesh *cmesh = new TPZCompMesh(geometry);
 
        // Material medio poroso
    TPZMatPoisson3d * mat = new TPZMatPoisson3d(1,dim);
    cmesh->InsertMaterialObject(mat);
        
    // Void material
    int matIdL2Proj = 2;
    TPZVec<STATE> sol(1,0.);
    TPZL2Projection *matl2proj = new TPZL2Projection(matIdL2Proj,dim,1,sol);
    cmesh->InsertMaterialObject(matl2proj);
        

    val2(0,0)=0.0;
    TPZMaterial * face = mat->CreateBC(mat,-1,1,val1,val2);
    cmesh->InsertMaterialObject(face);
    
    val2(0,0)=2.0;
    TPZMaterial * face_2 = mat->CreateBC(mat,-2,0,val1,val2);
    cmesh->InsertMaterialObject(face_2);
    
    val2(0,0)=0.0;
    TPZMaterial * face_3 = mat->CreateBC(mat,-3,1,val1,val2);
    cmesh->InsertMaterialObject(face_3);
    
    val2(0,0)=1.0;
    TPZMaterial * face_4 = mat->CreateBC(mat,-4,0,val1,val2);
    cmesh->InsertMaterialObject(face_4);
    
    
    
    // Setando L2
    cmesh->SetDimModel(dim);
    cmesh->SetDefaultOrder(0);
    
    cmesh->SetAllCreateFunctionsDiscontinuous();
    cmesh->AutoBuild();
    
    cmesh->AdjustBoundaryElements();
    cmesh->CleanUpUnconnectedNodes();
    
    
    return cmesh;
}
void CreateInterfaces(TPZMultiphysicsCompMesh *MPCompMesh){
    
   // TPZGeoMesh *fgmesh = MPCompMesh->ResetReference();
    MPCompMesh->LoadReferences();
    
    
    std::ofstream file1("flux_antes.vtk");
    TPZVTKGeoMesh::PrintCMeshVTK(MPCompMesh, file1);
    
    MPCompMesh->LoadReferences();
    
    // Creation of interface elements
    int nel = MPCompMesh->ElementVec().NElements();
    for(int el = 0; el < nel; el++)
    {
        TPZCompEl * compEl = MPCompMesh->ElementVec()[el];
        if(!compEl) continue;
        TPZGeoEl * gel = compEl->Reference();
        if(!gel) {continue;}
        if(gel->HasSubElement()) {continue;}
        int index = compEl ->Index();
        if(compEl->Dimension() == MPCompMesh->Dimension())
        {
            TPZMultiphysicsElement * InterpEl = dynamic_cast<TPZMultiphysicsElement *>(MPCompMesh->ElementVec()[index]);
            if(!InterpEl) continue;
            InterpEl->CreateInterfaces();
        }
    }
    std::ofstream file2("flux_despues.vtk");
    TPZVTKGeoMesh::PrintCMeshVTK(MPCompMesh, file2);
}
