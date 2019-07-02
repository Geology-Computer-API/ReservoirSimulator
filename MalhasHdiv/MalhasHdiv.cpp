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
#include "TPZHybridizeHDiv.h"
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
#include "TPZNullMaterial.h"
#include "mixedpoisson.h"
#include "pzbndcond.h"
#include "TPZSSpStructMatrix.h"
#include "pzskylstrmatrix.h"
#include "TPZSkylineNSymStructMatrix.h"

#include "pzanalysis.h"
#include "pzfstrmatrix.h"
#include "pzskylstrmatrix.h"
#include "pzstepsolver.h"
#include "pzsolve.h"
#include "TPZPersistenceManager.h"
#include "pzbuildmultiphysicsmesh.h"
#include "TPZMultiphysicsCompMesh.h"
#include "TPZCompMeshTools.h"

#ifdef _AUTODIFF
#include "tfad.h"
#include "fad.h"
#include "pzextractval.h"
#endif


//Creating geometric mesh
TPZGeoMesh * GenerateGmesh(int nx, int ny, double l, double h);

//Stablish de form fuction
void Ladoderecho (const TPZVec<REAL> &pt, TPZVec<STATE> &disp);

//Stablish an exact solution
void SolExact(const TPZVec<REAL> &ptx, TPZVec<STATE> &sol, TPZFMatrix<STATE> &flux);

//Creating computational pressure mesh
TPZCompMesh * GeneratePressureCmesh(TPZGeoMesh *Gmesh, int order_internal);

//Creating computational constant pressure mesh
TPZCompMesh * GenerateConstantCmesh(TPZGeoMesh *Gmesh, bool third_LM);

//Creating computational flux mesh
TPZCompMesh * GenerateFluxCmesh(TPZGeoMesh *Gmesh, int order_internal, int order_border);

//Creating computational mixed mesh
TPZMultiphysicsCompMesh * GenerateMixedCmesh(TPZVec<TPZCompMesh *> fvecmesh, int order);

//Test for the HDiv case
void HDivTest(int nx, int ny, int order_small, int order_high);

using namespace std;

int main(){
    
#ifdef LOG4CXX
    InitializePZLOG();
#endif
    HDivTest(2, 2, 1, 2);
}

/**
 * @brief Generates a Hdiv test
 * @param nx: number of partions on x
 * @param ny: number of partions on y
 * @param order_small: for the low order mesh
 * @param order_high: for the high order mesh
 * @return Two computational meshes with different orders
 */
void HDivTest(int nx, int ny, int order_small, int order_high){
    
    //Generating low order mesh
    TPZGeoMesh *gmesh = GenerateGmesh(nx, ny, 1, 1);
    TPZMultiphysicsCompMesh *MixedMesh_coarse = 0;
    TPZManVector<TPZCompMesh *> vecmesh(4); //vecmesh: Stands for vector coarse mesh
    {
        TPZCompMesh  *pmesh = GeneratePressureCmesh(gmesh, order_high);
        TPZCompMesh *qmesh = GenerateFluxCmesh(gmesh, order_high, order_small);
        TPZCompMesh *p0mesh = GenerateConstantCmesh(gmesh,false);
        TPZCompMesh *distmesh = GenerateConstantCmesh(gmesh,true);
        vecmesh[0] = qmesh;
        vecmesh[1] = pmesh;
        vecmesh[2] = p0mesh; //Mean pressure
        vecmesh[3] = distmesh;  //Distributed flux
        MixedMesh_coarse = GenerateMixedCmesh(vecmesh, 1);
    }
    // Created condensed elements for the elements that have internal nodes
    bool KeepOneLagrangian = true;
    bool KeepMatrix = false;
    TPZCompMeshTools::CreatedCondensedElements(MixedMesh_coarse, KeepOneLagrangian, KeepMatrix);

    TPZMultiphysicsCompMesh *MixedMesh_fine = 0;
    TPZManVector<TPZCompMesh *> vefmesh(4); //vefmesh: Stands for vector fine mesh
    {
        //Generating second mesh
        TPZCompMesh *qmesh_2 = GenerateFluxCmesh(gmesh, order_high, order_high);
        TPZCompMesh  *pmesh_2 = GeneratePressureCmesh(gmesh, order_high);
        TPZCompMesh *p0mesh = GenerateConstantCmesh(gmesh,false);
        TPZCompMesh *distmesh = GenerateConstantCmesh(gmesh,true);
        vefmesh[0] = qmesh_2;
        vefmesh[1] = pmesh_2;
        vefmesh[2] = p0mesh;
        vefmesh[3] = distmesh;
        MixedMesh_fine = GenerateMixedCmesh(vefmesh, 2);
    }

    // Created condensed elements for the elements that have internal nodes
    TPZCompMeshTools::CreatedCondensedElements(MixedMesh_fine, KeepOneLagrangian, KeepMatrix);

    //Solving the system:
    MixedMesh_coarse->InitializeBlock();
    MixedMesh_fine->InitializeBlock();
    bool must_opt_band_width_Q = false;
    int number_threads = 4;
    
    //Analysis
    TPZAnalysis *an_coarse = new TPZAnalysis(MixedMesh_coarse,must_opt_band_width_Q);
    TPZAnalysis *an_fine = new TPZAnalysis(MixedMesh_fine,must_opt_band_width_Q);
    TPZSkylineStructMatrix sparse_matrix_coarse(MixedMesh_coarse);
    TPZSkylineStructMatrix sparse_matrix_fine(MixedMesh_fine);
    TPZStepSolver<STATE> step;
    sparse_matrix_coarse.SetNumThreads(number_threads);
    sparse_matrix_fine.SetNumThreads(number_threads);
    step.SetDirect(ELDLt);
    an_coarse->SetStructuralMatrix(sparse_matrix_coarse);
    an_fine->SetStructuralMatrix(sparse_matrix_fine);
    an_coarse->SetSolver(step);
    an_fine->SetSolver(step);
    an_coarse->Assemble();
    an_fine->Assemble();
    {
        TPZAnalysis anloc(MixedMesh_coarse,false);
        std::string filename("Shape.vtk");
        TPZVec<int64_t> indices(20);
        for(int i=0; i<20; i++) indices[i] = i;
        anloc.ShowShape(filename, indices,1,"Flux");
    }
//    an_coarse->Solver().Matrix()->Print("k = ",std::cout,EMathematicaInput);
//    an_coarse->Rhs().Print("f = ",std::cout,EMathematicaInput);
    an_coarse->Solve();
//    an_coarse->Solution().Print("x = ",std::cout,EMathematicaInput);
    an_fine->Solve();
    
    TPZBuildMultiphysicsMesh::TransferFromMultiPhysics(vecmesh, MixedMesh_coarse);
    
//    vecmesh[2]->Solution().Print("g = ",std::cout,EMathematicaInput);
    
    //PostProcess
    TPZStack<std::string> scalar, vectors;
    TPZManVector<std::string,10> scalnames(4), vecnames(1);
    vecnames[0]  = "Flux";
    scalnames[0] = "Pressure";
    scalnames[1] = "Permeability";
    scalnames[2] = "g_average";
    scalnames[3] = "u_average";
    
    std::ofstream filePrint_coarse("MixedHdiv_coarse.txt");
    MixedMesh_coarse->Print(filePrint_coarse);
    std::string name_coarse = "MixedHdiv_coarse.vtk";
    
    std::ofstream filePrint_fine("MixedHdiv_fine.txt");
    MixedMesh_fine->Print(filePrint_fine);
    std::string name_fine = "MixedHdiv_fine.vtk";
    
    an_coarse->DefineGraphMesh(2, scalnames, vecnames, name_coarse);
    an_coarse->PostProcess(0,2);
    
    an_fine->DefineGraphMesh(2, scalnames, vecnames, name_fine);
    an_fine->PostProcess(0,2);
    
    int64_t target_index = 1;
    
    TPZVec<int64_t> equ_indexes(1);
    equ_indexes[0] = target_index;
    std::string name_phi = "MixedHdiv2_shape.vtk";
    std::string scal_name("Pressure");
    std::string vec_name("Flux");
//    TPZBuildMultiphysicsMesh::ShowShape(vefmesh, MixedMesh_fine, *an_fine, name_phi, equ_indexes);
    
}

/**
 * @brief Generates the geometric mesh
 * @param nx: number of partions on x
 * @param ny: number of partions on y
 * @param l: lenght
 * @param h: height
 * @return Geometric mesh
 */
TPZGeoMesh * GenerateGmesh(int nx, int ny, double l, double h){
    
    TPZGeoMesh *gmesh = new TPZGeoMesh;
    
    TPZVec<int> nels(3,0);
    nels[0]=nx;
    nels[1]=ny;
    
    TPZVec<REAL> x0(3,0.0);
    TPZVec<REAL> x1(3,l);
    x1[1]=h;
    x1[2]=0;
    
    //Setting boundary conditions (negative numbers to recognize them)
    TPZGenGrid gen(nels,x0,x1);
    gen.SetElementType(EQuadrilateral);
    gen.Read(gmesh);
    gen.SetBC(gmesh, 4, -1);
    gen.SetBC(gmesh, 5, -2);
    gen.SetBC(gmesh, 6, -3);
    gen.SetBC(gmesh, 7, -4);
    return gmesh;
}

/**
 * @brief Generates the pressure computational mesh
 * @param Geometric mesh
 * @param Order
 * @return Pressure computational mesh
 */
TPZCompMesh * GeneratePressureCmesh(TPZGeoMesh *Gmesh, int order){
    
    TPZCompMesh *Cmesh= new TPZCompMesh (Gmesh);
    
    Cmesh->SetDimModel(Gmesh->Dimension());
    Cmesh->SetDefaultOrder(order);
    Cmesh->SetAllCreateFunctionsDiscontinuous();
    Cmesh->ApproxSpace().CreateDisconnectedElements(true);
    
    //Add material to the mesh
    int dimen = Gmesh->Dimension();
    int MaterialId = 1;
    STATE Permeability=1;
    
    TPZMatPoisson3d *mat = new TPZMatPoisson3d(MaterialId, dimen);
    
    //No convection
    REAL conv=0;
    TPZVec<REAL> convdir(dimen, 0);
    mat->SetParameters(Permeability, conv, convdir);
    
    //Insert material to mesh
    Cmesh->InsertMaterialObject(mat);
    
    //Autobuild
    Cmesh->AutoBuild();
    
    int ncon = Cmesh->NConnects();
    for(int i=0; i<ncon; i++)
    {
        TPZConnect &newnod = Cmesh->ConnectVec()[i];
        newnod.SetLagrangeMultiplier(1);
    }
    return Cmesh;
}

//Creating pressure mesh
TPZCompMesh * GenerateConstantCmesh(TPZGeoMesh *Gmesh, bool third_LM)
{
    TPZCompMesh *Cmesh= new TPZCompMesh (Gmesh);
    
    Cmesh->SetDimModel(Gmesh->Dimension());
    Cmesh->SetDefaultOrder(0);
    Cmesh->SetAllCreateFunctionsDiscontinuous();
    Cmesh->ApproxSpace().CreateDisconnectedElements(true);
    
    //Add material to the mesh
    int dimen = Gmesh->Dimension();
    int MaterialId = 1;
    
    TPZNullMaterial *mat =new TPZNullMaterial(MaterialId);
    mat->SetDimension(dimen);
    mat->SetNStateVariables(1);

    //Insert material to mesh
    Cmesh->InsertMaterialObject(mat);
    
    //Autobuild
    Cmesh->AutoBuild();
    
    int ncon = Cmesh->NConnects();
    for(int i=0; i<ncon; i++)
    {
        TPZConnect &newnod = Cmesh->ConnectVec()[i];
        if (third_LM) {
            newnod.SetLagrangeMultiplier(3);
        }else{
            newnod.SetLagrangeMultiplier(2);
        }
    }
    return Cmesh;
}


/**
 * @brief Generates the flux computational mesh
 * @param mesh: Geometric mesh
 * @param order_internal: Order used for internal elements
 * @param order_border: Order used for border elements
 * @return Flux computational mesh
 */
TPZCompMesh * GenerateFluxCmesh(TPZGeoMesh *mesh, int order_internal, int order_border){
    
    int dimen = mesh->Dimension();
    TPZCompMesh *Cmesh = new TPZCompMesh(mesh);
    Cmesh->SetDimModel(dimen);
    Cmesh->SetDefaultOrder(order_border);
    
    //Definition of the approximation space
    int perm=1;
    REAL conv=0;
    REAL perme=1;
    TPZVec<REAL> convdir(dimen , 0.0);
    TPZMatPoisson3d *mat = new TPZMatPoisson3d(perm , dimen);
    mat->SetParameters(perme, conv, convdir);
    
    //Inserting volumetric materials objects
    Cmesh->InsertMaterialObject(mat);
    
    //Create H(div) functions
    Cmesh->SetAllCreateFunctionsHDiv();
    
    //Insert boundary conditions
    int D=0;
    int BC1=-1;
    TPZFMatrix<STATE> val1(1,1,0.0);
    TPZFMatrix<STATE> val2(1,1,0.0);
    TPZMaterial *bc1 = mat->CreateBC(mat, BC1, D, val1, val2);
    Cmesh->InsertMaterialObject(bc1);
    
    int BC2=-2;
    val2(0,0)=0;
    TPZMaterial *bc2 = mat->CreateBC(mat, BC2, D, val1, val2);
    Cmesh->InsertMaterialObject(bc2);
    
    int BC3=-3;
    val2(0,0)=0;
    TPZMaterial *bc3 = mat->CreateBC(mat, BC3, D, val1, val2);
    Cmesh->InsertMaterialObject(bc3);
    
    int BC4=-4;
    val2(0,0)=0;
    TPZMaterial *bc4 = mat->CreateBC(mat, BC4, D, val1, val2);
    Cmesh->InsertMaterialObject(bc4);
    
    Cmesh->AutoBuild();
    
    int64_t nel = Cmesh->NElements();
    for (int el=0; el<nel; el++) {
        TPZCompEl *cel = Cmesh->Element(el);
        if(!cel) continue;
        TPZInterpolatedElement *intel = dynamic_cast<TPZInterpolatedElement *>(cel);
        if(!intel) DebugStop();
        TPZGeoEl *gel = intel->Reference();
        intel->SetSideOrder(gel->NSides()-1, order_internal);
    }
    Cmesh->ExpandSolution();
    return Cmesh;
}

/**
 * @brief Generates the mixed computational mesh
 * @param fvecmesh: Vector thats contains flux and pressure computational mesh
 * @param order: Order
 * @return Mixed computational mesh
 */

TPZMultiphysicsCompMesh * GenerateMixedCmesh(TPZVec<TPZCompMesh *> fvecmesh, int order){
    TPZGeoMesh *gmesh = fvecmesh[1]->Reference();
    TPZMultiphysicsCompMesh *MixedMesh = new TPZMultiphysicsCompMesh(gmesh);
    
    //Definition of the approximation space
    int dimen= gmesh->Dimension();
    int matnum=1;
    REAL perm=1;
    REAL conv=0;
    TPZVec<REAL> convdir(dimen , 0.0);
    
    //Inserting material
    TPZMixedPoisson *mat = new TPZMixedPoisson(matnum, dimen);
    
    mat->SetPermeability(perm);
    mat->SetParameters(perm, conv, convdir);
    
    TPZAutoPointer<TPZFunction<STATE> > sourceterm = new TPZDummyFunction<STATE>(Ladoderecho, 5);
    mat->SetForcingFunction(sourceterm);
    
    //Inserting volumetric materials objects
    MixedMesh->InsertMaterialObject(mat);
    
    //Boundary conditions
    int D=0;
    int BC1=-1;
    TPZFMatrix<STATE> val1(1,1,0.0);
    TPZFMatrix<STATE> val2(1,1,0.0);
    TPZMaterial *bc1 = mat->CreateBC(mat, BC1, D, val1, val2);
    MixedMesh->InsertMaterialObject(bc1);
    
    int BC2=-2;
    val2(0,0)=0;
    TPZMaterial *bc2 = mat->CreateBC(mat, BC2, D, val1, val2);
    MixedMesh->InsertMaterialObject(bc2);
    
    int BC3=-3;
    val2(0,0)=0;
    TPZMaterial *bc3 = mat->CreateBC(mat, BC3, D, val1, val2);
    MixedMesh->InsertMaterialObject(bc3);
    
    int BC4=-4;
    val2(0,0)=0;
    TPZMaterial *bc4 = mat->CreateBC(mat, BC4, D, val1, val2);
    MixedMesh->InsertMaterialObject(bc4);
    
    MixedMesh->SetAllCreateFunctionsMultiphysicElem();
    MixedMesh->SetDimModel(dimen);
    
    //Autobuild
    TPZManVector<int,5> active_approx_spaces(4); /// 1 stands for an active approximation spaces
    active_approx_spaces[0] = 1;
    active_approx_spaces[1] = 1;
    active_approx_spaces[2] = 1;
    active_approx_spaces[3] = 1;
    MixedMesh->BuildMultiphysicsSpace(active_approx_spaces,fvecmesh);
    
    TPZBuildMultiphysicsMesh::AddElements(fvecmesh, MixedMesh);
    TPZBuildMultiphysicsMesh::AddConnects(fvecmesh,MixedMesh);
    TPZBuildMultiphysicsMesh::TransferFromMeshes(fvecmesh, MixedMesh);
    
    std::cout<<"n equ: "<<MixedMesh->NEquations()<<std::endl;
    std::cout<<"n equ: "<<fvecmesh[0]->NEquations()<<std::endl;
    std::cout<<"n equ: "<<fvecmesh[1]->NEquations()<<std::endl;
    std::cout<<"n equ: "<<fvecmesh[2]->NEquations()<<std::endl;
    std::cout<<"n equ: "<<fvecmesh[3]->NEquations()<<std::endl;
    std::cout<<"------------------------------"<<std::endl;
    return MixedMesh;
};

/**
 * @brief Generates the force function
 * @param pt: Points values
 * @param disp:
 * @return Force function value
 */
void Ladoderecho (const TPZVec<REAL> &pt, TPZVec<STATE> &disp){
    
    STATE x = pt[0];
    STATE y = pt[1];
    
    //Force function definition
    double fx= -2*M_PI*M_PI*sin(M_PI*x)*sin(M_PI*y);
    
    disp[0]=fx;
      disp[0]=0.0;
}
