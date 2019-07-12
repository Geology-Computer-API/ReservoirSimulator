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
#include "TPZMixedDarcyWithFourSpaces.h"
#include "pzcmesh.h"


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

//Creates index vector
void IndexVectorCoFi(TPZMultiphysicsCompMesh *Coarse_sol, TPZMultiphysicsCompMesh *Fine_sol, TPZVec<int64_t> & indexvec);

//Test for the HDiv 1D and 2D cases
void HDivTest(int nx, int ny, int order_small, int order_high, bool condense_equations_Q, int dimen);

//Test for the HDiv case One Dimension
void HDivTestOne(int nx, int order_small, int order_high, bool condense_equations_Q);
TPZGeoMesh * GenerateGmeshOne(int nx, double l);

//Transfer DOF from coarse mesh to fine mesh
void TransferDegreeOfFreedom(TPZFMatrix<STATE> & CoarseDoF, TPZFMatrix<STATE> & FineDoF, TPZVec<int64_t> & DoFIndexes);

void ConfigurateAnalyses(TPZCompMesh * cmesh_c, TPZCompMesh * cmesh_f, bool must_opt_band_width_Q, int number_threads, TPZAnalysis *an_c,TPZAnalysis *an_f, bool UsePardiso_Q);


using namespace std;

int main(){
    
#ifdef LOG4CXX
    InitializePZLOG();
#endif
    HDivTest(4, 4, 1, 4, false, 3);
//    HDivTestOne(50, 1, 2, true);
}

/**
 * @brief Generates a Hdiv test
 * @param nx: number of partions on x
 * @param ny: number of partions on y
 * @param order_small: for the low order mesh
 * @param order_high: for the high order mesh
 * @param condense_equations_Q: true or false wether you want or not condensation
 * @param dimen: 1 or 2 dimension problem
 * @return Two computational meshes with different orders
 */
void HDivTest(int nx, int ny, int order_small, int order_high, bool condense_equations_Q, int dimen){
    
    if (dimen!=1 || dimen !=2) {
        DebugStop();
    }
    
    // Created condensed elements for the elements that have internal nodes
    bool KeepOneLagrangian = false;
    bool KeepMatrix = false;
    bool render_shapes_Q = false;
    
    
    //Generating low order mesh
    
    TPZGeoMesh *gmesh = GenerateGmesh(nx, ny, 1, 1);
    TPZMultiphysicsCompMesh *MixedMesh_coarse = 0;
    TPZManVector<TPZCompMesh *> vecmesh_c(4); //vecmesh: Stands for vector coarse mesh
    {
        TPZCompMesh *q_cmesh = GenerateFluxCmesh(gmesh, order_high, order_small);
        TPZCompMesh *p_cmesh = GeneratePressureCmesh(gmesh, order_high);
        TPZCompMesh *gavg_cmesh = GenerateConstantCmesh(gmesh,false);
        TPZCompMesh *pavg_cmesh = GenerateConstantCmesh(gmesh,true);
        vecmesh_c[0] = q_cmesh;   //Flux
        vecmesh_c[1] = p_cmesh;   //Pressure
        vecmesh_c[2] = gavg_cmesh;    //Average distribute flux
        vecmesh_c[3] = pavg_cmesh;    //Average pressure
        
        MixedMesh_coarse = GenerateMixedCmesh(vecmesh_c, 1);
    }

    if (condense_equations_Q) {
        MixedMesh_coarse->ComputeNodElCon();
        int dim = MixedMesh_coarse->Dimension();
        int64_t nel = MixedMesh_coarse->NElements();
        for (int64_t el =0; el<nel; el++) {
            TPZCompEl *cel = MixedMesh_coarse->Element(el);
            if(!cel) continue;
            TPZGeoEl *gel = cel->Reference();
            if(!gel) continue;
            if(gel->Dimension() != dim) continue;
            int nc = cel->NConnects();
            cel->Connect(nc-1).IncrementElConnected();
        }
        
        TPZCompMeshTools::CreatedCondensedElements(MixedMesh_coarse, KeepOneLagrangian, KeepMatrix);
    }
    
    TPZMultiphysicsCompMesh * MixedMesh_fine = 0;
    TPZManVector<TPZCompMesh *> vecmesh_f(4); //vefmesh: Stands for vector fine mesh
    {
        
        TPZCompMesh *q_cmesh = GenerateFluxCmesh(gmesh, order_high, order_high);
        TPZCompMesh *p_cmesh = GeneratePressureCmesh(gmesh, order_high);
        TPZCompMesh *gavg_cmesh = GenerateConstantCmesh(gmesh,false);
        TPZCompMesh *pavg_cmesh = GenerateConstantCmesh(gmesh,true);
        vecmesh_f[0] = q_cmesh;   //Flux
        vecmesh_f[1] = p_cmesh;   //Pressure
        vecmesh_f[2] = gavg_cmesh;    //Average distribute flux
        vecmesh_f[3] = pavg_cmesh;    //Average pressure

        MixedMesh_fine = GenerateMixedCmesh(vecmesh_f, 2);
    }

    
    if (condense_equations_Q) {
        
        MixedMesh_fine->ComputeNodElCon();
        int dim = MixedMesh_fine->Dimension();
        int64_t nel = MixedMesh_fine->NElements();
        for (int64_t el =0; el<nel; el++) {
            TPZCompEl *cel = MixedMesh_fine->Element(el);
            if(!cel) continue;
            TPZGeoEl *gel = cel->Reference();
            if(!gel) continue;
            if(gel->Dimension() != dim) continue;
            int nc = cel->NConnects();
            cel->Connect(nc-1).IncrementElConnected();
        }

        // Created condensed elements for the elements that have internal nodes
        TPZCompMeshTools::CreatedCondensedElements(MixedMesh_fine, KeepOneLagrangian, KeepMatrix);
    }

    //Solving the system:
    MixedMesh_coarse->InitializeBlock();
    MixedMesh_fine->InitializeBlock();
    bool must_opt_band_width_Q = true;
    int number_threads = 0;
    TPZAnalysis *an_c = new TPZAnalysis;
    TPZAnalysis *an_f = new TPZAnalysis;
    ConfigurateAnalyses(MixedMesh_coarse, MixedMesh_fine, must_opt_band_width_Q, number_threads, an_c, an_f, true);
    
    if(render_shapes_Q){
        TPZAnalysis anloc(MixedMesh_coarse,false);
        std::string filename("Shape.vtk");
        TPZVec<int64_t> indices(20);
        for(int i=0; i<20; i++) indices[i] = i;
        anloc.ShowShape(filename, indices,1,"Flux");
    }
    
    // Solving and postprocessing problems separately
    if(1){
        
        an_c->Assemble();
        an_f->Assemble();
        an_c->Rhs() *= -1.0;
        an_f->Rhs() *= -1.0;
        
        an_c->Solve();
        an_f->Solve();
        
        TPZBuildMultiphysicsMesh::TransferFromMultiPhysics(vecmesh_c, MixedMesh_coarse);
        
        //PostProcess
        TPZStack<std::string> scalar, vectors;
        TPZManVector<std::string,10> scalnames(4), vecnames(1);
        vecnames[0]  = "q";
        scalnames[0] = "p";
        scalnames[1] = "kappa";
        scalnames[1] = "div_q";
        scalnames[2] = "g_average";
        scalnames[3] = "u_average";
        
        std::ofstream filePrint_coarse("MixedHdiv_coarse.txt");
        MixedMesh_coarse->Print(filePrint_coarse);
        std::string name_coarse = "MixedHdiv_coarse.vtk";
        
        
        std::ofstream filePrint_fine("MixedHdiv_fine.txt");
        MixedMesh_fine->Print(filePrint_fine);
        std::string name_fine = "MixedHdiv_fine.vtk";
        
        an_c->DefineGraphMesh(2, scalnames, vecnames, name_coarse);
        an_c->PostProcess(0,2);
        
        an_f->DefineGraphMesh(2, scalnames, vecnames, name_fine);
        an_f->PostProcess(0,2);
    }
    
    /// An iterative solution
    {
        // Resolver coarse
        an_c->Assemble();
        
        // Resolver coarse
        an_c->Solve();
        an_c->Solution().Print("xc = ",std::cout,EMathematicaInput);
        
//        // Transfer the solution
//        TPZVec<int64_t> indexvec;
//        IndexVectorCoFi(MixedMesh_coarse, MixedMesh_fine,indexvec);
//        std::cout<<"n equ: "<<indexvec<<std::endl;
//        TransferDegreeOfFreedom(MixedMesh_coarse->Solution(), MixedMesh_fine->Solution(), indexvec);
        
        
        
        
        // Residual
        an_f->Assemble();
        an_f->Rhs().Print("r = ",std::cout,EMathematicaInput);
        an_f->Solve();
        an_f->Solution().Print("xf = ",std::cout,EMathematicaInput);
        

    }
    
//
//    int64_t target_index = 1;
//
//    TPZVec<int64_t> equ_indexes(1);
//    equ_indexes[0] = target_index;
//    std::string name_phi = "MixedHdiv2_shape.vtk";
//    std::string scal_name("Pressure");
//    std::string vec_name("Flux");
////    TPZBuildMultiphysicsMesh::ShowShape(vefmesh, MixedMesh_fine, *an_fine, name_phi, equ_indexes);
    

}

void ConfigurateAnalyses(TPZCompMesh * cmesh_c, TPZCompMesh * cmesh_f, bool must_opt_band_width_Q, int number_threads, TPZAnalysis *an_c,TPZAnalysis *an_f, bool UsePardiso_Q){
    
    an_c->SetCompMesh(cmesh_c,must_opt_band_width_Q);
    an_f->SetCompMesh(cmesh_f,must_opt_band_width_Q);
    TPZStepSolver<STATE> step;
    if (UsePardiso_Q) {
        
        TPZSymetricSpStructMatrix sparse_matrix_coarse(cmesh_c);
        TPZSymetricSpStructMatrix sparse_matrix_fine(cmesh_f);
        sparse_matrix_coarse.SetNumThreads(number_threads);
        sparse_matrix_fine.SetNumThreads(number_threads);
        an_c->SetStructuralMatrix(sparse_matrix_coarse);
        an_f->SetStructuralMatrix(sparse_matrix_fine);
        
    }else{
        
        TPZSkylineStructMatrix sparse_matrix_coarse(cmesh_c);           //Problems with the
        TPZSkylineStructMatrix sparse_matrix_fine(cmesh_f);
        sparse_matrix_coarse.SetNumThreads(number_threads);
        sparse_matrix_fine.SetNumThreads(number_threads);
        an_c->SetStructuralMatrix(sparse_matrix_coarse);
        an_f->SetStructuralMatrix(sparse_matrix_fine);

    }
    
    step.SetDirect(ELDLt);
    an_c->SetSolver(step);
    an_f->SetSolver(step);
    


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
    
    //Inserting material
    TPZMixedDarcyWithFourSpaces * mat = new TPZMixedDarcyWithFourSpaces(matnum, dimen);
    mat->SetPermeability(perm);
    
    TPZAutoPointer<TPZFunction<STATE> > sourceterm = new TPZDummyFunction<STATE>(Ladoderecho, 10);
    mat->SetForcingFunction(sourceterm);
    
    //Inserting volumetric materials objects
    MixedMesh->InsertMaterialObject(mat);
    
    //Boundary conditions
    int D=0;
    int BC1=-1;
    TPZFMatrix<STATE> val1(1,1,0.0);
    TPZFMatrix<STATE> val2(1,1,0.0);
    
    val2(0,0)=0.0;
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
 * @param disp: Vector
 */
void Ladoderecho (const TPZVec<REAL> &pt, TPZVec<STATE> &disp){
    
    STATE x = pt[0];
//    STATE y = pt[1];
    
    //Force function definition
    double fx= 4*M_PI*M_PI*sin(2*M_PI*x);
//    double fx= 4*M_PI*M_PI*sin(2*M_PI*x)*sin(2*M_PI*y);
    
    disp[0]=fx;
//      disp[0]=0.0;
}

/**
 * @brief Generates a index vector which relates coarse and fine mesh
 * @param Coarse_sol: Coarse mesh
 * @param Fine_sol: Fine mesh
 * @param indexvec:
 */
void IndexVectorCoFi(TPZMultiphysicsCompMesh *Coarse_sol, TPZMultiphysicsCompMesh *Fine_sol, TPZVec<int64_t> & indexvec)
{
    int64_t maxcone = Coarse_sol->NConnects();
    
    int64_t indexvecsize = 0;
    for (int j=0; j<maxcone; j++) {
        bool is_condensed = Coarse_sol->ConnectVec()[j].IsCondensed();
        if (is_condensed == true ) continue;

        int64_t sequence_coarse = Coarse_sol->ConnectVec()[j].SequenceNumber();
        int blocksize_coarse = Coarse_sol->Block().Size(sequence_coarse);

        indexvecsize += blocksize_coarse;
    }
    indexvec.Resize(indexvecsize,-1);
    
    for (int j=0; j<maxcone; j++) {
        bool is_condensed = Coarse_sol->ConnectVec()[j].IsCondensed();
        if (is_condensed == true ) continue;

        int64_t sequence_coarse = Coarse_sol->ConnectVec()[j].SequenceNumber();
        int64_t sequence_fine = Fine_sol->ConnectVec()[j].SequenceNumber();
        int blocksize_coarse = Coarse_sol->Block().Size(sequence_coarse);
        int blocksize_fine = Fine_sol->Block().Size(sequence_fine);

        if (blocksize_coarse > blocksize_fine) {
            DebugStop();

        }
    
        for(int i=0; i<blocksize_coarse; i++){
            int64_t pos_coarse = Coarse_sol->Block().Position(sequence_coarse);
            int64_t pos_fine = Fine_sol->Block().Position(sequence_fine);
            indexvec[pos_coarse+i] = pos_fine+i;
        }
    }
}

/**
 * @brief Transfer DOF from coarse to fine mesh
 * @param CoarseDoF: Solution coarse matrix
 * @param FineDoF: Solution fine matrix
 * @param DoFIndexes: DOF index vector
 */
void TransferDegreeOfFreedom(TPZFMatrix<STATE> & CoarseDoF, TPZFMatrix<STATE> & FineDoF, TPZVec<int64_t> & DoFIndexes){
    
    int64_t n_data = DoFIndexes.size();
    for (int64_t i = 0 ; i < n_data; i++) {
        FineDoF(DoFIndexes[i],0) = CoarseDoF(i,0);
    }

}


void HDivTestOne(int nx, int order_small, int order_high, bool condense_equations_Q){
    
    // Created condensed elements for the elements that have internal nodes
    bool KeepOneLagrangian = false;
    bool KeepMatrix = false;
    bool render_shapes_Q = false;
    
    
    //Generating low order mesh
    TPZGeoMesh *gmesh = GenerateGmeshOne(nx, 1);
    TPZMultiphysicsCompMesh *MixedMesh_coarse = 0;
    TPZManVector<TPZCompMesh *> vecmesh_c(4); //vecmesh: Stands for vector coarse mesh
    {
        TPZCompMesh *q_cmesh = GenerateFluxCmesh(gmesh, order_high, order_small);
        TPZCompMesh *p_cmesh = GeneratePressureCmesh(gmesh, order_high);
        TPZCompMesh *gavg_cmesh = GenerateConstantCmesh(gmesh,true);
        TPZCompMesh *pavg_cmesh = GenerateConstantCmesh(gmesh,false);
        vecmesh_c[0] = q_cmesh;   //Flux
        vecmesh_c[1] = p_cmesh;   //Pressure
        vecmesh_c[2] = gavg_cmesh;    //Average distribute flux
        vecmesh_c[3] = pavg_cmesh;    //Average pressure
        
        MixedMesh_coarse = GenerateMixedCmesh(vecmesh_c, 1);
    }
    
    if (condense_equations_Q) {
        MixedMesh_coarse->ComputeNodElCon();
        int dim = MixedMesh_coarse->Dimension();
        int64_t nel = MixedMesh_coarse->NElements();
        for (int64_t el =0; el<nel; el++) {
            TPZCompEl *cel = MixedMesh_coarse->Element(el);
            if(!cel) continue;
            TPZGeoEl *gel = cel->Reference();
            if(!gel) continue;
            if(gel->Dimension() != dim) continue;
            int nc = cel->NConnects();
            cel->Connect(nc-1).IncrementElConnected();
        }
        
        TPZCompMeshTools::CreatedCondensedElements(MixedMesh_coarse, KeepOneLagrangian, KeepMatrix);
    }
    
    TPZMultiphysicsCompMesh * MixedMesh_fine = 0;
    TPZManVector<TPZCompMesh *> vecmesh_f(4); //vefmesh: Stands for vector fine mesh
    {
        
        TPZCompMesh *q_cmesh = GenerateFluxCmesh(gmesh, order_high, order_high);
        TPZCompMesh *p_cmesh = GeneratePressureCmesh(gmesh, order_high);
        TPZCompMesh *gavg_cmesh = GenerateConstantCmesh(gmesh,true);
        TPZCompMesh *pavg_cmesh = GenerateConstantCmesh(gmesh,false);
        vecmesh_f[0] = q_cmesh;   //Flux
        vecmesh_f[1] = p_cmesh;   //Pressure
        vecmesh_f[2] = gavg_cmesh;    //Average distribute flux
        vecmesh_f[3] = pavg_cmesh;    //Average pressure
        
        MixedMesh_fine = GenerateMixedCmesh(vecmesh_f, 2);
    }
    
    
    if (condense_equations_Q) {
        
        MixedMesh_fine->ComputeNodElCon();
        int dim = MixedMesh_fine->Dimension();
        int64_t nel = MixedMesh_fine->NElements();
        for (int64_t el =0; el<nel; el++) {
            TPZCompEl *cel = MixedMesh_fine->Element(el);
            if(!cel) continue;
            TPZGeoEl *gel = cel->Reference();
            if(!gel) continue;
            if(gel->Dimension() != dim) continue;
            int nc = cel->NConnects();
            cel->Connect(nc-1).IncrementElConnected();
        }
        
        // Created condensed elements for the elements that have internal nodes
        TPZCompMeshTools::CreatedCondensedElements(MixedMesh_fine, KeepOneLagrangian, KeepMatrix);
    }
    
    //Solving the system:
    MixedMesh_coarse->InitializeBlock();
    MixedMesh_fine->InitializeBlock();
    bool must_opt_band_width_Q = true;
    int number_threads = 0;
    TPZAnalysis *an_c = new TPZAnalysis;
    TPZAnalysis *an_f = new TPZAnalysis;
    ConfigurateAnalyses(MixedMesh_coarse, MixedMesh_fine, must_opt_band_width_Q, number_threads, an_c, an_f,true);
    
    if(render_shapes_Q){
        TPZAnalysis anloc(MixedMesh_coarse,false);
        std::string filename("Shape.vtk");
        TPZVec<int64_t> indices(20);
        for(int i=0; i<20; i++) indices[i] = i;
        anloc.ShowShape(filename, indices,1,"Flux");
    }
    
    // Solving and postprocessing problems separately
    if(1){
        
        an_c->Assemble();
        an_f->Assemble();
        
        an_c->Rhs() *= -1.0;
        an_f->Rhs() *= -1.0;
        
        an_c->Solve();
        an_f->Solve();
        
        TPZBuildMultiphysicsMesh::TransferFromMultiPhysics(vecmesh_c, MixedMesh_coarse);
        
        //PostProcess
        TPZStack<std::string> scalar, vectors;
        TPZManVector<std::string,10> scalnames(4), vecnames(1);
        vecnames[0]  = "q";
        scalnames[0] = "p";
        scalnames[1] = "kappa";
        scalnames[1] = "div_q";
        scalnames[2] = "g_average";
        scalnames[3] = "u_average";
        
        std::ofstream filePrint_coarse("MixedHdiv_coarse.txt");
        MixedMesh_coarse->Print(filePrint_coarse);
        std::string name_coarse = "MixedHdiv_coarse.vtk";
        
        
        std::ofstream filePrint_fine("MixedHdiv_fine.txt");
        MixedMesh_fine->Print(filePrint_fine);
        std::string name_fine = "MixedHdiv_fine.vtk";
        
        an_c->DefineGraphMesh(1, scalnames, vecnames, name_coarse);
        an_c->PostProcess(0,1);
        
        an_f->DefineGraphMesh(1, scalnames, vecnames, name_fine);
        an_f->PostProcess(0,1);
        return ;

    }
    
    // An iterative solution
    {
        // Resolver coarse
        an_c->Assemble();
        
        // Resolver coarse
        an_c->Solve();
//        an_c->Solution().Print("xc = ",std::cout,EMathematicaInput);
        
        //        // Transfer the solution
        //        TPZVec<int64_t> indexvec;
        //        IndexVectorCoFi(MixedMesh_coarse, MixedMesh_fine,indexvec);
        //        std::cout<<"n equ: "<<indexvec<<std::endl;
        //        TransferDegreeOfFreedom(MixedMesh_coarse->Solution(), MixedMesh_fine->Solution(), indexvec);
        
        
        
        
        // Residual
        an_f->Assemble();
        an_f->Rhs().Print("r = ",std::cout,EMathematicaInput);
        
        an_f->Solve();
        an_f->Solution().Print("xf = ",std::cout,EMathematicaInput);
        
    }
    
    //
    //    int64_t target_index = 1;
    //
    //    TPZVec<int64_t> equ_indexes(1);
    //    equ_indexes[0] = target_index;
    //    std::string name_phi = "MixedHdiv2_shape.vtk";
    //    std::string scal_name("Pressure");
    //    std::string vec_name("Flux");
    ////    TPZBuildMultiphysicsMesh::ShowShape(vefmesh, MixedMesh_fine, *an_fine, name_phi, equ_indexes);
    
    
}

/**
 * @brief Generates a geometric 1D mesh
 * @param nx: number of partions on x
 * @param l: lenght
 * @return Geometric 1D mesh
 */
TPZGeoMesh * GenerateGmeshOne(int nx, double l){
    //Creates vector nodes
    double h = l/nx;
    int Domain_Mat_Id = 1;
    int Inlet_bc_Id = -1;
    int Outletbc_Id = -2;
    TPZVec<REAL> xp(3,0.0);
//    cout << xp;
    TPZGeoMesh *gmesh = new TPZGeoMesh;
    gmesh->NodeVec().Resize(nx+1);
//    gmesh->NodeVec().NElements();
    for (int64_t i=0; i<nx+1; i++) {
        xp[0] =(i)*h;
            cout << xp<<endl;
        gmesh->NodeVec()[i]= TPZGeoNode(i, xp, *gmesh);
    }

//    for (auto &item: gmesh->NodeVec()) {
//        item.Print();
//    }

    //Creates elements
    TPZVec<int64_t> cornerindexes(2);
    for (int64_t iel=0; iel<nx; iel++) {
        cornerindexes[0]=iel;
        cornerindexes[1]=iel+1;
        gmesh->CreateGeoElement(EOned, cornerindexes, Domain_Mat_Id, iel);
    }
    //Set BC
    gmesh->Element(0)->CreateBCGeoEl(0, Inlet_bc_Id);
    gmesh->Element(nx-1)->CreateBCGeoEl(1, Outletbc_Id);
    gmesh->SetDimension(1);
    gmesh->BuildConnectivity();
    
    return gmesh;
}
