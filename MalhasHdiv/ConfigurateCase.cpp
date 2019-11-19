//
//  ConfigurateCase.h
//  ReservoirSimulation
//
//  Created by Jorge Paúl Ordóñez Andrade on 10/11/19.

#include "ConfigurateCase.h"
#include "TPZRefPattern.h"
#include "TPZRefPatternTools.h"
#include "pzrefquad.h"
#include "TPZMHMixedMesh4SpacesControl.h"

void ComputeCoarseIndices(TPZGeoMesh *gmesh, TPZVec<int64_t> &coarseindices);
void DivideMesh(TPZGeoMesh * gmesh);

ConfigurateCase::ConfigurateCase(){
    
}
void PrintPoint(TPZVec<double> point){
    std::cout<<point[0]<<std::endl;
    std::cout<<point[1]<<std::endl;
    std::cout<<point[2]<<std::endl;
   
}
TPZCompMesh *ConfigurateCase::HDivMesh(TPZGeoMesh * gmesh, int orderfine, int ordercoarse){
    
    int dimension = gmesh->Dimension();         //Gets mesh dimension
    int nvols = fsim_case.omega_ids.size();     //
    int nbound = fsim_case.gamma_ids.size();    //
    
    TPZCompMesh *cmesh = new TPZCompMesh(gmesh);
    
    cmesh->SetDefaultOrder(ordercoarse);
    TPZFMatrix<STATE> val1(dimension,dimension,0.0),val2(dimension,1,0.0);
    
    for (int ivol=0; ivol<nvols; ivol++) {
        TPZMixedDarcyFlow * volume = new TPZMixedDarcyFlow(fsim_case.omega_ids[ivol],fsim_case.omega_dim[ivol]);
        volume->SetPermeability(fsim_case.permeabilities[ivol]);
        
        
        cmesh->InsertMaterialObject(volume);
        if (ivol==0) {
            for (int ibound=0; ibound<nbound; ibound++) {
                
                val2(0,0)=fsim_case.vals[ibound];
                int condType=fsim_case.type[ibound];
                TPZMaterial * face = volume->CreateBC(volume,fsim_case.gamma_ids[ibound],condType,val1,val2);
                cmesh->InsertMaterialObject(face);
            }
        }
    }
    
    
    cmesh->SetAllCreateFunctionsHDiv();
    cmesh->AutoBuild();
    
        int64_t nel = cmesh->NElements();
        for (int el=0; el<nel; el++) {
            TPZCompEl *cel = cmesh->Element(el);
            if(!cel) continue;
            TPZInterpolatedElement *intel = dynamic_cast<TPZInterpolatedElement *>(cel);
            if(!intel) DebugStop();
            TPZGeoEl *gel = intel->Reference();
            intel->SetSideOrder(gel->NSides()-1, ordercoarse);
        }
        cmesh->ExpandSolution();
    
    cmesh->InitializeBlock();
    
#ifdef PZDEBUG
    std::stringstream file_name;
    file_name << "q_cmesh_raw" << ".txt";
    std::ofstream sout(file_name.str().c_str());
    cmesh->Print(sout);
#endif
    
    return cmesh;
    
};
TPZCompMesh *ConfigurateCase::DiscontinuousMesh(TPZGeoMesh * gmesh, int order, int lagrangemult){
    int dimension = gmesh->Dimension();
    int nvols = fsim_case.omega_ids.size();
    TPZCompMesh *cmesh = new TPZCompMesh(gmesh);
    
    TPZFMatrix<STATE> val1(dimension,dimension,0.0),val2(dimension,1,0.0);
    
    cmesh->SetDimModel(dimension);
    cmesh->SetDefaultOrder(order);
    cmesh->SetAllCreateFunctionsDiscontinuous();
    
    if (lagrangemult ==1) {
        cmesh->ApproxSpace().CreateDisconnectedElements(true);
    }
    
     std::set<int> matids;
    if (order == 0) {
       
        for (int ivol=0; ivol < nvols; ivol++) {
            TPZNullMaterial * volume = new TPZNullMaterial(fsim_case.omega_ids[ivol]);
            cmesh->InsertMaterialObject(volume);
            if (fsim_case.omega_dim[ivol] == dimension) {
                matids.insert(fsim_case.omega_ids[ivol]);
            }
        }
        
    }
    if (order != 0) {
        
        for (int ivol=0; ivol < nvols; ivol++) {
            TPZMixedDarcyFlow * volume = new TPZMixedDarcyFlow(fsim_case.omega_ids[ivol],dimension);
            volume->SetPermeability(fsim_case.permeabilities[ivol]);
            cmesh->InsertMaterialObject(volume);
            if (fsim_case.omega_dim[ivol] == dimension) {
                matids.insert(fsim_case.omega_ids[ivol]);
            }
        }
    }
    
    
    cmesh->AutoBuild(matids);
    cmesh->InitializeBlock();
    
    int ncon = cmesh->NConnects();
    for(int i=0; i<ncon; i++)
    {
        TPZConnect &newnod = cmesh->ConnectVec()[i];
        newnod.SetLagrangeMultiplier(lagrangemult);
    }
    
#ifdef PZDEBUG
    std::stringstream file_name;
    file_name << "p_cmesh_raw" << ".txt";
    std::ofstream sout(file_name.str().c_str());
    cmesh->Print(sout);
#endif
    
    return cmesh;
};

TPZMultiphysicsCompMesh *ConfigurateCase::CreateMultCompMesh(){
    TPZVec<TPZCompMesh *> meshvec(4);

//    int orderfine = fsim_case.order_qfine;
//    int ordercoarse = fsim_case.order_qcoarse;
    
    int orderfine = m_fineorder;
    int ordercoarse = m_coarseorder;
    //
    
    int dimension = m_gmesh->Dimension();
    int nvols = fsim_case.omega_ids.size();
    int nbound= fsim_case.gamma_ids.size();
    if (nvols<1) {
        std::cout<<"Error: Omega is not defined."<<std::endl;
        DebugStop();
    }
    if (nbound<1) {
        std::cout<<"Error: Gamma is not defined."<<std::endl;
        DebugStop();
    }
    
    TPZMultiphysicsCompMesh *cmesh = new TPZMultiphysicsCompMesh(m_gmesh);
    TPZFNMatrix<9,STATE> val1(dimension,dimension,0.0),val2(dimension,1,0.0);
    
    for (int ivol=0; ivol<nvols; ivol++) {
        
        TPZMixedDarcyWithFourSpaces * volume = new TPZMixedDarcyWithFourSpaces(fsim_case.omega_ids[ivol],dimension);
        volume->SetPermeability(fsim_case.permeabilities[ivol]);
        cmesh->InsertMaterialObject(volume);

        
        if (ivol==0) {
            for (int ibound=0; ibound<nbound; ibound++) {
                val2(0,0)=fsim_case.vals[ibound];
                int condType=fsim_case.type[ibound];
                TPZMaterial * face = volume->CreateBC(volume,fsim_case.gamma_ids[ibound],condType,val1,val2);
                cmesh->InsertMaterialObject(face);
            }
        }
    }
    
    cmesh->SetDimModel(dimension);
    
    
    //
    meshvec[0] = HDivMesh(m_gmesh, orderfine, ordercoarse);
    meshvec[1] = DiscontinuousMesh(m_gmesh, ordercoarse, 1);
    meshvec[2] = DiscontinuousMesh(m_gmesh, 0, 2);
    meshvec[3] = DiscontinuousMesh(m_gmesh, 0, 3);
    TPZManVector<int,5> active_approx_spaces(4,1); 
   
    cmesh->BuildMultiphysicsSpace(active_approx_spaces,meshvec);
    if (fsim_case.IsCondensedQ) {             //Asks if you want to condesate the problem 
        cmesh->ComputeNodElCon();
        int dim = cmesh->Dimension();
        int64_t nel = cmesh->NElements();
        for (int64_t el =0; el<nel; el++) {
            TPZCompEl *cel = cmesh->Element(el);
            if(!cel) continue;
            TPZGeoEl *gel = cel->Reference();
            if(!gel) continue;
            if(gel->Dimension() != dim) continue;
            int nc = cel->NConnects();
            cel->Connect(nc-1).IncrementElConnected();
            
        }
        
        TPZCompMeshTools::CreatedCondensedElements(cmesh, fsim_case.KeepOneLagrangianQ, fsim_case.KeepMatrixQ);
        
        std::ofstream filePrint_coarse("MixedHdiv_coarse.txt");
        cmesh->Print(filePrint_coarse);
    }
    cmesh->InitializeBlock();
    return cmesh;
    
};

TPZAnalysis *ConfigurateCase::CreateAnalysis(TPZMultiphysicsCompMesh *mcmesh){
    
    TPZAnalysis *an = new TPZAnalysis(mcmesh);
    
    TPZStepSolver<STATE> step;
    if (fsim_case.UsePardisoQ) {
        
        TPZSymetricSpStructMatrix sparse_matrix(mcmesh);
        sparse_matrix.SetNumThreads(fsim_case.n_threads);
        an->SetStructuralMatrix(sparse_matrix);
        
    }else{
        
        TPZSkylineStructMatrix sparse_matrix(mcmesh);
        sparse_matrix.SetNumThreads(fsim_case.n_threads);
        an->SetStructuralMatrix(sparse_matrix);
        
    }
    step.SetDirect(ELDLt);
    an->SetSolver(step);
    
    return an;
}

TPZGeoMesh * ConfigurateCase::CreateUniformMesh(int nx, REAL L, int ny, REAL h, int nz, REAL w){
    
    TPZVec<int> nels(3,0);
    nels[0]=nx;         //Elements over x
    nels[1]=ny;         //Elements over y
    
    TPZVec<REAL> x0(3,0.0);
    TPZVec<REAL> x1(3,0.0);
    x1[0]=L;
    
    if (ny!=0) {
        x1[1]=h;
    }
    
    TPZGeoMesh *gmesh = new TPZGeoMesh;
    TPZGenGrid gen(nels,x0,x1);
    gen.SetRefpatternElements(true);

//    if (ny!=0) {
//        <#statements#>
//    }
//    gen.SetElementType(EQuadrilateral);
    gen.Read(gmesh);
    
    if (nz!=0 ) {
        double var = w/nz;
        TPZExtendGridDimension extend(gmesh,var);
        extend.SetElType(1);
        gmesh = extend.ExtendedMesh(nz);
    }
    
    if (nz!=0) {
        for (auto gel:gmesh->ElementVec()) {
            TPZFMatrix<REAL> coordinates;
            gel->NodesCoordinates(coordinates);
            if(coordinates(2,0)==0){
                gel->CreateBCGeoEl(20, -1);
            }
            if(coordinates(2,4)==w){
                gel->CreateBCGeoEl(25, -2);
            }
            
            if(coordinates(0,0)==0.0 ){
                gel->CreateBCGeoEl(24, -3);
            }
            if(coordinates(1,0)==0.0 ){
                gel->CreateBCGeoEl(21, -4);
            }
            
            if(coordinates(0,1)== L ){
                gel->CreateBCGeoEl(22, -5);
            }
            if(coordinates(1,3)==h){
                gel->CreateBCGeoEl(23, -6);
            }
        };
        gmesh->SetDimension(3);
    }
    
    if (ny!=0 && nz==0) {
//        gen.SetBC(gmesh, 4, -1);
//        gen.SetBC(gmesh, 5, -2);
//        gen.SetBC(gmesh, 6, -3);
//        gen.SetBC(gmesh, 7, -4);
        gmesh->SetDimension(2);
    }
  
    if (ny==0 && nz==0) {
        double dh = L/nx;
        int Domain_Mat_Id = 1;
        int Inlet_bc_Id = -1;
        int Outletbc_Id = -2;
        TPZVec<REAL> xp(3,0.0);   //Creates vector nodes
        
        gmesh->NodeVec().Resize(nx+1);
        for (int64_t i=0; i<nx+1; i++) {
            xp[0] =(i)*dh;
            gmesh->NodeVec()[i]= TPZGeoNode(i, xp, *gmesh);
        }
        
        TPZVec<int64_t> cornerindexes(2);   //Creates elements
        for (int64_t iel=0; iel<nx; iel++) {
            cornerindexes[0]=iel;
            cornerindexes[1]=iel+1;
            gmesh->CreateGeoElement(EOned, cornerindexes, Domain_Mat_Id, iel);
        }
        gmesh->Element(0)->CreateBCGeoEl(0, Inlet_bc_Id);     //Sets BC
        gmesh->Element(nx-1)->CreateBCGeoEl(1, Outletbc_Id);
        gmesh->SetDimension(1);
        gmesh->BuildConnectivity();
    }
    
    gmesh->BuildConnectivity();
    m_gmesh = gmesh;
    return gmesh;
    
}
TPZAutoPointer<TPZMHMixedMesh4SpacesControl> ConfigurateCase::CreateMHMMixedMesh(){
    
    TPZGeoMesh *gmeshcoarse = CreateGeowithRefPattern();
    
    
    int interface_mat_id = 600;
    int flux_order = 1;
    int p_order = 1;
    
    TPZAutoPointer<TPZMHMixedMesh4SpacesControl> MHMixed;
    
    {
        TPZAutoPointer<TPZGeoMesh> gmeshauto = new TPZGeoMesh(*gmeshcoarse);
        {
            std::ofstream out("gmeshauto.txt");
            gmeshauto->Print(out);
        }
        TPZMHMixedMesh4SpacesControl *mhm = new TPZMHMixedMesh4SpacesControl(gmeshauto);
        
        TPZVec<int64_t> coarseindices;
        ComputeCoarseIndices(gmeshauto.operator->(), coarseindices);
        
      
        
        
        mhm->DefinePartitionbyCoarseIndices(coarseindices);
        // criam-se apenas elementos geometricos
        std::ofstream out2("gmeshauto2.txt");
        mhm->Print(out2);
        //        MHMMixedPref << "MHMixed";
        MHMixed = mhm;
        
        TPZMHMixedMesh4SpacesControl &meshcontrol = *mhm;
        {
            std::set<int> matids;
            matids.insert(1);
            matids.insert(2);
            mhm->fMaterialIds = matids;
            matids.clear();
            matids.insert(-1);
            matids.insert(-2);
            matids.insert(-3);
            matids.insert(-4);
            matids.insert(-5);
            matids.insert(-6);
            mhm->fMaterialBCIds = matids;
        }
        
        
        InsertMaterialObjects(*mhm);
        
#ifdef PZDEBUG
        if(1)
        {
            std::ofstream out("MixedMeshControlHDiv.txt");
            meshcontrol.Print(out);
        }
#endif
        
        meshcontrol.SetInternalPOrder(1);
        meshcontrol.SetSkeletonPOrder(1);
        
        std::ofstream filemesh("despuesInterfaces.txt");
        meshcontrol.Print(filemesh);
        //        TPZVTKGeoMesh::PrintGMeshVTK(meshcontrol.GMesh(), filemesh);
        meshcontrol.DivideSkeletonElements(0);
        meshcontrol.DivideBoundarySkeletonElements();
        
        std::ofstream file_geo("geometry.txt");
        meshcontrol.CMesh()->Reference()->Print(file_geo);
        //
        bool substructure = false;
        //        std::ofstream filee("Submesh.txt");
        //        meshcontrol.CMesh()->Print(filee);
        std::map<int,std::pair<TPZGeoElSide,TPZGeoElSide>> test;
        
    
        
        meshcontrol.BuildComputationalMesh(substructure);
        
        
        
#ifdef PZDEBUG
        if(1)
        {
            std::ofstream file("GMeshControlHDiv.vtk");
            TPZVTKGeoMesh::PrintGMeshVTK(meshcontrol.GMesh().operator->(), file);
        }
#endif
#ifdef PZDEBUG2
        if(1)
        {
            std::ofstream out("MixedMeshControlHDiv.txt");
            meshcontrol.Print(out);
        }
#endif
        
        std::cout << "MHM Hdiv Computational meshes created\n";
#ifdef PZDEBUG2
        if(1)
        {
            std::ofstream gfile("geometryMHMHdiv.txt");
            gmeshauto->Print(gfile);
            std::ofstream out_mhm("MHM_hdiv.txt");
            meshcontrol.CMesh()->Print(out_mhm);
            
        }
#endif
        
        std::cout << "Number of equations MHMixed " << MHMixed->CMesh()->NEquations() << std::endl;
            
        
    }
    
    TPZCompMesh *MixedMesh = MHMixed->CMesh().operator->();
    std::ofstream multcmesh("multcompmesh.txt");
    MixedMesh->Print(multcmesh);
 
    
    return MHMixed;
}
void ComputeCoarseIndices(TPZGeoMesh *gmesh, TPZVec<int64_t> &coarseindices)
    {
        //    {
        //        std::ofstream out("gmeshref.txt");
        //        gmesh->Print(out);
        //    }
        coarseindices.Resize(gmesh->NElements());
        int count = 0;
        for (int64_t el=0; el<gmesh->NElements(); el++) {
            TPZGeoEl *gel = gmesh->Element(el);
            if(!gel || gel->Dimension() != gmesh->Dimension()) continue;
            if(gel->Father()) continue;
            coarseindices[count] = el;
            count++;
        }
        coarseindices.Resize(count);
    }





TPZGeoMesh* ConfigurateCase::CreateGeowithRefPattern(){
    TPZGeoMesh *gmeshrefpattern = CreateUniformMesh(1, 1,1,1);
    TPZGeoMesh *gmesgline =new TPZGeoMesh;
    
    TPZGeoEl *gel = gmeshrefpattern->Element(0);
    
    TPZFMatrix<REAL> cooridnates;
    gel->NodesCoordinates(cooridnates);
    TPZVec<int64_t> nodeindex(4);
    TPZVec<int64_t> nodeindexline(2);
    int nnodes = gmeshrefpattern->NNodes();
    gmeshrefpattern->NodeVec().Resize(nnodes +2);
    TPZVec<REAL> x(3,0);
    x[0]= 0.0;
    x[1]= 0.5*(cooridnates(1,0)+cooridnates(1,2));
    gmeshrefpattern->NodeVec()[nnodes].Initialize(nnodes, x, *gmeshrefpattern);
    x[0]=cooridnates(0,1);
    x[1]=0.5*(cooridnates(1,1)+cooridnates(1,3));
    gmeshrefpattern->NodeVec()[nnodes+1].Initialize(nnodes+1, x, *gmeshrefpattern);
    //
    gmesgline->NodeVec().Resize(2);
    x[0]=0.0;
    x[1]=1.0;
    gmesgline->NodeVec()[0].Initialize(0, x, *gmesgline);
    x[0]=0.5;
    x[1]=1.0;
    gmesgline->NodeVec()[1].Initialize(1, x, *gmesgline);
    nodeindexline[0]=0;
    nodeindexline[1]=1;
    int64_t indexline=0;
    gmesgline->CreateGeoElement(EOned, nodeindexline, 1,indexline );
    indexline=1;
    gmesgline->CreateGeoElement(EOned, nodeindexline, 1,indexline );
    gmesgline->Element(1)->SetFather(gmesgline->Element(0));
    gmesgline->BuildConnectivity();
    
    TPZAutoPointer<TPZRefPattern> refpatline = new TPZRefPattern(*gmesgline);
    //
    
    int64_t index = gmeshrefpattern->NElements();
    nodeindex[0]=0;
    nodeindex[1]=1;
    nodeindex[2]=5;
    nodeindex[3]=4;
    gmeshrefpattern->CreateGeoElement(EQuadrilateral, nodeindex, 1,index );
    gmeshrefpattern->Element(index)->SetFather(gmeshrefpattern->Element(0));
    
    index = gmeshrefpattern->NElements();
    nodeindex[0]=4;
    nodeindex[1]=5;
    nodeindex[2]=3;
    nodeindex[3]=2;
    gmeshrefpattern->CreateGeoElement(EQuadrilateral, nodeindex, 1,index );
    gmeshrefpattern->Element(index)->SetFather(gmeshrefpattern->Element(0));
    
    
    
    gmeshrefpattern->BuildConnectivity();
    
    TPZAutoPointer<TPZRefPattern> refpat = new TPZRefPattern(*gmeshrefpattern);
    
    
    TPZGeoMesh *gmesh2 = CreateUniformMesh(2, 1,1,1);
    gmesh2->Element(0)->SetRefPattern(refpat);
    gmesh2->Element(1)->SetRefPattern(refpat);
    TPZVec<TPZGeoEl *> sons;
    gmesh2->Element(0)->Divide(sons);
    gmesh2->Element(1)->Divide(sons);
    //    gmesh2->Element(2)->Divide(sons);
    //    gmesh2->Element(3)->Divide(sons);
    //    gmesh2->Element(4)->Divide(sons);
    //    gmesh2->Element(5)->Divide(sons);
    //    gmesh2->Element(6)->Divide(sons);
    //    gmesh2->Element(7)->Divide(sons);
    
    //    gmesh2->Element(0)->CreateBCGeoEl(4, -1);
    //    gmesh2->Element(0)->CreateBCGeoEl(6, -3);
    gmesh2->BuildConnectivity();
    std::ofstream file2("WITHREFPATTERNinitial.vtk");
    TPZVTKGeoMesh::PrintGMeshVTK(gmesh2, file2);
    gmesh2->Element(0)->CreateBCGeoEl(4, -1);
    gmesh2->Element(0)->CreateBCGeoEl(6, -3);
    gmesh2->Element(0)->CreateBCGeoEl(7, -4);
    
    gmesh2->Element(1)->CreateBCGeoEl(4, -1);
    gmesh2->Element(1)->CreateBCGeoEl(5, -2);
    gmesh2->Element(1)->CreateBCGeoEl(6, -3);
    
    gmesh2->Element(8)->Divide(sons);
    gmesh2->Element(10)->Divide(sons);
    
    
    
    //    gmesh2->Element(2)->CreateBCGeoEl(4, -1);
    //    gmesh2->Element(2)->CreateBCGeoEl(7, -4);
    //
    //    gmesh2->Element(3)->CreateBCGeoEl(6, -3);
    //    gmesh2->Element(3)->CreateBCGeoEl(7, -4);
    //
    //    gmesh2->Element(4)->CreateBCGeoEl(4, -5);
    //    gmesh2->Element(4)->CreateBCGeoEl(5, -2);
    //    gmesh2->Element(5)->CreateBCGeoEl(5, -2);
    //    gmesh2->Element(5)->CreateBCGeoEl(6, -6);
    
    //    gmesh2->Element(8)->Divide(sons);
    //    gmesh2->Element(10)->Divide(sons);
    
    //    gmesh2->Element(13)->SetMaterialId(2);
    //    gmesh2->Element(14)->SetMaterialId(2);
    //    gmesh2->Element(19)->SetMaterialId(-5);
    //    gmesh2->Element(24)->SetMaterialId(-6);
    std::set<int> matids;
    matids.insert(-1);
    //    matids.insert(1);
    
    
    
    
    //    TPZRefPatternTools::RefineDirectional(gmesh2,matids,1 );
    
    
    std::ofstream file3("WITHREFPATTERNinitial2.vtk");
    TPZVTKGeoMesh::PrintGMeshVTK(gmesh2, file3);
    
    //    gmesh2->Element(3)->CreateBCGeoEl(7, -4);
    //    gmesh2->Element(2)->CreateBCGeoEl(7, -4);
    //
    //    gmesh2->Element(2)->CreateBCGeoEl(4, -1);
    //    gmesh2->Element(4)->CreateBCGeoEl(4, -1);
    //
    //    gmesh2->Element(4)->CreateBCGeoEl(5, -2);
    //    gmesh2->Element(5)->CreateBCGeoEl(5, -2);
    
    //    gmesh2->Element(1)->CreateBCGeoEl(6, -3);
    
    //    gmesh2->Element(8)->Divide(sons);
    //    gmesh2->Element(10)->Divide(sons);
    
    gmesh2->Element(6)->SetRefPattern(refpatline);
    gmesh2->Element(7)->SetRefPattern(refpatline);
    gmesh2->Element(9)->SetRefPattern(refpatline);
    gmesh2->Element(11)->SetRefPattern(refpatline);
    
    TPZVec<TPZGeoEl *> sons2;
    gmesh2->Element(6)->Divide(sons2);
    gmesh2->BuildConnectivity();
    gmesh2->Element(7)->Divide(sons2);
    gmesh2->BuildConnectivity();
    gmesh2->Element(9)->Divide(sons2);
    gmesh2->BuildConnectivity();
    gmesh2->Element(11)->Divide(sons2);
    gmesh2->BuildConnectivity();
    
    gmesh2->Element(18)->SetMaterialId(-5);
    gmesh2->Element(19)->SetMaterialId(-6);
    
    gmesh2->Element(4)->SetMaterialId(2);
    gmesh2->Element(5)->SetMaterialId(2);
    //    gmesh2->BuildConnectivity();
    //    gmesh2->Element(17)->SetFather(gmesh2->Element(7));
    //    gmesh2->Element(16)->SetFather(gmesh2->Element(6));
    //    gmesh2->Element(18)->SetFather(gmesh2->Element(9));
    //    gmesh2->Element(19)->SetFather(gmesh2->Element(11));
    //
    //    gmesh2->Element(7)->SetSubElement(0, gmesh2->Element(17));
    //    gmesh2->Element(6)->SetSubElement(0, gmesh2->Element(16));
    //    gmesh2->Element(9)->SetSubElement(0, gmesh2->Element(18));
    //    gmesh2->Element(11)->SetSubElement(0, gmesh2->Element(19));
    //
    //    gmesh2->Element(0)->SetRefPattern(NULL);
    //
    //    gmesh2->Element(18)->SetMaterialId(-5);
    //    gmesh2->Element(19)->SetMaterialId(-6);
    //    gmesh2->Element(9)->SetMaterialId(-5);
    //    gmesh2->Element(11)->SetMaterialId(-6);
    
    //    gmesh2->Element(9)->SetMaterialId(-5);
    //    gmesh2->Element(11)->SetMaterialId(-6);
    //    gmesh2->Element(18)->SetMaterialId(-5);
    //    gmesh2->Element(19)->SetMaterialId(-6);
    //    gmesh2->Element(7)->Divide(sons);
    //    gmesh2->Element(4)->Divide(sons);
    
    //    gmesh2->Element(6)->CreateBCGeoEl(2, -1);
    //    gmesh2->Element(16)->SetFather(gmesh2->Element(6));
    
    
    
    
    std::ofstream file("WITHREFPATTERN.vtk");
    TPZVTKGeoMesh::PrintGMeshVTK(gmesh2, file);
    
    std::ofstream filetxt("mallarefpattern.txt");
    gmesh2->Print();
   
    return gmesh2;
}
void ConfigurateCase::InsertMaterialObjects(TPZMHMeshControl &control)
{
    TPZCompMesh &cmesh = control.CMesh();
    
    TPZGeoMesh &gmesh = control.GMesh();
    const int typeFlux = 1, typePressure = 0;
    TPZFMatrix<STATE> val1(1,1,0.), val2Flux(1,1,0.), val2Pressure(1,1,10.);
    val2Pressure(0,0) = 1000.;
    
    int dim = gmesh.Dimension();
    cmesh.SetDimModel(dim);
    
    TPZCompMesh *MixedFluxPressureCmesh = &cmesh;
    
    // Material medio poroso
    TPZMixedDarcyWithFourSpaces * mat = new TPZMixedDarcyWithFourSpaces(1,dim);
//    mat->SetSymmetric();
    mat->SetPermeability(1.);
    //    mat->SetForcingFunction(One);
    MixedFluxPressureCmesh->InsertMaterialObject(mat);
    
    TPZMixedDarcyWithFourSpaces * mat_2 = new TPZMixedDarcyWithFourSpaces(2,dim);
//    mat_2->SetSymmetric();
    mat_2->SetPermeability(1.); //cambiar para tener contraste de perm
    //    mat->SetForcingFunction(One);
    MixedFluxPressureCmesh->InsertMaterialObject(mat_2);
    
    // Bc N
    TPZBndCond * bcN = mat->CreateBC(mat, -1, typeFlux, val1, val2Flux);
    //    bcN->SetForcingFunction(0, force);
    
    MixedFluxPressureCmesh->InsertMaterialObject(bcN);
    bcN = mat->CreateBC(mat, -3, typeFlux, val1, val2Flux);
    //    bcN->SetForcingFunction(0, force);
    
    MixedFluxPressureCmesh->InsertMaterialObject(bcN);
    
    // Bc S
    TPZBndCond * bcS = mat->CreateBC(mat, -2, typeFlux, val1, val2Flux);
    
    MixedFluxPressureCmesh->InsertMaterialObject(bcS);
    bcS = mat->CreateBC(mat, -4, typeFlux, val1, val2Flux);
    MixedFluxPressureCmesh->InsertMaterialObject(bcS);
    
    TPZBndCond * bcIn = mat->CreateBC(mat, -5, typePressure, val1, val2Pressure);
    
    MixedFluxPressureCmesh->InsertMaterialObject(bcIn);
    
    TPZBndCond * bcOut = mat->CreateBC(mat, -6, typePressure, val1, val2Flux);
    
    MixedFluxPressureCmesh->InsertMaterialObject(bcOut);
    
}
