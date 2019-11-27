//$Id: main.cc,v 1.21 2010-07-22 17:43:43 caju Exp $
#ifdef HAVE_CONFIG_H
#include <pz_config.h>
#endif

#include "ConfigurateCase.h"
#include "pzcondensedcompel.h"
#include "pzelementgroup.h"

SimulationCase SimulationCase2DMHM();
SimulationCase SimulationCase2D();

int main(){
    
    ConfigurateCase confi;
    confi.SetSimulationCase(SimulationCase2DMHM());
    confi.SetFineOrder(1);
    confi.SetCoarseOrder(1);
    
    TPZAutoPointer<TPZMHMixedMesh4SpacesControl> mhmcontrol = confi.CreateMHMMixedMesh();
    TPZCompMesh* MHMIxed= mhmcontrol->CMesh().operator->();
    
    {
        std::ofstream filefinal("MHM_Before.txt");
        MHMIxed->Print(filefinal);
        std::cout<<MHMIxed->NEquations()<<std::endl;
    }
    
    SimulationCase sim;
    bool IsCondensedQ = sim.IsCondensedQ;
//    bool IsCondensedQ = 0;
    if (1) {
        
        MHMIxed->ComputeNodElCon();
        int dim = MHMIxed->Dimension();
        int64_t nel = MHMIxed->NElements();
//        for (int64_t el =0; el<nel; el++) {
//            TPZCompEl *cel = MHMIxed->Element(el);
//            if(!cel) continue;
//            TPZMultiphysicsElement *mul = dynamic_cast<TPZMultiphysicsElement *>(cel);
//            TPZCondensedCompEl *cond = dynamic_cast<TPZCondensedCompEl *>((cel));
//            TPZSubCompMesh *subcmesh = dynamic_cast<TPZSubCompMesh *>((cel));
//            if (cond) {
//                int ok=0;
//            }
//            TPZChunkVector<TPZCompEl *> celsubvec;
//            if(subcmesh){
//                subcmesh->LoadReferences();
//                celsubvec = subcmesh->ElementVec();
//            }
//
//            int nsubel = celsubvec.NElements();
//            for (int ielsub =0 ; ielsub < nsubel; ielsub++) {
//                TPZCompEl *celsub = celsubvec[ielsub];
//                if(!celsub){
//                    continue;
//                }
//                TPZCondensedCompEl *conden = dynamic_cast<TPZCondensedCompEl *>(celsub);
//                TPZElementGroup *elgroup;
//
//                if(conden){
//                    mul =dynamic_cast<TPZMultiphysicsElement *>(conden->ReferenceCompEl());
//                    elgroup = dynamic_cast<TPZElementGroup *>(conden->ReferenceCompEl());
//                }
//                else{
//                    continue;
//                }
//
//                TPZStack<TPZCompEl *, 5> groupvec;
//                if (elgroup) {
//                    groupvec = elgroup->GetElGroup();
//                }
//                int nelgroup = groupvec.size();
//                for (int ielgr=0; ielgr<nelgroup; ielgr++) {
//                    TPZCompEl *celgr =groupvec[ielgr];
//                    if(!celgr){
//                        continue;
//                    }
//                    TPZMultiphysicsElement *mulsub = dynamic_cast<TPZMultiphysicsElement *>(celgr);
//                    TPZCompEl *celpaver = mulsub->Element(3);
//                    if (!celpaver) {
//                        continue;
//                    }
//                    TPZGeoEl *gel = celgr->Reference();
//                    if(!gel) continue;
//                    if(gel->Dimension() != dim) continue;
//                    int nc = mulsub->NConnects();
//                    mulsub->Connect(nc-1).Print(*MHMIxed);
////                    mulsub->Connect(nc-1).IncrementElConnected();
////                    mulsub->Connect(nc-1).IncrementElConnected();
//
//                    int ok=0;
//                }
//            }
//
//
//        }
        MHMIxed->ComputeNodElCon();
        int nconnects = MHMIxed->NConnects();
        for (int icon=0; icon<nconnects; icon++) {
            TPZConnect &connect = MHMIxed->ConnectVec()[icon];
            int lagrange = connect.LagrangeMultiplier();
            if (lagrange==3) {
                MHMIxed->ConnectVec()[icon].IncrementElConnected();
            }
            else{
//                MHMIxed->ConnectVec()[icon].ResetElConnected();
            }
        }
       
        std::ofstream filefinal("mhmfinalantes.txt");
        MHMIxed->Print(filefinal);
        
        // Created condensed elements for the elements that have internal nodes
        //                TPZCompMesh * cmeshaux = &MHMIxed.operator*();
        TPZCompMeshTools::CreatedCondensedElements(MHMIxed, false, false);
        
       
        std::cout<<MHMIxed->NEquations()<<std::endl;
    }
//        TPZMultiphysicsCompMesh * multcompmesh = dynamic_cast<TPZMultiphysicsCompMesh *>(MHMIxed);
    

    
    
    {
        std::ofstream filefinal("MHM_After.txt");
        MHMIxed->Print(filefinal);
        std::cout<<MHMIxed->NEquations()<<std::endl;
        std::ofstream fileflux("fluxmhm_After.txt");
        std::ofstream filepressure("pressuremhm_After.txt");
        std::ofstream filefluxavg("fluxmhmavg_After.txt");
        std::ofstream filepressureavg("pressuremhmavg_After.txt");
        std::ofstream filemhm("mhmmesh.txt");
        mhmcontrol->GetMeshes()[0]->Print(fileflux);
        mhmcontrol->GetMeshes()[1]->Print(filepressure);
        mhmcontrol->GetMeshes()[2]->Print(filefluxavg);
        mhmcontrol->GetMeshes()[3]->Print(filepressureavg);
    }

    bool shouldrenumber = true;
    TPZAnalysis an_coarse(MHMIxed,shouldrenumber);
    
    TPZSymetricSpStructMatrix strmat(MHMIxed);
    strmat.SetNumThreads(2);
    
    an_coarse.SetStructuralMatrix(strmat);
    TPZStepSolver<STATE> step;
    step.SetDirect(ELDLt);
    an_coarse.SetSolver(step);
    std::cout << "Assembling\n";
    an_coarse.Assemble();
    std::ofstream filemate("MatrixCoarse.txt");
    an_coarse.Solver().Matrix()->Print("EkRs",filemate,EMathematicaInput);
    
    std::cout << "Solving\n";
    an_coarse.Solve();
    std::cout << "Finished\n";
    an_coarse.LoadSolution(); // Compute internal dofs
    
    
//    TPZBuildMultiphysicsMesh::TransferFromMultiPhysics(compmeshes, cmesh);
    
    {
        std::ofstream fileflux("fluxmhm.txt");
        std::ofstream filepressure("pressuremhm.txt");
        std::ofstream filemhm("mhmmesh.txt");
        mhmcontrol->GetMeshes()[0]->Print(fileflux);
        mhmcontrol->GetMeshes()[1]->Print(filepressure);
        mhmcontrol->CMesh()->Print(filemhm);
    }
    
    // PostProcess
    TPZStack<std::string> scalar, vectors;
    TPZManVector<std::string,10> scalnames(4), vecnames(1);
    vecnames[0]  = "q";             //Flux
    scalnames[0] = "p";             //Pressure
    scalnames[1] = "kappa";         //Permeability
    scalnames[1] = "div_q";         //Flux Divergent
    scalnames[2] = "g_average";     //Distributed flux
    scalnames[3] = "u_average";     //Average pressure
    
    
    
    {
        std::ofstream out("gmesh.vtk");
        TPZVTKGeoMesh::PrintGMeshVTK(MHMIxed->Reference(), out, true);
    }
    
    std::ofstream filefinal("mhmfinal.txt");
    MHMIxed->Print(filefinal);
    std::cout<<MHMIxed->NEquations()<<std::endl;
    
    
    std::string name_coarse("Simulation_Results.vtk");
    
    an_coarse.DefineGraphMesh(2, scalnames, vecnames, name_coarse);
    an_coarse.PostProcess(0,2);
    
    //
    //    conf.CreateGeowithRefPattern();
    ////    conf.CreateRefPattern();
    //
    ////    conf.SetSimulationCase(SimulationCase2d());
    ////    TPZMultiphysicsCompMesh *cmesh_fine = conf.CreateMultCompMesh();
    ////    conf.SetFineOrder(1);
    ////    TPZMultiphysicsCompMesh *cmesh_coarse = conf.CreateMultCompMesh();
    ////    TPZAnalysis *an_coarse = conf.CreateAnalysis(cmesh_coarse);
    ////    TPZAnalysis *an_fine = conf.CreateAnalysis(cmesh_fine);
    ////    std::cout<<"numero de elementos"<<std::endl;
    ////    TPZCompEl *cel = cmesh_coarse->Element(0);
    //
    ////    conf.CreateMHMGeoMesh(0, 0, 0, 0);
    //
    ////    TPZAutoPointer<TPZMHMixedMesh4SpacesControl> mhm = conf.CreateMHMMixedMesh();
    ////   TPZAutoPointer<TPZCompMesh>   MHMIxed = mhm->CMesh();
    ////
    ////    if (1) {
    ////
    ////        MHMIxed->ComputeNodElCon();
    ////        int dim = MHMIxed->Dimension();
    ////        int64_t nel = MHMIxed->NElements();
    ////        for (int64_t el =0; el<nel; el++) {
    ////            TPZCompEl *cel = MHMIxed->Element(el);
    ////            if(!cel) continue;
    ////            TPZGeoEl *gel = cel->Reference();
    ////            if(!gel) continue;
    ////            if(gel->Dimension() != dim) continue;
    ////            int nc = cel->NConnects();
    ////            cel->Connect(nc-1).IncrementElConnected();
    ////        }
    ////
    ////        // Created condensed elements for the elements that have internal nodes
    ////        TPZCompMesh * cmeshaux = &MHMIxed.operator*();
    ////        TPZCompMeshTools::CreatedCondensedElements(cmeshaux, false, false);
    ////    }
    ////    std::ofstream filefinal("mhmfinal.txt");
    ////    MHMIxed->Print(filefinal);
    ////    std::cout<<MHMIxed->NEquations()<<std::endl;
    //////    TPZMultiphysicsCompMesh * multcompmesh = dynamic_cast<TPZMultiphysicsCompMesh *>(MHMIxed);
    ////
    ////    bool shouldrenumber = true;
    ////    TPZAnalysis an_coarse(MHMIxed,shouldrenumber);
    ////
    ////    TPZSymetricSpStructMatrix strmat(MHMIxed.operator->());
    ////    strmat.SetNumThreads(2);
    ////
    ////    an_coarse.SetStructuralMatrix(strmat);
    ////    TPZStepSolver<STATE> step;
    ////    step.SetDirect(ELDLt);
    ////    an_coarse.SetSolver(step);
    ////    std::cout << "Assembling\n";
    ////    an_coarse.Assemble();
    ////    std::ofstream filemate("MatrixCoarse.txt");
    ////    an_coarse.Solver().Matrix()->Print("EkRs",filemate,EMathematicaInput);
    ////
    ////    std::cout << "Solving\n";
    ////    an_coarse.Solve();
    ////    std::cout << "Finished\n";
    ////    an_coarse.LoadSolution(); // compute internal dofs
    ////
    ////    std::ofstream fileflux("fluxmhm.txt");
    ////    std::ofstream filepressure("pressuremhm.txt");
    ////    std::ofstream filemhm("mhmmesh.txt");
    ////
    ////    mhm->GetMeshes()[0]->Print(fileflux);
    ////    mhm->GetMeshes()[1]->Print(filepressure);
    ////    mhm->CMesh()->Print(filemhm);
    ////    TPZVec<TPZAutoPointer<TPZCompMesh>> compmeshes = mhm->GetMeshes();
    ////    TPZAutoPointer<TPZCompMesh> cmesh = mhm->CMesh();
    ////    TPZBuildMultiphysicsMesh::TransferFromMultiPhysics(compmeshes, cmesh);
    ////
    ////
    //////    PostProcess
    ////
    ////    TPZStack<std::string> scalar, vectors;
    ////    TPZManVector<std::string,10> scalnames(4), vecnames(1);
    ////    vecnames[0]  = "q";
    ////    scalnames[0] = "p";
    ////    scalnames[1] = "kappa";
    ////    scalnames[1] = "div_q";
    ////    scalnames[2] = "g_average";
    ////    scalnames[3] = "u_average";
    ////
    ////    std::string name_coarse("results.vtk");
    ////
    ////    an_coarse.DefineGraphMesh(2, scalnames, vecnames, name_coarse);
    ////    an_coarse.PostProcess(0,2);
    ////
    
    return 0;
};

SimulationCase SimulationCase2DMHM(){
    SimulationCase sim;
    sim.UsePardisoQ=true;
    sim.KeepMatrixQ =false;
    sim.KeepOneLagrangianQ = false;
    sim.IsCondensedQ = true;
    sim.n_threads = 0;
    sim.omega_ids.push_back(1);
    sim.omega_dim.push_back(2);
    sim.omega_ids.push_back(2);
    sim.omega_dim.push_back(2);
    
    int bc_non_flux = -1;
    int bc_non_flux2 = -3;
    int bc_non_flux3 = -2;
    int bc_non_flux4 = -4;
    int bc_inlet  = -5;
    int bc_outlet = -6;
    
    sim.gamma_ids.push_back(bc_non_flux);
    sim.gamma_dim.push_back(1);
    sim.gamma_ids.push_back(bc_non_flux2);
    sim.gamma_dim.push_back(1);
    sim.gamma_ids.push_back(bc_non_flux3);
    sim.gamma_dim.push_back(1);
    sim.gamma_ids.push_back(bc_non_flux4);
    sim.gamma_dim.push_back(1);
    
    sim.gamma_ids.push_back(bc_inlet);
    sim.gamma_dim.push_back(1);
    sim.gamma_ids.push_back(bc_outlet);
    sim.gamma_dim.push_back(1);
    
    int bc_type_D = 0;
    int bc_type_N = 1;
    REAL p_inlet  = 1000.0;
    REAL p_outlet = 10.0;
    REAL qn       = 0.0;
    
    sim.type.push_back(bc_type_N);
    sim.type.push_back(bc_type_N);
    sim.type.push_back(bc_type_N);
    sim.type.push_back(bc_type_N);
    
    sim.type.push_back(bc_type_D);
    sim.type.push_back(bc_type_D);
    
    sim.vals.push_back(qn);
    sim.vals.push_back(qn);
    sim.vals.push_back(qn);
    sim.vals.push_back(qn);
    sim.vals.push_back(p_inlet);
    sim.vals.push_back(p_outlet);
    
    return sim;
};


SimulationCase SimulationCase2D(){
    SimulationCase sim;
    sim.UsePardisoQ=true;
    sim.KeepMatrixQ =false;
    sim.KeepOneLagrangianQ = false;
    sim.IsCondensedQ = true;
    sim.n_threads = 24;
    sim.omega_ids.push_back(1);
    sim.omega_dim.push_back(2);
    sim.permeabilities.push_back(1.0);
    
    
    int bc_non_flux = -1;
    int bc_inlet  = -2;
    int bc_non_flux2 = -3;
    int bc_outlet = -4;
    
    sim.gamma_ids.push_back(bc_non_flux);
    sim.gamma_dim.push_back(1);
    sim.gamma_ids.push_back(bc_inlet);
    sim.gamma_dim.push_back(1);
    sim.gamma_ids.push_back(bc_non_flux2);
    sim.gamma_dim.push_back(1);
    sim.gamma_ids.push_back(bc_outlet);
    sim.gamma_dim.push_back(1);
    
    int bc_type_D = 0;    //    D = 0;
    int bc_type_N = 1;    //    N = 1;
    REAL p_inlet  = 0.0;
    REAL p_outlet = 0.0;
    REAL qn       = 0.0;
    
    sim.type.push_back(bc_type_N);
    sim.type.push_back(bc_type_D);
    sim.type.push_back(bc_type_N);
    sim.type.push_back(bc_type_D);
    
    sim.vals.push_back(qn);
    sim.vals.push_back(p_inlet);
    sim.vals.push_back(qn);
    sim.vals.push_back(p_outlet);
    
    return sim;
};
















//
//
//
//#include "ConfigurateCase.h"
//#include "TPZDistributedMeshAlgebra.h"
//#include "pzmultiphysicselement.h"
//
////Simulations cases configuration
//SimulationCase SimulationCase1D();
//SimulationCase SimulationCase2D();
//SimulationCase SimulationCase3D();
//SimulationCase SimulationCase2DMHM();
//
//
//int main(int argc, char **argv){
//
////#ifdef LOG4CXX
////    InitializePZLOG();
////#endif
//
//
//    ConfigurateCase conf;
////    TPZGeoMesh *gmesh = conf.CreateUniformMesh(2,1,2,1);
//    conf.SetGmesh(conf.CreateGeowithRefPattern());
//    conf.SetSimulationCase(SimulationCase2DMHM());
//    conf.SetFineOrder(1);
//
//
//    TPZAutoPointer<TPZMHMixedMesh4SpacesControl> mhmcontrol = conf.CreateMHMMixedMesh();
//    TPZCompMesh* MHMIxed= mhmcontrol->CMesh().operator->();
//
//
//
//    {
//        std::ofstream filefinal("mhmfinal.txt");
//        MHMIxed->Print(filefinal);
//        std::cout<<MHMIxed->NEquations()<<std::endl;
//    }
//
//    //    TPZMultiphysicsCompMesh * multcompmesh = dynamic_cast<TPZMultiphysicsCompMesh *>(MHMIxed);
//
//    bool shouldrenumber = true;
//    TPZAnalysis an_coarse(MHMIxed,shouldrenumber);
//
//    TPZSymetricSpStructMatrix strmat(MHMIxed);
//    strmat.SetNumThreads(2);
//
//    an_coarse.SetStructuralMatrix(strmat);
//    TPZStepSolver<STATE> step;
//    step.SetDirect(ELDLt);
//    an_coarse.SetSolver(step);
//    std::cout << "Assembling\n";
//    an_coarse.Assemble();
//    std::ofstream filemate("MatrixCoarse.txt");
//    an_coarse.Solver().Matrix()->Print("EkRs",filemate,EMathematicaInput);
//
//    std::cout << "Solving\n";
//    an_coarse.Solve();
//    std::cout << "Finished\n";
//    an_coarse.LoadSolution(); // Compute internal dofs
//
//    std::ofstream fileflux("fluxmhm.txt");
//    std::ofstream filepressure("pressuremhm.txt");
//    std::ofstream filemhm("mhmmesh.txt");
//
//    //    PostProcess
//
//    TPZStack<std::string> scalar, vectors;
//    TPZManVector<std::string,10> scalnames(4), vecnames(1);
//    vecnames[0]  = "q";             //Flux
//    scalnames[0] = "p";             //Pressure
//    scalnames[1] = "kappa";         //Permeability
//    scalnames[1] = "div_q";         //Flux Divergent
//    scalnames[2] = "g_average";     //Distributed flux
//    scalnames[3] = "u_average";     //Average pressure
//
//    std::string name_coarse("Simulation_Results.vtk");
//
//    an_coarse.DefineGraphMesh(2, scalnames, vecnames, name_coarse);
//    an_coarse.PostProcess(0,2);
//
//
//    conf.CreateGeowithRefPattern();
////    conf.CreateRefPattern();
//
////    conf.SetSimulationCase(SimulationCase2d());
////    TPZMultiphysicsCompMesh *cmesh_fine = conf.CreateMultCompMesh();
////    conf.SetFineOrder(1);
////    TPZMultiphysicsCompMesh *cmesh_coarse = conf.CreateMultCompMesh();
////    TPZAnalysis *an_coarse = conf.CreateAnalysis(cmesh_coarse);
////    TPZAnalysis *an_fine = conf.CreateAnalysis(cmesh_fine);
////    std::cout<<"numero de elementos"<<std::endl;
////    TPZCompEl *cel = cmesh_coarse->Element(0);
//
////    conf.CreateMHMGeoMesh(0, 0, 0, 0);
//
////    TPZAutoPointer<TPZMHMixedMesh4SpacesControl> mhm = conf.CreateMHMMixedMesh();
////   TPZAutoPointer<TPZCompMesh>   MHMIxed = mhm->CMesh();
////
////    if (1) {
////
////        MHMIxed->ComputeNodElCon();
////        int dim = MHMIxed->Dimension();
////        int64_t nel = MHMIxed->NElements();
////        for (int64_t el =0; el<nel; el++) {
////            TPZCompEl *cel = MHMIxed->Element(el);
////            if(!cel) continue;
////            TPZGeoEl *gel = cel->Reference();
////            if(!gel) continue;
////            if(gel->Dimension() != dim) continue;
////            int nc = cel->NConnects();
////            cel->Connect(nc-1).IncrementElConnected();
////        }
////
////        // Created condensed elements for the elements that have internal nodes
////        TPZCompMesh * cmeshaux = &MHMIxed.operator*();
////        TPZCompMeshTools::CreatedCondensedElements(cmeshaux, false, false);
////    }
////    std::ofstream filefinal("mhmfinal.txt");
////    MHMIxed->Print(filefinal);
////    std::cout<<MHMIxed->NEquations()<<std::endl;
//////    TPZMultiphysicsCompMesh * multcompmesh = dynamic_cast<TPZMultiphysicsCompMesh *>(MHMIxed);
////
////    bool shouldrenumber = true;
////    TPZAnalysis an_coarse(MHMIxed,shouldrenumber);
////
////    TPZSymetricSpStructMatrix strmat(MHMIxed.operator->());
////    strmat.SetNumThreads(2);
////
////    an_coarse.SetStructuralMatrix(strmat);
////    TPZStepSolver<STATE> step;
////    step.SetDirect(ELDLt);
////    an_coarse.SetSolver(step);
////    std::cout << "Assembling\n";
////    an_coarse.Assemble();
////    std::ofstream filemate("MatrixCoarse.txt");
////    an_coarse.Solver().Matrix()->Print("EkRs",filemate,EMathematicaInput);
////
////    std::cout << "Solving\n";
////    an_coarse.Solve();
////    std::cout << "Finished\n";
////    an_coarse.LoadSolution(); // compute internal dofs
////
////    std::ofstream fileflux("fluxmhm.txt");
////    std::ofstream filepressure("pressuremhm.txt");
////    std::ofstream filemhm("mhmmesh.txt");
////
////    mhm->GetMeshes()[0]->Print(fileflux);
////    mhm->GetMeshes()[1]->Print(filepressure);
////    mhm->CMesh()->Print(filemhm);
////    TPZVec<TPZAutoPointer<TPZCompMesh>> compmeshes = mhm->GetMeshes();
////    TPZAutoPointer<TPZCompMesh> cmesh = mhm->CMesh();
////    TPZBuildMultiphysicsMesh::TransferFromMultiPhysics(compmeshes, cmesh);
////
////
//////    PostProcess
////
////    TPZStack<std::string> scalar, vectors;
////    TPZManVector<std::string,10> scalnames(4), vecnames(1);
////    vecnames[0]  = "q";
////    scalnames[0] = "p";
////    scalnames[1] = "kappa";
////    scalnames[1] = "div_q";
////    scalnames[2] = "g_average";
////    scalnames[3] = "u_average";
////
////    std::string name_coarse("results.vtk");
////
////    an_coarse.DefineGraphMesh(2, scalnames, vecnames, name_coarse);
////    an_coarse.PostProcess(0,2);
////
//
//    return 0;
//
//}
//SimulationCase SimulationCase1D(){
//    std::cout<<"Please implement me"<<std::endl;
//    DebugStop();
//};
//
//SimulationCase SimulationCase2D(){
//    SimulationCase sim;
//    sim.UsePardisoQ=true;
//    sim.KeepMatrixQ =false;
//    sim.KeepOneLagrangianQ = false;
//    sim.IsCondensedQ = true;
//    sim.n_threads = 24;
//    sim.omega_ids.push_back(1);
//    sim.omega_dim.push_back(2);
//    sim.permeabilities.push_back(1.0);
//
//
//    int bc_non_flux = -1;
//    int bc_inlet  = -2;
//    int bc_non_flux2 = -3;
//    int bc_outlet = -4;
//
//    sim.gamma_ids.push_back(bc_non_flux);
//    sim.gamma_dim.push_back(1);
//    sim.gamma_ids.push_back(bc_inlet);
//    sim.gamma_dim.push_back(1);
//    sim.gamma_ids.push_back(bc_non_flux2);
//    sim.gamma_dim.push_back(1);
//    sim.gamma_ids.push_back(bc_outlet);
//    sim.gamma_dim.push_back(1);
//
//    int bc_type_D = 0;    //    D = 0;
//    int bc_type_N = 1;    //    N = 1;
//    REAL p_inlet  = 0.0;
//    REAL p_outlet = 0.0;
//    REAL qn       = 0.0;
//
//    sim.type.push_back(bc_type_N);
//    sim.type.push_back(bc_type_D);
//    sim.type.push_back(bc_type_N);
//    sim.type.push_back(bc_type_D);
//
//    sim.vals.push_back(qn);
//    sim.vals.push_back(p_inlet);
//    sim.vals.push_back(qn);
//    sim.vals.push_back(p_outlet);
//
//    return sim;
//};
//SimulationCase SimulationCase3D(){
//    std::cout<<"Please implement me"<<std::endl;
//    DebugStop();
//};
//SimulationCase SimulationCase2DMHM(){
//    SimulationCase sim;
//    sim.UsePardisoQ=true;
//    sim.KeepMatrixQ =false;
//    sim.KeepOneLagrangianQ = false;
//    sim.IsCondensedQ = true;
//    sim.n_threads = 0;
//    sim.omega_ids.push_back(1);
//    sim.omega_dim.push_back(2);
//    sim.permeabilities.push_back(1.0);
//    sim.omega_ids.push_back(2);
//    sim.omega_dim.push_back(2);
//    sim.permeabilities.push_back(1.0);
//
//    int bc_non_flux = -1;
//    int bc_non_flux2 = -3;
//    int bc_non_flux3 = -2;
//    int bc_non_flux4 = -4;
//    int bc_inlet  = -5;
//    int bc_outlet = -6;
//
//    sim.gamma_ids.push_back(bc_non_flux);
//    sim.gamma_dim.push_back(1);
//    sim.gamma_ids.push_back(bc_non_flux2);
//    sim.gamma_dim.push_back(1);
//    sim.gamma_ids.push_back(bc_non_flux3);
//    sim.gamma_dim.push_back(1);
//    sim.gamma_ids.push_back(bc_non_flux4);
//    sim.gamma_dim.push_back(1);
//
//    sim.gamma_ids.push_back(bc_inlet);
//    sim.gamma_dim.push_back(1);
//    sim.gamma_ids.push_back(bc_outlet);
//    sim.gamma_dim.push_back(1);
//
//    int bc_type_D = 0;
//    int bc_type_N = 1;
//    REAL p_inlet  = 1000.0;
//    REAL p_outlet = 0.0;
//    REAL qn       = 0.0;
//
//    sim.type.push_back(bc_type_N);
//    sim.type.push_back(bc_type_N);
//    sim.type.push_back(bc_type_N);
//    sim.type.push_back(bc_type_N);
//
//    sim.type.push_back(bc_type_D);
//    sim.type.push_back(bc_type_D);
//
//    sim.vals.push_back(qn);
//    sim.vals.push_back(qn);
//    sim.vals.push_back(qn);
//    sim.vals.push_back(qn);
//    sim.vals.push_back(p_inlet);
//    sim.vals.push_back(p_outlet);
//
//    return sim;
//};
