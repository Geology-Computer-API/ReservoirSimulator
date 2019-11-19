//$Id: main.cc,v 1.21 2010-07-22 17:43:43 caju Exp $
#ifdef HAVE_CONFIG_H
#include <pz_config.h>
#endif

//#include "pzgeoel.h"
//#include "pzgnode.h"
//#include "pzgmesh.h"
//#include "pzbiharmonic.h"

#include "ConfigurateCase.h"
#include "TPZDistributedMeshAlgebra.h"
#include "pzmultiphysicselement.h"
SimulationCase SimulationCase1d();
SimulationCase SimulationCase2d();
SimulationCase SimulationCase3d();
SimulationCase SimulationCase2dMHM();
#ifdef _AUTODIFF
//#include "tfad.h"
//#include "fad.h"
//#include "pzextractval.h"
#endif

std::ofstream log_file("Results.txt");

int main(int argc, char **argv){
    
//#ifdef LOG4CXX
//    InitializePZLOG();
//#endif
    ConfigurateCase conf;
//    TPZGeoMesh *gmesh = conf.CreateUniformMesh(2,1,2,1);
    conf.SetGmesh(conf.CreateGeowithRefPattern());
    conf.SetSimulationCase(SimulationCase2dMHM());
    conf.SetFineOrder(1);
    TPZMultiphysicsCompMesh *MHMIxed = conf.CreateMultCompMesh();
    
    std::ofstream filefinal("mhmfinal.txt");
    MHMIxed->Print(filefinal);
    std::cout<<MHMIxed->NEquations()<<std::endl;
    //    TPZMultiphysicsCompMesh * multcompmesh = dynamic_cast<TPZMultiphysicsCompMesh *>(MHMIxed);
    
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
    an_coarse.LoadSolution(); // compute internal dofs
    
    std::ofstream fileflux("fluxmhm.txt");
    std::ofstream filepressure("pressuremhm.txt");
    std::ofstream filemhm("mhmmesh.txt");
    
    
    
    
    //    PostProcess
    
    TPZStack<std::string> scalar, vectors;
    TPZManVector<std::string,10> scalnames(4), vecnames(1);
    vecnames[0]  = "q";
    scalnames[0] = "p";
    scalnames[1] = "kappa";
    scalnames[1] = "div_q";
    scalnames[2] = "g_average";
    scalnames[3] = "u_average";
    
    std::string name_coarse("results.vtk");
    
    an_coarse.DefineGraphMesh(2, scalnames, vecnames, name_coarse);
    an_coarse.PostProcess(0,2);
    
    
    
    conf.CreateGeowithRefPattern();
//    conf.CreateRefPattern();
    
//    conf.SetSimulationCase(SimulationCase2d());
//    TPZMultiphysicsCompMesh *cmesh_fine = conf.CreateMultCompMesh();
//    conf.SetFineOrder(1);
//    TPZMultiphysicsCompMesh *cmesh_coarse = conf.CreateMultCompMesh();
//    TPZAnalysis *an_coarse = conf.CreateAnalysis(cmesh_coarse);
//    TPZAnalysis *an_fine = conf.CreateAnalysis(cmesh_fine);
//    std::cout<<"numero de elementos"<<std::endl;
//    TPZCompEl *cel = cmesh_coarse->Element(0);
    
//    conf.CreateMHMGeoMesh(0, 0, 0, 0);
    
//    TPZAutoPointer<TPZMHMixedMesh4SpacesControl> mhm = conf.CreateMHMMixedMesh();
//   TPZAutoPointer<TPZCompMesh>   MHMIxed = mhm->CMesh();
//    
//    if (1) {
//        
//        MHMIxed->ComputeNodElCon();
//        int dim = MHMIxed->Dimension();
//        int64_t nel = MHMIxed->NElements();
//        for (int64_t el =0; el<nel; el++) {
//            TPZCompEl *cel = MHMIxed->Element(el);
//            if(!cel) continue;
//            TPZGeoEl *gel = cel->Reference();
//            if(!gel) continue;
//            if(gel->Dimension() != dim) continue;
//            int nc = cel->NConnects();
//            cel->Connect(nc-1).IncrementElConnected();
//        }
//        
//        // Created condensed elements for the elements that have internal nodes
//        TPZCompMesh * cmeshaux = &MHMIxed.operator*();
//        TPZCompMeshTools::CreatedCondensedElements(cmeshaux, false, false);
//    }
//    std::ofstream filefinal("mhmfinal.txt");
//    MHMIxed->Print(filefinal);
//    std::cout<<MHMIxed->NEquations()<<std::endl;
////    TPZMultiphysicsCompMesh * multcompmesh = dynamic_cast<TPZMultiphysicsCompMesh *>(MHMIxed);
//
//    bool shouldrenumber = true;
//    TPZAnalysis an_coarse(MHMIxed,shouldrenumber);
//
//    TPZSymetricSpStructMatrix strmat(MHMIxed.operator->());
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
//    an_coarse.LoadSolution(); // compute internal dofs
//
//    std::ofstream fileflux("fluxmhm.txt");
//    std::ofstream filepressure("pressuremhm.txt");
//    std::ofstream filemhm("mhmmesh.txt");
//
//    mhm->GetMeshes()[0]->Print(fileflux);
//    mhm->GetMeshes()[1]->Print(filepressure);
//    mhm->CMesh()->Print(filemhm);
//    TPZVec<TPZAutoPointer<TPZCompMesh>> compmeshes = mhm->GetMeshes();
//    TPZAutoPointer<TPZCompMesh> cmesh = mhm->CMesh();
//    TPZBuildMultiphysicsMesh::TransferFromMultiPhysics(compmeshes, cmesh);
//
//
////    PostProcess
//
//    TPZStack<std::string> scalar, vectors;
//    TPZManVector<std::string,10> scalnames(4), vecnames(1);
//    vecnames[0]  = "q";
//    scalnames[0] = "p";
//    scalnames[1] = "kappa";
//    scalnames[1] = "div_q";
//    scalnames[2] = "g_average";
//    scalnames[3] = "u_average";
//
//    std::string name_coarse("results.vtk");
//
//    an_coarse.DefineGraphMesh(2, scalnames, vecnames, name_coarse);
//    an_coarse.PostProcess(0,2);
//
    
    return 0;

}
SimulationCase SimulationCase1d(){
    
    
};
SimulationCase SimulationCase2d(){
    SimulationCase sim;
    sim.UsePardisoQ=true;
    sim.KeepMatrixQ =false;
    sim.KeepOneLagrangianQ = false;
    sim.IsCondensedQ = true;
//    sim.IsHybrid=true;
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
SimulationCase SimulationCase3dd(){
    
    
};
SimulationCase SimulationCase2dMHM(){
    SimulationCase sim;
    sim.UsePardisoQ=true;
    sim.KeepMatrixQ =false;
    sim.KeepOneLagrangianQ = false;
    sim.IsCondensedQ = true;
    //    sim.IsHybrid=true;
    sim.n_threads = 24;
    sim.omega_ids.push_back(1);
    sim.omega_dim.push_back(2);
    sim.permeabilities.push_back(1.0);
    sim.omega_ids.push_back(2);
    sim.omega_dim.push_back(2);
    sim.permeabilities.push_back(1.0);
    
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
    
    int bc_type_D = 0;    //    D = 0;
    int bc_type_N = 1;    //    N = 1;
    REAL p_inlet  = 1000.0;
    REAL p_outlet = 0.0;
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
