//
//  RSimulationCase.h
//  MonophasicTest
//
//  Created by Jose on 5/9/19.
//

#ifndef RSimulationCase_h
#define RSimulationCase_h

struct SimulationCase {
    bool            IsMHMQ;
    bool            IsCondensedQ;
    bool            UsePardisoQ;
    bool            UseFrontalQ;
    bool            KeepOneLagrangianQ;
    bool            KeepMatrixQ;
    int             elemen_type;
    int             n_h_levels;
    int             n_p_levels;
    int             n_acc_terms;
    int             order_qfine;
    int             order_qcoarse;
    int             order_p;
    int             n_threads;
    int             perturbation_type;
    std::string     mesh_type;
    std::string     domain_type;
    std::string     conv_summary;
    std::string     dump_folder;
    TPZStack<int>   omega_ids;
    TPZStack<int>   omega_dim;
    TPZStack<int>   gamma_ids;
    TPZStack<int>   gamma_dim;
    TPZStack<REAL>   permeabilities;
    TPZStack<REAL>   porosities;
    TPZStack<REAL>   type;
    TPZStack<REAL>   vals;
    REAL            c_inlet;
    
    SimulationCase() : IsMHMQ(false), UsePardisoQ(true), IsCondensedQ(true),UseFrontalQ(false), KeepOneLagrangianQ(true), KeepMatrixQ(true), elemen_type(0), n_h_levels(0), n_p_levels(1), n_acc_terms(0), order_qfine(2),order_qcoarse(1),order_p(1), n_threads(0),perturbation_type(0), mesh_type(""), domain_type(""),conv_summary(""),dump_folder(""),omega_ids(),omega_dim(),gamma_ids(), gamma_dim(), permeabilities(),porosities(), type(), vals(), c_inlet(0)
    {
        
    }
    
    SimulationCase(const SimulationCase &copy) :IsCondensedQ(copy.IsCondensedQ) ,IsMHMQ(copy.IsMHMQ), UsePardisoQ(copy.UsePardisoQ), UseFrontalQ(copy.UseFrontalQ),
    KeepOneLagrangianQ(copy.KeepOneLagrangianQ), KeepMatrixQ(copy.KeepMatrixQ), elemen_type(copy.elemen_type), n_h_levels(copy.n_h_levels), n_p_levels(copy.n_p_levels), n_acc_terms(copy.n_acc_terms), order_qfine(copy.order_qfine),order_qcoarse(copy.order_qcoarse),order_p(copy.order_p),n_threads(copy.n_threads),perturbation_type(copy.perturbation_type), mesh_type(copy.mesh_type), domain_type(copy.domain_type), conv_summary(copy.conv_summary),
    dump_folder(copy.dump_folder), omega_ids(copy.omega_ids), omega_dim(copy.omega_dim), gamma_ids(copy.gamma_ids), gamma_dim(copy.gamma_dim),
    permeabilities(copy.permeabilities),porosities(copy.porosities),
    type(copy.type),
    vals(copy.vals),c_inlet(copy.c_inlet)
    {
        
    }
    
    SimulationCase &operator=(const SimulationCase &copy)
    {
        /// check for self-assignment
        if(&copy == this){
            return *this;
        }
        
        IsMHMQ = copy.IsMHMQ;
        IsCondensedQ = copy.IsCondensedQ;
        UsePardisoQ = copy.UsePardisoQ;
        UseFrontalQ = copy.UseFrontalQ;
        KeepOneLagrangianQ = copy.KeepOneLagrangianQ;
        KeepMatrixQ = copy.KeepMatrixQ;
        elemen_type = copy.elemen_type;
        n_h_levels = copy.n_h_levels;
        n_p_levels = copy.n_p_levels;
        n_acc_terms = copy.n_acc_terms;
        order_qfine = copy.order_qfine;
        order_qcoarse = copy.order_qcoarse;
        order_p = copy.order_p;
        n_threads = copy.n_threads;
        perturbation_type = copy.perturbation_type;
        mesh_type = copy.mesh_type;
        domain_type = copy.domain_type;
        conv_summary = copy.conv_summary;
        dump_folder = copy.dump_folder;
        omega_ids = copy.omega_ids;
        omega_dim = copy.omega_dim;
        gamma_ids = copy.gamma_ids;
        gamma_dim = copy.gamma_dim;
        permeabilities=copy.permeabilities;
        porosities = copy.porosities;
        type=copy.type;
        vals=copy.vals;
        c_inlet=copy.c_inlet;
        return *this;
    }
};
#endif /* RSimulationCase_h */
