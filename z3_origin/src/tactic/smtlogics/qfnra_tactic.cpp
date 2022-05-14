/*++
Copyright (c) 2012 Microsoft Corporation

Module Name:

    qfnra_tactic.cpp

Abstract:

    Tactic for QF_NRA

Author:

    Leonardo (leonardo) 2012-02-28

Notes:

--*/
#include "tactic/tactical.h"
#include "tactic/core/simplify_tactic.h"
#include "tactic/core/propagate_values_tactic.h"
#include "tactic/arith/nla2bv_tactic.h"
#include "nlsat/tactic/qfnra_nlsat_tactic.h"
#include "tactic/smtlogics/smt_tactic.h"

static tactic * mk_qfnra_sat_solver(ast_manager& m, params_ref const& p, unsigned bv_size) {
    params_ref nra2sat_p = p;
    nra2sat_p.set_uint("nla2bv_max_bv_size", p.get_uint("nla2bv_max_bv_size", bv_size));   
    
    return and_then(mk_nla2bv_tactic(m, nra2sat_p),
                    mk_smt_tactic(m),
                    mk_fail_if_undecided_tactic());
}

tactic * mk_nlsat_mixed_solver(ast_manager& m, params_ref const& p) {
    ptr_vector<tactic> ts;
    {
        params_ref p_sc = p;
        p_sc.set_bool("linxi_simple_check", true);
        p_sc.set_uint("seed", 997);
        ts.push_back(try_for(and_then(mk_qfnra_nlsat_tactic(m, p_sc), mk_fail_if_undecided_tactic()), 20 * 1000));
    }
    {
        params_ref p_heuristic = p;
        p_heuristic.set_uint("seed", 233);
        ts.push_back(try_for(mk_qfnra_nlsat_tactic(m, p_heuristic), 5 * 1000));

        // params_ref p_order_3 = p;
        // p_order_3.set_uint("linxi_variable_ordering_strategy", 3);
        // p_order_3.set_uint("seed", 17);
        // ts.push_back(try_for(mk_qfnra_nlsat_tactic(m, p_order_3), 10 * 1000));

        // params_ref p_order_1 = p;
        // p_order_1.set_uint("linxi_variable_ordering_strategy", 1);
        // ts.push_back(try_for(mk_qfnra_nlsat_tactic(m, p_order_1), 15 * 1000));

        // params_ref p_order_2 = p;
        // p_order_2.set_uint("linxi_variable_ordering_strategy", 2);
        // ts.push_back(try_for(mk_qfnra_nlsat_tactic(m, p_order_2), 20 * 1000));
    }
    {
        params_ref p_l = p;
        p_l.set_bool("arith.greatest_error_pivot", true);
        ts.push_back(and_then(try_for(using_params(mk_smt_tactic(m), p_l), 400 * 1000), mk_fail_if_undecided_tactic()));
    }
    // for (unsigned i = 0; i < 100; ++i) { // 5s * 100 = 500s
    //     params_ref p_i = p;
    //     p_i.set_uint("seed", i);
    //     p_i.set_bool("shuffle_vars", true);
    //     ts.push_back(try_for(mk_qfnra_nlsat_tactic(m, p_i), 5 * 1000));
    // }
    {
        ts.push_back(mk_qfnra_nlsat_tactic(m, p));
    }
    return or_else(ts.size(), ts.data());
}


tactic * mk_qfnra_tactic(ast_manager & m, params_ref const& p) {
    params_ref p0 = p;
    p0.set_bool("inline_vars", true);
    params_ref p1 = p;    
    p1.set_uint("seed", 11);
    p1.set_bool("factor", false);
    params_ref p2 = p;
    p2.set_uint("seed", 13);
    p2.set_bool("factor", false);

    return and_then(mk_simplify_tactic(m, p), 
                    mk_propagate_values_tactic(m, p),

                    mk_nlsat_mixed_solver(m, p)
    );
}


