#include "nlsat/nlsat_equal.h"

namespace nlsat {

    struct equal_constraints::imp {
        pmanager & m_pm;
        polynomial::cache & m_cache;
        polynomial_ref_vector m_psc_tmp;

        polynomial_ref_vector m_choices;
        unsigned m_left;
        solver & m_solver;

        scoped_literal_vector * m_lemma;

        imp(solver & s, pmanager & pm, polynomial::cache & cache)
        : m_solver(s),
        m_pm(pm),
        m_cache(cache),
        m_psc_tmp(m_pm),
        m_choices(m_pm),
        m_lemma(nullptr)
        {
            m_left = UINT_MAX;
        }

        ~imp(){}

        /**
         * @brief Propagation for equational constraints
         * 
         * @param ps all equational constraints
         * @return Ek with max var xn
         */
        void heuristic_equal(polynomial_ref_vector const & ps, polynomial_ref_vector & res, var_vector const & m_vars, var_vector const & m_index){
            m_choices.reset();
            res.reset();
            unsigned m_num_vars = m_vars.size();
            res.resize(m_num_vars);
            if(ps.empty()){
                return;
            }
            m_left = m_num_vars;
            polynomial_ref ele(m_pm);
            TRACE("wzh", tout << "[heuristic EC] enter initial polys" << std::endl;);
            for(unsigned i = 0; i < ps.size(); i++){
                ele = ps.get(i);
                if(ele == nullptr){
                    TRACE("wzh", tout << "[heuristic EC] nullptr\n";);
                }
                m_choices.push_back(ele);
                var curr_max = m_index[m_pm.max_var(ele)];
                if(res.get(curr_max) == nullptr){
                    TRACE("wzh", tout << "[heuristic EC] add EC:\n";
                        m_pm.display(tout, ele);
                        tout << std::endl;
                    );
                    res[curr_max] = ele;
                    m_left --;
                    if(m_left == 0){
                        return;
                    }
                }
            }
            // fill in res with bfs
            psc_update(m_choices, res, m_index);
            
            TRACE("wzh", 
                for(auto ele: res){
                    if(ele != nullptr){
                        TRACE("wzh", tout << "[res loop] " << ele << std::endl;);
                    }
                }
            );
        }

        void psc_update(polynomial_ref_vector & choices, polynomial_ref_vector & res, var_vector const & m_index){
            unsigned sz = choices.size();
            polynomial_ref poly1(m_pm), poly2(m_pm);
            for(unsigned i = 0; i < sz; i++){
                poly1 = choices.get(i);
                for(unsigned j = i+1; j < sz; j++){
                    poly2 = choices.get(j);
                    polynomial_ref_vector & S = m_psc_tmp;
                    var v1 = m_pm.max_var(poly1), v2 = m_pm.max_var(poly2);
                    var v_max = v1 > v2 ? v1 : v2;
                    m_cache.psc_chain(poly1, poly2, v_max, S);
                    polynomial_ref ele(m_pm);
                    for(unsigned k = 0; k < S.size(); k++){
                        ele = S.get(k);
                        if(m_pm.is_const(ele)){
                            continue;
                        }
                        var pm_max = m_pm.max_var(ele);
                        if(!m_index.contains(pm_max)){
                            continue;
                        }
                        var curr_max = m_index[pm_max];
                        if(res.get(curr_max) == nullptr){
                            res[curr_max] = ele;
                            TRACE("wzh", tout << "[heuristic EC] add EC:\n";
                                m_pm.display(tout, ele);
                                tout << std::endl;
                            );
                            // poly1 = 0 && poly2 = 0 ==> ele = 0
                            // !(poly1 = 0) || !(poly2 = 0) || ele = 0
                            scoped_literal_vector m_lemma(m_solver);

                            poly * ptr1 = poly1.get(), * ptr2 = poly2.get(), * ptr3 = ele.get();
                            bool is_even1 = m_pm.degree(poly1, v1) % 2 == 0;
                            bool is_even2 = m_pm.degree(poly2, v2) % 2 == 0;
                            bool is_even3 = m_pm.degree(ele, pm_max) % 2 == 0;

                            literal lit1 = m_solver.mk_ineq_literal(atom::EQ, 1, &ptr1, &is_even1);
                            literal lit2 = m_solver.mk_ineq_literal(atom::EQ, 1, &ptr2, &is_even2);
                            literal lit3 = m_solver.mk_ineq_literal(atom::EQ, 1, &ptr3, &is_even3);
                            m_lemma.push_back(~lit1);
                            m_lemma.push_back(~lit2);
                            m_lemma.push_back(lit3);

                            // m_solver.add_equational_clause(m_lemma.size(), m_lemma.data(), true, nullptr);

                            m_left --;
                            if(m_left == 0){
                                return;
                            }
                        }
                    }
                }
            }
        }
    };

    equal_constraints::equal_constraints(solver & s, pmanager & pm, polynomial::cache & cache){
        m_imp = alloc(imp, s, pm, cache);
    }

    equal_constraints::~equal_constraints(){
        dealloc(m_imp);
    }

    void equal_constraints::heuristic_equal(polynomial_ref_vector const & ps, polynomial_ref_vector & res, var_vector const & m_vars, var_vector const & index){
        m_imp->heuristic_equal(ps, res, m_vars, index);
    }
};