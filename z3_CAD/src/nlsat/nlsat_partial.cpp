#include "nlsat/nlsat_partial.h"

namespace nlsat {
    struct projection_manager::imp {
        pmanager & m_pm;
        anum_manager & m_am;
        const assignment & m_assignment;
        solver & m_solver;
        
        unsigned m_levels;
        bool use_step;


        imp(pmanager & pm, anum_manager & am, assignment const & ass, solver & s, bool step)
        : m_pm(pm),
        m_am(am),
        m_assignment(ass),
        m_solver(s),
        use_step(step)
        {
            m_levels = UINT_MAX;
        }

        ~imp(){}

        var heuristic_partial_level(var_vector const & m_vars, var_vector const & m_index) const {
            // return m_vars.size() > 2 ? m_vars.size() - 2 : 0;
            return m_vars.size() > 100 ? 1 : 0;
        }

        // level means how many literals !(x = val(x))
        void mk_partial_lemma(var_vector const & m_vars, var_vector const & m_index, unsigned level, scoped_literal_vector & res, svector<char> & m_already_added_literal){
            for(unsigned i = 0; i < level; i++){
                var cur = m_vars[i];
                const anum & val = m_assignment.value(cur);
                SASSERT(m_assignment.is_all_rational());
                rational ra_val;
                m_am.to_rational(val, ra_val);
                rational num = numerator(ra_val);
                rational denom = denominator(ra_val);
                // !(x - val(x) == 0)
                // !(denum * x - num == 0)
                poly * p = m_pm.mk_linear(1, &denom, &cur, -num);
                bool is_even = false;
                literal l = ~m_solver.mk_ineq_literal(atom::EQ, 1, &p, &is_even);
                res.push_back(l);
                unsigned lidx = l.index();
                m_already_added_literal.setx(lidx, true, false);
            }
        }
    };

    projection_manager::projection_manager(pmanager & pm, anum_manager & am, assignment const & ass, solver & s, bool set_step){
        m_imp = alloc(imp, pm, am, ass, s, set_step);
    }

    projection_manager::~projection_manager(){
        dealloc(m_imp);
    }

    var projection_manager::heuristic_partial_level(var_vector const & m_vars, var_vector const & m_index) const {
        return m_imp->heuristic_partial_level(m_vars, m_index);
    }

    void projection_manager::mk_partial_lemma(var_vector const & m_vars, var_vector const & m_index, var level, scoped_literal_vector & res, svector<char> & m_already_added_literal){
        m_imp->mk_partial_lemma(m_vars, m_index, level, res, m_already_added_literal);
    }
};