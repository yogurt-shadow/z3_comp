#include "nlsat/nlsat_subtropical.h"

namespace nlsat {
    struct subtropical_processor::imp {
        anum_manager & m_am;
        pmanager & m_pm;
        atom_vector m_atoms;
        polynomial_ref_vector m_polys;
        vector<atom::kind> m_kinds;

        var_vector m_vars;
        var_vector m_index;
        var m_max_var;

        // frame: {{x1, x2, x3}, {y1, y2, y3}, ...}
        frame_vector m_frames;
        frame_vector m_pos_frames;
        frame_vector m_neg_frames;
        

        imp(anum_manager & am, pmanager & pm, atom_vector & atoms)
        : m_am(am),
        m_pm(pm),
        m_atoms(atoms),
        m_polys(pm)
        {
            m_vars.reset();
            m_index.reset();
            m_max_var = UINT_MAX;
            m_frames.reset();
            m_pos_frames.reset();
            m_neg_frames.reset();
        }

        ~imp(){}

        void collect(){
            m_polys.reset();
            for(atom * a: m_atoms){
                polynomial_ref curr(m_pm);
                if(a->is_ineq_atom()){
                    ineq_products(to_ineq_atom(a), curr);
                    m_kinds.push_back(a->get_kind());
                }
                else {
                    UNREACHABLE();
                }
                collect(curr);
                m_polys.push_back(curr);
            }
            mk_index();
        }

        void collect(poly * p){
            m_vars.reset();
            var_vector curr;
            m_pm.vars(p, curr);
            for(var v: curr){
                if(!m_vars.contains(v)){
                    m_vars.push_back(v);
                    if(v > m_max_var || m_max_var == UINT_MAX){
                        m_max_var = v;
                    }
                }
            }
        }
        
        void mk_index(){
            m_index.resize(m_max_var + 1);
            for(unsigned i = 0; i < m_vars.size(); i++){
                m_index[m_vars[i]] = i;
            }
        }

        void ineq_products(ineq_atom * a, polynomial_ref & res){
            res.reset();
            unsigned sz = a->size();
            if(sz == 0){
                return;
            }
            res = a->p(0);
            for(unsigned i = 1; i < sz; i++){
                res = m_pm.mul(res, a->p(i));
            }
        }

        void get_frame(polynomial_ref_vector & ps){
            m_frames.reset();
            m_pos_frames.reset();
            m_neg_frames.reset();
            polynomial_ref curr(m_pm);
            for(unsigned i = 0; i < ps.size(); i++){
                curr = ps.get(i);
                get_frame(curr);
            }
        }

        void get_frame(poly * p){
            frame res, pos, neg;
            unsigned size = m_pm.size(p);
            for(unsigned i = 0; i < size; i++){
                monomial * curr = m_pm.get_monomial(p, i);
                const numeral & curr_coeff = m_pm.coeff(p, i);
                numeral curr_num;
                m_pm.m().set(curr_num, 1);
                poly * mono_poly = m_pm.mk_polynomial(1, &curr_num, &curr);
                var_vector curr_vec;
                for(var v: m_vars){
                    unsigned deg = m_pm.degree(mono_poly, v);
                    curr_vec.push_back(deg);
                }
                res.push_back(curr_vec);
                if(m_pm.m().is_pos(curr_coeff)){
                    pos.push_back(curr_vec);
                }
                else if(m_pm.m().is_neg(curr_coeff)){
                    neg.push_back(curr_vec);
                }
                else{
                    UNREACHABLE();
                }
            }
            m_frames.push_back(res);
            m_pos_frames.push_back(pos);
            m_neg_frames.push_back(neg);
        }

        lbool process(){
            collect();
            get_frame(m_polys);
            
        }
    };

    subtropical_processor::subtropical_processor(anum_manager & am, pmanager & pm, atom_vector & atoms){
        m_imp = alloc(imp, am, pm, atoms);
    }

    subtropical_processor::~subtropical_processor(){
        dealloc(m_imp);
    }

    lbool subtropical_processor::process(){
        return m_imp->process();
    }
};