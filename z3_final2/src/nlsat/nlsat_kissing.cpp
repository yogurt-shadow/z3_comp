#include "nlsat/nlsat_kissing.h"

namespace nlsat {
    struct Kissing_Checker::imp {
        pmanager & pm;
        const clause_vector & clauses;
        const atom_vector & atoms;
        unsigned clause_num;
        unsigned n, m;
        unsigned num_square, num_diff;

        enum Kind {
            SQUARE_SUM, DIFF_SUM, OTHER
        };

        vector<Kind> clause_kind;

        // n: 1 2 3 4
        const unsigned KISSING_NUM[4] = {2, 6, 12, 24};

        imp(pmanager & _pm, clause_vector const & _clauses, atom_vector const & _atoms):
        pm(_pm), 
        clauses(_clauses),
        atoms(_atoms)
        {
            clause_num = clauses.size();
            clause_kind.resize(clause_num);
            n = 0, m = 0;
            num_square = 0, num_diff = 0;
        }

        lbool get_answer(){
            if(!is_kissing()){
                // std::cout << "no kissing\n";
                return l_undef;
            }
                // std::cout << "kissing\n";
            SASSERT(n - 1 >= 0 && n - 1 < 4);
            return m > KISSING_NUM[n-1] ? l_false : l_true;
        }

        bool is_kissing(){
            // std::cout << "clause num: " << clause_num << std::endl;
            m = get_m(clause_num);
            // std::cout << "m: " << m << std::endl;
            if(m == UINT_MAX){
                return false;
            }
            for(unsigned i = 0; i < clause_num; i++){
                const clause & curr = *clauses[i];
                clause_kind = getKind(curr);
                if(clause_kind == SQUARE_SUM){
                    num_square++;
                }
                else if(clause_kind == DIFF_SUM){
                    num_diff++;
                }
                else {
                    return false;
                }
            }
            return true;
        }

        Kind getKind(clause const & cls){
            if(cls.size() != 1){
                return OTHER;
            }
            literal l = cls[0];
            atom * at = atoms[l.var()];
            if(at == nullptr || !at->is_ineq_atom()){
                return OTHER;
            }
            ineq_atom * at2 = to_ineq_atom(at);
            unsigned poly_size = at2->size();
            if(poly_size != 1){
                return OTHER;
            }
            poly * curr_poly = at2->p(0);
            return getKind(curr_poly);
        }

        Kind getKind(poly * p){
            unsigned mono_size = pm.size(p);
            unsigned cross = 0, square = 0, one = 0;
            for(unsigned i = 0; i < mono_size; i++){
                monomial * curr_mono = pm.get_monomial(p, i);
                unsigned curr_deg = pm.total_degree(curr_mono);
                if(curr_deg == 0){
                    one++;
                    if(one > 1){
                        return OTHER;
                    }
                } 
                else if(curr_deg == 2){
                    numeral curr_num;
                    pm.m().set(curr_num, 1);
                    poly * curr_poly = pm.mk_polynomial(1, &curr_num, &curr_mono);
                    var_vector curr_vars;
                    pm.vars(curr_poly, curr_vars);
                    // x^2
                    if(curr_vars.size() == 1 && pm.degree(curr_poly, curr_vars[0]) == 2){
                        square++;
                    }
                    // x*y
                    else if(curr_vars.size() == 2 && pm.degree(curr_poly, curr_vars[0]) == 1
                                                  && pm.degree(curr_poly, curr_vars[1]) == 1){
                        cross++;
                    }
                    else {
                        return OTHER;
                    }
                }   
                else {
                    return OTHER;
                }
            }
            if(square > 0 && cross == 0){
                return update_check_n(square) ? SQUARE_SUM : OTHER;
            }
            if(cross > 0 && square == 2*cross){
                return update_check_n(cross) ? DIFF_SUM : OTHER;
            }
            return OTHER;
        }

        bool update_check_n(unsigned x){
            if(n == 0){
                n = x;
            }
            return n == x;
        }

        unsigned get_m(unsigned num){
            unsigned res = 1;
            while(true){
                unsigned sum1 = res;
                unsigned sum2 = res * (res - 1) / 2;
                unsigned sum = sum1 + sum2;
                if(sum == num){
                    return res;
                }
                if(sum > num){
                    return UINT_MAX;
                }
                res++;
            }
            return UINT_MAX;
        }
    };

    Kissing_Checker::Kissing_Checker(pmanager & _pm, clause_vector const & _clauses, atom_vector const & _atoms){
        m_imp = alloc(imp, _pm, _clauses, _atoms);
    }

    Kissing_Checker::~Kissing_Checker(){
        dealloc(m_imp);
    }

    lbool Kissing_Checker::get_answer() {
        return m_imp->get_answer();
    }
};