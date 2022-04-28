#include "nlsat/nlsat_cac.h"

namespace nlsat {
    cinterval::cinterval(bool lower_inf, bool upper_inf, bool lower_open, bool upper_open, 
                poly_vector const & p1, poly_vector const & p2, poly_vector const & p3, poly_vector const & p4)
    {
        m_lower_inf = lower_inf;
        m_upper_inf = upper_inf;
        m_lower_open = lower_open;
        m_upper_open = upper_open;
        for(auto ele: p1){
            p_main.push_back(ele);
        }
        for(auto ele: p2){
            p_down.push_back(ele);
        }
        for(auto ele: p3){
            p_lower.push_back(ele);
        }
        for(auto ele: p4){
            p_upper.push_back(ele);
        }
    }

    struct cac_manager::imp {
        anum_manager & m_am;
        pmanager & m_pm;
        solver & m_solver;

        evaluator & m_evaluator;
        interval_set_manager & m_ism;

        polynomial::cache & m_cache;
        polynomial_ref_vector m_psc_tmp;
        polynomial_ref_vector m_factors;
        scoped_anum_vector      m_roots;

        ProjectionMode m_mode;

        imp(solver & s, anum_manager & am, pmanager & pm, evaluator & eva, interval_set_manager & ism, polynomial::cache & cache, ProjectionMode mode)
        : m_am(am),
        m_pm(pm),
        m_solver(s),
        m_evaluator(eva),
        m_ism(ism),
        m_cache(cache),
        m_psc_tmp(pm),
        m_factors(pm),
        m_roots(am),
        m_mode(mode)
        {}

        ~imp(){}

        // main down lower upper
        cinterval mk_cinterval(bool lower_inf, bool upper_inf, bool lower_open, bool upper_open, 
                            anum const & lower, anum const & upper, 
                            poly_vector & p1, poly_vector & p2, poly_vector & p3, poly_vector & p4) const {
            cinterval res(lower_inf, upper_inf, lower_open, upper_open, p1, p2, p3, p4);
            for(auto ele: p1){
                m_pm.inc_ref(ele);
            }
            for(auto ele: p2){
                m_pm.inc_ref(ele);
            }       
            for(auto ele: p3){
                m_pm.inc_ref(ele);
            }       
            for(auto ele: p4){
                m_pm.inc_ref(ele);
            }
            if(!lower_inf){
                m_am.set(res.m_lower, lower);
            }
            else {
                m_am.set(res.m_lower, 0);
            }
            if(!upper_inf){
                m_am.set(res.m_upper, upper);
            }
            else {
                m_am.set(res.m_upper, 0);
            }
            return res;
        }

        void copy_cinterval(cinterval & lhs, cinterval const & rhs){
            lhs.m_lower_inf = rhs.m_lower_inf;
            lhs.m_lower_open = rhs.m_lower_open;
            lhs.m_upper_inf = rhs.m_upper_inf;
            lhs.m_upper_open = rhs.m_upper_open;
            lhs.p_down = rhs.p_down;
            lhs.p_main = rhs.p_main;
            lhs.p_lower = rhs.p_lower;
            lhs.p_upper = rhs.p_upper;
            m_am.set(lhs.m_lower, rhs.m_lower);
            m_am.set(lhs.m_upper, rhs.m_upper);
        }

        // main down lower upper
        cinterval mk_root_cinterval(anum const & w, poly_vector & p1, poly_vector & p2, poly_vector & p3, poly_vector & p4){
            return mk_cinterval(false, false, false, false, w, w, p1, p2, p3, p4);
        }

        bool is_subset(cinterval const & lhs, cinterval const & rhs) const {
            return contains(rhs, lhs);
        }

        bool contains(cinterval const & lhs, cinterval const & rhs) const {
            bool lower = lhs.m_lower_inf || 
                        (!rhs.m_lower_inf &&
                                (m_am.lt(lhs.m_lower, rhs.m_lower) 
                                        || (m_am.eq(lhs.m_lower, rhs.m_lower) && !(lhs.m_lower_open && !rhs.m_lower_open))
                                )
                        );
            bool upper = lhs.m_upper_inf ||
                        (!rhs.m_upper_inf &&
                                (m_am.gt(lhs.m_upper, rhs.m_upper)
                                        || (m_am.eq(lhs.m_upper, rhs.m_upper) && !(lhs.m_upper_open && !rhs.m_upper_open))
                                )
                        );
            return lower && upper;
        }

        // whether two intervals connect
        // we assume lhs < rhs
        bool interval_connect(cinterval const & lhs, cinterval const & rhs) const {
            // TODO: ASSERT lhs < rhs
            if(lhs.m_upper_inf || rhs.m_lower_inf){
                return true;
            }
            return m_am.eq(lhs.m_upper, rhs.m_lower) ? !(lhs.m_upper_open && rhs.m_lower_open) : m_am.gt(lhs.m_upper, rhs.m_lower);
        }

        lbool intervals_cover(cinterval const & a, cinterval const & b, cinterval const & rhs) const {
            if(!interval_connect(a, b)){
                return l_undef;
            }
            cinterval c = interval_union(a, b);
            return contains(c, rhs) ? l_true : l_false;
        }

        // assume lhs < rhs
        cinterval interval_union(cinterval const & lhs, cinterval const & rhs) const {
            poly_vector poly_empty;
            return mk_cinterval(lhs.m_lower_inf, rhs.m_upper_inf, lhs.m_lower_open, rhs.m_upper_open, lhs.m_lower, rhs.m_lower, poly_empty, poly_empty, poly_empty, poly_empty);
        }

        interval_set * mk_union(cinterval_vector const & vec) const {
            interval_set * res = nullptr;
            for(unsigned i = 0; i < vec.size(); i++){
                cinterval ele = vec[i];
                interval inter;
                inter.m_lower_inf = ele.m_lower_inf;
                inter.m_upper_inf = ele.m_upper_inf;
                inter.m_lower_open = ele.m_lower_open;
                inter.m_upper_open = ele.m_upper_open;
                m_am.set(inter.m_lower, ele.m_lower);
                m_am.set(inter.m_upper, ele.m_upper);

                interval_set * curr = m_ism.mk(inter.m_lower_open, inter.m_lower_inf, inter.m_lower, inter.m_upper_open, inter.m_upper_inf, inter.m_upper, null_literal, nullptr);
                res = m_ism.mk_union(res, curr);
                // break;
            }
            return res;
        }

        bool is_valid(cinterval const & inter) const {
            if(inter.m_lower_inf || inter.m_upper_inf){
                return true;
            }
            if(m_am.gt(inter.m_lower, inter.m_upper)){
                return false;
            }
            if(m_am.eq(inter.m_lower, inter.m_upper)){
                return !inter.m_lower_open && !inter.m_upper_open;
            }
            return true;
        }

        bool is_section(cinterval const & inter) const {
            return !inter.m_lower_inf && !inter.m_upper_inf 
                    && !inter.m_lower_open && !inter.m_upper_open
                    && m_am.eq(inter.m_lower, inter.m_upper);
        }

        bool is_sector(cinterval const & inter) const {
            return is_valid(inter) && !is_section(inter);
        }

        cinterval_vector get_unsat_cintervals(literal const & l, atom_vector const & m_atoms, var v) {
            atom * a = m_atoms[l.var()];
            SASSERT(a != 0);
            TRACE("wzh", tout << "[cdcac] enter infeasible intervals\n";);
            interval_set_ref curr_set = m_evaluator.infeasible_intervals(a, l.sign(), nullptr);
            TRACE("wzh", tout << "[cdcac] exit infeasible intervals\n";);
            TRACE("wzh", tout << "[cdcac] show infeasible intervals:\n";
                m_ism.display(tout << std::endl, curr_set);
            );
            if(m_ism.is_empty(curr_set)){
                cinterval_vector empty;
                return empty;
            }
            polynomial_ref curr_poly(m_pm);
            collect_literal_polys(l, curr_poly, m_atoms);
            return split_intervals(curr_set, curr_poly, v);
        }

        cinterval_vector get_unsat_cintervals(unsigned num, literal const * ls, atom_vector const & m_atoms, var v) {
            cinterval_vector res;
            for(unsigned i = 0; i < num; i++){
                const literal & l = ls[i];
                cinterval_vector curr = get_unsat_cintervals(l, m_atoms, v);
                for(auto ele: curr){
                    res.push_back(ele);
                }
            }
            prune_redundant_cintervals(res);
            return res;
        }

        cinterval_vector split_intervals(interval_set const * st, poly * p, var v) {
            cinterval_vector res;
            for(unsigned i = 0; i < st->m_num_intervals; i++){
                interval cur = st->m_intervals[i];
                add_cinterval_from_interval(cur, p, res, v);
            }
            return res;
        }

        // main down lower upper
        void add_cinterval_from_interval(interval const & inter, poly * p, cinterval_vector & res, var v){
            poly_vector poly_empty;
            poly_vector ps, p_down;

            m_factors.reset();
            m_cache.factor(p, m_factors);
            for(auto ele: m_factors){
                if(m_pm.max_var(ele) == v){
                    ps.push_back(ele);
                }
                else{
                    p_down.push_back(ele);
                }
            }

            TRACE("wzh", tout << "[cdcac] enter cinterval from interval\n";
                tout << "show poly:\n";
                m_pm.display(tout, p);
            );
            // [a, a]
            if(!inter.m_lower_inf && !inter.m_upper_inf && !inter.m_lower_open && !inter.m_upper_open
                    && m_am.eq(inter.m_lower, inter.m_upper)){
                res.push_back(mk_root_cinterval(inter.m_lower, ps, p_down, ps, ps));
                return;
            }
            // (-oo, +oo)
            if(inter.m_lower_inf && inter.m_upper_inf){
                res.push_back(mk_cinterval(true, true, true, true, inter.m_lower, inter.m_upper, ps, p_down, poly_empty, poly_empty));
                return;
            }
            // (-oo, a)
            if(inter.m_lower_inf && inter.m_upper_open){
                res.push_back(mk_cinterval(true, false, true, true, inter.m_lower, inter.m_upper, ps, p_down, poly_empty, ps));
                return;
            }
            // (-oo, a]
            if(inter.m_lower_inf && !inter.m_upper_open){
                res.push_back(mk_cinterval(true, false, true, true, inter.m_lower, inter.m_upper, ps, p_down, poly_empty, ps));
                res.push_back(mk_root_cinterval(inter.m_upper, ps, p_down, ps, ps));
                return;
            }
            // (a, +oo)
            if(inter.m_lower_open && inter.m_upper_inf){
                res.push_back(mk_cinterval(false, true, true, true, inter.m_lower, inter.m_upper, ps, p_down, ps, poly_empty));
                return;
            }
            // [a, +oo)
            if(!inter.m_lower_open && inter.m_upper_inf){
                res.push_back(mk_root_cinterval(inter.m_lower, ps, p_down, ps, ps));
                res.push_back(mk_cinterval(false, true, true, true, inter.m_lower, inter.m_upper, ps, p_down, ps, poly_empty));
                return;
            }
            // (a, b)
            if(inter.m_lower_open && inter.m_upper_open){
                res.push_back(mk_cinterval(false, false, true, true, inter.m_lower, inter.m_upper, ps, p_down, ps, ps));
                return;
            }
            // (a, b]
            if(inter.m_lower_open && !inter.m_upper_open){
                res.push_back(mk_cinterval(false, false, true, true, inter.m_lower, inter.m_upper, ps, p_down, ps, ps));
                res.push_back(mk_root_cinterval(inter.m_upper, ps, p_down, ps, ps));
                return;
            }
            // [a, b)
            if(!inter.m_lower_open && inter.m_upper_open){
                res.push_back(mk_root_cinterval(inter.m_lower, ps, p_down, ps, ps));
                res.push_back(mk_cinterval(false, false, true, true, inter.m_lower, inter.m_upper, ps, p_down, ps, ps));
                return;
            }
            // [a, b]
            if(!inter.m_lower_open && !inter.m_upper_open){
                res.push_back(mk_root_cinterval(inter.m_lower, ps, p_down, ps, ps));
                res.push_back(mk_cinterval(false, false, true, true, inter.m_lower, inter.m_upper, ps, p_down, ps, ps));
                res.push_back(mk_root_cinterval(inter.m_upper, ps, p_down, ps, ps));
                return;
            }
            UNREACHABLE();
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

        // In this version, we do not collect factors but product of all polynomials
        void collect_cac_polys(unsigned num, literal const * ls, polynomial_ref_vector & ps, atom_vector const & m_atoms) {
            ps.reset();
            for (unsigned i = 0; i < num; i++) {
                const literal & l = ls[i];
                polynomial_ref curr(m_pm);
                collect_literal_polys(l, curr, m_atoms);
                ps.push_back(curr);
            }
        }

        void collect_literal_polys(literal const & l, polynomial_ref & ps, atom_vector const & m_atoms){
            ps.reset();
            atom * a = m_atoms[l.var()];
            SASSERT(a != 0);
            if (a->is_ineq_atom()) {
                ineq_atom * a_ = to_ineq_atom(a);
                ineq_products(a_, ps);
            }
            else {
                ps = to_root_atom(a)->p();
            }
        }

        void prune_redundant_cintervals(cinterval_vector & vec){
            clean_cintervals(vec);
            remove_redundant_cintervals(vec);
        }

        void heap_sort(cinterval_vector & vec){
            int size = vec.size();
            if(size <= 1){
                return;
            }
            for(int i = (size - 1) / 2; i >= 0; i--){
                shiftDown(vec, size, i);
            }
            for(int i = size - 1; i > 0; i--){
                std::swap(vec[0], vec[i]);
                shiftDown(vec, i, 0);
            }
        }

        void shiftDown(cinterval_vector & vec, int size, int i){
            unsigned k = i;
            while(2 * k + 1 < size){
                unsigned j = 2 * k + 1;
                if((j + 1 < size) && clean_up_compare(vec[j], vec[j+1])){
                    j++;
                }
                if(clean_up_compare(vec[k], vec[j])){
                    std::swap(vec[k], vec[j]);
                    k = j;
                    continue;
                }
                break;
            }
        }

        /**
        * Induces an ordering on poly intervals that is suitable for redundancy
        * removal as implemented in clean_intervals.
        */
        bool clean_up_compare(cinterval const & lhs, cinterval const & rhs) const {
            // lower bound is smaller
            if(lhs.m_lower_inf && !rhs.m_lower_inf){
                return true;
            }
            // lower bound is larger
            if(!lhs.m_lower_inf && rhs.m_lower_inf){
                return false;
            }
            if(!lhs.m_lower_inf && !rhs.m_lower_inf){
                if(!m_am.eq(lhs.m_lower, rhs.m_lower)){
                    return m_am.lt(lhs.m_lower, rhs.m_lower);
                }
            }
            // lower bound type is smaller
            if(!lhs.m_lower_open && rhs.m_lower_open){
                return true;
            }
            // lower bound type is larger
            if(lhs.m_lower_open && !rhs.m_lower_open){
                return false;
            }

            // upper bound is larger
            if(lhs.m_upper_inf && !rhs.m_upper_inf){
                return true;
            }
            // upper bound is smaller
            if(!lhs.m_upper_inf && rhs.m_upper_inf){
                return false;
            }
            if(!lhs.m_upper_inf && !rhs.m_upper_inf){
                if(!m_am.eq(lhs.m_upper, rhs.m_upper)){
                    return m_am.gt(lhs.m_upper, rhs.m_upper);
                }
            }
            // upper bound type is larger
            if(!lhs.m_upper_open && rhs.m_upper_open){
                return true;
            }
            // upper bound type is smaller
            if(lhs.m_upper_open && !rhs.m_upper_open){
                return false;
            }

            // Identical
            return false;
        }

        void clean_cintervals(cinterval_vector & vec){
            if(vec.size() < 2){
                return;
            }
            TRACE("wzh", tout << "enter heap sort" << std::endl;);
            heap_sort(vec);
            TRACE("wzh", tout << "exit heap sort" << std::endl;);
            unsigned first = 0;
            while(first < vec.size() - 1){
                if(contains(vec[first], vec[first+1])){
                    break;
                }
                first++;
            }
            if(first < vec.size() - 1){
                for(unsigned i = first + 2; i < vec.size(); i++){
                    if(!contains(vec[first], vec[i])){
                        // not covered, move it to the front
                        first++;
                        copy_cinterval(vec[first], vec[i]);
                    }
                    // else interval is covered still
                }
                while(vec.size() > first + 1){
                    vec.pop_back();
                }
            }
        }

        void remove_redundant_cintervals(cinterval_vector & vec){
            unsigned mid = 1, right = 2, size = vec.size();
            while(right < size){
                bool found = false;
                for(unsigned r = right; r < size; r++){
                    const auto & below = vec[mid - 1];
                    const auto & middle = vec[mid];
                    const auto & above = vec[r];
                    if(intervals_cover(below, above, middle) == l_true){
                        found = true;
                        break;
                    }
                }
                if(found){
                    copy_cinterval(vec[mid], vec[right]);
                }
                else {
                    mid++;
                    if(mid < right){
                        copy_cinterval(vec[mid], vec[right]);
                    }
                }
                right++;
            }
            while(vec.size() > mid + 1){
                vec.pop_back();
            }
        }

        var construct_characterization(cinterval_vector & vec, polynomial_ref_vector & res, var max_x, assignment const & m_assignment){
            var next = UINT_MAX;
            TRACE("wzh", tout << "enter construct characterization" << std::endl;);
            TRACE("wzh", tout << "enter construct characterization" << std::endl;
                tout << "[nlcac] show max var: " << max_x << std::endl;
            );
            res.reset();

            for(unsigned i = 0; i < vec.size(); i++){
                // add all down polynomials
                for(const auto & ele: vec[i].p_down){
                    if(!res.contains(ele)){
                        res.push_back(ele);
                        var curr = m_pm.max_var(ele);
                        next = curr > next || next == UINT_MAX ? curr : next;
                    }
                }
                
                polynomial_ref ele(m_pm);
                polynomial_ref p_prime(m_pm);
                for(unsigned j = 0; j < vec[i].p_main.size(); j++){
                    ele = vec[i].p_main.get(j);

                    // add discriminants
                    p_prime = m_pm.derivative(ele, max_x);

                    TRACE("wzh", tout << "[nlcac] show derivative:\n";
                        tout << max_x << std::endl;
                        m_pm.display(tout, p_prime);
                    );

                    polynomial_ref_vector disc(m_pm);
                    psc_resultant(ele, p_prime, max_x, disc);
                    TRACE("wzh", 
                        tout << "show resultant: \n";
                        display(tout, disc);
                    );
                    for(auto disc_ele: disc){
                        if(!res.contains(disc_ele)){
                            res.push_back(disc_ele);
                            var curr = m_pm.max_var(disc_ele);
                            next = curr > next || next == UINT_MAX ? curr : next;
                        }
                    }
                    TRACE("wzh", tout << "[nlcac] after discriminant: \n";
                        display(tout, res);
                    );

                    // add required coefficients
                    for(auto ele_coeff: required_coefficients(ele, max_x, m_assignment)){
                        m_factors.reset();
                        m_cache.factor(ele_coeff, m_factors);
                        for(auto factor_ele: m_factors){
                            if(!res.contains(factor_ele)){
                                res.push_back(factor_ele);
                                var curr = m_pm.max_var(factor_ele);
                                next = curr > next || next == UINT_MAX ? curr : next;
                            }
                        }
                    }
                    TRACE("wzh", tout << "[nlcac] after coefficients: \n";
                        display(tout, res);
                    );

                    // lower polys
                    for(unsigned lower_index = 0; lower_index < vec[i].p_lower.size(); lower_index++){
                        polynomial_ref curr_lower(m_pm);
                        curr_lower = vec[i].p_lower.get(lower_index);
                        if(m_pm.eq(ele, curr_lower)){
                            continue;
                        }
                        if(!has_root_below(curr_lower, max_x, vec[i].m_lower, m_assignment)){
                            continue;
                        }
                        polynomial_ref_vector lower_res(m_pm);
                        psc_resultant(ele, curr_lower, max_x, lower_res);
                        for(auto res_ele: lower_res){
                            if(!res.contains(res_ele)){
                                res.push_back(res_ele);
                                var curr = m_pm.max_var(res_ele);
                                next = curr > next || next == UINT_MAX ? curr : next;
                            }
                        }
                    }

                    // upper polys
                    for(unsigned upper_index = 0; upper_index < vec[i].p_upper.size(); upper_index++){
                        polynomial_ref curr_upper(m_pm);
                        curr_upper = vec[i].p_upper.get(upper_index);
                        if(m_pm.eq(ele, curr_upper)){
                            continue;
                        }
                        if(!has_root_above(curr_upper, max_x, vec[i].m_upper, m_assignment)){
                            continue;
                        }
                        polynomial_ref_vector upper_res(m_pm);
                        psc_resultant(ele, curr_upper, max_x, upper_res);
                        for(auto res_ele: upper_res){
                            if(!res.contains(res_ele)){
                                res.push_back(res_ele);
                                var curr = m_pm.max_var(res_ele);
                                next = curr > next || next == UINT_MAX ? curr : next;
                            }
                        }
                    }
                }
            }

            // resultants of consecutive intervals
            for(int i = 0; i < vec.size() - 1; i++){
                polynomial_ref upper(m_pm), lower(m_pm);
                for(unsigned upper_i = 0; upper_i < vec[i].p_upper.size(); upper_i++){
                    upper = vec[i].p_upper.get(upper_i);
                    for(unsigned lower_i = 0; lower_i < vec[i+1].p_lower.size(); lower_i++){
                        lower = vec[i+1].p_lower.get(lower_i);
                        polynomial_ref_vector curr_res(m_pm);
                        psc_resultant(upper, lower, max_x, curr_res);
                        for(auto ele: curr_res){
                            if(!res.contains(ele)){
                                res.push_back(ele);
                                var curr = m_pm.max_var(ele);
                                next = curr > next || next == UINT_MAX ? curr : next;
                            }
                        }
                    }
                }
            }
            TRACE("wzh", tout << "exit construct characterization" << std::endl;);
            return next;
        }

        bool has_root_below(polynomial_ref const & p, var max_x, anum const & w, assignment const & m_assignment){
            m_roots.reset();
            m_am.isolate_roots(p, undef_var_assignment(m_assignment, max_x), m_roots);
            if(m_roots.size() == 0){
                return false;
            }
            return m_am.le(m_roots[0], w);
        }
        
        bool has_root_above(polynomial_ref const & p, var max_x, anum const & w, assignment const & m_assignment){
            m_roots.reset();
            m_am.isolate_roots(p, undef_var_assignment(m_assignment, max_x), m_roots);
            if(m_roots.size() == 0){
                return false;
            }
            return m_am.ge(m_roots.back(), w);
        }

        polynomial_ref_vector required_coefficients(polynomial_ref & p, var max_x, assignment const & m_assignment){
            switch(m_mode) {
                case MCCALLUM:
                    return required_coefficients_mccallum(p, max_x, m_assignment);
                case LAZARD:
                    return required_coefficients_lazard(p, max_x, m_assignment);
                case LAZARDMOD:
                    return required_coefficients_lazardmod(p, max_x, m_assignment);
                default:
                    UNREACHABLE();
            }
        }

        polynomial_ref_vector required_coefficients_mccallum(polynomial_ref & p, var max_x, assignment const & m_assignment){
            polynomial_ref_vector res(m_pm);
            
            TRACE("wzh", tout << "[poly debug] show poly:\n";
                m_pm.display(tout, p);
            );
            unsigned k = degree(p, max_x);
            SASSERT(k > 0);
            if (m_pm.nonzero_const_coeff(p, max_x, k)) {
                TRACE("nlsat_explain", tout << "constant coefficient, skipping...\n";);
                return res;
            }
            TRACE("wzh", tout << "[poly debug] show poly:\n";
                m_pm.display(tout, p);
            );

            for(int deg = k; deg >= 0; deg--){
                polynomial_ref lc(m_pm);
                lc = m_pm.coeff(p, max_x, deg);
                if(is_const(lc)){
                    break;
                }
                res.push_back(lc);
                if(m_am.eval_sign_at(lc, undef_var_assignment(m_assignment, max_x)) != 0){
                    break;
                }
                // break;
            }

            TRACE("wzh", tout << "[poly debug] show poly:\n";
                m_pm.display(tout, p);
            );
            return res;
        }

        polynomial_ref_vector required_coefficients_lazard(polynomial_ref & p, var max_x, assignment const & m_assignment){
            polynomial_ref_vector res(m_pm);
            unsigned k = degree(p, max_x);
            polynomial_ref lc(m_pm);
            lc = m_pm.coeff(p, max_x, k);
            if(is_const(lc)){
                return res;
            }
            res.push_back(lc);
            if(m_am.eval_sign_at(lc, undef_var_assignment(m_assignment, max_x)) != 0){
                return res;
            }
            polynomial_ref tc(m_pm);
            tc = m_pm.coeff(p, max_x, 0);
            if(is_const(tc)){
                return res;
            }
            res.push_back(tc);
            return res;
        }
        
        polynomial_ref_vector required_coefficients_lazardmod(polynomial_ref & p, var max_x, assignment const & m_assignment){
            polynomial_ref_vector res(m_pm);
            unsigned k = degree(p, max_x);
            polynomial_ref lc(m_pm);
            lc = m_pm.coeff(p, max_x, k);
            if(is_const(lc)){
                return res;
            }
            res.push_back(lc);
            polynomial_ref tc(m_pm);
            tc = m_pm.coeff(p, max_x, 0);
            if(is_const(tc)){
                return res;
            }
            if(m_am.eval_sign_at(lc, undef_var_assignment(m_assignment, max_x)) != 0){
                return res;
            }
            bool vanish = true;
            polynomial_ref curr_coeff(m_pm);
            for(int deg = k - 1; deg >= 0; deg--){
                curr_coeff = m_pm.coeff(p, max_x, deg);
                if(m_am.eval_sign_at(curr_coeff, undef_var_assignment(m_assignment, max_x)) != 0){
                    vanish = false;
                    break;
                }
            }
            if(!vanish){
                return res;
            }
            else {
                res.push_back(tc);
                return res;
            }
        }

        void psc_resultant(polynomial_ref & p, polynomial_ref & q, var max_x, polynomial_ref_vector & res){
            res.reset();
            m_psc_tmp.reset();
            polynomial_ref_vector & S = m_psc_tmp;
            S.reset();
            TRACE("wzh", tout << "enter psc chain\n";
                m_pm.display(tout, p);
                tout << std::endl;
                m_pm.display(tout, q);
                tout << std::endl;
                tout << "max var: " << max_x << std::endl;
                tout << "poly max var: " << m_pm.max_var(p) << std::endl;
            );
            m_cache.psc_chain(p, q, max_x, S);
            polynomial_ref curr(m_pm);
            for(unsigned i = 0; i < S.size(); i++){
                curr = S.get(i);
                TRACE("wzh", tout << "show poly:\n";
                    m_pm.display(tout, curr);
                );
                if (m_pm.is_zero(curr)) {
                    TRACE("nlsat_explain", tout << "skipping psc is the zero polynomial\n";);
                    continue;
                }
                if(m_pm.is_const(curr)){
                    return;
                }
                m_factors.reset();
                m_cache.factor(curr, m_factors);
                for(auto ele: m_factors){
                    res.push_back(ele);
                }
            }
        }

        bool all_univ(cinterval_vector const & vec, var x) const {
            polynomial_ref curr(m_pm);
            for(unsigned i = 0; i < vec.size(); i++){
                for(unsigned j = 0; j < vec[i].p_main.size(); j++){
                    curr = vec[i].p_main.get(j);
                    if(m_pm.max_var(curr) != x){
                        return false;
                    }
                    if(!m_pm.is_univariate(curr)){
                        return false;
                    }
                }
            }
            return true;
        }

        void add_cac_ineq_atom(atom::kind k, var x, anum const & w, scoped_literal_vector & res, svector<char> & already){
            rational ra_val;
            m_am.to_rational(w, ra_val);
            rational num = numerator(ra_val);
            rational denom = denominator(ra_val);
            poly * p = m_pm.mk_linear(1, &denom, &x, -num);
            bool is_even = false;
            literal l = ~m_solver.mk_ineq_literal(k, 1, &p, &is_even);
            res.push_back(l);
            unsigned lidx = l.index();
            already.setx(lidx, true, false);
        }

        cinterval interval_from_characterization(polynomial_ref_vector & ps, var x, anum const & y_val, assignment & m_assignment){
            poly_vector m, d, l, u;
            polynomial_ref curr(m_pm);
            for(unsigned i = 0; i < ps.size(); i++){
                curr = ps.get(i);
                if(m_pm.max_var(curr) == x){
                    m.push_back(curr);
                }
                else{
                    d.push_back(curr);
                }
            }
            // push_down(m, d, x);
            TRACE("wzh", tout << "[cdcac] enter interval from characterization\n";
                m_am.display(tout, y_val);
                tout << std::endl;
                tout << "var: " << x << std::endl;
                display(tout, ps);
            );
            bool lower_inf = true;
            bool upper_inf = true;
            m_roots.reset();
            scoped_anum_vector & roots = m_roots;
            scoped_anum lower(m_am);
            scoped_anum upper(m_am);

            polynomial_ref p(m_pm);
            unsigned sz = m.size();
            for (unsigned k = 0; k < sz; k++) {
                p = m.get(k);
                if (max_var(p) != x)
                    continue;
                roots.reset();
                // Variable y is assigned in m_assignment. We must temporarily unassign it.
                // Otherwise, the isolate_roots procedure will assume p is a constant polynomial.
                m_am.isolate_roots(p, undef_var_assignment(m_assignment, x), roots);
                unsigned num_roots = roots.size();
                for (unsigned i = 0; i < num_roots; i++) {
                    int s = m_am.compare(y_val, roots[i]);
                    TRACE("nlsat_explain", 
                          m_am.display_decimal(tout << "comparing root: ", roots[i]); tout << "\n";
                          m_am.display_decimal(tout << "with y_val:", y_val); 
                          tout << "\nsign " << s << "\n";
                          tout << "poly: " << p << "\n";
                          );
                    if (s == 0) {
                        // y_val == roots[i]
                        // add literal
                        // ! (y = root_i(p))
                        m_am.set(lower, y_val);
                        m_am.set(upper, y_val);
                        lower_inf = false;
                        upper_inf = false;
                        break;
                    }
                    else if (s < 0) {
                        // y_val < roots[i]
                        
                        // check if roots[i] is a better upper bound
                        if (upper_inf || m_am.lt(roots[i], upper)) {
                            upper_inf = false;
                            m_am.set(upper, roots[i]);
                        }
                    }
                    else if (s > 0) {
                        // roots[i] < y_val

                        // check if roots[i] is a better lower bound
                        if (lower_inf || m_am.lt(lower, roots[i])) {
                            lower_inf = false;
                            m_am.set(lower, roots[i]);
                        }
                    }
                }
            }
            if(!lower_inf){
                m_assignment.set(x, lower);
                for(auto ele: m){
                    scoped_anum cur_anum(m_am);
                    TRACE("wzh", tout << "[cdcac] show poly:\n";
                        m_pm.display(tout, ele);
                        tout << "show eval value:\n";
                        m_am.display(tout, cur_anum);
                    );
                    m_pm.eval(ele, m_assignment, cur_anum);
                    if(m_am.is_zero(cur_anum)){
                        l.push_back(ele);
                    }
                }
            }
            if(!upper_inf){
                m_assignment.set(x, upper);
                for(auto ele: m){
                    scoped_anum cur_anum(m_am);
                    TRACE("wzh", tout << "[cdcac] show poly:\n";
                        m_pm.display(tout, ele);
                        tout << "show eval value:\n";
                        m_am.display(tout, cur_anum);
                    );
                    m_pm.eval(ele, m_assignment, cur_anum);
                    if(m_am.is_zero(cur_anum)){
                        u.push_back(ele);
                    }
                }
            }
            m_assignment.reset(x);
            if(lower_inf && upper_inf){
                return mk_cinterval(true, true, true, true, lower, upper, m, d, l, u);
            }
            if(!lower_inf && upper_inf){
                return mk_cinterval(false, true, true, true, lower, upper, m, d, l, u);
            }
            if(lower_inf && !upper_inf){
                return mk_cinterval(true, false, true, true, lower, upper, m, d, l, u);
            }
            if(!lower_inf && !upper_inf){
                return m_am.eq(lower, upper) ? 
                mk_cinterval(false, false, false, false, lower, upper, m, d, l, u) : mk_cinterval(false, false, true, true, lower, upper, m, d, l, u);
            }
            UNREACHABLE();
        }

        void push_down(poly_vector & m, poly_vector & d, var x){
            for(auto ele: m){
                if(m_pm.max_var(ele) != x){
                    m.erase(ele);
                    d.push_back(ele);
                }
            }
        }

        std::ostream & display(std::ostream & out, cinterval const & curr) const {
            out << "[cinterval] interval: ";
            if (curr.m_lower_inf) {
                out << "(-oo, ";
            }
            else {
                if (curr.m_lower_open)
                        out << "(";
                else
                    out << "[";
                m_am.display_decimal(out, curr.m_lower);
                out << ", ";
            }
            if (curr.m_upper_inf) {
                out << "oo)";
            }
            else {
                m_am.display_decimal(out, curr.m_upper);
                if (curr.m_upper_open)
                    out << ")";
                else
                    out << "]";
            }
            out << std::endl;
            out << "[cinterval] main: ";
            display(out, curr.p_main);
            out << std::endl;
            out << "[cinterval] down: ";
            display(out, curr.p_down);
            out << std::endl;
            out << "[cinterval] upper: ";
            display(out, curr.p_upper);
            out << std::endl;
            out << "[cinterval] lower: ";
            display(out, curr.p_lower);
            out << std::endl;
            return out;
        }

        std::ostream & display(std::ostream & out, cinterval_vector const & vec) const {
            out << "[cinterval vector] size: " << vec.size() << std::endl;
            for(unsigned i = 0; i < vec.size(); i++){
                display(out, vec[i]);
            }
            out << std::endl;
            return out;
        }

        std::ostream & display(std::ostream & out, poly_vector const & ps) const {
            for(unsigned i = 0; i < ps.size(); i++){
                m_pm.display(out, ps.get(i));
                out << std::endl;
            }
            return out;
        }

        std::ostream & display(std::ostream & out, polynomial_ref_vector const & ps) const {
            for(unsigned i = 0; i < ps.size(); i++){
                m_pm.display(out, ps.get(i));
                out << std::endl;
            }
            return out;
        }
    };

    cac_manager::cac_manager(solver & s, anum_manager & am, pmanager & pm, evaluator & eva, interval_set_manager & ism, polynomial::cache & cache, ProjectionMode mode){
        m_imp = alloc(imp, s, am, pm, eva, ism, cache, mode);
    }

    cac_manager::~cac_manager(){
        dealloc(m_imp);
    }

    bool cac_manager::is_subset(cinterval const & lhs, cinterval const & rhs) const {
        return m_imp->is_subset(lhs, rhs);
    }

    bool cac_manager::contains(cinterval const & lhs, cinterval const & rhs) const {
        return m_imp->contains(lhs, rhs);
    }

    cinterval cac_manager::mk_cinterval(bool lower_inf, bool upper_inf, bool lower_open, bool upper_open, 
                            anum const & lower, anum const & upper, 
                            poly_vector & p1, poly_vector & p2, poly_vector & p3, poly_vector & p4){
        return m_imp->mk_cinterval(lower_inf, upper_inf, lower_open, upper_open, lower, upper, p1, p2, p3, p4);
    }

    cinterval cac_manager::mk_root_cinterval(anum const & w, poly_vector & p1, poly_vector & p2, poly_vector & p3, poly_vector & p4){
        return m_imp->mk_root_cinterval(w, p1, p2, p3, p4);
    }

    bool cac_manager::is_valid(cinterval const & c) const {
        return m_imp->is_valid(c);
    }

    bool cac_manager::is_sector(cinterval const & c) const {
        return m_imp->is_sector(c);
    }

    bool cac_manager::is_section(cinterval const & c) const {
        return m_imp->is_section(c);
    }

    cinterval_vector cac_manager::get_unsat_cintervals(literal const & l, atom_vector const & atoms, var v){
        return m_imp->get_unsat_cintervals(l, atoms, v);
    }

    cinterval_vector cac_manager::get_unsat_cintervals(unsigned num, literal const * l, atom_vector const & atoms, var v){
        return m_imp->get_unsat_cintervals(num, l, atoms, v);
    }

    cinterval_vector cac_manager::split_intervals(interval_set const * st, poly * p, var v){
        return m_imp->split_intervals(st, p, v);
    }

    void cac_manager::collect_cac_polys(unsigned num, literal const * ls, polynomial_ref_vector & ps, atom_vector const & atoms){
        m_imp->collect_cac_polys(num, ls, ps, atoms);
    }

    void cac_manager::prune_redundant_cintervals(cinterval_vector & vec){
        m_imp->prune_redundant_cintervals(vec);
    }

    interval_set * cac_manager::mk_union(cinterval_vector const & vec) const {
        return m_imp->mk_union(vec);
    }
    
    var cac_manager::construct_characterization(cinterval_vector & vec, polynomial_ref_vector & res, var x, assignment const & m_assignment){
        return m_imp->construct_characterization(vec, res, x, m_assignment);
    }

    bool cac_manager::all_univ(cinterval_vector const & vec, var x) const {
        return m_imp->all_univ(vec, x);
    }

    cinterval cac_manager::interval_from_characterization(polynomial_ref_vector & ps, var x, anum const & w, assignment & m_assignment){
        return m_imp->interval_from_characterization(ps, x, w, m_assignment);
    }

    void cac_manager::add_cac_ineq_atom(atom::kind k, var x, anum const & w, scoped_literal_vector & res, svector<char> & already){
        m_imp->add_cac_ineq_atom(k, x, w, res, already);
    }

    std::ostream & cac_manager::display(std::ostream & out, cinterval const & inter) const {
        return m_imp->display(out, inter);
    }

    std::ostream & cac_manager::display(std::ostream & out, cinterval_vector const & vec) const {
        return m_imp->display(out, vec);
    }

    std::ostream & cac_manager::display(std::ostream & out, polynomial_ref_vector const & ps) const {
        return m_imp->display(out, ps);
    }
};