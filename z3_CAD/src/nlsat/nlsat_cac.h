#pragma once

#include "math/polynomial/algebraic_numbers.h"
#include "math/polynomial/polynomial_cache.h"
#include "nlsat/nlsat_types.h"
#include "nlsat/nlsat_evaluator.h"
#include "nlsat/nlsat_interval_set.h"
#include "nlsat/nlsat_assignment.h"
#include "nlsat/nlsat_scoped_literal_vector.h"

namespace nlsat {
    // wzh nlcac
    enum ProjectionMode {
      MCCALLUM, LAZARD, LAZARDMOD
    };
    // hzw nlcac

    struct cinterval {
        bool m_lower_inf;
        bool m_upper_inf;
        bool m_lower_open;
        bool m_upper_open;
        anum m_lower;
        anum m_upper;
        poly_vector p_main;
        poly_vector p_down;
        poly_vector p_lower;
        poly_vector p_upper;

        cinterval(bool lower_inf, bool upper_inf, bool lower_open, bool upper_open, 
                poly_vector const & p1, poly_vector const & p2, poly_vector const & p3, poly_vector const & p4);
    };

    typedef vector<cinterval> cinterval_vector;

    class cac_manager {
    public:
        struct imp;
    private:
        imp * m_imp;
    public:
        cac_manager(solver &, anum_manager &, pmanager &, evaluator &, interval_set_manager & ism, polynomial::cache &, ProjectionMode);
        ~cac_manager();

        bool is_subset(cinterval const &, cinterval const &) const;
        bool contains(cinterval const &, cinterval const &) const;

        cinterval mk_cinterval(bool lower_inf, bool upper_inf, bool lower_open, bool upper_open,
                                anum const & lower, anum const & upper, 
                                poly_vector & p1, poly_vector & p2, poly_vector & p3, poly_vector & p4);

        cinterval mk_root_cinterval(anum const & w, poly_vector & p1, poly_vector & p2, poly_vector & p3, poly_vector & p4);

        bool is_valid(cinterval const &) const;
        // section means root
        bool is_section(cinterval const &) const;
        // sector means range
        bool is_sector(cinterval const &) const;
        cinterval_vector get_unsat_cintervals(literal const &, atom_vector const &, var);
        cinterval_vector get_unsat_cintervals(unsigned, literal const *, atom_vector const &, var);
        cinterval_vector split_intervals(interval_set const *, poly * p, var v);
        void collect_cac_polys(unsigned num, literal const * ls, polynomial_ref_vector & ps, atom_vector const &);
        void prune_redundant_cintervals(cinterval_vector &);
        // cac projection operator
        var construct_characterization(cinterval_vector &, polynomial_ref_vector &, var x, assignment const & m_assignment);
        // cac lifting operator
        bool all_univ(cinterval_vector const &, var) const;

        void add_cac_ineq_atom(atom::kind, var x, anum const &, scoped_literal_vector &, svector<char> &);

        cinterval interval_from_characterization(polynomial_ref_vector & ps, var x, anum const &, assignment &);

        interval_set * mk_union(cinterval_vector const &) const;

        std::ostream & display(std::ostream &, cinterval const &) const;
        std::ostream & display(std::ostream &, cinterval_vector const &) const;

        std::ostream & display(std::ostream & out, polynomial_ref_vector const & ps) const;
    };
};