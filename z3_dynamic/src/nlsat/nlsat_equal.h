#pragma once

#include "math/polynomial/polynomial_cache.h"
#include "math/polynomial/algebraic_numbers.h"
#include "nlsat/nlsat_scoped_literal_vector.h"
#include "nlsat/nlsat_solver.h"

namespace nlsat {
    class equal_constraints {
    public:
        struct imp;
    private:
        imp * m_imp;
    public:
        equal_constraints(solver &, pmanager &, polynomial::cache &);

        ~equal_constraints();

        void heuristic_equal(polynomial_ref_vector const &, polynomial_ref_vector &, var_vector const &, var_vector const &);
    };
};