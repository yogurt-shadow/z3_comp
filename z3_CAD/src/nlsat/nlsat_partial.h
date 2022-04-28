#pragma once

#include "math/polynomial/polynomial_cache.h"
#include "math/polynomial/algebraic_numbers.h"
#include "nlsat/nlsat_scoped_literal_vector.h"
#include "nlsat/nlsat_assignment.h"

namespace nlsat {
    class projection_manager {
    public:
        struct imp;
    private:
        imp * m_imp;
    public:
        projection_manager(pmanager & pm, anum_manager & am, assignment const & ass, solver &, bool set_step);
        ~projection_manager();

        var heuristic_partial_level(var_vector const & m_vars, var_vector const & m_index) const;
        void mk_partial_lemma(var_vector const & m_vars, var_vector const & m_index, var level, scoped_literal_vector & res, svector<char> & m_already_added_literal);
    };
};