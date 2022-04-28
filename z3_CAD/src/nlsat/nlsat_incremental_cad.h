#pragma once

#include "math/polynomial/polynomial_cache.h"
#include "math/polynomial/algebraic_numbers.h"
#include "nlsat/nlsat_scoped_literal_vector.h"
#include "nlsat/nlsat_types.h"
#include "nlsat/nlsat_assignment.h"

namespace nlsat {
    class incremental_cad_manager {
    public:
        struct imp;
    private: 
        imp * m_imp;
    public:
        incremental_cad_manager(anum_manager &, pmanager &, polynomial::cache &);
        ~incremental_cad_manager();
    };
};