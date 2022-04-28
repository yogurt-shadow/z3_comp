#pragma once

#include "math/polynomial/polynomial_cache.h"
#include "math/polynomial/algebraic_numbers.h"
#include "nlsat/nlsat_types.h"

namespace nlsat {
    class peo_calculator {
    public:
        struct imp;
    private:
        imp * m_imp;
    public:
        peo_calculator(pmanager &, anum_manager &);

        ~peo_calculator();

        var_vector get_peo_order(polynomial_ref_vector const &, var_vector const &, var_vector const &);
    };
};