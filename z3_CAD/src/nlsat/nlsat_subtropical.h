#pragma once

#include "math/polynomial/polynomial_cache.h"
#include "math/polynomial/algebraic_numbers.h"
#include "nlsat/nlsat_types.h"
#include "nlsat/nlsat_clause.h"


namespace nlsat {
    typedef vector<var_vector> frame;
    typedef vector<frame> frame_vector;

    class subtropical_processor {
    public:
        struct imp;
    private:
        imp * m_imp;
    public:
        subtropical_processor(anum_manager &, pmanager &, atom_vector &);
        ~subtropical_processor();

        lbool process();
    };
};