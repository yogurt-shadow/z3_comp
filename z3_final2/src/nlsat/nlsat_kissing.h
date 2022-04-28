#include "nlsat/nlsat_clause.h"

namespace nlsat {
    typedef polynomial::manager::scoped_numeral scoped_numeral;

    class Kissing_Checker {
        struct imp;
        imp *  m_imp;
    public:
        Kissing_Checker(pmanager &, clause_vector const &, atom_vector const &);
        ~Kissing_Checker();
        lbool get_answer();
    };
}