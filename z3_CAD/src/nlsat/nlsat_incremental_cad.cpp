#include "nlsat/nlsat_incremental_cad.h"

namespace nlsat {
    class LiftingTreeNode {
    private:
        LiftingTreeNode * m_parent;
        var m_var;
        anum m_value;
        unsigned m_level;
        // the set of polynomials that vanish at the current sample point
        poly_vector m_root_of;
        vector<LiftingTreeNode *> m_children;

    public:
        friend struct incremental_cad_manager::imp;
    };

    struct incremental_cad_manager::imp {
        anum_manager & m_am;
        pmanager & m_pm;
        polynomial::cache & m_cache;

        imp(anum_manager & am, pmanager & pm, polynomial::cache & cache)
        : m_am(am),
        m_pm(pm),
        m_cache(cache)
        {}

        ~imp(){}

        assignment get_assignment(LiftingTreeNode * node){
            assignment res(m_am);
            while(node){
                res.set(node->m_var, node->m_value);
                node = node->m_parent;
            }
            return res;
        }

    };

    incremental_cad_manager::incremental_cad_manager(anum_manager & am, pmanager & pm, polynomial::cache & cache){
        m_imp = alloc(imp, am, pm, cache);
    }

    incremental_cad_manager::~incremental_cad_manager(){
        dealloc(m_imp);
    }
};