#include "nlsat/nlsat_peo.h"

namespace nlsat {
   struct peo_calculator::imp {
       pmanager & m_pm;
       anum_manager & m_am;
       unsigned m_num_vars;

       vector<vector<unsigned>> graph;

       var_vector label;
       var_vector visited;

       imp(pmanager & pm, anum_manager & am)
       : m_pm(pm),
       m_am(am)
       {
           m_num_vars = 0;
       }

       ~imp(){}

       var_vector get_peo_order(polynomial_ref_vector const & ps, var_vector const & m_vars, var_vector const & m_index){
           m_num_vars = m_vars.size();
           build_graph(ps, m_index);
           var_vector res = MCS();
           return is_chordal(res) ? res : var_vector(0);
       }

       void build_graph(polynomial_ref_vector const & ps, var_vector const & m_index){
           graph.resize(m_num_vars, vector<unsigned>(m_num_vars, (unsigned) 0));
           polynomial_ref curr(m_pm);
           for(unsigned i = 0; i < ps.size(); i++){
               curr = ps.get(i);
               var_vector curr_vars;
               m_pm.vars(curr, curr_vars);
               for(unsigned i = 0; i < curr_vars.size(); i++){
                   var v1 = curr_vars[i];
                   for(unsigned j = i+1; j < curr_vars.size(); j++){
                       var v2 = curr_vars[j];
                       graph[m_index[v1]][m_index[v2]] = 1;
                       graph[m_index[v2]][m_index[v1]] = 1;
                   }
               }
           }
       }

       var_vector MCS(){
           var_vector res;
           label.resize(m_num_vars, 0);
           visited.resize(m_num_vars, 0);
           var cur = select_max_label();
           res.push_back(cur);
           for(unsigned i = 0; i < m_num_vars; i++){
               if(i != cur && graph[i][cur] == 1){
                   label[i]++;
               }
           }
           res.reverse();
           return res;
       }
       
       bool is_chordal(var_vector const & order){
           var_vector curr_set;   
           for(int i = m_num_vars - 1; i >= 0; i--){
               curr_set.reset();
               for(int j = m_num_vars - 1; j > i; j--){
                   if(graph[order[i]][order[j]]){
                       curr_set.push_back(order[j]);
                   }
               }

               for(unsigned p = 0; p < curr_set.size(); p++){
                   for(unsigned q = p+1; q < curr_set.size(); q++){
                       if(!graph[curr_set[p]][curr_set[q]]){
                           return false;
                       }
                   }
               }
           }
           return true;
       }

       var select_max_label(){
           var cur = UINT_MAX;
           var max_label = 0;
           for(unsigned i = 0; i < m_num_vars; i++){
               if(visited.contains(i)){
                   continue;
               }  
               if(label[i] >= max_label){
                   cur = i;
                   max_label = label[i];
               }
           }
           visited.push_back(cur);
           return cur;
       }
   };

   peo_calculator::peo_calculator(pmanager & pm, anum_manager & am){
       m_imp = alloc(imp, pm, am);
   }

   peo_calculator::~peo_calculator(){
       dealloc(m_imp);
   }

   var_vector peo_calculator::get_peo_order(polynomial_ref_vector const & ps, var_vector const & vars, var_vector const & index){
       return m_imp->get_peo_order(ps, vars, index);
   }
};