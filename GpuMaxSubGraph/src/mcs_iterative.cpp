//
// Created by davide on 5/3/24.
//

#include "test.hpp"
#include <vector>
#include <string>
#include <queue>
#include <algorithm>
#include <chrono>  // Include la libreria chrono

using namespace std;
    std::vector<std::vector<float>> gg0;
    std::vector<std::vector<float>> gg1;
    std::vector<float> edge_label;
    vector<pair<int,int>> m_final;


    struct queue_elem{
        vector<LabelClass> labels;
        vector<pair<int,int>> m_local;
    };
    vector<queue_elem> stack;


LabelClass *select_label(std::vector<LabelClass*>& label_classes, int map_size);

bool matching( int v, int w, LabelClass lc ) {
    std::vector<int> vector;
    vector.push_back(v);
    std::vector<int>  v_ring_atoms = {};
    v_ring_atoms = lc.get_ring_match_data(vector).at(0);

    if( !v_ring_atoms.empty() ) {
        for(int x : v_ring_atoms){

            if( x == -1 ) {
                return false;
            }
            if(x == w ){
                return true;
            }
        }
        return false;
    }
    return true;
}

vector<LabelClass> genLabels(int v, int w, const vector<LabelClass>& lcs) {
            std::vector<LabelClass> l_draft;

            for(LabelClass label : lcs){

                for(float edge_l : edge_label){
                    std::vector<int> v_conn;
                    std::vector<int> w_conn;
                    std::vector<std::vector<int> > v_c_rings;

                    for(int vtx : hood(v,gg0,edge_l)){
                        if( std::find(label.g.begin(),label.g.end(),vtx) != label.g.end() ){
                            v_conn.push_back(vtx);
                        }
                    }
                    v_c_rings = label.get_ring_match_data(v_conn);

                    for(int vtx : hood(w,gg1,edge_l)){
                        if(std::find(label.h.begin(),label.h.end(),vtx) != label.h.end() ){
                            w_conn.push_back(vtx);
                        }
                    }

                    int adj;
                    if(!v_conn.empty() && !w_conn.empty()){
                        if(edge_l != 0.0 || label.adj == 1){
                            adj = 1;
                        }else{
                            adj = 0;
                        }
                        LabelClass tmp(v_conn,w_conn,v_c_rings,adj=adj, label.label);
                        l_draft.push_back(tmp);
                    }
                }
            }
    return l_draft;
}




bool solve() {
        queue_elem elem =  stack.back();
        
        stack.pop_back();

        vector<LabelClass> lcs = elem.labels;

        vector<pair<int,int>> m_local = elem.m_local;

        std::vector<LabelClass*> label_class_pointers;

        label_class_pointers.reserve(lcs.size());
        
        for (LabelClass& item : lcs) {label_class_pointers.push_back(&item);}

        LabelClass *lcc;
        int v,w;
        queue_elem qel;
        pair<int,int> m_temp;

        lcc = select_label(label_class_pointers, m_local.size());


        if ( m_local.size() + calc_bound(lcs) <= m_final.size() || ( !lcc && !m_local.empty() )  ){ if( !stack.empty() ){ return true; } return false;}

        LabelClass lc = *lcc;
        vector<int> lcg = lc.g;

        while( !lcg.empty() ){
            v = select_vertex(lcg ,gg0);
            vector<int> lch = lc.h;
             while( !lch.empty()){
                w = select_vertex(lch, gg1);
                if( matching(v,w,lc ) )
                    {
                        m_temp.first = v;
                        m_temp.second = w;
                        m_local.push_back(m_temp);
                        qel.labels = genLabels(v,w,lcs);
                        qel.m_local = m_local;
                        stack.push_back(qel);
                        if ( m_local.size() > m_final.size() ) m_final = m_local; 
                        m_local.pop_back();
                    }
                auto it = std::find(lch.begin(), lch.end(), w);
                lch.erase(it);
                }
            auto it = std::find(lcg.begin(), lcg.end(), v);
            lcg.erase(it);
        }
    

    if( !stack.empty() ){ return true; } return false;
}

vector<pair<int,int>> mcs_iterative(const std::vector<std::vector<float>>& g00, const std::vector<std::vector<float>>& g11,
                                          const std::vector<std::string>& l0, const std::vector<std::string>& l1,
                                          std::vector<std::vector<int> >& ring_classes){


    gg0 = g00;
    gg1 = g11;
    edge_label = gen_bond_labels(gg0, gg1);
    int min = std::min(l0.size(), l1.size());
    std::vector<LabelClass> initial_label_classes = gen_initial_labels(l0, l1, ring_classes);

    vector<pair<int,int>> m_local={{}};
    queue_elem elem;
    int v,w;
    
    std::vector<LabelClass*> label_class_pointers;
        label_class_pointers.reserve(initial_label_classes.size());
    for (LabelClass& item : initial_label_classes) {label_class_pointers.push_back(&item);}
    
    LabelClass *lcc;
  

        lcc = select_label(label_class_pointers, 0);
        LabelClass lc = *lcc;
        vector<int> lcg = lc.g;
        while( !lcg.empty()  ){
            v = select_vertex(lcg,gg0);
            vector<int> lch = lc.h;
             while( !lch.empty()){
                w = select_vertex(lch, gg1);
                if( matching(v,w,lc ) )
                    {
                        m_local.at(0).first = v;
                        m_local.at(0).second = w;
                        elem.labels = genLabels(v,w,initial_label_classes);  
                        elem.m_local=m_local;
                        stack.push_back(elem);
                    }
                auto it = std::find(lch.begin(), lch.end(), w);
                lch.erase(it);

             }
            auto it = std::find(lcg.begin(), lcg.end(), v);
            lcg.erase(it);
            
        }
    
    bool flag;
    do{
        flag = solve();
     }while(flag);

    vector<pair<int,int>> best = m_final;


    m_final.clear();
    
    return best;
}
