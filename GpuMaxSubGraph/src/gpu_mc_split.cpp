//
// Created by davide on 5/3/24.
//
#include "gpu_header.hpp"


using namespace std;

std::vector<std::vector<float>> g0;
std::vector<std::vector<float>> g1;
std::vector<std::string> l0;
std::vector<std::string> l1;
int size_initial_label_classes ;
std::vector<float> edge_labels;
vector<pair<int,int>> m_best;
int max_first_len = 0;
int max_second_len = 0;


 
int max_depth = 0;
bool first_bb = false;


vector<queue_elem> Q;
vector<queue_elem> Q_gpu;
vector<queue_elem> Q_cpu;


int LIMIT_DEPTH = 4 ;


void sortLabels(std::vector<LabelClass>& labels){
    int index_max = 0;
    int max = 0;
    int i = 0;
    
    for(LabelClass lc : labels ){
        if(lc.g.size() > max){
            max = lc.g.size();
            index_max = i;
        }
        i++;
    }


    if(index_max != 0){
        LabelClass tmp = labels[0];
        labels[0] = labels[index_max];
        labels[index_max] = tmp;
    }


    index_max = 1;
    max = 0;
    for(int i = 1; i < labels.size(); i++){
        if(labels[i].g.size() > max){
            max = labels[i].g.size();
            index_max = i;
        }
    }

    if(index_max != 1){
        LabelClass tmp = labels[1];
        labels[1] = labels[index_max];
        labels[index_max] = tmp;
    }


}


LabelClass *select_label(std::vector<LabelClass*>& label_classes, int map_size);


bool matchable( int v, int w, LabelClass lc ) {
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

vector<LabelClass> genNewLabels(int v, int w, const vector<LabelClass>& lcs) {
            std::vector<LabelClass> l_draft;

            for(LabelClass label : lcs){

                for(float edge_l : edge_labels){
                    std::vector<int> v_conn;
                    std::vector<int> w_conn;
                    std::vector<std::vector<int> > v_c_rings;

                    for(int vtx : hood(v,g0,edge_l)){
                        if( std::find(label.g.begin(),label.g.end(),vtx) != label.g.end() ){
                            v_conn.push_back(vtx);
                        }
                    }
                    v_c_rings = label.get_ring_match_data(v_conn);

                    for(int vtx : hood(w,g1,edge_l)){
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




bool solve_mcs() {
        if(Q.empty()) return false;

        queue_elem elem =  Q.back();
        
        Q.pop_back();

        vector<LabelClass> lcs = elem.labels;

        vector<pair<int,int>> m_local = elem.m_local;

        std::vector<LabelClass*> label_class_pointers;

        label_class_pointers.reserve(lcs.size());
        
        for (LabelClass& item : lcs) {label_class_pointers.push_back(&item);}

        LabelClass *lcc;
        int v,w;
        queue_elem qel;
        pair<int,int> m_temp;

           
    // first brench&bound 
    if(first_bb == false){
        if ( m_local.size() >= max_depth) max_depth = m_local.size();
        else if( m_local.size() < max_depth){ 
            Q.push_back(elem);
            first_bb = true; 
            return false;
        }
    }

        lcc = select_label(label_class_pointers, m_local.size());


        if ( m_local.size() + calc_bound(lcs) <= m_best.size() || ( !lcc && !m_local.empty() )  ){ if( !Q.empty() ){ return true; } return false;}

        LabelClass lc = *lcc;
        vector<int> lcg = lc.g;

        while( !lcg.empty() ){
            v = select_vertex(lcg ,g0);
            vector<int> lch = lc.h;
             while( !lch.empty()){
                w = select_vertex(lch, g1);
                if( matchable(v,w,lc ) )
                    {
                        m_temp.first = v;
                        m_temp.second = w;
                        m_local.push_back(m_temp);
                        qel.labels = genNewLabels(v,w,lcs);
                        qel.m_local = m_local;
                        Q.push_back(qel);
                        if ( m_local.size() > m_best.size() ) m_best = m_local; 
                        m_local.pop_back();
                    }
                auto it = std::find(lch.begin(), lch.end(), w);
                lch.erase(it);
                }
            auto it = std::find(lcg.begin(), lcg.end(), v);
            lcg.erase(it);
        }
    

    if( !Q.empty() ){ return true; } return false;
}

// Functions identically to solve_mcs but operates with a different QUEUE ( Testing purpose )
bool solve_2(vector<queue_elem>& Q_tmp) {
        if(Q_tmp.empty()) return false;
        
        queue_elem elem =  Q_tmp.back();
        
        Q_tmp.pop_back();

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


        if ( m_local.size() + calc_bound(lcs) <= m_best.size() || ( !lcc && !m_local.empty() )  ){ if( !Q_tmp.empty() ){ return true; } return false;}

        LabelClass lc = *lcc;
        vector<int> lcg = lc.g;

        while( !lcg.empty() ){
            v = select_vertex(lcg ,g0);
            vector<int> lch = lc.h;
             while( !lch.empty()){
                w = select_vertex(lch, g1);
                if( matchable(v,w,lc ) )
                    {
                        m_temp.first = v;
                        m_temp.second = w;
                        m_local.push_back(m_temp);
                        qel.labels = genNewLabels(v,w,lcs);
                        qel.m_local = m_local;
                        Q_tmp.push_back(qel);
                        if ( m_local.size() > m_best.size() ) m_best = m_local; 
                        m_local.pop_back();
                    }
                auto it = std::find(lch.begin(), lch.end(), w);
                lch.erase(it);
                }
            auto it = std::find(lcg.begin(), lcg.end(), v);
            lcg.erase(it);
        }
    

    if( !Q_tmp.empty() ){ return true; } return false;
}


void filter_queue(vector<queue_elem> Q){
    int pos = 0;
    int max_Q_GPU = 0;
    for (queue_elem elem : Q){
            
        vector<LabelClass> lcs = elem.labels;

        std::vector<LabelClass*> label_class_pointers;

        label_class_pointers.reserve(lcs.size());
        
        for (LabelClass& item : lcs) {label_class_pointers.push_back(&item);}
        vector<pair<int,int >> m_local = elem.m_local;
        

        LabelClass *lcc = select_label(label_class_pointers, m_local.size());

        if ( m_local.size() + calc_bound(lcs) <= m_best.size() || ( !lcc && !m_local.empty() )  ){
            
        }else{
           
            if(m_local.size() <= LIMIT_DEPTH ){ 
                
                Q_gpu.push_back(elem); 
                Q.erase(next(Q.begin(), pos));
                if(Q_gpu.size() == 16){
                    vector<queue_elem> Q_tmp;
                    for(int i = 0; i < Q_gpu.size(); i++){
                        sortLabels(Q_gpu[i].labels);
                        Q_tmp.push_back(Q_gpu[i]);
                    }
                    cout << "\n-------> Kernel called with Q_gpu size : " << Q_tmp.size() << endl;
                    kernel(Q_tmp);
                    Q_gpu.clear();
                    return;
                }    
            }
            else{
                Q_cpu.push_back(elem);
                Q.erase(next(Q.begin(), pos));
            } 
        }
    
        pos++;
    }
}



vector<pair<int,int>> gpu_mc_split(const std::vector<std::vector<float>>& g00, const std::vector<std::vector<float>>& g11,
                                          const std::vector<std::string>& l00, const std::vector<std::string>& l11,
                                          std::vector<std::vector<int> >& ring_classes){

    m_best.clear();
    Q.clear();
    Q_cpu.clear();
    Q_gpu.clear();
    first_bb = false;
    max_depth = 0;
    g0 = g00;
    g1 = g11;
    l0 = l00;
    l1 = l11;
    edge_labels = gen_bond_labels(g0, g1);
    int min = std::min(l00.size(), l11.size());
    std::vector<LabelClass> initial_label_classes = gen_initial_labels(l00, l11, ring_classes);
    size_initial_label_classes = initial_label_classes.size();
    if(state_initialized){
        initialized(l00,l11,initial_label_classes);
        vector<pair<int,int>> null_map; 
        return null_map;
    } 

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
            v = select_vertex(lcg,g0);
            vector<int> lch = lc.h;
             while( !lch.empty()){
                w = select_vertex(lch, g1);
                if( matchable(v,w,lc ) )
                    {
                        m_local.at(0).first = v;
                        m_local.at(0).second = w;
                        elem.labels = genNewLabels(v,w,initial_label_classes);  
                        elem.m_local=m_local;
                        Q.push_back(elem);
                    }
                auto it = std::find(lch.begin(), lch.end(), w);
                lch.erase(it);

             }
            auto it = std::find(lcg.begin(), lcg.end(), v);
            lcg.erase(it);
            
        }
    
    bool flag;
    do{
        flag = solve_mcs();
        if(first_bb){
            filter_queue(Q);
            flag = solve_mcs();
        }
     }while(flag);
    

   if(Q_gpu.size() > 0) {
        flag = true;
        do
        {
            flag = solve_2(Q_gpu);
        } while (flag);
   }


    if(Q_cpu.size() > 0){
        flag = true;
        do
        {
            flag = solve_2(Q_cpu);
        } while (flag);
    }

    vector<pair<int,int>> best = m_best;


    m_best.clear();
    
    return best;
}
