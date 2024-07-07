//
// Created by davide on 5/3/24.
//


#include "gpu_header.hpp"

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

void printLabelClass(LabelClass lb) {
    if( true) {
        cout<< lb.label << " [ ";
        cout<< " G("<< lb.g.size() << "): ";
        if(!lb.g.empty()) {for ( int i : lb.g ) cout<<"["<<i<<"]";}
        cout<< " H("<< lb.h.size() << "): ";
        if(!lb.h.empty()) {for ( int i : lb.h ) cout<<"["<<i<<"]";}
        cout<< " RINGS("<< lb.rings_g.size() << "): [";
        for( vector<int> i : lb.rings_g ){cout<<"("<<i.size()<<")"<<"["; for( int j: i ) cout<<j<<", ";  cout<<" ]";}
        cout<<"]";
        cout<< " edge : " <<lb.adj<<" " ;
        cout<< lb.label << " ] "<<endl;
    }
}



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
    vector<LabelClass> l_draft;
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
                LabelClass tmp(v_conn,w_conn,v_c_rings,adj, label.label);
                l_draft.push_back(tmp);
            }
        }
    }
    return l_draft;
}


int iterazione = 0;

bool solve_mcs() {
    
    queue_elem elem =  Q.back();  

    Q.pop_back();   

    vector<LabelClass> lcs = elem.labels;

    std::vector<LabelClass*> label_class_pointers;

    label_class_pointers.reserve(lcs.size());
    
    for (LabelClass& item : lcs) {label_class_pointers.push_back(&item);}

    vector<pair<int,int >> m_local = elem.m_local;

    
    // first brench&bound 
    if ( m_local.size() >= max_depth) max_depth = m_local.size();
    else if( m_local.size() < max_depth){ 
        Q.push_back(elem);
        first_bb = true; 
        return false;
    }



    LabelClass *lcc = select_label(label_class_pointers, m_local.size());

    if ( m_local.size() + calc_bound(lcs) <= m_best.size() || ( !lcc && !m_local.empty() )  ){ if( !Q.empty() ){ return true; } return false;}

    queue_elem qel;

    LabelClass lc = *lcc;

    pair<int,int> m_temp;

    for( int v : lc.g )  {
        for ( int w : lc.h ) {
            if ( !matchable(v,w,lc) ) continue;
            m_temp.first = v;
            m_temp.second = w;
            m_local.push_back(m_temp);
            qel.labels = genNewLabels(v,w,lcs);
            qel.m_local = m_local;
            Q.push_back(qel);
            if ( m_local.size() > m_best.size() ) m_best = m_local;
            
            
            if(m_local.size() <= LIMIT_DEPTH ){ 
                if(Q.size() == 16) {
                    for(int i = 0; i < Q.size(); i++){
                        sortLabels(Q[i].labels);
                    }
                    cout << "\n Q_gpu size : " << Q.size() << endl;
                    kernel(Q);
                    Q.clear();
                } 
            }
            else Q_cpu.push_back(elem);

            m_local.pop_back();
        }
    }
    if( !Q.empty() ){ return true; } return false;
}


bool solve_2_mcs() {
    queue_elem elem =  Q_cpu.back();  

    Q_cpu.pop_back();   

    vector<LabelClass> lcs = elem.labels;

    std::vector<LabelClass*> label_class_pointers;

    label_class_pointers.reserve(lcs.size());
    
    for (LabelClass& item : lcs) {label_class_pointers.push_back(&item);}

    vector<pair<int,int >> m_local = elem.m_local;

    LabelClass *lcc = select_label(label_class_pointers, m_local.size());

    if ( m_local.size() + calc_bound(lcs) <= m_best.size() || ( !lcc && !m_local.empty() )  ){ if( !Q_cpu.empty() ){ return true; } return false;}

    queue_elem qel;

    LabelClass lc = *lcc;

    pair<int,int> m_temp;

    for( int v : lc.g )  {
        for ( int w : lc.h ) {
            if ( !matchable(v,w,lc) ) continue;
            m_temp.first = v;
            m_temp.second = w;
            m_local.push_back(m_temp);
            qel.labels = genNewLabels(v,w,lcs);
            qel.m_local = m_local;
            Q_cpu.push_back(qel);
            if ( m_local.size() > m_best.size() ) m_best = m_local;
            m_local.pop_back();
        }
    }
    if( !Q_cpu.empty() ){ return true; } return false;
}


int calcSize(vector<LabelClass> lcs) {
    int result=0;
    for( LabelClass lc : lcs) {
        result = result + (lc.g.size()+lc.h.size());
        for( vector<int> i : lc.rings_g ) {
            result += i.size();
        }
        result += 3;
    }
    return result;
}



void filter_queue(vector<queue_elem> Q){
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
            }
            else Q_cpu.push_back(elem);
        }
    

    }
}



vector<pair<int,int>> gpu_mc_split(const std::vector<std::vector<float>>& g00, const std::vector<std::vector<float>>& g11,
                                          const std::vector<std::string>& l00, const std::vector<std::string>& l11,
                                          std::vector<std::vector<int> >& ring_classes){

    int depth = 1;
    g0 = g00;
    g1 = g11;
    l0 = l00;
    l1 = l11;

    edge_labels = gen_bond_labels(g0, g1);
    int min = std::min(l0.size(), l1.size());
    std::vector<LabelClass> initial_label_classes = gen_initial_labels(l0, l1, ring_classes);

    if(state_initialized){
        initialized(l00,l11,initial_label_classes);
        vector<pair<int,int>> null_map; 
        return null_map;
    } 


    size_initial_label_classes = initial_label_classes.size();
    
    Q.reserve(32*32);
    for(auto & x : Q) {
        x.labels.reserve(size_initial_label_classes);
        x.m_local.reserve(min);
    }

 


    vector<pair<int,int>> m_local={{}};
    queue_elem elem;
    int v,w;
   for( LabelClass lc : initial_label_classes ) {
        v = select_vertex(lc.g,g0);
        w = select_vertex(lc.h,g1);
        if( !matchable(v,w,lc ) ) {
            
                for( int ww : lc.h ){
                    if(matchable(v,ww,lc)){
                        m_local.at(0).first = v;
                        m_local.at(0).second = ww;
                        elem.labels = genNewLabels(v,ww,initial_label_classes);
                        elem.m_local = m_local;
                        Q.push_back(elem);
                    }
                }
            
            
        }else{
            m_local.at(0).first = v;
            m_local.at(0).second = w;
            elem.labels = genNewLabels(v,w,initial_label_classes);  
            elem.m_local=m_local;
            Q.push_back(elem);
        }
    }

    bool flag;
    do{
        flag = solve_mcs();
            if(first_bb){
                cout << "\nm_best_size before kernel : " << m_best.size() << endl;
                filter_queue(Q);

                int num_iter = Q_gpu.size() / 16;
                int remaining = Q_gpu.size() % 16;

                vector<queue_elem> Q_tmp;

                for(int j = 0; j < num_iter; j++){
                    for(int i = j * 16; i < 16 * (j + 1); i++){
                        sortLabels(Q_gpu[i].labels);
                        Q_tmp.push_back(Q_gpu[i]);
                    } 
                    cout << "\n-------> Kernel called with Q_gpu size : " << Q_tmp.size() << endl;
                    kernel(Q_tmp);
                    Q_tmp.clear();
                }

                if(remaining > 0){
                    for(int i = num_iter * 16 ; i < (num_iter * 16) + remaining; i++){
                        sortLabels(Q_gpu[i].labels);
                        Q_tmp.push_back(Q_gpu[i]);
                    } 
                    cout << "\n-------->  Kernel called with Q_gpu size : " << Q_tmp.size() << endl;
                    kernel(Q_tmp);
                    Q_tmp.clear();
                }

                
                cout << "\nm_best_size after kernel : " << m_best.size() << endl;

                cout << "\n Q_cpu size : " << Q_cpu.size() << endl;
                iterazione = 0;
                if(Q_cpu.size() > 0 ){
                    do
                    {
                        flag = solve_2_mcs();
                    } while (flag);
                }
                return m_best;
            }


            if( flag==false && first_bb == false ) {
                cout << "\n Q_cpu size : " << Q_cpu.size() << endl;
                iterazione = 0;
                if(Q_cpu.size() > 0 ){
                    do
                    {
                        flag = solve_2_mcs();
                    } while (flag);
                }
                return m_best;
            }


        }while(flag);
    return m_best;
}
