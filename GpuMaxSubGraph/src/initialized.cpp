#include "test.hpp"


void initialized(const std::vector<std::string>& l00 ,const std::vector<std::string>& l11 , std::vector<LabelClass> initial_label_classes){
    if(l00.size() > max_l0_size) max_l0_size = l00.size();
    if(l11.size() > max_l1_size) max_l1_size = l11.size();
    if(initial_label_classes.size() > max_initial_label_size) max_initial_label_size = initial_label_classes.size();

    for( LabelClass lc : initial_label_classes ){
        if(lc.g.size() > max_first_len_initialized) max_first_len_initialized = lc.g.size();
    }
}