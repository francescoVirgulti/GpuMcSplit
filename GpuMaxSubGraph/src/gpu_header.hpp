#include <vector>
#include "test.hpp"

#include <string>
#include <algorithm>
using namespace std;



extern std::vector<std::vector<float>> g0;
extern std::vector<std::vector<float>> g1;
extern std::vector<float> edge_labels;
extern vector<pair<int,int>> m_best;


struct queue_elem{
    vector<LabelClass> labels;
    vector<pair<int,int>> m_local;
};



vector<pair<int,int>> first_solution(
                                const std::vector<std::string>& l0,
                                const std::vector<std::string>& l1,
                                vector<queue_elem> Q_filter,
                                int size_initial_label_classes ) ;