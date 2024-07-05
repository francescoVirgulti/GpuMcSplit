#include <vector>
#include "test.hpp"

#include <string>
#include <algorithm>
using namespace std;



extern std::vector<std::vector<float>> g0;
extern std::vector<std::vector<float>> g1;
extern std::vector<std::string> l0;
extern std::vector<std::string> l1;
extern int size_initial_label_classes ;
extern std::vector<float> edge_labels;
extern vector<pair<int,int>> m_best;

extern int max_first_len;
extern int max_second_len;

struct queue_elem{
    vector<LabelClass> labels;
    vector<pair<int,int>> m_local;
};



//struct 
typedef struct{
    int g_size;
    int h_size;
    int row_ring_size;
    int *col_ring_size;
    int *g;
    int *h;
    int adj;
    char label[4];
    int **rings_g;
}GpuLabelClass;


typedef struct{
    int first;
    int second;
}Pair;

extern Pair *m_best_solution;


typedef struct {
    int labels_size;
    int m_size;
    GpuLabelClass *labels;
    GpuLabelClass single_label;
    int *idxList;
    Pair *m_local;
}ThreadVar;

//auto is for autonomous

extern ThreadVar **thread_pool_list;


extern int *auto_pool_size;
extern Pair **auto_pool_m_best;
extern int *auto_pool_len_m_best;
extern ThreadVar *auto_pool_tmp;

extern vector<int> length_list;







extern bool malloc_done;













void kernel( vector<queue_elem> Q_gpu  );