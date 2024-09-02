//#include <string>
#include <vector>
#include "Label.cpp"

#include <rdkit/GraphMol/GraphMol.h>
#include <GraphMol/FileParsers/MolSupplier.h>

#include <Numerics/Matrix.h>


using namespace std;
using namespace RDKit;

extern int algorithm_used;

extern double malloc_elapsed_seconds;
extern int max_l0_size ;
extern int max_l1_size ; 
extern int max_first_len_initialized ;

extern int edge_label_size;
extern int max_initial_label_size ;
extern bool state_initialized;

extern bool malloc_done;

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

extern float *main_gpu_edge_labels;
extern float **main_gpu_g0;
extern float **main_gpu_g1;

std::vector<std::vector<int> > gen_rings_classes(
     std::vector<std::string> l0,
     std::vector<std::string> l1,
     std::vector<std::vector<int> > ring_info_m0,
     std::vector<std::vector<int> > ring_info_m1
);



std::vector<std::pair<std::string, int> > gen_rotations(const std::string& s);

int select_vertex(std::vector<int>& vtx_set, std::vector<std::vector<float> >& g);


std::vector<int> hood(int vtx, const std::vector<std::vector<float>>& g, float edge);

std::pair<int, int> pair_vertex(LabelClass lb,  std::vector<std::vector<float>> g0);

std::vector<std::pair<int,int>> gpu_mc_split(const std::vector<std::vector<float>>& g00, const std::vector<std::vector<float>>& g11,
                                          const std::vector<std::string>& l0, const std::vector<std::string>& l1,
                                          std::vector<std::vector<int> >& ring_classes);


int calc_bound( std::vector<LabelClass> label_classes);

std::vector<float> gen_bond_labels(const std::vector<std::vector<float> >& g0, const std::vector<std::vector<float> >& g1);

ROMol smiles_mcs( std::string& smile0,  std::string& smile1, int bond_match = 1, int ring_match = 1);

std::vector<LabelClass> gen_initial_labels(const std::vector<std::string>& l0, const std::vector<std::string>& l1,     std::vector<std::vector<int> >& ring_classes);

std::vector<std::vector<int>> gen_ring_classes(const RDKit::RWMol& mol0, const RDKit::RWMol& mol1);


std::vector<std::pair<int, int>> mc_split(const std::vector<std::vector<float>> g0, const std::vector<std::vector<float>> g1,
                                          const std::vector<std::string>& l0, const std::vector<std::string>& l1,
                                          std::vector<std::vector<int> >& ring_classes);

void search_mcs(std::vector<std::vector<float> > g0, std::vector<std::vector<float> > g1, std::vector<LabelClass>& label_classes, std::vector<float> edge_labels, const std::vector<std::pair<int, int> >& m);

RDKit::ROMol g2mol( std::vector<std::string>& labels,  std::vector<std::vector<float>>& adj) ;

ROMol mol_mcs(const RDKit::RWMol &mol0, const RDKit::RWMol &mol1, int bond_match=1, int ring_match=1, int return_map=0) ;

void initialized(const std::vector<std::string>& l00 ,const std::vector<std::string>& l11 , std::vector<LabelClass> initial_label_classes);

vector<pair<int,int>> mcs_iterative(const std::vector<std::vector<float>>& g00, const std::vector<std::vector<float>>& g11,
                                          const std::vector<std::string>& l0, const std::vector<std::string>& l1,
                                          std::vector<std::vector<int> >& ring_classes);