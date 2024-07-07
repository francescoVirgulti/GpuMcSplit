
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



__device__ __managed__ float *gpu_edge_labels;
__device__ __managed__ int size_edge_labels;


__device__ __managed__ float **gpu_g0;
__device__ __managed__ int size_gpu_g0_row;
__device__ __managed__ int size_gpu_g0_col;

__device__ __managed__ float **gpu_g1;
__device__ __managed__ int size_gpu_g1_row;
__device__ __managed__ int size_gpu_g1_col;





extern bool malloc_done;










