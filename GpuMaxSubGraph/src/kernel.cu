//
// Created by davide on 5/3/24.
//

#include "gpu_header.hpp"
#include "cuda_header.h"
#include <vector>
#include <string.h>
#include <string>
#include <queue>
#include <algorithm>
#include <unordered_set>
#include <chrono>  // Include la libreria chrono
using namespace std;

double malloc_elapsed_seconds;




// Initializes GPU pointer variables with corresponding main GPU pointers.
void pointer_initialized(){
    gpu_edge_labels = main_gpu_edge_labels;
    gpu_g0 = main_gpu_g0;
    gpu_g1 = main_gpu_g1;
}

__device__ void copyIntArray(int *a, int *b, int sizeb){
    for ( int i = 0 ; i < sizeb ; i++){
        a[i] = b[i];
    }
}

__device__ void copyIntMatrix(int **a, int **b, int rowsize, int *colsize )
{

    for( int i = 0 ; i < rowsize ; i++){
        for( int j = 0 ; j < colsize[i] ; j++){
            a[i][j] = b[i][j];
        }
    }

}

// Copies the contents of one GpuLabelClass object to another on the GPU, including arrays and matrix data.
__device__ void cpyGpuLabelClass(GpuLabelClass *l1, GpuLabelClass l2){
    l1->adj = l2.adj;
    l1->row_ring_size = l2.row_ring_size;
    l1->g_size = l2.g_size;
    l1->h_size = l2.h_size;
    for(int i = 0; i < 4 ; i++) {
        l1->label[i] = l2.label[i];
    }
    copyIntArray( l1->g , l2.g, l2.g_size);
    copyIntArray( l1->h, l2.h, l2.h_size);
    copyIntArray( l1->col_ring_size, l2.col_ring_size , l2.row_ring_size);
    copyIntMatrix( l1->rings_g, l2.rings_g, l2.row_ring_size, l2.col_ring_size );/**/
}

// Copies data from a vector of edge labels to a float pointer array, or sets the pointer to nullptr if the vector is empty.
void vectorToPointerEdge(float *gpu_edge_labels){
    if(edge_labels.size() == 0){
        gpu_edge_labels = nullptr;
        return;
    }
    int size = 0;
    for(float edg : edge_labels){
        gpu_edge_labels[size] = edg;
        size++;
    }

    return;
}

// Converts a 2D vector to a 2D pointer matrix, copying data from the vector to the allocated matrix.
void vectorToPointerMatrix(const std::vector<std::vector<float>>& g,float** gpu_g) {
    // Get dimensions of the vector
    int numRows = g.size();
    if (numRows == 0) {
        // Empty vector, set pointers to nullptr
        gpu_g = nullptr;
        return;
    }

    int numCol = g[0].size();

    for (int i = 0; i < numRows; ++i) {
        for (int j = 0; j < numCol; ++j) {
            gpu_g[i][j] = g[i][j];
        }
    }
}

// Converts data from a vector of LabelClass objects (CPU) to an array of GpuLabelClass objects (GPU).
void LabelFromCpuToGpu(GpuLabelClass *new_label, const vector<LabelClass>& old_label ){

    for (int idx = 0 ; idx < old_label.size() ; ++idx ){
        new_label[idx].g_size = old_label.at(idx).g.size();
        new_label[idx].h_size = old_label.at(idx).h.size();
        new_label[idx].row_ring_size = old_label.at(idx).rings_g.size();
        new_label[idx].adj = old_label.at(idx).adj;
        strcpy(new_label[idx].label, old_label.at(idx).label.c_str() );

        for ( int j = 0 ; j < old_label.at(idx).g.size() ; ++j ){
            new_label[idx].g[j] = old_label.at(idx).g.at(j);
        }
        for ( int j = 0 ; j < old_label.at(idx).h.size() ; ++j ){
            new_label[idx].h[j] = old_label.at(idx).h.at(j);}

        for ( int row = 0 ; row < old_label.at(idx).rings_g.size() ; ++row ){
            new_label[idx].col_ring_size[row] = old_label.at(idx).rings_g.at(row).size();
            for ( int col = 0 ; col < old_label.at(idx).rings_g.at(row).size() ; ++col ){
                new_label[idx].rings_g[row][col] = old_label.at(idx).rings_g.at(row).at(col);
            }
        }

    }
}

// vtx_set: selected label class
// g: selected graph
__device__ bool device_contains(int value, int *arr, int size) {
    for (int i = 0; i < size; ++i) {
        if (arr[i] == value) {
            return true;
        }
    }
    return false;
}

//      puts into a 2D array the data regarding indexes of rings related to the array of elements
//      2D-array that will contain the result that will be modified
//      1D-array containing idxList
//      1D-array of elements
//      int size of elements
__device__ int device_get_ring_match_data(int *dim_col, int **result, int *idxList, int *elems, int elem_size, GpuLabelClass *lc){
    int index;
    int idx_list_size= 0 ;


    for( int i = 0; i < elem_size ; ++i){
        index = 0;
        for( int j = 0 ; j < lc->g_size ; ++j ){
            if( lc->g[j] == elems[i] ) {idxList[idx_list_size] = index; idx_list_size++;}
            index++;
        }
    }


    for ( int i = 0  ; i < idx_list_size ; ++i  ){
        dim_col[i] = lc->col_ring_size[idxList[i]];
        for ( int j = 0 ; j < dim_col[i] ; ++j ){
            result[i][j] = lc->rings_g[idxList[i]][j];
        }
    }



    return idx_list_size;
}

// return the best select label given an array of labels
__device__ void device_select_label(GpuLabelClass *label , GpuLabelClass *lcs, int map_size, int lcs_size){
    int min = 999;
    int max ;
    for( int i = 0 ; i < lcs_size ; ++i ){
        //printf("LABEL CLASSES INTERNE[%d]\n", i);
        //printLabelClass(lcs[i]);
        if( lcs[i].adj == 1 || map_size == 0 ){
            if( lcs[i].g_size > lcs[i].h_size ) max = lcs[i].g_size;
            else max = lcs[i].h_size;
            //printf("\nMAX : %d\n", max);
            if( max < min ){
                
                min = max;
                cpyGpuLabelClass(label, lcs[i] );
            }
        }
    }
    return;
}

// compute the bound given a 1D array of struct GpuLabelClass and its size
__device__ int device_calc_bound(GpuLabelClass *lcs, int lc_size) {
    int bound = 0;
    for( int i = 0 ; i < lc_size ; ++i){
        if ( lcs[i].g_size > lcs[i].h_size ) bound = bound + lcs[i].h_size;
        else bound = bound + lcs[i].g_size;
    }
    return bound;
}

//return = size of the friends
//friend is the OUTPUT
__device__ int device_hoodG(int *friends,int vtx, float edge, float **g, int size_g) {
    int size = 0;
    
    for (int i = 0; i < size_g; i++) {
        if ( g[i][vtx] == edge && vtx != i) {
            friends[size] = i;
            size++;
        }
    }

    return size;
}

// result == size of generated label
// output is l_draft
// input : v
__device__ void device_resize(int *array, int size_arr, int place_availabel){
    int count = 0;

    for( int i = 0 ; i < size_arr && place_availabel > 0 ; i++){
        if( array[i] != -1 ){
            place_availabel--;
            array[count] = array[i];
            count ++;
        }
    }
}

__device__ int device_gen_new_labels(GpuLabelClass *l_draft ,  int v, int w, GpuLabelClass *lcs, int lcs_size, int *idxList) {
    int vs,ws, draft_size = 0;
    int dim_row;
    //int count = 0;
    for ( int i = 0 ; i < lcs_size ; ++i ){
      // printf("\ndevice_gen_new_labels iterazione num : %d", i );
        for ( int j = 0 ; j < size_edge_labels ; ++j ){
            int friendsize;
            friendsize = device_hoodG( l_draft[draft_size].g , v , gpu_edge_labels[j], gpu_g0 , size_gpu_g0_row);
        
            vs = 0;
            for ( int k = 0; k < friendsize ; ++k ){
                if( device_contains(l_draft[draft_size].g[k] , lcs[i].g, lcs[i].g_size) ){ vs++;  }
                else{ l_draft[draft_size].g[k] = -1;}
            }
        
            device_resize(l_draft[draft_size].g, friendsize, vs );

            dim_row = device_get_ring_match_data(l_draft[draft_size].col_ring_size, l_draft[draft_size].rings_g, idxList ,l_draft[draft_size].g, vs, &lcs[i] );
            

            friendsize = device_hoodG(l_draft[draft_size].h, w, gpu_edge_labels[j], gpu_g1, size_gpu_g1_row );
            //printf("\n esco da hood 2");
            ws = 0;
            for ( int k = 0 ; k < friendsize ; ++k ){
                
                if( device_contains(l_draft[draft_size].h[k], lcs[i].h, lcs[i].h_size) ){  ws++; }
                else {
                    l_draft[draft_size].h[k] = -1;
                }
            }
            device_resize(l_draft[draft_size].h, friendsize, ws );
    
            int adj;
            if ( ws > 0 && vs > 0 ){
                if( gpu_edge_labels[j] != 0.0 || lcs[i].adj == 1 ) {adj = 1;}
                else { adj = 0; }

                l_draft[draft_size].g_size = vs;
                l_draft[draft_size].h_size = ws;
                l_draft[draft_size].row_ring_size = dim_row;
                l_draft[draft_size].adj = adj;
                for( int c = 0 ; c < 4 ; c++){
                    l_draft[draft_size].label[c] = lcs[i].label[c];
                }
                draft_size++;
            }
        }
    }
    return draft_size;
}

//given two atoms from the same label, return true if they are matchable, false otherwise
// based on how their rings matches
__device__ bool device_matchable(int **v_ring_atoms, int v, int w, GpuLabelClass *lc, int *idxList) {

    device_get_ring_match_data(lc->col_ring_size, v_ring_atoms, idxList , &v, 1 ,lc);
    if( lc->col_ring_size[idxList[0]] > 0  ){
        for(int i = 0; i < lc->col_ring_size[idxList[0]]  ; i++){
            if( v_ring_atoms[0][i] == -1 )return false;
            if( v_ring_atoms[0][i] == w ) return true;
        }
        return false;
    }
    return true;
}

// vtx_set: selected label class
// g: selected graph
__device__ int device_select_vertex(int result, int *result_pos, int *vtx_set, int vtx_size, float **g, int num_row, int num_column) {
    int max_deg = -1;
    
    //printf("vtx_size %d  i:   %d    nr: %d   nc:  %d ", vtx_size, i, num_row, num_column);
    for(int i = 0; i < vtx_size; i++){
       
        int deg = 0;
        for(int j = 0; j < num_column; j++){
            int consider = g[vtx_set[i]][j];
            if(consider != 0){
                deg++;
            }
        }

        if(deg>max_deg){
            max_deg = deg;
            result = vtx_set[i];
            *result_pos = i;
        }
    }
    return result;
}

//copy var1 in var2
__device__ void copy_single_ThreadVar(ThreadVar *var2, ThreadVar var1){

    var2->labels_size = var1.labels_size;
    var2->m_size = var1.m_size;

    for (int i = 0; i < var1.labels_size ; i++){
         var2->labels[i].adj =  var1.labels[i].adj;
        var2->labels[i].g_size = var1.labels[i].g_size ;
        var2->labels[i].h_size = var1.labels[i].h_size ;
        var2->labels[i].row_ring_size = var1.labels[i].row_ring_size ; 
        copyIntArray(var2->labels[i].g, var1.labels[i].g,var1.labels[i].g_size );
        copyIntArray(var2->labels[i].h, var1.labels[i].h,var1.labels[i].h_size );
        copyIntArray(var2->labels[i].col_ring_size, var1.labels[i].col_ring_size,var1.labels[i].row_ring_size );

        for(int c = 0; c < 4 ; c++) {
            var2->labels[i].label[c] = var1.labels[i].label[c] ;
        }

        copyIntMatrix(var2->labels[i].rings_g,
                      var1.labels[i].rings_g, 
                      var1.labels[i].row_ring_size,
                      var1.labels[i].col_ring_size ) ;
    }

    for(int i = 0; i < var1.m_size; i++) {
        var2->m_local[i].first =  var1.m_local[i].first;
        var2->m_local[i].second =  var1.m_local[i].second;
    }
}

__global__ void autonomouslySolve(ThreadVar **thread_pool_list, int *queue_size_list, int* m_best_size,Pair **auto_pool_m_best_list, ThreadVar *tmp, int Q_size ){
    
    int globalIdx = blockIdx.x * blockDim.x + threadIdx.x;
    if(globalIdx >= Q_size) return;

    ThreadVar *thread_pool = thread_pool_list[globalIdx];
    int queue_size = queue_size_list[globalIdx];
    int m_best_size_index = globalIdx;
    Pair *auto_pool_m_best = auto_pool_m_best_list[globalIdx];
    tmp = tmp + globalIdx;



    int flag = 0;
    int z ;
    int iterazione;
    int index;

    int max_legth_queue = 1;


    while(queue_size > 0) {
        
        if(queue_size > max_legth_queue) {
            max_legth_queue = queue_size;
        }


        queue_size -- ;
        
        copy_single_ThreadVar(tmp , thread_pool[queue_size] );
        ThreadVar TMP = *tmp;
        GpuLabelClass *label = &TMP.single_label;
        device_select_label(label, TMP.labels, TMP.m_size, TMP.labels_size);



        index = 0;
        flag = 0;
        int v_real_size = label->g_size;
        int w_real_size = label->h_size;
        int v;
        int w; 
        int pos_v;
        int pos_w;
        int tmp;

        if((TMP.m_size + device_calc_bound(TMP.labels, TMP.labels_size) < m_best_size[m_best_size_index]) || !label) {flag = 1;}

        if( flag == 0 ){
            
            index = queue_size;
            while( v_real_size > 0 ){
                v = device_select_vertex(v , &pos_v, label->g, v_real_size, gpu_g0, size_gpu_g0_row, size_gpu_g0_col );
                v_real_size = v_real_size -1;
                tmp = label->g[v_real_size];
                label->g[v_real_size] = label->g[pos_v];
                label->g[pos_v] = tmp;
                w_real_size = label->h_size;
                while( w_real_size > 0 ){
                    w = device_select_vertex(w , &pos_w, label->h, w_real_size, gpu_g1, size_gpu_g1_row, size_gpu_g1_col );
                    w_real_size = w_real_size-1;
                    tmp = label->h[w_real_size];
                    label->h[w_real_size] = label->h[pos_w];
                    label->h[pos_w] = tmp;
                    if( !device_matchable(label->rings_g, v, w, label, TMP.idxList ) ) continue;
                    

                    for(z = 0; z < TMP.m_size; z ++){
                        thread_pool[index].m_local[z].first = TMP.m_local[z].first;
                        thread_pool[index].m_local[z].second = TMP.m_local[z].second;
                    } 
                    thread_pool[index].m_size = TMP.m_size +1;
                    thread_pool[index].m_local[z].first = v;
                    thread_pool[index].m_local[z].second =w;
                    
                    int l_s = device_gen_new_labels( 
                        thread_pool[index].labels, 
                        v, 
                        w , 
                        TMP.labels,  
                        TMP.labels_size, 
                        TMP.idxList );


                    thread_pool[index].labels_size = l_s;

                    if(thread_pool[index].m_size > m_best_size[m_best_size_index] ){
                        m_best_size[m_best_size_index] = thread_pool[index].m_size;
                        for( int z = 0 ; z < thread_pool[index].m_size; z++ ){
                            auto_pool_m_best[z].first = thread_pool[index].m_local[z].first;
                            auto_pool_m_best[z].second = thread_pool[index].m_local[z].second;
                        }
                    }
                    
                    index ++;
                }
            }
        }


        queue_size = index;
        iterazione ++;
    }

   return ;
}

void checkError(int iterazione, int line, cudaError_t r) {
    if (r != cudaSuccess) {
        printf("CUDA error on line %d - iterazione %d : %s\n", line, iterazione,  cudaGetErrorString(r));
        exit(0);
    }
}

void malloc(){
    size_edge_labels = edge_label_size;

    int Q_GPU_size = 16;
    int min_mol_size;

    if(max_l0_size > max_l1_size) min_mol_size = max_l1_size;
    else  min_mol_size = max_l0_size;
    
    /*
    cout << "size_edge_labels : " << size_edge_labels;
    cout << "min_mol_size : " << min_mol_size <<endl;
    cout << "max_l0_size : " << max_l0_size <<endl;
    */

    clock_t start = clock();                          
    //cuda Mallocs
    cudaMallocManaged(&m_best_solution , sizeof(Pair)* max_l1_size);
    //cuda malloc edge labels
    checkError(0, __LINE__ , cudaMallocManaged(&main_gpu_edge_labels, sizeof(float) * size_edge_labels));
    //cuda malloc adj matrix mol 0
    checkError(0, __LINE__ , cudaMallocManaged((void **) &main_gpu_g0, max_l0_size * sizeof(float *)));
    for (int i = 0; i < max_l0_size; ++i) { checkError(0, __LINE__ , cudaMallocManaged((void **) &(main_gpu_g0[i]), max_l0_size * sizeof(float))); }
    //cuda malloc adj matrix mol 1
    checkError(0, __LINE__ , cudaMallocManaged((void **) &main_gpu_g1, max_l1_size * sizeof(float *)));
    for (int i = 0; i < max_l1_size; ++i) { checkError(0, __LINE__ , cudaMallocManaged((void **) &(main_gpu_g1[i]), max_l1_size * sizeof(float))); }


    //lista di code per ogni singolo thread
    checkError(0, __LINE__ , cudaMallocManaged((void **) &thread_pool_list, (Q_GPU_size)  * sizeof(ThreadVar *) ));
    checkError(0, __LINE__ , cudaMallocManaged(&auto_pool_size, (Q_GPU_size) * sizeof(int *)));
    checkError(0, __LINE__ , cudaMallocManaged(&auto_pool_len_m_best, (Q_GPU_size) * sizeof(int *)));
    checkError(0, __LINE__ , cudaMallocManaged((void **) &auto_pool_m_best, (Q_GPU_size)  * sizeof(Pair *) ));
    checkError(0, __LINE__ , cudaMallocManaged( &auto_pool_tmp, (Q_GPU_size)  * sizeof(ThreadVar ) ));

    int common_queue_element_size = 3;
    
    for(int i = 0; i < (Q_GPU_size); i++) {
        // Esempi di allocazioni, assicurati di gestire gli errori per ciascuna
            checkError(i, __LINE__, cudaMallocManaged(&(auto_pool_tmp[i].single_label.col_ring_size), sizeof(int) * min_mol_size));
            checkError(i, __LINE__, cudaMallocManaged(&(auto_pool_tmp[i].single_label.g), sizeof(int) * max_first_len_initialized));
            checkError(i, __LINE__, cudaMallocManaged(&(auto_pool_tmp[i].single_label.h), sizeof(int) * max_first_len_initialized));
            
            // Allocazione per array di puntatori e iterazione su di essi
            checkError(i, __LINE__, cudaMallocManaged(&(auto_pool_tmp[i].single_label.rings_g), sizeof(int *) * max_first_len_initialized));
            for (int h = 0; h < max_l0_size; ++h) {
                checkError(i, __LINE__, cudaMallocManaged(&(auto_pool_tmp[i].single_label.rings_g[h]), sizeof(int) * max_first_len_initialized));
            }

            // Allocazioni successive
            checkError(i, __LINE__, cudaMallocManaged(&(auto_pool_m_best[i]), sizeof(Pair) * min_mol_size));
            checkError(i, __LINE__, cudaMallocManaged(&(auto_pool_tmp[i].m_local), sizeof(Pair) * min_mol_size));
            checkError(i, __LINE__, cudaMallocManaged(&(auto_pool_tmp[i].idxList), sizeof(int) * min_mol_size));


            checkError(i, __LINE__, cudaMallocManaged(&(auto_pool_tmp[i].labels), ( max_initial_label_size + (2 * size_edge_labels)) * sizeof(GpuLabelClass)));

            //first + edge_labels
            for(int s = 0; s < size_edge_labels ; s++){
                checkError(i, __LINE__, cudaMallocManaged(&(auto_pool_tmp[i].labels[s].col_ring_size), sizeof(int) * max_first_len_initialized));
                checkError(i, __LINE__, cudaMallocManaged(&(auto_pool_tmp[i].labels[s].g), sizeof(int) * max_first_len_initialized));
                checkError(i, __LINE__, cudaMallocManaged(&(auto_pool_tmp[i].labels[s].h), sizeof(int) * max_first_len_initialized));

                checkError(i, __LINE__, cudaMallocManaged(&(auto_pool_tmp[i].labels[s].rings_g), sizeof(int *) * max_first_len_initialized));
                for (int h = 0; h < max_first_len_initialized ; ++h) {
                    checkError(i, __LINE__, cudaMallocManaged(&(auto_pool_tmp[i].labels[s].rings_g[h]), sizeof(int) * max_first_len_initialized));
                }
            }

            //second + edge_labels
            for(int s = size_edge_labels; s < (2 * size_edge_labels); s++ ){
                checkError(i, __LINE__, cudaMallocManaged(&(auto_pool_tmp[i].labels[s].col_ring_size), sizeof(int) * max_first_len_initialized));
                checkError(i, __LINE__, cudaMallocManaged(&(auto_pool_tmp[i].labels[s].g), sizeof(int) * max_first_len_initialized));
                checkError(i, __LINE__, cudaMallocManaged(&(auto_pool_tmp[i].labels[s].h), sizeof(int) * max_first_len_initialized));

                checkError(i, __LINE__, cudaMallocManaged(&(auto_pool_tmp[i].labels[s].rings_g), sizeof(int *) * max_first_len_initialized));
                for (int h = 0; h < max_first_len_initialized ; ++h) {
                    checkError(i, __LINE__, cudaMallocManaged(&(auto_pool_tmp[i].labels[s].rings_g[h]), sizeof(int) * max_first_len_initialized));
                }
            }


            // 5 is an indicative number

            for (int k = (2 * size_edge_labels); k < max_initial_label_size + (2 * size_edge_labels); ++k) {
                checkError(i, __LINE__, cudaMallocManaged(&(auto_pool_tmp[i].labels[k].col_ring_size), sizeof(int) * common_queue_element_size));
                checkError(i, __LINE__, cudaMallocManaged(&(auto_pool_tmp[i].labels[k].g), sizeof(int) * common_queue_element_size));
                checkError(i, __LINE__, cudaMallocManaged(&(auto_pool_tmp[i].labels[k].h), sizeof(int) * common_queue_element_size));

                checkError(i, __LINE__, cudaMallocManaged(&(auto_pool_tmp[i].labels[k].rings_g), sizeof(int *) * common_queue_element_size));
                for (int h = 0; h < 5; ++h) {
                    checkError(i, __LINE__, cudaMallocManaged(&(auto_pool_tmp[i].labels[k].rings_g[h]), sizeof(int) * common_queue_element_size));
                }
            }

    }


    //bool flag = false;
    
    for (int f = 0; f <  Q_GPU_size ; f++) {
        int length = (min_mol_size) / 2  ;
        length_list.push_back(length);


        checkError(f, __LINE__, cudaMallocManaged((void **)&thread_pool_list[f], sizeof(ThreadVar) * length));

        // Allocazione per le strutture dati all'interno di ciascun ThreadVar
        // prova
        //length = 15;

        for (int j = 0; j < length; ++j) {
            // Esempi di allocazioni, assicurati di gestire gli errori per ciascuna
            checkError(f, __LINE__, cudaMallocManaged((void **)&(thread_pool_list[f][j].single_label.col_ring_size), sizeof(int) * max_first_len_initialized));
            checkError(f, __LINE__, cudaMallocManaged((void **)&(thread_pool_list[f][j].single_label.g), sizeof(int) * max_first_len_initialized));
            checkError(f, __LINE__, cudaMallocManaged((void **)&(thread_pool_list[f][j].single_label.h), sizeof(int) * max_first_len_initialized));
            
            // Allocazione per array di puntatori e iterazione su di essi
            checkError(f, __LINE__, cudaMallocManaged((void **)&(thread_pool_list[f][j].single_label.rings_g), sizeof(int *) * max_first_len_initialized));
            for (int h = 0; h < max_first_len_initialized; ++h) {
                checkError(f, __LINE__, cudaMallocManaged((void **)&(thread_pool_list[f][j].single_label.rings_g[h]), sizeof(int) * max_first_len_initialized));
            }

            // Allocazioni successive
            checkError(f, __LINE__, cudaMallocManaged((void **)&(auto_pool_m_best[f]), sizeof(Pair) * min_mol_size));
            checkError(f, __LINE__, cudaMallocManaged((void **)&(thread_pool_list[f][j].m_local), sizeof(Pair) * min_mol_size));
            checkError(f, __LINE__, cudaMallocManaged((void **)&(thread_pool_list[f][j].idxList), sizeof(int) * min_mol_size));


            checkError(f, __LINE__, cudaMallocManaged((void **)&(thread_pool_list[f][j].labels), (max_initial_label_size  + (2 * size_edge_labels) ) * sizeof(GpuLabelClass)));
            
            //first + edge_label
            for(int s = 0; s < size_edge_labels ; s++){
                checkError(f, __LINE__, cudaMallocManaged((void **)&(thread_pool_list[f][j].labels[s].col_ring_size), sizeof(int) * max_first_len_initialized));
                checkError(f, __LINE__, cudaMallocManaged((void **)&(thread_pool_list[f][j].labels[s].g), sizeof(int) * max_first_len_initialized));
                checkError(f, __LINE__, cudaMallocManaged((void **)&(thread_pool_list[f][j].labels[s].h), sizeof(int) * max_first_len_initialized));

                checkError(f, __LINE__, cudaMallocManaged((void **)&(thread_pool_list[f][j].labels[s].rings_g), sizeof(int *) * max_first_len_initialized));
                for (int h = 0; h < max_first_len_initialized; ++h) {
                    checkError(f, __LINE__, cudaMallocManaged((void **)&(thread_pool_list[f][j].labels[s].rings_g[h]), sizeof(int) * max_first_len_initialized));
                }
            }

            //second + edge_labels
            for( int s = size_edge_labels ; s < (2 * size_edge_labels); s++){
                checkError(f, __LINE__, cudaMallocManaged((void **)&(thread_pool_list[f][j].labels[s].col_ring_size), sizeof(int) * max_first_len_initialized));
                checkError(f, __LINE__, cudaMallocManaged((void **)&(thread_pool_list[f][j].labels[s].g), sizeof(int) * max_first_len_initialized));
                checkError(f, __LINE__, cudaMallocManaged((void **)&(thread_pool_list[f][j].labels[s].h), sizeof(int) * max_first_len_initialized));

                checkError(f, __LINE__, cudaMallocManaged((void **)&(thread_pool_list[f][j].labels[s].rings_g), sizeof(int *) * max_first_len_initialized));
                for (int h = 0; h < max_first_len_initialized; ++h) {
                    checkError(f, __LINE__, cudaMallocManaged((void **)&(thread_pool_list[f][j].labels[s].rings_g[h]), sizeof(int) * max_first_len_initialized));
                }
            }

            for (int k = ((2 * size_edge_labels)); k < max_initial_label_size + (2 * size_edge_labels) ; ++k) {
                checkError(f, __LINE__, cudaMallocManaged((void **)&(thread_pool_list[f][j].labels[k].col_ring_size), sizeof(int) * common_queue_element_size));
                checkError(f, __LINE__, cudaMallocManaged((void **)&(thread_pool_list[f][j].labels[k].g), sizeof(int) * common_queue_element_size));
                checkError(f, __LINE__, cudaMallocManaged((void **)&(thread_pool_list[f][j].labels[k].h), sizeof(int) * common_queue_element_size));

                checkError(f, __LINE__, cudaMallocManaged((void **)&(thread_pool_list[f][j].labels[k].rings_g), sizeof(int *) * common_queue_element_size));
                for (int h = 0; h < 5; ++h) {
                    checkError(f, __LINE__, cudaMallocManaged((void **)&(thread_pool_list[f][j].labels[k].rings_g[h]), sizeof(int) * common_queue_element_size));
                }
            }
        }
    }
    clock_t end = clock();
     // Calculate elapsed time in seconds
    double elapsed_seconds = (double)(end - start) / CLOCKS_PER_SEC;
    malloc_elapsed_seconds = elapsed_seconds;


}

void kernel( vector<queue_elem> Q_filter  ) {

    //cout << "Kernel START" << endl;

    if(malloc_done == false){
        malloc();
        malloc_done = true;
    //    cout << "Kernel END OF MALLOC" << endl;
    }

    pointer_initialized();
  


    //initialize
    //init edge labels
    vectorToPointerEdge(gpu_edge_labels);
    //init adj matrix mol0
    vectorToPointerMatrix(g0, gpu_g0);
    size_gpu_g0_row = g0.size();
    size_gpu_g0_col = g0[0].size();
    //init adj matrix mol 1
    vectorToPointerMatrix(g1, gpu_g1);
    size_gpu_g1_row = g1.size();
    size_gpu_g1_col = g1[0].size();





    //copy the element of the Q_filter inside each QUEUE
    for(int i = 0; i < Q_filter.size(); i++) {
        auto_pool_size[i] = 1;
        auto_pool_len_m_best[i] = m_best.size();
        
        LabelFromCpuToGpu(thread_pool_list[i][0].labels, Q_filter[i].labels );
        thread_pool_list[i][0].labels_size = Q_filter[i].labels.size();
        thread_pool_list[i][0].m_size = Q_filter[i].m_local.size();
        for(int h = 0; h < Q_filter[i].m_local.size(); h++){
            thread_pool_list[i][0].m_local[h].first = Q_filter[i].m_local[h].first;
            thread_pool_list[i][0].m_local[h].second = Q_filter[i].m_local[h].second;
        }
    }

    
    autonomouslySolve<<<1,32>>>(thread_pool_list, auto_pool_size ,auto_pool_len_m_best, auto_pool_m_best, auto_pool_tmp ,Q_filter.size());
    // Synchronize threads
    cudaDeviceSynchronize();

    
    for(int i = 0; i < Q_filter.size(); i++) {
       if(auto_pool_len_m_best[i] > m_best.size()){
        m_best.clear();
        pair<int,int> tmp;
        for(int j = 0; j < auto_pool_len_m_best[i] ; j++){
            tmp.first = auto_pool_m_best[i][j].first;
            tmp.second = auto_pool_m_best[i][j].second;
            m_best.push_back(tmp);
        }
       }
    }

    return ;

}

