//
// Created by davide on 4/19/24.
//

#include <fstream>
#include <iostream>
#include <vector>
#include <string>
#include "test.hpp"
#include <ctime>
#include <chrono>  // Include la libreria chrono

using namespace std;
using namespace RDKit;

int max_l0_size = 0;
int max_l1_size = 0; 
int max_first_len_initialized = 0;

int edge_label_size = 4;
int max_initial_label_size = 0;
bool state_initialized = true;




int main()
{

    int *pointer_tmp;
    cudaMallocManaged(&pointer_tmp, sizeof(int) * 2);

    
  string s0 = "CC(C)[C@@H]1CC[C@@H](C)C[C@H]1OC(=O)c2ccccc2c3c(C)cccc3CN(C)C4CCCCC4";
    string s1 = "O=C(OC1CCCCC1)c2ccccc2c3ccccc3CNC4CCCCC4";
 ROMol result = smiles_mcs(s0, s1 );
 state_initialized = false;

    clock_t start = clock();
     result = smiles_mcs(s0, s1 );
    clock_t end = clock();

    
    // Calculate elapsed time in seconds
    double elapsed_seconds = (double)(end - start) / CLOCKS_PER_SEC;

    std::vector<std::string> result_string;
    for (const auto &atom : result.atoms()) {
        result_string.push_back(atom->getSymbol());
    }


    cout << "[";
    for ( int idx = 0; idx < result_string.size(); idx++ ){
        if(idx == result_string.size()-1 ){
            cout <<"'"<<result_string.at(idx)<<"']"<<endl;
        }
        else cout <<"'"<<result_string.at(idx)<<"', ";
    }
    cout<<"done\n\n";

     // Print the elapsed time in seconds
    std::cout << "\nElapsed time: " << elapsed_seconds << " seconds" << std::endl;
    std::cout << "\nElapsed time [WITHOUT MALLOC]: " << elapsed_seconds - malloc_elapsed_seconds << " seconds" << std::endl;

    return 0;
}
