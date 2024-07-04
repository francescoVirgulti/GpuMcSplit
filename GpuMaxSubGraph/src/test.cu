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



int main()
{

    int *pointer_tmp;
    cudaMallocManaged(&pointer_tmp, sizeof(int) * 2);

    
    string s0 = "Fc1cc(NCCN2C=Cn3nc(cc3C2=O)c4occc4)ccc1C(=O)CC5CCC5";
    string s1 = "O=C(CC1CCC1)c2ccc(NCCN3C=Cn4nc(cc4C3=O)c5occc5)cc2";


    clock_t start = clock();
    ROMol result = smiles_mcs(s0, s1 );
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
