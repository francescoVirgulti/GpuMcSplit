//
// Created by davide on 4/19/24.
//
#include <rdkit/GraphMol/GraphMol.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <fstream>
#include <iostream>
#include <vector>
#include <string>
#include "test.hpp"
#include "cuda_header.h"
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


//auto is for autonomous
ThreadVar **thread_pool_list; 
int *auto_pool_size;
Pair **auto_pool_m_best;
int *auto_pool_len_m_best;
ThreadVar *auto_pool_tmp;

vector<int> length_list;

Pair *m_best_solution;

float *main_gpu_edge_labels;
float **main_gpu_g0;
float **main_gpu_g1;

bool malloc_done = false;


int main()
{

    int *pointer_tmp;
    cudaMallocManaged(&pointer_tmp, sizeof(int) * 2);

      std::vector<std::pair<std::string, std::string>> molecules;

        std::string filename = "input.txt";
        std::string outputFilename = "output.txt";
        // Open the output file for writing
        std::ofstream outputFile(outputFilename);
        if (!outputFile) {
            std::cerr << "Error opening output file: " << outputFilename << std::endl;
            return 1;
        }
        //std::vector<std::pair<std::string, std::string>> smiles;
        std::ifstream file(filename);
        int skip;
        std::string first, second;
        ROMol result;
        int i = 0;


        if (file.is_open()) {
            std::string line;
            while (std::getline(file, line) ) {

                pair<string,string> tmp;
                std::istringstream iss(line);
                
                // Skip the first int
                iss >> skip;
                // Read the first string
                iss >> first;
                // Skip the next four ints
                for (int i = 0; i < 4; ++i) {
                    iss >> skip;
                }
                // Read the second string
                iss >> second; 
                clock_t start = clock();
                // Store the pair of strings
                tmp.first = first;
                tmp.second = second;
                molecules.push_back(tmp);
                
                i++;
            }
        } else {
            std::cerr << "Error opening file: " << filename << std::endl;
        }

        cout << "molecules size : " << molecules.size() <<endl;


        clock_t start = clock();
        for(pair mol_pair : molecules ){
            smiles_mcs(mol_pair.first, mol_pair.second);
        }
        clock_t end = clock();
        state_initialized = false;
         double elapsed_seconds = (double)(end - start) / CLOCKS_PER_SEC;


        std::vector<std::string> result_string;
        int index = 0;
        for(pair mol_pair : molecules ){
            cout << "__________________________MOLECULES PAIR NUM : " << index << endl;
            index++;
            result.clear();
            RWMol mol0 = *SmilesToMol(mol_pair.first);
            RWMol mol1 = *SmilesToMol(mol_pair.second); 
            clock_t start = clock();
            result = mol_mcs(mol0, mol1, 1,1,0);
            clock_t end = clock();
            elapsed_seconds = elapsed_seconds + (double)(end - start) / CLOCKS_PER_SEC;
            result_string.clear();
            for (const auto &atom : result.atoms()) {
                result_string.push_back(atom->getSymbol());
            }

                outputFile << "[";
                for ( int idx = 0; idx < result_string.size(); idx++ ){
                    if(idx == result_string.size()-1 ){
                        outputFile <<"'"<<result_string.at(idx)<<"']" << endl;;
                    }
                    else outputFile <<"'"<<result_string.at(idx)<<"', ";
                }
        }

    




    
    
    // Calculate elapsed time in seconds
    //double elapsed_seconds = (double)(end - start) / CLOCKS_PER_SEC;

    



     // Print the elapsed time in seconds
    std::cout << "\nElapsed time: " << elapsed_seconds << " seconds" << std::endl;
    std::cout << "\n MALLOC Elapsed time: " << malloc_elapsed_seconds << " seconds" << std::endl;
    std::cout << "\nElapsed time [WITHOUT MALLOC]: " << elapsed_seconds - malloc_elapsed_seconds << " seconds" << std::endl;

    return 0;
}
