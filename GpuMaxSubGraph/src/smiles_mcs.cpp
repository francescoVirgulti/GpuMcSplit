#include <rdkit/GraphMol/GraphMol.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <RDGeneral/test.h>
#include "main.hpp"
using namespace RDKit;

//translates the SMILES into Mol and call the main function to compute MCS
ROMol smiles_mcs( std::string& smile0,  std::string& smile1, int bond_match , int ring_match ) {

    RWMol mol0 = *SmilesToMol(smile0);
    RWMol mol1 = *SmilesToMol(smile1);

    return mol_mcs(mol0, mol1, 1,1,0);
}