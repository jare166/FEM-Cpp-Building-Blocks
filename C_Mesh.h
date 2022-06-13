#ifndef C_MESH_H
#define C_MESH_H

#include <math.h>
#include <vector>
#include <iostream> // Read/Write to terminal
#include <fstream>  // Read/Write to file
#include <sstream>  // Read/Write from string

#include "f_MiscellaneousFunctions.h"

// DECLARATIONS
template<typename T> std::vector<T> linspace(int num_in, T start_in, T end_in);

//! This class is used to store data describing a mesh of 1, 2, or 3 dimensions.
class C_Mesh{
    public:
        std::vector< std::vector<double> > nodes;
        std::vector< std::vector<int> >    elements;
        int num_NPE = 0; //! Nodes per element
        int num_Nd  = 0; //! Total number of nodes
        int num_El  = 0; //! Total number of elements
        int dim = 2;     //! Dimension of mesh

    C_Mesh() { }

    void write_connectivity(std::string filePath) {
        std::fstream fID;
        fID.open(filePath);

        // i. List Number of Lines
        fID << "num_lines," << num_El << "\n";

        // ii. List Connectivity
        for (int ii = 0; ii < num_El; ii++) {  fID << elements[ii][0] << "," << elements[ii][1] << "\n"; }

        fID.close();
    }

    void write_nodes(std::string filePath) {
        std::fstream fID;
        fID.open(filePath);

        // i. List Number of Lines
        fID << "num_lines," << num_Nd << "\n";

        for (int ii = 0; ii < num_Nd; ii++) { 
            // ii. List Node Number
            fID << ii << ", ";

            // iii. List Nodes
            for (int jj = 0; jj < dim; jj++) { fID << nodes[ii][jj] << ","; } 
            fID << "\n";
        }

        fID.close();
    }
};

#endif