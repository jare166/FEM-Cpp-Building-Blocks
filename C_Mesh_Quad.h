#ifndef C_MESH_QUAD_H
#define C_MESH_QUAD_H

#include <math.h>
#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>

#include "f_MiscellaneousFunctions.h"
#include "C_Mesh.h"

// DECLARATIONS
template<typename T> std::vector<T> linspace(int num_in, T start_in, T end_in);

//! This class is used to store data describing a mesh of 1, 2, or 3 dimensions.
class C_Mesh_Quad : public C_Mesh{
    /* 
    INHERITED PROPERTIES:

    public:
        std::vector< std::vector<double> > nodes;
        std::vector< std::vector<int> >    elements;
        int num_NPE = 2; //! Nodes per element
        int num_Nd  = 0; //! Total number of nodes
        int num_El  = 0; //! Total number of elements
        int dim = 2;     //! Dimension of mesh
    */

    public:
    
    C_Mesh_Quad() {  }
};

#endif