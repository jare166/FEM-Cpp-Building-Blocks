#ifndef C_MESH_1D_H
#define C_MESH_1D_H

#include <math.h>
#include <vector>

#include "f_MiscellaneousFunctions.h"

// DECLARATIONS
template<typename T> std::vector<T> linspace(int num_in, T start_in, T end_in);

//! This class is used to store data describing a 1D mesh.
class C_Mesh_1D{
    public:
        std::vector<double> nodes;
        std::vector< std::vector<int> > elements;
        int num_Nd = 0;
        int num_El = 0;

    C_Mesh_1D(double l, int num_Nd_in){
        num_Nd = num_Nd_in;
        num_El = num_Nd-1;

        nodes = linspace(num_Nd, 0.0, l);
        elements.resize(num_Nd, std::vector<int>(2));
        
        for (int ii = 0; ii < num_El; ii++) { elements[ii] = {ii, ii+1}; }
    }
};

#endif