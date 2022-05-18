#ifndef C_FEM_GAUSSPOINT_1D_H
#define C_FEM_GAUSSPOINT_1D_H

#include <vector>
#include <math.h>
#include <iostream>

class C_GaussPoint_1D{
    public:
        int num_GP; //!< Number of Gauss Points
        std::vector<double> pt; //!< Vector of Gauss Points
        std::vector<double> wt; //!< Corresponding Gauss Weights
    
    C_GaussPoint_1D (int n_in){
        num_GP = n_in;

        if (num_GP == 1){
            pt = {0.0};
            wt = {2.0};
        }
        else if (num_GP == 2){
            pt = {0.577350269189626, -0.577350269189626};
            wt = {1.0,                1.0};
        }
        else if (num_GP == 3){
            pt = {-0.774596669241483, 0,   0.774596669241483};
            wt = { 5/9,               8/9, 5/9};
        }
        else if (num_GP == 4){
            pt = {-0.86113631, -0.33998104, 0.33998104, 0.86113631};
            wt = { 0.34785485,  0.65214515, 0.65214515, 0.34785485};
        }
        else {
            pt = {-0.906180, -0.538469, 0.0,      0.538469, 0.906180};
            wt = { 0.236927,  0.478629, 0.568889, 0.478629, 0.236927};
            num_GP = 5;
        }
    }
};

#endif