#ifndef C_FEM_BASISFUNCTION_1D_H
#define C_FEM_BASISFUNCTION_1D_H

#include <vector>
#include <math.h>
#include <iostream>
#include <tuple>

#include "C_FEM_GaussPoint_1D.h" 
#include "C_Matrix_Dense.h"


//! This class produces Hermite basis functions at Gauss Points for an element.
class C_HermiteBasis{
    public:    
        C_Matrix_Dense sp;   //!< Interpolation functions
        C_Matrix_Dense dsp;  //!< 1st Derivative of Interpolation functions
        C_Matrix_Dense ddsp; //!< 2nd Derivative of Interpolations Functions
        int num_sp;

    //! Class Constructor
    /*!
        \param order the order of the requested interpolation functions. Not used here, included for parallelism.
        \param GP_Data object containing Gauss Point information, such as coordinate data.
        \param le length of element.
    */
    C_HermiteBasis (int order, C_GaussPoint_1D GP_Data, double le){
        num_sp = 4;
        
        // Assign new, correctly-sized dense matrices to member variables
        this -> sp = C_Matrix_Dense(GP_Data.num_GP, num_sp);
        this -> dsp = C_Matrix_Dense(GP_Data.num_GP, num_sp);
        this -> ddsp = C_Matrix_Dense(GP_Data.num_GP, num_sp);

        for (int ii = 0; ii < GP_Data.num_GP; ii++){
            double x = GP_Data.pt[ii]; // Evaluate x at Gauss Points

            // (v1, th1, v2, th2)
            sp(ii,0) = 0.25*pow((1-x),2)*(2+x);
            sp(ii,1) = le/8*pow((1-x),2)*(x-1);
            sp(ii,2) = 0.25*pow((1+x),2)*(2+x);
            sp(ii,3) = le/8*pow((1+x),2)*(x+1);

            ddsp(ii,0) = (1/le)*(-6*x/le); 
            ddsp(ii,1) = (1/le)*(3*x-1);
            ddsp(ii,2) = (1/le)*(6*x/le);
            ddsp(ii,3) = (1/le)*(3*x+1);
        }
    }
};


//! This class produces 1st- and 2nd-order Lagrange basis functions at Gauss Points for an element.
class C_LagrangeBasis{
    public:
        C_Matrix_Dense sp;  //!< Interpolation functions
        C_Matrix_Dense dsp; //!< 1st Derivative of Interpolation functions
        int num_sp;

    //! Class Constructor
    /*!
        \param order the order of the requested interpolation functions. 
        \param GP_Data object containing Gauss Point information, such as coordinate data.
    */
    C_LagrangeBasis (int order, C_GaussPoint_1D GP_Data){
        if (order == 1) {
            // Linear
            num_sp = 2;

            // Assign new, correctly-sized dense matrices to member variables
            this -> sp = C_Matrix_Dense(GP_Data.num_GP, num_sp);
            this -> dsp = C_Matrix_Dense(GP_Data.num_GP, num_sp);

            for (int ii = 0; ii < GP_Data.num_GP; ii++){
                double x = GP_Data.pt[ii]; // Evaluate x at Gauss Points

                sp(ii,0) = 0.5*(1-x);
                sp(ii,1) = 0.5*(1+x);

                dsp(ii,0) =-0.5;
                dsp(ii,0) = 0.5;
            }
        }
        else {
            // Quadratic
            num_sp = 3;
            
            // Assign new, correctly-sized dense matrices to member variables
            this -> sp = C_Matrix_Dense(GP_Data.num_GP, num_sp);
            this -> dsp = C_Matrix_Dense(GP_Data.num_GP, num_sp);

            for (int ii = 0; ii < GP_Data.num_GP; ii++){
                double x = GP_Data.pt[ii]; // Evaluate x at Gauss Points

                sp(ii,0) = 0.5*(x - 1)*x;
                sp(ii,1) = 1 - x*x;
                sp(ii,2) = 0.5*(x + 1)*x;

                dsp(ii,0) = x - 0.5;
                dsp(ii,1) = -2*x;
                dsp(ii,2) = x + 0.5;  
            }
        }
    }
};

#endif