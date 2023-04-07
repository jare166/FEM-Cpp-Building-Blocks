#ifndef C_FEM_BASISFUNCTION_2D_H
#define C_FEM_BASISFUNCTION_2D_H

#include <math.h>
#include <iostream>

#include "C_FEM_GaussPoint_1D.h" 
#include "C_Matrix_Dense.h"


//! This class produces 1st- and 2nd-order Lagrange basis functions at Gauss Points for a 1-dimensional element.
class C_LagrangeBasis_2D{
    public:
        C_Matrix_Dense sp;  //!< Interpolation functions
        C_Matrix_Dense dsp; //!< 1st Derivative of Interpolation functions
        int num_sp;

    //! Class Constructor
    /*!
        \param nodes_per_elem indirectly specifies the order of the requested interpolation functions. 
        \param GP_Data_x object containing Gauss Point x-coordinate data.
        \param GP_Data_y object containing Gauss Point y-coordinate data.
    */
    C_LagrangeBasis_2D (int nodes_per_elem, C_GaussPoint_1D GP_Data_x, C_GaussPoint_1D GP_Data_y){

        // Assign new, correctly-sized dense matrices to member variables
        this -> sp  = C_Matrix_Dense(GP_Data_x.num_GP, 1);
        this -> dsp = C_Matrix_Dense(GP_Data_x.num_GP, 2);

        std::vector<int> sgn_x = {-1, 1, 1,-1, 0, 1, 0,-1, 0};
        std::vector<int> sgn_y = {-1,-1, 1, 1,-1, 0, 1, 0, 0};

        if (nodes_per_elem == 4) {
            // 1. Linear, 4 nodes per element
            for (int ii = 0; ii < nodes_per_elem; ii++){
                double x = GP_Data_x.pt[ii]; // Evaluate x at Gauss Points
                double y = GP_Data_y.pt[ii]; // Evaluate y at Gauss Points

                sp(ii,0)  = 0.25*(1 + x*sgn_x[ii])*(1 + y*sgn_y[ii]);

                dsp(ii,0) = 0.25*(1 + y*sgn_y[ii])*x;
                dsp(ii,1) = 0.25*(1 + x*sgn_x[ii])*y;
            }
        }
        else if (nodes_per_elem == 8) {
            // 2. Quadratic, 8 nodes per element
            for (int ii = 0; ii < nodes_per_elem; ii++){
                double x = GP_Data_x.pt[ii]; // Evaluate x at Gauss Points
                double y = GP_Data_y.pt[ii]; // Evaluate y at Gauss Points

                double xp_1, yp_1, xp_2, yp_2, XI2, ETA2;
                xp_1 = 1.0 + x*sgn_x[ii];
                yp_1 = 1.0 + y*sgn_y[ii];
                xp_2 = 1.0 - pow(x,2);
                yp_2 = 1.0 - pow(y,2);

                if(ii < 4) {
                    sp(ii,0)  = 0.25*xp_1*yp_1*(x*sgn_x[ii] + y*sgn_y[ii] - 1.0);
                    dsp(ii,0) = 0.25*yp_1*sgn_x[ii]*(2.0*x*sgn_x[ii] + y*sgn_y[ii]);
                    dsp(ii,1) = 0.25*xp_1*sgn_y[ii]*(2.0*y*sgn_y[ii] + x*sgn_x[ii]);
                } 
                else if (ii < 6) {
                    sp(ii,0)  = 0.5*xp_2*yp_1;
                    dsp(ii,0) = -x*yp_1;
                    dsp(ii,1) = 0.5*sgn_y[ii]*xp_2;
                } 
                else {
                    sp(ii,0)  = 0.5*yp_2*xp_1;
                    dsp(ii,0) = 0.5*sgn_x[ii]*yp_2;
                    dsp(ii,1) = -y*xp_1;
                }

            }
        }
        else {
            // 3. Quadratic, 9 nodes per element (serendipity element)
            for (int ii = 0; ii < nodes_per_elem; ii++){
                double x = GP_Data_x.pt[ii]; // Evaluate x at Gauss Points
                double y = GP_Data_y.pt[ii]; // Evaluate y at Gauss Points

                double xp_1, yp_1, xp_2, yp_2, XI2;
                xp_1 = 1.0 + x*sgn_x[ii];
                yp_1 = 1.0 + y*sgn_y[ii];
                xp_2 = 1.0 - pow(x,2);
                yp_2 = 1.0 - pow(y,2);

                if (ii < 4) {
                    sp(ii,0)  = 0.25*sgn_x[ii]*sgn_y[ii]*x*y*xp_1*yp_1;
                    dsp(ii,0) = 0.25*sgn_x[ii]*sgn_y[ii]*y*yp_1*(1.0+2.0*sgn_x[ii]*x);
                    dsp(ii,1) = 0.25*sgn_y[ii]*sgn_x[ii]*x*xp_1*(1.0+2.0*sgn_y[ii]*y);
                } 
                else if(ii < 6) {
                    sp(ii,0)  = 0.5*xp_2*yp_1*sgn_y[ii]*y;
                    dsp(ii,0) = -x*sgn_y[ii]*y*yp_1;
                    dsp(ii,1) = 0.5*xp_2*sgn_y[ii]*(1.0+2.0*sgn_y[ii]*y);
                } 
                else if(ii < 8) {
                    sp(ii,0)  = 0.5*yp_2*xp_1*sgn_x[ii]*x;
                    dsp(ii,0) = -y*sgn_x[ii]*x*xp_1;
                    dsp(ii,1) = 0.5*yp_2*sgn_x[ii]*(1.0+2.0*sgn_x[ii]*x);
                } 
                else {
                    sp(ii,0)  = (xp_2*yp_2);
                    dsp(ii,0) = -2.0*x*yp_2;
                    dsp(ii,1) = -2.0*y*xp_2;
                }  
            }
        }
    }
};

#endif