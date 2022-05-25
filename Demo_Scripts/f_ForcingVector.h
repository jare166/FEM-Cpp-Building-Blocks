#ifndef F_FORCINGVECTOR_H
#define F_FORCINGVECTOR_H

#include <math.h>
#include <iostream>
#include <tuple>

#include "C_FEM_BasisFunction_1D.h"
#include "C_FEM_GaussPoint_1D.h" 
#include "C_Mesh_1D.h"
#include "C_Matrix_Dense.h"

#include "f_MiscellaneousFunctions.h"

//! Construct/Assemble 1D Axial Stiffness Matrix
/*!
This function computes the element matrices and positions them in the global matrix, kGlobal, for
the linear elastic stiffness.

    \param flags object containing material flags
    \param material object containing material data
    \param kGlobal Pass-by-reference global sparse matrix. Initialized in COO form.
    \param fGlobal Pass-by-reference forcing vector.
    \param itEl Current element index which is being initialized.
    \param elNodes Vector containing start and end locations (in x) of element nodes.
    \param elDOF Degree-of-Freedom element matrix.
    \param GP_Data Object containing Gauss-Points.

Author: Dominic Jarecki
Date: 4-14-2022
*/
void stiffnessMatrix_AxialPlusBending_1D(int itEl, C_Mesh_1D mesh, C_Flags flags, C_Material material, 
    C_Matrix_Sparse& kGlobal, Eigen::VectorXd& fGlobal, C_LagrangeBasis* feL, C_GaussPoint_1D GP_Data){
    
    std::vector<int>      elDOF = mesh.elements[itEl];
    std::vector<double> elNodes = {mesh.nodes[elDOF[0]], mesh.nodes[elDOF[1]]};

    // Total Lengths:
    double le = elNodes[1] - elNodes[0];  // Axial Element

    // Recompute Hermite Shape Functions
    C_HermiteBasis feH(flags.interpOrder,  GP_Data, le);
    C_Matrix_Dense dsp_L, sp_H, ddsp_H;

    double GJ_d  = (2/le);
    double GJ_dd = pow((2/le),2);

    for (int itGp = 0; itGp < GP_Data.num_GP; itGp++) {
        // Number of Shape Functions associated w/ basis (for indexing)
        int ns_L = (*feL).num_sp;
        int ns_H = feH.num_sp;
        
        // Extract Rows at each Gauss Point
        dsp_L = (*feL).dsp(itGp, intspace(0,ns_L));

        sp_H   = feH.sp(itGp, intspace(0,ns_H));
        ddsp_H = feH.ddsp(itGp, intspace(0,ns_H));

        // Transformation to axial element
        dsp_L = GJ_d*dsp_L;

        double xL, xG_L, xG_H;
        double AX_L, FX;

        // Gaussian Quadrature
        double JxW = 0.5*GP_Data.wt[itGp]*le;
        
        AX_L = material.EA; // AX0 + AX1*xG_L;

        // i. Global Forcing Vector
        for (int i = 0; i < (*feL).num_sp; i++) {
            xL   = GP_Data.pt[i];                  // Local  Position
            xG_L = elNodes[0] + 0.5*(1.0 + xL)*le; // Global Position, Axial Element
            FX   = flags.FX0 + flags.FX1*xG_L;

            fGlobal(3*elDOF[i])   += FX*dsp_L(0, i)*JxW;    // Axial DOF
        }

        // ii. Global Stiffness Matrix
        // Products of Interpolation Functions
        C_Matrix_Dense S11 = dsp_L.T()*dsp_L*JxW;   // Axial DOF

        // Active Indices
        std::vector<int> ind_L = {3*elDOF[0],   3*elDOF[1]};   // Axial DOF

        // Add to Global Coefficient Matrix
        kGlobal.add_matr(AX_L*S11, ind_L, ind_L); // Axial DOF
        
    }
}


#endif