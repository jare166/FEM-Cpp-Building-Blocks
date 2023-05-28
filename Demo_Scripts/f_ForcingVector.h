#ifndef F_FORCINGVECTOR_H
#define F_FORCINGVECTOR_H

#include <math.h>
#include <iostream>
#include <tuple>

#include "../C_FEM_BasisFunction_1D.h"
#include "../C_FEM_GaussPoint_1D.h" 
#include "../C_Mesh_Frame.h"
#include "../C_Matrix_Dense.h"

#include "../f_MiscellaneousFunctions.h"

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
class C_Material{
    public:
        // Elastic Constants
        double EI = 1;
        double EA = 1;
        
        // Elastic-Plastic Constants
        double my = 15.66; // Plastic Yield Moment
        double a  = 10.10; // Sharpness Parameter

    C_Material(){}
};

void stiffnessMatrix_AxialBar(int itEl, C_Mesh_Frame& mesh, C_Material& mat, 
    C_Matrix_Sparse& kGlobal, C_LagrangeBasis_1D& feL, C_GaussPoint_1D& GP_Data){
    
    int elDOF_0 = mesh.elements(itEl,0);
    int elDOF_1 = mesh.elements(itEl,1);
    std::vector<double> elNodes = {mesh.nodes(elDOF_0, 0), mesh.nodes(elDOF_1, 0)};

    double le = elNodes[1] - elNodes[0];

    C_Matrix_Dense dsp_L;
    for (int itGp = 0; itGp < GP_Data.num_GP; itGp++) {

        int ns_L = feL.num_sp;
        
        // Extract Rows at each Gauss Point
        dsp_L = (2/le)*feL.dsp(itGp, intspace(0,ns_L));
        double JxW = 0.5*GP_Data.wt[itGp]*le;
        C_Matrix_Dense S11 = JxW*mat.EA*(dsp_L.T()*dsp_L);

        // Add to Global Coefficient Matrix
        kGlobal.add_matr(S11, {elDOF_0, elDOF_1}, {elDOF_0, elDOF_1});
    }
}


#endif