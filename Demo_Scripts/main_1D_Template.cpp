#include <eigen-3.4.0/Eigen/Dense>  // Linear Solver Libraries
#include <eigen-3.4.0/Eigen/Sparse>

#include "C_FEM_BasisFunction_1D.h"
#include "C_FEM_GaussPoint_1D.h" 
#include "C_Mesh_1D.h"
#include "C_Matrix_Sparse.h"
#include "C_Matrix_Dense.h"

#include "f_ForcingVector.h"
#include "f_SolverInterfaces.h"


int main()
{
    C_Material mat;
    C_GaussPoint_1D GP_Data(2);
    C_Mesh_1D       mesh(10.0, 11);
    C_LagrangeBasis feL(1, GP_Data);
    int totDof=mesh.num_Nd;
    C_Matrix_Sparse kGlobal;
    Eigen::VectorXd fGlobal(totDof);
    fGlobal.setZero();

    for (int itEl = 0; itEl < mesh.num_El; itEl++) {
        stiffnessMatrix_AxialBar(itEl, mesh, mat, kGlobal, feL, GP_Data);
    }

    // i. Neumann Boundary Conditions
    int fDof=totDof-1;
    fGlobal(fDof)=10.0;

    // ii. Dirichlet Boundary Conditions
    int bcDof=0;
    kGlobal.col_NonSparseAssign(0.0, bcDof);
    kGlobal.row_NonSparseAssign(0.0, bcDof);
    kGlobal(bcDof,bcDof) = 1;
    fGlobal(bcDof) = 0.0; 

    
    // 5. Solve, using Cholesky Factorization of kGlob
    Eigen::SparseMatrix<double> kG_eigen(totDof, totDof);
    convert_to_Eigen(kGlobal, kG_eigen);

    Eigen::SimplicialLDLT< Eigen::SparseMatrix<double> > chol;
    chol.compute(kG_eigen);  
    Eigen::VectorXd sol = chol.solve(fGlobal);
    std::cout << sol;
    std::cout << "\n";

    return 0;
}