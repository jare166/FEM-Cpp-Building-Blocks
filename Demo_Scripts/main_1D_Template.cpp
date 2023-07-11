#include <eigen-3.4.0/Eigen/Dense>  // Linear Solver Libraries
#include <eigen-3.4.0/Eigen/Sparse>

#include "../C_FEM_BasisFunction_1D.h"
#include "../C_FEM_GaussPoint_1D.h" 
#include "../C_Mesh_Frame.h"
#include "../C_Matrix_Sparse.h"
#include "../C_Matrix_Dense.h"

#include "f_ForcingVector.h"
#include "../Solver_Interfaces/f_SolverInterfaces.h"


int main()
{
    // 0. Initialize Governing Objects
    C_Material mat;
    C_Mesh_Frame mesh(3,2,5);

    int kk_n = 0, kk_l = 0;
    mesh.construct_elems( 0.0, 0.0, 0.0,  1.0, 0.0, 0.0,  1, 2,  5, kk_n, kk_l);

    C_GaussPoint_1D GP_Data(2);
    C_LagrangeBasis_1D feL(1, GP_Data);

    int totDof = mesh.num_Nd;
    C_Matrix_Sparse kGlobal(totDof);
    Eigen::VectorXd fGlobal(totDof);
    fGlobal.setZero();

    // 1. Assemble sub-matrices
    for (int itEl = 0; itEl < mesh.num_El; itEl++) { stiffnessMatrix_AxialBar(itEl, mesh, mat, kGlobal, feL, GP_Data); }

    // Assign Boundary Conditions
    //  i. Neumann BCs
    int fDof=totDof-1;
    fGlobal(fDof)=10.0;

    //  ii. Dirichlet BCs
    int bcDof=0;
    kGlobal.col_NonSparseAssign(0.0, bcDof);
    kGlobal.row_NonSparseAssign(0.0, bcDof);
    kGlobal(bcDof,bcDof) = 1;
    fGlobal(bcDof) = 0.0; 

    // 2. Solve, using Cholesky Factorization of kGlob
    Eigen::SparseMatrix<double> kG_eigen(totDof, totDof);
    convert_to_Eigen(kGlobal, kG_eigen);

    Eigen::SimplicialLDLT< Eigen::SparseMatrix<double> > chol;
    chol.compute(kG_eigen);  
    Eigen::VectorXd sol = chol.solve(fGlobal);
    std::cout << sol;
    std::cout << "\n";

    return 0;
}