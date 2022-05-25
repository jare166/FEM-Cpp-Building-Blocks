#include <eigen-3.4.0/Eigen/Dense>  // Linear Solver Libraries
#include <eigen-3.4.0/Eigen/Sparse>
#include <iostream> // Read/Write to terminal
#include <iomanip>  // Set width of output to terminal
#include <fstream>  // Read/Write to file
#include <sstream>  // Read/Write from string
#include <vector> 
#include <ctime>    // Auto filenames

#include "C_FEM_BasisFunction_1D.h"
#include "C_FEM_GaussPoint_1D.h" 
#include "C_Mesh_1D.h"
#include "C_Matrix_Sparse.h"
#include "C_Matrix_Dense.h"

#include "f_ForcingVector.h"
#include "f_SolverInterfaces.h"

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

class C_Flags{
    public:
        int num_Elements = 1;    //!< Number of Elements in Mesh
        int num_GaussPoints = 2; //!< Number of Gauss Points per element
        int interpOrder = 1;     //!< Interpolation Order used for approximation functions
        int totDof = 1;          //!< Total number of degrees of freedom associated w/ FEM problem

        double length_Mesh = 1;  //!< Total Length of Mesh
        //@{
        //! Finite-Element Constants (1D Problem)
        /*! The following constants are associated w/ the forcing vector to a 2nd-order PDE:
            \f$ f = (f_{x0} + f_{x1}*x) \f$
        */
        double FX0 = 0;
        double FX1 = 0;
        //@}

        /*!
        Boundary Conditions specified as:
        DBC,     node, DOF, value
             or
        NBC,     node, DOF, value
        */
        std::vector<int> DBC_nods; //!< Vector containing nodal values at which Dirichlet BCs are specified
        std::vector<int> NBC_nods; //!< Vector containing nodal values at which Neumann BCs are specified

        std::vector<int> DBC_DOFs; //!< Vector containing corresponding DOF 
        std::vector<int> NBC_DOFs; //!< Vector containing corresponding DOF 

        std::vector<double> DBC_vals; //!< Vector containing corresponding values
        std::vector<double> NBC_vals; //!< Vector containing corresponding values

        //! Miscellaneous Preferences
        int show_matrix_preBC  = 0;
        int show_matrix_postBC = 0;

    C_Flags(){}
};


int main()
{
    // 1. User Selections
    C_Flags flags;
    C_Material material;

    // 2. Initialize Governing Objects
    C_GaussPoint_1D GP_Data(flags.num_GaussPoints);
    C_Mesh_1D       mesh(flags.length_Mesh, flags.num_Elements+1);
    C_LagrangeBasis feL(flags.interpOrder, GP_Data);

    
    // 3. Populate Global Sparse Matrix
    C_Matrix_Sparse kGlobal;
    Eigen::VectorXd fGlobal(flags.totDof);
    for (int ii = 0; ii < flags.totDof; ii++) { fGlobal(ii) = 0.0; }

    for (int itEl = 0; itEl < mesh.num_El; itEl++) {
        stiffnessMatrix_AxialPlusBending_1D(itEl, mesh, flags, material, kGlobal, fGlobal, &feL, GP_Data);
    }

    // [OPTIONAL] Display Sparse Matrix (pre-BC)
    if (flags.show_matrix_preBC == 1) {
        std::cout << "\n Matrix Values (pre-BC): \n";
        std::cout << kGlobal;
    }


    // 4. Apply Boundary Conditions
    // i. Neumann
    for (int ii = 0; ii < flags.NBC_vals.size(); ii++) { fGlobal(3*flags.NBC_nods[ii] + (flags.NBC_DOFs[ii]-1)) = flags.NBC_vals[ii]; }

    // ii. Dirichlet
    int it;
    for (int ii = 0; ii < flags.DBC_vals.size(); ii++) { 
        it = 3*flags.DBC_nods[ii] + (flags.DBC_DOFs[ii]-1);

        kGlobal.col_NonSparseAssign(0.0, it);
        kGlobal.row_NonSparseAssign(0.0, it);
        kGlobal(it,it) = 1;

        fGlobal(it) = flags.DBC_vals[ii]; 
    }

    
    // 5. Solve, using Cholesky Factorization of kGlob
    Eigen::SparseMatrix<double> kG_eigen(kGlobal.row_size+1, kGlobal.col_size+1);
    convert_to_Eigen(kGlobal, kG_eigen);

    // [OPTIONAL] Display Eigen Sparse Matrix (post-BC)
    if (flags.show_matrix_postBC == 1) {
        std::cout << "\n Matrix Values (post-BC): \n";
        std::cout << "\n" << Eigen::MatrixXd(kG_eigen) << "\n";
    }

    Eigen::SimplicialLDLT< Eigen::SparseMatrix<double> > chol;
    chol.compute(kG_eigen);  
    Eigen::VectorXd sol = chol.solve(fGlobal);

    return 0;
}