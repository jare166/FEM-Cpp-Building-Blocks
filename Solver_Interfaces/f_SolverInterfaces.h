#ifndef f_SOLVERINTERFACES_H
#define f_SOLVERINTERFACES_H

#include <eigen-3.4.0/Eigen/Dense>
#include <eigen-3.4.0/Eigen/Sparse>

#include <math.h>
#include <vector>
#include <iostream>
#include <iterator>
#include <algorithm>
#include <forward_list>

#include "../C_Matrix_Sparse.h"
#include "../C_Matrix_Dense.h"

void convert_to_Eigen(C_Matrix_Sparse& kG, Eigen::SparseMatrix<double>& kG_eigen) {
    /*!
    Initialize an array suitable for use with Eigen linear solver.

    DIJ (4-19-22)
    */

    std::vector< Eigen::Triplet<double> > temp_COO(kG.NNZ);

    // Iterate over rows
    int kk = 0;
    for (int jj = 0; jj <= kG.row_size; jj++) {
        std::forward_list<double>::iterator it_v = kG.value_list[jj].begin();
        std::forward_list<int>::iterator    it_c = kG.col_ind_list[jj].begin();

        int init_length = std::distance(kG.value_list[jj].begin(), kG.value_list[jj].end());

        // Iterate over columns
        for (int ii = 0; ii < init_length; ii++) {
            // Store as (row, column, value) for each element of Triplet Vector
            temp_COO[kk++] = Eigen::Triplet<double> (jj, *it_c, *it_v);
            ++it_v; ++it_c;
        }
    }

    // Store as Eigen SparseMatrix
    kG_eigen.setFromTriplets(temp_COO.begin(), temp_COO.end());
    // kG_eigen.makeCompressed();
}

#endif