#ifndef C_MATRIX_SPARSE_H
#define C_MATRIX_SPARSE_H

#include <eigen-3.4.0/Eigen/Dense>
#include <eigen-3.4.0/Eigen/Sparse>
#include <math.h>
#include <vector>
#include <iostream>
#include <iterator>
#include <algorithm>
#include <forward_list>

#include "C_Matrix_Dense.h"

//! This code initializes sparse matrix objects for use w/ various linear solvers.
/*!
When initializing, the data is stored as Coordinate Lists (COO).

Once initialized, convert to CSR for use w/ needed algorithms:
    C_Matrix_Sparse::convert_to_CSR()

Author: Dominic Jarecki
Date: 3-24-2022
*/
class C_Matrix_Sparse{   
    public:
        //! CONSTRUCTION:
        //@{
        /*! 
        COO Coordinate List
        Vector of linked lists; using this format allows fast access of elements by row. 
        */
        std::vector<std::forward_list<double>> value_list;
        std::vector<std::forward_list<int>>    col_ind_list;
        //@}

        //! MATRIX CALCULATIONS:
        //@{
        /*! 
        CSR (Compressed Sparse Row)
        */
        std::vector<double> value;
        std::vector<int>    col_ind, row_ind;
        //@}

        // Miscellaneous
        int row_size = -1, col_size = -1; //! Track size as object grows
        int NNZ = 0;                      //! Number of non-zero elements

        int active_flag = 1; //! Active data format; 1: COO, 2: CSR, 3: CSC

    // I. Constructor
    C_Matrix_Sparse(){}

    C_Matrix_Sparse(int row_size_in, int col_size_in) {
        row_size = row_size_in;
        col_size = col_size_in;
        NNZ = 0;

        // Add elements to vectors of lists
        for (int ii = 0; ii < row_size; ii++) {
            col_ind_list.push_back(std::forward_list<int>()); 
            value_list.push_back(std::forward_list<double>()); 
        }
    }

    // II. Access Contents
    double& operator() (int r_i, int c_i) { 
        //! Access Contents of List
        /*! 
        This function returns elements contained in C_Matrix_Sparse. Currently, they
        are accessed through COO only, not through CSR.

        DIJ (4-07-22)
        */
        // COO       
        // Ensure an element is present at queried location
        add_elem(r_i, c_i, 0.0); 

        // Iterate through list
        std::forward_list<double>::iterator it_v = value_list[r_i].begin();
        std::forward_list<int>::iterator    it_c = col_ind_list[r_i].begin();

        int init_length = std::distance(value_list[r_i].begin(), value_list[r_i].end());

        for (int ii = 0; ii < init_length; ii++){
            // Check row
            if (*it_c == c_i) { break; }

            ++it_v;
            ++it_c;
        }

        // Access value assigned
        double& access_val = *it_v;
        return access_val;
        }

    //  i. row and column query and/or vector assignment
    // std::vector<double>& row(int r_i) {}
    // std::vector<double>& col(int c_i) {}

    //  ii. row and column scalar assignment
    void row_NonSparseAssign(int r_i, double val) {
        //! This function assigns all NON-SPARSE vals in row r_i to be val.
        /*!
        SPARSE VALUES ARE NOT MODIFIED.
        NOTE: This is a quick operation.
        DIJ (4-13-22)
        */
        // Iterate through values in row r_i
        std::forward_list<double>::iterator it_v = value_list[r_i].begin();

        int init_length = std::distance(value_list[r_i].begin(), value_list[r_i].end());

        for (int ii = 0; ii < init_length; ii++){
            *it_v = val;
            ++it_v;
        }
    }
    
    void col_NonSparseAssign(int c_i, double val) {
        //! This function assigns all NON-SPARSE vals in col c_i to be val.
        /*!
        SPARSE VALUES ARE NOT MODIFIED.
        NOTE: THIS IS A VERY COSTLY OPERATION.     
        DIJ (4-13-22)
        */
        // Iterate through list
        for (int r_i = 0; r_i < row_size; r_i++) {
            std::forward_list<double>::iterator it_v = value_list[r_i].begin();
            std::forward_list<int>::iterator    it_c = col_ind_list[r_i].begin();

            int init_length = std::distance(value_list[r_i].begin(), value_list[r_i].end());

            for (int ii = 0; ii < init_length; ii++){
                // Check row
                if (*it_c == c_i) { *it_v = val; break; }

                ++it_v;
                ++it_c;
            }
        }
    }

    void add_elem(int r_i, int c_i, double val) {
        /*!
        This function adds elements to the COO representation of the sparse matrix.
        When this happens, the CSR or CSC representation (if previously calculated) will 
        be out of date, and will need to be updated. 
        (This is taken care of automatically.)

        Repeat indices are automatically handled, thanks to the vector-of-lists structure
        used.

        DIJ (3-29-22)
        */

        // Update size parameters if expansion occurs
        if(r_i > row_size) { 
            // Add elements to vectors of lists
            for (int ii = 0; ii < (r_i - row_size); ii++) {
                col_ind_list.push_back(std::forward_list<int>()); 
                value_list.push_back(std::forward_list<double>()); 
                }
            row_size = r_i;                    // i.  ROWS
            }
        if(c_i > col_size) { col_size = c_i; } // ii. COLUMNS


        // Iterate through list, inserting new (column, value) pair at the correct
        // location.
        int init_length = std::distance(value_list[r_i].begin(), value_list[r_i].end());

        if (init_length == 0) {
            value_list[r_i].emplace_front(val);
            col_ind_list[r_i].emplace_front(c_i);

            NNZ++;
        }
        else {
            std::forward_list<double>::iterator it_v_Last = value_list[r_i].before_begin();
            std::forward_list<int>::iterator    it_c_Last = col_ind_list[r_i].before_begin();

            std::forward_list<double>::iterator it_v = value_list[r_i].begin();
            std::forward_list<int>::iterator    it_c = col_ind_list[r_i].begin();

            std::forward_list<double>::iterator it_v_Next = value_list[r_i].begin();
            std::forward_list<int>::iterator    it_c_Next = col_ind_list[r_i].begin();
            it_v_Next++; it_c_Next++;

            for(int ii = 0; ii < init_length; ii++){
                // Check row
                if (c_i < *it_c) { 
                    value_list[r_i].insert_after(it_v_Last,   val);
                    col_ind_list[r_i].insert_after(it_c_Last, c_i);

                    NNZ++;
                    break;
                    }
                else if (c_i == *it_c) {
                    // If data is already stored at this location, add new data to original
                    *it_v = *it_v + val;
                    break;
                    }
                else if (c_i > *it_c) {
                    if (ii == (init_length - 1)) {
                        // Final Element, Insert
                        value_list[r_i].insert_after(it_v,   val);
                        col_ind_list[r_i].insert_after(it_c, c_i);

                        NNZ++;
                        break;
                    } else if (c_i < *it_c_Next) {
                        // Next Element NOT c_i == *it_c, Insert
                        value_list[r_i].insert_after(it_v,   val);
                        col_ind_list[r_i].insert_after(it_c, c_i);

                        NNZ++;
                        break;
                    }
                    }
                
                // Update Iterators
                it_v_Last = it_v;
                it_c_Last = it_c;

                ++it_v;
                ++it_c;
            }
        }

        // Update CSR Representation Automatically
        if (active_flag == 2) { convert_to_CSR(); }
    }

    void add_matr(C_Matrix_Dense mat, std::vector<int> row_v, std::vector<int> col_v) {
        //! Add Dense Matrix at Sliced Locations
        /*!
        Stores complete matrix mat at (row, pair) locations given by vectors a and b,
        ADDING to the previously stored value.
        \param mat Dense matrix to be inserted.
        \param row_v vector of row values at which dense matrix should be stored in global sparse matrix
        \param col_v vector of column values "" "".
        */
        for (int ii = 0; ii < row_v.size(); ii++) {
            for (int jj = 0; jj < col_v.size(); jj++) { add_elem(row_v[ii], col_v[jj], mat(ii,jj)); }
        }
    }

    void set_elem(int r_i, int c_i, double val) {
        /*!
        This function sets elements in the COO representation of the sparse matrix,
        OVERWRITING the previously contained value, if present.
        When this happens, the CSR or CSC representation (if previously calculated) will 
        be out of date, and will need to be updated. 
        (This is taken care of automatically.)

        Repeat indices are automatically handled, thanks to the vector-of-lists structure
        used.

        DIJ (3-29-22)
        */

        // Update size parameters if expansion occurs
        if(r_i > row_size) { 
            // Add elements to vectors of lists
            for (int ii = 0; ii < (r_i - row_size); ii++) {
                col_ind_list.push_back(std::forward_list<int>()); 
                value_list.push_back(std::forward_list<double>()); 
                }
            row_size = r_i;                    // i.  ROWS
            }
        if(c_i > col_size) { col_size = c_i; } // ii. COLUMNS


        // Iterate through list, inserting new (column, value) pair at the correct
        // location.
        int init_length = std::distance(value_list[r_i].begin(), value_list[r_i].end());

        if (init_length == 0) {
            value_list[r_i].emplace_front(val);
            col_ind_list[r_i].emplace_front(c_i);

            NNZ++;
        }
        else {
            std::forward_list<double>::iterator it_v_Last = value_list[r_i].before_begin();
            std::forward_list<int>::iterator    it_c_Last = col_ind_list[r_i].before_begin();

            std::forward_list<double>::iterator it_v = value_list[r_i].begin();
            std::forward_list<int>::iterator    it_c = col_ind_list[r_i].begin();

            std::forward_list<double>::iterator it_v_Next = value_list[r_i].begin();
            std::forward_list<int>::iterator    it_c_Next = col_ind_list[r_i].begin();
            it_v_Next++; it_c_Next++;

            for(int ii = 0; ii < init_length; ii++){
                // Check row
                if (c_i < *it_c) { 
                    value_list[r_i].insert_after(it_v_Last,   val);
                    col_ind_list[r_i].insert_after(it_c_Last, c_i);

                    NNZ++;
                    break;
                    }
                else if (c_i == *it_c) {
                    // If data is already stored at this location, overwrite.
                    *it_v = val;
                    break;
                    }
                else if (c_i > *it_c) {
                    if (ii == (init_length - 1)) {
                        // Final Element, Insert
                        value_list[r_i].insert_after(it_v,   val);
                        col_ind_list[r_i].insert_after(it_c, c_i);

                        NNZ++;
                        break;
                    } else if (c_i < *it_c_Next) {
                        // Next Element NOT c_i == *it_c, Insert
                        value_list[r_i].insert_after(it_v,   val);
                        col_ind_list[r_i].insert_after(it_c, c_i);

                        NNZ++;
                        break;
                    }
                    }
                
                // Update Iterators
                it_v_Last = it_v;
                it_c_Last = it_c;

                ++it_v;
                ++it_c;
            }
        }

        // Update CSR Representation Automatically
        if (active_flag == 2) { convert_to_CSR(); }
    }

    void set_matr(C_Matrix_Dense mat, std::vector<int> row_v, std::vector<int> col_v) {
        //! Set Dense Matrix at Sliced Locations
        /*!
        Stores complete matrix mat at (row, pair) locations given by vectors a and b,
        OVERWRITING the previously contained value.
        \param mat Dense matrix to be inserted.
        \param row_v vector of row values at which dense matrix should be stored in global sparse matrix
        \param col_v vector of column values "" "".
        */
        for (int ii = 0; ii < row_v.size(); ii++) {
            for (int jj = 0; jj < col_v.size(); jj++) { set_elem(row_v[ii], col_v[jj], mat(ii,jj)); }
        }
    }


    // III. Conversions
    void convert_to_CSR() {
        /*!
        Convert array to CSR format.

        DIJ (4-15-22)
        */
        active_flag = 2; // CSR Format

        value.resize(NNZ,   0); 
        row_ind.resize(row_size+1, 0); 
        col_ind.resize(NNZ, 0);

        int kk = 0;
        for (int ii = 0; ii <= row_size; ii++){
            std::forward_list<double>::iterator it_v = value_list[ii].begin();
            std::forward_list<int>::iterator    it_c = col_ind_list[ii].begin();
            
            int init_length = std::distance(value_list[ii].begin(), value_list[ii].end());

            for(int jj = 0; jj < init_length; jj++){
                value[kk]   = *it_v; 
                col_ind[kk] = *it_c;

                row_ind[ii + 1]++;

                kk += 1;
                ++it_v;
                ++it_c;
            }

            row_ind[ii + 1] += row_ind[ii]; 
        }
    }

    void convert_to_Eigen(Eigen::SparseMatrix<double>& kG_eigen) {
        /*!
        Initialize an array suitable for use with Eigen linear solver.

        DIJ (4-19-22)
        */

        std::vector< Eigen::Triplet<double> > temp_COO(NNZ);

        // Iterate over rows
        int kk = 0;
        for (int jj = 0; jj <= row_size; jj++) {
            std::forward_list<double>::iterator it_v = value_list[jj].begin();
            std::forward_list<int>::iterator    it_c = col_ind_list[jj].begin();

            int init_length = std::distance(value_list[jj].begin(), value_list[jj].end());

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


    // IV. Display
    void print_contents() {
        if (active_flag == 1){
            // List
            std::cout << "\n";
            std::cout <<  row_size << " x " << col_size;
            std::cout << " Sparse Matrix stored in COO Form\n";
            std::cout << "(Row, Column, Value)\n";

            // Iterate through list
            for (int ii = 0; ii <= row_size; ii++) {
                std::forward_list<double>::iterator it_v = value_list[ii].begin();
                std::forward_list<int>::iterator    it_c = col_ind_list[ii].begin();
                
                // Skip, if no elements are contained in this row
                int init_length = std::distance(value_list[ii].begin(), value_list[ii].end());
                if (init_length < 1) { continue; }

                for(int jj = 0; jj < init_length; jj++){
                    std::cout << "(" << ii << "," << *it_c << "," << *it_v << ")\n";
                    ++it_v;
                    ++it_c;
                }
            }

        } else if (active_flag == 2){
            // CSR
            std::cout << "\n";
            std::cout <<  row_size << " x " << col_size;
            std::cout << " Sparse Matrix stored in CSR Form\n";
            std::cout << "(Row, Column, Value)\n";
            for (int ii = 0; ii <= row_size; ii++){
                int row_start = row_ind[ii];
                int row_end   = row_ind[ii + 1];

                std::vector<int>    col_slice = std::vector<int>(col_ind.begin()  + row_start, col_ind.begin() + row_end);
                std::vector<double> val_slice = std::vector<double>(value.begin() + row_start, value.begin()   + row_end);

                for (int jj = 0; jj < col_slice.size(); jj++) { 
                    std::cout << "(" << ii << "," << col_slice[jj] << ",";
                    std::cout << val_slice[jj] << ")\n";
                    }
            }
        } 
        
        }

};

//! EXTERNAL FUNCTIONS
//! i. Operator Overload for Output 
/*! 
This overload is performed exterior to the C_Matrix_... class so that it can be 
accessed by the std::ostream class and standard syntax can be employed.
*/
std::ostream& operator<<(std::ostream& os, C_Matrix_Sparse& obj) {
    double output_width = 6;

    os << "\n";
    os << obj.row_size << " x " << obj.col_size;
    os << " Sparse Matrix\n";

    // Iterate through list
    for (int ii = 0; ii < obj.row_size; ii++) {
        std::forward_list<double>::iterator it_v = obj.value_list[ii].begin();
        std::forward_list<int>::iterator    it_c = obj.col_ind_list[ii].begin();
        
        // Skip, if no elements are contained in this row
        int init_length = std::distance(obj.value_list[ii].begin(), obj.value_list[ii].end());
        int num_val_row = 0;

        for(int jj = 0; jj < obj.col_size; jj++){
            if (num_val_row < init_length) {
                // 1. Iterator contains additional elements.
                if ( *it_c == jj ) {
                    // i. Element is contained at this position
                    os << std::setw(output_width) << *it_v << " "; 

                    ++it_c; ++it_v;
                    num_val_row++;
                }
                else {
                    // ii. Zero element here, empty
                    os << std::setw(output_width) << 0.0 << " "; 
                }
            }
            else {
                // 2. No remaining elements in iterator
                // Zero element here, empty
               os << std::setw(output_width) << 0.0 << " "; 
            }
        }
        os << "\n";
    }

    return os;
}


#endif