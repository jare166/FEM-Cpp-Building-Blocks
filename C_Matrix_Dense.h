#ifndef C_MATRIX_DENSE_H
#define C_MATRIX_DENSE_H

#include <iostream>
#include <iomanip>
#include <exception>

#include "f_MiscellaneousFunctions.h"

// Report File and Line Number
#define throw_line(arg) throw my_exception(arg, __FILE__, __LINE__);


//! This code initializes dense matrix objects for use in constructing global FEM matrices.
/*!

Contents are stored in row-major format, so the default indexing scheme for pulling from the main
vector of values is, for (i,j)
    NOT TRANSPOSED: i*col_size + j;
    TRANSPOSED:     i*row_size + j;

Because the parentheses operator () has been overloaded to take care of this indexing scheme 
automatically, all code (where possible) is written in terms of this.

Note that when accessing elements from within the object itself, a slightly ambiguous syntax is used;
namely:
    (*this)(ii,jj)

Author: Dominic Jarecki
Date: 5-03-2022
*/
class C_Matrix_Dense {
    public:
        std::vector<double> values;        //! All matrix values stored sequentially in vector
        int  row_size = -1, col_size = -1; //! Fixed size, not allowed to grow.
        int  NNZ;                          //! Total number of elements
        bool transpose = false;            //! Object currently transposed?

    // I. Constructor
    C_Matrix_Dense() {} 

    C_Matrix_Dense(int row_size_in, int col_size_in) {
        row_size = row_size_in;
        col_size = col_size_in;
        NNZ = row_size*col_size;

        values.resize(NNZ);
        for (int ii = 0; ii < NNZ; ii++) { values[ii] = 0.0; }
    }


    // II. Operator Overloads
    //! 1. Basic Access Operation
    double& operator() (int r_i, int c_i) { 
        //! Access Contents of Matrix
        /*! 
        This function returns elements contained in C_Matrix_Dense. Passed by
        reference so that they can be modified through this method also.

        Contents are stored in row-major format!

        DIJ (5-04-22)
        */

        int row_1 = (*this).row_size; // first matrix
        int col_1 = (*this).col_size;
        int ij; 

        // EXCEPTION: Check for out-of-range access
        check_size_access((*this), r_i, c_i);

        // Determine access index for vector of values
        if (!(*this).transpose) { ij = r_i*col_1 + c_i; }
        else                    { ij = c_i*row_1 + r_i; }

        // Access value assigned
        double& access_val = (*this).values[ij];
        return access_val;
    }
    
    //! 2. Slice ACCESS Operations --> NOT PASSED BY REFERENCE
    //!     i.   Slice: (row, col_range)
    C_Matrix_Dense operator() (int r_i, std::vector<int> c_i) { 
        //! Access Contents of Matrix
        /*! 
        This function returns a slice of elements contained in C_Matrix_Dense. NOT passed by
        reference, so the values within the matrix itself cannot be modified directly through this method.

        Contents are stored in row-major format!

        DIJ (5-17-22)
        */

        int row_1 = (*this).row_size; // first matrix
        int col_1 = (*this).col_size;

        int ij; 

        // EXCEPTION: Check for out-of-range access
        check_size_access((*this), r_i, c_i);

        // Initialize Output to Correct Size
        int num_rows = 1;
        int num_cols = c_i.size();
        C_Matrix_Dense access_slice(num_rows, num_cols);

        // Determine access index for vector of values
        for (int jj = 0; jj < num_cols; jj++) {
            if (!(*this).transpose) { ij = r_i*col_1 + c_i[jj]; }
            else                    { ij = c_i[jj]*row_1 + r_i; }

            access_slice(0,jj) = (*this).values[ij];
        }

        // Access value assigned
        return access_slice;
    }
    //!     ii.  Slice: (row_range, col)
    C_Matrix_Dense operator() (std::vector<int> r_i, int c_i) { 
        //! Access Contents of Matrix
        /*! 
        This function returns a slice of elements contained in C_Matrix_Dense. NOT passed by
        reference, so the values within the matrix itself cannot be modified directly through this method.

        Contents are stored in row-major format!

        DIJ (5-17-22)
        */

        int row_1 = (*this).row_size; // first matrix
        int col_1 = (*this).col_size;

        int ij; 

        // EXCEPTION: Check for out-of-range access
        check_size_access((*this), r_i, c_i);

        // Initialize Output to Correct Size
        int num_rows = r_i.size();
        int num_cols = 1;
        C_Matrix_Dense access_slice(num_rows, num_cols);

        // Determine access index for vector of values
        for (int ii = 0; ii < num_rows; ii++) {
            if (!(*this).transpose) { ij = r_i[ii]*col_1 + c_i; }
            else                    { ij = c_i*row_1 + r_i[ii]; }

            access_slice(ii,0) = (*this).values[ij];
        }

        // Access value assigned
        return access_slice;
    }
    //!     iii. Slice: (row_range, col_range)
    C_Matrix_Dense operator() (std::vector<int> r_i, std::vector<int> c_i) { 
        //! Access Contents of Matrix
        /*! 
        This function returns a slice of elements contained in C_Matrix_Dense. NOT passed by
        reference, so the values within the matrix itself cannot be modified directly through this method.

        Contents are stored in row-major format!

        DIJ (5-17-22)
        */

        int row_1 = (*this).row_size; // first matrix
        int col_1 = (*this).col_size;

        int ij; 

        // EXCEPTION: Check for out-of-range access
        check_size_access((*this), r_i, c_i);

        // Initialize Output to Correct Size
        int num_rows = r_i.size();
        int num_cols = c_i.size();
        C_Matrix_Dense access_slice(num_rows, num_cols);

        // Determine access index for vector of values
        for (int ii = 0; ii < num_rows; ii++) {
            for (int jj = 0; jj < num_cols; jj++) {
                if (!(*this).transpose) { ij = r_i[ii]*col_1 + c_i[jj]; }
                else                    { ij = c_i[jj]*row_1 + r_i[ii]; }

                access_slice(ii,jj) = (*this).values[ij];
            }
        }

        // Access value assigned
        return access_slice;
    }
    
    //! 3. Slice MODIFICATION Operations --> NO VALUE RETURNED
    //! add(): ADD TO PREVIOUSLY-STORED DATA
    //!     i. Add element at index
    void add_elem(double val, int r_i, int c_i) {
        /*!
        This function adds elements to the dense matrix; operations here are 
        redundant w/ parentheses operator ().

        \param val value which should be stored
        \param r_i row value at which value should be stored in matrix
        \param c_i column value at which value should be stored in matrix

        DIJ (5-18-22)
        */

        int row_1 = (*this).row_size; // first matrix
        int col_1 = (*this).col_size;

        int ij; 

        // EXCEPTION: Check for out-of-range access
        check_size_access((*this), r_i, c_i);

        // Determine access index for vector of values
        if (!(*this).transpose) { ij = r_i*col_1 + c_i; }
        else                    { ij = c_i*row_1 + r_i; }

        (*this).values[ij] += val;
    }
    //!     ii. Add element at slices
    void add_matr(double mat_val, std::vector<int> r_i, std::vector<int> c_i) {
        //! Add Dense Matrix at Sliced Locations
        /*!
        Stores complete matrix mat at (row, pair) locations given by vectors a and b,
        ADDING to the previously stored value.

        \param mat_val Dense matrix value to be inserted.
        \param r_i vector of row values at which dense matrix should be stored in global sparse matrix
        \param c_i vector of column values "" "".

        DIJ (5-18-22)
        */

        // EXCEPTION: Check for out-of-range access
        // Checked in sub-function call

        // Initialize Output to Correct Size
        int num_rows = r_i.size();
        int num_cols = c_i.size();

        // Determine access index for vector of values
        for (int ii = 0; ii < num_rows; ii++) {
            for (int jj = 0; jj < num_cols; jj++) { add_elem(mat_val, r_i[ii], c_i[jj]); }
        }

    }
    //!     iii. Add element over whole matrix
    void add_matr(double mat_val) {
        //! Add double to entire matrix
        /*!
        Augments contents of matrix with double value mat_val at all points in matrix. 
        Will throw an error if used with an empty matrix.

        \param mat_val Double value to be added.

        DIJ (7-12-22)
        */
        
        // EXCEPTION: Check for empty object
        if ( (*this).row_size == -1 ) { 
            std::string error_message = "Assignment failed. Empty dense matrix object; unable to determine size.";
            error_message = error_message + ". Error in: " + __FILE__ + ", at line " + std::to_string(__LINE__) + ".";

            throw std::length_error(error_message); 
        }

        for (int ii = 0; ii < (*this).NNZ; ii++) { (*this).values[ii] += mat_val; }
    }
    //!     iv. Add matrix at slices
    void add_matr(C_Matrix_Dense mat, std::vector<int> r_i, std::vector<int> c_i) {
        //! Add Dense Matrix at Sliced Locations
        /*!
        Stores complete matrix mat at (row, pair) locations given by vectors a and b,
        ADDING to the previously stored value.

        \param mat Dense matrix to be inserted.
        \param r_i vector of row values at which dense matrix should be stored in global sparse matrix
        \param c_i vector of column values "" "".

        DIJ (5-18-22)
        */

        // EXCEPTION: Check for out-of-range access
        // Checked in sub-function call

        // Initialize Output to Correct Size
        int num_rows = r_i.size();
        int num_cols = c_i.size();

        // Determine access index for vector of values
        for (int ii = 0; ii < num_rows; ii++) {
            for (int jj = 0; jj < num_cols; jj++) { add_elem(mat(ii, jj), r_i[ii], c_i[jj]); }
        }
        
    }
    
    //! set(): OVER-WRITE PREVIOUSLY-STORED DATA
    //!     i. Set element at index
    void set_elem(const double val, int r_i, int c_i) {
        /*!
        This function sets elements in the dense matrix,
        OVERWRITING the previously contained value, if present.

        \param val value which should be stored
        \param r_i row value at which value should be stored in matrix
        \param c_i column value at which value should be stored in matrix

        DIJ (5-18-22)
        */

        int row_1 = (*this).row_size; // first matrix
        int col_1 = (*this).col_size;

        int ij; 

        // EXCEPTION: Check for out-of-range access
        check_size_access((*this), r_i, c_i);

        // Determine access index for vector of values
        if (!(*this).transpose) { ij = r_i*col_1 + c_i; }
        else                    { ij = c_i*row_1 + r_i; }

        (*this).values[ij] = val;
    }
    //!     ii. Set element at slices
    void set_matr(double mat_val, std::vector<int> r_i, std::vector<int> c_i) {
        //! Set Double at Sliced Locations
        /*!
        Stores double mat_val at (row, pair) locations given by vectors a and b,
        OVERWRITING the previously contained value.

        \param mat_val Double value to be inserted.
        \param r_i vector of row values at which dense matrix should be stored in global sparse matrix
        \param c_i vector of column values "" "".

        DIJ (5-18-22)
        */
        
        // EXCEPTION: Check for out-of-range access
        // Checked in sub-function call

        // Initialize Output to Correct Size
        int num_rows = r_i.size();
        int num_cols = c_i.size();

        // Determine access index for vector of values
        for (int ii = 0; ii < num_rows; ii++) {
            for (int jj = 0; jj < num_cols; jj++) { set_elem(mat_val, r_i[ii], c_i[jj]); }
        }

    }
    //!     iii. Set element over whole matrix
    void set_matr(double mat_val) {
        //! Set double over entire matrix
        /*!
        Overwrites contents of matrix with double value mat_val at all points in matrix. 
        Will throw an error if used with an empty matrix.

        \param mat_val Double value to be inserted.

        DIJ (7-12-22)
        */
        
        // EXCEPTION: Check for empty object
        if ( (*this).row_size == -1 ) { 
            std::string error_message = "Assignment failed. Empty dense matrix object; unable to determine size.";
            error_message = error_message + ". Error in: " + __FILE__ + ", at line " + std::to_string(__LINE__) + ".";

            throw std::length_error(error_message); 
        }

        for (int ii = 0; ii < (*this).NNZ; ii++) { (*this).values[ii] = mat_val; }
    }
    //!     iv. Set matrix at slices
    void set_matr(C_Matrix_Dense& mat, std::vector<int> r_i, std::vector<int> c_i) {
        //! Set Dense Matrix at Sliced Locations
        /*!
        Stores complete matrix mat at (row, pair) locations given by vectors a and b,
        OVERWRITING the previously contained value.

        \param mat Dense matrix to be inserted.
        \param r_i vector of row values at which dense matrix should be stored in global sparse matrix
        \param c_i vector of column values "" "".

        DIJ (5-18-22)
        */

        // EXCEPTION: Check for out-of-range access
        // Checked in sub-function call

        // Initialize Output to Correct Size
        int num_rows = r_i.size();
        int num_cols = c_i.size();

        // Determine access index for vector of values
        for (int ii = 0; ii < num_rows; ii++) {
            for (int jj = 0; jj < num_cols; jj++) { set_elem(mat(ii, jj), r_i[ii], c_i[jj]); }
        }
    }
    //!     v. Set matrix/vector to zero values, no copy
    void setZero() {
        // 1. Assignment to Empty Object: Unable to determine size
        if ( (*this).row_size == -1 ) { 
            std::string error_message = "Assignment failed. Empty dense matrix object; unable to determine size.";
            error_message = error_message + ". Error in: " + __FILE__ + ", at line " + std::to_string(__LINE__) + ".";

            throw std::length_error(error_message); 
        }

        // 2. Assign Values
        for (int ii = 0; ii < (*this).NNZ; ii++) { (*this).values[ii] = 0.0; }
    }

    //! Assignment Operator, from Dense Matrix
    C_Matrix_Dense operator=(C_Matrix_Dense obj2) {
        // 1. Assignment to Empty Object
        if ( (*this).row_size == -1 ) { 
            // Update Size
            (*this).row_size = obj2.row_size;
            (*this).col_size = obj2.col_size;

            (*this).NNZ = obj2.NNZ;
            (*this).values.resize(NNZ);
        }

        // 2. EXCEPTION: Object Size Mismatch
        check_size_addition((*this), obj2, "Assignment");

        // 2. Object Sizes Match, Assign Values
        for (int ii = 0; ii < (*this).row_size; ii++) {
            for (int jj = 0; jj < (*this).col_size; jj++) { (*this)(ii,jj) = obj2(ii,jj); }
        }
        return (*this);
    }
    //! Assignment Operator, from Array
    C_Matrix_Dense operator=(double *obj2) {
        // 1. Assignment to Empty Object: Unable to determine size
        if ( (*this).row_size == -1 ) { 
            std::string error_message = "Assignment failed. Empty dense matrix object; unable to determine size.";
            error_message = error_message + ". Error in: " + __FILE__ + ", at line " + std::to_string(__LINE__) + ".";

            throw std::length_error(error_message); 
            }

        // 2. Assume Object Sizes Match, Assign Values
        for (int ii = 0; ii < (*this).NNZ; ii++) { (*this).values[ii] = obj2[ii]; }
        
        return (*this);
    }
    //! Assignment Operator, from Vector
    C_Matrix_Dense operator=(std::vector<double> obj2) {
        // 1. Assignment to Empty Object: Unable to determine size
        if ( (*this).row_size == -1 ) { 
            std::string error_message = "Assignment failed. Empty dense matrix object; unable to determine size.";
            error_message = error_message + ". Error in: " + __FILE__ + ", at line " + std::to_string(__LINE__) + ".";
            
            throw std::length_error(error_message); 
        }

        // 2. Assume Object Sizes Match, Assign Values
        for (int ii = 0; ii < (*this).NNZ; ii++) { (*this).values[ii] = obj2[ii]; }
        
        return (*this);
    }

    //! Matrix Multiplication
    C_Matrix_Dense operator*(C_Matrix_Dense& obj2) {
        /*!
        Overloads * operator to allow two matrices to be multiplied.

        Based on Dr. Srinivasa's code.

        DIJ (5-04-22)
        */

        // i. Check if inner dimensions match
        int ii_M, jj_M, kk_M; // Bounds for Loops     
        check_size_multiplication((*this), obj2, "Multiplication", ii_M, jj_M, kk_M);

        // Initialize object to appropriate size
        C_Matrix_Dense obj_out(ii_M, jj_M); 

        // ii. Otherwise, compute matrix product
        double sum;
        int ij;

        for (int ii = 0; ii < ii_M; ii++) {
            for (int jj = 0; jj < jj_M; jj++) {
                sum=0;
                for(int kk = 0; kk < kk_M; kk++) { sum += (*this)(ii,kk)*obj2(kk,jj); }
                ij = ii*jj_M + jj;
                obj_out.values[ij] = sum;
            }
        }
        return obj_out;
    }
    
    //! Matrix Addition
    C_Matrix_Dense operator+(C_Matrix_Dense& obj2) {
        // Check: Ensure both matrices are the same size
        check_size_addition((*this), obj2, "Addition");

        C_Matrix_Dense obj_out((*this).row_size, (*this).col_size);
        for (int ii = 0; ii < (*this).row_size; ii++) {
            for (int jj = 0; jj < (*this).col_size; jj++) { obj_out(ii,jj) = (*this)(ii,jj) + obj2(ii,jj); }
        }
        return obj_out;
    }
    //! Matrix Subtraction
    C_Matrix_Dense operator-(C_Matrix_Dense& obj2) {
        // Check: Ensure both matrices are the same size
        check_size_addition((*this), obj2, "Subtraction");

        C_Matrix_Dense obj_out((*this).row_size, (*this).col_size);
        for (int ii = 0; ii < (*this).row_size; ii++) {
            for (int jj = 0; jj < (*this).col_size; jj++) { obj_out(ii,jj) = (*this)(ii,jj) - obj2(ii,jj); }
        }
        return obj_out;
    }
    
    //! Return copy of object containing matrix negation
    C_Matrix_Dense operator-() {
        C_Matrix_Dense obj_out = (*this);
        for (int ii = 0; ii < (*this).NNZ; ii++) { obj_out.values[ii] = -obj_out.values[ii]; }
        return obj_out;
    }
    //! Return copy of object flagged as Transposed
    C_Matrix_Dense T () {
        // Copy Matrix
        C_Matrix_Dense obj_out = (*this);

        // Redefine Indices
        obj_out.row_size = (*this).col_size;
        obj_out.col_size = (*this).row_size;
        obj_out.transpose = true;

        return obj_out;
    }

    //! Inner Product
    double inner_product(const C_Matrix_Dense& obj2) {
        if ((*this).NNZ != obj2.NNZ) { 
            std::string error_message = "Inner product failed. NNZ elements mismatch.";
            error_message = error_message + ". Error in: " + __FILE__ + ", at line " + std::to_string(__LINE__) + ".";
            
            throw std::length_error(error_message);
        }

        double in_prod = 0.0;
        for (int ii = 0; ii < (*this).NNZ; ii++) { in_prod += obj2.values[ii] * (*this).values[ii];  }

        return in_prod;
    }
    //! L2 Vector Norm
    double norm() {
        return std::sqrt( inner_product((*this)) );
    }

    // III. Utility Functions
    //! ERROR-CHECK UTILITIES
    //! 1. Check Access
    //!     i. Check Access (input: index, index)
    void check_size_access(const C_Matrix_Dense& obj1, int r_i, int c_i) {
        //! Check that matrix sizes match before addition, etc.
        /*!
        This function checks matrix having its value accessed to ensure no
        out-of-bounds accesses are permitted.

        \param obj1 Matrix to have value accessed
        \param r_i row index or set of indices ("slice")
        \param c_i column index or set of indices ("slice")

        \author Dominic Jarecki
        \date 5-18-22 
        */

        // Check row access
        if ( r_i >= obj1.row_size ) { 
            std::string error_message = "Matrix access failed. Requested row index ";
            error_message = error_message + std::to_string(r_i) + " > total row size " + std::to_string(obj1.row_size-1);
            error_message = error_message + ". Error in: " + __FILE__ + ", at line " + std::to_string(__LINE__) + ".";

            throw std::out_of_range(error_message); 
        }
        
        // Check column access
        if ( c_i >= obj1.col_size ){
            std::string error_message = "Matrix access failed. Requested column index ";
            error_message = error_message + std::to_string(c_i) + " > total column size " + std::to_string(obj1.col_size-1);
            error_message = error_message + ". Error in: " + __FILE__ + ", at line " + std::to_string(__LINE__) + ".";

            throw std::out_of_range(error_message); 
        }
        
    }
    //!     ii. Check Access (input: index, slice)
    void check_size_access(const C_Matrix_Dense& obj1, int r_i, std::vector<int> c_i) {
        //! Check that matrix sizes match before addition, etc.
        /*!
        This function checks matrix having its value accessed to ensure no
        out-of-bounds accesses are permitted.

        \param obj1 Matrix to have value accessed
        \param r_i row index or set of indices ("slice")
        \param c_i column index or set of indices ("slice")

        \author Dominic Jarecki
        \date 5-18-22 
        */

        // Check row access
        if ( r_i >= obj1.row_size ) { 
            std::string error_message = "Matrix access failed. Requested row index ";
            error_message = error_message + std::to_string(r_i) + " > total row size " + std::to_string(obj1.row_size-1);
            error_message = error_message + ". Error in: " + __FILE__ + ", at line " + std::to_string(__LINE__) + ".";

            throw std::out_of_range(error_message); 
        }

        // Check column access
        for (int ii = 0; ii < c_i.size(); ii++) {
            if ( c_i[ii] >= obj1.col_size ){
                std::string error_message = "Matrix access failed. Requested column index ";
                error_message = error_message + std::to_string(c_i[ii]) + " > total column size " + std::to_string(obj1.col_size-1);
                error_message = error_message + ". Error in: " + __FILE__ + ", at line " + std::to_string(__LINE__) + ".";

                throw std::out_of_range(error_message); 
            }
        }
    }
    //!     iii. Check Access (input: slice, index)
    void check_size_access(const C_Matrix_Dense& obj1, std::vector<int> r_i, int c_i) {
        //! Check that matrix sizes match before addition, etc.
        /*!
        This function checks matrix having its value accessed to ensure no
        out-of-bounds accesses are permitted.

        \param obj1 Matrix to have value accessed
        \param r_i row index or set of indices ("slice")
        \param c_i column index or set of indices ("slice")

        \author Dominic Jarecki
        \date 5-18-22 
        */

        // Check row access
        for (int ii = 0; ii < r_i.size(); ii++) {
            if ( r_i[ii] >= obj1.row_size ) { 
                std::string error_message = "Matrix access failed. Requested row index ";
                error_message = error_message + std::to_string(r_i[ii]) + " > total row size " + std::to_string(obj1.row_size-1);
                error_message = error_message + ". Error in: " + __FILE__ + ", at line " + std::to_string(__LINE__) + ".";

                throw std::out_of_range(error_message); 
            }
        }

        // Check column access
        if ( c_i >= obj1.col_size ){
            std::string error_message = "Matrix access failed. Requested column index ";
            error_message = error_message + std::to_string(c_i) + " > total column size " + std::to_string(obj1.col_size-1);
            error_message = error_message + ". Error in: " + __FILE__ + ", at line " + std::to_string(__LINE__) + ".";

            throw std::out_of_range(error_message); 
        }
        
    }
    //!     iv. Check Access (input: slice, slice)
    void check_size_access(const C_Matrix_Dense& obj1, std::vector<int> r_i, std::vector<int> c_i) {
        //! Check that matrix sizes match before addition, etc.
        /*!
        This function checks matrix having its value accessed to ensure no
        out-of-bounds accesses are permitted.

        \param obj1 Matrix to have value accessed
        \param r_i row index or set of indices ("slice")
        \param c_i column index or set of indices ("slice")

        \author Dominic Jarecki
        \date 5-18-22 
        */

        // Check row access
        for (int ii = 0; ii < r_i.size(); ii++) {
            if ( r_i[ii] >= obj1.row_size ) { 
                std::string error_message = "Matrix access failed. Requested row index ";
                error_message = error_message + std::to_string(r_i[ii]) + " > total row size " + std::to_string(obj1.row_size-1);
                error_message = error_message + ". Error in: " + __FILE__ + ", at line " + std::to_string(__LINE__) + ".";

                throw std::out_of_range(error_message); 
            }
        }

        // Check column access
        for (int ii = 0; ii < c_i.size(); ii++) {
            if ( c_i[ii] >= obj1.col_size ){
                std::string error_message = "Matrix access failed. Requested column index ";
                error_message = error_message + std::to_string(c_i[ii]) + " > total column size " + std::to_string(obj1.col_size-1);
                error_message = error_message + ". Error in: " + __FILE__ + ", at line " + std::to_string(__LINE__) + ".";

                throw std::out_of_range(error_message); 
            }
        }
    }

    //! 2. Check dimensions for addition, subtration, and asignment
    void check_size_addition(const C_Matrix_Dense& obj1, const C_Matrix_Dense& obj2, std::string op_string) {
        //! Check that matrix sizes match before addition, etc.
        /*!
        This function checks size of two matrices intended for operations such as 
        addition, subtraction, and assignment to ensure that a result can be computed.

        \param obj1 First matrix to be compared
        \param obj2 Second matrix to be compared
        \param op_string String defining operation performed; e.g. "Addition"

        \author Dominic Jarecki
        \date 5-17-22 
        */

        if (!((*this).row_size == obj2.row_size)) {
            std::string error_message = op_string + " failed. Row size mismatch.";
            error_message = error_message + " Error in: " + __FILE__ + ", at line " + std::to_string(__LINE__) + ".";
            
            throw std::length_error(error_message);
        }
        if (!((*this).col_size == obj2.col_size)) {
            std::string error_message = op_string + " failed. Column size mismatch.";
            error_message = error_message + " Error in: " + __FILE__ + ", at line " + std::to_string(__LINE__) + ".";
            
            throw std::length_error(error_message);
        }
    }

    //! 3. Check dimensions for multiplication
    void check_size_multiplication(const C_Matrix_Dense& obj1, const C_Matrix_Dense& obj2, std::string op_string, int& ii_M, int& jj_M, int& kk_M) {
        //! Check that matrix inner dimensions match before multiplication
        /*!
        This function checks size of two matrices intended for a multiplication
        to ensure that a result can be computed.

        Additionally, this function will initialize indices ii_M, jj_M, and kk_M 
        which can be used for multiplication, etc. (This will save on the number 
        of conditional statements that need to be invoked.)

        \param obj1 First matrix to be compared
        \param obj2 Second matrix to be compared
        \param op_string String defining operation performed; e.g. "Multiplication"
        \param ii_M First index to be iterated over in multiplication, etc.
        \param jj_M Second index to be iterated over in multiplication, etc.
        \param kk_M Third index to be iterated over in multiplication, etc.

        \author Dominic Jarecki
        \date 5-10-22 
        */

        int row_1 = (*this).row_size; // first matrix
        int col_1 = (*this).col_size;
        int row_2 = obj2.row_size;    // second matrix
        int col_2 = obj2.col_size;

        // i. Check if inner dimensions match
        // int ii_M, jj_M, kk_M; // Bounds for Loops

        if (col_1 != row_2) {
            std::string error_message = op_string + " failed. Inner dimensions do not match.";
            error_message = error_message + ". Error in: " + __FILE__ + ", at line " + std::to_string(__LINE__) + ".";
            
            throw std::length_error(error_message); 
        }

        ii_M = row_1; // Left  Outer Dimension
        jj_M = col_2; // Right Outer Dimension
        kk_M = col_1; // Inner Dimension
    }
};


//! EXTERNAL FUNCTIONS
//! Operator Overload for Output 
/*! 
This overload is performed exterior to the C_Matrix_... class so that it can be 
accessed by the std::ostream class and standard syntax can be employed.
*/
std::ostream& operator<<(std::ostream& os, C_Matrix_Dense obj) {
    double output_width = 6;

    os << "\n";
    os << obj.row_size << " x " << obj.col_size;
    os << " Dense Matrix\n";

    for (int ii = 0; ii < obj.row_size; ii++) {
        for (int jj = 0; jj < obj.col_size; jj++) { os << std::setw(output_width) << obj(ii,jj) << " "; }
        os << "\n";
    }
    return os;
}

//! Generate Identity Matrix
C_Matrix_Dense identity(int dim) {
    C_Matrix_Dense obj_out(dim, dim);
    for (int ii = 0; ii < dim; ii++) { obj_out.values[(ii*(1+dim))] = 1; }
    return obj_out;
}

//! (Pre-) Scalar Multiplication
C_Matrix_Dense operator* (double scalar, const C_Matrix_Dense& obj2) {
    C_Matrix_Dense obj_out(obj2.row_size, obj2.col_size);
    for (int ii = 0; ii < obj2.NNZ; ii++) { obj_out.values[ii] = scalar*obj2.values[ii]; }
    return obj_out;
}

//! (Post-) Scalar Multiplication
C_Matrix_Dense operator* (const C_Matrix_Dense& obj2, double scalar) {
    C_Matrix_Dense obj_out(obj2.row_size, obj2.col_size);
    for (int ii = 0; ii < obj2.NNZ; ii++) { obj_out.values[ii] = scalar*obj2.values[ii]; }
    return obj_out;
}

//! Division by a Scalar
C_Matrix_Dense operator/ (const C_Matrix_Dense& obj2, double scalar) {
    C_Matrix_Dense obj_out(obj2.row_size, obj2.col_size);
    for (int ii = 0; ii < obj2.NNZ; ii++) { obj_out.values[ii] = obj2.values[ii]/scalar; }
    return obj_out;
}

#endif