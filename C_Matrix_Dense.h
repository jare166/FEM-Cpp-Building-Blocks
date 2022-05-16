#ifndef C_MATRIX_DENSE_H
#define C_MATRIX_DENSE_H

#include <iostream>
#include <iomanip>
#include <exception>

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
class C_Matrix_Dense{
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
    //! Access Operation
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
        if ( r_i > row_1 ) { 
            std::string error_message = "Requested row index ";
            error_message = error_message + std::to_string(r_i) + " > total row size " + std::to_string(row_1);
            error_message = error_message + ". Error in: " + __FILE__ + ", at line " + std::to_string(__LINE__) + ".";

            throw std::out_of_range(error_message); 
            }
        if ( c_i > col_1 ) {
            std::string error_message = "Requested column index ";
            error_message = error_message + std::to_string(c_i) + " > total column size " + std::to_string(col_1);
            error_message = error_message + ". Error in: " + __FILE__ + ", at line " + std::to_string(__LINE__) + ".";
            
            throw std::out_of_range(error_message); 
        }

        // Determine access index for vector of values
        if (!(*this).transpose) { ij = r_i*col_1 + c_i; }
        else                    { ij = c_i*row_1 + r_i; }

        // Access value assigned
        double& access_val = (*this).values[ij];
        return access_val;
        }

    //! Assignment Operator, from Dense Matrix
    C_Matrix_Dense operator=(C_Matrix_Dense obj2) {
        // 1. Assignment to Empty Object
        if ( (*this).row_size == -1 ) { 
            // Update Size
            if (!obj2.transpose) {
                // i. Input Object not Transposed
                (*this).row_size = obj2.row_size;
                (*this).col_size = obj2.col_size;
            }
            else {
                // ii. Input Object Transposed
                (*this).row_size = obj2.col_size;
                (*this).col_size = obj2.row_size;
            }

            (*this).NNZ = obj2.NNZ;
            (*this).values.resize(NNZ);
        }

        // 2. EXCEPTION: Object Size Mismatch
        if (!check_size((*this), obj2)) { 
            std::string error_message = "Assignment failed. Dense matrix size mismatch.";
            error_message = error_message + ". Error in: " + __FILE__ + ", at line " + std::to_string(__LINE__) + ".";

            throw std::length_error(error_message); 
        }

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
    C_Matrix_Dense operator*(C_Matrix_Dense obj2) {
        /*!
        Overloads * operator to allow two matrices to be multiplied.

        Based on Dr. Srinivasa's code.

        DIJ (5-04-22)
        */

        // i. Check if inner dimensions match
        int ii_M, jj_M, kk_M; // Bounds for Loops
        if (!check_size((*this), obj2, ii_M, jj_M, kk_M)) { 
            std::string error_message = "Multiplication failed. Inner dimensions do not match.";
            error_message = error_message + ". Error in: " + __FILE__ + ", at line " + std::to_string(__LINE__) + ".";
            
            throw std::length_error(error_message); 
            }

        // Initialize object to appropriate size
        C_Matrix_Dense obj_out(ii_M, jj_M); 

        // ii. Otherwise, compute matrix product
        double sum;
        int ij;

        for (int ii = 0; ii < ii_M; ii++) {
            for (int jj = 0; jj < jj_M; jj++) {
                sum=0;
                for(int kk = 0; kk < kk_M; kk++) {
                    sum += (*this)(ii,kk)*obj2(kk,jj);
                }
                ij = ii*jj_M + jj;
                obj_out.values[ij] = sum;
            }
        }
        return obj_out;
    }
    
    //! Matrix Addition
    C_Matrix_Dense operator+(C_Matrix_Dense obj2) {
        // Check: Ensure both matrices are the same size
        if (!((*this).row_size == obj2.row_size)) {
            std::string error_message = "Addition failed. Row size mismatch.";
            error_message = error_message + ". Error in: " + __FILE__ + ", at line " + std::to_string(__LINE__) + ".";
            
            throw std::length_error(error_message);
            }
        if (!((*this).col_size == obj2.col_size)) {
            std::string error_message = "Addition failed. Column size mismatch.";
            error_message = error_message + ". Error in: " + __FILE__ + ", at line " + std::to_string(__LINE__) + ".";
            
            throw std::length_error(error_message);
            }

        C_Matrix_Dense obj_out((*this).row_size, (*this).col_size);
        for (int ii = 0; ii < (*this).row_size; ii++) {
            for (int jj = 0; jj < (*this).col_size; jj++) { obj_out(ii,jj) = (*this)(ii,jj) + obj2(ii,jj); }
        }
        return obj_out;
    }
    //! Matrix Subtraction
    C_Matrix_Dense operator-(C_Matrix_Dense obj2) {
        // Check: Ensure both matrices are the same size
        if (!((*this).row_size == obj2.row_size)) { throw std::length_error("Addition failed. Row size mismatch."); }
        if (!((*this).col_size == obj2.col_size)) { throw std::length_error("Addition failed. Column size mismatch."); }

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
    double inner_product(C_Matrix_Dense obj2) {
        if ((*this).NNZ != obj2.NNZ) { 
            std::string error_message = "Inner product failed. NNZ elements mismatch.";
            error_message = error_message + ". Error in: " + __FILE__ + ", at line " + std::to_string(__LINE__) + ".";
            
            throw std::length_error(error_message);
            }

        double in_prod = 0.0;
        for (int ii = 0; ii < (*this).row_size; ii++) {
            for (int jj = 0; jj < (*this).col_size; jj++) { in_prod += obj2(ii,jj) * (*this)(ii,jj); }
        }

        return in_prod;
    }


    // III. Utility Functions
    bool check_size(C_Matrix_Dense obj1, C_Matrix_Dense obj2) {
        //! Check that matrix sizes match before operation
        /*!
        This function checks size of two matrices intended for an operation
        (such as +, -, or *) to ensure that a result can be computed.

        If one or more of the objects is transposed, then the dimensions to be 
        compared must be changed accordingly; this is because the transpose 
        operation DOES NOT rearrange the contents of the matrix, but merely 
        triggers a flag indicating that the object should be treated differently.

        \param obj1 First matrix to be compared
        \param obj2 Second matrix to be compared

        \author Dominic Jarecki
        \date 5-10-22 
        */

        int row_1 = (*this).row_size; // first matrix
        int col_1 = (*this).col_size;
        int row_2 = obj2.row_size;    // second matrix
        int col_2 = obj2.col_size;

        if (col_1 != row_2) { return false; }

        return true;
    }
    
    bool check_size(C_Matrix_Dense obj1, C_Matrix_Dense obj2, int& ii_M, int& jj_M, int& kk_M) {
        //! Check that matrix sizes match before operation, AND return indices for use in sizing/multiplication
        /*!
        This function checks size of two matrices intended for an operation
        (such as +, -, or *) to ensure that a result can be computed.

        If one or more of the objects is transposed, then the dimensions to be 
        compared must be changed accordingly; this is because the transpose 
        operation DOES NOT rearrange the contents of the matrix, but merely 
        triggers a flag indicating that the object should be treated differently.

        Additionally, this function will initialize indices ii_M, jj_M, and kk_M 
        which can be used for multiplication, etc. (This will save on the number 
        of conditional statements that need to be invoked.)

        \param obj1 First matrix to be compared
        \param obj2 Second matrix to be compared
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

        if (col_1 != row_2) { return false; }

        ii_M = row_1; // Left  Outer Dimension
        jj_M = col_2; // Right Outer Dimension
        kk_M = col_1; // Inner Dimension

        return true;
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

#endif