# FEM-CPP-Building-Blocks


Introduction
------------
This repository contains C++ classes intended to allow rapid building and construction of Finite Element code. Implementation is split into two main classes:
 * **C_Matrix_Dense** - Stores a matrix as a vector in Row-Major format. Supports all standard matrix manipulations, including addition, subtraction, matrix and scalar multiplication and transposition. Size is determined at initialization/assignment; dynamic resize is *not allowed*. Use `C_Matrix_Dense` to construct individual elements.
 * **C_Matrix_Sparse** - Stores a matrix as a vector of linked lists, retaining only non-zero (and accessed) elements. By default, stores in COO format; use this for initialization of the object. If needed, convert to CSR representation, by calling `C_Matrix_Sparse::convert_to_CSR()`. Use `C_Matrix_Sparse` to assemble element matrices into a lightweight sparse matrix in formats suitable for use with different solver libraries.

Overloaded output-stream (<<) operator allows user-friendly viewing of matrix contents for each class. `C_Matrix_Sparse::print_contents()` displays non-zero elements only.

Installation
------------
`FEM-CPP-Building-Blocks` is a header-only library. Installation should be as simple as cloning the repository and using a suitable compiler. Currently, this repository is tested and known to work with `g++17`.
`FEM-CPP-Building-Blocks` is designed to be used as for construction of solver-agnostic objects, flexible enough to grow with a finite element research project. [Solver interfaces][Solver_LINK] are stored separately

Implementation
------------


[Solver_LINK]: 
