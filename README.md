# FEM-CPP-Building-Blocks


Introduction
------------
This repository contains C++ classes intended to allow rapid building and construction of Finite Element code. Implementation is split into two main classes:
 * **C_Matrix_Dense** - Stores a matrix as a vector in [Row-Major format][RMF_LINK]. Supports all standard matrix manipulations, including addition, subtraction, matrix and scalar multiplication and transposition. Size is determined at initialization/assignment; dynamic resize is *not allowed*. Use `C_Matrix_Dense` to construct individual elements.
 * **C_Matrix_Sparse** - Stores a matrix as a vector of linked lists, retaining only non-zero (and accessed) elements. By default, stores in [COO][COO_LINK] format; use this for initialization of the object. If needed, convert to [CSR][CSR_LINK] representation, by calling `C_Matrix_Sparse::convert_to_CSR()`. Use `C_Matrix_Sparse` to assemble element matrices into a lightweight sparse matrix in formats suitable for use with different solver libraries.

Overloaded output-stream (<<) operator allows user-friendly viewing of matrix contents for each class. `C_Matrix_Sparse::print_contents()` displays non-zero elements only.

Installation
------------
`FEM-CPP-Building-Blocks` is a header-only library. Installation should be as simple as cloning the repository and using a suitable compiler. Currently, this repository is tested and known to work with `g++17`.

`FEM-CPP-Building-Blocks` is designed to be used as for construction of solver-agnostic objects, flexible enough to grow with a finite element research project. [Solver interfaces][Solver_LINK] are stored separately so that base classes do not have dependencies on external libraries unless they are used.

Implementation
------------
[Demo scripts][Demo_LINK] are included that show how to construct matrices and display to the terminal. The syntax is simple, and is designed to be MATLAB-like where possible.
>Basic Assignment to Matrices

```c++
C_Matrix_Dense A(2,2) = {1, 2, 3, 4}; // Stored as { {1, 2}; {3, 4} } --> Enter in Row-Major ordering

C_Matrix_Dense B(2,2);
B = -A.T(); // Set B equal to negative of A-transpose.
```

>Write to Sparse Matrix

```c++
C_Matrix_Sparse C(6,6); // Default NNZ = 0;
C.add_matr(B, {3,4}, {2, 3}); // Store Contents of Dense Matrix on rows 4:5 and columns 3:5 (C++ is zero-indexed)
```

>Display to Terminal

```c++
std::cout << A;
std::cout << B;
std::cout << C;
```

Contributors
------------
**Dominic Jarecki**

**Bensingh Dhas**


Acknowledgements
------------
In association with:
[CONTINUUM MECHANICS @ TAMU][CM_LINK]
[Center of Innovation in Mechanics for Design and Manufacturing (CIMDM)][CIMDM_LINK]




[RMF_LINK]: https://en.wikipedia.org/wiki/Row-_and_column-major_order
[COO_LINK]: https://en.wikipedia.org/wiki/Sparse_matrix#Coordinate_list_(COO)
[CSR_LINK]: https://en.wikipedia.org/wiki/Sparse_matrix#Compressed_sparse_row_(CSR,_CRS_or_Yale_format)
[Solver_LINK]: https://github.com/jare166/FEM-Cpp-Building-Blocks/tree/main/Solver_Interfaces
[Demo_LINK]: https://github.com/jare166/FEM-Cpp-Building-Blocks/tree/main/Demo_Scripts
[CM_LINK]: https://continuummechanics.tamu.edu/home/
[CIMDM_LINK]: https://cimdm.tamu.edu/
