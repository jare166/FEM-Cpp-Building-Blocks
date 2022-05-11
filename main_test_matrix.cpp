#include <math.h>
#include <vector>
#include <iostream>
#include <iterator>
#include <list>

#include "C_Matrix_Sparse.h"
#include "C_Matrix_Dense.h"

int main()
{
    // Test 1: Dense Matrix Operations
    C_Matrix_Dense A1(4,4), B1(4,4), C1(4,6);
    C_Matrix_Dense A2(4,3), B2(3,4), C2(3,3);

    // Assignments
    // i. Part I
    A1(0,0) = 4;   A1(0,2) = 3; 
    A1(1,1) = 2.5; A1(1,2) = 4.5;
    A1(2,2) = 3;   A1(2,3) = 1; 
    A1(3,0) = 1.7; A1(3,2) = 0.9;

    B1(0,1) = 4;   B1(0,2) = 2; 
    B1(1,1) = 4.5; B1(1,2) = 1.5; B1(1,3) = 7;
    B1(2,2) = 1;  
    B1(3,0) = 1.2; B1(3,2) = 0.6;

    C1(0,0) = 3;   C1(0,3) = 4; 
    C1(1,1) = 4.5; C1(1,2) = 1.5; C1(1,3) = 2;
    C1(2,3) = 1;    
    C1(3,1) = 1.2; C1(3,2) = 0.6; C1(3,4) = 6.7; C1(3,5) = 2; 
    
    // ii. Part II
    A2(0,0) = 4;   A2(0,2) = 3; 
    A2(1,1) = 2.5; A2(1,2) = 4.5;
    A2(2,2) = 3;   A2(2,1) = 1; 
    A2(3,0) = 1.7; A2(3,2) = 0.9;

    B2(0,1) = 4;   B2(0,2) = 2; 
    B2(1,1) = 4.5; B2(1,2) = 1.5;
    B2(2,0) = 1.2; B2(2,2) = 0.6; B2(2,3) = 9;

    C2(0,0) = 3;   C2(0,2) = 4; 
    C2(1,1) = 4.5; C2(1,2) = 1.5;
    C2(2,1) = 1.2; C2(2,2) = 0.6;

    std::cout << "\nBefore Multiplication: \n";
    std::cout << A1;
    std::cout << B1;
    std::cout << C1;
    std::cout << A2;
    std::cout << B2;
    std::cout << B2.T();
    std::cout << C2;


    // 2. Multiplications
    C_Matrix_Dense D1 = A1*B1;
    C_Matrix_Dense D2 = A1*C1;
    C_Matrix_Dense D3 = A1*B2.T();

    C_Matrix_Dense D4 = B1*B2.T();
    C_Matrix_Dense D5 = C2.T()*A2.T();
    C_Matrix_Dense D6 = C2*A2.T();

    C_Matrix_Dense D7 = B1*(B2.T() + A2);
    C_Matrix_Dense D8 = B1*(B2.T() - A2*5.0);
    C_Matrix_Dense D9 = C2*6.0*A2.T();

    C_Matrix_Dense D10 = C1*C1.T(); 
    C_Matrix_Dense D11 = C1.T()*C1;

    std::cout << "\nAfter Multiplication: \n";
    std::cout << D1;
    std::cout << D2;
    std::cout << D3;
    std::cout << D4;
    std::cout << D5;
    std::cout << D6;
    std::cout << D7;
    std::cout << D8;
    std::cout << D9;
    std::cout << D10;
    std::cout << D11;


    // 3. Trivial Multiplications
    C_Matrix_Dense A(3,3), B(3,3);

    A(0,0) = 3; A(0,1) = 4; A(0,2) = 5;
    A(1,0) = 3; A(1,1) = 4; A(1,2) = 5;
    A(2,0) = 3; A(2,1) = 4; A(2,2) = 5;

    B(0,0) = 1; B(0,1) = 2; B(0,2) = 3;
    B(1,0) = 2; B(1,1) = 3; B(1,2) = 4;
    B(2,0) = 3; B(2,1) = 4; B(2,2) = 5;

    std::cout << "\nTrivial Multiplications: \n";
    std::cout << A;
    std::cout << B;
    std::cout << A*B;
    std::cout << B*A;
    std::cout << "\nPre-Transpose: \n";
    std::cout << A.T()*B;
    std::cout << B.T()*A;
    std::cout << "\nPost-Transpose: \n";
    std::cout << A*B.T();
    std::cout << B*A.T();
    std::cout << "\nTranspose: \n";
    std::cout << A.T()*B.T();
    std::cout << B.T()*A.T();


    // 4. Inner Product and Outer Product
    C_Matrix_Dense a(1,4), b(4, 1);

    a(0,0) = 1; a(0,1) = 2; a(0,2) = 3; a(0,3) = 4;
    b(0,0) = 1; b(1,0) = 2; b(2,0) = 3; b(3,0) = 4;

    std::cout << "\n Inner Product \n";
    std::cout << a*b;
    std::cout << "\n Outer Product \n";
    std::cout << b*a;

    /*
    // 5. Assign to Sparse Matrix
    C_Matrix_Sparse Z1(6,6);
    C_Matrix_Sparse Z2(6,6);

    std::vector<int> c = {1, 2, 3, 4};
    std::vector<int> d = {1, 2, 3, 4};

    std::cout << "\n Assignment to Sparse Matrix, 1: \n";
    Z1.add_matr(A1, c, d);
    std::cout << A1;
    std::cout << Z1;

    std::vector<int> e = {0, 1, 4, 5};
    std::vector<int> f = {0, 1, 4, 5};

    std::cout << "\n Assignment to Sparse Matrix, 2: \n";
    Z2.add_matr(A1, e, f);
    std::cout << A1;
    std::cout << Z2;
    */

    return 0;
}