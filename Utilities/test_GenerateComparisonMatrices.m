% This code replicates the C++ code used to test dense matrix
% implementation, displaying the same results to the terminal for the
% purposes of comparison.
%
% DIJ (5-9-22)



%% 2. Multiplications
A1 = zeros(4,4); B1 = zeros(4,4); C1 = zeros(4,6);
A2 = zeros(4,3); B2 = zeros(3,4); C2 = zeros(3,3);

% i. Part I
A1(1,1) = 4; A1(1,3) = 3; A1(2,2) = 2.5; A1(2,3) = 4.5;
A1(3,3) = 3; A1(3,4) = 1; A1(4,1) = 1.7; A1(4,3) = 0.9;

B1(1,2) = 4; B1(1,3) = 2; B1(2,2) = 4.5; B1(2,3) = 1.5;
B1(3,3) = 1; B1(2,4) = 7; B1(4,1) = 1.2; B1(4,3) = 0.6;

C1(1,1) = 3; C1(1,4) = 4; C1(2,2) = 4.5; C1(2,3) = 1.5;
C1(3,4) = 1; C1(2,4) = 2; C1(4,2) = 1.2; C1(4,3) = 0.6;
C1(4,6) = 2; C1(4,5) = 6.7;

% ii. Part II
A2(1,1) = 4; A2(1,3) = 3; A2(2,2) = 2.5; A2(2,3) = 4.5;
A2(3,3) = 3; A2(3,2) = 1; A2(4,1) = 1.7; A2(4,3) = 0.9;

B2(1,2) = 4; B2(1,3) = 2; B2(2,2) = 4.5; B2(2,3) = 1.5;
B2(3,1) = 1.2; B2(3,3) = 0.6; B2(3,4) = 9;

C2(1,1) = 3; C2(1,3) = 4; C2(2,2) = 4.5; C2(2,3) = 1.5;
C2(3,2) = 1.2; C2(3,3) = 0.6;

disp('Before Multiplication');
disp(A1);
disp(B1);
disp(C1);
disp(A2);
disp(B2);
disp(C2);

disp('After Multiplication');
disp(A1*B1);
disp(A1*C1);
disp(A1*B2');
disp(B1*B2');
disp(C2'*A2');
disp(C2*A2');
disp(B1*(B2' + A2));
disp(B1*(B2' - A2*5));
disp(C2*6*A2');
disp(C1*C1');
disp(C1'*C1);


%% 3. Trivial Multiplication
A = [3, 4, 5; 3, 4, 5; 3, 4, 5];
B = [1, 2, 3; 2, 3, 4; 3, 4, 5];

disp('Trivial Multiplication');
disp(A);
disp(B);
disp(A*B);
disp(B*A);

disp('Pre-Transpose');
disp(A'*B);
disp(B'*A);

disp('Post-Transpose');
disp(A*B');
disp(B*A');

disp('Transpose');
disp(A'*B');
disp(B'*A');


%% 4. Inner Product and Outer Product
a = [1, 2, 3, 4];
b = a';

disp('Inner Product');
disp(a*b);
disp('Outer Product');
disp(b*a);


%% 7. Sachin Tests
K1 = C2 + 5*(B'*B);
K2 = C2 + 5*B'*B;   %#ok<MHERM>

K3 = B1 + (B2'*B2);
K4 = B1 + B2'*B2;

K5 = B1 + 5*(B2'*B2);
K6 = B1 + 5*B2'*B2; %#ok<MHERM>

disp("Square Matrices")
disp(K1)
disp(K2)

disp("Non-Square Matrices, No Constant")
disp(K3)
disp(K4)

disp("Non-Square Matrices, Multiplied by Constant")
disp(K5)
disp(K6)


