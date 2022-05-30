eigenPath= $(CURDIR)/
solverIntef= $(CURDIR)/Solver_Interfaces/
bar_localDirEigen:
	g++ -std=c++17 -g -c */main_1D_Template.cpp -I$(eigenPath) -I$(solverIntef)
	g++ -std=c++17 -o main.exe main_1D_Template.o
	./main.exe

bar_includePathEigen:
	g++ -std=c++17 -g -c */main_1D_Template.cpp
	g++ -std=c++17 -o main.exe main_1D_Template.o
	./main.exe

matrix_test_1:
	g++ -std=c++17 -g -c */main_MatrixDemo_1.cpp
	g++ -std=c++17 -o main.exe main_MatrixDemo_1.o

clean_MT_1:
	rm -v main_MatrixDemo_1.o
	rm -v main.exe
