INC_DIR_1 = ../
CPPFLAGS = -I $(INC_DIR_1)
VPATH = ./Demo_Scripts/main_MatrixDemo_1.cpp

sparse_test:
	g++ -std=c++17 -g -c */main_MatrixDemo_1.cpp $(CPPFLAGS)
	g++ -std=c++17 -o main.exe main_MatrixDemo_1.o

clean_sparse:
	rm -v main_MatrixDemo_1.o
	rm -v main.exe
