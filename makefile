EigenDir = /home/staff/b/bensinghdhas/tools/eigen-3.4.0/
CPPFLAGS = -I$(EigenDir)

matrix_test_1:
	g++ -std=c++17 -g -c */main_MatrixDemo_1.cpp -I$(CPPFLAGS)
	g++ -std=c++17 -o main.exe main_MatrixDemo_1.o

clean_MT_1:
	rm -v main_MatrixDemo_1.o
	rm -v main.exe
