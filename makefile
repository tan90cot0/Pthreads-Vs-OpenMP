serial:
	g++ LU_serial.cpp -O3 -std=c++2a -o LU_serial.o

run_serial: serial
	./LU_serial.o 5 1

openmp:
	g++ LU_openmp.cpp -O3 -std=c++2a -fopenmp -o LU_openmp.o

run_openmp:
	./LU_openmp.o 5 2

pthreads:
	g++ LU_pthreads.cpp -O3 -std=c++2a -pthread -o  LU_pthreads.o

run_pthreads: pthreads
	./LU_pthreads.o 5 1

check_hari:
	g++ checker.cpp -o checker.o

run_hari:
	./checker.o