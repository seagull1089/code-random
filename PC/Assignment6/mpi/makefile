mpi_method:mpi_method.c utils.o
	 mpicc -I ../ -o mpi_method mpi_method.c  utils.o -g 
utils.o:../utils.h ../utils.c
	gcc -c ../utils.c -g 
clean:
	rm -rf *.o mpi_method
	
