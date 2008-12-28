mpi_method:mpi_method.c utils.o 
	 mpicc -o mpi_method mpi_method.c  utils.o   -g 
utils.o:utils.c utils.h
	gcc -c utils.c	-g 
clean:
	rm -rf *.o mpi_method
	
