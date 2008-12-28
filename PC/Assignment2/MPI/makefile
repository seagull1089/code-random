mpi_method:mpi_method_with_array_distribution.c random.o
	mpicc -o mpi_method random.o mpi_method_with_array_distribution.c
random.o:../random.c
	gcc -c ../random.c
clean:
	rm -rf *.o mpi_method	
