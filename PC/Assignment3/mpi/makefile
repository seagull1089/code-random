lcs_mpi:lcs_mpi.c common.o
	mpicc -o lcs_mpi common.o lcs_mpi.c
common.o:../common.c ../common.h
	gcc -c ../common.c 
clean:
	rm -rf *.o lcs_mpi	
