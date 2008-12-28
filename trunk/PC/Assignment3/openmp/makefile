lcs_openmp:lcs_openmp.c common.o
	gcc-4.2   -fopenmp lcs_openmp.c common.o -o lcs_openmp -g 
common.o:../common.h ../common.c
	gcc-4.2  -c ../common.c
clean:
	rm -rf *.o lcs_openmp
