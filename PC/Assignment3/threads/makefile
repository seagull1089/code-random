lcs_thread:lcs_threads.c common.o
	gcc  lcs_threads.c common.o -o lcs_thread -g -lpthread
common.o:../common.h ../common.c
	gcc  -c ../common.c -g 
clean:
	rm -rf *.o lcs_thread
