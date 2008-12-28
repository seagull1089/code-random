thread_method:thread_method.c	utils.o
	cc -I ../ thread_method.c utils.o -lpthread -o thread_method -g 
utils.o: ../utils.c ../utils.h 
	cc  -c ../utils.c -g 
clean:
	rm -rf *.o thread_method
	
