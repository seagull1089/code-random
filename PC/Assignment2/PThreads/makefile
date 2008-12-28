thread_method:thread_method.c random.o
	gcc -g thread_method.c random.o -lpthread -o thread_method 
random.o:../random.c
	gcc -c ../random.c
clean:
	rm *.o thread_method 
