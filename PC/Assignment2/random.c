#include<stdlib.h>
#include<stdio.h>

int get_random(int lower,int higher,int N,int entropy){
	int r_num = -1;
//	printf("DEBUG lower = %d higher = %d \n",lower,higher);
	unsigned int x = (entropy)*time(NULL);
	if(lower > higher || N==0) {
		return r_num;
	}
	rand_r(&x);
	while(r_num < lower ||  r_num >higher){
	//printf("DEBUG2 random = %d\n",r_num);	
		r_num = rand_r(&x)%N;
	}
//	printf("DEBUG random = %d\n",r_num);	
	return r_num;
}
