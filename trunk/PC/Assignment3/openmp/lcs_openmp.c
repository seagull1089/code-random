#include <omp.h> 
#include <stdio.h>
#include"../common.h"
#include<pthread.h>

int **weights;
int N,M;
pthread_mutex_t **locks;
void set_value(int id,int row,int col,int value);
int get_value(int id,int row,int col);
int average_chunk;
int max_threads;

int main(int argc, char **argv)
{

        int   myid;

	if(argc<3) {
		printf("Usage:progname filename1 filename2 num_threads\n"); 
		exit(1);
	}
        max_threads = atoi(argv[3]);
        char *string1 = read_data(argv[1]);
        char *string2 = read_data(argv[2]);

        omp_set_num_threads(max_threads);
      //  int N,M;
        N= strlen(string1)+1;
        M= strlen(string2)+1;

        //int **weights;
        weights = (int **) malloc(M*sizeof(int *));
        int i=0,j=0;
        for(i=0;i<M;i++){
		weights[i] = (int *)malloc(N*sizeof(int));
		for(j=0;j<N;j++){
			weights[i][j]= 0;
		}
        }



	locks = (pthread_mutex_t **)malloc((max_threads-1)*sizeof(pthread_mutex_t *));
	
	for(i=0;i<max_threads-1;i++){
		locks[i] = (pthread_mutex_t *)malloc(M*sizeof(pthread_mutex_t));
		for(j=0;j<M;j++){
			//printf(" i = %d , j = %d \n",i,j);		
			pthread_mutex_init(&locks[i][j],NULL);
		}
	}

//#pragma omp parallel
//{
//int i =0,j=0;
//#pragma omp parallel for shared(locks,max_threads)
	for(i=0;i<max_threads-1;i++){
		for(j=0;j<M;j++){
		//printf(" i = %d , j = %d \n",i,j);
			pthread_mutex_lock(&locks[i][j]);
		}
	}

//}
	average_chunk = (N/max_threads)+1 ;
	//printf("Average_chunk = %d N=%d....\n",average_chunk,N);
struct timeval t1,t2;
gettimeofday(&t1,NULL);
#pragma omp parallel  shared(weights,string1,string2,average_chunk) 
{
	int i=0,j=0;
	for(i=1;i<M;i++){
#pragma omp parallel for schedule(dynamic,average_chunk) shared(weights,string1,string2,average_chunk) num_threads(max_threads)
      		for(j=1;j<N;j++){
//#pragma omp parallel for schedule(dynamic,4) shared(weights,string1,string2) 
			if(string2[i-1]==string1[j-1]){
				int myid = omp_get_thread_num();
				set_value(myid,i,j,get_value(myid,i-1,j-1)+1);
			}else{
				set_value(myid,i,j,max(get_value(myid,i-1,j),get_value(myid,i,j-1)));
			}
		}
	}
}

	printf("The final weight = %d\n",weights[M-1][N-1]);
gettimeofday(&t2,NULL);
	int total_seconds = (t2.tv_sec - t1.tv_sec)*1000000 +
			    (t2.tv_usec - t1.tv_usec);

	printf("total time taken to calculate matrix = %f\n",total_seconds/1000000.0);
	//print_matrix(weights,M,N);
	output_lcs(weights,string1,string2);
	free(weights);
	free(locks);
	free(string1);
	free(string2);
}


int get_value(int id,int row,int col){
	int value = -1;
	/* Identify if the asked value  belongs to a critical column or not */
	if((id !=0) && (row > 0) &&((id)*average_chunk == col)){
		/* belongs to a critical section */
		//printf("waiting for (id-1) = %d row = %d col = %d \n",id-1,row,col);
		pthread_mutex_lock(&locks[id-1][row]);
		value = weights[row][col];
		pthread_mutex_unlock(&locks[id-1][row]);
		//printf("finished for (id-1) = %d row = %d col = %d \n",id-1,row,col);
	}else{
		value = weights[row][col];
	}
	return value;
}

void set_value(int id,int row,int col,int value){
	weights[row][col] = value;
	if((id!=max_threads-1) && (id*average_chunk == col) && (row>0)){
		/* this is a critical section and there is already a lock on it*/
		/* just write the value and unlock the lock */
		pthread_mutex_unlock(&locks[id][row]);	
		//printf("unlocked for id = %d row = %d col = %d\n",id,row,col); 
	}

}
