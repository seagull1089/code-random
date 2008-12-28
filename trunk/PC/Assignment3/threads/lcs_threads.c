#include"../common.h"
#include<pthread.h>

int **weights;
char *string1,*string2;
int N, M;

pthread_mutex_t **locks;
struct thread_params{
	int id;
	int s_index;
	int l_index;
	int col_height;
};
int max_threads;
struct thread_params *thread_info;
void *initialize_locks(void *thread_info);
void *compute_weights(void *thread_info);
int get_value(int id,int col,int row);
void set_value(int id,int col,int row,int value);

int main(int argc,char **argv){

if(argc < 4) {
	printf("Usage : progname file1 file2 num_threads");
	exit(1);
}

max_threads = atoi(argv[3]);
if(max_threads == 0){
	printf(" 0 threads are not allowed. Must be some positive number of threads \n");
	exit(1);
}
string1 = read_data(argv[1]);
string2 = read_data(argv[2]);

N = strlen(string1)+1;
M = strlen(string2)+1;


int i,j;

/* Initialize the locks. No of locks per thread = M (column height) */
locks = (pthread_mutex_t **)malloc((max_threads-1)*sizeof(pthread_mutex_t *));
for(i=0;i<max_threads-1;i++){
	locks[i] = (pthread_mutex_t *)malloc(M*sizeof(pthread_mutex_t));
	for(j=0;j<M;j++){
		pthread_mutex_init(&locks[i][j],NULL);
	}
}

/* Create the thread messages having an id,start index and last index */
thread_info = (struct thread_params *)malloc(max_threads*(sizeof(struct thread_params)));
int chunk = (N-1)/max_threads;
for(i=0;i<max_threads;i++){
	/* divide the columns equally to all threads */
	thread_info[i].id = i;
	thread_info[i].s_index=(i)*chunk+1;
	thread_info[i].col_height=M;
	if(i==max_threads-1){
		thread_info[i].l_index= N;
	}else{
		thread_info[i].l_index=(i+1)*chunk;
	}

}

/* creat a matrix of MxN M rows N columns */
/* Initialize the weights */
weights = (int **) malloc((M)*sizeof(int *));
for(i=0;i<M;i++){
	weights[i]  = (int *) malloc((N)*sizeof(int));
	for(j=0;j<N;j++){
		weights[i][j] = 0;
	}
}

pthread_t *tids = (pthread_t *)malloc(max_threads*sizeof(pthread_t));
void *tret;

for(i=0;i<max_threads-1;i++){
	pthread_create(&tids[i],NULL,initialize_locks,&thread_info[i]);
}

for(i=0;i<max_threads-1;i++){
	pthread_join(tids[i],&tret);
}

struct timeval t1,t2;
gettimeofday(&t1,NULL);

for(i=0;i<max_threads;i++){
	pthread_create(&tids[i],NULL,compute_weights,&thread_info[i]);
}

for(i=0;i<max_threads;i++){
	pthread_join(tids[i],&tret);
}


//print_matrix(weights,M,N);

printf("No of matching chars = %d \n",weights[M-1][N-1]);

gettimeofday(&t2,NULL);
	int total_seconds = (t2.tv_sec - t1.tv_sec)*1000000 +
			    (t2.tv_usec - t1.tv_usec);

	printf("total time taken to calculate matrix = %f\n",total_seconds/1000000.0);

output_lcs(weights,string1,string2);
free(tids);
free(thread_info);
free(locks);
free(string1);
free(string2);
free(weights);
}

int get_value(int id,int row,int col){
	int value = -1;
	/* Identify if the asked value  belongs to a critical column or not */
	if((id !=0) && (row > 0) &&(thread_info[id-1].l_index == col)){
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
	if((id!=max_threads-1) && (thread_info[id].l_index == col) && (row>0)){
		/* this is a critical section and there is already a lock on it*/
		/* just write the value and unlock the lock */
		pthread_mutex_unlock(&locks[id][row]);	
		//printf("unlocked for id = %d row = %d col = %d\n",id,row,col); 
	}

}

void *initialize_locks(void *thread_info){
	struct thread_params  *t_info = (struct thread_params *)thread_info;	
	int i =0;
	for(i=1;i<t_info->col_height;i++){
		pthread_mutex_lock(&locks[t_info->id][i]);
	}	
}

void *compute_weights(void *thread_info){
	struct thread_params  *t_info = (struct thread_params *)thread_info;	
//	printf("%d %d %d %d\n",t_info->id,t_info->s_index,t_info->l_index,t_info->col_height);
	int i =0,j=0;
	for(i=1;i<t_info->col_height;i++){
		for(j=t_info->s_index;j<=t_info->l_index;j++){
			//printf("Indices id = %d %d %d \n",t_info->id, i,j);
			if(string2[i-1]==string1[j-1]){
				set_value(t_info->id,i,j,get_value(t_info->id,i-1,j-1)+1);
			}else{
				set_value(t_info->id,i,j,
					max(get_value(t_info->id,i,j-1),get_value(t_info->id,i-1,j)));
			}
		}

	}	
}

