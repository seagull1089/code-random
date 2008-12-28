#include <mpi.h>  // remember to include the mpi header
#include <stdio.h>
#include<stdlib.h>
#include<sys/time.h>

extern int get_random(int lower,int higher,int N,int entropy);
void set_intitial_partitions(int indices[],int N,int num_procs,int chunksize);
int partition(float A[], int p, int r,float pivot);
int get_local_ranks_sum(int a[],int num_procs);
void set_partitions(int indices[],int local_ranks[],int num_procs,int left);
int  get_new_pivot(int indices[],int num_procs,int chunksize,int N);
int has_converged(int local_ranks_sum,int final_rank);
int main(int argc, char *argv[])
{

  int num_procs, myid, name_len,N,partition_size;
  char proc_name[MPI_MAX_PROCESSOR_NAME];
  char *fileName;
  // Initialize MPI
	float *array,var;
	int partitition_size;
  if(argc<2){
	printf("Usage : progname fileName\n");
	exit(1);
  }
  fileName = argv[1];
  MPI_Status status; 
  int tag = 100;
  MPI_Init(&argc, &argv);

  // Obtain the number of processes
  MPI_Comm_size(MPI_COMM_WORLD, &num_procs);

  // Obtain the process id
  MPI_Comm_rank(MPI_COMM_WORLD, &myid);

  // Obtain name of machine process is executing on
  MPI_Get_processor_name(proc_name, &name_len);
	
//  printf("Hello World from processor %d out of %d, executing on %s filename = %s \n",
//  	      myid, num_procs, proc_name,fileName);
int chunksize;
int i =0,final_rank;

    if(myid == 0){
	FILE *fp;
	fp = fopen(fileName,"r");
	fscanf(fp,"%d",&N);
 	 final_rank = N/2;

	chunksize = (N/num_procs);
	for(i=1;i<num_procs;i++){
		MPI_Send(&chunksize,1,MPI_INT,i,tag,MPI_COMM_WORLD);
		if(i == num_procs-1){
			partition_size = N - chunksize*(num_procs-1);
		}else {
			partition_size = chunksize;
		}
		MPI_Send(&partition_size,1,MPI_INT,i,tag,MPI_COMM_WORLD);
	}

	/* for processor 0 */ 
	partition_size = chunksize;
	//printf("read in processid  %d  partition_size = %d chunksize =%d\n" 
	//	, myid , partition_size,chunksize);

	int scan_count=0;
	/* read the elements for processor 0 */	
	array = (float *)malloc(partition_size*sizeof(float));
	while(scan_count < partition_size){
			fscanf(fp,"%f",&var);
		  	array[scan_count] = var;
			scan_count++;	
	}

	for(i=1;i<num_procs;i++){
		float *tmparray;
		int elements_to_read = chunksize;
		if(i==num_procs-1){
			elements_to_read =  N - chunksize*(num_procs-1); 
		}
		tmparray = (float *)malloc(elements_to_read*sizeof(float));
		scan_count=0;
		while(scan_count < elements_to_read){
			fscanf(fp,"%f",&var);
			tmparray[scan_count] = var;
			scan_count++;
		}
		MPI_Send(tmparray,scan_count,MPI_FLOAT,i,tag,MPI_COMM_WORLD);
		free(tmparray);	
	}
		
	fclose(fp);
	//  end of loading for processor 0.
	}else{
		/* rest of the processors just recieve 
		the input from the processor 0 */

		MPI_Recv(&chunksize,1,MPI_INT,0,tag,MPI_COMM_WORLD,&status);
		MPI_Recv(&partition_size,1,MPI_INT,0,tag,MPI_COMM_WORLD,&status);

	 	array = (float *)malloc(partition_size * sizeof(float));	
	//	printf("read in processid  %d  partition_size = %d chunksize =%d\n" 
	//	, myid , partition_size,chunksize);
		MPI_Recv(array,partition_size,MPI_FLOAT,0,tag,MPI_COMM_WORLD,&status);
		
	} 

  
     /* 
	for(i=0;i<num_procs;i++){
		if(myid == i){
			int j = 0;
			printf(" processor id  = %d , name = %s partition_size=%d\n",
					myid,proc_name,partition_size);
			for(j=0;j<partition_size;j++){
				printf("%f  ",array[j]);
			}
			printf("\n");
		}
	}
	*/	
/* end of initialising arrays */

/* main logic begins here. process id = 0 is the root process */
	if(myid ==0){
		double t1,t2;
		t1= MPI_Wtime();
	
		int counter = 0;
		int pivot_index = get_random(0,N,N,(++counter));
		int i = 0,debug=0;
		int *partition_indices = (int *)malloc(num_procs*2*sizeof(int));
		int *local_ranks = (int *) malloc(num_procs*sizeof(int));
		set_intitial_partitions(partition_indices,N,num_procs,chunksize);

		int local_ranks_sum = -1;
		int run_flag = 1;
		float median = -1;
		
		while(run_flag){
		for(i=1;i<num_procs;i++){
			MPI_Send(&run_flag,1,MPI_INT,i,tag,MPI_COMM_WORLD);
			MPI_Send(partition_indices+2*i,2,MPI_INT,i,tag,MPI_COMM_WORLD);
		}
		
		for(i=1;i<num_procs;i++){
			MPI_Send(&pivot_index,1,MPI_INT,i,tag,MPI_COMM_WORLD);
		}
		float pivot_value,tmp_pivot;
		for(i=0;i<num_procs;i++){
			if(i==0){
				tmp_pivot = -1;
				if(pivot_index < partition_size){
					tmp_pivot = array[pivot_index];	
				}

			}else{
				MPI_Recv(&tmp_pivot,1,MPI_FLOAT,i,tag,MPI_COMM_WORLD,&status);
			}
			if(tmp_pivot >=0){
				pivot_value = tmp_pivot;
				if(debug)
				printf("Got the pivot value as -- %f index = %d pid = %d\n",pivot_value,pivot_index,i);
				tmp_pivot =-1;
			}
		}

		/* send back the pivot value to all 
		processes so that they start partition */
		for(i=1;i<num_procs;i++){
			MPI_Send(&pivot_value,1,MPI_FLOAT,i,tag,MPI_COMM_WORLD);
		}
		int l_rank=-1;
		for(i=0;i<num_procs;i++){
			if(i==0){
				l_rank = partition(array,partition_indices[0],partition_indices[1],pivot_value) + 1;
			}else{
				MPI_Recv(&l_rank,1,MPI_INT,i,tag,MPI_COMM_WORLD,&status);
			}
			local_ranks[i]=l_rank;	
		//	printf(" id = %d , local rank = %d \n",i,l_rank);
		}		
			
		local_ranks_sum = get_local_ranks_sum(local_ranks,num_procs);

		if(debug)
		printf("Got local sum ranks = %d \n" , local_ranks_sum);

		if((local_ranks_sum == final_rank) || has_converged(local_ranks_sum,final_rank)){
			run_flag = 0;
			median = pivot_value;	
		}else if(local_ranks_sum < final_rank){
				set_partitions(partition_indices,local_ranks,num_procs,0);
		}else if(local_ranks_sum > final_rank){
				set_partitions(partition_indices,local_ranks,num_procs,1);
		}
			pivot_index = get_new_pivot(partition_indices,num_procs,chunksize,N);
		
			
		} // end of while loop

		/* say bye bye */	
		for(i=1;i<num_procs;i++){
			MPI_Send(&run_flag,1,MPI_INT,i,tag,MPI_COMM_WORLD);
		}
		
		t2 = MPI_Wtime();	
		printf("total Execution time  = %f\n" ,(t2-t1)/1000);
		
		printf("Median value = %f \n" , median);
	
	}else{

		/* to be executed by other processors other than myid =0 */
		int recv_pivot_index;
		int indices[2];
		int run_flag;
		int debug =0;	
		MPI_Recv(&run_flag,1,MPI_INT,0,tag,MPI_COMM_WORLD,&status);
		while(run_flag==1){
	
		MPI_Recv(indices,2,MPI_INT,0,tag,MPI_COMM_WORLD,&status);
		if(debug){	
		printf("recieved run_flag = %d indices = %d %d id = %d \n",
			run_flag,indices[0],indices[1],myid);
		}
		MPI_Recv(&recv_pivot_index,1,MPI_INT,0,tag,MPI_COMM_WORLD,&status);

		float send_pivot_value = -1;
		if((recv_pivot_index >= myid*chunksize) && 
		  (recv_pivot_index < (myid*chunksize) + partition_size)){
			send_pivot_value = array[recv_pivot_index - (myid*chunksize)];
		}
		if(debug){
		printf("chunksize = %d , partition_size = %d , send_pivot_value = %f , recv_pivot_index =%d myid = %d\n",
			chunksize,partition_size,send_pivot_value,recv_pivot_index,myid);
		}
		MPI_Send(&send_pivot_value,1,MPI_FLOAT,0,tag,MPI_COMM_WORLD);
		
		/* recieve the finalised pivot value from the processor 0 */
		float pivot_value;
		int local_rank;
		MPI_Recv(&pivot_value,1,MPI_FLOAT,0,tag,MPI_COMM_WORLD,&status);
	//	printf(" myid = %d recieved pivot value = %f \n",myid,pivot_value);
		local_rank = partition(array,indices[0],indices[1],pivot_value) + 1;
		MPI_Send(&local_rank,1,MPI_INT,0,tag,MPI_COMM_WORLD);
	//	printf("processor = %d , local_rank = %d \n",myid,local_rank);
		MPI_Recv(&run_flag,1,MPI_INT,0,tag,MPI_COMM_WORLD,&status);
		}
	} 
	
 // Last call to MPI (REQUIRED)
  MPI_Finalize();
}
int  get_new_pivot(int indices[],int num_procs,int chunksize,int N){
	static int counter = 0;
	int pivot_index;
	int pid = get_random(0,num_procs-1,num_procs,(++counter));	
	
	pivot_index = get_random((chunksize*pid) + indices[2*pid] , (chunksize*pid)+indices[2*pid+1],N,(++counter));
//	printf(" pid = %d , pivot_index = %d \n",pid,pivot_index);
	return pivot_index;	
}
void set_partitions(int indices[],int local_ranks[],int num_procs,int left){
	int i = 0;
	for(i=0;i<num_procs;i++){
		if(left==1){
			if(local_ranks[i]-1 > indices[2*i])
			indices[2*i+1] = local_ranks[i]-1;
		}else{
			if(local_ranks[i] < indices[2*i+1])
			indices[2*i] = local_ranks[i];
		}
	}
}
int get_local_ranks_sum(int a[],int num_procs){
	int i =0;
	int sum=0;
	for(i=0;i<num_procs;i++){
		sum = sum+a[i];
	}	
	return sum+1;
}
void set_intitial_partitions(int indices[],int N,int num_procs,int chunksize){
	int i =0;
	for(i=0;i<num_procs-1;i++){
		indices[2*i] = 0;
		indices[2*i+1] = chunksize-1;
	}
	indices[2*(num_procs-1)] = 0;
	indices[2*(num_procs-1)+1]=N-((num_procs-1)*chunksize) -1 ;
/*	for(i=0;i<num_procs;i++){
		printf(" i = %d indices = %d , %d \n",i,indices[2*i],indices[2*i+1]);
	}
*/
	return;
}

int has_converged(int rank,int required_rank){
	static int ranks[10]={0};
	static int maxrank = 0,pointer=-1,iter=0;;
	if(rank <=required_rank && rank >=maxrank){
		maxrank = rank;
		ranks[(++pointer)%10] = maxrank;
		iter++;
	}
	
	int i =0, converged = 0;
	if(iter>10){
		converged =1;
	for(i=0;i<10;i++){
		if(ranks[i] != maxrank){
			converged = 0;
			break;
		}
	}
	}
	return converged;
}


int partition(float A[], int p, int r,float pivot)
{
	/* should you changed the intial value of i?*/
	int             j,i = p-1;
	float           tmp;
	 /* loop runs from p to r */
	for (j = p; j <= r; j++) {
		/* modified the condition of < to <=*/ 
		if (A[j] < pivot) {
			i = i + 1;
			tmp = A[i];
			A[i] = A[j];
			A[j] = tmp;
		}
	}
	/* just return the last position */
	return i ;
}
