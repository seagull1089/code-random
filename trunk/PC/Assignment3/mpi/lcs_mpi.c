#include"../common.h"
#include<mpi.h>

int main(int argc,char **argv){

  int myid,name_len,num_procs;
  char proc_name[MPI_MAX_PROCESSOR_NAME];
  

  int tag=100;   

  MPI_Init(&argc, &argv);

  // Obtain the number of processes
  MPI_Comm_size(MPI_COMM_WORLD, &num_procs);

  // Obtain the process id
  MPI_Comm_rank(MPI_COMM_WORLD, &myid);

  // Obtain name of machine process is executing on
  MPI_Get_processor_name(proc_name, &name_len);

  MPI_Status status; 
  int chunk;
  int N,M;
 //pass the strings to each of the processorsa

 char *string1=NULL,*string2=NULL,*dummy1=NULL,*dummy2=NULL,*string1_original;
  if(myid == 0){
	
  	dummy1 = read_data(argv[1]); //will become N
 	dummy2 = read_data(argv[2]); // will form M

	string1 = (char *)malloc(strlen(dummy1)*sizeof(char));
	string1_original = (char *)malloc(strlen(dummy1)*sizeof(char));
	string2 = (char *)malloc(strlen(dummy2)*sizeof(char));

	if(dummy1[0]==3 && dummy2[0]==3){
		strcpy(string1,dummy1+1);
		strcpy(string2,dummy2+1);
	}else{
		strcpy(string1,dummy1);
		strcpy(string2,dummy2);
	}
	
	strcpy(string1_original,string1);

	free(dummy1);
	free(dummy2);
	 N = strlen(string1);
       	 M = strlen(string2);
	//printf(" strings are ->%s<- len = %d ->%s<- len = %d \n",string1,N,string2,M);

	chunk = N/num_procs;
	int partition_size = chunk;
	int i =0;
/*	for(i=0;i<N;i++){
		printf(" i = %d , string1[i] = %d\n",i,string1[i]);
	}
*/
	for(i=1;i<num_procs;i++){
		if(i==num_procs-1){
			 partition_size= N-chunk*i;
		}
		MPI_Send(&partition_size,1,MPI_INT,i,tag,MPI_COMM_WORLD);
		MPI_Send(&M,1,MPI_INT,i,tag,MPI_COMM_WORLD);
	}

	partition_size = chunk;
	for(i=1;i<num_procs;i++){
		if(i==num_procs-1){
			partition_size = N-chunk*i;
		}
	//	printf(" Sending to id = %d , start = %d ,partition_size=%d\n",i,i*chunk,partition_size);
		MPI_Send(string1+i*chunk,partition_size,MPI_CHAR,i,tag,MPI_COMM_WORLD);
		MPI_Send(string2,M,MPI_CHAR,i,tag,MPI_COMM_WORLD);
	}
	chunk = N/num_procs;
	
	string1[chunk]='\0';
	N = strlen(string1);
	
  }else{
			
	MPI_Recv(&N,1,MPI_INT,0,tag,MPI_COMM_WORLD,&status);
	string1 = (char *)malloc((N+1)*sizeof(char));
	MPI_Recv(&M,1,MPI_INT,0,tag,MPI_COMM_WORLD,&status);
	string2 = (char *)malloc((M+1)*sizeof(char));
//	printf(" id = %d , N = %d , M = %d \n", myid,N,M);
	MPI_Recv(string1,N,MPI_CHAR,0,tag,MPI_COMM_WORLD,&status);
	string1[N]='\0';
	int i =0;

	MPI_Recv(string2,M,MPI_CHAR,0,tag,MPI_COMM_WORLD,&status);
	string2[M]='\0';
	//printf("id = %d string1 %s len = %d string2 %s \n",myid,string1,strlen(string1),string2); 	
  }

 // Now compute the weights matrix.
 // Everyone deoes it on their own.
	int **weights;
	weights = (int **)malloc( (M+1)*sizeof(int *));
//	printf("Adress of weights = %d id = %d\n",weights,myid);
	int i =0,j=0;
	for(i=0;i<M+1;i++){
		weights[i] = (int *)malloc((N+1)*sizeof(int));
		for(j=0;j<N+1;j++){
			weights[i][j]=0;
		}
	}
	
	//printf(" myid = %d , N = %d M = %d , string1 = %s string2 = %s \n",myid,N,M,string1,string2);

	double t1;
	if(myid ==0){
		t1=MPI_Wtime();
	}	
	for(i=1;i<M+1;i++){
		for(j=1;j<N+1;j++){

			if(myid!=0 && j==1){
				MPI_Recv(&weights[i][j-1],1,MPI_INT,myid-1,tag,MPI_COMM_WORLD,&status);
			}	
			

			if(string1[j-1]==string2[i-1]){
				weights[i][j]=weights[i-1][j-1]+1;
			}else{
				weights[i][j]= max(weights[i][j-1]
						,weights[i-1][j]);
							
			}

			if((myid!= num_procs-1) && j==N){
					MPI_Send(&weights[i][j],1,MPI_INT,myid+1,tag,MPI_COMM_WORLD);
				}	
		}
	}


	if(myid ==0){
		int final_weight = 0;
		if(num_procs-1>0){
			MPI_Recv(&final_weight,1,MPI_INT,num_procs-1,tag,MPI_COMM_WORLD,&status);
		}else{
			final_weight = weights[M][N];
		}
		printf("No of matching chars  = %d \n",final_weight);
		double t2 = MPI_Wtime();
		printf("Total time taken to compute the matrix = %f\n",(t2-t1));
	
	}else{
		if(myid == num_procs -1){
			MPI_Send(&weights[M][N],1,MPI_INT,0,tag,MPI_COMM_WORLD);
		}
	}


	i = strlen(string2);
	j = strlen(string1);
	int last_index=0;
	char *aggregate_lcs;
	char *previous_proc_lcs;
	int prev_proc_lcs_length=0;
	int aggregate_lcs_length=0;
	if(myid != num_procs-1){
		MPI_Recv(&i,1,MPI_INT,myid+1,tag,MPI_COMM_WORLD,&status);
	//	printf(" i = %d j = %d myid = %d \n",i,j,myid);
		MPI_Recv(&prev_proc_lcs_length,1,MPI_INT,myid+1,tag,MPI_COMM_WORLD,&status);
		aggregate_lcs_length = j+prev_proc_lcs_length+1;
		previous_proc_lcs= (char *)malloc(prev_proc_lcs_length*sizeof(char));
		MPI_Recv(previous_proc_lcs,prev_proc_lcs_length,MPI_CHAR,myid+1,tag,MPI_COMM_WORLD,&status);
					
	}else{
		previous_proc_lcs=malloc(sizeof(char));
		previous_proc_lcs[0]='\0';
		aggregate_lcs_length=strlen(string1)+1;
	}
	aggregate_lcs = (char *)malloc(aggregate_lcs_length*sizeof(char));
	int l =0;
	for(l=0;l<strlen(string1);l++){
		aggregate_lcs[l]='_';
	}
	aggregate_lcs[l]='\0';	
	int k =strlen(string1)-1;
	while(i>0 && j>0){
		if(string2[i-1]==string1[j-1]){
			aggregate_lcs[k]=string2[i-1];	
			k--;i--;j--;
		}else{
			if(weights[i-1][j]>weights[i][j-1]){
				i--;
			}else{
				j--;
				k--;		
			}		
		}
	}		
//	aggregate_lcs[k]='\0';	
	last_index = i;
	strcat(aggregate_lcs,previous_proc_lcs);
	if(myid != 0 && (myid-1)>=0){
		MPI_Send(&last_index,1,MPI_INT,myid-1,tag,MPI_COMM_WORLD);
		MPI_Send(&aggregate_lcs_length,1,MPI_INT,myid-1,tag,MPI_COMM_WORLD);
		MPI_Send(aggregate_lcs,aggregate_lcs_length,MPI_CHAR,myid-1,tag,MPI_COMM_WORLD);
	}else{
		printf("The LCS found is : %s\n",aggregate_lcs);		
	}

	free(aggregate_lcs);
	free(previous_proc_lcs);
	free(weights);
	free(string1);
	free(string2);	


  //printf("Hello World from processor %d name = %s\n",myid,proc_name);
  // Last call to MPI (REQUIRED)
  MPI_Finalize();

}
