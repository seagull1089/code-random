#include"utils.h"
#include<mpi.h>

int compare(const void  *a,const void  *b){
	return (*(int *)a - *(int *)b);
}
void sample_sort(int A[],int length,int (*cmp)(const void *, const void *));
int main(int argc,char **argv){

  int num_procs, myid, name_len;
  
  char proc_name[MPI_MAX_PROCESSOR_NAME];

  // Initialize MPI
  MPI_Init(&argc, &argv);

  // Obtain the number of processes
  MPI_Comm_size(MPI_COMM_WORLD, &num_procs);

  // Obtain the process id
  MPI_Comm_rank(MPI_COMM_WORLD, &myid);

  // Obtain name of machine process is executing on
  MPI_Get_processor_name(proc_name, &name_len);

  MPI_Status status;
  
  int tag=100;

  char *filename = argv[1];
  
  int i,N,chunksize,last_proc_chunksize;
  int *mychunk;

  if(myid == 0 ) {
	FILE *fp = fopen(filename,"r");
	fscanf(fp,"%d ",&N);
//	printf("Total number of elements = %d\n",N);
	chunksize = N/num_procs;
	last_proc_chunksize = N - ((num_procs-1)*chunksize);
//	printf(" chunksize  = %d , last_proc_chunksize = %d \n",chunksize,last_proc_chunksize);
 	mychunk = readfromFile(fp,chunksize);


	int proc_chunk_size = chunksize;

	for(i=1;i<num_procs;i++){
		if(i == num_procs -1) {
			proc_chunk_size = last_proc_chunksize;
		}
	
		MPI_Send(&proc_chunk_size,1,MPI_INT,i,tag,MPI_COMM_WORLD);
	
		int *proc_chunk = readfromFile(fp,proc_chunk_size);
		MPI_Send(proc_chunk,proc_chunk_size,MPI_INT,i,tag,MPI_COMM_WORLD);
		free(proc_chunk);
	
	}
	fclose(fp);
 
  }else{
		
		MPI_Recv(&chunksize,1,MPI_INT,0,tag,MPI_COMM_WORLD,&status);
		mychunk = malloc(chunksize*sizeof(int));
			
		MPI_Recv(mychunk,chunksize,MPI_INT,0,tag,MPI_COMM_WORLD,&status);
		

   }

  MPI_Bcast(&N,1,MPI_INT,0,MPI_COMM_WORLD);
//  printf("received values : N = %d chunksize = %d myid = %d\n",N,chunksize,myid);

  MPI_Barrier(MPI_COMM_WORLD);
  double T1,T2;
  if(myid == 0) {
	T1 = MPI_Wtime();
  }  

 	
  sample_sort(mychunk,chunksize,compare);
  MPI_Barrier(MPI_COMM_WORLD);

  if(myid == 0){
	T2 = MPI_Wtime();
	printf(" Time taken for sample sort  for %d processors = %lf \n",num_procs,T2-T1); 
	FILE *fout = fopen("out.txt","w");
	for(i=0;i<chunksize;i++){
		fprintf(fout,"%d\n",mychunk[i]);

	}
	int recv_size = chunksize;
	for(i=1;i<num_procs;i++){
		if(i==num_procs-1) {
			recv_size = last_proc_chunksize;
		}
		int *recvbuf = malloc(recv_size*sizeof(int));
		MPI_Recv(recvbuf,recv_size,MPI_INT,i,tag,MPI_COMM_WORLD,&status);
		int k =0;
		for(k=0;k<recv_size;k++){
			fprintf(fout,"%d\n",recvbuf[k]);
		}
		free(recvbuf);
	}
	fclose(fout);
  }else {
	MPI_Send(mychunk,chunksize,MPI_INT,0,tag,MPI_COMM_WORLD);
  }	

  free(mychunk);
  MPI_Finalize();
}


void sample_sort(int A[],int length,int (*cmp)(const void *, const void *)){
	qsort(A,length,sizeof(int),cmp);
  		
	int num_procs,myid;

	MPI_Comm_size(MPI_COMM_WORLD, &num_procs);
  	MPI_Comm_rank(MPI_COMM_WORLD, &myid);
  	MPI_Status status;
	int tag = 100;
	int i = 0;
	// Send the p-1 elements now. 
	int segsize = length/num_procs;
	int *sample_elems = malloc((num_procs-1)*sizeof(int));
	for(i=0;i<num_procs-1;i++){
		sample_elems[i]= A[segsize*(i+1)];
	}
/*	
	for(i=0;i<length;i++){
		printf(" %d " , A[i]);
	}
	printf(" my id = %d \n",myid);
*/	
	int *all_procs_sample ;
	if(myid == 0){
		all_procs_sample = malloc(num_procs*(num_procs-1)*sizeof(int));	
	}
	MPI_Gather(sample_elems,num_procs-1, MPI_INT,
            all_procs_sample, num_procs-1, MPI_INT, 0,
            MPI_COMM_WORLD);

	free(sample_elems);

	int *bucket_limits = malloc((num_procs-1)*sizeof(int));
	if(myid==0){
		qsort(all_procs_sample,num_procs*(num_procs-1),sizeof(int),cmp);
		for(i=0;i<num_procs-1;i++){
			bucket_limits[i]=all_procs_sample[(num_procs-1)*(i+1)];
		}
/*
		for(i=0;i<num_procs-1;i++){
			printf(" bucket limit = %d ",bucket_limits[i]);
		}
		printf("\n");
*/
		free(all_procs_sample);
	}
	//BroadCast the bucket limits 
	MPI_Bcast(bucket_limits,num_procs-1, MPI_INT, 0, MPI_COMM_WORLD);
	int *partition_indices = get_partition_indices(A,length,
					bucket_limits,num_procs-1);

	free(bucket_limits);
/*
	for(i=0;i<num_procs-1;i++){
		printf(" b=%d -- %d ",bucket_limits[i],partition_indices[i]);
	}
	printf("; id =  %d\n",myid);
*/
	int *num_data_tobe_sent = malloc(num_procs*sizeof(int));
	for(i=0;i<num_procs;i++){
		int n_elems;
		if(i==0){
			n_elems = partition_indices[0] -0;
			n_elems = (n_elems >0) ? n_elems: 0;
		}else if(i==num_procs-1){

			if(partition_indices[num_procs-2] >=0){
				n_elems = length - partition_indices[num_procs-2];
			}else{
				n_elems = 0;
			}
		}else if(partition_indices[i-1]>=0 && partition_indices[i]>=0){
			n_elems=partition_indices[i]-partition_indices[i-1];
		}else{
			n_elems =0;
		}
		num_data_tobe_sent[i]=n_elems;
	}

	int *num_data_tobe_recv = malloc(num_procs*sizeof(int));

	MPI_Alltoall(num_data_tobe_sent,1, MPI_INT, 
		     num_data_tobe_recv, 1, MPI_INT, 
			MPI_COMM_WORLD);
	
	int *recv_disp_vec = malloc(num_procs*sizeof(int));
	int *send_disp_vec = malloc(num_procs*sizeof(int));
		recv_disp_vec[0]=0;
		send_disp_vec[0]=0;
	int total_recv_elems =0;
		total_recv_elems += num_data_tobe_recv[0];
	for(i=1;i<num_procs;i++){
		total_recv_elems += num_data_tobe_recv[i];
		recv_disp_vec[i]=recv_disp_vec[i-1]+num_data_tobe_recv[i-1];
		send_disp_vec[i]=send_disp_vec[i-1]+num_data_tobe_sent[i-1];
	}	


	int *recv_vec = malloc(sizeof(int)*total_recv_elems);
	MPI_Alltoallv(A, num_data_tobe_sent,
            send_disp_vec, MPI_INT,
            recv_vec, num_data_tobe_recv,
            recv_disp_vec, MPI_INT,MPI_COMM_WORLD);

	free(num_data_tobe_sent);
	free(num_data_tobe_recv);
	free(send_disp_vec);
	free(recv_disp_vec);

	
	qsort(recv_vec,total_recv_elems,sizeof(int),cmp);


	int *deficit = malloc(num_procs*sizeof(int));
	int *cumulative_deficit = malloc(num_procs*sizeof(int));

	int my_deficit = total_recv_elems - length;

	MPI_Allgather(&my_deficit, 1, MPI_INT, 
			deficit, 1, MPI_INT, MPI_COMM_WORLD);

	cumulative_deficit[0] = deficit[0];
	//printf(" %d ",deficit[0]); 
	for(i=1;i<num_procs;i++){
		cumulative_deficit[i]= cumulative_deficit[i-1]+deficit[i];
	//	printf(" %d ",deficit[i]);
	}
//	printf(" | myid = %d\n",myid);

	MPI_Request ls_request,lr_request,rs_request,rr_request;
	int lsr=-1,rsr = -1;
	int left_send=0,my_excess=0;
	if(myid > 0){
		left_send = cumulative_deficit[myid-1];
		if(left_send < 0){
			MPI_Isend(recv_vec,(-1)*left_send, MPI_INT, myid-1, tag, MPI_COMM_WORLD, &ls_request);
			lsr = 0;
		}else if(left_send >0){
			//MPI_Irecv(left_recv, left_send, MPI_INT, myid-1, tag, MPI_COMM_WORLD,&lr_request);
			MPI_Irecv(A, left_send, MPI_INT, myid-1, tag, MPI_COMM_WORLD,&lr_request);
			lsr = 1;
		}
	}
	if(myid != num_procs-1){
		my_excess  = cumulative_deficit[myid];
		if(my_excess > 0){
			MPI_Isend(recv_vec + total_recv_elems - my_excess,my_excess,MPI_INT, myid+1, tag, MPI_COMM_WORLD,&rs_request);
			rsr = 0;
		}else if(my_excess <0){
		//	MPI_Irecv(right_recv,(-1)*my_excess,MPI_INT, myid+1, tag, MPI_COMM_WORLD, &rr_request);
			MPI_Irecv(A+length+my_excess,(-1)*my_excess,MPI_INT, myid+1, tag, MPI_COMM_WORLD, &rr_request);
			rsr = 1;
		}
	}

	if(lsr == 0){
		MPI_Wait(&ls_request,&status);
	}else if (lsr == 1) {
		MPI_Wait(&lr_request,&status);
	}
	if(rsr == 0) {
		MPI_Wait(&rs_request,&status);
	}else if(rsr ==1){
		MPI_Wait(&rr_request,&status);
	}


	int t_start= 0,a_start=0;
	if(lsr == 0){
		t_start = -1*left_send;
	}else if(lsr ==1){
		a_start = left_send;
	}
	int t_end = total_recv_elems ,a_end = length;
	if(rsr == 0){
		t_end = total_recv_elems - my_excess;
	}else if(rsr ==1){
		a_end = length+ my_excess;
	}

	int j =0;
	for(i=t_start;i<t_end;i++){
		A[a_start+j] = recv_vec[i];
		j++;
	}
	free(deficit);
	free(cumulative_deficit);
	free(recv_vec);
}


