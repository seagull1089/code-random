#include <mpi.h> 
#include "utils.h" 

double get_B_Value(double B[],double Missing_B[],
		int M_B_ColInfo[],int M_Bsize,
			int start,int end,int row_num,int *last_B_index);

int  binSearch(int M_B_ColInfo[],int start , int end, int col) ;
int main(int argc, char *argv[])
{

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

  char *A_file = argv[1];
  char *B_file = argv[2];;


  int tag=100;
  int N,rows=0;
 

/* this has the vectorB values for the corresponding processor */ 
  double *vec_B;	
/* this has the corresponding row values for the corresponding processor */ 
  struct matrixChunk *myrows;
  
  if(myid == 0) {
	N = countlines(B_file);
	int rows_per_process = N/num_procs;
	int i;
//	printf(" lines = %d , rank = %d rows_per_process = %d \n",N,myid,rows_per_process);
	rows = rows_per_process;
	FILE *bv = fopen(B_file,"r");
	vec_B = readValues(rows,bv);
	for(i = 1; i<num_procs;i++){
		if(i==num_procs-1){
			rows_per_process = N-((num_procs-1)*(N/num_procs));
		}
		MPI_Send(&rows_per_process,1,MPI_INT,i,tag,MPI_COMM_WORLD);		
		double *tmp_vecB  = readValues(rows_per_process,bv);
		MPI_Send(tmp_vecB,rows_per_process,MPI_DOUBLE,i,tag,MPI_COMM_WORLD);
		free(tmp_vecB);
	}
	fclose(bv);
	
	FILE *fp = fopen(A_file,"r");

	myrows = readMatrixRows(rows,fp);
	rows_per_process = N/num_procs;
	for(i=1;i<num_procs;i++){
		if(i==num_procs-1){
			rows_per_process = N-((num_procs-1)*(N/num_procs));
			//printf("rows_per_process = %d \n",rows_per_process);
		}
		struct matrixChunk *mchunk = readMatrixRows(rows_per_process,fp);
		//Send out the structure information 

		MPI_Send(&mchunk->row_ptr,1,MPI_INT,i,tag,MPI_COMM_WORLD);
		MPI_Send(&mchunk->cv_ptr,1,MPI_INT,i,tag,MPI_COMM_WORLD);

		MPI_Send(mchunk->row_info,mchunk->row_ptr,MPI_INT,i,tag,MPI_COMM_WORLD);
		MPI_Send(mchunk->col_info,mchunk->cv_ptr,MPI_INT,i,tag,MPI_COMM_WORLD);
		MPI_Send(mchunk->values,mchunk->cv_ptr,MPI_DOUBLE,i,tag,MPI_COMM_WORLD);
		
		free_struct(mchunk);
	}
	fclose(fp);
	
  }else{
		MPI_Recv(&rows,1,MPI_INT,0,tag,MPI_COMM_WORLD,&status);	
		vec_B = malloc(rows*sizeof(double));
		MPI_Recv(vec_B,rows,MPI_DOUBLE,0,tag,MPI_COMM_WORLD,&status);

		myrows = malloc(sizeof(struct matrixChunk));
		MPI_Recv(&myrows->row_ptr,1,MPI_INT,0,tag,MPI_COMM_WORLD,&status);
		MPI_Recv(&myrows->cv_ptr ,1,MPI_INT,0,tag,MPI_COMM_WORLD,&status);
		myrows->row_info=malloc(myrows->row_ptr*sizeof(int));
		myrows->col_info=malloc(myrows->cv_ptr*sizeof(int));
		myrows->values  =malloc(myrows->cv_ptr*sizeof(double));
		
		MPI_Recv(myrows->row_info,myrows->row_ptr,MPI_INT,0,tag,MPI_COMM_WORLD,&status);
		MPI_Recv(myrows->col_info,myrows->cv_ptr, MPI_INT,0,tag,MPI_COMM_WORLD,&status);
		MPI_Recv(myrows->values,myrows->cv_ptr, MPI_DOUBLE,0,tag,MPI_COMM_WORLD,&status);
//		printf(" Recieved all info  \n");		

	}

// Broadcast the total size to every one.
	MPI_Bcast(&N,1,MPI_INT,0,MPI_COMM_WORLD);
//	printStruct(myrows);
// END  of distribution of matrices .

double T1,T2,T3;
// take the time here. Put a barrier and then take time . 
MPI_Barrier(MPI_COMM_WORLD);
if(myid == 0){
	T1  =  MPI_Wtime();
	//printf("Start of step 2 : time = %lf \n",T1); 
}

// Now Identify what is not there for computation in B.
	
	int start_row_index,last_row_index;
	start_row_index = (N/num_procs)*(myid); //myrows->row_info[0];
	last_row_index  = start_row_index + rows-1;//myrows->row_info[2*(rows-1)];

/* identify the unique required B indices 
   that we need from other processors.
   Get them in sorted order. We need the values in sorted fashion 
   to make further calculations simpler. 
*/

	int num_required_B_elems;
	int *required_B_indices = 
		get_non_existing_indices(myrows->col_info,
					myrows->cv_ptr,
					start_row_index,last_row_index,
					&num_required_B_elems,N);	

//	printf("Got the required indices vector , id = %d \n",myid);

// Identify which process would be having the missing elements
// Count how many elements do you need from the each of the processors.

/*
	The variables here are not so intuitive.
	Have to bear with that.  Don't get confused.
*/
	int i;
	int *rcnt_per_process  = malloc(num_procs*sizeof(int));
	int *recvbuf = malloc(num_procs*sizeof(int));
	int *sdispls = malloc(num_procs*sizeof(int));
	int *rdispls = malloc(num_procs*sizeof(int));

	for(i=0;i<num_procs;i++){
		rcnt_per_process[i]=0;
		recvbuf[i]=0;
		sdispls[i]=0;
		rdispls[i]=0;
	}

	int num_data_tobe_sent=0,num_data_expected=0;
	int segsize = N/num_procs;
	for(i=0;i<num_required_B_elems;i++){
		int related_pid = required_B_indices[i]/segsize;
		if(related_pid >=num_procs){
			related_pid = num_procs-1;
		}
		rcnt_per_process[related_pid]++;	
		num_data_expected++;
	}
/*
	if(myid ==1) {
		for(i=0;i<num_required_B_elems;i++){
			printf(" --%d ",required_B_indices[i]);
		}
		printf("id  =  %d\n",myid); 
	}
*/


// tell others how much data you are expecting from others.
	MPI_Alltoall(rcnt_per_process, 1,
            MPI_INT, recvbuf,1,
            MPI_INT, MPI_COMM_WORLD);
	
	for(i=0;i<num_procs;i++){
		num_data_tobe_sent +=recvbuf[i];
//		printf("%d | %d ",recvbuf[i],rcnt_per_process[i]);
	}
//	printf("-->id = %d\n",myid);

	int *to_be_sent_cols_info = malloc(num_data_tobe_sent*sizeof(int));

	for(i=1;i<num_procs;i++){
		sdispls[i]=sdispls[i-1]+rcnt_per_process[i-1];	
		rdispls[i]=rdispls[i-1]+recvbuf[i-1];
//		printf("%d %d --", sdispls[i],rdispls[i]);	
	}	


/* 
	Send the required indices to other processors.
		rcnt_per_process : has the info on which how many elements a processor has to get 
		sdispls : has the info on from which location of required_B_indices to start reading from 
			 to send the data to a particular process.
		recvbuf : has the info on how much the processor is expecting the data from other processors. 
		rdipls : has the info of displacement in to_be_sent_cols_info to store the data.
*/

	MPI_Alltoallv(required_B_indices, rcnt_per_process,
            sdispls, MPI_INT,
            to_be_sent_cols_info, recvbuf,
            rdispls, MPI_INT,MPI_COMM_WORLD);

MPI_Barrier(MPI_COMM_WORLD);
if(myid == 0){
	T2  =  MPI_Wtime();
//	printf("Start of step 5 : time = %lf \n",T2); 
}

	double *to_be_sent_Bvals = malloc(sizeof(double)*num_data_tobe_sent);
	double *missing_Bvals = malloc(sizeof(double)*num_required_B_elems);

	for(i=0;i<num_data_tobe_sent;i++){
			to_be_sent_Bvals[i]=vec_B[to_be_sent_cols_info[i]-start_row_index];
			//printf(" %d -%lf ,",to_be_sent_cols_info[i],to_be_sent_Bvals[i]); 
	}
//	printf("\n");

/*
	Send the back the values which other processors are asking for. and recieve the elements which we have requested earlier.
	to_be_sent_Bvals : has the B vector values which others are asking for 
	recvbuf : has the info on which processor has to recieve how many elements.
	rdispls : has the displacement info to tell from where in to_be_sent_Bvals to start reading from 
		  for a particular processor.
	missing_Bvals : is the buffer where we store all the incoming data from other processors. This is the B values 
			which we have earlier requested for.
	rcnt_per_process : how many from each processor we are expecting 
	sdispls: displacement in missing_Bvals from where the incoming elements of a processor are stored.
*/
	MPI_Alltoallv(to_be_sent_Bvals, recvbuf,
            rdispls, MPI_DOUBLE,
            missing_Bvals, rcnt_per_process,
            sdispls, MPI_DOUBLE,MPI_COMM_WORLD);

/*	
	for(i=0;i<num_required_B_elems;i++){

		printf(" %d - %lf \n",required_B_indices[i],missing_Bvals[i]); 
	}
	printf("\n");

*/
MPI_Barrier(MPI_COMM_WORLD);
if(myid == 0){
	T3  =  MPI_Wtime();
//	printf("Start of step 6 : time = %lf \n",T3); 
}


	double *result = malloc(sizeof(double)*rows);
	for(i=0;i<rows;i++){
		result[i]=0;
	}
	
//	printf(" start index = %d , last index = %d ,myid = %d \n",start_row_index,last_row_index,myid);
	int j=0,last_j_index=0,r_index=0;
	for(i=0;i<myrows->row_ptr;i=i+2){
		int i_cols = myrows->row_info[i+1];
		int r_index = myrows->row_info[i]-start_row_index;
	//	printf(" i = %d , myid = %d \n",r_index,myid);
		
		// The columns indexes  of a given row are in increasing order. 
		// so while calculating the product , we any way have to start from lower col index 
		// and we don't have to access the lower col index values. so just keep track of where we
                // had looked for a particular j value. Next time , you just have to start from this index. 
		// You don't have to scan the entire missing_Bvals vector for the next value.
		
		//The order of accessing B values is O(1) in this case of matrix multiplicarion.
		//

		int b_val_index = 0;
		for(j=last_j_index;j<last_j_index+i_cols;j++){
			int j_col = myrows->col_info[j];
			result[r_index] += myrows->values[j]*
					get_B_Value(vec_B,missing_Bvals,
						required_B_indices,num_required_B_elems,
						start_row_index,last_row_index,j_col,&b_val_index);		
		}
		last_j_index +=i_cols;

	//	r_index++;
	}
MPI_Barrier(MPI_COMM_WORLD);
if(myid == 0){
	double T4  =  MPI_Wtime();
	printf("End of calculations\n");
	printf("Start of step 6 : time = %lf \n",T4); 
	printf("Time taken to perfrom steps 2 thru 6 = %lf  \n",(T4-T1));
	printf("Time taken to perfrom step  5 thru 6 = %lf  \n",(T4-T2));
//	printf("Time taken to perfrom step  6 = %lf  \n",(T4-T3));
	

}


	if(myid == 0){
		/*
			Aggregate all the results and write them down to disk. 
			No more high logic. 
		*/
		FILE *fres = fopen("out.txt","w");
		for(i=0;i<rows;i++){
			fprintf(fres,"%lf\n",result[i]);
			//printf(" final -- %lf id = 0\n",result[i]);
		}
		int recv_size = rows;
		for(i=1;i<num_procs;i++){
			if(i==num_procs-1){
				recv_size = N - (num_procs-1)*rows;
			}
			double *tmp_result = malloc(recv_size*sizeof(double));
			MPI_Recv(tmp_result,recv_size,MPI_DOUBLE,i,tag,MPI_COMM_WORLD,&status);
			int k =0;
			for(k=0;k<recv_size;k++){
				fprintf(fres,"%lf\n",tmp_result[k]);
				//printf(" final -- %lf id = %d \n",result[k],k);
			}
			free(tmp_result);
		}

		fclose(fres);

	}else{
		MPI_Send(result,rows,MPI_DOUBLE,0,tag,MPI_COMM_WORLD);		
	}

	//MPI_Gather(&myid,1,MPI_INT,tmp,1,MPI_INT,0,MPI_COMM_WORLD);
	//printVec(vec_B,rows);
	/*
		free all the mallocs
	*/
	free(to_be_sent_cols_info);
	free(to_be_sent_Bvals);
	free(required_B_indices);
	free(missing_Bvals);
	free(result);
	free(vec_B);
	free(rcnt_per_process);
	free(recvbuf);
	free(sdispls);
	free(rdispls);
	free_struct(myrows);
  // Last call to MPI (REQUIRED)
  MPI_Finalize();
}


double get_B_Value(double B[],double Missing_B[],
		int M_B_ColInfo[],int M_Bsize,
			int start,int end,int row_num,int *last_B_index){

	/*
		last_B_index is the last index we have searched in Missing_B vector. 
		the next element we are searching for would be either at this index or 
		after this index only.
	*/	
	double res = 0;
	int i=0;

	
	if(row_num >= start && row_num <=end){
		res= B[row_num - start];
	}else{
		/*
		for(i=*last_B_index;i<M_Bsize;i++){
			if(M_B_ColInfo[i]==row_num){
				res = Missing_B[i];
				*last_B_index = i;
				//printf(" %d %lf \n",row_num,res);
				break;	
			}
		}	
		*/
		int k  = binSearch(M_B_ColInfo,0,M_Bsize-1,row_num);
		res = (k>=0) ? Missing_B[k] : 0;
	}
	if(res == 0) {
		printf(" res = 0.000 for j = %d\n",row_num);
	}
	return res;
	
}

int  binSearch(int M_B_ColInfo[],int start , int end, int col) {
	int i = (start + end)/2;
	if(M_B_ColInfo[i] == col) {
		return i;
	}else if(start >= end) {
		return -1;
	}
	else if(M_B_ColInfo[i] >  col ) {
		return binSearch(M_B_ColInfo,start,i-1,col);
	} else if(M_B_ColInfo[i] < col){
		return binSearch(M_B_ColInfo,i+1,end,col);
	}
	
	return -1;
}
