#include <mpi.h> 
#include "utils.h" 
#include <libgen.h>


int compare(const void  *a,const void  *b){
        return (*(int *)a - *(int *)b);
}
void computeMIS(int *nodes, int num_nodes, int num_direct_nodes,int *direct_adjnodes, int *direct_node_lengths,int *MIS,MPI_Comm comm) ;
MPI_Status status;
int tag=100;
int N,chunksize,last_proc_chunksize;

int  binSearch(int M_B_ColInfo[],int start , int end, int col) ;
char *bitarray;
int  N_bits;
unsigned int rand_seed ;

int * get_N_randomNumbers(int num){
	int i,*r_vector = malloc(num*sizeof(int));
	for(i=0;i<num;i++){
                  int rand_x = (rand_r(&rand_seed)*rand_r(&rand_seed))%N_bits;
                  while(BITTEST(bitarray,rand_x)){
                          rand_x = (rand_r(&rand_seed)*rand_r(&rand_seed))%N_bits;
                  }
                  BITSET(bitarray,rand_x);
                  r_vector[i] = rand_x;
         }
	return r_vector;
}
int *NL_C_vector,*C_vector,*random_vector ;  
int isMinimumInAdjacent(int num_nodes, int *nodes,int *direct_adjnodes_lengths, int *direct_adjnodes,
			int *indirect_nodes,int *indirect_node_lengths,
			int u_count,int *required_indices,int *required_values,int index);

int *  removeAdjacentCandidateNodes(int num_nodes, int *nodes,int *direct_adjnodes_lengths, int *direct_adjnodes,
			int *indirect_nodes,int *indirect_node_lengths,
			int u_count,int *required_indices,int *required_values,int index,int *removed_nodes,int *r_capacity,int *r_counter);

int main(int argc, char *argv[])
{

  int i,num_procs, myid, name_len;
  char proc_name[MPI_MAX_PROCESSOR_NAME];

  // Initialize MPI
  MPI_Init(&argc, &argv);

  // Obtain the number of processes
  MPI_Comm_size(MPI_COMM_WORLD, &num_procs);

  // Obtain the process id
  MPI_Comm_rank(MPI_COMM_WORLD, &myid);

  // Obtain name of machine process is executing on
  MPI_Get_processor_name(proc_name, &name_len);

  struct AdjacentNodesInfo *mydirect_adjnodes;
  struct IndirectNodesInfo *myindirect_adjnodes;


  char *filename  = argv[1];

	if(myid == 0){
		FILE *fp = fopen(filename,"r");
		N = maxNodes(fp);
		rewind(fp);
		N_bits = N*2;
		bitarray = malloc(BITNSLOTS(N_bits)*sizeof(char));
		memset(bitarray,0,BITNSLOTS(N_bits));
		rand_seed  = time(NULL)*N;
		

		chunksize = N/num_procs;
		last_proc_chunksize = N - ((num_procs-1)*chunksize);
		mydirect_adjnodes = readData(fp,1,chunksize);
		printf(" chunksize = %d \n",chunksize);
		int last_index,first_index,r_size;
		struct AdjacentNodesInfo *tmp_adj;
		r_size = chunksize;
		random_vector = get_N_randomNumbers(r_size);

		for(i=1;i<num_procs;i++){
			first_index = i*chunksize+1;
			last_index = (i+1)*chunksize;
			if(i== num_procs-1){
				last_index = N;
				r_size = last_proc_chunksize;
			}
			tmp_adj = readData(fp,first_index,last_index);
			MPI_Send(&tmp_adj->num_nodes,1,MPI_INT,i,tag,MPI_COMM_WORLD);
			MPI_Send(&tmp_adj->num_direct_nodes,1,MPI_INT,i,tag,MPI_COMM_WORLD);
			MPI_Send(tmp_adj->nodes,tmp_adj->num_nodes,MPI_INT,i,tag,MPI_COMM_WORLD);
			MPI_Send(tmp_adj->direct_adjnodes_lengths,tmp_adj->num_nodes+1,MPI_INT,i,tag,MPI_COMM_WORLD);
			
			MPI_Send(tmp_adj->direct_adjnodes,tmp_adj->num_direct_nodes,MPI_INT,i,tag,MPI_COMM_WORLD);
			int k =0;
			for(k=0;k<tmp_adj->num_nodes;k++){
				//printf(" node : %d range = %d %d \n",tmp_adj->nodes[k],tmp_adj->direct_adjnodes_lengths[k],tmp_adj->direct_adjnodes_lengths[k+1]);

			}	
			free(tmp_adj);
			int *r_vector = get_N_randomNumbers(r_size);
			MPI_Send(r_vector,r_size,MPI_INT,i,tag,MPI_COMM_WORLD);
			free(r_vector);
		}
			

	}else{

		mydirect_adjnodes = malloc(sizeof(struct AdjacentNodesInfo));
		MPI_Recv(&mydirect_adjnodes->num_nodes,1,MPI_INT,0,tag,MPI_COMM_WORLD,&status);
		MPI_Recv(&mydirect_adjnodes->num_direct_nodes,1,MPI_INT,0,tag,MPI_COMM_WORLD,&status);
		mydirect_adjnodes->nodes = malloc(mydirect_adjnodes->num_nodes*sizeof(int));
		mydirect_adjnodes->direct_adjnodes_lengths = malloc((mydirect_adjnodes->num_nodes+1)*sizeof(int));
		mydirect_adjnodes->direct_adjnodes = malloc(mydirect_adjnodes->num_direct_nodes*sizeof(int));
		MPI_Recv(mydirect_adjnodes->nodes,mydirect_adjnodes->num_nodes,MPI_INT,0,tag,MPI_COMM_WORLD,&status);
		MPI_Recv(mydirect_adjnodes->direct_adjnodes_lengths,mydirect_adjnodes->num_nodes+1,MPI_INT,0,tag,MPI_COMM_WORLD,&status);
		MPI_Recv(mydirect_adjnodes->direct_adjnodes,mydirect_adjnodes->num_direct_nodes,MPI_INT,0,tag,MPI_COMM_WORLD,&status);
		random_vector = malloc(sizeof(int)*mydirect_adjnodes->num_nodes);
		MPI_Recv(random_vector,mydirect_adjnodes->num_nodes,MPI_INT,0,tag,MPI_COMM_WORLD,&status);


	}


  	MPI_Bcast(&N,1,MPI_INT,0,MPI_COMM_WORLD);
  	MPI_Bcast(&chunksize,1,MPI_INT,0,MPI_COMM_WORLD);
	
	/* calculate the indirect edges and populate them */

	/* create a temporary adjodes vector to be sorted out */
	int *MIS = malloc(sizeof(int)*mydirect_adjnodes->num_nodes);
	for(i=0;i<mydirect_adjnodes->num_nodes;i++){
		MIS[i]=0;
	}


	MPI_Barrier(MPI_COMM_WORLD);
	double t1,t2;
	t1 = MPI_Wtime();
	
	computeMIS(mydirect_adjnodes->nodes,mydirect_adjnodes->num_nodes,
			mydirect_adjnodes->num_direct_nodes,mydirect_adjnodes->direct_adjnodes,
			mydirect_adjnodes->direct_adjnodes_lengths,MIS,MPI_COMM_WORLD);	
	MPI_Barrier(MPI_COMM_WORLD);
	t2 = MPI_Wtime();

	free(C_vector);
	free(NL_C_vector);		
	if(myid == 0){

		printf("time taken to compute MIS = %f\n",t2-t1);

		char *outfilename = malloc(sizeof(char)*(strlen(filename)+3));
		sprintf(outfilename,"%s.mis",filename);
		outfilename = basename(outfilename);
		FILE *ofp = fopen(outfilename,"w");
		int min_node = 1;
		for(i=0;i<chunksize;i++){
			if(MIS[i] == 1){
			fprintf(ofp,"%d\n",i+min_node);
			}
		}
		min_node += chunksize;
		free(MIS);
		int j,recv_size = chunksize;
		for(j=1;j<num_procs;j++){
			if(j== num_procs -1){
				recv_size = last_proc_chunksize;
			}
			int *tmp_mis = malloc(sizeof(int) * recv_size);
			//MPI_Recv(random_vector,mydirect_adjnodes->num_nodes,MPI_INT,0,tag,MPI_COMM_WORLD,&status);
			MPI_Recv(tmp_mis,last_proc_chunksize,MPI_INT,j,tag,MPI_COMM_WORLD,&status);
			for(i=0;i<recv_size;i++){
				if(tmp_mis[i] == 1){
				fprintf(ofp,"%d\n",i+min_node);
				}
			}
			min_node += recv_size;
			free(tmp_mis);	
		}
		fclose(ofp);

	}else{
		MPI_Send(MIS,mydirect_adjnodes->num_nodes,MPI_INT,0,tag,MPI_COMM_WORLD);
		free(MIS);
	}

	free(mydirect_adjnodes->nodes);
	free(mydirect_adjnodes->direct_adjnodes);
	free(mydirect_adjnodes->direct_adjnodes_lengths);
	free(mydirect_adjnodes);
 
	//printf("Hello world from proceesor %d %d \n",myid,N);
  // Last call to MPI (REQUIRED)
  MPI_Finalize();
}

void computeMIS(int *nodes, int num_nodes, int num_direct_nodes,int *direct_adjnodes, int *direct_node_lengths,int *MIS,MPI_Comm comm) {


	

	int num_procs,myid;
  	MPI_Comm_size(comm, &num_procs);
        MPI_Comm_rank(comm, &myid);
	//printf(" Calling from proc %d \n",myid);

	int min_node = nodes[0];	
	int i,*tmp_adjnodes = malloc(num_direct_nodes*sizeof(int));
	for(i=0;i<num_direct_nodes;i++){
		tmp_adjnodes[i]= direct_adjnodes[i];
	}
	qsort(tmp_adjnodes,num_direct_nodes,sizeof(int),compare);
	int unique =-1;
	int u_count = 0;
	for(i=0;i<num_direct_nodes;i++){
		if(tmp_adjnodes[i] != unique){
			u_count++;
			unique = tmp_adjnodes[i];
		}
	}


	int *u_vertices = malloc(u_count*sizeof(int));
	int *u_lengths 	= malloc(u_count*sizeof(int));
	unique = -1;
	int tot_u_count = u_count;
	u_count =0;
	int v_count =0;
	

	//printf(" tot_u_count = %d myrange = %d to %d and myid = %d \n", tot_u_count,nodes[0],nodes[num_nodes-1],myid);
	for(i=0;i<num_direct_nodes;i++){
		if(tmp_adjnodes[i] != unique){
			u_vertices[u_count] = tmp_adjnodes[i];
			if(u_count != 0){
				u_lengths[u_count-1] = v_count;
				v_count = 0;
			}
			unique = tmp_adjnodes[i];
		//	printf(" -- ucount = %d\n",u_count);
			u_count++;
		}
		v_count++;
	}
	u_lengths[u_count-1] = v_count;

	int j,*rev_lookup = malloc(sizeof(int)*num_direct_nodes);
	int *i_displs 	= malloc(u_count*sizeof(int));
	i_displs[0]=0;
	//printf("\n Unique indices for %d = ",myid);
	for(i=1;i<tot_u_count;i++){
		i_displs[i] = i_displs[i-1]+ u_lengths[i-1];
		//printf(" -- %d",u_vertices[i]);
	}	
	//printf("\n");

	for(i=0;i<num_nodes;i++){
		//printf(" node = %d , indexes = %d , %d \n",nodes[i],direct_node_lengths[i],direct_node_lengths[i+1]);
		for(j=direct_node_lengths[i];j<direct_node_lengths[i+1];j++){
			int i_node = direct_adjnodes[j];
			int i_index = binSearch(u_vertices,0,tot_u_count,i_node);
		//	printf(" searching for  = %d  got index = %d \n" ,i_node,i_index);
			rev_lookup[i_displs[i_index]]=nodes[i];
			i_displs[i_index]++;			
		}
	}	




	// calculate who needs to recieve how many 
	int *num_elems_tobesent = malloc(num_procs*sizeof(int));
	int *num_elems_toberecv = malloc(num_procs*sizeof(int));
	int *send_indices_displs = malloc(num_procs*sizeof(int));	
	int *recv_indices_displs = malloc(num_procs*sizeof(int));	
	int *sdispls_vec = malloc(num_procs*sizeof(int));	
	int *rdispls_vec = malloc(num_procs*sizeof(int));	
	for(i=0;i<num_procs;i++){
		num_elems_tobesent[i] = 0;
		num_elems_toberecv[i] = 0;
		send_indices_displs[i]=0;
		recv_indices_displs[i]=0;
		sdispls_vec[i]=0;
		rdispls_vec[i]=0;
	}


	for(i =0;i<tot_u_count;i++){
		int belongs_to = (u_vertices[i]-1)/chunksize;
		if(belongs_to >=num_procs-1){
			belongs_to = num_procs-1;
		}		
		send_indices_displs[belongs_to]++;
		num_elems_tobesent[belongs_to] += u_lengths[i];
	}

	MPI_Alltoall(send_indices_displs,1,MPI_INT,
			recv_indices_displs,1,MPI_INT,
				comm);

	/*
		recv_buff_size holds the total number of indices 
		recieved from other procs.
	*/

	int recv_buff_size = recv_indices_displs[0];
	for(i=1;i<num_procs;i++){
		sdispls_vec[i]= sdispls_vec[i-1]+send_indices_displs[i-1];
		rdispls_vec[i] = rdispls_vec[i-1]+recv_indices_displs[i-1];
		recv_buff_size += recv_indices_displs[i];
	}

	int *recieved_indices = malloc(recv_buff_size*sizeof(int));
	int *recieved_indices_lengths = malloc(recv_buff_size*sizeof(int));
   	
	MPI_Alltoallv(u_vertices, send_indices_displs,
             sdispls_vec, MPI_INT,
            recieved_indices, recv_indices_displs,
             rdispls_vec, MPI_INT,comm);
		
	MPI_Alltoallv(u_lengths, send_indices_displs,
             sdispls_vec, MPI_INT,
            recieved_indices_lengths, recv_indices_displs,
             rdispls_vec, MPI_INT,comm);
	
	free(u_lengths);


	/*
	printf(" Received indices = ");
	for(i=0;i<recv_buff_size;i++){
		printf(" %d - %d ",recieved_indices[i],recieved_indices_lengths[i]);
	}
	*/

        MPI_Alltoall(num_elems_tobesent,1, MPI_INT,
                      num_elems_toberecv, 1, MPI_INT,
                         comm);

	int tot_recv_len=num_elems_toberecv[0];	
	rdispls_vec[0] = 0;
	sdispls_vec[0] =0;
	for(i=1;i<num_procs;i++){
		tot_recv_len += num_elems_toberecv[i];
		rdispls_vec[i] = rdispls_vec[i-1]+num_elems_toberecv[i-1];
		sdispls_vec[i] = sdispls_vec[i-1]+num_elems_tobesent[i-1];
	}

	int *indirect_recv_vector = malloc(sizeof(int)*tot_recv_len);
	
	MPI_Alltoallv(rev_lookup, num_elems_tobesent,
             sdispls_vec, MPI_INT,
            indirect_recv_vector, num_elems_toberecv,
             rdispls_vec, MPI_INT,comm);

	/* freee */
	free(rev_lookup);

	int *indirect_node_lengths = malloc(sizeof(int)*(num_nodes+1));
	int *tmp_indirect_lengths = malloc(sizeof(int)*(num_nodes));
	int *fill_displs = malloc(sizeof(int)*(num_nodes));
	int *indirect_nodes = malloc(sizeof(int)*tot_recv_len);
	for(i=0;i<num_nodes;i++){
			indirect_node_lengths[i] = 0;
			tmp_indirect_lengths[i]=0;
			fill_displs[i]=0;
	}
	for(i=0;i<recv_buff_size;i++){
		int tmp_node_index  = recieved_indices[i]-min_node;
		tmp_indirect_lengths[tmp_node_index] += recieved_indices_lengths[i];
	}
	
	for(i=1;i<num_nodes;i++){
		indirect_node_lengths[i]= indirect_node_lengths[i-1]+tmp_indirect_lengths[i-1];
		fill_displs[i] = indirect_node_lengths[i];	
	}
	indirect_node_lengths[num_nodes] = tot_recv_len;

	int curr_displs = 0 ;
	
	/*printf("\n indirect recv vector = ");
	for(i=0;i<tot_recv_len;i++){
		printf(" %d ",indirect_recv_vector[i]);
	} 
	printf("\n");
	*/
	for(i=0;i<recv_buff_size;i++){
		int tmp_node_index  = recieved_indices[i]-min_node;
		for(j=curr_displs;j<curr_displs+recieved_indices_lengths[i];j++){
		//	printf(" ##  %d ",indirect_recv_vector[j]);
			indirect_nodes[fill_displs[tmp_node_index]] = indirect_recv_vector[j];
			fill_displs[tmp_node_index]++;
		}	
		curr_displs = curr_displs+recieved_indices_lengths[i];
	}

	free(recieved_indices);
	free(recieved_indices_lengths);
	free(fill_displs);
	free(tmp_indirect_lengths);
	/*	
	printf(" indirect nodes = ");
	for(i=0;i<tot_recv_len;i++){
		printf(" %d ",indirect_nodes[i]);
	} 
	*/
	
	qsort(indirect_recv_vector,tot_recv_len,sizeof(int),compare);
	u_count=0;
	unique = -1;

//	printf(" last element = %d \n",indirect_recv_vector[tot_recv_len-1]);

	for(i=0;i<tot_recv_len;i++){
		if(indirect_recv_vector[i] != unique){
			u_count++;
			unique = indirect_recv_vector[i];
		}
	}

	u_vertices = realloc(u_vertices,sizeof(int)*(tot_u_count+u_count));

	unique=-1;
	
	for(i=0;i<tot_recv_len;i++){
		if(indirect_recv_vector[i] != unique){
			unique = indirect_recv_vector[i];
			u_vertices[tot_u_count] = unique; 
			tot_u_count++;
		}
	}
	free(indirect_recv_vector);	
	qsort(u_vertices,tot_u_count,sizeof(int),compare);
	u_count = 0;
	unique=-1;
	for(i=0;i<tot_u_count;i++){
		if(u_vertices[i] != unique){
			unique=u_vertices[i];
			u_count++;
		}
	} 

	int *required_indices = malloc(sizeof(int)*u_count);
	u_count = 0;
	unique=-1;
	for(i=0;i<tot_u_count;i++){
		if(u_vertices[i] != unique){
			unique=u_vertices[i];
			required_indices[u_count] = unique;
			u_count++;
		}
	} 
	free(u_vertices);
	
	for(i=0;i<num_procs;i++){
	//	num_elems_tobesent[i] = 0;
	//	num_elems_toberecv[i] = 0;
		send_indices_displs[i]=0;
		recv_indices_displs[i]=0;
		sdispls_vec[i]=0;
		rdispls_vec[i]=0;
	}

	for(i =0;i<u_count;i++){
		int belongs_to = (required_indices[i]-1)/chunksize;
		if(belongs_to >=num_procs-1){
			belongs_to = num_procs-1;
		}		
		send_indices_displs[belongs_to]++;
	//	num_elems_tobesent[belongs_to] += u_lengths[i];
	}

	MPI_Alltoall(send_indices_displs,1,MPI_INT,
			recv_indices_displs,1,MPI_INT,
				comm);

	int request_length = recv_indices_displs[0];
	for(i=1;i<num_procs;i++){
		sdispls_vec[i]= sdispls_vec[i-1]+send_indices_displs[i-1];
		rdispls_vec[i] = rdispls_vec[i-1]+recv_indices_displs[i-1];
		request_length += recv_indices_displs[i];
	}
	
	int *requested_indices = malloc(sizeof(int)*request_length);
	int *requested_values = malloc(sizeof(int)*request_length);
	int *required_values = malloc(sizeof(int)*u_count);
	MPI_Alltoallv(required_indices, send_indices_displs,
             sdispls_vec, MPI_INT,
            requested_indices, recv_indices_displs,
             rdispls_vec, MPI_INT,comm);

	
	for(i=0;i<request_length;i++){
		int req_index = binSearch(nodes,0,num_nodes,requested_indices[i]);
		requested_values[i] = random_vector[req_index];
	}


	MPI_Alltoallv(requested_values, recv_indices_displs,
             rdispls_vec, MPI_INT,
            required_values, send_indices_displs,
             sdispls_vec, MPI_INT,comm);

	free(requested_indices);
	free(requested_values);
/*
	printf("\nCurrent Random Numbers :");	
	for(i=0;i<num_nodes;i++){
		printf(" %d - %d ",nodes[i],random_vector[i]);
	}
	printf("\n Gotten values :");
	for(i=0;i<u_count;i++){
		printf(" %d - %d",required_indices[i],required_values[i]);
	}
	printf("\n");
*/
	/*
	for(i=0;i<num_nodes;i++){
		printf( "\n node = %d  ", nodes[i]);
		for(j=indirect_node_lengths[i];j<indirect_node_lengths[i+1];j++){
			printf(" -> %d ",indirect_nodes[j]);
		}
		printf("\n");
	}
	*/

	C_vector = malloc(sizeof(int)*num_nodes);
	for(i=0;i<num_nodes;i++){
		C_vector[i] = 1;
	}

	NL_C_vector = malloc(sizeof(int)*u_count);
	for(i=0;i<u_count;i++){
		NL_C_vector[i] = 1;
	}

	int iteration = 0;
	while(1){
	iteration++;
	int mis_found = 0;	
	int *removed_nodes,r_capacity =1,r_counter=0;
	removed_nodes = malloc(sizeof(int)*r_capacity);
	for(i=0;i<num_nodes;i++){
	
		if(C_vector[i] ==1){
			if(isMinimumInAdjacent(num_nodes,nodes,direct_node_lengths,direct_adjnodes,
						indirect_nodes,indirect_node_lengths,u_count,required_indices,required_values,i)){
				MIS[i] =1;
				C_vector[i] = 0;
				removed_nodes = removeAdjacentCandidateNodes( num_nodes,nodes,direct_node_lengths,direct_adjnodes,
						indirect_nodes,indirect_node_lengths,u_count,required_indices,required_values,i,removed_nodes,&r_capacity,&r_counter);

				//printf("MIS =  %d ",nodes[i]);
				mis_found++;
			}
		}
	
	}

//	printf(" r_capacity = %d , r_counter = %d \n", r_capacity,r_counter);

	qsort(removed_nodes,r_counter,sizeof(int),compare);
	int r_count = 0; 
	unique = -1;
	for(i=0;i<r_counter;i++){
		if(removed_nodes[i] != unique){
			r_count++;
			unique = removed_nodes[i];
		}
	}

	int min_size = (r_count>0)? r_count : 1;
	int *send_remove = malloc(min_size*sizeof(int));
	r_count = 0;
	unique = -1;

	for(i=0;i<r_counter;i++){
		if(removed_nodes[i] != unique){
			unique = removed_nodes[i];
			send_remove[r_count] = unique;
			r_count++;
		}
	}
	free(removed_nodes);

	int *remove_sizes = malloc(num_procs*sizeof(int));

	 MPI_Allgather(&r_count, 1,
            MPI_INT, remove_sizes, 1,
             MPI_INT, comm);

	
	int *remove_displs = malloc(num_procs*sizeof(int));
	remove_displs[0] =0;
	int r_vector_size = remove_sizes[0];
	for(i=1;i<num_procs;i++){
		remove_displs[i] = remove_displs[i-1] + remove_sizes[i-1];
		r_vector_size += remove_sizes[i];
	} 
	int *all_remove_vector = malloc(r_vector_size*sizeof(int));
	
	MPI_Allgatherv(send_remove,r_count,
            MPI_INT, all_remove_vector,remove_sizes,
            remove_displs, MPI_INT, comm);
	
	free(send_remove);

	for(i=0;i<r_vector_size;i++){
		int tmp_node = all_remove_vector[i];
		if(tmp_node >= nodes[0] && tmp_node <= nodes[num_nodes-1]) {
			if(C_vector[tmp_node - nodes[0]] == 1) {
				C_vector[tmp_node -nodes[0]] = 0;
			}	
		}

		
		int tmp_index = binSearch(required_indices,0,u_count,tmp_node);
		if(tmp_index >= 0 && NL_C_vector[tmp_index] == 1) {
			NL_C_vector[tmp_index] = 0;
		}

	}	
	free(remove_displs);
	free(all_remove_vector);
	int cnt_remaining=0;

	
//	int *all_procs_remaining = malloc(sizeof(int)*num_procs);
	for(i=0;i<num_nodes;i++){
		if(C_vector[i] == 1) {
			cnt_remaining++;
		}	
	}

//	printf(" proc = %d , remaing = %d \n",myid,cnt_remaining);
	int total_remaining =0;
/*
	MPI_Allreduce(&cnt_remaining, 1,
            MPI_INT, total_remaining, 1,
             MPI_INT, comm);
*/	
	MPI_Allreduce(&cnt_remaining, &total_remaining,1,
            MPI_INT, MPI_SUM, comm);



	int all_mis_found_count=0; 
	MPI_Reduce(&mis_found, &all_mis_found_count, 1,
            MPI_INT, MPI_SUM, 0, comm);
	if(myid == 0){
		printf("iteration = %d  number of mis nodes found = %d \n",iteration, all_mis_found_count);
	}
	
	if(total_remaining > 0){
	//	printf(" still %d nodes remaing overall \n",total_status);
	}else{
		
		//printf(" Over and Out  %d  \n",total_status);
		break; 
	}

	}// end of while loop .
	/* start your iterations here */


}

int *  addtoRemovedNode(int *removed_nodes,int *r_capacity,int *r_counter,int node){
	if(*r_capacity <= *r_counter){
		*r_capacity = *r_capacity + 1000;
		removed_nodes = realloc(removed_nodes,(*r_capacity)*sizeof(int));
	}
	removed_nodes[*r_counter] = node;
	*r_counter  = *r_counter+1;
	return removed_nodes;

}



int *  removeAdjacentCandidateNodes(int num_nodes, int *nodes,int *direct_adjnodes_lengths, int *direct_adjnodes,
			int *indirect_nodes,int *indirect_node_lengths,
			int u_count,int *required_indices,int *required_values,int index,int *removed_nodes,int *r_capacity,int *r_counter){
		
	int i,j;
	for(i=direct_adjnodes_lengths[index];i<direct_adjnodes_lengths[index+1];i++){
		int tmp_node = direct_adjnodes[i];
		if(tmp_node >= nodes[0] && tmp_node <= nodes[num_nodes-1] && C_vector[tmp_node-nodes[0]] == 1){
			C_vector[tmp_node-nodes[0]] = 0;
			removed_nodes = addtoRemovedNode(removed_nodes,r_capacity,r_counter,tmp_node);
		}else{
			int tmp_index = binSearch(required_indices,0,u_count,tmp_node);
			if(NL_C_vector[tmp_index] == 1) {
				NL_C_vector[tmp_index] = 0;
			removed_nodes	= addtoRemovedNode(removed_nodes,r_capacity,r_counter,tmp_node);
			}
		}
		
		
	}

	for(i= indirect_node_lengths[index];i<indirect_node_lengths[index+1];i++){
		int tmp_node = indirect_nodes[i];
		if(tmp_node >= nodes[0] && tmp_node <= nodes[num_nodes-1] && C_vector[tmp_node-nodes[0]]==1){
			C_vector[tmp_node-nodes[0]] = 0;
			removed_nodes = addtoRemovedNode(removed_nodes,r_capacity,r_counter,tmp_node);
		}else{
			int tmp_index = binSearch(required_indices,0,u_count,tmp_node);
			if(NL_C_vector[tmp_index] == 1) {
				NL_C_vector[tmp_index] = 0;
			removed_nodes	= addtoRemovedNode(removed_nodes,r_capacity,r_counter,tmp_node);
			}
		}
	}
	
	return removed_nodes;

}
int isMinimumInAdjacent(int num_nodes, int *nodes,int *direct_adjnodes_lengths, int *direct_adjnodes,
			int *indirect_nodes,int *indirect_node_lengths,
			int u_count,int *required_indices,int *required_values,int index){
		
	int i,j;
	int val = random_vector[index];
	for(i=direct_adjnodes_lengths[index];i<direct_adjnodes_lengths[index+1];i++){
		int tmp_node = direct_adjnodes[i];
		if(tmp_node >= nodes[0] && tmp_node <= nodes[num_nodes-1] && C_vector[tmp_node-nodes[0]] == 1){
			if(random_vector[tmp_node- nodes[0]] < val){
				return 0;
			}
		}else{
			int tmp_index = binSearch(required_indices,0,u_count,tmp_node);
			if(NL_C_vector[tmp_index] == 1 && required_values[tmp_index] < val){
				return 0;
			}
		}
	}

	for(i= indirect_node_lengths[index];i<indirect_node_lengths[index+1];i++){
		int tmp_node = indirect_nodes[i];
		if(tmp_node >= nodes[0] && tmp_node <= nodes[num_nodes-1] && C_vector[tmp_node-nodes[0]]==1){
			if(random_vector[tmp_node- nodes[0]] < val){
				return 0;
			}
		}else{
			int tmp_index = binSearch(required_indices,0,u_count,tmp_node);
			if(NL_C_vector[tmp_index] == 1 && required_values[tmp_index] < val ){
				return 0;
			}

		}

	}
	
	return 1;

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
