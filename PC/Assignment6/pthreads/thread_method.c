#include"utils.h"
#include<pthread.h>
#include <libgen.h>

int max_node,num_threads;
struct AdjacentNodesInfo *adjnodesinfo;
struct IndirectNodesInfo *indirectnodesinfo; 
int *C_Vector;
int *MIS; 
int *rand_vector; 
int            *thread_ids;
pthread_t      *tids;


void  *thread_function(void *arg);
int num_threads,max_node;
int chunksize; 
void markIndepdendet(int index,int val);
int isMinimumInAdjacent(int index,int val);
void synchronizeIteration();
int CandidateListLength(){
	int i,cnt =0;
	for(i=0;i<max_node;i++){
		if(C_Vector[i] ==1){
			return 1;
		}
	}
	return 0;
}
int main(int argc,char **argv){
	if(argc < 3) {
		printf("Usage : progname  filename num_threads\n");
		exit(1);
	}

	char *filename = argv[1];
	char *outfilename = malloc(sizeof(char)*(strlen(filename)+3));
	sprintf(outfilename,"%s.mis",filename);
	outfilename = basename(outfilename);
	num_threads = atoi(argv[2]);
	FILE *fp = fopen(filename,"r");
	max_node = maxNodes(fp);
	rewind(fp);
		
	adjnodesinfo = readData(fp,1,max_node); 
	fclose(fp);
	
	chunksize = max_node/num_threads;
	indirectnodesinfo = getIndirectNodes(adjnodesinfo);
	
	printf(" max nodes = %d \n",max_node); 
	MIS = malloc(max_node*sizeof(int));
	C_Vector = malloc(max_node*sizeof(int)); 


	int i=0,j=0;
	for(i=0;i<max_node;i++){
		MIS[i]= 0;
		C_Vector[i] =1;
	} 	
	rand_vector = getUniqueRandomNumbers(max_node);



	thread_ids = (int *) malloc(num_threads * sizeof(int));
	tids = (pthread_t *) malloc(num_threads * sizeof(pthread_t));
        
	for (i = 0; i < num_threads; i++) {
		 thread_ids[i]=i;
                 pthread_create(&tids[i], NULL,
                                       thread_function, &thread_ids[i]);
        }

	int mis_marked = 0;
	int iteration_no = 0; 

	struct timeval t1,t2;
	gettimeofday(&t1,NULL);

	while(1){
		for (i = 0; i < num_threads; i++) {
        	         pthread_create(&tids[i], NULL,
                                       thread_function, &thread_ids[i]);
        	}

    	   	 for (i = 0; i < num_threads; i++) {
        	          (void) pthread_join(tids[i], NULL);
         	}
		iteration_no++;
		synchronizeIteration();
		int curr_mis_marked_cnt = 0; 
		for(i=0;i<max_node;i++){
			if(MIS[i] == 1) {
				curr_mis_marked_cnt++;
			}
		}

		printf("iteration no = %d , no of mis nodes marked = %d \n",iteration_no,curr_mis_marked_cnt - mis_marked);
		mis_marked = curr_mis_marked_cnt;
		if(CandidateListLength() == 0){
			break; 
		}
	
	}

	gettimeofday(&t2,NULL);
        int total_seconds = (t2.tv_sec - t1.tv_sec)*1000000 +
                            (t2.tv_usec - t1.tv_usec);

        printf("total time taken to calculate MIS  = %f\n",total_seconds/1000000.0);


	FILE *ofp = fopen(outfilename,"w");
//	printf(" MIS nodes are : ");
	for(i=0;i<max_node;i++){
		if(MIS[i] == 1) {
		//	printf(" %d \n", adjnodesinfo->nodes[i]);
			fprintf(ofp,"%d\n",adjnodesinfo->nodes[i]);
		}
	}
//	printf("\n");
	fclose(ofp);
	
	free(thread_ids);
	if(rand_vector)	
	free(rand_vector);
	if(adjnodesinfo){
	if(adjnodesinfo->nodes)
	free(adjnodesinfo->nodes);
	if(adjnodesinfo->direct_adjnodes_lengths)
	free(adjnodesinfo->direct_adjnodes_lengths);
	if(adjnodesinfo->direct_adjnodes)
	free(adjnodesinfo->direct_adjnodes);
	free(adjnodesinfo);
	}
	if(indirectnodesinfo){
	if(indirectnodesinfo->indirect_adjnodes)
	free(indirectnodesinfo->indirect_adjnodes);
	if(indirectnodesinfo->indirect_adjnodes_lengths)
	free(indirectnodesinfo->indirect_adjnodes_lengths);
	free(indirectnodesinfo);
	}
	if(MIS)
	free(MIS);
	if(C_Vector)
	free(C_Vector);
}

int isMinimumInAdjacent(int index,int val){
	int i = 0;
//	printf(" Given node = %d , val = %d , ", index+1,val);
	for(i=adjnodesinfo->direct_adjnodes_lengths[index];
		i<adjnodesinfo->direct_adjnodes_lengths[index+1];i++){
		int node = adjnodesinfo->direct_adjnodes[i];
	
		if((C_Vector[node-1] == 1 && rand_vector[node-1] < val)   ||  (MIS[node-1] == 1)) {
	//		printf(" adj node for %d :  %d , value = %d \n",index+1,adjnodesinfo->nodes[node-1],rand_vector[node-1]);
			return 0;
		}
	}
	for(i=indirectnodesinfo->indirect_adjnodes_lengths[index];
		i<indirectnodesinfo->indirect_adjnodes_lengths[index+1];i++){
		int node = indirectnodesinfo->indirect_adjnodes[i];
	
		if((C_Vector[node-1] == 1 && rand_vector[node-1] < val)   || (  MIS[node-1] == 1)) {
//			printf(" adj node for %d :  %d , value = %d \n",index+1,adjnodesinfo->nodes[node-1],rand_vector[node-1]);
			return 0;
		}
	}
	
	return 1;
} 

void markIndepdendet(int index,int val){

	int i = 0;
	C_Vector[index] =0;
	for(i=adjnodesinfo->direct_adjnodes_lengths[index];
		i<adjnodesinfo->direct_adjnodes_lengths[index+1];i++){
		int node = adjnodesinfo->direct_adjnodes[i];
		if(C_Vector[node-1] == 1 )
			C_Vector[node-1] = 0;
	}
	for(i=indirectnodesinfo->indirect_adjnodes_lengths[index];
		i<indirectnodesinfo->indirect_adjnodes_lengths[index+1];i++){
		int node = indirectnodesinfo->indirect_adjnodes[i];
		if(C_Vector[node-1] == 1 )
			C_Vector[node-1] = 0;
	}
	

}

void  *thread_function(void *arg)
{
        int  *myid = (int *) arg;
	int l_index = (*myid)*chunksize;
	int r_index = (*myid+1)*chunksize;
	if(*myid == num_threads -1) {
		r_index = max_node;
	}
//	printf("\n This is from thread %d min_index = %d max_index = %d  max_node = %d \n",*myid,l_index,r_index,max_node);	
	
	int i=0,j=0;
	
	for(i=l_index,j=0; i< r_index; i++,j++){
		if(C_Vector[i] == 1){
			if( isMinimumInAdjacent(i,rand_vector[i])){
				//printf(" Marking %d as independent ..\n",i+1);
				MIS[i] = 1;
			}
		} 
	} 
	
        return (void *) 0;
}

void  synchronizeIteration() { 

	int i = 0; 
	int count = 0;
	for(i=0;i<max_node;i++){
		if(MIS[i] == 1  &&  C_Vector[i] == 1) {	
			markIndepdendet(i,rand_vector[i]);		
			count++;		
		}
	}

//	printf("Marked %d nodes as independent  \n");


} 
