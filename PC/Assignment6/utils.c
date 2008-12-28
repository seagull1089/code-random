#include"utils.h"

int * getUniqueRandomNumbers(int num){

	int *random_numbers = malloc(sizeof(int)*num);
	int n_bits = num*2; 
	char *bitarray = malloc(BITNSLOTS(n_bits)*sizeof(char));
        memset(bitarray,0,BITNSLOTS(n_bits));
	int i = 0; 
	unsigned int x = time(NULL)*num;
	for(i=0;i<num;i++){
		
		//int rand_x = (rand_r(&x)*rand_r(&x))%n_bits;
		int rand_x = ((int) rand_r(&x))%n_bits;
		while(BITTEST(bitarray,rand_x)){
			rand_x = (rand_r(&x))%n_bits;
		}
		BITSET(bitarray,rand_x);
		random_numbers[i] = rand_x;
	}

/*
			if(!BITTEST(bitarray,indices[i])){
                                BITSET(bitarray,indices[i]);
                                u_indices_count++;
                        }
*/
	return random_numbers;
} 

int max(int a , int b) { 
	return (a>b)?a:b;
} 
int maxNodes(FILE *fp) {
	int u,v,max_nodes=0;
	while(fscanf(fp,"%d %d",&u,&v)>0){
		max_nodes = max(max_nodes,max(u,v));
	}
	return max_nodes;
}
struct AdjacentNodesInfo *readData(FILE *fp, int min_v, int max_v){
	
	struct AdjacentNodesInfo *adjnodesinfo;
	adjnodesinfo = malloc(sizeof(struct AdjacentNodesInfo));
	int u,v;
	int  dadj_capacity = 100;
 	int i=0;

	//printf(" min vertex = %d , max_vertex = %d \n",min_v,max_v);

	int num_vertices = max_v-min_v+1;
	adjnodesinfo->nodes = malloc(num_vertices*sizeof(int));

	for(i=0;i<num_vertices;i++){
		adjnodesinfo->nodes[i]= min_v+i;
	} 
	
	adjnodesinfo->direct_adjnodes_lengths = malloc((num_vertices+1)*sizeof(int));
	adjnodesinfo->direct_adjnodes = malloc(dadj_capacity*sizeof(int));

	adjnodesinfo->num_direct_nodes = 0; 
	int  num_lines= 0;
	int last_vertex_id = min_v-1; 
	int file_offset = ftell(fp);

	adjnodesinfo->direct_adjnodes_lengths[0]=0;	

	while(fscanf(fp,"%d %d",&u,&v)>0){
		if(u != last_vertex_id){
	//		printf(" last vertex = %d , curr = %d , num lines = %d \n",last_vertex_id,u,num_lines);
			for(i=last_vertex_id;i<u;i++){
	//			printf("assinging %d  with  %d \n",i-min_v+1,num_lines);
				adjnodesinfo->direct_adjnodes_lengths[i-min_v+1] = num_lines;
			}
			last_vertex_id = u; 
			if(last_vertex_id > max_v) {
				fseek(fp,file_offset,SEEK_SET);
				break; 
			} 

		} 		

		//printf("%d %d\n",u , v);
		if(num_lines >= dadj_capacity -1){
			dadj_capacity = dadj_capacity+1000;
			adjnodesinfo->direct_adjnodes = realloc(adjnodesinfo->direct_adjnodes,dadj_capacity*sizeof(int));
		}
		file_offset = ftell(fp);

		adjnodesinfo->direct_adjnodes[num_lines] = v;
		num_lines++;
	}

	adjnodesinfo->num_nodes = num_vertices;
	// i<u .. 
//	printf(" last_vertex_id = %d , minv-1 = %d , numlines = %d \n",last_vertex_id,min_v,num_lines);
	for(i=last_vertex_id;i<max_v+1;i++){
		adjnodesinfo->direct_adjnodes_lengths[i-min_v+1] = num_lines;
	}
	//adjnodesinfo->direct_adjnodes_lengths[num_vertices] = num_lines;
	adjnodesinfo->num_direct_nodes = num_lines;
	return adjnodesinfo;
}

struct IndirectNodesInfo *getIndirectNodes(struct AdjacentNodesInfo *adjnodesinfo){
	struct IndirectNodesInfo *indirectnodes; 
	indirectnodes = malloc(sizeof(struct IndirectNodesInfo));
	indirectnodes->num_indirect_nodes = adjnodesinfo->num_direct_nodes;
	int num_indirect_nodes = adjnodesinfo->num_direct_nodes;
	int num_nodes = adjnodesinfo->num_nodes;

	indirectnodes->indirect_adjnodes = malloc(sizeof(int)*num_indirect_nodes);
	int *i_counters = malloc(adjnodesinfo->num_nodes* sizeof(int));
	int *i_displs = malloc(adjnodesinfo->num_nodes*sizeof(int));
	int *i_curdispls = malloc(adjnodesinfo->num_nodes*sizeof(int));
	int j=0,i = 0; 
	for(i=0;i<adjnodesinfo->num_nodes;i++){
		i_counters[i]=0;
		i_displs[i] = 0;
		i_curdispls[i]=0;	
	}

	
	for(i=0;i<num_nodes; i++){
		for(j=adjnodesinfo->direct_adjnodes_lengths[i];j<adjnodesinfo->direct_adjnodes_lengths[i+1];j++){
			i_counters[adjnodesinfo->direct_adjnodes[j]-1]++;
		}
	}
	for(i=1;i<num_nodes+1;i++){
		i_displs[i] = i_counters[i-1]+i_displs[i-1];
		i_curdispls[i] = i_displs[i];
	}
	for(i=0;i<adjnodesinfo->num_nodes; i++){
		for(j=adjnodesinfo->direct_adjnodes_lengths[i];j<adjnodesinfo->direct_adjnodes_lengths[i+1];j++){
			indirectnodes->indirect_adjnodes[i_curdispls[adjnodesinfo->direct_adjnodes[j]-1]] = adjnodesinfo->nodes[i];
			i_curdispls[adjnodesinfo->direct_adjnodes[j]-1]++;
		}
	}

	indirectnodes->indirect_adjnodes_lengths = i_displs;
	free(i_counters);
	free(i_curdispls);	
	
	return indirectnodes;

}


