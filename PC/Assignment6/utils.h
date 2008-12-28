#include<stdio.h>
#include<stdlib.h>
#include<limits.h>
#include<string.h>
#define BITMASK(b) (1<<((b)%CHAR_BIT))
#define BITSLOT(b) ((b)/CHAR_BIT)
#define BITSET(a,b) ((a)[BITSLOT(b)] |= BITMASK(b))
#define BITCLEAR(a,b) ((a)[BITSLOT(b)] &= ~BITMASK(b))
#define BITTEST(a,b) ((a)[BITSLOT(b)] & BITMASK(b))
#define BITNSLOTS(nb) ((nb + CHAR_BIT -1 )/ CHAR_BIT)



struct AdjacentNodesInfo {

	int *nodes;
	int num_nodes; 
	int num_direct_nodes; 
	int *direct_adjnodes_lengths;
	int *direct_adjnodes;
	/*
	int num_indirect_nodes; 
	int *indirect_adjnodes_lengths; 
	int *indirect_adjnodes; 
	*/
};


struct IndirectNodesInfo{
	int num_indirect_nodes;
	int *indirect_adjnodes_lengths;
	int *indirect_adjnodes;

};

struct MinRandIndex{
	int min_index;
	int min_rand; 
};
int * getUniqueRandomNumbers(int num);
int maxNodes(FILE *fp);
int max(int a , int b) ; 
struct AdjacentNodesInfo *readData(FILE *fp, int min_v,int max_v);
struct IndirectNodesInfo *getIndirectNodes(struct AdjacentNodesInfo *adjnodesinfo);
