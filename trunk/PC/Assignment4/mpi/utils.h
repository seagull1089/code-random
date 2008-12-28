#include<stdlib.h>
#include<stdio.h>

#define STEPSIZE 10

int countlines(char *filename);

struct matrixChunk{

	int *row_info;
	int *col_info;
	double *values;
	int row_ptr;
	int cv_ptr;
	int last_access;

} ;


int * get_non_existing_indices(int indices[],int size,
		int s_index,int l_index,int *result_size,int max_index);

void sort(int [],int ,int);
int partition(int [],int ,int);

void printVec(double *d,int N);
double * readValues(int rows,FILE *fp);
void printStruct(struct matrixChunk *mchunk);
struct matrixChunk * readMatrixRows(int rows,FILE *fp);
void free_struct(struct matrixChunk *mchunk);
