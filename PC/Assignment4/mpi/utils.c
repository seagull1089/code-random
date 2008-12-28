#include"utils.h"
#include<string.h>

/*
	Bit array manipulation macros taken from 
	http://c-faq.com/misc/bitsets.html
*/

#include<limits.h>
#define BITMASK(b) (1<<((b)%CHAR_BIT))
#define BITSLOT(b) ((b)/CHAR_BIT)
#define BITSET(a,b) ((a)[BITSLOT(b)] |= BITMASK(b))
#define BITCLEAR(a,b) ((a)[BITSLOT(b)] &= ~BITMASK(b))
#define BITTEST(a,b) ((a)[BITSLOT(b)] & BITMASK(b))
#define BITNSLOTS(nb) ((nb + CHAR_BIT -1 )/ CHAR_BIT)


int countlines(char *filename){
	FILE *fp;
	double df;
	int count=0;
	fp = fopen(filename,"r");
	while(fscanf(fp,"%lf",&df)>0){
//	printf("-- %lf\n",df);
	count++;
	}
	fclose(fp);
	return count;
}

void printVec(double *d,int N){
	int i;
	for(i=0;i<N;i++){
		printf(" %lf , ",d[i]);
	}
	printf("\n");
}


int * get_non_existing_indices(int indices[],int size,
				int s_index,int l_index,
					int *result_size,int max_index){

	int i=0,j=0;
	int u_indices_count =0;
	/*
	int *tmp_indices = malloc(max_index*sizeof(int));
	for(i=0;i<max_index;i++){
		tmp_indices[i]=0;
	}
	*/
	char *bitarray = malloc(BITNSLOTS(max_index)*sizeof(char));
	memset(bitarray,0,BITNSLOTS(max_index));
	for(i=0;i<size;i++){
		if(indices[i] < s_index ||  indices[i]>l_index){
			if(!BITTEST(bitarray,indices[i])){
				BITSET(bitarray,indices[i]);
				u_indices_count++;
			}
			/*
			if(tmp_indices[indices[i]] ==0){
				u_indices_count++;
				tmp_indices[indices[i]]=1;
			}
			*/		
		}

	}
	//printf("\n");
	int *u_indices = malloc(u_indices_count*sizeof(int));
	*result_size = u_indices_count;
	
	j = 0;
	for(i=0;i<max_index;i++){
		if(BITTEST(bitarray,i)){		
//		if(tmp_indices[i]==1){
			u_indices[j]=i;
			j++;
		}
	}
	free(bitarray);
	//free(tmp_indices);
//	printf(" Final  size = %d\n",u_indices_count); 
//	sort(u_indices,0,u_indices_count-1);
	return u_indices;
}


void sort(int A[],int beg, int end){
	if(beg<end){
		int k;
		k = partition(A,beg,end);
		sort(A,beg,k-1);
		sort(A,k+1,end);
	}
}

int 
partition(int A[], int p, int r)
{
        int             x = A[r];
        int             i = p - 1;
        int             j, tmp;
        for (j = p; j <= r - 1; j++) {
                if (A[j] < x) {
                        i = i + 1;
                        tmp = A[i];
                        A[i] = A[j];
                        A[j] = tmp;
                }   
        }   

        A[r] = A[i + 1]; 
        A[i + 1] = x;
        return i + 1;
}


double * readValues(int rows,FILE *fp){
	int count=0;
	double *values = malloc(rows*sizeof(double));
	while(count <rows){
		fscanf(fp,"%lf",&values[count]);
		count++;
	}
	return values;
}

void free_struct(struct matrixChunk *mchunk){
	if(mchunk !=NULL){
		if(mchunk->row_info!=NULL){
			free(mchunk->row_info);
		}
		if(mchunk->col_info!=NULL) {
			free(mchunk->col_info);
		}
		if(mchunk->values!=NULL){
			free(mchunk->values);
		}
		free(mchunk);
	}
}
struct matrixChunk *  readMatrixRows(int rows,FILE *fp){
/* Assume some value to initially allocate the memory */
	struct matrixChunk *mchunk;
	mchunk = malloc(sizeof(struct matrixChunk));
	int r_count=0,r_start;
	int __start = 0;
	int i = 0, j=0;
	double val=0;
	int last_row =0,r_elements_count=0;
 
	int current_file_pos = ftell(fp);
	
	int curr_vec_size = 10;//STEPSIZE;

	mchunk->row_info = (int *)malloc(rows*2*sizeof(int));
	mchunk->values = (double *)malloc(curr_vec_size*sizeof(double)); 
	mchunk->col_info = (int *)malloc(curr_vec_size*sizeof(int));
	mchunk->row_ptr=0;
	mchunk->cv_ptr=0;
	while(fscanf(fp,"%d %d %lf",&i,&j,&val)>0){
		if(__start == 0){
			__start =1;
			r_start = i;
			last_row = i;
		}

		if(last_row != i) {
			mchunk->row_info[mchunk->row_ptr]=last_row;
			mchunk->row_info[mchunk->row_ptr+1]=r_elements_count;
			last_row=i;
			mchunk->row_ptr +=2;
			r_elements_count=0;
		}
		

		
		if(i-r_start>= rows){
			int tmp_pos = ftell(fp);
			fseek(fp,current_file_pos - tmp_pos,SEEK_CUR);
			break;
		}

		current_file_pos = ftell(fp);
/* double the size of the vector each time you fill out the vector */
		if(mchunk->cv_ptr>=curr_vec_size-1){
			curr_vec_size *=2;
			mchunk->col_info = 
				(int *)realloc(mchunk->col_info,curr_vec_size*sizeof(int));
			mchunk->values = 
				(double *)realloc(mchunk->values,curr_vec_size*sizeof(double));

		}	
					
		r_elements_count++;
		mchunk->col_info[mchunk->cv_ptr]=j;
		mchunk->values[mchunk->cv_ptr]=val;
		mchunk->cv_ptr++;
	}

	if(feof(fp)){
		mchunk->row_info[mchunk->row_ptr]=last_row;
		mchunk->row_info[mchunk->row_ptr+1]=r_elements_count;
		mchunk->row_ptr +=2;
	}
	return mchunk;
}


void printStruct(struct matrixChunk *mchunk){

	printf("row ptr = %d , cv_ptr = %d\n",mchunk->row_ptr,mchunk->cv_ptr);
	int i =0;
	printf(" row values : \n");
	for(i=0;i<mchunk->row_ptr;i=i+2){
		printf("%d - %d ",mchunk->row_info[i],mchunk->row_info[i+1]); 
	}
		
	printf("\n j /values :\n");
	for(i=0;i<mchunk->cv_ptr;i++){
	printf("j= %d val =%lf |",mchunk->col_info[i],mchunk->values[i]);
	}
	printf("\n");
}
