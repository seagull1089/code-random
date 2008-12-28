#include"utils.h"

int * readfromFile(FILE *fp, int N){

	int *data = malloc(N*sizeof(int));
	int var=0;
	int cntr = 0;
	while((cntr <N) && (fscanf(fp,"%d ",&var)>0)){
		data[cntr]=var;
		cntr++;
        }
	return data;
}

int * get_partition_indices(int A[],int length,int bucket_limits[],int b_length){

	int *p_indices = malloc(b_length*sizeof(int));
	int j=0,i = 0;
	for(i=0;i<b_length;i++){
		p_indices[i]=-1;
	}
	
	for(i=0;i<length;i++){
		
		if(A[i]<bucket_limits[j]){
			continue;
		}else if(A[i]>bucket_limits[j]){
			p_indices[j]=i;
			while(j<(b_length) && A[i]>bucket_limits[j]){
			//	p_indices[j]=i;
				j++;
			}
			p_indices[j-1]=i;
			if(j==b_length){
				break;
			}
		}	
		
	}
	if(i==length){
		p_indices[j]=length;
	}
	return p_indices;

}
		/*
		if(A[i]<bucket_limits[0]){
			if(j!=0){
				j=0;
				p_indices[j]=i;
			}
		}else if(A[i]){

		}else if(A[i]>bucket_limits[b_length-1]){
			if(j!=b_length-1){
				j=b_length-1;
				p_indices[j]=i;
			}

		}
		*/
