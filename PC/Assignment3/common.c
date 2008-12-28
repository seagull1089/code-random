#include"common.h"

int __MAXLEN=100;
int max(int a,int b){
	return ((a>b)?a:b);
}

int min(int a,int b){
	return ((a<b)?a:b);
}
void output_lcs(int **A,char *string1,char *string2){
	int i= strlen(string2);
	int j= strlen(string1);
	int lcs_len = max(i,j);
	int  __follow_i= 0;
	if(lcs_len == i){
		__follow_i= 1;
	}
	char *tmp_str = (char *)malloc((lcs_len+1)*sizeof(char));
	int l =0;
	for(l=0;l<lcs_len;l++){
		tmp_str[l] = '_';
	}
	tmp_str[l]='\0';	
	int k =lcs_len-1;
	while(i>0 && j>0){
		if(string2[i-1] == string1[j-1]){
			tmp_str[k] = string1[j-1];
			//printf("%c",string1[j-1]);	
			i--;
			j--;
			k--;
		}else{
			if(A[i-1][j] > A[i][j-1]){
				i--;
				if(__follow_i){
				k--;
				}	
			}else{
				j--;
				if(!__follow_i){
				k--;
				}
			}
		}
	}
      //  printf("\nThe Large sequence is : %s\n",string1);	
       // printf("The Short sequence is : %s\n",string2);	
	printf("The outputString is   : %s\n",tmp_str);
	free(tmp_str);
}


void print_matrix(int **A,int M,int N){
	int i=0,j=0;
	for(i=0;i<M;i++){
		for(j=0;j<N;j++){
			printf("%d ",A[i][j]);	
		}
		printf("\n");
	}
}


char * read_data(char *file_name){

	char dummy[__MAXLEN];
	FILE *fp = fopen(file_name,"r");
	fseek(fp,0,SEEK_END);
 	int length = ftell(fp) +1;
	char *string = NULL;
	string = (char *)malloc(length*sizeof(char));
	memset(string,length,'\0');
 	//printf("Length of the file = %d \n",length);
	rewind(fp);
	int i =0;
/*
	for(i=0;i<strlen(string);i++){
		printf("DEBUG1 i = %d , string[i] = %c \n",i,string[i]);
	}
*/
	memset(dummy,__MAXLEN,'\0');
	while(fgets(dummy,__MAXLEN-1,fp)!=NULL){
		if(dummy[strlen(dummy)-1]=='\n'){
			dummy[strlen(dummy)-1]='\0';
		}
		strcat(string,dummy);
		memset(dummy,__MAXLEN,'\0');	
 	}
/*
	for(i=0;i<strlen(string);i++){
		printf("DEBUG i = %d , string[i] = %c \n",i,string[i]);
	}
*/
	string[strlen(string)]='\0';
	fclose(fp);
	return string;
}
