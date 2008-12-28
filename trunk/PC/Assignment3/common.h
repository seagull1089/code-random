#include<stdio.h>
#include<stdlib.h>
#include<string.h>


int min(int,int);
int max(int,int);
char * read_data(char *file_name);
void print_matrix(int **A,int N,int M);
void output_lcs(int **A,char *string1,char *string2);
