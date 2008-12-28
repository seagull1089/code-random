#include<stdio.h>
#include<pthread.h>
#include<stdlib.h>
#include<time.h>

pthread_t      *tids;
float          *array, pivot;
int             N;
int             t_count;
int             rand_index, rank;

int            *local_ranks;
int            *partition_indices;

extern          get_random(int, int, int, int);
void           *thread_function(void *);
int             get_my_partition(int *, int);
int             get_local_ranks_sum();
float           get_new_pivot();
void            set_partitions(int);
void            set_initial_partitions();
int             partition(float A[], int p, int r, float pivot);
int             has_converged(int rank, int required_rank);
int             debug = 0;


int 
main(int argc, char **argv)
{

	if (argc < 3) {
		printf("Usage : progname filename thread_count [rank]\n");
		exit(1);

	}
	char           *filename = argv[1];
	t_count = atoi(argv[2]);

	int             i = 0;
	int            *thread_ids;
	float           var;
	int             local_ranks_sum, median = -1;
	int             loop = 0;
	clock_t         clk_start, clk_end;
	FILE           *fp = fopen(filename, "r");
	fscanf(fp, "%d", &N);

	if (argc == 4) {
		rank = atoi(argv[3]);
	} else {
		rank = N / 2;
	}
	if (debug)
		printf(" total numbers = %d total thrreads = %d median rank = %d\n", N, t_count, rank);

	array = (float *) malloc(N * sizeof(float));
	while ((fscanf(fp, "%f ", &var) > 0) && i < N) {
		array[i] = var;
		i++;
	}

	fclose(fp);

	clk_start = clock();

	thread_ids = (int *) malloc(t_count * sizeof(int));
	tids = (pthread_t *) malloc(t_count * sizeof(pthread_t));
	local_ranks = (int *) malloc(t_count * sizeof(int));
	partition_indices = (int *) malloc(2 * t_count * sizeof(int));

	for (i = 0; i < t_count; i++) {
		local_ranks[i] = -1;
		thread_ids[i] = i;

	}
	set_initial_partitions();
	int             rand_index = get_random(0, N - 1, N, getpid());
	pivot = array[rand_index];

	if (debug)
		printf("initial random value = %d\
		       	pivot value = %f\n",
		       rand_index, pivot);

	while (1) {

		for (i = 0; i < t_count; i++) {
			pthread_create(&tids[i], NULL,
				       thread_function, &thread_ids[i]);
		}

		for (i = 0; i < t_count; i++) {
			(void) pthread_join(tids[i], NULL);
		}
		local_ranks_sum = get_local_ranks_sum();

		//printf("INFO  local_ranks sum = %d pivot value = %f \n", local_ranks_sum, pivot);

		if ((local_ranks_sum == rank) || has_converged(local_ranks_sum, rank)) {
			median = pivot;
			break;
		} else {
			if (local_ranks_sum < rank) {
				/* search on the right side */
				set_partitions(0);
			} else if (local_ranks_sum > rank) {
				/* search on left side */
				set_partitions(1);
			}
			pivot = get_new_pivot();

		}
		loop++;

	}
	clk_end = clock();
	printf("Total Execution time = %f \n", (double) (clk_end - clk_start) / CLOCKS_PER_SEC);
	printf("median = %d \n", median);

	free(array);
	free(thread_ids);
	free(tids);
	free(partition_indices);
	free(local_ranks);

}
int 
has_converged(int rank, int required_rank)
{
	static int      ranks[10] = {0};
	static int      maxrank = 0, pointer = -1, iter = 0;;
	if (rank <= required_rank && rank >= maxrank) {
		maxrank = rank;
		ranks[(++pointer) % 10] = maxrank;
		iter++;
	}
	int             i = 0, converged = 0;
	if (iter > 10) {
		converged = 1;
		for (i = 0; i < 10; i++) {
			if (ranks[i] != maxrank) {
				converged = 0;
				break;
			}
		}
	}
	return converged;
}


float 
get_new_pivot()
{
	float           __tmp = -1;
	int             part_num;
	long            __index;
	static int      count = 0;
	while (__tmp < 0) {
		/*
		 * slect a partition and then chose random element in that
		 * partition
		 */
		part_num = get_random(0, t_count - 1, t_count, ++count);

		if (partition_indices[2 * part_num] == partition_indices[2 * part_num + 1]) {
			__index = partition_indices[2 * part_num];
		} else {
			if (debug)
				printf(" indices : %d %d ", partition_indices[2 * part_num], partition_indices[2 * part_num + 1]);
			__index = get_random(partition_indices[2 * part_num],
					     partition_indices[2 * part_num + 1], N, (++count) * getpid());
		}
		__tmp = array[__index];

	}
	if (debug)
		printf("DEBUG new pivot value = %f index = %d part_num = %d  pivot = %f \n", __tmp, __index, part_num, array[__index]);
	return __tmp;
}
void 
set_partitions(int left)
{
	int             i = 0;
	int             __tmp;
	int             partition_size = N / t_count;
	for (i = 0; i < t_count; i++) {
		if (left == 1) {
			__tmp = i * partition_size + local_ranks[i] - 1;
			if (__tmp > partition_indices[2 * i]) {
				partition_indices[2 * i + 1] = __tmp;
			}
			//partition_indices[2 * i + 1] = i * partition_size + local_ranks[i] - 1;
		} else {
			__tmp = i * partition_size + local_ranks[i];
			if (__tmp < partition_indices[2 * i + 1]) {
				partition_indices[2 * i] = __tmp;;
			}
			/*
			 * else{ partition_indices[2*i] =
			 * partition_indices[2*i+1]; }
			 */
		}

	}
}
void 
set_initial_partitions()
{
	int             i;
	int             partition_size = N / t_count;
	for (i = 0; i < t_count; i++) {
		partition_indices[2 * i] = i * partition_size;
		if (i == t_count - 1) {
			partition_indices[2 * i + 1] = N - 1;
		} else {
			partition_indices[2 * i + 1] = (i + 1) * partition_size - 1;
		}
	}

}
int 
get_local_ranks_sum()
{
	int             i = 0, sum = 0;
	for (i = 0; i < t_count; i++) {
		sum += local_ranks[i];
	}
	return sum + 1;
}

int 
get_my_partition(int *ind, int id)
{

	int             partition_size = N / t_count;
	ind[0] = id * partition_size;
	if (id == t_count - 1) {
		ind[1] = N - 1;
	} else {
		ind[1] = (id + 1) * partition_size - 1;
	}
	return 0;

}
void           *
thread_function(void *arg)
{

	int            *myid = (int *) arg;
	int             local_rank;

	int             lower = partition_indices[2 * (*myid)];
	int             higher = partition_indices[2 * (*myid) + 1];
	int             partition_size = N / t_count;


	pthread_t       pid = pthread_self();
	local_rank = partition(array, lower,
			    higher, pivot) - ((*myid) * partition_size) + 1;
	local_ranks[*myid] = local_rank;
	if (debug)
		printf("DEBUG , id = %d  pthread_t pid = %d \
		       	lower = %d higher = %d local_rank =%d\n",
		       *myid, pid, lower, higher, local_ranks[*myid]);

	return (void *) 0;
}



int 
partition(float A[], int p, int r, float pivot)
{
	/* should you changed the intial value of i? */
	int             j, i = p - 1;
	float           tmp;
	/* loop runs from p to r */
	for (j = p; j <= r; j++) {
		/* modified the condition of < to <= */
		if (A[j] < pivot) {
			i = i + 1;
			tmp = A[i];
			A[i] = A[j];
			A[j] = tmp;
		}
	}
	/* just return the last position */
	return i;
}
