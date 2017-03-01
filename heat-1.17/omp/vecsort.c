/*
 * vecsort.c
 *
 *  Created on: Mar 1, 2017
 *      Author: erwin
 */



#include <stdlib.h>
#include <unistd.h>
#include <getopt.h>
#include <stdio.h>
#include <ctype.h>
#include <omp.h>

/* Ordering of the vector */
typedef enum Ordering {ASCENDING, DESCENDING, RANDOM} Order;

int debug = 0;

/* Sort vector v of l elements using mergesort */
/*
 * main.cpp
 *
 *  Created on: Feb 25, 2017
 *      Author: erwin
 */

#define run_Testroutine 1

#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>

void printArray(int a[], int N)
{
	for (size_t i = 0; i < N; i++)
	{
		printf("%d,", a[i]);
	}
	printf("\n");
}

void printArrayft(int a[], int from, int to)
{
	for (size_t i = from; i < to; i++)
	{
		printf("%d,", a[i]);
	}
	printf("\n");
}

char checkArray(int a[], int N)
{
	unsigned int i;
	for (i = 1; i < N; i++)
	{
		if (!(a[i-1] <= a[i])) {
			//printArray(a,N);
			return 0;
		}
	}

	return 1;
}

void merge(int a[], int start, int end, int tmp[])
{
	int mid =  (start + end) / 2;
	int i, j;
	int leftCounter = start;
	int rightCounter = mid + 1;

//	printArrayft(a, start, end);

	for (i = start; (leftCounter<=mid) && (rightCounter<=end); i++)
	{
		if (a[leftCounter] < a[rightCounter])
		{
			tmp[i]= a[leftCounter++];
		}
		else
		{
			tmp[i] = a[rightCounter++];
		}
	}

	//merge leftovers in the left array
	for ( ; leftCounter<=mid; i++, leftCounter++)
	{
		tmp[i] = a[leftCounter];
	}

	//merge leftovers in the right array
	for ( ; rightCounter<=end; i++, rightCounter++)
	{
		tmp[i] = a[rightCounter];
	}

	//update the main array
	for ( i = start, j = start; i <= end; i++, j++) //Because temp always starts from 0 till local end, while the main array might have separate sections needed to be updated.
	{
		a[i] = tmp[j];
	}

//	printArrayft(a, start, end);
}

void splitsp(int a[], int start, int end, int tmp[])
{
	if (start < end)
	{
		int mid = (start + end) / 2;
#pragma omp task
		splitsp(a, start, mid, tmp);
#pragma omp task
		splitsp(a, mid + 1, end, tmp);
#pragma omp taskwait
		merge(a, start, end, tmp);
	}
}


void splits(int a[], int start, int end, int tmp[])
{
	if (start < end)
	{
		int mid = (start + end) / 2;
		splits(a, start, mid, tmp);
		splits(a, mid + 1, end, tmp);
		merge(a, start, end, tmp);
	}
}

void vecsort(int **a, int *lengths, unsigned long size) {
	int i,j;

	for(i=0; i<size; i++) {
		int * tmp;
		if (posix_memalign(&tmp, 16, sizeof(int)*lengths[i])) {
			fprintf(stderr, "Error while allocating tmp memory.");
		}
		splits(a[i], 0, lengths[i]-1, tmp);
		free(tmp);
	}
}

void vecsort_loop(int **a, int *lengths, unsigned long size) {
	int i,j;

#pragma omp parallel for schedule(guided) shared(a,lengths, size) private(i)
	for(i=0; i<size; i++) {
		int * tmp;
		if (posix_memalign(&tmp, 16, sizeof(int)*lengths[i])) {
			fprintf(stderr, "Error while allocating tmp memory.");
		}
		splits(a[i], 0, lengths[i]-1, tmp);
		free(tmp);
	}
}

void vecsort_task(int **a, int *lengths, unsigned long size) {
	int i,j;
#pragma omp parallel sections
	{
#pragma omp single
		{
			for(i=0; i<size; i++) {
				int * tmp;
				if (posix_memalign(&tmp, 16, sizeof(int)*lengths[i])) {
					fprintf(stderr, "Error while allocating tmp memory.");
				}
				splitsp(a[i], 0, lengths[i]-1, tmp);
				free(tmp);
			}
		}
	}
}

void vecsort_fullpar(int **a, int *lengths, unsigned long size) {
	int i,j;
#pragma omp parallel for schedule(guided) shared(a,lengths, size) private(i)
	for(i=0; i<size; i++) {
		int * tmp;
		if (posix_memalign(&tmp, 16, sizeof(int)*lengths[i])) {
			fprintf(stderr, "Error while allocating tmp memory.");
		}
		splitsp(a[i], 0, lengths[i]-1, tmp);
		free(tmp);
	}
}

#define MB (1000000/sizeof(int))
#define sec 1.0e-6

int ** generateVec(int * lengths, unsigned long N, long seed) {
	srand(seed);

	int ** vecvec;
	vecvec = malloc(sizeof(int*)*N);

	unsigned long i,j;
	for(i=0; i<N; i++) {
		if (posix_memalign(&vecvec[i], 16, sizeof(int)*lengths[i])) {
			exit(666);
		}

		for(j=0; j<lengths[i]; j++) {
			vecvec[i][j] = rand();
		}
	}

	return vecvec;
}

double* runTest(int *lengths, unsigned long N, long seed, unsigned int repeat) {
	//data
	double *results;
	int result_iter = 0;
	results = (double*)malloc(sizeof(double)*4);

	//printf("Elements: %u\n", N);

	//gen
	int i,h;

	//printf("Input: ");
	//printArray(a, N);

	struct timeval start, end;
	double t1, t2, t3, t4;
	int **a;
	t1=t2=t3=t4=0;

	for(h=0; h<repeat; h++) {
		a = generateVec(lengths, N, seed);
		gettimeofday(&start, 0);
		{
			//seq
			vecsort(a,lengths,N);
		}
		gettimeofday(&end, 0);

		t1 += (end.tv_sec - start.tv_sec)*1e6 + (end.tv_usec - start.tv_usec);

		a = generateVec(lengths, N, seed);
		gettimeofday(&start, 0);
		{
			//Parallel loop
			vecsort_loop(a,lengths,N);
		}
		gettimeofday(&end, 0);

		t2 += (end.tv_sec - start.tv_sec)*1e6 + (end.tv_usec - start.tv_usec);

		a = generateVec(lengths, N, seed);
		gettimeofday(&start, 0);
		{
			//Parallel task
			vecsort_task(a,lengths,N);
		}
		gettimeofday(&end, 0);

		t3 += (end.tv_sec - start.tv_sec)*1e6 + (end.tv_usec - start.tv_usec);

		a = generateVec(lengths, N, seed);
		gettimeofday(&start, 0);
		{
			//Parallel task + loop
			vecsort_fullpar(a,lengths,N);
		}
		gettimeofday(&end, 0);

		t4 += (end.tv_sec - start.tv_sec)*1e6 + (end.tv_usec - start.tv_usec);


	}
	//Average
	t1 /= (double)repeat;
	t2 /= (double)repeat;
	t3 /= (double)repeat;
	t4 /= (double)repeat;
	printf("Sequential: %.5f\n", t1*sec);
	printf("looppar   : %.5f\n", t2*sec);
	printf("taskpar   : %.5f\n", t3*sec);
	printf("allpar    : %.5f\n", t4*sec);
	results[result_iter++] = t1;
	results[result_iter++] = t2;
	results[result_iter++] = t3;
	results[result_iter++] = t4;

	return results;
}

void msort(int *vector, size_t length) {
	int * tmp;
	tmp = (int*) malloc(sizeof(int)*length);

	splits(vector, 0, length-1, tmp);
}


//int testroutine(unsigned int * elem, unsigned short n, unsigned int repeat, long seed) {
//	//Declare
//	int *a, *tmp, i,j;
//	a=tmp=0;
//	int N;
//
//
//	double **results;
//	results = malloc(sizeof(double*)*n);
//
//	for(i=0; i<n; i++) {
//		N = elem[i];
//		//Notify
//		printf("Array size: %i (%.5f MB)\n", N, (double)N/(double)MB);
//		//alloc
//		int res = posix_memalign((void**)&tmp, 16, sizeof(int)*N);
//		res |= posix_memalign((void**)&a, 16, sizeof(int)*N);
//		if (res)
//			return res;
//		//Run
//		results[i] = runTest(a,tmp,elem[i], seed, repeat);
//
//		//free
//		free(a);
//		free(tmp);
//
//		if (results[i]==0) {
//			for(j=i; j>=0; j--)
//				free(results[j]);
//			free(results);
//			printf("Invalid result in test\n");
//			return 1;
//		}
//	}
//
//
//	//Out results
//	FILE* f = fopen("Results.csv", "w");
//
//	for(i=0; i<n; i++) {
//		for(j=0; j<repeat*3; j++) {
//			fprintf(f, "%.5f; ", results[i][j]);
//		}
//		fprintf(f,"\n");
//	}
//
//	fclose(f);
//
//	return 0;
//}


void print_v(int *v, long l) {
  printf("\n");
  for(long i = 0; i < l; i++) {
    if(i != 0 && (i % 10 == 0)) {
      printf("\n");
    }
    printf("%d ", v[i]);
  }
  printf("\n");
}

int main(int argc, char **argv) {

	if(argc < 3) {
		printf("Invalid usage: vecsort seed length maxveclength\n");
		return 1;
	}

	int * lengths;
	unsigned long long size = 0;
	unsigned long i, N = atoll(argv[2]);
	long seed = atoll(argv[1]);
	int max = atoi(argv[3]);

	srand(seed);

	lengths = malloc(sizeof(int)*N);
	for(i=0; i<N; i++) {
		lengths[i] = ((double)rand()/(double)RAND_MAX)*(double)max;
		size += lengths[i];
	}

	if (size/MB>4000) {
		printf("Expected size %llu MB, terminating...\n", size/MB);
		exit(2);
	}
	double * res = runTest(lengths, N, seed, 3);

	FILE *fp;
//	fp = fopen("/var/scratch/ppp1620/vecsort.csv", "a");
	fp = fopen("vecsort.csv", "a");

	if (ftell(fp)==0) {
		fprintf(fp, "seed; length; seq; looppar; taskpar; allpar;\n");
	}

	fprintf(fp, "%l;%lu;%.5f;%.5f;%.5f;%.5f;\n",seed, N, res[0],res[1],res[2],res[3]);

	fflush(fp);
	fclose(fp);
	free(res);

  return 0;
}

