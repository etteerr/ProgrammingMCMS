
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
		{
		splitsp(a, start, mid, tmp);
		splitsp(a, mid + 1, end, tmp);
		merge(a, start, end, tmp);
		}
#pragma omp taskwait
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


#define MB (1000000/sizeof(int))
#define sec 1.0e-6

double* runTest(int* a, int * tmp, unsigned int N, long seed, unsigned int repeat) {
	//data
	double *results;
	int result_iter = 0;
	results = (double*)malloc(sizeof(double)*3);

	//printf("Elements: %u\n", N);

	//gen
	int i,h;

	//printf("Input: ");
	//printArray(a, N);

	struct timeval start, end;
	double t1, t2;
	t1=t2=0;
	for(h=0; h<repeat; h++) {
		//Start seq
		//sleep(1); //let the core's chill out
		//Generate
		srand(seed);
		for(i = 0; i < N; i++)
			a[i] = N-i;//rand();

		gettimeofday(&start, 0);
		{
			splits(a, 0, N-1, tmp);
		}
		gettimeofday(&end, 0);

		if (!checkArray(a,N-1)) {
			printf("Sequential fail\n");
			return 0;
		}

		t1 += (end.tv_sec - start.tv_sec)*1e6 + (end.tv_usec - start.tv_usec);

		//Generate again
		srand(seed);
		for(i = 0; i < N; i++)
			a[i] = N-i;//rand();

		//sleep(1);
		gettimeofday(&start, 0);
		//Start paralllel
		// Define every splits call as a task
		#pragma omp parallel
		{
			#pragma omp single
			splitsp(a, 0, N-1, tmp);
		}
		gettimeofday(&end, 0);

//		printArray(a,N);

		if (!checkArray(a,N-1)) {
			printf("Parallel fail\n");
			return 0;
		}

		t2 += (end.tv_sec - start.tv_sec)*1e6 + (end.tv_usec - start.tv_usec);


		//printf("Output: ");
		//printArray(a, N);


	}
	//Average
	t1 /= (double)repeat;
	t2 /= (double)repeat;
	printf("Sequential: %.5f\nTask Parallel: %.5f (%.5f%%)\n", t1*sec, t2*sec, (t1/t2)*100.0);
	results[result_iter++] = t1;
	results[result_iter++] = t2;
	results[result_iter++] = (t1/t2)*100.0;

	return results;
}

void msort(int *vector, size_t length) {
	int * tmp;
	tmp = (int*) malloc(sizeof(int)*length);

	splits(vector, 0, length-1, tmp);
}


int testroutine(unsigned int * elem, unsigned short n, unsigned int repeat, long seed) {
	//Declare
	int *a, *tmp, i,j;
	a=tmp=0;
	int N;


	double **results;
	results = malloc(sizeof(double*)*n);

	for(i=0; i<n; i++) {
		N = elem[i];
		//Notify
		printf("Array size: %i (%.5f MB)\n", N, (double)N/(double)MB);
		//alloc
		int res = posix_memalign((void**)&tmp, 16, sizeof(int)*N);
		res |= posix_memalign((void**)&a, 16, sizeof(int)*N);
		if (res)
			return res;
		//Run
		results[i] = runTest(a,tmp,elem[i], seed, repeat);

		//free
		free(a);
		free(tmp);

		if (results[i]==0) {
			for(j=i; j>=0; j--)
				free(results[j]);
			free(results);
			printf("Invalid result in test\n");
			return 1;
		}
	}


	//Out results
	FILE* f = fopen("Results.csv", "w");

	for(i=0; i<n; i++) {
		for(j=0; j<repeat*3; j++) {
			fprintf(f, "%.5f; ", results[i][j]);
		}
		fprintf(f,"\n");
	}

	fclose(f);

	return 0;
}


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

  int c;
  int seed = 42;
  long length = 1e4;
  Order order = ASCENDING;
  int *vector;

  /* Read command-line options. */
  while((c = getopt(argc, argv, "adrgl:s:")) != -1) {
    switch(c) {
      case 'a':
        order = ASCENDING;
        break;
      case 'd':
        order = DESCENDING;
        break;
      case 'r':
        order = RANDOM;
        break;
      case 'l':
        length = atol(optarg);
        break;
      case 'g':
        debug = 1;
        break;
      case 's':
        seed = atoi(optarg);
	break;
      case '?':
        if(optopt == 'l' || optopt == 's') {
          fprintf(stderr, "Option -%c requires an argument.\n", optopt);
        }
        else if(isprint(optopt)) {
          fprintf(stderr, "Unknown option '-%c'.\n", optopt);
        }
        else {
          fprintf(stderr, "Unknown option character '\\x%x'.\n", optopt);
        }
        return -1;
      default:
        return -1;
      }
  }

  /* Seed such that we can always reproduce the same random vector */
  srand(seed);

  /* Allocate vector. */
  vector = (int*)malloc(length*sizeof(int));
  if(vector == NULL) {
    fprintf(stderr, "Malloc failed...\n");
    return -1;
  }

  /* Fill vector. */
  switch(order){
    case ASCENDING:
      for(long i = 0; i < length; i++) {
        vector[i] = (int)i;
      }
      break;
    case DESCENDING:
      for(long i = 0; i < length; i++) {
        vector[i] = (int)(length - i);
      } 
      break;
    case RANDOM:
      for(long i = 0; i < length; i++) {
        vector[i] = rand();
      }
      break;
  }

  if(debug) {
    print_v(vector, length);
  }


  if(run_Testroutine) {
		int size = 5;
		unsigned int arr[size];
		int i;
		for(i=0; i<size; i++) {
			arr[i] = i*i*i*MB+3;
		}
		return testroutine(arr, size, 2, seed);
  }else
	  /* Sort */
	    msort(vector, length);

  if(debug) {
    print_v(vector, length);
  }

  return 0;
}

