/*
 * pipeline.c
 *
 *  Created on: Mar 12, 2017
 *      Author: erwin
 */
#include <pthread.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/time.h>

#define printarr 0

#define null 0
#define kilobyte 1000
#define megabyte 1000*kilobyte
#define gigabyte 1000*megabyte

#define endsig 0xFAFAFAFAFAFAFAFA //mask is 64 bits while rand numbers are 32 bits and this is unique

enum order {
	ascending,
	decending,
	orandom
};

//thread data (shared)
struct thread_data {
	volatile long int number; //8
	long int rank; //8
	volatile int flag; //4 mail postal service flag
	pthread_attr_t * pAttr;
	long int * arr;
	char pad[64-36]; //pad to page
};

void * thread(struct thread_data * data);

void printArray(long int *arr, unsigned long size);
char checkArray(long int *arr, unsigned long size);

int main(int nargs, char ** args) {
	register long int i = 1;
	enum order o = orandom;
	unsigned long nElem = 400;
	unsigned int seed = 0;

	//Process cmdline opts
	char *cmd;
	while(i < nargs) {
		cmd = args[i++];
		if (strcmp(cmd , "-a")==0)
			o = ascending;

		if (strcmp(cmd , "-d")==0)
			o = decending;

		if (strcmp(cmd , "-r")==0)
			o = orandom;

		if (strcmp(cmd , "-n")==0)
			nElem = strtol(args[i++],NULL, 10);

		if (strcmp(cmd , "-s")==0)
			nElem = atoi(args[i++]);
	}

	//set seed
	srand(seed);

	//Check
	if (nElem == 0) {
		printf("nElem must be greater than 0 or requires valid input.\n");
		return 1;
	}

	//Create elements
	long int * elems;
	i = posix_memalign((void**)&elems, 128/8, sizeof(long int) * nElem);
	if (i) {
		printf("Allocation failed with code %li.\n", i);
		return 2;
	}

	//Fill elements
	if (o==ascending)
		for(i=0; i<nElem;i++) {
			elems[i] = i;
		}

	if (o==decending)
		for(i=nElem; i>0;i--) {
			elems[i] = i;
		}

	if (o==orandom)
		for(i=0; i<nElem;i++) {
			elems[i] = rand();
		}

	//Init pthreads attr
	pthread_attr_t pAttr;
	pthread_attr_init(&pAttr);
	pthread_attr_setstacksize(&pAttr, 32*kilobyte); //default
	pthread_attr_setdetachstate(&pAttr, PTHREAD_CREATE_DETACHED);

	//Allocate structs for threads
	struct thread_data * tds;
	i = posix_memalign((void**)&tds, 128/8, sizeof(struct thread_data)*(nElem+1));
	if (i) {
		printf("Allocation of thread_data failed with code %li.\n", i);
		return 3;
	}

	//init mem
	for(i=0; i<(nElem+1);i++) {
		tds[i].flag = 0;
		tds[i].number = 0;
		tds[i].rank = i;
		tds[i].pAttr = &pAttr;
		tds[i].arr = elems;
	}

	//start posix
	pthread_t first_thread;
	pthread_create(&first_thread, &pAttr, thread, (void*) &tds[0]);

	//print array
	if (printarr) printArray(elems, nElem);

	//Get wallclock things
	struct timeval start, end;
	//feed data
	gettimeofday(&start, 0);
	for(i=0; i<(nElem); i++) {

		while(1) {
			if (tds[0].flag==0) {
				tds[0].number = elems[i];
				tds[0].flag = 1;
				elems[i] = 0;
				break;
			}
			sched_yield();
		}
	}


	//send end
	while(1) {
		if (tds[0].flag==0) {
			tds[0].number = endsig;
			tds[0].flag = 1;
			break;
		}
		sched_yield();
	}

	while(elems[nElem-1] == 0) sched_yield();
	//All is done
	gettimeofday(&end, 0);

	if (printarr) printArray(elems, nElem);

	if (checkArray(elems, nElem)) {
		printf("Incorrect result\n");
	}

	//print time
	double time = (end.tv_sec - start.tv_sec) + (end.tv_usec - start.tv_usec)/1.0e6;
	printf("Time: %.9f\n", time);


	//store results
	FILE *f;
	f = fopen("/var/scratch/ppp1620/pipeline.csv", "a");

	if (!f)
		return 1;

	if (ftell(f)==0)
		fprintf(f, "Size; order (asc=0 desc=1 rand=2); time;\n");

	fprintf(f, "%lu; %i; %f;\n", nElem, (int)o, time);
	fclose(f);


	return 0;

}


//To pass on send function
/**
 * passon_send
 * 	Sends num to the next thread in line
 * 	Assumes thread exists (and its thread_data)
 */
static inline void passon_send(long int num, struct thread_data * data) {
	while(1) {
		if ((data+1)->flag==0) {
			(data+1)->number = num;
			(data+1)->flag = 1;
			break;
		}
		sched_yield();
	}
}

void * thread(struct thread_data * data) {
	register long int my_number;
	pthread_t child = 0;
	register long int pass_on_buffer;
	register long int numbuf;
	register int state = 0;

	//initial recv
	while(1) {
		if (data->flag) {
			my_number = data->number;
			data->flag = 0;
			break;
		}
		sched_yield();
	}

	//If we are the collector
	if (my_number == endsig) {
		register long unsigned int i = 0;
		while(1) {
			if (data->flag) {
				//accept data
				numbuf = data->number;
				data->flag = 0;
				__sync_synchronize();

				if (numbuf==endsig)
					break;

				//output data
				data->arr[i++] = numbuf;

			}
			sched_yield();
		}

		return 0;
	}

	//comparator / endsig
	while(1) {
		if (data->flag) {
			//get number & reset flag
			numbuf = data->number;
			data->flag = 0;

			//Check for sig
			if (numbuf == endsig) { //sig
				if (state == 0) {
					if (!child)
						pthread_create(&child, data->pAttr, thread, (void*) (data+1));
					passon_send(endsig, data);
					passon_send(my_number, data);
					state++;
				}else {
					passon_send(endsig, data);
					return 0;
				}
			}else { //No sig
				//Comppare if state is not set, other wise just forward
				if (!state && my_number < numbuf) {
					pass_on_buffer = my_number;
					my_number = numbuf;
				}else
					pass_on_buffer = numbuf;

				//Create child
				if (!child)
					pthread_create(&child, data->pAttr, thread, (void*) (data+1));

				//Wait for target buffer to be free'ed and pass on number
				passon_send(pass_on_buffer, data);

			}

			//Memory barrier
			__sync_synchronize();
		}
		sched_yield();
	}
}

void printArray(long int *arr, unsigned long size) {
	register unsigned long i;
	for(i=0 ; i<size; i++)
			printf("%i ", arr[i]);
		printf("\n");
}
char checkArray(long int *arr, unsigned long size) {
	register unsigned long i;
	char fail = 0;

	for(i=1 ; i<size; i++) {
		if (arr[i-1] > arr[i]) {
			fail = 1;
			break;
		}
	}

	return fail;
}


