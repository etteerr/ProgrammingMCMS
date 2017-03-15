/*
 ============================================================================
 Name        : histo.c
 Author      : Erwin Diepgrond
 Version     :
 Copyright   : 
 Description : Histo bin counting algorithm
 ============================================================================
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <pthread.h>
#include <semaphore.h>
#include <sys/time.h>
#include <unistd.h>

void die(char * msg) {
	printf("%s\n", msg);
	exit(3);
}

static void readpgm_uchar(const char *fname,
                          size_t height, size_t width, unsigned char *data)
{
    char format[2];
    FILE *f;
    unsigned imgw, imgh, maxv, v;
    size_t i;

    printf("Reading PGM data from %s...\n", fname);

    if (!(f = fopen(fname, "r"))) die("fopen");

    fscanf(f, "%2s", format);
    if (format[0] != 'P' || format[1] != '2') die("only ASCII PGM input is supported");

    if (fscanf(f, "%u", &imgw) != 1 ||
        fscanf(f, "%u", &imgh) != 1 ||
        fscanf(f, "%u", &maxv) != 1) die("invalid input");

    if (imgw != width || imgh != height) {
        fprintf(stderr, "input data size (%ux%u) does not match cylinder size (%zux%zu)\n",
                imgw, imgh, width, height);
        die("invalid input");
    }

    for (i = 0; i < width * height; ++i)
    {
        if (fscanf(f, "%u", &v) != 1) die("invalid data");
        v %= 256;
        data[i] = v;//dmin + (double)v * (dmax - dmin) / maxv;
    }

    fclose(f);
}

void printhelp() {
	printf("histo [options] -f [file]\n"
			"Options:\n"
			"\t-w width\n"
			"\t-h height\n"
			"\t--help prints this message.\n");
}

void parseInput(int nargs, char **args, unsigned char ** data, unsigned long * length, int * nThreads) {
	//Check input pointers
	if (!(length && data)) {
		printf("length or data is a null pointer.\n");
		exit(1);
	}

	//Set defaults
	char file_default[128] = "plasma_100x150.pgm";
	char * file = 0;
	unsigned int width = 100;
	unsigned int height = 150;
	*nThreads = 2;

	//Get inputs
	register int i = 1;
	while(i < nargs) {

		char *cmd = args[i++];

		if (!strcmp(cmd, "-f"))
			file = args[i++];

		if (!strcmp(cmd, "-h"))
			height = atoi(args[i++]);

		if (!strcmp(cmd, "-w"))
			width = atoi(args[i++]);

		if (!strcmp(cmd, "--help"))
			printhelp();

		if (!strcmp(cmd, "-t"))
			*nThreads = atoi(args[i++]);
	}
	if (file==0) {
		file = file_default;
	}

	//Set length
	*length = width * height;

	//print arguments
	printf("hist, starting with:\n"
			"\t-h %i\n"
			"\t-w %i\n"
			"\t file: %s\n"
			"\t threads: %i\n",
			height,
			width,
			file,
			*nThreads);

	//Alloc
	int res = posix_memalign((void **)data, 128/8, sizeof(char)*width*height);
	if (res) {printf("Allocation error (%i)\n", res); exit(5);}

	//load
	readpgm_uchar(file, height, width, *data);

}

typedef enum {
	init,
	work,
	saving,
	done
}thread_flag;

typedef struct {
	pthread_t pid;
	unsigned int id;
	unsigned long from;
	unsigned long amount;
	unsigned int bin[256];
	unsigned int * shared_bin;
	pthread_mutex_t * lock;
	unsigned char * data;
	volatile thread_flag *flag;
	sem_t * semaphore;
}thread_data_mutex_sync;

void * threader(thread_data_mutex_sync * data) {
	register unsigned long length, from, i;
	register unsigned char * img;
	img = data->data;

	//Init registers
	from = data->from;
	length = data->amount;

	//run
	*data->flag = work;

	for(i=from; i<from+length; i++) {
		//printf("tf: %i %i\n", img[i], data->bin[img[i]]); fflush(stdout);
		data->bin[img[i]]++;
	}

	//dump
	*data->flag = saving;
	pthread_mutex_lock(data->lock);
	for(i=0; i<256;i++) {
		data->shared_bin[i] += data->bin[i];
	}
	pthread_mutex_unlock(data->lock);

	//done
	*data->flag = done;

	return 0;

}

double run_sync_mutex(unsigned int *bin, unsigned char* data, unsigned long length, int nThreads) {
	register int i;

	//init pthreads
	thread_data_mutex_sync * tds;
	pthread_attr_t tAttr;
	pthread_mutex_t mLock;
	pthread_mutexattr_t mAttr;
	volatile thread_flag flag = init;
	struct timeval start, stop;

	//Init memory
	register int res = posix_memalign(&tds, 128/8, sizeof(thread_data_mutex_sync) * nThreads);
	if (res) die("Allocation error for thread_data_mutex_sync");
	memset(tds, 0, sizeof(thread_data_mutex_sync) * nThreads);
	pthread_attr_init(&tAttr);
	pthread_attr_setdetachstate(&tAttr, PTHREAD_CREATE_JOINABLE);
	pthread_mutexattr_init(&mAttr);
	pthread_mutex_init(&mLock, &mAttr);

	//determine stride
	unsigned long stride, left;
	stride = length / nThreads;
	left = length % nThreads;


	//Start time, because after this forloop, the threads have already done some work...
	gettimeofday(&start, 0);

	//Create threads
	for(i=nThreads-1; i>=0; i--) {
		tds[i].id = i;
		tds[i].flag = &flag;
		tds[i].data = data;
		tds[i].from = i*stride;
		tds[i].shared_bin = bin;
		tds[i].lock = &mLock;

		//The last thread gets the extra
		//Note: Inverse order, so last thread is the first :D
		if (i==nThreads-1)
			tds[i].amount = stride + left;
		else
			tds[i].amount = stride;

		//Create & run
		pthread_create(&tds[i].pid, &tAttr, threader, &tds[i]);
	}


	//join
	for(i=0; i<nThreads; i++)
		pthread_join(tds[i].pid, 0);

	//done
	gettimeofday(&stop, 0);

	//free
	free(tds);
	pthread_attr_destroy(&tAttr);
	pthread_mutex_destroy(&mLock);
	pthread_mutexattr_destroy(&mAttr);

	//return time
	return (double)(stop.tv_sec-start.tv_sec) + (stop.tv_usec-start.tv_usec)/1.0e6;
}

void * threader_fullmutex(thread_data_mutex_sync * data) {
	register unsigned long length, from, i;
	register unsigned char * img;
	img = data->data;

	//Init registers
	from = data->from;
	length = data->amount;

	//run
	*data->flag = work;

	for(i=from; i<from+length; i++) {
		pthread_mutex_lock(data->lock);
		data->shared_bin[img[i]]++;
		pthread_mutex_unlock(data->lock);
	}

	//done
	*data->flag = done;

	return 0;

}

double run_full_mutex(unsigned int *bin, unsigned char* data, unsigned long length, int nThreads) {
	register int i;

	//init pthreads
	thread_data_mutex_sync * tds;
	pthread_attr_t tAttr;
	pthread_mutex_t mLock;
	pthread_mutexattr_t mAttr;
	volatile thread_flag flag = init;
	struct timeval start, stop;

	//Init memory
	register int res = posix_memalign(&tds, 128/8, sizeof(thread_data_mutex_sync) * nThreads);
	if (res) die("Allocation error for thread_data_mutex_sync");
	pthread_attr_init(&tAttr);
	pthread_attr_setdetachstate(&tAttr, PTHREAD_CREATE_JOINABLE);
	pthread_mutexattr_init(&mAttr);
	pthread_mutex_init(&mLock, &mAttr);

	//determine stride
	unsigned long stride, left;
	stride = length / nThreads;
	left = length % nThreads;


	//Start time, because after this forloop, the threads have already done some work...
	gettimeofday(&start, 0);

	//Create threads
	for(i=nThreads-1; i>=0; i--) {
		tds[i].id = i;
		tds[i].flag = &flag;
		tds[i].data = data;
		tds[i].from = i*stride;
		tds[i].shared_bin = bin;
		tds[i].lock = &mLock;

		//The last thread gets the extra
		//Note: Inverse order, so last thread is the first :D
		if (i==nThreads-1)
			tds[i].amount = stride + left;
		else
			tds[i].amount = stride;

		//Create & run
		pthread_create(&tds[i].pid, &tAttr, threader_fullmutex, &tds[i]);
	}


	//join
	for(i=0; i<nThreads; i++)
		pthread_join(tds[i].pid, 0);

	//done
	gettimeofday(&stop, 0);

	//free
	free(tds);
	pthread_attr_destroy(&tAttr);
	pthread_mutex_destroy(&mLock);
	pthread_mutexattr_destroy(&mAttr);

	//return time
	return (double)(stop.tv_sec-start.tv_sec) + (stop.tv_usec-start.tv_usec)/1.0e6;
}

void * threader_atomic(thread_data_mutex_sync * data) {
	register unsigned long length, from, i;
	register unsigned char * img;
	img = data->data;
	volatile unsigned int * bin = data->shared_bin;

	//Init registers
	from = data->from;
	length = data->amount;

	//run
	*data->flag = work;

	for(i=from; i<from+length; i++) {
		__sync_fetch_and_add(&bin[img[i]], 1);
	}

	//done
	*data->flag = done;

	return 0;

}

double run_atomic(unsigned int *bin, unsigned char* data, unsigned long length, int nThreads) {
	register int i;

	//init pthreads
	thread_data_mutex_sync * tds;
	pthread_attr_t tAttr;
	pthread_mutex_t mLock;
	pthread_mutexattr_t mAttr;
	volatile thread_flag flag = init;
	struct timeval start, stop;

	//Init memory
	register int res = posix_memalign(&tds, 128/8, sizeof(thread_data_mutex_sync) * nThreads);
	if (res) die("Allocation error for thread_data_mutex_sync");
	pthread_attr_init(&tAttr);
	pthread_attr_setdetachstate(&tAttr, PTHREAD_CREATE_JOINABLE);
	pthread_mutexattr_init(&mAttr);
	pthread_mutex_init(&mLock, &mAttr);

	//determine stride
	unsigned long stride, left;
	stride = length / nThreads;
	left = length % nThreads;


	//Start time, because after this forloop, the threads have already done some work...
	gettimeofday(&start, 0);

	//Create threads
	for(i=nThreads-1; i>=0; i--) {
		tds[i].id = i;
		tds[i].flag = &flag;
		tds[i].data = data;
		tds[i].from = i*stride;
		tds[i].shared_bin = bin;
		tds[i].lock = &mLock;

		//The last thread gets the extra
		//Note: Inverse order, so last thread is the first :D
		if (i==nThreads-1)
			tds[i].amount = stride + left;
		else
			tds[i].amount = stride;

		//Create & run
		pthread_create(&tds[i].pid, &tAttr, threader_atomic, &tds[i]);
	}


	//join
	for(i=0; i<nThreads; i++)
		pthread_join(tds[i].pid, 0);

	//done
	gettimeofday(&stop, 0);

	//free
	free(tds);
	pthread_attr_destroy(&tAttr);
	pthread_mutex_destroy(&mLock);
	pthread_mutexattr_destroy(&mAttr);

	//return time
	return (double)(stop.tv_sec-start.tv_sec) + (stop.tv_usec-start.tv_usec)/1.0e6;
}

//void * threader_stm(thread_data_mutex_sync * data) {
//	register unsigned long length, from, i;
//	register unsigned char * img;
//	img = data->data;
//
//	//Init registers
//	from = data->from;
//	length = data->amount;
//
//	//run
//	*data->flag = work;
//
//	for(i=from; i<from+length; i++) {
//		__transaction_atomic { data->shared_bin[img[i]]++; }
//	}
//
//	//done
//	*data->flag = done;
//
//	return 0;
//
//}
//
//double run_stm(unsigned int *bin, unsigned char* data, unsigned long length, int nThreads) {
//	register int i;
//
//	//init pthreads
//	thread_data_mutex_sync * tds;
//	pthread_attr_t tAttr;
//	pthread_mutex_t mLock;
//	pthread_mutexattr_t mAttr;
//	volatile thread_flag flag = init;
//	struct timeval start, stop;
//
//	//Init memory
//	register int res = posix_memalign(&tds, 128/8, sizeof(thread_data_mutex_sync) * nThreads);
//	if (res) die("Allocation error for thread_data_mutex_sync");
//	pthread_attr_init(&tAttr);
//	pthread_attr_setdetachstate(&tAttr, PTHREAD_CREATE_JOINABLE);
//	pthread_mutexattr_init(&mAttr);
//	pthread_mutex_init(&mLock, &mAttr);
//
//	//determine stride
//	unsigned long stride, left;
//	stride = length / nThreads;
//	left = length % nThreads;
//
//
//	//Start time, because after this forloop, the threads have already done some work...
//	gettimeofday(&start, 0);
//
//	//Create threads
//	for(i=nThreads-1; i>=0; i--) {
//		tds[i].id = i;
//		tds[i].flag = &flag;
//		tds[i].data = data;
//		tds[i].from = i*stride;
//		tds[i].shared_bin = bin;
//		tds[i].lock = &mLock;
//
//		//The last thread gets the extra
//		//Note: Inverse order, so last thread is the first :D
//		if (i==nThreads-1)
//			tds[i].amount = stride + left;
//		else
//			tds[i].amount = stride;
//
//		//Create & run
//		pthread_create(&tds[i].pid, &tAttr, threader_stm, &tds[i]);
//	}
//
//
//	//join
//	for(i=0; i<nThreads; i++)
//		pthread_join(tds[i].pid, 0);
//
//	//done
//	gettimeofday(&stop, 0);
//
//	//free
//	free(tds);
//	pthread_attr_destroy(&tAttr);
//	pthread_mutex_destroy(&mLock);
//	pthread_mutexattr_destroy(&mAttr);
//
//	//return time
//	return (double)(stop.tv_sec-start.tv_sec) + (stop.tv_usec-start.tv_usec)/1.0e6;
//}


void * threader_sema(thread_data_mutex_sync * data) {
	register unsigned long length, from, i;
	register unsigned char * img;
	img = data->data;

	//Init registers
	from = data->from;
	length = data->amount;

	//run
	*data->flag = work;

	for(i=from; i<from+length; i++) {
		sem_wait(data->semaphore);
		data->shared_bin[img[i]]++;
		sem_post(data->semaphore);
	}

	//done
	*data->flag = done;

	return 0;

}

double run_sema(unsigned int *bin, unsigned char* data, unsigned long length, int nThreads) {
	register int i;

	//init pthreads
	thread_data_mutex_sync * tds;
	pthread_attr_t tAttr;
	sem_t semap;
	volatile thread_flag flag = init;
	struct timeval start, stop;

	//Init memory
	register int res = posix_memalign(&tds, 128/8, sizeof(thread_data_mutex_sync) * nThreads);
	if (res) die("Allocation error for thread_data_mutex_sync");
	pthread_attr_init(&tAttr);
	pthread_attr_setdetachstate(&tAttr, PTHREAD_CREATE_JOINABLE);
	sem_init(&semap, 0, 1);

	//determine stride
	unsigned long stride, left;
	stride = length / nThreads;
	left = length % nThreads;


	//Start time, because after this forloop, the threads have already done some work...
	gettimeofday(&start, 0);

	//Create threads
	for(i=nThreads-1; i>=0; i--) {
		tds[i].id = i;
		tds[i].flag = &flag;
		tds[i].data = data;
		tds[i].from = i*stride;
		tds[i].shared_bin = bin;
		tds[i].semaphore = &semap;

		//The last thread gets the extra
		//Note: Inverse order, so last thread is the first :D
		if (i==nThreads-1)
			tds[i].amount = stride + left;
		else
			tds[i].amount = stride;

		//Create & run
		pthread_create(&tds[i].pid, &tAttr, threader_sema, &tds[i]);
	}


	//join
	for(i=0; i<nThreads; i++)
		pthread_join(tds[i].pid, 0);

	//done
	gettimeofday(&stop, 0);

	//free
	free(tds);
	pthread_attr_destroy(&tAttr);
	sem_destroy(&semap);

	//return time
	return (double)(stop.tv_sec-start.tv_sec) + (stop.tv_usec-start.tv_usec)/1.0e6;
}

void * threader_sync_atomic(thread_data_mutex_sync * data) {
	register unsigned long length, from, i;
	register unsigned char * img;
	img = data->data;



	//Init registers
	from = data->from;
	length = data->amount;

	//run
	*data->flag = work;

	for(i=from; i<from+length; i++) {
		data->bin[img[i]]++;
	}

	//dump
	*data->flag = saving;
	for(i=0; i<256;i++) {
		__sync_fetch_and_add(&data->shared_bin[i] , data->bin[i]);
	}

	//done
	*data->flag = done;

	return 0;

}

double run_sync_atomic(unsigned int *bin, unsigned char* data, unsigned long length, int nThreads) {
	register int i;

	//init pthreads
	thread_data_mutex_sync * tds;
	pthread_attr_t tAttr;
	volatile thread_flag flag = init;
	struct timeval start, stop;

	//Init memory
	register int res = posix_memalign(&tds, 128/8, sizeof(thread_data_mutex_sync) * nThreads);
	if (res) die("Allocation error for thread_data_mutex_sync");
	memset(tds, 0, sizeof(thread_data_mutex_sync) * nThreads);
	pthread_attr_init(&tAttr);
	pthread_attr_setdetachstate(&tAttr, PTHREAD_CREATE_JOINABLE);

	//determine stride
	unsigned long stride, left;
	stride = length / nThreads;
	left = length % nThreads;


	//Start time, because after this forloop, the threads have already done some work...
	gettimeofday(&start, 0);

	//Create threads
	for(i=nThreads-1; i>=0; i--) {
		tds[i].id = i;
		tds[i].flag = &flag;
		tds[i].data = data;
		tds[i].from = i*stride;
		tds[i].shared_bin = bin;

		//The last thread gets the extra
		//Note: Inverse order, so last thread is the first :D
		if (i==nThreads-1)
			tds[i].amount = stride + left;
		else
			tds[i].amount = stride;

		//Create & run
		pthread_create(&tds[i].pid, &tAttr, threader_sync_atomic, &tds[i]);
	}


	//join
	for(i=0; i<nThreads; i++)
		pthread_join(tds[i].pid, 0);

	//done
	gettimeofday(&stop, 0);

	//free
	free(tds);
	pthread_attr_destroy(&tAttr);

	//return time
	return (double)(stop.tv_sec-start.tv_sec) + (stop.tv_usec-start.tv_usec)/1.0e6;
}

unsigned int * getTrueBin(unsigned char *data, unsigned long length, double * time) {
	struct timeval start, stop;
	unsigned int * bin;
	register int res = posix_memalign((void**)&bin, 128/8, 256*sizeof(unsigned int));
	if (res)
		die("Allocation error for true bin");

	memset(bin, 0, 256*sizeof(unsigned int));

	int i;
	gettimeofday(&start, 0);
	for (i=0; i<length; i++)
		bin[data[i]]++;
	gettimeofday(&stop, 0);

	if(time)
		*time=(double)(stop.tv_sec-start.tv_sec) + (stop.tv_usec-start.tv_usec)/1.0e6;

	return bin;
}

/**
 * returns true upon success
 */
char bincorrect(unsigned int * bin, unsigned int * bin2) {
	int i;
	char res = 1;
	for(i=0; i<255; i++) {
		if (!(bin[i] == bin2[i])) {
			printf("fail on %i: %i != %i\n", i, bin[i], bin2[i]);
			res = 0;
		}
	}

	return res;
}

int main(int nargs, char ** args) {
	unsigned char * data;
	unsigned long length;
	unsigned int * bin;
	unsigned int * truebin;
	double seqTime;
	int nThreads;

	//Init data & lenght
	parseInput(nargs, args, &data, &length, &nThreads);

	//calculate norm
	truebin = getTrueBin(data, length, &seqTime);
	printf("seq:      \t%f\n", seqTime);

	//Init bin
	//Note: bin is a case of extreme vector shit
	register int res = posix_memalign((void**)&bin, 128/8, 256*sizeof(unsigned int));
	if (res)
		die("Allocation error for bin");

	int i;
	for(i=0; i < 256;) { //Vectorize!
		bin[i++] = 0;
		bin[i++] = 0;
		bin[i++] = 0;
		bin[i++] = 0;
	}


	//run_exp
	double time;
	time = run_sync_mutex(bin, data, length, nThreads);
	if (!bincorrect(bin, truebin))
		die("Bin failure\n");
	printf("run_sync_mutex:\t%f\n", time);

	//reset bin
	for(i=0; i < 256;) { //Vectorize!
		bin[i++] = 0;
		bin[i++] = 0;
		bin[i++] = 0;
		bin[i++] = 0;
	}

	//run_exp2
	double time2;
	time2 = run_full_mutex(bin, data, length, nThreads);
	if (!bincorrect(bin, truebin))
		die("Bin failure\n");
	printf("run_full_mutex:\t%f\n", time2);

	//reset bin
	for(i=0; i < 256;) { //Vectorize!
		bin[i++] = 0;
		bin[i++] = 0;
		bin[i++] = 0;
		bin[i++] = 0;
	}

	//run_exp2
	double time3;
	time3 = run_atomic(bin, data, length, nThreads);
	if (!bincorrect(bin, truebin))
		die("Bin failure\n");
	printf("run_atomic:\t%f\n", time3);

	//reset bin
		for(i=0; i < 256;) { //Vectorize!
			bin[i++] = 0;
			bin[i++] = 0;
			bin[i++] = 0;
			bin[i++] = 0;
		}

	//run_exp2
	double time4;
	time4 = run_sema(bin, data, length, nThreads);
	if (!bincorrect(bin, truebin))
		die("Bin failure\n");
	printf("run_semap:\t%f\n", time4);


	//reset bin
		for(i=0; i < 256;) { //Vectorize!
			bin[i++] = 0;
			bin[i++] = 0;
			bin[i++] = 0;
			bin[i++] = 0;
		}

	//run_exp5
	double time5;
	time5 = run_sync_atomic(bin, data, length, nThreads);
	if (!bincorrect(bin, truebin))
		die("Bin failure\n");
	printf("run_sync_atomic:\t%f\n", time5);

	//store results
	FILE *f;
	f = fopen("/var/scratch/ppp1620/histo.csv", "a");

	if (!f)
		return 0;

	if (ftell(f)==0)
		fprintf(f, "nThreads; Size; seq; sync_mutex; full_mutex; atomic; semap; sync_mutex\n");

	fprintf(f, "%i; %lu; %f; %f ;%f;%f;%f;%f\n", nThreads, length, seqTime, time, time2, time3, time4, time5);
	fclose(f);


	//free
	free(bin);
	free(truebin);


	return 0;
}
