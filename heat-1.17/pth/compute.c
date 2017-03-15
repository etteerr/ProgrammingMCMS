#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <pthread.h>
//#include <emmintrin.h>
#include <x86intrin.h>
#include <unistd.h>
#include <sys/time.h>
//#include <semaphore.h>
#include "../src/compute.h"
#include "../src/output.h"

#define FPOPS_PER_POINT_PER_ITERATION (                 \
        1     /* current point 1 mul */ +               \
        3 + 1 /* direct neighbors 3 adds + 1 mul */ +   \
        3 + 1 /* diagonal neighbors 3 adds + 1 mul */ + \
        2     /* final add */ +                         \
        1     /* difference old/new */                  \
        )

//Pthreads
//If commented, used env variable ED_NUM_THREADS
//TODO: Remove mutex from work After report ;)
#define ED_NUM_THREADS 32
#define Verbose_work_request 0
#define Sleep_work_request 0


#define checkSanity 0
#define sanityPGM 0
#define sideCoef 0.58578643762

//Static pointers to data
double * hm_data = 0;
double * hm_prevData = 0;
double * hm_coeff = 0;

//Fixed global dimentions
size_t hm_row = 0;
size_t hm_column = 0;

//defines
void do_compute_select(double*hm_data, double *hm_prevData, unsigned int lowerbound, unsigned int upperbound);

//work queue
struct work {
	unsigned int from;
	unsigned int amount;
	unsigned int step;

	//work lock
	pthread_mutex_t lock;
};


struct work_queue {
	unsigned int buf_ptr;
	unsigned int buf_size;

	struct work * queue_buffer;
	pthread_spinlock_t slock;
};

inline struct work * get_work(struct work_queue *wq, unsigned int max_step, unsigned long id) {
	//Lock workqueue
	pthread_spin_lock(&wq->slock);
	if (Verbose_work_request) printf("Request work from %lu (%i) ", id, wq->buf_ptr);

	//Lock workk
	if (pthread_mutex_trylock(&wq->queue_buffer[wq->buf_ptr].lock)) {
		//lock failed. already someone working on this job?
		pthread_spin_unlock(&wq->slock);
		if (Verbose_work_request) printf("failed: work locked.\n");
		if (Sleep_work_request) sleep(1);
		return 0;
	}

	//Check if work is needed
	if (wq->queue_buffer[wq->buf_ptr].step >= max_step) {
		//Release locks and return 0
		pthread_mutex_unlock(&wq->queue_buffer[wq->buf_ptr].lock);
		pthread_spin_unlock(&wq->slock);
		if (Verbose_work_request) printf("failed: work exceeds step limit (%u >= %u).\n", wq->queue_buffer[wq->buf_ptr].step, max_step);
		if (Sleep_work_request) sleep(1);
		return 0;
	}

	struct work * nw = &wq->queue_buffer[wq->buf_ptr];
	//inc work pointer
	wq->buf_ptr = (wq->buf_ptr+1)%wq->buf_size;
	if (Verbose_work_request) printf("Succeeded! New buf_ptr: %u\n", wq->buf_ptr);
	if (Sleep_work_request) sleep(1);

	//unlock workqueue
	pthread_spin_unlock(&wq->slock);

	return nw;
}

inline void release_work(struct work * w) {
	w->step++;
	pthread_mutex_unlock(&w->lock);
}

//thread struct info
struct thread_info {
	pthread_t thread_handle; //ul
	pthread_t thread_id; //0 to X

	//data
	double * hm_data;
	double * hm_prevData;

	//queue
	struct work_queue * q;

	//private vars
	unsigned long cStep;
	volatile unsigned long target_step; //Use __sync_fetch_and_add ? yes

	//Notify
	unsigned int * working_threads; //number

};


#define atomic_read(POINTER) __sync_fetch_and_add(POINTER, 0)
#define atomic_add(POINTER,VAL) __sync_fetch_and_add(POINTER,VAL)
#define atomic_sub(POINTER,VAL) __sync_fetch_and_sub(POINTER,VAL)


//Getters and setters
// Inline for performance
#define hm_get(R,C) hm_prevData[hm_column * (R+1) + C]
#define hm_geta(R,C) (hm_prevData +hm_column * (R+1) + C)
//__attribute__((assume_aligned(16)))
//inline double hm_get(int row, int column) {
//	if (checkSanity) {
//		if (row > (int)hm_row || column > hm_column || row < -1 || column < 0)
//			printf("Invalid coördinate requested (%i,%i)\n", row, column);
//	}
//	return hm_prevData[hm_column * (row+1) + column];
//}

#define hm_get_current(R,C) hm_data[hm_column * (R+1) + C]
//__attribute__((assume_aligned(16)))
//inline double hm_get_current(size_t  row, size_t  column) {
//	if (checkSanity) {
//		if (row > hm_row || column > hm_column || row < 0 || column < 0)
//			printf("Invalid coördinate requested (current) (%i,%i)\n", row, column);
//	}
//	return hm_data[hm_column * (row+1) + column];
//}

#define hm_set(R,C,V) hm_data[hm_column * (R+1) + C] = V
//inline void hm_set(size_t row, size_t column, double val) {
//	if (checkSanity) {
//		if (row > hm_row || column > hm_column || row < 0 || column < 0)
//			printf("Invalid coördinate set (%i,%i)\n", row, column);
//	}
//	*(hm_data + hm_column * (row+1) + column) = val;
//}

#define hm_get_coeff(R,C) hm_coeff[hm_column * R + C]
#define hm_get_coeffa(R,C) (hm_coeff + hm_column * R + C)
//__attribute__((assume_aligned(16)))
//inline double hm_get_coeff(size_t row, size_t column) {
//	return hm_coeff[hm_column * row + column];
//}

#define hm_set_coeff(R,C,V) hm_coeff[hm_column * R + C] = V
//__attribute__((assume_aligned(16)))
//inline void hm_set_coeff(size_t row, size_t column, double val) {
//	hm_coeff[hm_column * row + column] = val;
//}

//Inline is faster in this case
//#define hm_swap { double * tmp = hm_data; hm_data = hm_prevData; hm_prevData = tmp; }
inline void hm_swap() {
	double * tmp = hm_data;
	hm_data = hm_prevData;
	hm_prevData = tmp;
}


inline void renderImageMemory(size_t key, double min, double max) {
	begin_picture(key, hm_column, (hm_row+2)*2, min, max);

	for(int i = 0; i < hm_column; i++) {
		for (int j =0; j < (hm_row+2)*2; j++)
			draw_point(i,j,hm_data[i*hm_column + j]);
	}
	end_picture();
}

void calcStatistics(struct parameters* p,struct results* r) {
    	//convergence
    	r->maxdiff = 0;
    	r->tmax = 0;
    	r->tmin = 1e9;
    	r->tavg = 0;
    	double sum = 0.0;
    	for(int i = 0; i < hm_row; i++)
    		for (int j =0; j < hm_column; j++) {
				//max
				if (hm_get_current(i,j) > r->tmax) r->tmax = hm_get_current(i,j);
				//min
				if (hm_get_current(i,j) < r->tmin) r->tmin = hm_get_current(i,j);
				//maxdiff
				if (fabs(hm_get_current(i,j)-hm_get(i,j)) > r->maxdiff) r->maxdiff = fabs(hm_get_current(i,j)-hm_get(i,j));
				//searchme
				sum += hm_get_current(i,j);
    		}

    	r->tavg += sum/(hm_column*hm_row);
}

void checkAnswer(unsigned int i, unsigned int j, double cn, double cd) {
	double cnext = (1-hm_get_coeff(i,j))*sideCoef;
	double cdiag = 1-hm_get_coeff(i,j)-cnext;
	cnext *= .25;
	cdiag *= .25;
	static char first = 1;
	double expected =
			//SUm neighbours (direct)
			((hm_get(i-1,j) + hm_get(i+1,j) + hm_get(i,j-1) + hm_get(i,j+1)))*cnext +
			//Sum diag
			((hm_get(i-1,j-1) + hm_get(i+1,j+1) + hm_get(i+1,j-1) + hm_get(i-1,j+1)))*cdiag +
			//Add current
			hm_get(i,j) * hm_get_coeff(i,j);

	if (expected != hm_get_current(i,j)/* || cnext!=cn || cdiag!=cd */) {
		if  (first) {
			first = 0;
			fprintf(stderr,"unexpected value at (i, j) Value (diff)\t& cnext, cdiag (diff, diff)\n");
		}
		fprintf(stderr,"unexpected value at (%u, %u) %.8f (%.8f)\t& %.8f, %.8f (%.8f, %.8f)\n",i,j, hm_get_current(i,j), hm_get_current(i,j) - expected, cn,cd, cn-cnext, cd-cdiag);
//		fprintf(stderr,"unexpected value at (%u, %u) %.8f (%.8f)\t& %.8f, %.8f (%.8f, %.8f)\n",i,j, hm_get_current(i,j), expected, cn,cd, cnext, cdiag);
		sleep(1);
	}
}

void * threader(struct thread_info *args) {
	//init local data (mainly for use in registery of compiler magic does its jobs)
	char working = 0;
	double *hmd, *hmp, *hmt;
	struct work * w;

	hmd = args->hm_data;
	hmp = args->hm_prevData;
	unsigned long cStep = args->cStep;

	while(1) {

		//get work
		w = get_work(args->q, args->target_step, args->thread_id);

		//If we got work, do it
		if (w!=0) {
			//report
//			printf("Thread %i:\n"
//					"\tcStep: %lu\n"
//					"\ttStep: %lu\n"
//					"\tfrom:  %lu\n", args->thread_id,
//					args->cStep,
//					args->target_step,
//					w->from);
			if (!working) {
				working = 1;
				atomic_add(args->working_threads,1);
			}

			//inc local step and swap buffers accordingly
			if (w->step != cStep) {
				while(w->step > cStep) {
					cStep++;

					//Swap
					hmt = hmd;
					hmd = hmp;
					hmp = hmt;
				}
			}
			do_compute_select(hmd, hmp,w->from, w->from+w->amount);

			release_work(w);
		}else {
			//No work left
			//Notify main thread
			if (working) {
				working = 0;
				atomic_sub(args->working_threads, 1);
			}

			//Update when we have the time (from registery to pointed memory)
			args->cStep = cStep;

			sched_yield();
		}
	}

	return 0;
}


void do_compute(const struct parameters *p, struct results *r) {
	//Def
	struct timeval start,end;
    int i,j;
    unsigned int num_threads;
    struct thread_info * thread_handles;

    //get thread num
#ifdef ED_NUM_THREADS
    num_threads = ED_NUM_THREADS;
#else
    num_threads = atoi(getenv("ED_NUM_THREADS"));
#endif

    //quit if invalid num_threads
    if (num_threads <= 0) {
    	die("Invalid num_threads, run export ED_NUM_THREADS=#");
    }

    //Mutex attr
    pthread_mutexattr_t mAttr;
    pthread_mutexattr_init(&mAttr);

    //Init workqueue
    unsigned int jobs_per_step = 4*num_threads;
    unsigned int job_size = hm_row/jobs_per_step;
    if (job_size*jobs_per_step < hm_row) jobs_per_step++;
    struct work_queue workq;
    workq.buf_ptr = 0;
    workq.buf_size = jobs_per_step;
    workq.queue_buffer = malloc(sizeof(struct work)*workq.buf_size);
    pthread_spin_init(&workq.slock, 0);

    //fill the queue
    for(i = 0; i < jobs_per_step; i++) {
    	//From
    	workq.queue_buffer[i].from = i*job_size;

    	//To (amount)
    	if ((i+1)*job_size >= hm_row)
    		workq.queue_buffer[i].amount = hm_row-(i*job_size);
    	else
    		workq.queue_buffer[i].amount = job_size;

    	//Init step & finally lock
    	workq.queue_buffer[i].step = 0;
    	pthread_mutex_init(&workq.queue_buffer[i].lock, &mAttr);
    }

    //Init pthreads
    thread_handles = malloc(num_threads*sizeof(struct thread_info));

    if (thread_handles==0)
    	exit(100);

    //Pthread attr
    pthread_attr_t pAttr;
    pthread_attr_init(&pAttr);
    pthread_attr_setscope(&pAttr, PTHREAD_SCOPE_SYSTEM);
    pthread_attr_setdetachstate(&pAttr, PTHREAD_CREATE_DETACHED);
//    pthread_attr_setstacksize(&pAttr, 512*1000); //512 kbytes

    //Create theads & their shared condition var
    unsigned int working_threads = 0;
    pthread_condattr_t cAttr;
    pthread_cond_t prCond;
    pthread_condattr_init(&cAttr);
    pthread_cond_init(&prCond, &cAttr);

    for(i=0; i<num_threads; i++) {
    	//Create worker info
    	thread_handles[i].thread_id = i;
    	thread_handles[i].hm_data = hm_data;
    	thread_handles[i].hm_prevData = hm_prevData;
    	thread_handles[i].q = &workq;
    	thread_handles[i].cStep = 0;
    	thread_handles[i].target_step = 0;
    	thread_handles[i].working_threads = &working_threads;


    	if (pthread_create(&thread_handles[i].thread_handle, &pAttr, (void*)&threader, (void*)&thread_handles[i]))
    		exit(110 + i);
    }
    unsigned int period = p->period;
    //main calc loop
    while(r->niter<p->maxiter) {

    	//Calc & time
    	gettimeofday(&start, 0);

    	if (period + r->niter > p->maxiter)
    		period = p->maxiter - r->niter;

    	//Set state (all working)
//    	atomic_add(&working_threads, num_threads);

    	//Increment (by period) all thread max iters
    	for (i=0; i<num_threads; i++) {
    		atomic_add(&thread_handles[i].target_step, period);
    	}

    	//Wait for atleast one thread to be working
    	while(!atomic_read(&working_threads)) sched_yield();

    	//Wait for working threads to be 0
		while(atomic_read(&working_threads)) sched_yield();

		r->niter += period;

    	gettimeofday(&end, 0);

    	//Calc statistics
    	calcStatistics(p, r);

    	//calc time diff and store in p
    	r->time += ((end.tv_sec-start.tv_sec)*1e6 + (end.tv_usec - start.tv_usec))/1e6;

    	//report results
    	report_results(p,r);

    	//check stop condition
    	if (r->maxdiff < p->threshold) break;
    }
    //**** Cleanup *****
    //Threads
    for(i = 0; i <num_threads; i++) {
    	pthread_cancel(thread_handles[i].thread_handle);
    	pthread_mutex_destroy(&workq.queue_buffer[i].lock);
    }

    //pthread structures
    pthread_mutexattr_destroy(&mAttr);
    pthread_attr_destroy(&pAttr);
    pthread_spin_destroy(&workq.slock);
    pthread_condattr_destroy(&cAttr);
    pthread_cond_destroy(&prCond);

    //Mallocs
    free(workq.queue_buffer);
    free(thread_handles);

	//store results
	FILE *f;
	f = fopen("/var/scratch/ppp1620/pth_heat.csv", "a");

	if (!f)
		return;

	if (ftell(f)==0)
		fprintf(f, "nThreads; width; height; iters; time; flops\n");

	fprintf(f, "%i; %lu; %lu; %lu; %f; %f\n", num_threads, p->M, p->N, r->niter, r->time,
	(double) p->N * (double)p->M *
	(double)(r->niter * FPOPS_PER_POINT_PER_ITERATION +
	(double) r->niter / p->period) / r->time);

	fclose(f);


}

void do_compute_select(double*hm_data, double *hm_prevData, unsigned int lowerbound, unsigned int upperbound)
{
	//Intresting way of declaring privates
	double cdiag, cnext;
	unsigned int i, j;
	// Calculate the inner sums

	for ( i = lowerbound; i < upperbound ; i++)
	{

		//left
		cnext = (1-hm_get_coeff(i,0))*sideCoef;
		cdiag = 1-hm_get_coeff(i,0)-cnext;
		hm_set(i,0,
				//SUm neighbours (direct)
				((hm_get(i-1,0) + hm_get(i+1,0) + hm_get(i,hm_column-1) + hm_get(i,1))*.25)*cnext +
				//Sum diag
				((hm_get(i-1,hm_column-1) + hm_get(i+1,1) + hm_get(i+1,hm_column-1) + hm_get(i-1,1))*.25)*cdiag +
				//Add current
				hm_get(i,0) * hm_get_coeff(i,0)
		);

		//Middle
		for ( j = 1; j < hm_column-1; j++)
		{
			cnext = (1-hm_get_coeff(i,j))*sideCoef;
			cdiag = 1-hm_get_coeff(i,j)-cnext;
			hm_set(i,j,
					//SUm neighbours (direct)
					((hm_get(i-1,j) + hm_get(i+1,j) + hm_get(i,j-1) + hm_get(i,j+1))*.25)*cnext +
					//Sum diag
					((hm_get(i-1,j-1) + hm_get(i+1,j+1) + hm_get(i+1,j-1) + hm_get(i-1,j+1))*.25)*cdiag +
					//Add current
					hm_get(i,j) * hm_get_coeff(i,j)
			);
		}

		//right
		cnext = (1-hm_get_coeff(i,(hm_column-1)))*sideCoef;
		cdiag = 1-hm_get_coeff(i,(hm_column-1))-cnext;
		hm_set(i,(hm_column-1),
				//SUm neighbours (direct)
				((hm_get(i-1,(hm_column-1)) + hm_get(i+1,(hm_column-1)) + hm_get(i,(hm_column-1)-1) + hm_get(i,0))*.25)*cnext +
				//Sum diag
				((hm_get(i-1,(hm_column-1)-1) + hm_get(i+1,0) + hm_get(i+1,(hm_column-1)-1) + hm_get(i-1,0))*.25)*cdiag +
				//Add current
				hm_get(i,(hm_column-1)) * hm_get_coeff(i,(hm_column-1))
		);
	}
}

/**
 * hm_init_map
 *  returns 0 on succes
 *
 *  Note that the allocation is padded with two extra rows per map
 *  Making x(-1) pointers valid as well as x(size).
 *
 */
int hm_init_map(struct parameters * p) {
	//Set vars
	size_t row = p->N;
	size_t column = p->M;

	//Init tmps
	void * tmp1 = 0;
	void * tmp3 = 0;

	//Allocate alligned memory
	int result1 = posix_memalign(&tmp1, 128/8, (sizeof(double)*(row+2)*(column))*2);
	int result2 = posix_memalign(&tmp3, 128/8, sizeof(double)*row*column);

	if (result1 || result2) {
		hm_data = 0;
		if (HM_VERBOSE) printf("Error while allocating memory in hm_init_map: %i, %i", result1, result2);
	}else {
		hm_data     =   (double*) (tmp1);
		hm_prevData =   (double*) (tmp1+column*(row+2)*sizeof(double)); //skip two ghost rows and a map.
		hm_coeff    =   (double*) tmp3;
		hm_row     =   row;
		hm_column  =   column;

		//Copy data
		// Heatmap
		memcpy((void*)(hm_prevData+column), (void*)p->tinit, sizeof(double)*row*column);
		memcpy((void*)(hm_data+column), (void*)p->tinit, sizeof(double)*row*column); // Fixes strange problems
		// Prevdata
		memcpy((void*)hm_coeff, (void*)p->conductivity, sizeof(double)*row*column);

		//Ghost rows
		// prevData
		memcpy((void*)(hm_prevData), (void*)p->tinit, sizeof(double)*column);
		memcpy((void*)(hm_prevData+(row+1)*column), (void*)(p->tinit+row*column-column), sizeof(double)*column);

		// current map
		memcpy((void*)(hm_data), (void*)p->tinit, sizeof(double)*column);
		memcpy((void*)(hm_data+(row+1)*column), (void*)p->tinit, sizeof(double)*column);

		//print test
//		begin_picture(666, column, row, p->io_tmin, p->io_tmax);
//			for (int i = 0; i < row; i++)
//				for(int j = 0; j < column; j++)
//					draw_point(i,j,p->tinit[i*column + j]);
//
//		end_picture();
//		begin_picture(667, column, row, p->io_tmin, p->io_tmax);
//			for (int i = 1; i < row+1; i++)
//				for(int j = 0; j < column; j++)
//					draw_point(i,j,hm_get(i,j));
//
//		end_picture();

		if (checkSanity) {
	//		Test map
			for (int i = 0; i < (column*(row+2))*2; i++)
				hm_data[i] = 1.0;// ((double)rand()/RAND_MAX)*(double)rand();
		}
	}//end of allocate succes if
	return result1 | result2;
}
