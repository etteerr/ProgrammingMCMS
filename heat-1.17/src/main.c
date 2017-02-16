#include "compute.h"
#include <sys/time.h>
#include <unistd.h>

#define convThreshold 1.0e-3

int main(int argc, char **argv)
{
    struct parameters p;
    struct results r;
    struct timeval start,end;


    read_parameters(&p, argc, argv);

    //Allocate& init heatmap
    if (hm_init_map(&p))
        return 1;

    //set time
    r.time = 0;

    //Main loop
    int j;
    while(r.niter<p.maxiter) {

    	//Calc & time
    	gettimeofday(&start, 0);
    	for (j=0; j<p.period; j++ )
    		do_compute(&p, &r);
    	gettimeofday(&end, 0);

    	//Calc statistics
    	calcStatistics(&p, &r);

    	//calc time diff and store in p
    	r.time += ((end.tv_sec-start.tv_sec)*1e6 + (end.tv_usec - start.tv_usec))/1e6;

    	//report results
    	report_results(&p,&r);

    	//check stop condition
    	if (r.maxdiff < p.threshold) break;
    }

    return 0;
}
