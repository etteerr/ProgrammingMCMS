#include "compute.h"
#include <sys/time.h>
#include <unistd.h>
#include <stdio.h>

#define convThreshold 1.0e-3

#define FPOPS_PER_POINT_PER_ITERATION (                 \
        1     /* current point 1 mul */ +               \
        3 + 1 /* direct neighbors 3 adds + 1 mul */ +   \
        3 + 1 /* diagonal neighbors 3 adds + 1 mul */ + \
        2     /* final add */ +                         \
        1     /* difference old/new */                  \
        )

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
    r.niter = 0;

    //wait a sec or 2
    sleep(2);

    printf("%lu, %lu\n", p.maxiter, r.niter);

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

    //Save results appending to file
    FILE *fp;
    fp = fopen("heatResult.csv","a");
    printf("%li", ftell(fp));

    //Check filestream
    if (fp==0) {
    	fprintf(stderr, "Failed to write results to file.\n");
    	return 1;
    }

    //If its a new file, print header
    if (ftell(fp)==0) {
    	fprintf(fp, "Exe;threads;File;Rows;Columns; Iterations; FLOP/s;\n");
    }

    //Print data
    fprintf(fp,"%s;%s;%lu;%lu;%lu;% .6e\n", argv[0], p.filename, p.N, p.M, r.niter,
    		(double) p.N * (double)p.M *
			(double)(r.niter * FPOPS_PER_POINT_PER_ITERATION +
			(double) r.niter / p.period) / r.time);

    fflush(fp);
    fclose(fp);


    return 0;
}
