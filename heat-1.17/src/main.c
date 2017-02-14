#include "compute.h"

int main(int argc, char **argv)
{
    struct parameters p;
    struct results r;

    read_parameters(&p, argc, argv);

    //Allocate& init heatmap
    if (hm_init_map(&p))
        return 1;

    printf("Ok");

    do_compute(&p, &r);

    report_results(&p, &r);

    return 0;
}
