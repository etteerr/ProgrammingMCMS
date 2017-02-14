#include <sys/time.h>
#include <math.h>
#include <stdlib.h>
#include "../src/compute.h"

//Static pointers to data
double * hm_data = 0;
double * hm_prevData = 0;
double * hm_coeff = 0;

//Fixed global dimentions
size_t hm_rows = 0;
size_t hm_columns = 0;

void do_compute(const struct parameters* p, struct results *r)
{
    
}


/**
 * hm_init_map
 *  returns 0 on succes
 */
int hm_init_map(struct parameters * p) {
    //Set vars
    size_t rows = p->N;
    size_t columns = p->M;

    //Init tmps
    void * tmp1 = 0;
    void * tmp2 = 0;
    void * tmp3 = 0;

    //Allocate alligned memory
    int result1 = posix_memalign(&tmp1, 128/8, sizeof(double)*rows*columns);
    int result2 = posix_memalign(&tmp2, 128/8, sizeof(double)*rows*columns);
    int result3 = posix_memalign(&tmp3, 128/8, sizeof(double)*rows*columns);

    if (result1 || result2 || result2) {
        hm_data = 0;
        if (HM_VERBOSE) printf("Error while allocating memory in hm_init_map: %i, %i, %i", result1, result2, result3);
    }else {
        hm_data     =   (double*) tmp1;
        hm_prevData =   (double*) tmp2;
        hm_coeff    =   (double*) tmp3;
        hm_rows     =   rows;
        hm_columns  =   columns;

        memcpy((void*)hm_prevData, (void*)p->tinit, sizeof(double)*rows*columns);
        memcpy((void*)hm_coeff, (void*)p->conductivity, sizeof(double)*rows*columns);

        free((void*)p->tinit);
        free((void*)p->conductivity);
        
    }

    return result1 | result2 | result3;
}

//Getters and setters
// Inline for performance
inline double hm_get(size_t rows, size_t columns) {
    return hm_prevData[hm_columns*rows + columns];
}

inline double hm_get_current(size_t rows, size_t columns) {
    return hm_data[hm_columns*rows + columns];
}

inline void hm_set(size_t rows, size_t columns, double val) {
    hm_data[hm_columns*rows + columns] = val;
}

inline double hm_get_coeff(size_t rows, size_t columns) {
    return hm_coeff[hm_columns*rows + columns];
}

inline void hm_set_coeff(size_t rows, size_t columns, double val) {
    hm_coeff[hm_columns*rows + columns] = val;
}

inline void hm_swap() {
    double * tmp = hm_data;
    hm_data = hm_prevData;
    hm_prevData = tmp;
}
