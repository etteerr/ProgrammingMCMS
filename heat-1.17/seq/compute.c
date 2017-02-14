#include <sys/time.h>
#include <math.h>
#include <stdlib.h>
#include "../src/compute.h"

//Static pointers to data
double * hm_data = 0;
double * hm_prevData = 0;
double * hm_coeff = 0;

//Fixed global dimentions
size_t hm_row = 0;
size_t hm_column = 0;


//Getters and setters
// Inline for performance
inline double hm_get(size_t row, size_t column) {
    return hm_prevData[hm_column*row + column];
}

inline double hm_get_current(size_t row, size_t column) {
    return hm_data[hm_column*row + column];
}

inline void hm_set(size_t row, size_t column, double val) {
    hm_data[hm_column*row + column] = val;
}

inline double hm_get_coeff(size_t row, size_t column) {
    return hm_coeff[hm_column*row + column];
}

inline void hm_set_coeff(size_t row, size_t column, double val) {
    hm_coeff[hm_column*row + column] = val;
}

inline void hm_swap() {
    double * tmp = hm_data;
    hm_data = hm_prevData;
    hm_prevData = tmp;
}

void do_compute(const struct parameters* p, struct results *r)
{
   int i, j;
   int row = p->N;
   int column = p->M;
   double tmp1 =0 ;
   double tmp2 = 0;
   double tmp3 = 0;
   double tmp4 = 0;
   double coeff = 0;
   double sum = 0;
   // Calculate the inner sums       
   for ( i = 1; i < row - 1; i++)
   {
        for ( j = 1; j < column - 1; j++)
        {
            coeff = hm_get_coeff(i,j);
            // Get 4 direct neighbours
            tmp1 = hm_get(i - 1, j);
            tmp2 = hm_get(i + 1, j);
            tmp3 = hm_get(i, j - 1);
            tmp4 = hm_get(i, j + 1);
            sum = 
            // Get 4 diaginal neighbours
            tmp1 = hm_get(i - 1, j - 1);
            tmp2 = hm_get(i + 1, j + 1);
            tmp3 = hm_get(i - 1, j + 1);
            tmp4 = hm_get(i + 1, j - 1);
        }   
   }

   // Calculate the boundary cases, row wise
   //  Top Row
   for ( i = 1; i < column - 1; i++)
   {
        coeff = hm_get_coeff(0,i);
        // Get 4 direct neighbours
        tmp1 = hm_get(1, i);
        tmp3 = hm_get(0, i - 1);
        tmp4 = hm_get(0, i + 1);
        sum = 
        // Get 4 diaginal neighbours
        tmp1 = hm_get(1, i - 1);
        tmp2 = hm_get(1, i + 1);
   }

   //  Bottom Row
   for ( i = 1; i < column - 1; i++)
   {
        coeff = hm_get_coeff(0,i);
        // Get 4 direct neighbours
        tmp1 = hm_get(row - 2, i);
        tmp3 = hm_get(row - 1, i - 1);
        tmp4 = hm_get(row - 1, i + 1);
        sum = 
        // Get 4 diaginal neighbours
        tmp1 = hm_get(row - 2, i - 1);
        tmp2 = hm_get(row - 2, i + 1);
   }
   // Calculate the boundary cases, col wise
   //  Top Left corner
    // Get 4 direct neighbours
    coeff = hm_get_coeff(0, 0);
    tmp1 = hm_get(0, column - 1);
    tmp2 = hm_get(0, 1);
    tmp3 = hm_get(1, 0);
    sum = 
    // Get 4 diaginal neighbours
    tmp1 = hm_get(1, column - 1);
    tmp2 = hm_get(1, 1);

   //  Bottom Left corner
    // Get 4 direct neighbours
    coeff = hm_get_coeff(row - 1, 0);
    tmp1 = hm_get(row - 1, column - 1);
    tmp2 = hm_get(row - 1, 1);
    tmp3 = hm_get(row - 2, 0);
    sum = 
    // Get 4 diaginal neighbours
    tmp1 = hm_get(row - 2, column - 1);
    tmp2 = hm_get(row - 1, 1);
   //  Left side
   for ( i = 1; i < row - 1; i++)
   {
        coeff = hm_get_coeff(i, 0);
        tmp1 = hm_get(i - 1, 0);
        tmp2 = hm_get(i + 1, 0);
        tmp3 = hm_get(i, column - 1);
        tmp4 = hm_get(i, 1);
        sum = 
        tmp1 = hm_get(i - 1, 1);
        tmp2 = hm_get(i - 1, column - 1);
        tmp3 = hm_get(i + 1, 1);
        tmp4 = hm_get(i + 1, column - 1);
   }

   //  right side
   for ( i = 1; i < row - 1; i++)
   {
        coeff = hm_get_coeff(i, column - 1);
        tmp1 = hm_get(i - 1, column - 1);
        tmp2 = hm_get(i + 1, column - 1);
        tmp3 = hm_get(i, column - 2);
        tmp4 = hm_get(i, 0);
        sum = 
        tmp1 = hm_get(i - 1, column - 2);
        tmp2 = hm_get(i - 1, 0);
        tmp3 = hm_get(i + 1, 0);
        tmp4 = hm_get(i + 1, column - 2);
   }
   //  Top right corner
    // Get 4 direct neighbours
    coeff = hm_get_coeff(0, column - 1);
    tmp1 = hm_get(0, column - 2);
    tmp2 = hm_get(0, 0);
    tmp3 = hm_get(1, column - 1);
    sum = 
    // Get 4 diaginal neighbours
    tmp1 = hm_get(1, column - 2);
    tmp2 = hm_get(1, 0);

   //  Bottom right corner
    // Get 4 direct neighbours
    coeff = hm_get_coeff(row - 1, column - 1);
    tmp1 = hm_get(row - 2, column - 1);
    tmp2 = hm_get(row - 2, column -2);
    tmp3 = hm_get(row - 1, 0);
    sum = 
    // Get 4 diaginal neighbours
    tmp1 = hm_get(row - 2, column - 2);
    tmp2 = hm_get(row - 2, 0);
}


/**
 * hm_init_map
 *  returns 0 on succes
 */
int hm_init_map(struct parameters * p) {
    //Set vars
    size_t row = p->N;
    size_t column = p->M;

    //Init tmps
    void * tmp1 = 0;
    void * tmp2 = 0;
    void * tmp3 = 0;

    //Allocate alligned memory
    int result1 = posix_memalign(&tmp1, 128/8, sizeof(double)*row*column);
    int result2 = posix_memalign(&tmp2, 128/8, sizeof(double)*row*column);
    int result3 = posix_memalign(&tmp3, 128/8, sizeof(double)*row*column);

    if (result1 || result2 || result2) {
        hm_data = 0;
        if (HM_VERBOSE) printf("Error while allocating memory in hm_init_map: %i, %i, %i", result1, result2, result3);
    }else {
        hm_data     =   (double*) tmp1;
        hm_prevData =   (double*) tmp2;
        hm_coeff    =   (double*) tmp3;
        hm_row     =   row;
        hm_column  =   column;

        memcpy((void*)hm_prevData, (void*)p->tinit, sizeof(double)*row*column);
        memcpy((void*)hm_coeff, (void*)p->conductivity, sizeof(double)*row*column);

        
    }

    return result1 | result2 | result3;
}