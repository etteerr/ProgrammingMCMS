#include <sys/time.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include "../src/compute.h"

#define sideCoef 0.58578643762

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
	return hm_prevData[hm_column * row + column];
}

inline double hm_get_current(size_t row, size_t column) {
	return hm_data[hm_column * row + column];
}

inline void hm_set(size_t row, size_t column, double val) {
	hm_data[hm_column * row + column] = val;
}

inline double hm_get_coeff(size_t row, size_t column) {
	return hm_coeff[hm_column * row + column];
}

inline void hm_set_coeff(size_t row, size_t column, double val) {
	hm_coeff[hm_column * row + column] = val;
}

inline void hm_swap() {
	double * tmp = hm_data;
	hm_data = hm_prevData;
	hm_prevData = tmp;
}

inline void do_cell(unsigned int i, unsigned int j) {


}

void do_compute(const struct parameters* p, struct results *r)
{
	unsigned int i, j;
	double cdiag, cnext;

	// Calculate the inner sums
	for ( i = 1; i < hm_row - 1; i++)
	{
		for ( j = 1; j < hm_column - 1; j++)
		{
			cnext = (1-hm_get_coeff(i,j))*sideCoef;
			cdiag = 1-hm_get_coeff(i,j)-cnext;
			hm_set(i,j,
					//SUm neighbours (direct)
					((hm_get(i-1,j) + hm_get(i+1,j) + hm_get(i,j-1) + hm_get(i,j+1))/4.0)*cnext +
					//Sum diag
					((hm_get(i-1,j-1) + hm_get(i+1,j+1) + hm_get(i+1,j-1) + hm_get(i-1,j+1))/4.0)*cdiag +
					//Add current
					hm_get(i,j) * hm_get_coeff(i,j)
			);
//			printf("%.5f -> %.5f\n", hm_get(i,j), hm_get_current(i,j));fflush(stdout);
//			printf("%.5f -> %.5f\n", cdiag, cnext);fflush(stdout);
//			printf("%.16f\n", ((hm_get(i-1,j) + hm_get(i+1,j) + hm_get(i,j-1) + hm_get(i,j+1))/4.0)*cnext);
		}
	}

	//left side
	for ( i = 1; i < hm_row - 1; i++)
	{

		cnext = (1-hm_get_coeff(i,j))*sideCoef;
		cdiag = 1-hm_get_coeff(i,j)-cnext;
		hm_set(i,j,
				//SUm neighbours (direct)
				((hm_get(i-1,j) + hm_get(i+1,j) + hm_get(i,hm_column-1) + hm_get(i,j+1))/4.0)*cnext +
				//Sum diag
				((hm_get(i-1,hm_column-1) + hm_get(i+1,j+1) + hm_get(i+1,hm_column-1) + hm_get(i-1,j+1))/4.0)*cdiag +
				//Add current
				hm_get(i,j) * hm_get_coeff(i,j)
		);
	}

	//  right side
	j=hm_column-1;
	for ( i = 1; i < hm_row - 1; i++)
	{
		cnext = (1-hm_get_coeff(i,j))*sideCoef;
		cdiag = 1-hm_get_coeff(i,j)-cnext;
		hm_set(i,j,
				//SUm neighbours (direct)
				((hm_get(i-1,j) + hm_get(i+1,j) + hm_get(i,j-1) + hm_get(i,0))/4.0)*cnext +
				//Sum diag
				((hm_get(i-1,j-1) + hm_get(i+1,0) + hm_get(i+1,j-1) + hm_get(i-1,0))/4.0)*cdiag +
				//Add current
				hm_get(i,j) * hm_get_coeff(i,j)
		);
	}



	// add inter
	r->niter++;

	//convergence
	r->maxdiff = 0;
	r->tmax = 0;
	r->tmin = p->io_tmax;
	r->tavg = 0;
	for(i = hm_column; i < hm_row * hm_column-hm_column; i++) {
		//max
		if (hm_data[i] > r->tmax) r->tmax = hm_data[i];
		//min
		if (hm_data[i] < r->tmin) r->tmin = hm_data[i];
		//maxdiff
		if (abs(hm_data[i]-hm_prevData[i]) > r->maxdiff) r->maxdiff = abs(hm_data[i]-hm_prevData[i]);
		//searchme
		r->tavg += hm_data[i]/(hm_column*hm_row);
	}

	report_results(p,r);

	hm_swap();

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
	int result1 = posix_memalign(&tmp1, 128/8, sizeof(double)*row*column+2*sizeof(double)*column);
	int result2 = posix_memalign(&tmp2, 128/8, sizeof(double)*row*column+2*sizeof(double)*column);
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

		memcpy((void*)(hm_prevData+column), (void*)p->tinit, sizeof(double)*row*column);
		memcpy((void*)hm_coeff, (void*)p->conductivity, sizeof(double)*row*column);

		//Ghost rows
		memcpy((void*)hm_prevData, (void*)p->tinit, sizeof(double)*column);
		memcpy((void*)(hm_prevData+(row+1)*column), (void*)(p->tinit+row*column-column), sizeof(double)*column);

		//fix row
		hm_row--;

		//print test
		begin_picture(666, column, row, p->io_tmin, p->io_tmax);
			for (int i = 0; i < row; i++)
				for(int j = 0; j < column; j++)
					draw_point(i,j,p->tinit[i*column + j]);

		end_picture();
		begin_picture(667, column, row, p->io_tmin, p->io_tmax);
			for (int i = 0; i < row; i++)
				for(int j = 0; j < column; j++)
					draw_point(i,j,hm_get(i,j));

		end_picture();

	}

	return result1 | result2 | result3;
}
