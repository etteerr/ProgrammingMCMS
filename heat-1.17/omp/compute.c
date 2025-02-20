#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include "../src/compute.h"
#include "../src/output.h"

#define checkSanity 0
#define sanityPGM 0
#define sideCoef 0.58578643762
#define _E_parallel 1
#define ompsection parallel default(none) \
	shared(hm_data, hm_prevData, hm_coeff, hm_row, hm_column, r) \
	if((hm_row*hm_column > 1) && _E_parallel)

#define ompforloop2d for nowait schedule(static)
//#define ompforloop1d for nowait schedule(guided)

//Static pointers to data
double * hm_data = 0;
double * hm_prevData = 0;
double * hm_coeff = 0;

//Fixed global dimentions
size_t hm_row = 0;
size_t hm_column = 0;


//Getters and setters
// Inline for performance
#define hm_get_l(R,C) hm_prevData[hm_column * (R+1) + C]
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

inline void do_cell(unsigned int i, unsigned int j) {


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
				if (fabs(hm_get_current(i,j)-hm_get_l(i,j)) > r->maxdiff) r->maxdiff = fabs(hm_get_current(i,j)-hm_get_l(i,j));
				//searchme
				sum += hm_get_current(i,j);
    		}

    	r->tavg += sum/(hm_column*hm_row);
}

void do_compute(const struct parameters* p, struct results *r)
{
	#pragma omp ompsection
	{
		//Intresting way of declaring privates
		double cdiag, cnext;
		unsigned int i, j;
		// Calculate the inner sums
		#pragma omp ompforloop2d
		for ( i = 0; i < hm_row ; i++)
		{

			//left
			cnext = (1-hm_get_coeff(i,0))*sideCoef;
			cdiag = 1-hm_get_coeff(i,0)-cnext;
			hm_set(i,0,
					//SUm neighbours (direct)
					((hm_get_l(i-1,0) + hm_get_l(i+1,0) + hm_get_l(i,hm_column-1) + hm_get_l(i,1))*.25)*cnext +
					//Sum diag
					((hm_get_l(i-1,hm_column-1) + hm_get_l(i+1,1) + hm_get_l(i+1,hm_column-1) + hm_get_l(i-1,1))*.25)*cdiag +
					//Add current
					hm_get_l(i,0) * hm_get_coeff(i,0)
			);

			//Middle
			for ( j = 1; j < hm_column-1; j++)
			{
				cnext = (1-hm_get_coeff(i,j))*sideCoef;
				cdiag = 1-hm_get_coeff(i,j)-cnext;
				hm_set(i,j,
						//SUm neighbours (direct)
						((hm_get_l(i-1,j) + hm_get_l(i+1,j) + hm_get_l(i,j-1) + hm_get_l(i,j+1))*.25)*cnext +
						//Sum diag
						((hm_get_l(i-1,j-1) + hm_get_l(i+1,j+1) + hm_get_l(i+1,j-1) + hm_get_l(i-1,j+1))*.25)*cdiag +
						//Add current
						hm_get_l(i,j) * hm_get_coeff(i,j)
				);
			}

			//right
			cnext = (1-hm_get_coeff(i,(hm_column-1)))*sideCoef;
			cdiag = 1-hm_get_coeff(i,(hm_column-1))-cnext;
			hm_set(i,(hm_column-1),
					//SUm neighbours (direct)
					((hm_get_l(i-1,(hm_column-1)) + hm_get_l(i+1,(hm_column-1)) + hm_get_l(i,(hm_column-1)-1) + hm_get_l(i,0))*.25)*cnext +
					//Sum diag
					((hm_get_l(i-1,(hm_column-1)-1) + hm_get_l(i+1,0) + hm_get_l(i+1,(hm_column-1)-1) + hm_get_l(i-1,0))*.25)*cdiag +
					//Add current
					hm_get_l(i,(hm_column-1)) * hm_get_coeff(i,(hm_column-1))
			);
		}
//		#pragma omp barrier

		#pragma omp single
		r->niter++;

		#pragma omp barrier
		#pragma omp single
		hm_swap();

	}//end omp

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
