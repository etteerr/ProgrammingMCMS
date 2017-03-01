#include <math.h>
#include <stdlib.h>
#include <stdio.h>
//#include <emmintrin.h>
#include <x86intrin.h>
#include <unistd.h>
#include "../src/compute.h"
#include "../src/output.h"

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


void do_compute(const struct parameters* p, struct results *r)
{
	//Intresting way of declaring privates
	double cdiag, cnext;
	unsigned int i, j;
	// Calculate the inner sums

	for ( i = 0; i < hm_row ; i++)
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
//		#pragma omp barrier


	r->niter++;



	hm_swap();



}

void do_compute_sse(const struct parameters* p, struct results *r)
{
		//Intresting way of declaring privates
		double cdiag, cnext;
		unsigned int i, j;
		// Calculate the inner sums
		for ( i = 0; i < hm_row ; i++)
		{
			double diag, next;
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
			for ( j = 1; j < hm_column-1; j+=2)
			{
				//Calculate both cnext and cdiag for i,j and i,j+1
				__m128d cnextd, cdiagd;
//				cnext = (1-hm_get_coeff(i,j))*sideCoef;
//				cdiag = 1-hm_get_coeff(i,j)-cnext;
				cdiagd = cnextd = 1 - _mm_loadu_pd(hm_get_coeffa(i,j));
				cnextd *= sideCoef;
				cdiagd -= cnextd;
				cnextd *= .25;
				cdiagd *= .25;

				//Assigned calculated cnext & cdiag
				cnext = (double)cnextd[0];
				cdiag = (double)cdiagd[0];

				//i,j
				__m128d top2, top1, mid2, mid1, bot2, bot1, diag1, diag2, next1, next2; //10/16 vars
//				╔════════════════════════════════════════╗
//				║                Load (1)                ║
//				╠═════════╦═════════╦═════════╦══════════╣
//				║ top2[0] ║ top2[1] ║ top1[0] ║ top1[1]  ║
//				╠═════════╬═════════╬═════════╬══════════╣
//				║ mid2[0] ║ mid2[1] ║ mid1[0] ║ mid1[1]  ║
//				╠═════════╬═════════╬═════════╬══════════╣
//				║ bot2[0] ║ bot2[1] ║ bot1[0] ║ bot1[1]  ║
//				╠═════════╩═════════╩═════════╩══════════╣

				//Load stage 1
				top2 = _mm_load_pd(hm_geta(i-1,j-1));//upper  left + 1
				top1 = _mm_load_pd(hm_geta(i-1,j+1));//upper  right + 1
				mid2 = _mm_load_pd(hm_geta(i  ,j-1));//middle left + middle
				mid1 = _mm_load_pd(hm_geta(i  ,j+1));//middle right + 1
				bot2 = _mm_load_pd(hm_geta(i+1,j-1));//bot    left + 1
				bot1 = _mm_load_pd(hm_geta(i+1,j+1));//bot    right + 1

				//Sort
				//Note: _MM_SHUFFLE2(SecondVar, firstVar)
				diag1 = _mm_shuffle_pd(top2, top1, _MM_SHUFFLE2(0,0)); //First of top2 and first of top1
				diag2 = _mm_shuffle_pd(bot2, bot1, _MM_SHUFFLE2(0,0)); //First of bot2 and first of bot1
				next2 = _mm_shuffle_pd(mid2, mid1, _MM_SHUFFLE2(0,0)); //First and first
				next1 = _mm_shuffle_pd(top2, bot2, _MM_SHUFFLE2(1,1)); //second of top2 and bot2

				//Calculate weighted
				diag2 += diag1;
				diag = (double)diag2[0] + (double)diag2[1];
				next2 += next1;
				next = (double)next2[0] + (double)next2[1];

				next *= cnext;
				diag *= cdiag;

				//Sum
				// hm_set instead of _mm_storeu, as this method is alot faster (almost x2)
				hm_set(i,j,next + diag + hm_get(i,j)*hm_get_coeff(i,j));

				//checkAnswer(i,j,cnext, cdiag);

				//([-1 (0][1)2])[3 4]

				//Assigned calculated cnext & cdiag
				cnext = (double)cnextd[1];
				cdiag = (double)cdiagd[1];

				//Sort
				diag1 = _mm_shuffle_pd(top1, top2, _MM_SHUFFLE2(1,1));
				diag2 = _mm_shuffle_pd(bot1, bot2, _MM_SHUFFLE2(1,1));
				next1 = _mm_shuffle_pd(mid1, mid2, _MM_SHUFFLE2(1,1));
				next2 = _mm_shuffle_pd(top1, bot1, _MM_SHUFFLE2(0,0));


				//Calculate weighted
				diag2 += diag1;
				diag = (double)diag2[0] + (double)diag2[1];
				next2 += next1;
				next = (double)next2[0] + (double)next2[1];

				next *= cnext;
				diag *= cdiag;

				//Sum
				hm_set(i,j+1,next + diag + hm_get(i,j+1) * hm_get_coeff(i,j+1));
				//checkAnswer(i,j+1, cnext, cdiag);

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

		r->niter++;
		hm_swap();

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
