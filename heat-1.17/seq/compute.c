#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <emmintrin.h>
#include <time.h>
#include "../src/compute.h"
#include "../src/output.h"

#define checkSanity 1
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
__attribute__((assume_aligned(16)))
inline double hm_get(int row, int column) {
	if (checkSanity) {
		if (row > (int)hm_row || column > hm_column || row < -1 || column < 0)
			printf("Invalid coördinate requested (%i,%i)\n", row, column);
	}
	return hm_prevData[hm_column * (row+1) + column];
}

inline double * hm_get_p(int row, int column) {
	return hm_prevData + hm_column * (row+1) + column;
}

__attribute__((assume_aligned(16)))
inline double hm_get_current(size_t  row, size_t  column) {
	if (checkSanity) {
		if (row > hm_row || column > hm_column || row < 0 || column < 0)
			printf("Invalid coördinate requested (current) (%i,%i)\n", row, column);
	}
	return hm_data[hm_column * (row+1) + column];
}

inline void hm_set(size_t row, size_t column, double val) {
	if (checkSanity) {
		if (row > hm_row || column > hm_column || row < 0 || column < 0)
			printf("Invalid coördinate set (%i,%i)\n", row, column);
	}
	*(hm_data + hm_column * (row+1) + column) = val;
}
__attribute__((assume_aligned(16)))
inline double hm_get_coeff(size_t row, size_t column) {
	return hm_coeff[hm_column * row + column];
}

__attribute__((assume_aligned(16)))
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
    	for(int i = 0; i < hm_row; i++)
    		for (int j =0; j < hm_column; j++) {
				//max
				if (hm_get_current(i,j) > r->tmax) r->tmax = hm_get_current(i,j);
				//min
				if (hm_get_current(i,j) < r->tmin) r->tmin = hm_get_current(i,j);
				//maxdiff
				if (abs(hm_get_current(i,j)-hm_get(i,j)) > r->maxdiff) r->maxdiff = abs(hm_get_current(i,j)-hm_get(i,j));
				//searchme
				r->tavg += hm_get_current(i,j)/(hm_column*hm_row);
    		}
}
const int a = _MM_SHUFFLE2(1,0); //Shuffle operator 'a' [0 (0][1) 1]=> [0 1]
const int b = _MM_SHUFFLE2(0,0); //Shuffle operator 'b' [0 (0][1) 1]=> [1 0]
#define _ESWAPA(A,B,OP1, OP2){ \
	__v2df __ESWAPA_TMP; \
	__ESWAPA_TMP = _mm_shuffle_pd(A,B,OP1); \
	B   = _mm_shuffle_pd(A,B,OP2); \
	A   = __ESWAPA_TMP; \
	}

void do_compute(const struct parameters* p, struct results *r)
{
	size_t i, j;
	double cdiag, cnext;
	double * p2v2; //pointer 2 vector 2
	__v2df top2, bot2, mid2, mid1, top1, bot1, next1, res; // 1+ 6

	// Calculate the inner sums
	for ( i = 0; i < hm_row ; i++)
	{
		for ( j = 1; j < hm_column-1-1; j+=2) //Extra -1 for prevent j overflow in heatmap
		{
			//Init
			next1 = _mm_set1_pd(0.0);

			//Load + shuffle top
			top2 = _mm_load_pd(hm_get_p(i-1, j-1));
			top1 = _mm_load_pd(hm_get_p(i-1, j+1));
			_ESWAPA(top2, top1, a, b); //top2 is diag top, top1 is middle

			//Load + shuffle bot
			bot2 = _mm_load_pd(hm_get_p(i+1, j-1));
			bot1 = _mm_load_pd(hm_get_p(i+1, j+1));
			printf("%.5f %.5f\n", ((double *)&bot2)[0], ((double *)&bot2)[1]);
			printf("%.5f %.5f\n", ((double *)&bot1)[0], ((double *)&bot1)[1]);
			_ESWAPA(bot2, bot1, a, b); //bot2 is diag bot, bot1 is middle
			printf("%.5f %.5f\n", ((double *)&bot2)[0], ((double *)&bot2)[1]);
			printf("%.5f %.5f\n", ((double *)&bot1)[0], ((double *)&bot1)[1]);
			//shuffle  bot1&top1
			_ESWAPA(top1, bot1, a, b); //Shuffles 0 x 1 x to 01xx, bot1 now next side info

			//Load  + shuffle mid
			mid2 = _mm_load_pd(hm_get_p(i+0, j-1));
			mid1 = _mm_load_pd(hm_get_p(i+0, j+1));
			_ESWAPA(mid2, mid1, a, b); //mid2 is next mid, mid1 is middle (+ mid next)

			//Shuffle

			_ESWAPA(mid1, next1,a,b); //replace mid1[1] with 0 and set next[0] to the second
			//((double *)&mid1)[1] = 0.0;

			//Calc modifiers
			cnext = (1-hm_get_coeff(i,j))*sideCoef;
			cdiag = 1-hm_get_coeff(i,j)-cnext;
			cnext /= 4;
			cdiag /= 4;

			//Calc with weight and sum
			// Note that filling xmmX is left to compiler  (*= instead of _mm_mul_pd)
			res = (top2 * cdiag) + (bot2 * cdiag) + (mid2 * cnext) + (top1 * cnext) + (mid1 * hm_get_coeff(i,j)); //mid1[1] == 0

			//sum final
			p2v2 = (double*)&res; //reinterpret as double[2]
			p2v2[0] += p2v2[1];    //Add double[1] to double[0]

			//Save
			hm_set(i,j,p2v2[0]); //Store double[0]

			//Next cell

			//Shuffle
			_ESWAPA(mid1, next1, a, b);

			//Calc modifiers
			cnext = (1-hm_get_coeff(i,j+1))*sideCoef;
			cdiag = 1-hm_get_coeff(i,j+1)-cnext;
			cnext /= 4;
			cdiag /= 4;

			//calculate
			top2 *= cnext;
			mid1 *= cnext;
			top1 *= cdiag;
			bot1 *= cdiag;
			mid2 *= hm_get_coeff(i,j+1);

			//sum
			top2 += mid1 + top2 + bot1;
			p2v2 = (double *)&mid2; //Our mid value is in mid2[1]
			top2 += p2v2[1];
			p2v2 = (double *)&top2; //Result is stored in top2[0:1]
			p2v2[0] += p2v2[1];

			//store
			hm_set(i,j+1, p2v2[0]);

//			hm_set(i,j,
//					//SUm neighbours (direct)
//					((hm_get(i-1,j) + hm_get(i+1,j) + hm_get(i,j-1) + hm_get(i,j+1))/4.0)*cnext +
//					//Sum diag
//					((hm_get(i-1,j-1) + hm_get(i+1,j+1) + hm_get(i+1,j-1) + hm_get(i-1,j+1))/4.0)*cdiag +
//					//Add current
//					hm_get(i,j) * hm_get_coeff(i,j)
//			);

			if (checkSanity) {
				if (hm_get(i,j) != 1.0 || hm_get_current(i,j) != 1.0) {
					printf("(%i,%i) %.5f -> %.5f\n",i,j, hm_get(i,j), hm_get_current(i,j));fflush(stdout);
					printf("(%i,%i) %.5f -> %.5f\n",i,j+1, hm_get(i,j+1), hm_get_current(i,j+1));fflush(stdout);
					sleep(1);
				}
			}
		}

		//Uneven case
		if (hm_column % 2 == 1) {

			//set j
			j = hm_column - 2;

			//Calc modifiers
			cnext = (1-hm_get_coeff(i,j+1))*sideCoef;
			cdiag = 1-hm_get_coeff(i,j+1)-cnext;
			//Calc
			hm_set(i,j,
					//SUm neighbours (direct)
					((hm_get(i-1,j) + hm_get(i+1,j) + hm_get(i,j-1) + hm_get(i,j+1))/4.0)*cnext +
					//Sum diag
					((hm_get(i-1,j-1) + hm_get(i+1,j+1) + hm_get(i+1,j-1) + hm_get(i-1,j+1))/4.0)*cdiag +
					//Add current
					hm_get(i,j) * hm_get_coeff(i,j)
			);
		}
	}


	//left side
	for ( i = 0; i < hm_row; i++)
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

		if (checkSanity) {
			if (hm_get(i,j) != 1.0 || hm_get_current(i,j) != 1.0) {
				printf("(%i,%i) %.5f -> %.5f\n",i,j, hm_get(i,j), hm_get_current(i,j));fflush(stdout);
			}
		}
	}

	//  right side
	j=hm_column-1;
	for ( i = 0; i < hm_row; i++)
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

	if (checkSanity && sanityPGM) renderImageMemory(r->niter, -100.0, +100.0);

	hm_swap();

	// add inter
	r->niter++;

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
		// coeff
		memcpy((void*)hm_coeff, (void*)p->conductivity, sizeof(double)*row*column);

		//Ghost rows
		// prevData
		memcpy((void*)(hm_prevData), (void*)p->tinit, sizeof(double)*column);
		memcpy((void*)(hm_prevData+(row+1)*column), (void*)(p->tinit+row*column-column), sizeof(double)*column);

		// current map
		memcpy((void*)(hm_data), (void*)p->tinit, sizeof(double)*column);
		memcpy((void*)(hm_data+(row+1)*column), (void*)p->tinit, sizeof(double)*column);

		//print test
		begin_picture(666, column, row, p->io_tmin, p->io_tmax);
			for (int i = 0; i < row; i++)
				for(int j = 0; j < column; j++)
					draw_point(i,j,p->tinit[i*column + j]);

		end_picture();
		begin_picture(667, column, row, p->io_tmin, p->io_tmax);
			for (int i = 1; i < row+1; i++)
				for(int j = 0; j < column; j++)
					draw_point(i,j,hm_get(i,j));

		end_picture();
		begin_picture(670, column, row, p->io_tmin, p->io_tmax);
			for (int i = 1; i < row+1; i++)
				for(int j = 0; j < column; j++)
					draw_point(i,j,hm_get_coeff(i,j));

		end_picture();

		renderImageMemory(668, p->io_tmin, p->io_tmax);

		if (checkSanity) {
	//		Test map
			for (int i = 0; i < (column*(row+2))*2; i++)
				hm_data[i] = i;// ((double)rand()/RAND_MAX)*(double)rand();
		}
	}//end of allocate succes if
	return result1 | result2;
}
