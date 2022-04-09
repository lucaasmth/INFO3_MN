#include <stdio.h>

#include "../include/mnblas.h"
#include "../include/complexe.h"

#include "flop.h"

#define VECSIZE    65536

#define NB_FOIS    10

typedef float vfloat [VECSIZE] ;
typedef double vdouble [VECSIZE] ;
typedef complexe_float_t vcomplexe_float [VECSIZE] ;
typedef complexe_double_t vcomplexe_double [VECSIZE] ;

vfloat vfloat1, vfloat2 ;
vdouble vdouble1, vdouble2 ;
vcomplexe_float vcomplexe_float1, vcomplexe_float2 ;
vcomplexe_double vcomplexe_double1, vcomplexe_double2 ;

void vfloat_init (vfloat V, float x)
{
	register unsigned int i ;

	for (i = 0; i < VECSIZE; i++)
		V [i] = x ;

	return ;
}

void vdouble_init (vdouble V, double x)
{
	register unsigned int i ;

	for (i = 0; i < VECSIZE; i++)
		V [i] = x ;

	return ;
}

void vcomplexe_float_init (vcomplexe_float V, float x, float y)
{
	register unsigned int i ;

	for (i = 0; i < VECSIZE; i++) {
		V [i].real = x ;
		V [i].imaginary = y ;
	}

	return ;
}

void vcomplexe_double_init (vcomplexe_double V, float x, float y)
{
	register unsigned int i ;

	for (i = 0; i < VECSIZE; i++) {
		V [i].real = x ;
		V [i].imaginary = y ;
	}

	return ;
}

void vector_print (vfloat V)
{
	register unsigned int i ;

	for (i = 0; i < VECSIZE; i++)
		printf ("%f ", V[i]) ;
	printf ("\n") ;
	
	return ;
}

int main (int argc, char **argv)
{
	struct timeval start, end ;
	unsigned long long int start_nano, end_nano ;

	// Test nano de sgemm
	init_flop_nano () ;
	for (int i = 0 ; i < NB_FOIS; i++) {
		vfloat_init (vfloat1, 1.0) ;
		vfloat_init (vfloat2, 2.0) ;

		start_nano = _rdtsc () ;
        mncblas_sgemm (MNCblasRowMajor, MNCblasNoTrans, MNCblasNoTrans, 10, 10, 10, 3.0, vfloat1, 10, vfloat2, 10, 1.0, vfloat1, 10) ;
		end_nano = _rdtsc () ;

		calcul_flop_nano ("sgemm nano ", 2 * VECSIZE, end_nano-start_nano) ;
	}
	printf ("==========================================================\n") ;

	// Test nano de dgemm
	init_flop_nano () ;
	for (int i = 0 ; i < NB_FOIS; i++) {
		vdouble_init (vdouble1, 1.0) ;
		vdouble_init (vdouble2, 2.0) ;

		start_nano = _rdtsc () ;
        mncblas_dgemm (MNCblasRowMajor, MNCblasNoTrans, MNCblasNoTrans, 10, 10, 10, 3.0, vdouble1, 10, vdouble2, 10, 1.0, vdouble1, 10) ;
		end_nano = _rdtsc () ;

		calcul_flop_nano ("dgemm nano ", 2 * VECSIZE, end_nano-start_nano) ;
	}
	printf ("==========================================================\n") ;

	// Test nano de cgemm
	init_flop_nano () ;
	for (int i = 0 ; i < NB_FOIS; i++) {
		vcomplexe_float_init (vcomplexe_float1, 1.0, 1.0) ;
		vcomplexe_float_init (vcomplexe_float2, 2.0, 2.0) ;

		start_nano = _rdtsc () ;
        mncblas_cgemm (MNCblasRowMajor, MNCblasNoTrans, MNCblasNoTrans, 10, 10, 10, vcomplexe_float1, vcomplexe_float1, 10, vcomplexe_float2, 10, vcomplexe_float2, vcomplexe_float1, 10) ;
		end_nano = _rdtsc () ;

		calcul_flop_nano ("cgemm nano ", 2 * VECSIZE, end_nano-start_nano) ;
	}
	printf ("==========================================================\n") ;

	// Test nano de zgemm
	init_flop_nano () ;
	for (int i = 0 ; i < NB_FOIS; i++) {
		vcomplexe_double_init (vcomplexe_double1, 1.0, 1.0) ;
		vcomplexe_double_init (vcomplexe_double2, 2.0, 2.0) ;

		start_nano = _rdtsc () ;
        mncblas_zgemm (MNCblasRowMajor, MNCblasNoTrans, MNCblasNoTrans, 10, 10, 10, vcomplexe_double1, vcomplexe_double1, 10, vcomplexe_double2, 10, vcomplexe_double2, vcomplexe_double1, 10) ;
		end_nano = _rdtsc () ;

		calcul_flop_nano ("zgemm nano ", 2 * VECSIZE, end_nano-start_nano) ;
	}
	printf ("==========================================================\n") ;
 
	// Test micro de sgemm
	init_flop_micro () ;
	for (int i = 0 ; i < NB_FOIS; i++) {
		vfloat_init (vfloat1, 1.0) ;
		vfloat_init (vfloat2, 2.0) ;

		TOP(start) ;
		mncblas_sgemm (MNCblasRowMajor, MNCblasNoTrans, MNCblasNoTrans, 10, 10, 10, 3.0, vfloat1, 10, vfloat2, 10, 1.0, vfloat1, 10) ;
		TOP(end) ;

		calcul_flop_micro ("sgemm micro", 2 * VECSIZE, tdiff_micro (&start, &end)) ;
	}
 	printf ("==========================================================\n") ;

	// Test micro de dgemm
	init_flop_micro () ;
	for (int i = 0 ; i < NB_FOIS; i++) {
		vdouble_init (vdouble1, 1.0) ;
		vdouble_init (vdouble2, 2.0) ;

		TOP(start) ;
        mncblas_dgemm (MNCblasRowMajor, MNCblasNoTrans, MNCblasNoTrans, 10, 10, 10, 3.0, vdouble1, 10, vdouble2, 10, 1.0, vdouble1, 10) ;
		TOP(end) ;

		calcul_flop_micro ("dgemm micro", 2 * VECSIZE, tdiff_micro (&start, &end)) ;
	}
 	printf ("==========================================================\n") ;

	// Test micro de cgemm
	init_flop_micro () ;
	for (int i = 0 ; i < NB_FOIS; i++) {
		vcomplexe_float_init (vcomplexe_float1, 1.0, 1.0) ;
		vcomplexe_float_init (vcomplexe_float2, 2.0, 2.0) ;

		TOP(start) ;
        mncblas_cgemm (MNCblasRowMajor, MNCblasNoTrans, MNCblasNoTrans, 10, 10, 10, vcomplexe_float1, vcomplexe_float1, 10, vcomplexe_float2, 10, vcomplexe_float2, vcomplexe_float1, 10) ;
		TOP(end) ;

		calcul_flop_micro ("cgemm micro", 2 * VECSIZE, tdiff_micro (&start, &end)) ;
	}
 	printf ("==========================================================\n") ;

	// Test micro de zgemm
	init_flop_micro () ;
	for (int i = 0 ; i < NB_FOIS; i++) {
		vcomplexe_double_init (vcomplexe_double1, 1.0, 1.0) ;
		vcomplexe_double_init (vcomplexe_double2, 2.0, 2.0) ;

		TOP(start) ;
        mncblas_zgemm (MNCblasRowMajor, MNCblasNoTrans, MNCblasNoTrans, 10, 10, 10, vcomplexe_double1, vcomplexe_double1, 10, vcomplexe_double2, 10, vcomplexe_double2, vcomplexe_double1, 10) ;
		TOP(end) ;

		calcul_flop_micro ("zgemm micro", 2 * VECSIZE, tdiff_micro (&start, &end)) ;
	}
 	printf ("==========================================================\n") ;
}
