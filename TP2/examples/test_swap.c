#include <stdio.h>

#include "mnblas.h"
#include "complexe.h"

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

int main (int argc, char **argv) {
//	struct timeval start, end ;
    unsigned long long int start_nano, end_nano;

    // Test nano de sswap
    init_flop_nano();
    for (int i = 0; i < NB_FOIS; i++) {
        vfloat_init(vfloat1, 1.0);
        vfloat_init(vfloat2, 2.0);

        start_nano = _rdtsc ();
        mncblas_sswap(VECSIZE, vfloat1, 1, vfloat1, 1);
        end_nano = _rdtsc ();

        calcul_byte("sswap", VECSIZE * sizeof(float) * 3, end_nano - start_nano);
    }
    printf("==========================================================\n");

    // Test nano de dswap
    init_flop_nano();
    for (int i = 0; i < NB_FOIS; i++) {
        vdouble_init(vdouble1, 1.0);
        vdouble_init(vdouble2, 2.0);

        start_nano = _rdtsc ();
        mncblas_dswap(VECSIZE, vdouble1, 1, vdouble2, 1);
        end_nano = _rdtsc ();

        calcul_byte("dswap", VECSIZE * sizeof(double) * 3, end_nano - start_nano);

    }
    printf("==========================================================\n");

    // Test nano de cswap
    init_flop_nano();
    for (int i = 0; i < NB_FOIS; i++) {
        vcomplexe_float_init(vcomplexe_float1, 1.0, 1.0);
        vcomplexe_float_init(vcomplexe_float2, 2.0, 2.0);

        start_nano = _rdtsc ();
        mncblas_cswap(VECSIZE, vcomplexe_float1, 1, vcomplexe_float2, 1);
        end_nano = _rdtsc ();

        calcul_byte("cswap", VECSIZE * sizeof(float) * 6, end_nano - start_nano);

    }
    printf("==========================================================\n");

    // Test nano de zswap
    init_flop_nano();
    for (int i = 0; i < NB_FOIS; i++) {
        vcomplexe_double_init(vcomplexe_double1, 1.0, 1.0);
        vcomplexe_double_init(vcomplexe_double2, 2.0, 2.0);

        start_nano = _rdtsc ();
        mncblas_zswap(VECSIZE, vcomplexe_double1, 1, vcomplexe_double2, 1);
        end_nano = _rdtsc ();

        calcul_byte("zswap", VECSIZE * sizeof(double) * 6, end_nano - start_nano);

    }
}
