#include <stdio.h>

#include "mnblas.h"
#include "complexe.h"

#include "flop.h"

#define VECSIZE    65536

#define NB_FOIS    1000

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
    float sdot_nano = 0, ddot_nano = 0, cdotu_nano = 0, zdotu_nano = 0, cdotc_nano = 0, zdotc_nano = 0, sdot_micro = 0, ddot_micro = 0, cdotu_micro = 0, zdotu_micro = 0, cdotc_micro = 0, zdotc_micro = 0;
    void *res = NULL;
    // Test nano de saxpy
    init_flop_nano () ;
    for (int i = 0 ; i < NB_FOIS; i++) {
        vfloat_init (vfloat1, 1.0) ;
        vfloat_init (vfloat2, 2.0) ;

        start_nano = _rdtsc () ;
        mncblas_sdot (VECSIZE, vfloat1, 1, vfloat2, 1) ;
        end_nano = _rdtsc () ;

        sdot_nano+=calcul_flop_nano ("saxpy nano ", 2 * VECSIZE, end_nano-start_nano) ;

    }

    printf ("==========================================================\n") ;

    // Test nano de daxpy
    init_flop_nano () ;
    for (int i = 0 ; i < NB_FOIS; i++) {
        vdouble_init (vdouble1, 1.0) ;
        vdouble_init (vdouble2, 2.0) ;

        start_nano = _rdtsc () ;
        mncblas_ddot (VECSIZE, vdouble1, 1, vdouble2, 1) ;
        end_nano = _rdtsc () ;

        ddot_nano+=calcul_flop_nano ("daxpy nano ", 2 * VECSIZE, end_nano-start_nano) ;

    }

    printf ("==========================================================\n") ;

    // Test nano de caxpy
    init_flop_nano () ;
    for (int i = 0 ; i < NB_FOIS; i++) {
        vcomplexe_float_init (vcomplexe_float1, 1.0, 1.0) ;
        vcomplexe_float_init (vcomplexe_float2, 2.0, 2.0) ;

        start_nano = _rdtsc () ;
        mncblas_cdotu_sub (VECSIZE, vcomplexe_float1, 1, vcomplexe_float2, 1, res) ;
        res = NULL;
        end_nano = _rdtsc () ;

        cdotu_nano+=calcul_flop_nano ("caxpy nano ", 2 * VECSIZE, end_nano-start_nano) ;
    }

//    printf ("==========================================================\n") ;
//
//    // Test nano de zaxpy
//    init_flop_nano () ;
//    for (int i = 0 ; i < NB_FOIS; i++) {
//        vcomplexe_double_init (vcomplexe_double1, 1.0, 1.0) ;
//        vcomplexe_double_init (vcomplexe_double2, 2.0, 2.0) ;
//
//        start_nano = _rdtsc () ;
//        mncblas_zdotu_sub (VECSIZE, vcomplexe_double1, 1, vcomplexe_double2, 1, res) ;
//        end_nano = _rdtsc () ;
//
//        zdotu_nano+=calcul_flop_nano ("zaxpy nano ", 2 * VECSIZE, end_nano-start_nano) ;
//    }
//
//    printf ("==========================================================\n") ;
//    // Test nano de zaxpy
//    init_flop_nano () ;
//    for (int i = 0 ; i < NB_FOIS; i++) {
//        vcomplexe_double_init (vcomplexe_double1, 1.0, 1.0) ;
//        vcomplexe_double_init (vcomplexe_double2, 2.0, 2.0) ;
//
//        start_nano = _rdtsc () ;
//        mncblas_zdotc_sub (VECSIZE,  vcomplexe_double1, 1, vcomplexe_double2, 1, res) ;
//        end_nano = _rdtsc () ;
//
//        zdotc_nano+=calcul_flop_nano ("zaxpy nano ", 2 * VECSIZE, end_nano-start_nano) ;
//    }
//
//    printf ("==========================================================\n") ;
//    // Test nano de zaxpy
//    init_flop_nano () ;
//    for (int i = 0 ; i < NB_FOIS; i++) {
//        vcomplexe_double_init (vcomplexe_double1, 1.0, 1.0) ;
//        vcomplexe_double_init (vcomplexe_double2, 2.0, 2.0) ;
//
//        start_nano = _rdtsc () ;
//        mncblas_cdotc_sub (VECSIZE, vcomplexe_double1, 1, vcomplexe_double2, 1, res) ;
//        end_nano = _rdtsc () ;
//
//        cdotc_nano+=calcul_flop_nano ("zaxpy nano ", 2 * VECSIZE, end_nano-start_nano) ;
//    }
//
//    printf ("==========================================================\n") ;
//
//    // Test micro de saxpy
//    init_flop_micro () ;
//    for (int i = 0 ; i < NB_FOIS; i++) {
//        vfloat_init (vfloat1, 1.0) ;
//        vfloat_init (vfloat2, 2.0) ;
//
//        TOP(start) ;
//        mncblas_sdot (VECSIZE, vfloat1, 1, vfloat2, 1) ;
//        TOP(end) ;
//
//        sdot_micro+=calcul_flop_micro ("saxpy micro", 2 * VECSIZE, tdiff_micro (&start, &end)) ;
//    }
//
//    printf ("==========================================================\n") ;
//
//    // Test micro de daxpy
//    init_flop_micro () ;
//    for (int i = 0 ; i < NB_FOIS; i++) {
//        vdouble_init (vdouble1, 1.0) ;
//        vdouble_init (vdouble2, 2.0) ;
//
//        TOP(start) ;
//        mncblas_ddot (VECSIZE, vdouble1, 1, vdouble2, 1) ;
//        TOP(end) ;
//
//        ddot_micro+=calcul_flop_micro ("daxpy micro", 2 * VECSIZE, tdiff_micro (&start, &end)) ;
//    }
//
//    printf ("==========================================================\n") ;
//
//    // Test micro de caxpy
//    init_flop_micro () ;
//    for (int i = 0 ; i < NB_FOIS; i++) {
//        vcomplexe_float_init (vcomplexe_float1, 1.0, 1.0) ;
//        vcomplexe_float_init (vcomplexe_float2, 2.0, 2.0) ;
//
//        TOP(start) ;
//        mncblas_cdotu_sub (VECSIZE, vcomplexe_float1, 1, vcomplexe_float2, 1, res) ;
//        TOP(end) ;
//
//        cdotu_micro+=calcul_flop_micro ("caxpy micro", 2 * VECSIZE, tdiff_micro (&start, &end)) ;
//    }
//
//    printf ("==========================================================\n") ;
//
//    // Test micro de cdotc_sub
//    init_flop_micro () ;
//    for (int i = 0 ; i < NB_FOIS; i++) {
//        vcomplexe_double_init (vcomplexe_double1, 1.0, 1.0) ;
//        vcomplexe_double_init (vcomplexe_double2, 2.0, 2.0) ;
//
//        TOP(start) ;
//        mncblas_cdotc_sub (VECSIZE, vcomplexe_double1, 1, vcomplexe_double2, 1, res) ;
//        TOP(end) ;
//
//        cdotc_micro+=calcul_flop_micro ("zaxpy micro", 2 * VECSIZE, tdiff_micro (&start, &end)) ;
//    }
//
//    printf ("==========================================================\n") ;
//
//    // Test micro de zaxpy
//    init_flop_micro () ;
//    for (int i = 0 ; i < NB_FOIS; i++) {
//        vcomplexe_double_init (vcomplexe_double1, 1.0, 1.0) ;
//        vcomplexe_double_init (vcomplexe_double2, 2.0, 2.0) ;
//
//        TOP(start) ;
//        mncblas_zdotu_sub (VECSIZE, vcomplexe_double1, 1, vcomplexe_double2, 1, res) ;
//        TOP(end) ;
//
//        zdotu_micro+=calcul_flop_micro ("zaxpy micro", 2 * VECSIZE, tdiff_micro (&start, &end)) ;
//    }
//    printf ("==========================================================\n") ;
//    // Test micro de zaxpy
//    init_flop_micro () ;
//    for (int i = 0 ; i < NB_FOIS; i++) {
//        vcomplexe_double_init (vcomplexe_double1, 1.0, 1.0) ;
//        vcomplexe_double_init (vcomplexe_double2, 2.0, 2.0) ;
//
//        TOP(start) ;
//        mncblas_zdotc_sub (VECSIZE, vcomplexe_double1, 1, vcomplexe_double2, 1, res) ;
//        TOP(end) ;
//
//        zdotc_micro+=calcul_flop_micro ("zaxpy micro", 2 * VECSIZE, tdiff_micro (&start, &end)) ;
//    }
//    printf ("==========================================================\n") ;

    printf("sdot_nano %5.3f\n", sdot_nano / NB_FOIS);
    printf("ddot_nano %5.3f\n", ddot_nano / NB_FOIS);
    printf("cdotu_nano %5.3f\n", cdotu_nano / NB_FOIS);
    printf("cdotc_nano %5.3f\n", cdotc_nano / NB_FOIS);
    printf("zdotu_nano %5.3f\n", zdotu_nano / NB_FOIS);
    printf("zdotc_nano %5.3f\n", zdotc_nano / NB_FOIS);
    printf("sdot_micro %5.3f\n", sdot_micro / NB_FOIS);
    printf("ddot_micro %5.3f\n", ddot_micro / NB_FOIS);
    printf("cdotu_micro %5.3f\n", cdotu_micro / NB_FOIS);
    printf("cdotc_micro %5.3f\n", cdotc_micro / NB_FOIS);
    printf("zdotu_micro %5.3f\n", zdotu_micro / NB_FOIS);
    printf("zdotc_micro %5.3f\n", zdotc_micro / NB_FOIS);
}
