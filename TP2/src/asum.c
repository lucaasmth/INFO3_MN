#include "../include/mnblas.h"
#include "../include/complexe2.h"

#define abs(a) ((a)>0?(a):(-a))


float mncblas_sasum (const int N, const float *X, const int incX)
{
    register float res = 0;
    register unsigned int i = 0 ;

    for (i = 0 ; i < N ; i += incX)
    {
        res += abs(X[i]);
    }
    return res;
}

double mncblas_dasum (const int N, const double *X, const int incX)
{
    register double res = 0;
    register unsigned int i = 0 ;

    for (i = 0 ; i < N ; i += incX)
    {
        res += abs(X[i]);
    }
    return res;
}

float mncblas_scasum (const int N, const void *X, const int incX)
{
    register float res = 0;
    register unsigned int i = 0 ;

    for (i = 0 ; i < N ; i += incX){
        complexe_float_t elem = *((complexe_float_t *) X+i);
        res += abs(elem.imaginary);
        res += abs(elem.real);
    }
    return res;
}


double mncblas_dzasum (const int N, const void *X, const int incX)
{
    register double res = 0;
    register unsigned int i = 0 ;

    for (i = 0 ; i < N ; i += incX){
        complexe_double_t elem = *((complexe_double_t *) X+i);
        res += abs(elem.imaginary);
        res += abs(elem.real);
    }
    return res;
}
