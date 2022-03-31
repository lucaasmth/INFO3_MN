#include "../include/mnblas.h"
#include "../include/complexe2.h"

#define abs(a) ((a)>0?(a):(-a))

CBLAS_INDEX mncblas_isamax(const int N, const float  *X, const int incX) {
    if( N < 0 )
        return 0;
    CBLAS_INDEX res = 0;
    float largest = -1;
    register unsigned int i = 1 ;

    for (i = 0 ; i < N ; i += incX)
    {
        if(X[i] > largest) {
            res = i;
            largest = X[i];
        }
    }
    return res;
}

CBLAS_INDEX mncblas_idamax(const int N, const double *X, const int incX) {
    if( N < 0 )
        return 0;
    CBLAS_INDEX res = 0;
    double largest = X[0];
    register unsigned int i = 1 ;

    for (i = 0 ; i < N ; i += incX)
    {
        if(X[i] > largest) {
            res = i;
            largest = X[i];
        }
    }
    return res;
}

CBLAS_INDEX mncblas_icamax(const int N, const void   *X, const int incX) {
    if( N < 0 )
        return 0;
    CBLAS_INDEX res = 0;
    float largest = -1;
    register unsigned int i = 1;

    for (i = 0 ; i < N ; i += incX)
    {
        complexe_float_t elem = *((complexe_float_t *) X + i);
        float calc = abs(elem.imaginary) + abs(elem.real);
        if(calc > largest) {
            res = i;
            largest = calc;
        }
    }
    return res;
}

CBLAS_INDEX mncblas_izamax(const int N, const void   *X, const int incX) {
    if( N < 0 )
        return 0;
    CBLAS_INDEX res = 0;
    double largest = -1;
    register unsigned int i = 1;

    for (i = 0 ; i < N ; i += incX)
    {
        complexe_double_t elem = *((complexe_double_t *) X + i);
        double calc = abs(elem.imaginary) + abs(elem.real);
        if(calc > largest) {
            res = i;
            largest = calc;
        }
    }
    return res;
}