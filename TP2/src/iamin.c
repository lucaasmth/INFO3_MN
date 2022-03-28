#include "../include/mnblas.h"
#include "../include/complexe2.h"

#define abs(a) ((a)>0?(a):(-a))

CBLAS_INDEX mncblas_isamin(const int N, const float  *X, const int incX) {
    if( N < 0 )
        return 0;
    CBLAS_INDEX res = 0;
    float smallest = X[0];
    register unsigned int i = 1 ;

    for (i = 0 ; i < N ; i += incX)
    {
        if(X[i] < smallest) {
            res = i;
            smallest = X[i];
        }
    }
    return res;
}

CBLAS_INDEX mncblas_idamin(const int N, const double *X, const int incX) {
    if( N < 0 )
        return 0;
    CBLAS_INDEX res = 0;
    double smallest = X[0];
    register unsigned int i = 1 ;

    for (i = 0 ; i < N ; i += incX)
    {
        if(X[i] < smallest) {
            res = i;
            smallest = X[i];
        }
    }
    return res;
}

CBLAS_INDEX mncblas_icamin(const int N, const void   *X, const int incX) {
    if( N < 0 )
        return 0;
    CBLAS_INDEX res = 0;
    complexe_float_t first_elem = *((complexe_float_t *) X + 0);
    float smallest = abs(first_elem.imaginary) + abs(first_elem.real);
    register unsigned int i = 1;

    for (i = 0 ; i < N ; i += incX)
    {
        complexe_float_t elem = *((complexe_float_t *) X + i);
        float calc = abs(elem.imaginary) + abs(elem.real);
        if(calc < smallest) {
            res = i;
            smallest = calc;
        }
    }
    return res;
}

CBLAS_INDEX mncblas_izamin(const int N, const void   *X, const int incX) {
    if( N < 0 )
        return 0;
    CBLAS_INDEX res = 0;
    complexe_double_t first_elem = *((complexe_double_t *) X + 0);
    double smallest = abs(first_elem.imaginary) + abs(first_elem.real);
    register unsigned int i = 1;

    for (i = 0 ; i < N ; i += incX)
    {
        complexe_double_t elem = *((complexe_double_t *) X + i);
        double calc = abs(elem.imaginary) + abs(elem.real);
        if(calc < smallest) {
            res = i;
            smallest = calc;
        }
    }
    return res;
}