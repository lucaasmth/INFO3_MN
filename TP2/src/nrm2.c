#include "../include/mnblas.h"
#include "../include/complexe2.h"
#include <math.h>

float mnblas_snrm2(const int N, const float *X, const int incX){
    float res = 0;
    register unsigned int i = 0;
    for (i = 0 ; i < N ; i += incX)
    {
        res+= powf(X[i], 2);
    }
    return sqrtf(res);
}

double mnblas_dnrm2(const int N, const double *X, const int incX){
    double res = 0;
    register unsigned int i = 0;
    for (i = 0 ; i < N ; i += incX)
    {
        res+= pow(X[i], 2);
    }
    return sqrt(res);
}

float mnblas_scnrm2(const int N, const void *X, const int incX){
    float res = 0;
    register unsigned int i = 0;
    for (i = 0 ; i < N ; i += incX)
    {
        complexe_float_t elem = *((complexe_float_t *) X + i);
        res += sqrtf(powf(elem.real, 2) + powf(elem.imaginary, 2));
    }
    return sqrtf(res);
}

double mnblas_dznrm2(const int N, const void *X, const int incX){
    double res = 0;
    register unsigned int i = 0;
    for (i = 0 ; i < N ; i += incX)
    {
        complexe_double_t elem = *((complexe_double_t *) X + i);
        res += sqrt(pow(elem.real, 2) + pow(elem.imaginary, 2));
    }
    return sqrtf(res);
}