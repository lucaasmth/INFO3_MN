#include "mnblas.h"
#include "../include/complexe2.h"

void mncblas_scopy(const int N, const float *X, const int incX, 
                 float *Y, const int incY)
{
  register unsigned int i = 0 ;

#pragma omp parallel for
  for (i=0; i < N ; i++)
    {
      Y[i] = X[i];
    }

  return ;
}

void mncblas_dcopy(const int N, const double *X, const int incX, 
                 double *Y, const int incY)
{
    register unsigned int i = 0 ;

#pragma omp parallel for
    for (i=0; i < N ; i++)
    {
        Y[i] = X[i];
    }

    return ;
}

void mncblas_ccopy(const int N, const void *X, const int incX, 
		                    void *Y, const int incY)
{

    register unsigned int i = 0 ;

    complexe_float_t* X_copy = (complexe_float_t *) X;
    complexe_float_t* Y_copy = (complexe_float_t *) Y;

#pragma omp parallel for
    for (i=0; i < N ; i++)
    {
        Y_copy[i] = X_copy[i];
    }

    X = (void *) X_copy;
    Y = (void *) Y_copy;

    return ;

}

void mncblas_zcopy(const int N, const void *X, const int incX, 
		                    void *Y, const int incY)
{
    register unsigned int i = 0 ;

    complexe_double_t* X_copy = (complexe_double_t *) X;
    complexe_double_t* Y_copy = (complexe_double_t *) Y;

#pragma omp parallel for
    for (i=0; i < N ; i++)
    {
        Y_copy[i] = X_copy[i];
    }

    X = (void *) X_copy;
    Y = (void *) Y_copy;

    return ;
}

