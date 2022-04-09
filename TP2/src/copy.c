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

#pragma omp parallel for
    for (i=0; i < N ; i++)
    {
        *((complexe_float_t*)(Y+i)) = *((complexe_float_t*)(X+i)) ;
    }
    return ;

}

void mncblas_zcopy(const int N, const void *X, const int incX, 
		                    void *Y, const int incY)
{
    register unsigned int i = 0 ;

#pragma omp parallel for
    for (i=0; i < N ; i++)
    {
        *((complexe_double_t*)(Y+i)) = *((complexe_double_t*)(X+i)) ;
    }
    return ;
}

