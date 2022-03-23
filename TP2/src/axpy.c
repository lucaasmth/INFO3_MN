#include "mnblas.h"
#include <stdio.h>
#include "../include/mnblas.h"
#include "../include/complexe2.h"


void mncblas_saxpy(const int N, const float a, const float *X, const int incX, 
                 float *Y, const int incY)
{
  register unsigned int i = 0 ;
  register unsigned int j = 0 ;
  
  for (i = 0 ; i < N ; i += incX)
    {
      Y[j] += a * X [i] * Y [j] ;
      j+=incY ;
    }

}

void mncblas_daxpy(const int N, const double a, const double *X, const int incX, 
                  double *Y, const int incY)
{
    register unsigned int i = 0 ;
    register unsigned int j = 0 ;
  
    for (i = 0 ; i < N ; i += incX){
        Y[j] += a * X [i] * Y [j] ;
        j+=incY ;
    }
}

void   mncblas_caxpy(const int N, const void *a, const void *X, const int incX,
                        void *Y, const int incY, void *dotu)
{
    register unsigned int i = 0 ;
    register unsigned int j = 0 ;
  
    for (i = 0 ; i < N ; i += incX){
        complexe_float_t multiplication = mult_complexe_float(*(complexe_float_t *)a, ((complexe_float_t *)X)[i]);
        ((complexe_float_t *)Y)[j] = add_complexe_float(multiplication, ((complexe_float_t *)Y)[j]);
        j += incY;
    }
}

  
void   mncblas_zaxpy(const int N, const void *a, const void *X, const int incX,
                        void *Y, const int incY, void *dotc)
{
    register unsigned int i = 0 ;
    register unsigned int j = 0 ;
  
    for (i = 0 ; i < N ; i += incX){
        complexe_double_t multiplication = mult_complexe_double(*(complexe_double_t *)a, ((complexe_double_t *)X)[i]);
        ((complexe_double_t *)Y)[j] = add_complexe_double(multiplication, ((complexe_double_t *)Y)[j]);
        j += incY;
    }
}
