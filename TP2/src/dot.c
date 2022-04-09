#include "../include/mnblas.h"
#include <stdio.h>
#include "../include/complexe.h"

float mncblas_sdot(const int N, const float *X, const int incX,
                   const float *Y, const int incY) {
    register unsigned int i = 0;
    register unsigned int j = 0;
    float dot = 0.0;

#pragma omp parallel for
    for (i = 0; i < N; i += incX) {
        dot += X[i] * Y[j];
        j += incY;
    }

    return dot;
}

double mncblas_ddot(const int N, const double *X, const int incX,
                    const double *Y, const int incY) {
    register unsigned int i = 0;
    register unsigned int j = 0;
    double dot = 0.0;

#pragma omp parallel for
    for (i = 0; i < N; i += incX) {
        dot += X[i] * Y[j];
        j += incY;
    }

    return dot;
}

void mncblas_cdotu_sub(const int N, const void *X, const int incX,
                       const void *Y, const int incY, void *dotu) {
    register unsigned int i = 0;
    register unsigned int j = 0;
    register complexe_float_t dot;
    register complexe_float_t multiplication;
#pragma omp parallel for
    for (i = 0; i < N; i += incX) {
        multiplication = mult_complexe_float(((complexe_float_t *) X)[i], ((complexe_float_t *) Y)[j]);
        dot = add_complexe_float(dot, multiplication);
        j += incY;
    }
    *(complexe_float_t *) dotu = dot;
}

void mncblas_cdotc_sub(const int N, const void *X, const int incX,
                       const void *Y, const int incY, void *dotc) {
    register unsigned int i = 0;
    register unsigned int j = 0;
    register complexe_float_t dot;
    register complexe_float_t multiplication;
    register complexe_float_t conjugue_X;
#pragma omp parallel for
    for (i = 0; i < N; i += incX) {
        conjugue_X.real = ((complexe_float_t *) X)[i].real;
        conjugue_X.imaginary = -((complexe_float_t *) X)[i].imaginary;
        multiplication = mult_complexe_float(conjugue_X, ((complexe_float_t *) Y)[j]);
        dot = add_complexe_float(dot, multiplication);
        j += incY;
    }
    *(complexe_float_t *) dotc = dot;
}

void mncblas_zdotu_sub(const int N, const void *X, const int incX,
                       const void *Y, const int incY, void *dotu) {
    register unsigned int i = 0;
    register unsigned int j = 0;
    register complexe_double_t dot;
    register complexe_double_t multiplication;
#pragma omp parallel for
    for (i = 0; i < N; i += incX) {
        multiplication = mult_complexe_double(((complexe_double_t *) X)[i], ((complexe_double_t *) Y)[j]);
        dot = add_complexe_double(dot, multiplication);
        j += incY;
    }
    *(complexe_double_t *) dotu = dot;
}

void mncblas_zdotc_sub(const int N, const void *X, const int incX,
                       const void *Y, const int incY, void *dotc) {
    register unsigned int i = 0;
    register unsigned int j = 0;
    register complexe_double_t dot;
    register complexe_double_t multiplication;
    register complexe_double_t conjugue_X;
#pragma omp parallel for
    for (i = 0; i < N; i += incX) {
        conjugue_X.real = ((complexe_double_t *) X)[i].real;
        conjugue_X.imaginary = -((complexe_double_t *) X)[i].imaginary;
        multiplication = mult_complexe_double(conjugue_X, ((complexe_double_t *) Y)[j]);
        dot = add_complexe_double(dot, multiplication);
        j += incY;
    }
    *(complexe_double_t *) dotc = dot;
}




