#include "../include/mnblas.h"
#include "../include/complexe2.h"
#include <stdlib.h>
#include <stdio.h>

/*
Layout : indique si trié par colonne (MNCblasColMajor) ou par ligne (MNCblasRowMajor)

TransA :
    indique l'opération:
    if trans = CblasNoTrans, then y:= alpha*A*x+ beta*y;
    if trans = CblasTrans, then y:= alpha*A'*x+ beta*y ;
    if trans = CblasConjTrans , then y:= alpha*conjg(A')*x + beta*y.


M : nombre de lignes
N : nombre de colonnes

alpha : scalaire
beta : scalaire

A :  matrice de taille

lda :  taille du tableau en une dimension
    si CblasColMajor : lda vaut au entre 1 et m
    si CblasRowMajor : lda vaut au entre 1 et n
    
x : vecteur
incX :incrément de X

y : vecteur
incY : incrément de Y
*/

void mncblas_sgemv(const MNCBLAS_LAYOUT layout, const MNCBLAS_TRANSPOSE TransA, const int M, const int N, const float alpha,
                   const float *A, const int lda, const float *X, const int incX, const float beta, float *Y, const int incY)
{

    if(layout == MNCblasRowMajor) {
        if(TransA == MNCblasNoTrans) {
            register float sommeligne ;
            for(register unsigned int i = 0, j; i < M ; i+= incY) {
                sommeligne = 0;
                for(j = 0; j < N ; j+= incX) {
                    sommeligne += *(A+i*N+j)*X[j];
                }
                Y[i] = alpha*sommeligne + beta*Y[i];
            }
        } else if (TransA == MNCblasTrans) {
            register float sommeligne ;
            for(register unsigned int i = 0; i < N ; i+= incY) {
                sommeligne = 0;
                for(register int j = 0; j < M ; j+= incX) {
                    sommeligne += *(A+j*N+i)*X[j];
                }
                Y[i] = alpha*sommeligne + beta*Y[i];
            }
        }
    } else if (layout == MNCblasColMajor){
        if(TransA == MNCblasNoTrans) {
            register float sommecolonne ;
            for(register unsigned int i = 0, j; i < M ; i+= incY) {
                sommecolonne = 0;
                for(j = 0; j < N ; j+= incX) {
                    sommecolonne += *(A+j*M+i)*X[j];
                }
                Y[i] = alpha*sommecolonne + beta*Y[i];
            }
        } else if (TransA == MNCblasTrans) {
            register float sommecolonne ;
            for(register unsigned int i = 0, j; i < N ; i+= incY) {
                sommecolonne = 0;
                for(j = 0; j < M ; j+= incX) {
                    sommecolonne += *(A+i*M+j)*X[j];
                }
                Y[i] = alpha*sommecolonne + beta*Y[i];
            }
        }
    }
}

void mncblas_dgemv(MNCBLAS_LAYOUT layout, MNCBLAS_TRANSPOSE TransA, const int M, const int N, const double alpha, const double *A, const int lda,
                   const double *X, const int incX, const double beta, double *Y, const int incY)
{
 if(layout == MNCblasRowMajor) {
        if(TransA == MNCblasNoTrans) {
            register double sommeligne ;
            for(register unsigned int i = 0, j; i < M ; i+= incY) {
                sommeligne = 0;
                for(j = 0; j < N ; j+= incX) {
                    sommeligne += *(A+i*N+j)*X[j];
                }
                Y[i] = alpha*sommeligne + beta*Y[i];
            }
        } else if (TransA == MNCblasTrans) {
            register double sommeligne ;
            for(register unsigned int i = 0, j; i < N ; i+= incY) {
                sommeligne = 0;
                for(j = 0; j < M ; j+= incX) {
                    sommeligne += *(A+j*N+i)*X[j];
                }
                Y[i] = alpha*sommeligne + beta*Y[i];
            }
        }
    } else if (layout == MNCblasColMajor){
        if(TransA == MNCblasNoTrans) {
            register double sommecolonne ;
            for(register unsigned int i = 0, j; i < M ; i+= incY) {
                sommecolonne = 0;
                for(j = 0; j < N ; j+= incX) {
                    sommecolonne += *(A+j*M+i)*X[j];
                }
                Y[i] = alpha*sommecolonne + beta*Y[i];
            }
        } else if (TransA == MNCblasTrans) {
            register double sommecolonne ;
            for(register unsigned int i = 0, j; i < N ; i+= incY) {
                sommecolonne = 0;
                for(j = 0; j < M ; j+= incX) {
                    sommecolonne += *(A+i*M+j)*X[j];
                }
                Y[i] = alpha*sommecolonne + beta*Y[i];
            }
        }
    }
}

void mncblas_cgemv(MNCBLAS_LAYOUT layout, MNCBLAS_TRANSPOSE TransA,
                     const int M, const int N,
                     const void *alpha, const void *A, const int lda, const void *X, const int incX, const void *beta, void *Y, const int incY) { // NB OPE FLOTANTE = 8*M*N + 14*M
    if(layout == MNCblasRowMajor) {
        if(TransA == MNCblasNoTrans) {
            register complexe_float_t sommeligne ;
            for(register unsigned int i = 0, j; i < M ; i+= incY) {
                sommeligne.real = 0; sommeligne.imaginary = 0;
                for(j = 0; j < N ; j+= incX) {
                    sommeligne = add_complexe_float(sommeligne,mult_complexe_float(*(((complexe_float_t*)A)+i*N+j),*(((complexe_float_t*)X)+j)));
                }
                *(((complexe_float_t*)Y)+i) = add_complexe_float(mult_complexe_float(*((complexe_float_t*)alpha),sommeligne),mult_complexe_float(*((complexe_float_t*)beta),*(((complexe_float_t*)Y)+i)));
            }
        } else if (TransA == MNCblasTrans) {
            register complexe_float_t sommeligne ;
            for(register unsigned int i = 0, j; i < N ; i+= incY) {
                sommeligne.real = 0; sommeligne.imaginary = 0;
                for(j = 0; j < M ; j+= incX) {
                    sommeligne = add_complexe_float(sommeligne,mult_complexe_float(*(((complexe_float_t*)A)+j*N+i),*(((complexe_float_t*)X)+j)));
                }
                *(((complexe_float_t*)Y)+i) = add_complexe_float(mult_complexe_float(*((complexe_float_t*)alpha),sommeligne),mult_complexe_float(*((complexe_float_t*)beta),*(((complexe_float_t*)Y)+i)));
            }
        } else if (TransA == MNCblasConjTrans) {
            register complexe_float_t sommeligne ;
            for(register unsigned int i = 0, j; i < N ; i+= incY) {
                sommeligne.real = 0; sommeligne.imaginary = 0;
                for(j = 0; j < M ; j+= incX) {
                    (((complexe_float_t*)A)+j*N+i)->imaginary *= -1;
                    sommeligne = add_complexe_float(sommeligne,mult_complexe_float(*(((complexe_float_t*)A)+j*N+i),*(((complexe_float_t*)X)+j)));
                }
                *(((complexe_float_t*)Y)+i) = add_complexe_float(mult_complexe_float(*((complexe_float_t*)alpha),sommeligne),mult_complexe_float(*((complexe_float_t*)beta),*(((complexe_float_t*)Y)+i)));
            }
        }
    } else if (layout == MNCblasColMajor){
        if(TransA == MNCblasNoTrans) {
            register complexe_float_t sommecolonne ;
            for(register unsigned int i = 0, j; i < M ; i+= incY) {
                sommecolonne.real = 0; sommecolonne.imaginary = 0;
                for(j = 0; j < N ; j+= incX) {
                    sommecolonne = add_complexe_float(sommecolonne,mult_complexe_float(*(((complexe_float_t*)A)+j*M+i),*(((complexe_float_t*)X)+j)));
                }
                *(((complexe_float_t*)Y)+i) = add_complexe_float(mult_complexe_float(*((complexe_float_t*)alpha),sommecolonne),mult_complexe_float(*((complexe_float_t*)beta),*(((complexe_float_t*)Y)+i)));
            }
        } else if (TransA == MNCblasTrans) {
            register complexe_float_t sommecolonne ;
            for(register unsigned int i = 0, j; i < N ; i+= incY) {
                sommecolonne.real = 0; sommecolonne.imaginary = 0;
                for(j = 0; j < M ; j+= incX) {
                    sommecolonne = add_complexe_float(sommecolonne,mult_complexe_float(*(((complexe_float_t*)A)+i*M+j),*(((complexe_float_t*)X)+j)));
                }
                *(((complexe_float_t*)Y)+i) = add_complexe_float(mult_complexe_float(*((complexe_float_t*)alpha),sommecolonne),mult_complexe_float(*((complexe_float_t*)beta),*(((complexe_float_t*)Y)+i)));
            }
        }
        else if (TransA == MNCblasConjTrans) {
            register complexe_float_t sommecolonne ;
            for(register unsigned int i = 0, j; i < N ; i+= incY) {
                sommecolonne.real = 0; sommecolonne.imaginary = 0;
                for(j = 0; j < M ; j+= incX) {
                    (((complexe_float_t*)A)+i*M+j)->imaginary *= -1;
                    sommecolonne = add_complexe_float(sommecolonne,mult_complexe_float(*(((complexe_float_t*)A)+i*M+j),*(((complexe_float_t*)X)+j)));
                }
                *(((complexe_float_t*)Y)+i) = add_complexe_float(mult_complexe_float(*((complexe_float_t*)alpha),sommecolonne),mult_complexe_float(*((complexe_float_t*)beta),*(((complexe_float_t*)Y)+i)));
            }
        }
    }
}

void mncblas_zgemv(MNCBLAS_LAYOUT layout, MNCBLAS_TRANSPOSE TransA,
                     const int M, const int N,
                     const void *alpha, const void *A, const int lda, const void *X, const int incX, const void *beta, void *Y, const int incY){ 
    if(layout == MNCblasRowMajor) {
        if(TransA == MNCblasNoTrans) {
            register complexe_double_t sommeligne ;
            for(register unsigned int i = 0, j; i < M ; i+= incY) {
                sommeligne.real = 0; sommeligne.imaginary = 0;
                for(j = 0; j < N ; j+= incX) {
                    sommeligne = add_complexe_double(sommeligne,mult_complexe_double(*(((complexe_double_t*)A)+i*N+j),*(((complexe_double_t*)X)+j)));
                }
                *(((complexe_double_t*)Y)+i) = add_complexe_double(mult_complexe_double(*((complexe_double_t*)alpha),sommeligne),mult_complexe_double(*((complexe_double_t*)beta),*(((complexe_double_t*)Y)+i)));
            }
        } else if (TransA == MNCblasTrans) {
            register complexe_double_t sommeligne ;
            for(register unsigned int i = 0, j; i < N ; i+= incY) {
                sommeligne.real = 0; sommeligne.imaginary = 0;
                for(j = 0; j < M ; j+= incX) {
                    sommeligne = add_complexe_double(sommeligne,mult_complexe_double(*(((complexe_double_t*)A)+j*N+i),*(((complexe_double_t*)X)+j)));
                }
                *(((complexe_double_t*)Y)+i) = add_complexe_double(mult_complexe_double(*((complexe_double_t*)alpha),sommeligne),mult_complexe_double(*((complexe_double_t*)beta),*(((complexe_double_t*)Y)+i)));
            }
        } else if (TransA == MNCblasConjTrans) {
            register complexe_double_t sommeligne ;
            for(register unsigned int i = 0, j; i < N ; i+= incY) {
                sommeligne.real = 0; sommeligne.imaginary = 0;
                for(j = 0; j < M ; j+= incX) {
                    (((complexe_double_t*)A)+j*N+i)->imaginary *= -1;
                    sommeligne = add_complexe_double(sommeligne,mult_complexe_double(*(((complexe_double_t*)A)+j*N+i),*(((complexe_double_t*)X)+j)));
                }
                *(((complexe_double_t*)Y)+i) = add_complexe_double(mult_complexe_double(*((complexe_double_t*)alpha),sommeligne),mult_complexe_double(*((complexe_double_t*)beta),*(((complexe_double_t*)Y)+i)));
            }
        }
    } else if (layout == MNCblasColMajor){
        if(TransA == MNCblasNoTrans) {
            register complexe_double_t sommecolonne ;
            for(register unsigned int i = 0, j; i < M ; i+= incY) {
                sommecolonne.real = 0; sommecolonne.imaginary = 0;
                for(j = 0; j < N ; j+= incX) {
                    sommecolonne = add_complexe_double(sommecolonne,mult_complexe_double(*(((complexe_double_t*)A)+j*M+i),*(((complexe_double_t*)X)+j)));
                }
                *(((complexe_double_t*)Y)+i) = add_complexe_double(mult_complexe_double(*((complexe_double_t*)alpha),sommecolonne),mult_complexe_double(*((complexe_double_t*)beta),*(((complexe_double_t*)Y)+i)));
            }
        } else if (TransA == MNCblasTrans) {
            register complexe_double_t sommecolonne ;
            for(register unsigned int i = 0, j; i < N ; i+= incY) {
                sommecolonne.real = 0; sommecolonne.imaginary = 0;
                for(j = 0; j < M ; j+= incX) {
                    sommecolonne = add_complexe_double(sommecolonne,mult_complexe_double(*(((complexe_double_t*)A)+i*M+j),*(((complexe_double_t*)X)+j)));
                }
                *(((complexe_double_t*)Y)+i) = add_complexe_double(mult_complexe_double(*((complexe_double_t*)alpha),sommecolonne),mult_complexe_double(*((complexe_double_t*)beta),*(((complexe_double_t*)Y)+i)));
            }
        } else if (TransA == MNCblasConjTrans) {
            register complexe_double_t sommecolonne ;
            for(register unsigned int i = 0, j; i < N ; i+= incY) {
                sommecolonne.real = 0; sommecolonne.imaginary = 0;
                for(j = 0; j < M ; j+= incX) {
                    (((complexe_double_t*)A)+i*M+j)->imaginary *= -1;
                    sommecolonne = add_complexe_double(sommecolonne,mult_complexe_double(*(((complexe_double_t*)A)+i*M+j),*(((complexe_double_t*)X)+j)));
                }
                *(((complexe_double_t*)Y)+i) = add_complexe_double(mult_complexe_double(*((complexe_double_t*)alpha),sommecolonne),mult_complexe_double(*((complexe_double_t*)beta),*(((complexe_double_t*)Y)+i)));
            }
        }
    }
}
