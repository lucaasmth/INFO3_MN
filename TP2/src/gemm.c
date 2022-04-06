#include "../include/mnblas.h"
#include "../include/complexe2.h"

void mncblas_sgemm(MNCBLAS_LAYOUT layout, MNCBLAS_TRANSPOSE TransA, MNCBLAS_TRANSPOSE TransB, const int M, const int N, const int K, const float alpha, const float *A, const int lda, const float *B, const int ldb, const float beta, float *C, const int ldc) {
    if(layout == MNCblasRowMajor) {
        if(TransA == MNCblasNoTrans) {
            if(TransB == MNCblasNoTrans) {
                float somme;
                for(int i = 0; i < M; i++) {
                    for(int j = 0; j < N; j++) {
                        somme = 0;
                        for(int k = 0; k < K; k++) {
                            somme += A[i * lda + k] * B[k * ldb + j];
                        }
                        C[i * ldc + j] = alpha * somme + beta * C[i * ldc + j];
                    }
                }
            } else {
                float somme;
                for(int i = 0; i < M; i++) {
                    for(int j = 0; j < N; j++) {
                        somme = 0;
                        for(int k = 0; k < K; k++) {
                            somme += A[i * lda + k] * B[j * ldb + k];
                        }
                        C[i * ldc + j] = alpha * somme + beta * C[i * ldc + j];
                    }
                }
            }
        } else {
            if(TransB == MNCblasNoTrans) {
                float somme;
                for(int i = 0; i < M; i++) {
                    for(int j = 0; j < N; j++) {
                        somme = 0;
                        for(int k = 0; k < K; k++) {
                            somme += A[k * lda + i] * B[k * ldb + j];
                        }
                        C[i * ldc + j] = alpha * somme + beta * C[i * ldc + j];
                    }
                }
            } else {
                float somme;
                for(int i = 0; i < M; i++) {
                    for(int j = 0; j < N; j++) {
                        somme = 0;
                        for(int k = 0; k < K; k++) {
                            somme += A[k * lda + i] * B[j * ldb + k];
                        }
                        C[i * ldc + j] = alpha * somme + beta * C[i * ldc + j];
                    }
                }
            }
        }
    } else {
        if(TransA == MNCblasNoTrans) {
            if(TransB == MNCblasNoTrans) {
                mncblas_sgemm(layout, MNCblasTrans, MNCblasTrans, N, M, K, alpha, B, ldb, A, lda, beta, C, ldc);
            } else {
                mncblas_sgemm(layout, MNCblasTrans, MNCblasNoTrans, N, M, K, alpha, B, ldb, A, lda, beta, C, ldc);
            }
        } else {
            if(TransB == MNCblasNoTrans) {
                mncblas_sgemm(layout, MNCblasNoTrans, MNCblasTrans, M, N, K, alpha, A, lda, B, ldb, beta, C, ldc);
            } else {
                mncblas_sgemm(layout, MNCblasNoTrans, MNCblasNoTrans, M, N, K, alpha, A, lda, B, ldb, beta, C, ldc);
            }
        }
    }
}

void mncblas_dgemm(MNCBLAS_LAYOUT layout, MNCBLAS_TRANSPOSE TransA, MNCBLAS_TRANSPOSE TransB, const int M, const int N, const int K, const double alpha, const double *A, const int lda, const double *B, const int ldb, const double beta, double *C, const int ldc) {
    if(layout == MNCblasRowMajor) {
        if(TransA == MNCblasNoTrans) {
            if(TransB == MNCblasNoTrans) {
                double somme;
                for(int i = 0; i < M; i++) {
                    for(int j = 0; j < N; j++) {
                        somme = 0;
                        for(int k = 0; k < K; k++) {
                            somme += A[i * lda + k] * B[k * ldb + j];
                        }
                        C[i * ldc + j] = alpha * somme + beta * C[i * ldc + j];
                    }
                }
            } else {
                double somme;
                for(int i = 0; i < M; i++) {
                    for(int j = 0; j < N; j++) {
                        somme = 0;
                        for(int k = 0; k < K; k++) {
                            somme += A[i * lda + k] * B[j * ldb + k];
                        }
                        C[i * ldc + j] = alpha * somme + beta * C[i * ldc + j];
                    }
                }
            }
        } else {
            if(TransB == MNCblasNoTrans) {
                double somme;
                for(int i = 0; i < M; i++) {
                    for(int j = 0; j < N; j++) {
                        somme = 0;
                        for(int k = 0; k < K; k++) {
                            somme += A[k * lda + i] * B[k * ldb + j];
                        }
                        C[i * ldc + j] = alpha * somme + beta * C[i * ldc + j];
                    }
                }
            } else {
                double somme;
                for(int i = 0; i < M; i++) {
                    for(int j = 0; j < N; j++) {
                        somme = 0;
                        for(int k = 0; k < K; k++) {
                            somme += A[k * lda + i] * B[j * ldb + k];
                        }
                        C[i * ldc + j] = alpha * somme + beta * C[i * ldc + j];
                    }
                }
            }
        }
    } else {
        if(TransA == MNCblasNoTrans) {
            if(TransB == MNCblasNoTrans) {
                mncblas_dgemm(layout, MNCblasTrans, MNCblasTrans, N, M, K, alpha, B, ldb, A, lda, beta, C, ldc);
            } else {
                mncblas_dgemm(layout, MNCblasTrans, MNCblasNoTrans, N, M, K, alpha, B, ldb, A, lda, beta, C, ldc);
            }
        } else {
            if(TransB == MNCblasNoTrans) {
                mncblas_dgemm(layout, MNCblasNoTrans, MNCblasTrans, M, N, K, alpha, A, lda, B, ldb, beta, C, ldc);
            } else {
                mncblas_dgemm(layout, MNCblasNoTrans, MNCblasNoTrans, M, N, K, alpha, A, lda, B, ldb, beta, C, ldc);
            }
        }
    }
}

void mncblas_cgemm(MNCBLAS_LAYOUT layout, MNCBLAS_TRANSPOSE TransA, MNCBLAS_TRANSPOSE TransB, const int M, const int N, const int K, const void *alpha, const void *A, const int lda, const void *B, const int ldb, const void *beta, void *C, const int ldc) {
    if(layout == MNCblasRowMajor) {
        if(TransA == MNCblasNoTrans) {
            if(TransB == MNCblasNoTrans) {
                complexe_double_t somme;
                for(int i = 0; i < M; i++) {
                    for(int j = 0; j < N; j++) {
                        somme.real = 0;
                        somme.imaginary = 0;
                        for(int k = 0; k < K; k++) {
                            somme = add_complexe_double(somme, mult_complexe_double(((complexe_double_t *)A)[i * lda + k], ((complexe_double_t *)B)[k * ldb + j]));
                        }
                        ((complexe_double_t *)C)[i * ldc + j].real = ((complexe_float_t*)alpha)->real * somme.real + ((complexe_float_t*)alpha)->imaginary * somme.imaginary + ((complexe_float_t*)beta)->real * ((complexe_double_t *)C)[i * ldc + j].real - ((complexe_float_t*)beta)->imaginary * ((complexe_double_t *)C)[i * ldc + j].imaginary;
                        ((complexe_double_t *)C)[i * ldc + j].imaginary = ((complexe_float_t*)alpha)->real * somme.imaginary + ((complexe_float_t*)alpha)->imaginary * somme.real + ((complexe_float_t*)beta)->real * ((complexe_double_t *)C)[i * ldc + j].imaginary + ((complexe_float_t*)beta)->imaginary * ((complexe_double_t *)C)[i * ldc + j].real;
                    }
                }
            } else if(TransB == MNCblasTrans) {
                complexe_double_t somme;
                for(int i = 0; i < M; i++) {
                    for(int j = 0; j < N; j++) {
                        somme.real = 0;
                        somme.imaginary = 0;
                        for(int k = 0; k < K; k++) {
                            somme = add_complexe_double(somme, mult_complexe_double(((complexe_double_t *)A)[i * lda + k], ((complexe_double_t *)B)[j * ldb + k]));
                        }
                        ((complexe_double_t *)C)[i * ldc + j].real = ((complexe_float_t*)alpha)->real * somme.real + ((complexe_float_t*)alpha)->imaginary * somme.imaginary + ((complexe_float_t*)beta)->real * ((complexe_double_t *)C)[i * ldc + j].real - ((complexe_float_t*)beta)->imaginary * ((complexe_double_t *)C)[i * ldc + j].imaginary;
                        ((complexe_double_t *)C)[i * ldc + j].imaginary = ((complexe_float_t*)alpha)->real * somme.imaginary + ((complexe_float_t*)alpha)->imaginary * somme.real + ((complexe_float_t*)beta)->real * ((complexe_double_t *)C)[i * ldc + j].imaginary + ((complexe_float_t*)beta)->imaginary * ((complexe_double_t *)C)[i * ldc + j].real;
                    }
                }
            } else {
                complexe_double_t somme;
                for(int i = 0; i < M; i++) {
                    for(int j = 0; j < N; j++) {
                        somme.real = 0;
                        somme.imaginary = 0;
                        for(int k = 0; k < K; k++) {
                            somme = add_complexe_double(somme, mult_complexe_double(((complexe_double_t *)A)[i * lda + k], ((complexe_double_t *)B)[j * ldb + k]));
                        }
                        ((complexe_double_t *)C)[i * ldc + j].real = ((complexe_float_t*)alpha)->real * somme.real - ((complexe_float_t*)alpha)->imaginary * somme.imaginary + ((complexe_float_t*)beta)->real * ((complexe_double_t *)C)[i * ldc + j].real - ((complexe_float_t*)beta)->imaginary * ((complexe_double_t *)C)[i * ldc + j].imaginary;
                        ((complexe_double_t *)C)[i * ldc + j].imaginary = ((complexe_float_t*)alpha)->real * somme.imaginary + ((complexe_float_t*)alpha)->imaginary * somme.real + ((complexe_float_t*)beta)->real * ((complexe_double_t *)C)[i * ldc + j].imaginary + ((complexe_float_t*)beta)->imaginary * ((complexe_double_t *)C)[i * ldc + j].real;
                    }
                }
            }
        } else if (TransA == MNCblasTrans) {
            if(TransB == MNCblasNoTrans) {
                complexe_double_t somme;
                for(int i = 0; i < M; i++) {
                    for(int j = 0; j < N; j++) {
                        somme.real = 0;
                        somme.imaginary = 0;
                        for(int k = 0; k < K; k++) {
                            somme = add_complexe_double(somme, mult_complexe_double(((complexe_double_t *)A)[k * lda + i], ((complexe_double_t *)B)[k * ldb + j]));
                        }
                        ((complexe_double_t *)C)[i * ldc + j].real = ((complexe_float_t*)alpha)->real * somme.real + ((complexe_float_t*)alpha)->imaginary * somme.imaginary + ((complexe_float_t*)beta)->real * ((complexe_double_t *)C)[i * ldc + j].real - ((complexe_float_t*)beta)->imaginary * ((complexe_double_t *)C)[i * ldc + j].imaginary;
                        ((complexe_double_t *)C)[i * ldc + j].imaginary = ((complexe_float_t*)alpha)->real * somme.imaginary - ((complexe_float_t*)alpha)->imaginary * somme.real + ((complexe_float_t*)beta)->real * ((complexe_double_t *)C)[i * ldc + j].imaginary + ((complexe_float_t*)beta)->imaginary * ((complexe_double_t *)C)[i * ldc + j].real;
                    }
                }
            } else if(TransB == MNCblasTrans) {
                complexe_double_t somme;
                for(int i = 0; i < M; i++) {
                    for(int j = 0; j < N; j++) {
                        somme.real = 0;
                        somme.imaginary = 0;
                        for(int k = 0; k < K; k++) {
                            somme = add_complexe_double(somme, mult_complexe_double(((complexe_double_t *)A)[k * lda + i], ((complexe_double_t *)B)[j * ldb + k]));
                        }
                        ((complexe_double_t *)C)[i * ldc + j].real = ((complexe_float_t*)alpha)->real * somme.real + ((complexe_float_t*)alpha)->imaginary * somme.imaginary + ((complexe_float_t*)beta)->real * ((complexe_double_t *)C)[i * ldc + j].real - ((complexe_float_t*)beta)->imaginary * ((complexe_double_t *)C)[i * ldc + j].imaginary;
                        ((complexe_double_t *)C)[i * ldc + j].imaginary = ((complexe_float_t*)alpha)->real * somme.imaginary + ((complexe_float_t*)alpha)->imaginary * somme.real + ((complexe_float_t*)beta)->real * ((complexe_double_t *)C)[i * ldc + j].imaginary + ((complexe_float_t*)beta)->imaginary * ((complexe_double_t *)C)[i * ldc + j].real;
                    }
                }
            } else {
                complexe_double_t somme;
                for(int i = 0; i < M; i++) {
                    for(int j = 0; j < N; j++) {
                        somme.real = 0;
                        somme.imaginary = 0;
                        for(int k = 0; k < K; k++) {
                            somme = add_complexe_double(somme, mult_complexe_double(((complexe_double_t *)A)[k * lda + i], ((complexe_double_t *)B)[j * ldb + k]));
                        }
                        ((complexe_double_t *)C)[i * ldc + j].real = ((complexe_float_t*)alpha)->real * somme.real - ((complexe_float_t*)alpha)->imaginary * somme.imaginary + ((complexe_float_t*)beta)->real * ((complexe_double_t *)C)[i * ldc + j].real - ((complexe_float_t*)beta)->imaginary * ((complexe_double_t *)C)[i * ldc + j].imaginary;
                        ((complexe_double_t *)C)[i * ldc + j].imaginary = ((complexe_float_t*)alpha)->real * somme.imaginary + ((complexe_float_t*)alpha)->imaginary * somme.real + ((complexe_float_t*)beta)->real * ((complexe_double_t *)C)[i * ldc + j].imaginary + ((complexe_float_t*)beta)->imaginary * ((complexe_double_t *)C)[i * ldc + j].real;
                    }
                }
            }
        }
    } else {
        if(TransA == MNCblasNoTrans) {
            if(TransB == MNCblasNoTrans) {
                mncblas_cgemm(layout, MNCblasTrans, MNCblasTrans, N, M, K, alpha, B, ldb, A, lda, beta, C, ldc);
            } else {
                mncblas_cgemm(layout, MNCblasTrans, MNCblasNoTrans, N, M, K, alpha, B, ldb, A, lda, beta, C, ldc);
            }
        } else {
            if(TransB == MNCblasNoTrans) {
                mncblas_cgemm(layout, MNCblasNoTrans, MNCblasTrans, M, N, K, alpha, A, lda, B, ldb, beta, C, ldc);
            } else {
                mncblas_cgemm(layout, MNCblasNoTrans, MNCblasNoTrans, M, N, K, alpha, A, lda, B, ldb, beta, C, ldc);
            }
        }
    }
}

void mncblas_zgemm(MNCBLAS_LAYOUT layout, MNCBLAS_TRANSPOSE TransA, MNCBLAS_TRANSPOSE TransB, const int M, const int N, const int K, const void *alpha, const void *A, const int lda, const void *B, const int ldb, const void *beta, void *C, const int ldc) {
    if(layout == MNCblasRowMajor) {
        if(TransA == MNCblasNoTrans) {
            if(TransB == MNCblasNoTrans) {
                complexe_float_t somme;
                for(int i = 0; i < M; i++) {
                    for(int j = 0; j < N; j++) {
                        somme.real = 0;
                        somme.imaginary = 0;
                        for(int k = 0; k < K; k++) {
                            somme = add_complexe_float(somme, mult_complexe_float(((complexe_float_t *)A)[i * lda + k], ((complexe_float_t *)B)[k * ldb + j]));
                        }
                        ((complexe_float_t *)C)[i * ldc + j].real = ((complexe_float_t*)alpha)->real * somme.real + ((complexe_float_t*)alpha)->imaginary * somme.imaginary + ((complexe_float_t*)beta)->real * ((complexe_float_t *)C)[i * ldc + j].real - ((complexe_float_t*)beta)->imaginary * ((complexe_float_t *)C)[i * ldc + j].imaginary;
                        ((complexe_float_t *)C)[i * ldc + j].imaginary = ((complexe_float_t*)alpha)->real * somme.imaginary + ((complexe_float_t*)alpha)->imaginary * somme.real + ((complexe_float_t*)beta)->real * ((complexe_float_t *)C)[i * ldc + j].imaginary + ((complexe_float_t*)beta)->imaginary * ((complexe_float_t *)C)[i * ldc + j].real;
                    }
                }
            } else if(TransB == MNCblasTrans) {
                complexe_float_t somme;
                for(int i = 0; i < M; i++) {
                    for(int j = 0; j < N; j++) {
                        somme.real = 0;
                        somme.imaginary = 0;
                        for(int k = 0; k < K; k++) {
                            somme = add_complexe_float(somme, mult_complexe_float(((complexe_float_t *)A)[i * lda + k], ((complexe_float_t *)B)[j * ldb + k]));
                        }
                        ((complexe_float_t *)C)[i * ldc + j].real = ((complexe_float_t*)alpha)->real * somme.real + ((complexe_float_t*)alpha)->imaginary * somme.imaginary + ((complexe_float_t*)beta)->real * ((complexe_float_t *)C)[i * ldc + j].real - ((complexe_float_t*)beta)->imaginary * ((complexe_float_t *)C)[i * ldc + j].imaginary;
                        ((complexe_float_t *)C)[i * ldc + j].imaginary = ((complexe_float_t*)alpha)->real * somme.imaginary + ((complexe_float_t*)alpha)->imaginary * somme.real + ((complexe_float_t*)beta)->real * ((complexe_float_t *)C)[i * ldc + j].imaginary + ((complexe_float_t*)beta)->imaginary * ((complexe_float_t *)C)[i * ldc + j].real;
                    }
                }
            } else {
                complexe_float_t somme;
                for(int i = 0; i < M; i++) {
                    for(int j = 0; j < N; j++) {
                        somme.real = 0;
                        somme.imaginary = 0;
                        for(int k = 0; k < K; k++) {
                            somme = add_complexe_float(somme, mult_complexe_float(((complexe_float_t *)A)[i * lda + k], ((complexe_float_t *)B)[j * ldb + k]));
                        }
                        ((complexe_float_t *)C)[i * ldc + j].real = ((complexe_float_t*)alpha)->real * somme.real - ((complexe_float_t*)alpha)->imaginary * somme.imaginary + ((complexe_float_t*)beta)->real * ((complexe_float_t *)C)[i * ldc + j].real - ((complexe_float_t*)beta)->imaginary * ((complexe_float_t *)C)[i * ldc + j].imaginary;
                        ((complexe_float_t *)C)[i * ldc + j].imaginary = ((complexe_float_t*)alpha)->real * somme.imaginary + ((complexe_float_t*)alpha)->imaginary * somme.real + ((complexe_float_t*)beta)->real * ((complexe_float_t *)C)[i * ldc + j].imaginary + ((complexe_float_t*)beta)->imaginary * ((complexe_float_t *)C)[i * ldc + j].real;
                    }
                }
            }
        } else if (TransA == MNCblasTrans) {
            if(TransB == MNCblasNoTrans) {
                complexe_float_t somme;
                for(int i = 0; i < M; i++) {
                    for(int j = 0; j < N; j++) {
                        somme.real = 0;
                        somme.imaginary = 0;
                        for(int k = 0; k < K; k++) {
                            somme = add_complexe_float(somme, mult_complexe_float(((complexe_float_t *)A)[k * lda + i], ((complexe_float_t *)B)[k * ldb + j]));
                        }
                        ((complexe_float_t *)C)[i * ldc + j].real = ((complexe_float_t*)alpha)->real * somme.real + ((complexe_float_t*)alpha)->imaginary * somme.imaginary + ((complexe_float_t*)beta)->real * ((complexe_float_t *)C)[i * ldc + j].real - ((complexe_float_t*)beta)->imaginary * ((complexe_float_t *)C)[i * ldc + j].imaginary;
                        ((complexe_float_t *)C)[i * ldc + j].imaginary = ((complexe_float_t*)alpha)->real * somme.imaginary - ((complexe_float_t*)alpha)->imaginary * somme.real + ((complexe_float_t*)beta)->real * ((complexe_float_t *)C)[i * ldc + j].imaginary + ((complexe_float_t*)beta)->imaginary * ((complexe_float_t *)C)[i * ldc + j].real;
                    }
                }
            } else if(TransB == MNCblasTrans) {
                complexe_float_t somme;
                for(int i = 0; i < M; i++) {
                    for(int j = 0; j < N; j++) {
                        somme.real = 0;
                        somme.imaginary = 0;
                        for(int k = 0; k < K; k++) {
                            somme = add_complexe_float(somme, mult_complexe_float(((complexe_float_t *)A)[k * lda + i], ((complexe_float_t *)B)[j * ldb + k]));
                        }
                        ((complexe_float_t *)C)[i * ldc + j].real = ((complexe_float_t*)alpha)->real * somme.real + ((complexe_float_t*)alpha)->imaginary * somme.imaginary + ((complexe_float_t*)beta)->real * ((complexe_float_t *)C)[i * ldc + j].real - ((complexe_float_t*)beta)->imaginary * ((complexe_float_t *)C)[i * ldc + j].imaginary;
                        ((complexe_float_t *)C)[i * ldc + j].imaginary = ((complexe_float_t*)alpha)->real * somme.imaginary + ((complexe_float_t*)alpha)->imaginary * somme.real + ((complexe_float_t*)beta)->real * ((complexe_float_t *)C)[i * ldc + j].imaginary + ((complexe_float_t*)beta)->imaginary * ((complexe_float_t *)C)[i * ldc + j].real;
                    }
                }
            } else {
                complexe_float_t somme;
                for(int i = 0; i < M; i++) {
                    for(int j = 0; j < N; j++) {
                        somme.real = 0;
                        somme.imaginary = 0;
                        for(int k = 0; k < K; k++) {
                            somme = add_complexe_float(somme, mult_complexe_float(((complexe_float_t *)A)[k * lda + i], ((complexe_float_t *)B)[j * ldb + k]));
                        }
                        ((complexe_float_t *)C)[i * ldc + j].real = ((complexe_float_t*)alpha)->real * somme.real - ((complexe_float_t*)alpha)->imaginary * somme.imaginary + ((complexe_float_t*)beta)->real * ((complexe_float_t *)C)[i * ldc + j].real - ((complexe_float_t*)beta)->imaginary * ((complexe_float_t *)C)[i * ldc + j].imaginary;
                        ((complexe_float_t *)C)[i * ldc + j].imaginary = ((complexe_float_t*)alpha)->real * somme.imaginary + ((complexe_float_t*)alpha)->imaginary * somme.real + ((complexe_float_t*)beta)->real * ((complexe_float_t *)C)[i * ldc + j].imaginary + ((complexe_float_t*)beta)->imaginary * ((complexe_float_t *)C)[i * ldc + j].real;
                    }
                }
            }
        }
    } else {
        if(TransA == MNCblasNoTrans) {
            if(TransB == MNCblasNoTrans) {
                mncblas_zgemm(layout, MNCblasTrans, MNCblasTrans, N, M, K, alpha, B, ldb, A, lda, beta, C, ldc);
            } else {
                mncblas_zgemm(layout, MNCblasTrans, MNCblasNoTrans, N, M, K, alpha, B, ldb, A, lda, beta, C, ldc);
            }
        } else {
            if(TransB == MNCblasNoTrans) {
                mncblas_zgemm(layout, MNCblasNoTrans, MNCblasTrans, M, N, K, alpha, A, lda, B, ldb, beta, C, ldc);
            } else {
                mncblas_zgemm(layout, MNCblasNoTrans, MNCblasNoTrans, M, N, K, alpha, A, lda, B, ldb, beta, C, ldc);
            }
        }
    }
}