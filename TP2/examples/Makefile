all: test_complexe3 test_complexe2 test_complexe test_dot test_dot_dyn test_copy test_copy_dyn test_swap test_axpy test_asum test_iamin test_iamax test_nrm2 test_gemv test_gemm

INC_DIR=../include

LIB_DIR=../lib

LIBST=-lmnblas
LIBDYN=-lmnblasdyn

OPTIONS_COMPIL  =-Wall -O2 -fPIC -I$(INC_DIR) -lm

OPTIONS_LINK_STATIC  =  -fopenmp -L$(LIB_DIR) $(LIBST)
OPTIONS_LINK_DYN  = -fopenmp -L$(LIB_DIR) $(LIBDYN)  

test_dot: test_dot.o flop.o $(LIB_DIR)/libmnblas.a
	gcc -o test_dot test_dot.o flop.o $(OPTIONS_LINK_STATIC)

test_swap: test_swap.o flop.o $(LIB_DIR)/libmnblas.a
	gcc -o test_swap test_swap.o flop.o $(OPTIONS_LINK_STATIC)

test_axpy: test_axpy.o flop.o $(LIB_DIR)/libmnblas.a
	gcc -o test_axpy test_axpy.o flop.o $(OPTIONS_LINK_STATIC)

test_asum: test_asum.o flop.o $(LIB_DIR)/libmnblas.a
	gcc -o test_asum test_asum.o flop.o $(OPTIONS_LINK_STATIC)

test_iamin: test_iamin.o flop.o $(LIB_DIR)/libmnblas.a
	gcc -o test_iamin test_iamin.o flop.o $(OPTIONS_LINK_STATIC)

test_iamax: test_iamax.o flop.o $(LIB_DIR)/libmnblas.a
	gcc -o test_iamax test_iamax.o flop.o $(OPTIONS_LINK_STATIC)

test_nrm2: test_nrm2.o flop.o $(LIB_DIR)/libmnblas.a
	gcc -o test_nrm2 test_nrm2.o flop.o $(OPTIONS_LINK_STATIC) -lm

test_copy: test_copy.o flop.o $(LIB_DIR)/libmnblas.a
	gcc -o test_copy test_copy.o flop.o $(OPTIONS_LINK_STATIC)

test_gemv: test_gemv.o flop.o $(LIB_DIR)/libmnblas.a
	gcc -o test_gemv test_gemv.o flop.o $(OPTIONS_LINK_STATIC)

test_gemm: test_gemm.o flop.o $(LIB_DIR)/libmnblas.a
	gcc -o test_gemm test_gemm.o flop.o $(OPTIONS_LINK_STATIC)

flop.o: flop.c flop.h
	gcc $(OPTIONS_COMPIL) -c flop.c

test_dot_dyn: test_dot.o flop.o
	gcc -o test_dot_dyn flop.o test_dot.o $(OPTIONS_LINK_DYN) -lm

test_dot.o: test_dot.c $(INC_DIR)/mnblas.h
	gcc $(OPTIONS_COMPIL) -c test_dot.c

test_copy_dyn: test_copy.o flop.o
	gcc -o test_copy_dyn flop.o test_copy.o $(OPTIONS_LINK_DYN) -lm

test_copy.o: test_copy.c $(INC_DIR)/mnblas.h
	gcc $(OPTIONS_COMPIL) -c test_copy.c

test_swap.o: test_swap.c $(INC_DIR)/mnblas.h
	gcc $(OPTIONS_COMPIL) -c test_swap.c

test_axpy.o: test_axpy.c $(INC_DIR)/mnblas.h
	gcc $(OPTIONS_COMPIL) -c test_axpy.c

test_iamin.o: test_iamin.c $(INC_DIR)/mnblas.h
	gcc $(OPTIONS_COMPIL) -c test_iamin.c

test_iamax.o: test_iamax.c $(INC_DIR)/mnblas.h
	gcc $(OPTIONS_COMPIL) -c test_iamax.c

test_asum.o: test_asum.c $(INC_DIR)/mnblas.h
	gcc $(OPTIONS_COMPIL) -c test_asum.c

test_nrm2.o: test_nrm2.c $(INC_DIR)/mnblas.h
	gcc $(OPTIONS_COMPIL) -c test_nrm2.c

test_gemv.o: test_gemv.c $(INC_DIR)/mnblas.h
	gcc $(OPTIONS_COMPIL) -c test_gemv.c

test_gemm.o: test_gemm.c $(INC_DIR)/mnblas.h
	gcc $(OPTIONS_COMPIL) -c test_gemm.c

test_complexe.o: test_complexe.c  $(INC_DIR)/complexe.h
	gcc $(OPTIONS_COMPIL) -c test_complexe.c

test_complexe: test_complexe.o flop.o
	gcc -o test_complexe test_complexe.o flop.o $(OPTIONS_LINK_STATIC)

test_complexe2.o: test_complexe2.c  $(INC_DIR)/complexe2.h
	gcc $(OPTIONS_COMPIL) -c test_complexe2.c

test_complexe2: test_complexe2.o flop.o
	gcc -o test_complexe2 test_complexe2.o flop.o $(OPTIONS_LINK_STATIC)


test_complexe3.o: test_complexe3.c  $(INC_DIR)/complexe2.h
	gcc $(OPTIONS_COMPIL) -c test_complexe3.c

test_complexe3: test_complexe3.o flop.o
	gcc -o test_complexe3 test_complexe3.o flop.o $(OPTIONS_LINK_STATIC)


clean:
	rm -f *.o test_dot test_dot_dyn test_complexe test_complexe2 test_complexe3 *~
