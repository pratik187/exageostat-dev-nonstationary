/**
 *
 * @file cuda_zherfb.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2021 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon cuda_zherfb GPU kernel
 *
 * @version 1.0.0
 * @author Florent Pruvost
 * @author Mathieu Faverge
 * @date 2020-03-03
 * @precisions normal z -> c d s
 *
 */
#include "cudablas.h"

int
CUDA_zherfb( cham_uplo_t uplo, int n,
             int k, int ib, int nb,
             const cuDoubleComplex *A, int lda,
             const cuDoubleComplex *T, int ldt,
             cuDoubleComplex *C, int ldc,
             cuDoubleComplex *WORK, int ldwork,
             CUBLAS_STREAM_PARAM )
{
    /* Check input arguments */
    if ((uplo != ChamUpper) && (uplo != ChamLower)) {
        cudablas_error(1, "Illegal value of uplo");
        return -1;
    }
    if (n < 0) {
        cudablas_error(2, "Illegal value of n");
        return -2;
    }
    if (k < 0) {
        cudablas_error(3, "Illegal value of k");
        return -3;
    }
    if (ib < 0) {
        cudablas_error(4, "Illegal value of ib");
        return -4;
    }
    if (nb < 0) {
        cudablas_error(5, "Illegal value of nb");
        return -5;
    }
    if ( (lda < chameleon_max(1,n)) && (n > 0) ) {
        cudablas_error(7, "Illegal value of lda");
        return -7;
    }
    if ( (ldt < chameleon_max(1,ib)) && (ib > 0) ) {
        cudablas_error(9, "Illegal value of ldt");
        return -9;
    }
    if ( (ldc < chameleon_max(1,n)) && (n > 0) ) {
        cudablas_error(11, "Illegal value of ldc");
        return -11;
    }

    if (uplo == ChamLower) {
        /* Left */
        CUDA_zunmqrt( ChamLeft, ChamConjTrans, n, n, k, ib,
                      A, lda, T, ldt, C, ldc, WORK, ldwork,
                      CUBLAS_STREAM_VALUE );
        /* Right */
        CUDA_zunmqrt( ChamRight, ChamNoTrans, n, n, k, ib,
                      A, lda, T, ldt, C, ldc, WORK, ldwork,
                      CUBLAS_STREAM_VALUE );
    }
    else {
        /* Right */
        CUDA_zunmlqt( ChamRight, ChamConjTrans, n, n, k, ib,
                      A, lda, T, ldt, C, ldc, WORK, ldwork,
                      CUBLAS_STREAM_VALUE );
        /* Left */
        CUDA_zunmlqt( ChamLeft, ChamNoTrans, n, n, k, ib,
                      A, lda, T, ldt, C, ldc, WORK, ldwork,
                      CUBLAS_STREAM_VALUE );
    }
    return 0;
}
