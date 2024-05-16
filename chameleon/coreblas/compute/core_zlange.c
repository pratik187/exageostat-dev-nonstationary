/**
 *
 * @file core_zlange.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2021 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon core_zlange CPU kernel
 *
 * @version 1.0.0
 * @comment This file has been automatically generated
 *          from Plasma 2.6.0 for CHAMELEON 0.9.2
 * @author Julien Langou
 * @author Henricus Bouwmeester
 * @author Mathieu Faverge
 * @author Florent Pruvost
 * @date 2020-03-03
 * @precisions normal z -> c d s
 *
 */
#include "coreblas/lapacke.h"
#include "coreblas.h"

/**
 *
 * @ingroup CORE_CHAMELEON_Complex64_t
 *
 *  CORE_PLASMA_zlange returns the value
 *
 *     zlange = ( max(abs(A(i,j))), NORM = ChamMaxNorm
 *              (
 *              ( norm1(A),         NORM = ChamOneNorm
 *              (
 *              ( normI(A),         NORM = ChamInfNorm
 *              (
 *              ( normF(A),         NORM = ChamFrobeniusNorm
 *
 *  where norm1 denotes the one norm of a matrix (maximum column sum),
 *  normI denotes the infinity norm of a matrix (maximum row sum) and
 *  normF denotes the Frobenius norm of a matrix (square root of sum
 *  of squares). Note that max(abs(A(i,j))) is not a consistent matrix
 *  norm.
 *
 *******************************************************************************
 *
 * @param[in] norm
 *          = ChamMaxNorm: Max norm
 *          = ChamOneNorm: One norm
 *          = ChamInfNorm: Infinity norm
 *          = ChamFrobeniusNorm: Frobenius norm
 *
 * @param[in] M
 *          The number of rows of the matrix A. M >= 0.
 *
 * @param[in] N
 *          The number of columns of the matrix A. N >= 0.
 *
 * @param[in] A
 *          The M-by-N matrix A.
 *
 * @param[in] LDA
 *          The leading dimension of the array A. LDA >= max(1,M).
 *
 * @param[in,out] work
 *          Array of dimension (MAX(1,LWORK)), where LWORK >= M when NORM =
 *          ChamInfNorm; otherwise, WORK is not referenced.
 *
 * @param[out] normA
 *          On exit, normA is the norm of matrix A.
 *
 */

void CORE_zlange(cham_normtype_t norm, int M, int N,
                 const CHAMELEON_Complex64_t *A, int LDA,
                 double *work, double *normA)
{
    *normA = LAPACKE_zlange_work(
        LAPACK_COL_MAJOR,
        chameleon_lapack_const(norm),
        M, N, A, LDA, work);
}
