/**
 *
 * @file openmp/codelet_zsyr2k.c
 *
 * @copyright 2012-2021 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zsyr2k OpenMP codelet
 *
 * @version 1.0.0
 * @author Philippe Virouleau
 * @author Mathieu Faverge
 * @date 2020-03-03
 * @precisions normal z -> c d s
 *
 */
#include "chameleon_openmp.h"
#include "chameleon/tasks_z.h"
#include "coreblas/coreblas_ztile.h"

void INSERT_TASK_zsyr2k( const RUNTIME_option_t *options,
                       cham_uplo_t uplo, cham_trans_t trans,
                       int n, int k, int nb,
                       CHAMELEON_Complex64_t alpha, const CHAM_desc_t *A, int Am, int An,
                       const CHAM_desc_t *B, int Bm, int Bn,
                       CHAMELEON_Complex64_t beta, const CHAM_desc_t *C, int Cm, int Cn )
{
    ( void )nb;
    CHAM_tile_t *tileA = A->get_blktile( A, Am, An );
    CHAM_tile_t *tileB = B->get_blktile( B, Bm, Bn );
    CHAM_tile_t *tileC = C->get_blktile( C, Cm, Cn );
#pragma omp task firstprivate( uplo, trans, n, k, alpha, tileA, tileB, beta, tileC ) depend( in:tileA[0], tileB[0] ) depend( inout:tileC[0] )
    TCORE_zsyr2k( uplo, trans,
                 n, k, alpha, tileA, tileB, beta, tileC );
}
