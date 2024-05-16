/**
 *
 * @file openmp/codelet_zsyrk.c
 *
 * @copyright 2012-2021 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zsyrk OpenMP codelet
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

void INSERT_TASK_zsyrk( const RUNTIME_option_t *options,
                      cham_uplo_t uplo, cham_trans_t trans,
                      int n, int k, int nb,
                      CHAMELEON_Complex64_t alpha, const CHAM_desc_t *A, int Am, int An,
                      CHAMELEON_Complex64_t beta, const CHAM_desc_t *C, int Cm, int Cn )
{
    ( void )nb;
    CHAM_tile_t *tileA = A->get_blktile( A, Am, An );
    CHAM_tile_t *tileC = C->get_blktile( C, Cm, Cn );
#pragma omp task firstprivate( uplo, trans, n, k, alpha, tileA, beta, tileC ) depend( in:tileA[0] ) depend( inout:tileC[0] )
    TCORE_zsyrk( uplo, trans,
        n, k,
        alpha, tileA,
        beta, tileC );
}
