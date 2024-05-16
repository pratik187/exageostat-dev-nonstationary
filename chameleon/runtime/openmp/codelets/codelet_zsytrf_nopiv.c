/**
 *
 * @file openmp/codelet_zsytrf_nopiv.c
 *
 * @copyright 2012-2021 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zsytrf_nopiv OpenMP codelet
 *
 * @version 1.0.0
 * @author Philippe Virouleau
 * @author Mathieu Faverge
 * @date 2020-03-03
 * @precisions normal z -> c
 *
 */
#include "chameleon_openmp.h"
#include "chameleon/tasks_z.h"
#include "coreblas/coreblas_ztile.h"

void INSERT_TASK_zsytrf_nopiv( const RUNTIME_option_t *options,
                               cham_uplo_t uplo, int n, int nb,
                               const CHAM_desc_t *A, int Am, int An,
                               int iinfo )
{
    CHAM_tile_t *tileA = A->get_blktile( A, Am, An );
#pragma omp task firstprivate( uplo, n, tileA ) depend( inout:tileA[0] )
    TCORE_zsytf2_nopiv( uplo, n, tileA );
}
