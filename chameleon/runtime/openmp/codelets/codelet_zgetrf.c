/**
 *
 * @file openmp/codelet_zgetrf.c
 *
 * @copyright 2012-2021 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zgetrf OpenMP codelet
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

void INSERT_TASK_zgetrf( const RUNTIME_option_t *options,
                       int m, int n, int nb,
                       const CHAM_desc_t *A, int Am, int An,
                       int *IPIV,
                       cham_bool_t check_info, int iinfo )
{
    CHAM_tile_t *tileA = A->get_blktile( A, Am, An );
    int info = 0;
#pragma omp task firstprivate( m, n, tileA, IPIV, info ) depend( out:IPIV[0] ) depend( inout:tileA[0] )
    TCORE_zgetrf( m, n, tileA, IPIV, &info );
}
