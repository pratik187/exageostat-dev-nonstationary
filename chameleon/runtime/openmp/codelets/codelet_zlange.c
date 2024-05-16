/**
 *
 * @file openmp/codelet_zlange.c
 *
 * @copyright 2012-2021 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zlange OpenMP codelet
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

void INSERT_TASK_zlange( const RUNTIME_option_t *options,
                         cham_normtype_t norm, int M, int N, int NB,
                         const CHAM_desc_t *A, int Am, int An,
                         const CHAM_desc_t *B, int Bm, int Bn )
{
    CHAM_tile_t *tileA = A->get_blktile( A, Am, An );
    CHAM_tile_t *tileB = B->get_blktile( B, Bm, Bn );
    int ws_size = options->ws_wsize;
#pragma omp task firstprivate( ws_size, M, N, tileA, tileB, options ) depend( in:tileA[0] ) depend( inout:tileB[0] )
    {
        double work[ws_size];
        TCORE_zlange( norm, M, N, tileA, work, tileB->mat );
    }
}

void INSERT_TASK_zlange_max( const RUNTIME_option_t *options,
                           const CHAM_desc_t *A, int Am, int An,
                           const CHAM_desc_t *B, int Bm, int Bn )
{
    CHAM_tile_t *tileA = A->get_blktile( A, Am, An );
    CHAM_tile_t *tileB = B->get_blktile( B, Bm, Bn );

#pragma omp task firstprivate( tileA, tileB ) depend( in:tileA[0] ) depend( inout:tileB[0] )
    {
        double *A, *B;

        A = tileA->mat;
        B = tileB->mat;

        if ( A[0] > B[0] ) {
            B[0] = A[0];
        }
    }
}
