/**
 *
 * @file openmp/codelet_ztpmqrt.c
 *
 * @copyright 2012-2021 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * @brief Chameleon ztpmqrt OpenMP codelet
 *
 * @version 1.0.0
 * @author Philippe Virouleau
 * @author Mathieu Faverge
 * @date 2020-03-03
 * @precisions normal z -> s d c
 *
 */
#include "chameleon_openmp.h"
#include "chameleon/tasks_z.h"
#include "coreblas/coreblas_ztile.h"

void INSERT_TASK_ztpmqrt( const RUNTIME_option_t *options,
                          cham_side_t side, cham_trans_t trans,
                          int M, int N, int K, int L, int ib, int nb,
                          const CHAM_desc_t *V, int Vm, int Vn,
                          const CHAM_desc_t *T, int Tm, int Tn,
                          const CHAM_desc_t *A, int Am, int An,
                          const CHAM_desc_t *B, int Bm, int Bn )
{
    CHAM_tile_t *tileA = A->get_blktile( A, Am, An );
    CHAM_tile_t *tileB = B->get_blktile( B, Bm, Bn );
    CHAM_tile_t *tileT = T->get_blktile( T, Tm, Tn );
    CHAM_tile_t *tileV = V->get_blktile( V, Vm, Vn );
    int ws_size = options->ws_wsize;

#pragma omp task firstprivate( ws_size, side, trans, M, N, K, L, ib, nb, tileV, tileT, tileA, tileB ) depend( in:tileV[0], tileT[0] ) depend( inout:tileA[0], tileB[0] )
    {
        CHAMELEON_Complex64_t tmp[ws_size];
        TCORE_ztpmqrt( side, trans, M, N, K, L, ib,
                      tileV, tileT, tileA, tileB, tmp );
    }
}
