/**
 *
 * @file parsec/codelet_ztrsm.c
 *
 * @copyright 2009-2015 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2021 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon ztrsm PaRSEC codelet
 *
 * @version 1.0.0
 * @author Reazul Hoque
 * @author Mathieu Faverge
 * @date 2020-03-03
 * @precisions normal z -> c d s
 *
 */
#include "chameleon_parsec.h"
#include "chameleon/tasks_z.h"
#include "coreblas/coreblas_z.h"

static inline int
CORE_ztrsm_parsec( parsec_execution_stream_t *context,
                    parsec_task_t             *this_task )
{
    cham_side_t side, uplo, trans, diag;
    int tempmm, nb, ldak, ldam;
    CHAMELEON_Complex64_t alpha;
    CHAMELEON_Complex64_t *T;
    CHAMELEON_Complex64_t *C;

    parsec_dtd_unpack_args(
        this_task, &side, &uplo, &trans, &diag, &tempmm, &nb, &alpha, &T, &ldak, &C, &ldam );

    CORE_ztrsm( side, uplo, trans, diag,
                tempmm, nb, alpha,
                T, ldak, C, ldam );

    (void)context;
    return PARSEC_HOOK_RETURN_DONE;
}

void INSERT_TASK_ztrsm(const RUNTIME_option_t *options,
                      cham_side_t side, cham_uplo_t uplo, cham_trans_t transA, cham_diag_t diag,
                      int m, int n, int nb,
                      CHAMELEON_Complex64_t alpha, const CHAM_desc_t *A, int Am, int An,
                      const CHAM_desc_t *B, int Bm, int Bn)
{
    parsec_taskpool_t* PARSEC_dtd_taskpool = (parsec_taskpool_t *)(options->sequence->schedopt);
    CHAM_tile_t *tileA = A->get_blktile( A, Am, An );
    CHAM_tile_t *tileB = B->get_blktile( B, Bm, Bn );

    parsec_dtd_taskpool_insert_task(
        PARSEC_dtd_taskpool, CORE_ztrsm_parsec, options->priority, "Trsm",
        sizeof(int),    &side,                     VALUE,
        sizeof(int),    &uplo,                     VALUE,
        sizeof(int),    &transA,                   VALUE,
        sizeof(int),    &diag,                     VALUE,
        sizeof(int),           &m,                        VALUE,
        sizeof(int),           &n,                        VALUE,
        sizeof(CHAMELEON_Complex64_t),           &alpha,      VALUE,
        PASSED_BY_REF,     RTBLKADDR( A, CHAMELEON_Complex64_t, Am, An ), chameleon_parsec_get_arena_index( A ) | INPUT,
        sizeof(int), &(tileA->ld), VALUE,
        PASSED_BY_REF,     RTBLKADDR( B, CHAMELEON_Complex64_t, Bm, Bn ), chameleon_parsec_get_arena_index( B ) | INOUT | AFFINITY,
        sizeof(int), &(tileB->ld), VALUE,
        PARSEC_DTD_ARG_END );

    (void)nb;
}
