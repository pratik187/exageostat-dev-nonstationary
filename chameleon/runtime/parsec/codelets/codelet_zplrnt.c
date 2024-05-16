/**
 *
 * @file parsec/codelet_zplrnt.c
 *
 * @copyright 2009-2015 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2021 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zplrnt PaRSEC codelet
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
CORE_zplrnt_parsec( parsec_execution_stream_t *context,
                    parsec_task_t             *this_task )
{
    int m;
    int n;
    CHAMELEON_Complex64_t *A;
    int lda;
    int bigM;
    int m0;
    int n0;
    unsigned long long int seed;

    parsec_dtd_unpack_args(
        this_task, &m, &n, &A, &lda, &bigM, &m0, &n0, &seed );

    CORE_zplrnt( m, n, A, lda, bigM, m0, n0, seed );

    (void)context;
    return PARSEC_HOOK_RETURN_DONE;
}

void INSERT_TASK_zplrnt( const RUNTIME_option_t *options,
                        int m, int n, const CHAM_desc_t *A, int Am, int An,
                        int bigM, int m0, int n0, unsigned long long int seed )
{
    parsec_taskpool_t* PARSEC_dtd_taskpool = (parsec_taskpool_t *)(options->sequence->schedopt);
    CHAM_tile_t *tileA = A->get_blktile( A, Am, An );

    parsec_dtd_taskpool_insert_task(
        PARSEC_dtd_taskpool, CORE_zplrnt_parsec, options->priority, "zplrnt",
        sizeof(int),       &m,                          VALUE,
        sizeof(int),       &n,                          VALUE,
        PASSED_BY_REF,     RTBLKADDR( A, CHAMELEON_Complex64_t, Am, An ), chameleon_parsec_get_arena_index( A ) | OUTPUT | AFFINITY,
        sizeof(int), &(tileA->ld), VALUE,
        sizeof(int),       &bigM,                       VALUE,
        sizeof(int),       &m0,                         VALUE,
        sizeof(int),       &n0,                         VALUE,
        sizeof(unsigned long long int),       &seed,    VALUE,
        PARSEC_DTD_ARG_END );
}
