/**
 *
 * @file parsec/codelet_map.c
 *
 * @copyright 2018-2021 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon map PaRSEC codelet
 *
 * @version 1.0.0
 * @author Mathieu Faverge
 * @date 2020-03-03
 *
 */
#include "chameleon_parsec.h"
#include "chameleon/tasks.h"

static inline int
CORE_map_parsec( parsec_execution_stream_t *context,
                 parsec_task_t             *this_task )
{
    const CHAM_desc_t *desc;
    cham_uplo_t uplo;
    int m;
    int n;
    void *data;
    cham_unary_operator_t op_fct;
    void *op_args;

    parsec_dtd_unpack_args(
        this_task, &desc, &uplo, &m, &n, &data, &op_fct, &op_args );
    op_fct( desc, uplo, m, n, data, op_args );

    (void)context;
    return PARSEC_HOOK_RETURN_DONE;
}

void INSERT_TASK_map( const RUNTIME_option_t *options,
                      cham_uplo_t uplo, const CHAM_desc_t *A, int Am, int An,
                      cham_unary_operator_t op_fct, void *op_args )
{
    parsec_taskpool_t* PARSEC_dtd_taskpool = (parsec_taskpool_t *)(options->sequence->schedopt);

    parsec_dtd_taskpool_insert_task(
        PARSEC_dtd_taskpool, CORE_map_parsec, options->priority, "map",
        sizeof(CHAM_desc_t*),             &A,    VALUE,
        sizeof(cham_uplo_t),              &uplo, VALUE,
        sizeof(int),                      &Am,   VALUE,
        sizeof(int),                      &An,   VALUE,
        PASSED_BY_REF, RTBLKADDR(A, void, Am, An), chameleon_parsec_get_arena_index( A ) | INOUT,
        sizeof(cham_unary_operator_t),    &op_fct,  VALUE,
        sizeof(void*),                    &op_args, VALUE,
        PARSEC_DTD_ARG_END );
}
