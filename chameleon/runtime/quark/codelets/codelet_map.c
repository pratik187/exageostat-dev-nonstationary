/**
 *
 * @file quark/codelet_map.c
 *
 * @copyright 2018-2021 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon map Quark codelet
 *
 * @version 1.0.0
 * @author Mathieu Faverge
 * @date 2020-03-03
 *
 */
#include "chameleon_quark.h"
#include "chameleon/tasks.h"

void CORE_map_quark(Quark *quark)
{
    const CHAM_desc_t *desc;
    cham_uplo_t uplo;
    int m;
    int n;
    CHAM_tile_t *tile;
    cham_unary_operator_t op_fct;
    void *op_args;

    quark_unpack_args_7( quark, desc, uplo, m, n, tile, op_fct, op_args );
    op_fct( desc, uplo, m, n, tile, op_args );
}

void INSERT_TASK_map( const RUNTIME_option_t *options,
                      cham_uplo_t uplo, const CHAM_desc_t *A, int Am, int An,
                      cham_unary_operator_t op_fct, void *op_args )
{
    quark_option_t *opt = (quark_option_t*)(options->schedopt);

    QUARK_Insert_Task(
        opt->quark, CORE_map_quark, (Quark_Task_Flags*)opt,
        sizeof(CHAM_desc_t*),             &A,    VALUE,
        sizeof(cham_uplo_t),              &uplo, VALUE,
        sizeof(int),                      &Am,   VALUE,
        sizeof(int),                      &An,   VALUE,
        sizeof(void*), RTBLKADDR(A, void, Am, An), INOUT,
        sizeof(cham_unary_operator_t),    &op_fct,  VALUE,
        sizeof(void*),                    &op_args, VALUE,
        0);
}
