/**
 *
 * @file quark/codelet_zgetrf_nopiv.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2021 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zgetrf_nopiv Quark codelet
 *
 * @version 1.0.0
 * @author Omar Zenati
 * @author Mathieu Faverge
 * @author Emmanuel Agullo
 * @author Cedric Castagnede
 * @author Florent Pruvost
 * @date 2020-03-03
 * @precisions normal z -> c d s
 *
 */
#include "chameleon_quark.h"
#include "chameleon/tasks_z.h"
#include "coreblas/coreblas_ztile.h"

void CORE_zgetrf_nopiv_quark(Quark *quark)
{
    int m;
    int n;
    int ib;
    CHAM_tile_t *tileA;
    RUNTIME_sequence_t *sequence;
    RUNTIME_request_t *request;
    int iinfo;
    int info;

    quark_unpack_args_7(quark, m, n, ib, tileA, sequence, request, iinfo);
    TCORE_zgetrf_nopiv(m, n, ib, tileA, &info);
    if ( info != CHAMELEON_SUCCESS ) {
        RUNTIME_sequence_flush( (CHAM_context_t*)quark, sequence, request, iinfo+info );
    }
}

void INSERT_TASK_zgetrf_nopiv(const RUNTIME_option_t *options,
                             int m, int n, int ib, int nb,
                             const CHAM_desc_t *A, int Am, int An,
                             int iinfo)
{
    quark_option_t *opt = (quark_option_t*)(options->schedopt);
    DAG_CORE_GETRF;
    QUARK_Insert_Task(
        opt->quark, CORE_zgetrf_nopiv_quark, (Quark_Task_Flags*)opt,
        sizeof(int),                        &m,             VALUE,
        sizeof(int),                        &n,             VALUE,
        sizeof(int),                        &ib,            VALUE,
        sizeof(void*), RTBLKADDR(A, CHAMELEON_Complex64_t, Am, An),                     INOUT,
        sizeof(RUNTIME_sequence_t*),           &(options->sequence),      VALUE,
        sizeof(RUNTIME_request_t*),            &(options->request),       VALUE,
        sizeof(int),                        &iinfo,         VALUE,
        0);
}
