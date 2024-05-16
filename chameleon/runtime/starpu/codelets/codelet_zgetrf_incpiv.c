/**
 *
 * @file starpu/codelet_zgetrf_incpiv.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2021 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zgetrf_incpiv StarPU codelet
 *
 * @version 1.1.0
 * @comment This file has been automatically generated
 *          from Plasma 2.5.0 for CHAMELEON 0.9.2
 * @author Hatem Ltaief
 * @author Jakub Kurzak
 * @author Mathieu Faverge
 * @author Emmanuel Agullo
 * @author Cedric Castagnede
 * @author Lucas Barros de Assis
 * @author Florent Pruvost
 * @author Samuel Thibault
 * @date 2021-01-11
 * @precisions normal z -> c d s
 *
 */
#include "chameleon_starpu.h"
#include "runtime_codelet_z.h"

#if !defined(CHAMELEON_SIMULATION)
static void cl_zgetrf_incpiv_cpu_func(void *descr[], void *cl_arg)
{
    CHAMELEON_starpu_ws_t *h_work;
    int m;
    int n;
    int ib;
    CHAM_tile_t *tileA;
    int *IPIV;
    cham_bool_t check_info;
    int iinfo;
    RUNTIME_sequence_t *sequence;
    RUNTIME_request_t *request;
    int info = 0;

    tileA = cti_interface_get(descr[0]);

    starpu_codelet_unpack_args(cl_arg, &m, &n, &ib, &IPIV, &check_info, &iinfo, &h_work, &sequence, &request);
    TCORE_zgetrf_incpiv(m, n, ib, tileA, IPIV, &info);

    if ( (sequence->status == CHAMELEON_SUCCESS) && (info != 0) ) {
        RUNTIME_sequence_flush( NULL, sequence, request, iinfo+info );
    }
}
#endif /* !defined(CHAMELEON_SIMULATION) */

/*
 * Codelet definition
 */
CODELETS_CPU(zgetrf_incpiv, cl_zgetrf_incpiv_cpu_func)

void INSERT_TASK_zgetrf_incpiv(const RUNTIME_option_t *options,
                              int m, int n, int ib, int nb,
                              const CHAM_desc_t *A, int Am, int An,
                              const CHAM_desc_t *L, int Lm, int Ln,
                              int *IPIV,
                              cham_bool_t check_info, int iinfo)
{
    (void)nb;
    struct starpu_codelet *codelet = &cl_zgetrf_incpiv;
    void (*callback)(void*) = options->profiling ? cl_zgetrf_incpiv_callback : NULL;
    starpu_option_request_t* schedopt = (starpu_option_request_t *)(options->request->schedopt);
    int workerid = (schedopt == NULL) ? -1 : schedopt->workerid;

    CHAMELEON_starpu_ws_t *h_work = (CHAMELEON_starpu_ws_t*)(options->ws_host);

    CHAMELEON_BEGIN_ACCESS_DECLARATION;
    CHAMELEON_ACCESS_RW(A, Am, An);
    CHAMELEON_ACCESS_W(L, Lm, Ln);
    CHAMELEON_END_ACCESS_DECLARATION;

    rt_starpu_insert_task(
        codelet,
        STARPU_VALUE,    &m,                 sizeof(int),
        STARPU_VALUE,    &n,                 sizeof(int),
        STARPU_VALUE,    &ib,                sizeof(int),
        STARPU_RW,        RTBLKADDR(A, CHAMELEON_Complex64_t, Am, An),
        STARPU_W,         RTBLKADDR(L, CHAMELEON_Complex64_t, Lm, Ln),
        STARPU_VALUE,    &IPIV,              sizeof(int*),
        STARPU_VALUE,    &check_info,        sizeof(cham_bool_t),
        STARPU_VALUE,    &iinfo,             sizeof(int),
        STARPU_SCRATCH,   options->ws_worker,
        STARPU_VALUE,    &h_work,            sizeof(CHAMELEON_starpu_ws_t *),
        STARPU_VALUE,    &(options->sequence),       sizeof(RUNTIME_sequence_t*),
        STARPU_VALUE,    &(options->request),        sizeof(RUNTIME_request_t*),
        STARPU_PRIORITY,  options->priority,
        STARPU_CALLBACK,  callback,
        STARPU_EXECUTE_ON_WORKER, workerid,
#if defined(CHAMELEON_CODELETS_HAVE_NAME)
        STARPU_NAME, "zgetrf_incpiv",
#endif
        0);
}
