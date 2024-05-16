/**
 *
 * @file starpu/codelet_zunmqr.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2021 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zunmqr StarPU codelet
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
static void cl_zunmqr_cpu_func(void *descr[], void *cl_arg)
{
    cham_side_t side;
    cham_trans_t trans;
    int m;
    int n;
    int k;
    int ib;
    CHAM_tile_t *tileA;
    CHAM_tile_t *tileT;
    CHAM_tile_t *tileC;
    CHAM_tile_t *tileW;
    int ldW;

    tileA = cti_interface_get(descr[0]);
    tileT = cti_interface_get(descr[1]);
    tileC = cti_interface_get(descr[2]);
    tileW = cti_interface_get(descr[3]); /* ib * nb */

    starpu_codelet_unpack_args( cl_arg, &side, &trans, &m, &n, &k, &ib, &ldW );

    TCORE_zunmqr( side, trans, m, n, k, ib,
                  tileA, tileT, tileC, tileW->mat, ldW );
}

#if defined(CHAMELEON_USE_CUDA)
static void cl_zunmqr_cuda_func(void *descr[], void *cl_arg)
{
    cham_side_t side;
    cham_trans_t trans;
    int m;
    int n;
    int k;
    int ib;
    CHAM_tile_t *tileA;
    CHAM_tile_t *tileT;
    CHAM_tile_t *tileC;
    CHAM_tile_t *tileW;
    int ldW;

    tileA = cti_interface_get(descr[0]);
    tileT = cti_interface_get(descr[1]);
    tileC = cti_interface_get(descr[2]);
    tileW = cti_interface_get(descr[3]); /* ib * nb */

    starpu_codelet_unpack_args( cl_arg, &side, &trans, &m, &n, &k, &ib, &ldW );

    RUNTIME_getStream(stream);

    CUDA_zunmqrt(
            side, trans, m, n, k, ib,
            tileA->mat, tileA->ld,
            tileT->mat, tileT->ld,
            tileC->mat, tileC->ld,
            tileW->mat, ldW, stream );

#ifndef STARPU_CUDA_ASYNC
    cudaStreamSynchronize( stream );
#endif
}
#endif /* defined(CHAMELEON_USE_CUDA) */
#endif /* !defined(CHAMELEON_SIMULATION) */

/*
 * Codelet definition
 */
CODELETS(zunmqr, cl_zunmqr_cpu_func, cl_zunmqr_cuda_func, STARPU_CUDA_ASYNC)

void INSERT_TASK_zunmqr( const RUNTIME_option_t *options,
                         cham_side_t side, cham_trans_t trans,
                         int m, int n, int k, int ib, int nb,
                         const CHAM_desc_t *A, int Am, int An,
                         const CHAM_desc_t *T, int Tm, int Tn,
                         const CHAM_desc_t *C, int Cm, int Cn )
{
    struct starpu_codelet *codelet = &cl_zunmqr;
    void (*callback)(void*) = options->profiling ? cl_zunmqr_callback : NULL;
    starpu_option_request_t* schedopt = (starpu_option_request_t *)(options->request->schedopt);
    int workerid = (schedopt == NULL) ? -1 : schedopt->workerid;

    CHAMELEON_BEGIN_ACCESS_DECLARATION;
    CHAMELEON_ACCESS_R(A, Am, An);
    CHAMELEON_ACCESS_R(T, Tm, Tn);
    CHAMELEON_ACCESS_RW(C, Cm, Cn);
    CHAMELEON_END_ACCESS_DECLARATION;

    rt_starpu_insert_task(
        codelet,
        STARPU_VALUE,    &side,              sizeof(int),
        STARPU_VALUE,    &trans,             sizeof(int),
        STARPU_VALUE,    &m,                 sizeof(int),
        STARPU_VALUE,    &n,                 sizeof(int),
        STARPU_VALUE,    &k,                 sizeof(int),
        STARPU_VALUE,    &ib,                sizeof(int),
        STARPU_R,         RTBLKADDR(A, CHAMELEON_Complex64_t, Am, An),
        STARPU_R,         RTBLKADDR(T, CHAMELEON_Complex64_t, Tm, Tn),
        STARPU_RW,        RTBLKADDR(C, CHAMELEON_Complex64_t, Cm, Cn),
        /* ib * nb */
        STARPU_SCRATCH,   options->ws_worker,
        STARPU_VALUE,    &nb,                sizeof(int),
        STARPU_PRIORITY,  options->priority,
        STARPU_CALLBACK,  callback,
        STARPU_EXECUTE_ON_WORKER, workerid,
#if defined(CHAMELEON_CODELETS_HAVE_NAME)
        STARPU_NAME, "zunmqr",
#endif
        0);
}
