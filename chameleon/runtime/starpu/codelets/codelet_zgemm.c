/**
 *
 * @file starpu/codelet_zgemm.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2021 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zgemm StarPU codelet
 *
 * @version 1.1.0
 * @author Hatem Ltaief
 * @author Jakub Kurzak
 * @author Mathieu Faverge
 * @author Emmanuel Agullo
 * @author Cedric Castagnede
 * @author Lucas Barros de Assis
 * @author Florent Pruvost
 * @date 2021-03-16
 * @precisions normal z -> c d s
 *
 */
#include "chameleon_starpu.h"
#include "runtime_codelet_z.h"

struct cl_zgemm_args_s {
    cham_trans_t transA;
    cham_trans_t transB;
    int m;
    int n;
    int k;
    CHAMELEON_Complex64_t alpha;
    CHAM_tile_t *tileA;
    CHAM_tile_t *tileB;
    CHAMELEON_Complex64_t beta;
    CHAM_tile_t *tileC;
};

#if !defined(CHAMELEON_SIMULATION)
static void
cl_zgemm_cpu_func( void *descr[], void *cl_arg )
{
    struct cl_zgemm_args_s clargs;
    CHAM_tile_t *tileA;
    CHAM_tile_t *tileB;
    CHAM_tile_t *tileC;

    tileA = cti_interface_get(descr[0]);
    tileB = cti_interface_get(descr[1]);
    tileC = cti_interface_get(descr[2]);

    starpu_codelet_unpack_args( cl_arg, &clargs );
    TCORE_zgemm( clargs.transA, clargs.transB,
                 clargs.m, clargs.n, clargs.k,
                 clargs.alpha, tileA, tileB,
                 clargs.beta,  tileC );
}

#ifdef CHAMELEON_USE_CUDA
static void
cl_zgemm_cuda_func( void *descr[], void *_cl_arg )
{
    struct cl_zgemm_args_s clargs;
    CHAM_tile_t *tileA;
    CHAM_tile_t *tileB;
    CHAM_tile_t *tileC;

    tileA = cti_interface_get(descr[0]);
    tileB = cti_interface_get(descr[1]);
    tileC = cti_interface_get(descr[2]);

    starpu_codelet_unpack_args( _cl_arg, &clargs );

    RUNTIME_getStream( stream );

    CUDA_zgemm(
        clargs.transA, clargs.transB,
        clargs.m, clargs.n, clargs.k,
        (cuDoubleComplex*)&(clargs.alpha),
        tileA->mat, tileA->ld,
        tileB->mat, tileB->ld,
        (cuDoubleComplex*)&(clargs.beta),
        tileC->mat, tileC->ld,
        stream );

#ifndef STARPU_CUDA_ASYNC
    cudaStreamSynchronize( stream );
#endif

    return;
}
#endif /* defined(CHAMELEON_USE_CUDA) */
#endif /* !defined(CHAMELEON_SIMULATION) */

/*
 * Codelet definition
 */
CODELETS( zgemm, cl_zgemm_cpu_func, cl_zgemm_cuda_func, STARPU_CUDA_ASYNC )

void INSERT_TASK_zgemm( const RUNTIME_option_t *options,
                        cham_trans_t transA, cham_trans_t transB,
                        int m, int n, int k, int nb,
                        CHAMELEON_Complex64_t alpha, const CHAM_desc_t *A, int Am, int An,
                                                     const CHAM_desc_t *B, int Bm, int Bn,
                        CHAMELEON_Complex64_t beta,  const CHAM_desc_t *C, int Cm, int Cn )
{
    if ( alpha == 0. ) {
        return INSERT_TASK_zlascal( options, ChamUpperLower, m, n, nb,
                                    beta, C, Cm, Cn );
    }

    struct cl_zgemm_args_s clargs = {
        .transA = transA,
        .transB = transB,
        .m      = m,
        .n      = n,
        .k      = k,
        .alpha  = alpha,
        .tileA  = A->get_blktile( A, Am, An ),
        .tileB  = B->get_blktile( B, Bm, Bn ),
        .beta   = beta,
        .tileC  = C->get_blktile( C, Cm, Cn )
    };
    void (*callback)(void*);
    RUNTIME_request_t       *request  = options->request;
    starpu_option_request_t *schedopt = (starpu_option_request_t *)(request->schedopt);
    int                      workerid, accessC;
    char                    *cl_name = "zgemm";

    /* Handle cache */
    CHAMELEON_BEGIN_ACCESS_DECLARATION;
    CHAMELEON_ACCESS_R(A, Am, An);
    CHAMELEON_ACCESS_R(B, Bm, Bn);
    CHAMELEON_ACCESS_RW(C, Cm, Cn);
    CHAMELEON_END_ACCESS_DECLARATION;

    /* Callback fro profiling information */
    callback = options->profiling ? cl_zgemm_callback : NULL;

    /* Fix the worker id */
    workerid = (schedopt == NULL) ? -1 : schedopt->workerid;

    /* Reduce the C access if needed */
    accessC = ( beta == 0. ) ? STARPU_W : STARPU_RW;

    /* Insert the task */
    rt_starpu_insert_task(
        &cl_zgemm,
        /* Task codelet arguments */
        STARPU_VALUE, &clargs, sizeof(struct cl_zgemm_args_s),
        STARPU_R,      RTBLKADDR(A, CHAMELEON_Complex64_t, Am, An),
        STARPU_R,      RTBLKADDR(B, CHAMELEON_Complex64_t, Bm, Bn),
        accessC,       RTBLKADDR(C, CHAMELEON_Complex64_t, Cm, Cn),

        /* Common task arguments */
        STARPU_PRIORITY,          options->priority,
        STARPU_CALLBACK,          callback,
        STARPU_EXECUTE_ON_WORKER, workerid,
#if defined(CHAMELEON_CODELETS_HAVE_NAME)
        STARPU_NAME,              cl_name,
#endif

        0 );

    (void)nb;
}
