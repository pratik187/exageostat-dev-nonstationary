/**
 *
 * @file starpu/codelet_zher2k.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2021 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zher2k StarPU codelet
 *
 * @version 1.1.0
 * @author Hatem Ltaief
 * @author Jakub Kurzak
 * @author Mathieu Faverge
 * @author Emmanuel Agullo
 * @author Cedric Castagnede
 * @author Lucas Barros de Assis
 * @author Florent Pruvost
 * @date 2021-01-11
 * @precisions normal z -> c
 *
 */
#include "chameleon_starpu.h"
#include "runtime_codelet_z.h"

#if !defined(CHAMELEON_SIMULATION)
static void cl_zher2k_cpu_func(void *descr[], void *cl_arg)
{
    cham_uplo_t uplo;
    cham_trans_t trans;
    int n;
    int k;
    CHAMELEON_Complex64_t alpha;
    CHAM_tile_t *tileA;
    CHAM_tile_t *tileB;
    double beta;
    CHAM_tile_t *tileC;

    tileA = cti_interface_get(descr[0]);
    tileB = cti_interface_get(descr[1]);
    tileC = cti_interface_get(descr[2]);

    starpu_codelet_unpack_args(cl_arg, &uplo, &trans, &n, &k, &alpha, &beta);
    TCORE_zher2k(uplo, trans,
                n, k, alpha, tileA, tileB, beta, tileC);
}

#ifdef CHAMELEON_USE_CUDA
static void cl_zher2k_cuda_func(void *descr[], void *cl_arg)
{
    cham_uplo_t uplo;
    cham_trans_t trans;
    int n;
    int k;
    cuDoubleComplex alpha;
    CHAM_tile_t *tileA;
    CHAM_tile_t *tileB;
    double beta;
    CHAM_tile_t *tileC;

    tileA = cti_interface_get(descr[0]);
    tileB = cti_interface_get(descr[1]);
    tileC = cti_interface_get(descr[2]);

    starpu_codelet_unpack_args(cl_arg, &uplo, &trans, &n, &k, &alpha, &beta);

    RUNTIME_getStream(stream);

    CUDA_zher2k( uplo, trans,
                 n, k,
                 &alpha, tileA->mat, tileA->ld,
                         tileB->mat, tileB->ld,
                 &beta,  tileC->mat, tileC->ld,
                 stream);

#ifndef STARPU_CUDA_ASYNC
    cudaStreamSynchronize( stream );
#endif

    return;
}
#endif /* CHAMELEON_USE_CUDA */
#endif /* !defined(CHAMELEON_SIMULATION) */

/*
 * Codelet definition
 */
CODELETS(zher2k, cl_zher2k_cpu_func, cl_zher2k_cuda_func, STARPU_CUDA_ASYNC)

/**
 *
 * @ingroup INSERT_TASK_Complex64_t
 *
 */
void
INSERT_TASK_zher2k( const RUNTIME_option_t *options,
                    cham_uplo_t uplo, cham_trans_t trans,
                    int n, int k, int nb,
                    CHAMELEON_Complex64_t alpha, const CHAM_desc_t *A, int Am, int An,
                                                 const CHAM_desc_t *B, int Bm, int Bn,
                    double beta,                 const CHAM_desc_t *C, int Cm, int Cn )
{
    if ( alpha == 0. ) {
        return INSERT_TASK_zlascal( options, uplo, n, n, nb,
                                    beta, C, Cm, Cn );
    }

    (void)nb;
    struct starpu_codelet *codelet = &cl_zher2k;
    void (*callback)(void*) = options->profiling ? cl_zher2k_callback : NULL;
    starpu_option_request_t* schedopt = (starpu_option_request_t *)(options->request->schedopt);
    int workerid = (schedopt == NULL) ? -1 : schedopt->workerid;
    int accessC = ( beta == 0. ) ? STARPU_W : STARPU_RW;

    CHAMELEON_BEGIN_ACCESS_DECLARATION;
    CHAMELEON_ACCESS_R(A, Am, An);
    CHAMELEON_ACCESS_R(B, Bm, Bn);
    CHAMELEON_ACCESS_RW(C, Cm, Cn);
    CHAMELEON_END_ACCESS_DECLARATION;

    rt_starpu_insert_task(
        codelet,
        STARPU_VALUE,      &uplo,                sizeof(int),
        STARPU_VALUE,     &trans,                sizeof(int),
        STARPU_VALUE,         &n,                        sizeof(int),
        STARPU_VALUE,         &k,                        sizeof(int),
        STARPU_VALUE,     &alpha,         sizeof(CHAMELEON_Complex64_t),
        STARPU_R,                 RTBLKADDR(A, CHAMELEON_Complex64_t, Am, An),
        STARPU_R,                 RTBLKADDR(B, CHAMELEON_Complex64_t, Bm, Bn),
        STARPU_VALUE,      &beta,                     sizeof(double),
        accessC,                  RTBLKADDR(C, CHAMELEON_Complex64_t, Cm, Cn),
        STARPU_PRIORITY,    options->priority,
        STARPU_CALLBACK,    callback,
        STARPU_EXECUTE_ON_WORKER, workerid,
#if defined(CHAMELEON_CODELETS_HAVE_NAME)
        STARPU_NAME, "zher2k",
#endif
        0);
}
