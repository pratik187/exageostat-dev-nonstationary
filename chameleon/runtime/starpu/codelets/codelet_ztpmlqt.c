/**
 *
 * @file starpu/codelet_ztpmlqt.c
 *
 * @copyright 2009-2016 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2021 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * @brief Chameleon ztpmlqt StarPU codelet
 *
 * @version 1.1.0
 * @author Mathieu Faverge
 * @author Lucas Barros de Assis
 * @author Florent Pruvost
 * @date 2021-01-11
 * @precisions normal z -> s d c
 *
 */
#include "chameleon_starpu.h"
#include "runtime_codelet_z.h"

#if !defined(CHAMELEON_SIMULATION)
static void cl_ztpmlqt_cpu_func(void *descr[], void *cl_arg)
{
    cham_side_t side;
    cham_trans_t trans;
    int M;
    int N;
    int K;
    int L;
    int ib;
    size_t lwork;
    CHAM_tile_t *tileV;
    CHAM_tile_t *tileT;
    CHAM_tile_t *tileA;
    CHAM_tile_t *tileB;
    CHAM_tile_t *tileW;

    tileV = cti_interface_get(descr[0]);
    tileT = cti_interface_get(descr[1]);
    tileA = cti_interface_get(descr[2]);
    tileB = cti_interface_get(descr[3]);
    tileW = cti_interface_get(descr[4]); /* ib * nb */
    starpu_codelet_unpack_args( cl_arg, &side, &trans, &M, &N, &K, &L, &ib, &lwork );

    TCORE_ztpmlqt( side, trans, M, N, K, L, ib,
                   tileV, tileT, tileA, tileB, tileW->mat );
}

#if defined(CHAMELEON_USE_CUDA)
static void cl_ztpmlqt_cuda_func(void *descr[], void *cl_arg)
{
    cham_side_t side;
    cham_trans_t trans;
    int M;
    int N;
    int K;
    int L;
    int ib;
    size_t lwork;
    CHAM_tile_t *tileV;
    CHAM_tile_t *tileT;
    CHAM_tile_t *tileA;
    CHAM_tile_t *tileB;
    CHAM_tile_t *tileW;

    tileV = cti_interface_get(descr[0]);
    tileT = cti_interface_get(descr[1]);
    tileA = cti_interface_get(descr[2]);
    tileB = cti_interface_get(descr[3]);
    tileW = cti_interface_get(descr[4]); /* 3*ib*nb */

    starpu_codelet_unpack_args( cl_arg, &side, &trans, &M, &N, &K, &L, &ib, &lwork );

    RUNTIME_getStream(stream);

    CUDA_ztpmlqt(
            side, trans, M, N, K, L, ib,
            tileV->mat, tileV->ld,
            tileT->mat, tileT->ld,
            tileA->mat, tileA->ld,
            tileB->mat, tileB->ld,
            tileW->mat, lwork, stream );

#ifndef STARPU_CUDA_ASYNC
    cudaStreamSynchronize( stream );
#endif
}
#endif /* defined(CHAMELEON_USE_CUDA) */
#endif /* !defined(CHAMELEON_SIMULATION) */

/*
 * Codelet definition
 */
CODELETS(ztpmlqt, cl_ztpmlqt_cpu_func, cl_ztpmlqt_cuda_func, STARPU_CUDA_ASYNC)

void INSERT_TASK_ztpmlqt( const RUNTIME_option_t *options,
                          cham_side_t side, cham_trans_t trans,
                          int M, int N, int K, int L, int ib, int nb,
                          const CHAM_desc_t *V, int Vm, int Vn,
                          const CHAM_desc_t *T, int Tm, int Tn,
                          const CHAM_desc_t *A, int Am, int An,
                          const CHAM_desc_t *B, int Bm, int Bn )
{
    struct starpu_codelet *codelet = &cl_ztpmlqt;
    void (*callback)(void*) = options->profiling ? cl_ztpmlqt_callback : NULL;
    starpu_option_request_t* schedopt = (starpu_option_request_t *)(options->request->schedopt);
    int workerid = (schedopt == NULL) ? -1 : schedopt->workerid;

    CHAMELEON_BEGIN_ACCESS_DECLARATION;
    CHAMELEON_ACCESS_R(V, Vm, Vn);
    CHAMELEON_ACCESS_R(T, Tm, Tn);
    CHAMELEON_ACCESS_RW(A, Am, An);
    CHAMELEON_ACCESS_RW(B, Bm, Bn);
    CHAMELEON_END_ACCESS_DECLARATION;

    rt_starpu_insert_task(
        codelet,
        STARPU_VALUE, &side,  sizeof(int),
        STARPU_VALUE, &trans, sizeof(int),
        STARPU_VALUE, &M,     sizeof(int),
        STARPU_VALUE, &N,     sizeof(int),
        STARPU_VALUE, &K,     sizeof(int),
        STARPU_VALUE, &L,     sizeof(int),
        STARPU_VALUE, &ib,     sizeof(int),
        STARPU_VALUE, &(options->ws_wsize), sizeof(size_t),
        STARPU_R,      RTBLKADDR(V, CHAMELEON_Complex64_t, Vm, Vn),
        STARPU_R,      RTBLKADDR(T, CHAMELEON_Complex64_t, Tm, Tn),
        STARPU_RW,     RTBLKADDR(A, CHAMELEON_Complex64_t, Am, An),
        STARPU_RW,     RTBLKADDR(B, CHAMELEON_Complex64_t, Bm, Bn),
        /* Other options */
        STARPU_SCRATCH,   options->ws_worker,
        STARPU_PRIORITY,  options->priority,
        STARPU_CALLBACK,  callback,
        STARPU_EXECUTE_ON_WORKER, workerid,
#if defined(CHAMELEON_USE_MPI)
        STARPU_EXECUTE_ON_NODE, B->get_rankof(B, Bm, Bn),
#endif
#if defined(CHAMELEON_CODELETS_HAVE_NAME)
        STARPU_NAME, (( L == 0 ) ? "ztsmlq" : "ztpmlqt"),
#endif
        0);

    (void)ib; (void)nb;
}
