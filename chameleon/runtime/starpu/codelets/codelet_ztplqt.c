/**
 *
 * @file starpu/codelet_ztplqt.c
 *
 * @copyright 2009-2016 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2021 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon ztplqt StarPU codelet
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
static void cl_ztplqt_cpu_func(void *descr[], void *cl_arg)
{
    int M;
    int N;
    int L;
    int ib;
    CHAM_tile_t *tileA;
    CHAM_tile_t *tileB;
    CHAM_tile_t *tileT;
    CHAM_tile_t *tileWORK;

    tileA    = cti_interface_get(descr[0]);
    tileB    = cti_interface_get(descr[1]);
    tileT    = cti_interface_get(descr[2]);
    tileWORK = cti_interface_get(descr[3]); /* ib * nb */
    starpu_codelet_unpack_args( cl_arg, &M, &N, &L, &ib );

    TCORE_zlaset( ChamUpperLower, ib, M, 0., 0., tileT );
    TCORE_ztplqt( M, N, L, ib,
                 tileA, tileB, tileT, tileWORK->mat );
}
#endif /* !defined(CHAMELEON_SIMULATION) */

/*
 * Codelet definition
 */
CODELETS_CPU(ztplqt, cl_ztplqt_cpu_func)

void INSERT_TASK_ztplqt( const RUNTIME_option_t *options,
                         int M, int N, int L, int ib, int nb,
                         const CHAM_desc_t *A, int Am, int An,
                         const CHAM_desc_t *B, int Bm, int Bn,
                         const CHAM_desc_t *T, int Tm, int Tn )
{
    struct starpu_codelet *codelet = &cl_ztplqt;
    void (*callback)(void*) = options->profiling ? cl_ztplqt_callback : NULL;
    starpu_option_request_t* schedopt = (starpu_option_request_t *)(options->request->schedopt);
    int workerid = (schedopt == NULL) ? -1 : schedopt->workerid;

    CHAMELEON_BEGIN_ACCESS_DECLARATION;
    CHAMELEON_ACCESS_RW(A, Am, An);
    CHAMELEON_ACCESS_RW(B, Bm, Bn);
    CHAMELEON_ACCESS_W(T, Tm, Tn);
    CHAMELEON_END_ACCESS_DECLARATION;

    rt_starpu_insert_task(
        codelet,
        STARPU_VALUE, &M,     sizeof(int),
        STARPU_VALUE, &N,     sizeof(int),
        STARPU_VALUE, &L,     sizeof(int),
        STARPU_VALUE, &ib,    sizeof(int),
        STARPU_RW,     RTBLKADDR(A, CHAMELEON_Complex64_t, Am, An),
        STARPU_RW,     RTBLKADDR(B, CHAMELEON_Complex64_t, Bm, Bn),
        STARPU_W,      RTBLKADDR(T, CHAMELEON_Complex64_t, Tm, Tn),
        /* Other options */
        STARPU_SCRATCH,   options->ws_worker,
        STARPU_PRIORITY,  options->priority,
        STARPU_CALLBACK,  callback,
        STARPU_EXECUTE_ON_WORKER, workerid,
#if defined(CHAMELEON_USE_MPI)
        STARPU_EXECUTE_ON_NODE, B->get_rankof(B, Bm, Bn),
#endif
#if defined(CHAMELEON_CODELETS_HAVE_NAME)
        STARPU_NAME, (L == 0) ? "ztplqs" : "ztplqt",
#endif
        0);

    (void)ib; (void)nb;
}
