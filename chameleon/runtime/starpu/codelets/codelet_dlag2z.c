/**
 *
 * @file starpu/codelet_dlag2z.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2021 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon dlag2z StarPU codelet
 *
 * @version 1.1.0
 * @author Mathieu Faverge
 * @date 2021-01-11
 * @precisions normal z -> c
 *
 */
#include "chameleon_starpu.h"
#include "runtime_codelet_z.h"

#if !defined(CHAMELEON_SIMULATION)
static void cl_dlag2z_cpu_func(void *descr[], void *cl_arg)
{
    cham_uplo_t uplo;
    int m;
    int n;
    CHAM_tile_t *tileA;
    CHAM_tile_t *tileB;

    tileA = cti_interface_get(descr[0]);
    tileB = cti_interface_get(descr[1]);

    starpu_codelet_unpack_args(cl_arg, &uplo, &m, &n);
    TCORE_dlag2z( uplo, m, n, tileA, tileB );
}
#endif /* !defined(CHAMELEON_SIMULATION) */

/*
 * Codelet definition
 */
CODELETS_CPU(dlag2z, cl_dlag2z_cpu_func)

/**
 *
 * @ingroup INSERT_TASK_Complex64_t
 *
 */
void INSERT_TASK_dlag2z( const RUNTIME_option_t *options,
                         cham_uplo_t uplo, int m, int n,
                         const CHAM_desc_t *A, int Am, int An,
                         const CHAM_desc_t *B, int Bm, int Bn )
{
    struct starpu_codelet *codelet = &cl_dlag2z;
    void (*callback)(void*) = options->profiling ? cl_dlag2z_callback : NULL;
    starpu_option_request_t* schedopt = (starpu_option_request_t *)(options->request->schedopt);
    int workerid = (schedopt == NULL) ? -1 : schedopt->workerid;

    CHAMELEON_BEGIN_ACCESS_DECLARATION;
    CHAMELEON_ACCESS_R(A, Am, An);
    CHAMELEON_ACCESS_W(B, Bm, Bn);
    CHAMELEON_END_ACCESS_DECLARATION;

    rt_starpu_insert_task(
        codelet,
        STARPU_VALUE,    &uplo,              sizeof(uplo),
        STARPU_VALUE,    &m,                 sizeof(int),
        STARPU_VALUE,    &n,                 sizeof(int),
        STARPU_R,         RTBLKADDR(A, double, Am, An),
        STARPU_W,         RTBLKADDR(B, CHAMELEON_Complex64_t, Bm, Bn),
        STARPU_PRIORITY,  options->priority,
        STARPU_CALLBACK,  callback,
        STARPU_EXECUTE_ON_WORKER, workerid,
#if defined(CHAMELEON_CODELETS_HAVE_NAME)
        STARPU_NAME, "dlag2z",
#endif
        0);
}
