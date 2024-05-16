/**
 *
 * @file starpu/codelet_zaxpy.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2021 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zaxpy StarPU codelet
 *
 * @version 1.1.0
 * @author Florent Pruvost
 * @author Mathieu Faverge
 * @author Samuel Thibault
 * @date 2021-01-11
 * @precisions normal z -> c d s
 *
 */
#include "chameleon_starpu.h"
#include "runtime_codelet_z.h"

#if !defined(CHAMELEON_SIMULATION)
static void cl_zaxpy_cpu_func(void *descr[], void *cl_arg)
{
    int M;
    CHAMELEON_Complex64_t alpha;
    CHAM_tile_t *tileA;
    int incA;
    CHAM_tile_t *tileB;
    int incB;

    tileA = cti_interface_get(descr[0]);
    tileB = cti_interface_get(descr[1]);
    starpu_codelet_unpack_args(cl_arg, &M, &alpha, &incA, &incB);
    TCORE_zaxpy(M, alpha, tileA, incA, tileB, incB);
}
#endif /* !defined(CHAMELEON_SIMULATION) */

/*
 * Codelet definition
 */
CODELETS_CPU(zaxpy, cl_zaxpy_cpu_func)

void INSERT_TASK_zaxpy( const RUNTIME_option_t *options,
                        int M, CHAMELEON_Complex64_t alpha,
                        const CHAM_desc_t *A, int Am, int An, int incA,
                        const CHAM_desc_t *B, int Bm, int Bn, int incB )
{
    struct starpu_codelet *codelet = &cl_zaxpy;
    void (*callback)(void*) = options->profiling ? cl_zaxpy_callback : NULL;
    starpu_option_request_t* schedopt = (starpu_option_request_t *)(options->request->schedopt);
    int workerid = (schedopt == NULL) ? -1 : schedopt->workerid;

    if ( alpha == 0. ) {
        return;
    }

    CHAMELEON_BEGIN_ACCESS_DECLARATION;
    CHAMELEON_ACCESS_R(A, Am, An);
    CHAMELEON_ACCESS_RW(B, Bm, Bn);
    CHAMELEON_END_ACCESS_DECLARATION;

    rt_starpu_insert_task(
            codelet,
            STARPU_VALUE,    &M,                           sizeof(int),
            STARPU_VALUE,    &alpha,                       sizeof(CHAMELEON_Complex64_t),
            STARPU_R,        RTBLKADDR(A, CHAMELEON_Complex64_t, Am, An),
            STARPU_VALUE,    &incA,                        sizeof(int),
            STARPU_RW,       RTBLKADDR(B, CHAMELEON_Complex64_t, Bm, Bn),
            STARPU_VALUE,    &incB,                        sizeof(int),
            STARPU_PRIORITY, options->priority,
            STARPU_CALLBACK, callback,
            STARPU_EXECUTE_ON_WORKER, workerid,
#if defined(CHAMELEON_CODELETS_HAVE_NAME)
            STARPU_NAME, "zaxpy",
#endif
            0);
}
