/**
 *
 * @file quark/runtime_control.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2021 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon Quark control routines
 *
 * @version 1.1.0
 * @author Vijay Joshi
 * @author Cedric Castagnede
 * @author Florent Pruvost
 * @author Mathieu Faverge
 * @author Samuel Thibault
 * @date 2020-04-22
 *
 */
#include <stdio.h>
#include <stdlib.h>
#include "chameleon_quark.h"

/**
 *
 */
int RUNTIME_init( CHAM_context_t *chamctxt,
                  int ncpus,
                  int ncudas,
                  int nthreads_per_worker )
{
    int hres = -1;
    if ( ncudas > 0 ) {
        chameleon_warning( "RUNTIME_init_scheduler(quark)", "GPUs are not supported for now");
    }

    if ( nthreads_per_worker > 0 ) {
        chameleon_warning( "RUNTIME_init_scheduler(quark)", "Multi-threaded kernels are not supported for now");
    }

    chamctxt->schedopt = (void*)QUARK_New( ncpus );

    if(NULL != chamctxt->schedopt) {
        chamctxt->nworkers = ncpus;
        chamctxt->nthreads_per_worker = nthreads_per_worker;
        hres = 0;
    }

    return hres;
}

/**
 *
 */
void RUNTIME_finalize( CHAM_context_t *chamctxt )
{
    QUARK_Delete((Quark*)(chamctxt->schedopt));
    return;
}

/**
 *  To suspend the processing of new tasks by workers
 */
void RUNTIME_pause( CHAM_context_t *chamctxt )
{
    (void)chamctxt;
    return;
}

/**
 *  This is the symmetrical call to RUNTIME_pause,
 *  used to resume the workers polling for new tasks.
 */
void RUNTIME_resume( CHAM_context_t *chamctxt )
{
    (void)chamctxt;
    return;
}

/**
 *  Busy-waiting barrier
 */
void RUNTIME_barrier( CHAM_context_t *chamctxt )
{
    QUARK_Barrier((Quark*)(chamctxt->schedopt));
}

/**
 *  Display a progress information when executing the tasks
 */
void RUNTIME_progress( CHAM_context_t *chamctxt )
{
    (void)chamctxt;
    return;
}

/**
 * Thread rank.
 */
int RUNTIME_thread_rank( CHAM_context_t *chamctxt )
{
    return QUARK_Thread_Rank((Quark*)(chamctxt->schedopt));
}

/**
 * Number of threads.
 */
int RUNTIME_thread_size( CHAM_context_t *chamctxt )
{
    (void)chamctxt;
    /*
     * TODO: should add a function to Quark to get the number of thread from the
     * data structure and not from the system function
     */
    return quark_get_numthreads();
}

/**
 *  The process rank
 */
int RUNTIME_comm_rank( CHAM_context_t *chamctxt )
{
    (void)chamctxt;
    return 0;
}

/**
 *  This returns the size of the distributed computation
 */
int RUNTIME_comm_size( CHAM_context_t *chamctxt )
{
    (void)chamctxt;
    return 1;
}
