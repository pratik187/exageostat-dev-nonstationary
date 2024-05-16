/**
 *
 * @file starpu/runtime_async.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2021 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon StarPU asynchronous routines
 *
 * @version 1.1.0
 * @author Mathieu Faverge
 * @author Cedric Castagnede
 * @author Florent Pruvost
 * @author Samuel Thibault
 * @date 2020-04-22
 *
 */
#include <stdlib.h>
#include "chameleon_starpu.h"

/**
 *  Create a sequence
 */
int RUNTIME_sequence_create( CHAM_context_t  *chamctxt,
                             RUNTIME_sequence_t *sequence )
{
    (void)chamctxt;
    (void)sequence;
    return CHAMELEON_SUCCESS;
}

/**
 *  Destroy a sequence
 */
int RUNTIME_sequence_destroy( CHAM_context_t  *chamctxt,
                              RUNTIME_sequence_t *sequence )
{
    (void)chamctxt;
    (void)sequence;
    return CHAMELEON_SUCCESS;
}

/**
 *  Wait for the completion of a sequence
 */
int RUNTIME_sequence_wait( CHAM_context_t  *chamctxt,
                           RUNTIME_sequence_t *sequence )
{
    (void)chamctxt;
    (void)sequence;

    if (chamctxt->progress_enabled) {
        RUNTIME_progress(chamctxt);
    }

#if defined(CHAMELEON_USE_MPI)
#  if defined(HAVE_STARPU_MPI_WAIT_FOR_ALL)
    starpu_mpi_wait_for_all(MPI_COMM_WORLD);
#  else
    starpu_task_wait_for_all();
    starpu_mpi_barrier(MPI_COMM_WORLD);
#  endif
#else
    starpu_task_wait_for_all();
#endif
    return CHAMELEON_SUCCESS;
}

/**
 *  Terminate a sequence
 */
void RUNTIME_sequence_flush( CHAM_context_t  *chamctxt,
                             RUNTIME_sequence_t *sequence,
                             RUNTIME_request_t  *request,
                             int status )
{
    (void)chamctxt;
    sequence->request = request;
    sequence->status = status;
    request->status = status;
    return;
}

/**
 *  Create a request
 */
int RUNTIME_request_create( CHAM_context_t  *chamctxt,
                            RUNTIME_request_t *request )
{
    (void)chamctxt;
    /* allocate schedopt */
    request->schedopt = (starpu_option_request_t*)malloc(sizeof(starpu_option_request_t));
    /* initialize schedopt */
    starpu_option_request_t* schedopt = (starpu_option_request_t *)(request->schedopt);
    /* default is to not use "execute_on_a_specific_worker" i.e. -1 */
    schedopt->workerid = -1;
    request->status = CHAMELEON_SUCCESS;
    return CHAMELEON_SUCCESS;
}

/**
 *  Destroy a request
 */
int RUNTIME_request_destroy( CHAM_context_t  *chamctxt,
                             RUNTIME_request_t *request )
{
    (void)chamctxt;
    free(request->schedopt);
    return CHAMELEON_SUCCESS;
}

/**
 *  Set runtime parameter for a request
 */
int RUNTIME_request_set( CHAM_context_t  *chamctxt,
                         RUNTIME_request_t *request,
                         int param, int value )
{
    if ( request->schedopt == NULL ) {
        chameleon_error("RUNTIME_request_set", "request not initialized");
        return CHAMELEON_ERR_NOT_INITIALIZED;
    }
    starpu_option_request_t* schedopt = (starpu_option_request_t *)(request->schedopt);

    switch ( param ) {
        case CHAMELEON_REQUEST_WORKERID:
            if ( (value < -1) || (value >= chamctxt->nworkers) ) {
                chameleon_error("RUNTIME_request_set", "workerid should be in [-1, NCPUS-1]");
                return CHAMELEON_ERR_ILLEGAL_VALUE;
            }
            schedopt->workerid = value;
            break;
        default:
            chameleon_error("RUNTIME_request_set", "unknown parameter");
            return CHAMELEON_ERR_ILLEGAL_VALUE;
    }

    return CHAMELEON_SUCCESS;
}