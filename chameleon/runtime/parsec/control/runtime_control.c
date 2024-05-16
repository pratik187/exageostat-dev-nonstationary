/**
 *
 * @file parsec/runtime_control.c
 *
 * @copyright 2012-2017 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2021 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon PaRSEC control routines
 *
 * @version 1.0.0
 * @author Reazul Hoque
 * @author Mathieu Faverge
 * @author Samuel Thibault
 * @date 2020-03-03
 *
 */
#include <stdio.h>
#include <stdlib.h>
#include "chameleon_parsec.h"

#if defined(CHAMELEON_USE_MPI)
#include <mpi.h>
#endif

/**
 * Initialize CHAMELEON
 */
int RUNTIME_init( CHAM_context_t *chamctxt,
                  int ncpus,
                  int ncudas,
                  int nthreads_per_worker )
{
    int hres = -1, default_ncores = -1;
    int *argc = (int *)malloc(sizeof(int));
    *argc = 0;

    /* Initializing parsec context */
    if( 0 < ncpus ) {
        default_ncores = ncpus;
    }
    chamctxt->parallel_enabled = CHAMELEON_TRUE;
    chamctxt->schedopt = (void *)parsec_init(default_ncores, argc, NULL);

    if(NULL != chamctxt->schedopt) {
        chamctxt->nworkers = ncpus;
        chamctxt->nthreads_per_worker = nthreads_per_worker;
        hres = 0;
    }

    free(argc);

    (void)ncudas;
    return hres;
}

/**
 * Finalize CHAMELEON
 */
void RUNTIME_finalize( CHAM_context_t *chamctxt )
{
    parsec_context_t *parsec = (parsec_context_t*)chamctxt->schedopt;
    parsec_fini(&parsec);
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
 * Barrier CHAMELEON.
 */
void RUNTIME_barrier( CHAM_context_t *chamctxt )
{
    parsec_context_t *parsec = (parsec_context_t*)(chamctxt->schedopt);
    // This will be a problem with the fake tasks inserted to detect end of DTD algorithms
    parsec_context_wait( parsec );
    return;
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
    (void)chamctxt;
    return 0;
}

/**
 * Thread rank.
 */
int RUNTIME_thread_size( CHAM_context_t *chamctxt )
{
    // TODO: fixme
    //return vpmap_get_nb_total_threads();
    (void)chamctxt;
    return 1;
}

/**
 *  This returns the rank of this process
 */
int RUNTIME_comm_rank( CHAM_context_t *chamctxt )
{
    int rank = 0;
#if defined(CHAMELEON_USE_MPI)
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif

    (void)chamctxt;
    return rank;
}

/**
 *  This returns the size of the distributed computation
 */
int RUNTIME_comm_size( CHAM_context_t *chamctxt )
{
    int size = 1;
#if defined(CHAMELEON_USE_MPI)
    MPI_Comm_size(MPI_COMM_WORLD, &size);
#endif

    (void)chamctxt;
    return size;
}
