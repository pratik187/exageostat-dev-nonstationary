/**
 *
 * @file workspace.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2021 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon workspace routines
 *
 * @version 1.0.0
 * @author Mathieu Faverge
 * @author Cedric Castagnede
 * @author Florent Pruvost
 * @date 2020-03-03
 *
 ***
 *
 * @defgroup Workspace
 * @brief Group routines exposed to users about specific workspaces management
 *
 */

#include <stdlib.h>
#include "control/common.h"
#include "control/auxiliary.h"
#include "control/workspace.h"

/**
 *
 */
int chameleon_alloc_ibnb_tile(int M, int N, cham_tasktype_t func, int type, CHAM_desc_t **desc, int p, int q)
{
    int status;
    int IB, NB, MT, NT;
    int64_t lm, ln;
    CHAM_context_t *chamctxt;

    chamctxt = chameleon_context_self();
    if (chamctxt == NULL) {
        chameleon_fatal_error("chameleon_alloc_ibnb_tile", "CHAMELEON not initialized");
        return CHAMELEON_ERR_NOT_INITIALIZED;
    }

    /* Tune NB & IB depending on M & N; Set IBNBSIZE */
    status = chameleon_tune(func, M, N, 0);
    if (status != CHAMELEON_SUCCESS) {
        chameleon_error("chameleon_alloc_ibnb_tile", "chameleon_tune() failed");
        return CHAMELEON_ERR_UNEXPECTED;
    }

    /* Set MT & NT & allocate */
    NB = CHAMELEON_NB;
    IB = CHAMELEON_IB;
    MT = (M%NB==0) ? (M/NB) : (M/NB+1);
    NT = (N%NB==0) ? (N/NB) : (N/NB+1);

    /* Size is doubled for RH QR to store the reduction T */
    if ((chamctxt->householder == ChamTreeHouseholder) &&
        ((func == CHAMELEON_FUNC_SGELS)  ||
         (func == CHAMELEON_FUNC_DGELS)  ||
         (func == CHAMELEON_FUNC_CGELS)  ||
         (func == CHAMELEON_FUNC_ZGELS)  ||
         (func == CHAMELEON_FUNC_SGESVD) ||
         (func == CHAMELEON_FUNC_DGESVD) ||
         (func == CHAMELEON_FUNC_CGESVD) ||
         (func == CHAMELEON_FUNC_ZGESVD)))
    {
        NT *= 2;
    }

    lm = (int64_t)IB * (int64_t)MT;
    ln = (int64_t)NB * (int64_t)NT;

    return CHAMELEON_Desc_Create( desc, CHAMELEON_MAT_ALLOC_TILE, type, IB, NB, IB*NB,
                                  lm, ln, 0, 0, lm, ln, p, q );
}

/**
 *
 */
int chameleon_alloc_ipiv(int M, int N, cham_tasktype_t func, int type, CHAM_desc_t **desc, void **IPIV, int p, int q)
{
    int status;
    int NB, IB, MT, NT;
    int64_t lm, ln;
    size_t size;
    CHAM_context_t *chamctxt;

    chamctxt = chameleon_context_self();
    if (chamctxt == NULL) {
        chameleon_fatal_error("chameleon_alloc_ipiv", "CHAMELEON not initialized");
        return CHAMELEON_ERR_NOT_INITIALIZED;
    }

    /* Tune NB & IB depending on M & N; Set IBNBSIZE */
    status = chameleon_tune(func, M, N, 0);
    if (status != CHAMELEON_SUCCESS) {
        chameleon_error("chameleon_alloc_ipiv", "chameleon_tune() failed");
        return CHAMELEON_ERR_UNEXPECTED;
    }

    /* Set MT & NT & allocate */
    NB = CHAMELEON_NB;
    IB = CHAMELEON_IB;

    NT = (N%NB==0) ? (N/NB) : ((N/NB)+1);
    MT = (M%NB==0) ? (M/NB) : ((M/NB)+1);

    lm = (int64_t)IB * (int64_t)MT;
    ln = (int64_t)NB * (int64_t)NT;

    size = (size_t)(chameleon_min(MT, NT) * NB * NT * sizeof(int));
    if (size == 0) {
        *IPIV = NULL;
        return CHAMELEON_SUCCESS;
    }
    /* TODO: Fix the distribution for IPIV */
    *IPIV = (int*)malloc( size );

    return CHAMELEON_Desc_Create( desc, CHAMELEON_MAT_ALLOC_TILE, type, IB, NB, IB*NB,
                                  lm, ln, 0, 0, lm, ln, p, q );
}

/**
 *
 * @ingroup Workspace
 *
 *  CHAMELEON_Dealloc_Worksapce - Deallocate workspace descriptor allocated by
 *                            any workspace allocation routine.
 *
 *******************************************************************************
 *
 * @param[in] desc
 *          Workspace descriptor
 *
 *******************************************************************************
 *
 * @retval CHAMELEON_SUCCESS successful exit
 *
 */
int CHAMELEON_Dealloc_Workspace(CHAM_desc_t **desc)
{
    CHAM_context_t *chamctxt;

    chamctxt = chameleon_context_self();
    if (chamctxt == NULL) {
        chameleon_fatal_error("CHAMELEON_Dealloc_Workspace", "CHAMELEON not initialized");
        return CHAMELEON_ERR_NOT_INITIALIZED;
    }
    if (*desc == NULL) {
        chameleon_error("CHAMELEON_Dealloc_Workspace", "attempting to deallocate a NULL descriptor");
        return CHAMELEON_ERR_UNALLOCATED;
    }
    if ((*desc)->mat == NULL && (*desc)->use_mat == 1) {
        chameleon_error("CHAMELEON_Dealloc_Worspace", "attempting to deallocate a NULL pointer");
        return CHAMELEON_ERR_UNALLOCATED;
    }
    chameleon_desc_destroy( *desc );

    free(*desc);
    *desc = NULL;
    return CHAMELEON_SUCCESS;
}
