/**
 *
 * @file zgepdf_qr.c
 *
 * @copyright 2012-2021 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 * @copyright 2016-2020 KAUST. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zgepdf_qr wrapper.
 *
 * This wrapper exists only for testing purpose to ease the creation of
 * timer/testing fot this subroutine of the QDWH/Zolo algorithms. That is why
 * only the Tile version is available.
 *
 * @version 1.1.0
 * @author Mathieu Faverge
 * @date 2020-12-01
 * @precisions normal z -> s d c
 *
 */
#include "control/common.h"

/**
 ******************************************************************************
 *
 * @ingroup CHAMELEON_Complex64_t
 *
 * @brief Perfom a QR factorization of a matrix A and generate the associated Q.
 *
 * A is of the form [ A1 ] and Q is generated in a similar form [ Q1 ].
 *                  [ A2 ]                                      [ Q2 ]
 *
 * @warning This function is NOT a USER routine and should not be called directly.
 *
 *******************************************************************************
 *
 * @param[in] doqr
 *            Specify if A1 is already factorized or if it needs to be done within this call.
 *
 * @param[in] optid
 *            Specify if A2 is an identity structure matrix and optimization is performed or not.
 *
 * @param[in] qrtreeT
 *            Describe the reduction tree for the factorization of A1.
 *
 * @param[in] qrtreeB
 *            Describe the reduction tree for the factorization of A2.
 *
 * @param[in] A1
 *
 * @param[in] TS1
 *
 * @param[in] TT1
 *
 * @param[in] Q1
 *
 * @param[in] A2
 *
 * @param[in] TS2
 *
 * @param[in] TT2
 *
 * @param[in] Q2
 *
 *******************************************************************************
 *
 * @return
 *          \retval CHAMELEON_SUCCESS successful exit
 *          \retval <0 if -i, the i-th argument had an illegal value
 *
 *******************************************************************************
 *
 * @sa CHAMELEON_cgepdf_qr_Tile
 * @sa CHAMELEON_dgepdf_qr_Tile
 * @sa CHAMELEON_sgepdf_qr_Tile
 *
 */
int CHAMELEON_zgepdf_qr_Tile( int doqr, int optid,
                              const libhqr_tree_t *qrtreeT, const libhqr_tree_t *qrtreeB,
                              CHAM_desc_t *A1, CHAM_desc_t *TS1, CHAM_desc_t *TT1, CHAM_desc_t *Q1,
                              CHAM_desc_t *A2, CHAM_desc_t *TS2, CHAM_desc_t *TT2, CHAM_desc_t *Q2 )
{
    CHAM_context_t *chamctxt;
    RUNTIME_sequence_t *sequence = NULL;
    RUNTIME_request_t request = RUNTIME_REQUEST_INITIALIZER;
    int status;

    chamctxt = chameleon_context_self();
    if (chamctxt == NULL) {
        chameleon_fatal_error("CHAMELEON_zgepdf_qr_Tile", "CHAMELEON not initialized");
        return CHAMELEON_ERR_NOT_INITIALIZED;
    }
    chameleon_sequence_create( chamctxt, &sequence );

    chameleon_pzgepdf_qr( 1, doqr, optid, qrtreeT, qrtreeB,
                          A1, TS1, TT1, NULL, Q1,
                          A2, TS2, TT2, NULL, Q2,
                          sequence, &request );

    CHAMELEON_Desc_Flush( Q1, sequence );
    CHAMELEON_Desc_Flush( Q2, sequence );

    chameleon_sequence_wait( chamctxt, sequence );
    status = sequence->status;
    chameleon_sequence_destroy( chamctxt, sequence );
    return status;
}
