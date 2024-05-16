/**
 *
 * @file pzgepdf_qdwh.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2021 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon QDWH Polar Decomposition parallel algorithm
 *
 * @version 1.1.0
 * @author Mathieu Faverge
 * @author Hatem Ltaief
 * @date 2021-03-17
 * @precisions normal z -> s d c
 *
 */
#include "control/common.h"
#include "chameleon/flops.h"
#include <stdlib.h>
#include <stdio.h>

#if !defined(CHAMELEON_SIMULATION)
#include "coreblas/lapacke.h"
#endif

#if defined(CHAMELEON_SCHED_STARPU) || defined(CHAMELEON_USE_CUDA)
static int _zgepdf_qdwh_opt_genD = 1;
#else
static int _zgepdf_qdwh_opt_genD = 0;
#endif

static int _zgepdf_qdwh_opt_qr = 1;
static int _zgepdf_qdwh_opt_id = 1;
static int _zgepdf_qdwh_verbose = 0;

/**
 * @brief Compute the parameters of the QDWH iteration.
 *
 * @param[in] Li
 *        The estimated condition number on input.
 *
 * @param[out] params
 *        Array of size 3 that holds the parameters from the coming iteration
 *        computed from Li on exit.
 *
 * @retval The new estimated condition number after the iteration
 *
 */
static inline double
chameleon_pzgepdf_parameters( double Li, double *params )
{
    double L2, dd, sqd, a1;
    double a, b, c;

    L2  = Li * Li;
    /* Add this check because I add L2 that was slowly moving to 1.0 + eps, and generating NaN in the pow call */
    L2  = (L2 < 1.) ? L2 : 1.;
    dd  = pow( 4. * (1. - L2 ) / (L2 * L2), 1./3. );
    sqd = sqrt(1. + dd);
    a1  = sqd + sqrt( 8. - 4. * dd + 8. * (2. - L2) / (L2 * sqd) ) / 2.;
    a  = a1;
    b  = (a - 1.) * (a - 1.) / 4.;
    c  = a + b - 1.;

    params[0] = a;
    params[1] = b;
    params[2] = c;

    /* Return the new Li */
    return Li * (a + b * L2) / (1. + c * L2);
}

/**
 * @brief Initialize all the temporary descriptors and trees required.
 *
 * @param[in] U
 *        The input matrix U to determine the size and distribution of the top part.
 *
 * @param[in] H
 *        The input matrix H to determine the size and distribution of the bottom part.
 *
 * @param[out] qrtreeT, qrtreeB
 *        On exit, contains the reductions trees to perform QR factorization
 *        respectively on the Top and Bottom parts.
 *
 * @param[out] A
 *        A workspace to store a backup of A in order to compute the polar factor.
 *
 * @param[out] Ut
 *        A workspace to store a transposed U in the Cholesky iteration.
 *
 * @param[out] B1, TS1, TT1, Q1, D1
 *        Set of workspace that match the size of U for the QR iteration
 *
 * @param[out] B2, TS2, TT2, Q2, D2
 *        Set of workspace that match the size of H for the QR iteration
 *
 * @param[out] gemm_ws
 *        Pointer to the workspace data structure used by the GEMM operation
 *
 */
static inline void
chameleon_pzgepdf_qdwh_init( const CHAM_desc_t *U, const CHAM_desc_t *H,
                             libhqr_tree_t *qrtreeT, libhqr_tree_t *qrtreeB,
                             CHAM_desc_t *A,  CHAM_desc_t *Ut,
                             CHAM_desc_t *B1, CHAM_desc_t *TS1, CHAM_desc_t *TT1, CHAM_desc_t *Q1, CHAM_desc_t *D1,
                             CHAM_desc_t *B2, CHAM_desc_t *TS2, CHAM_desc_t *TT2, CHAM_desc_t *Q2, CHAM_desc_t *D2,
                             void **gemm_ws )
{
    CHAM_context_t *chamctxt;
    int ib, nb = U->nb;

    chamctxt = chameleon_context_self();
    ib = CHAMELEON_IB;

    /*
     * Let's initialize all the descriptors required for the computations:
     *     - the two Bs where U is copied for factorization
     *     - the two Qs where the Q is generated for QR iterations
     *     - the two Ds to create copies on the fly of the diagonal tiles in QR for faster computations on the GPU
     */
    chameleon_zdesc_copy_and_restrict( U, B1, U->m, U->n );
    chameleon_zdesc_copy_and_restrict( H, B2, U->n, U->n );

    chameleon_zdesc_copy_and_restrict( U, Q1, U->m, U->n );
    chameleon_zdesc_copy_and_restrict( H, Q2, U->n, U->n );

    if ( _zgepdf_qdwh_opt_genD ) {
        chameleon_zdesc_copy_and_restrict( U, D1, U->m, U->n );
        chameleon_zdesc_copy_and_restrict( H, D2, U->n, U->n );
    }

    /*
     * Let's add the workspace required to store the small T matrices from the QR factorization.
     * TS for the first QR factorization, TT for the level of reduction in the trees.
     * Note that while TS is lower trapezoidal, TT is a tile sparse matrix and
     * only tiles involved in the higher levels of the reduction tree are
     * generated.
     */
    chameleon_desc_init( TS1, CHAMELEON_MAT_ALLOC_TILE,
                         ChamComplexDouble, ib, nb, ib * nb,
                         ib * U->mt, nb * U->nt, 0, 0,
                         ib * U->mt, nb * U->nt, U->p, U->q,
                         NULL, NULL, NULL );
    chameleon_desc_init( TT1, CHAMELEON_MAT_ALLOC_TILE,
                         ChamComplexDouble, ib, nb, ib * nb,
                         ib * U->mt, nb * U->nt, 0, 0,
                         ib * U->mt, nb * U->nt, H->p, H->q,
                         NULL, NULL, NULL );

    chameleon_desc_init( TS2, CHAMELEON_MAT_ALLOC_TILE,
                         ChamComplexDouble, ib, nb, ib * nb,
                         ib * H->mt, nb * H->nt, 0, 0,
                         ib * H->mt, nb * H->nt, U->p, U->q,
                         NULL, NULL, NULL );
    chameleon_desc_init( TT2, CHAMELEON_MAT_ALLOC_TILE,
                         ChamComplexDouble, ib, nb, ib * nb,
                         ib * H->mt, nb * H->nt, 0, 0,
                         ib * H->mt, nb * H->nt, H->p, H->q,
                         NULL, NULL, NULL );

    /*
     * Create the adapted trees to perform the QR factorizations
     */
    {
        libhqr_matrix_t mat = {
            .mt    = B1->mt,
            .nt    = B1->nt,
            .nodes = B1->p * B1-> q,
            .p     = B1->p,
        };

        /* Tree for the top matrix */
        libhqr_init_hqr( qrtreeT, LIBHQR_QR, &mat,
                         -1,    /*low level tree   */
                         -1,    /* high level tree */
                         -1,    /* TS tree size    */
                         B1->p, /* High level size */
                         -1,    /* Domino */
                         0      /* TSRR (unstable) */ );

        /* Tree for the bottom matrix */
        mat.mt = B2->mt;
        mat.nt = B2->nt;
        libhqr_init_tphqr( qrtreeB, LIBHQR_TSQR,
                           mat.mt, _zgepdf_qdwh_opt_id ? (mat.nt-1) : 0, &mat,
                           /* high level tree (Could be greedy, but flat should reduce the volume of comm) */
                           LIBHQR_FLAT_TREE,
                           -1,   /* TS tree size    */
                           B2->p /* High level size */ );
    }

    /*
     * Let's create a backup of A, to be able to compute H in the end.
     */
    chameleon_zdesc_copy_and_restrict( U, A, U->m, U->n );

    /*
     * Let's create a transposed version of U for the solve step in the Cholesky iteration.
     */
    chameleon_desc_init( Ut, CHAMELEON_MAT_ALLOC_TILE,
                         ChamComplexDouble, U->mb, U->nb, U->mb * U->nb,
                         U->n, U->m, 0, 0,
                         U->n, U->m, U->p, U->q,
                         NULL, NULL, NULL );

    /*
     * Allocate the data descriptors for the lookahead if needed
     */
    *gemm_ws = CHAMELEON_zgemm_WS_Alloc( ChamNoTrans, ChamNoTrans, NULL, NULL, U );

    return;
}

/**
 * @brief Cleanup all the temporary data.
 *
 * @param[in,out] qrtreeT, qrtreeB
 *         On exit, the trees are destroyed.
 *
 * @param[in,out] A, Ut
 *        On exit, the workspace are destroyed.
 *
 * @param[in,out] B1, TS1, TT1, Q1, D1
 *        On exit, the workspace are destroyed.
 *
 * @param[in,out] B2, TS2, TT2, Q2, D2
 *
 */
static inline void
chameleon_pzgepdf_qdwh_fini( libhqr_tree_t *qrtreeT, libhqr_tree_t *qrtreeB,
                             CHAM_desc_t *A, CHAM_desc_t *Ut,
                             CHAM_desc_t *B1, CHAM_desc_t *TS1, CHAM_desc_t *TT1, CHAM_desc_t *Q1, CHAM_desc_t *D1,
                             CHAM_desc_t *B2, CHAM_desc_t *TS2, CHAM_desc_t *TT2, CHAM_desc_t *Q2, CHAM_desc_t *D2 )
{
    chameleon_desc_destroy( B1 );
    chameleon_desc_destroy( B2 );

    chameleon_desc_destroy( Q1 );
    chameleon_desc_destroy( Q2 );

    if ( _zgepdf_qdwh_opt_genD ) {
        chameleon_desc_destroy( D1 );
        chameleon_desc_destroy( D2 );
    }

    chameleon_desc_destroy( TS1 );
    chameleon_desc_destroy( TT1 );

    chameleon_desc_destroy( TS2 );
    chameleon_desc_destroy( TT2 );

    libhqr_finalize( qrtreeT );
    libhqr_finalize( qrtreeB );

    chameleon_desc_destroy( A );
    chameleon_desc_destroy( Ut );

    return;
}

/**
 * @brief Compute an estimation of the condition number with a QR factorization.
 *
 * @param[in] chamctxt
 *        The chameleon context.
 *
 * @param[in] Anorm
 *        The One norm of the scaled matrix A.
 *
 * @param[in] normest
 *        The estimated Two-norm of the matrix A.
 *
 * @param[in] qrtree
 *        The reduction tree to be used by the QR factorization of A.
 *
 * @param[in,out] A
 *        On entry, the M-by-N matrix A for which the condition number must be estimated.
 *        On exit, contains the QR factorization of A.
 *
 * @param[in,out] TS, TT, D
 *        The workspace required to performed a QR factorization on the matrix A.
 *
 * @param[in,out] W
 *        A workspace of size N-by-N to store the upper triangular R used in the estimation.
 *
 * @param[in,out] sequence
 *        The sequence to which the tasks should be submitted.
 *
 * @param[in,out] request
 *        The request to which the tasks should be submitted.
 *
 * @retval The estimator of the condition number of A.
 *
 */
static inline double
chameleon_pzgeqdwh_condest_qr( CHAM_context_t *chamctxt,
                               double Anorm, double normest, const libhqr_tree_t *qrtree,
                               CHAM_desc_t *A, CHAM_desc_t *TS, CHAM_desc_t *TT,
                               CHAM_desc_t *D, CHAM_desc_t *W,
                               RUNTIME_sequence_t *sequence, RUNTIME_request_t *request )
{
    CHAM_desc_t *upperA;
    double       Ainvnorm, Li;

    /*
     * Use QR for cond-est
     */
    chameleon_pzgeqrf_param( _zgepdf_qdwh_opt_genD, A->nt, qrtree, A, TS, TT, D,
                             sequence, request );

    /* Flush unused data in coming algorithms */
    CHAMELEON_Desc_Flush( TS, sequence );
    CHAMELEON_Desc_Flush( TT, sequence );
    if ( D != NULL ) {
        CHAMELEON_Desc_Flush( D, sequence );
    }

    /* Copy the R part */
    upperA = chameleon_desc_submatrix( A, 0, 0, A->n, A->n );
    chameleon_pzlacpy( ChamUpper, upperA, W, sequence, request );
    chameleon_pztrtri( ChamUpper, ChamNonUnit, W, sequence, request );
    chameleon_pzlange_generic( ChamOneNorm, ChamUpper, ChamNonUnit,
                               W, &Ainvnorm, sequence, request );

    CHAMELEON_Desc_Flush( A, sequence );
    CHAMELEON_Desc_Flush( W, sequence );
    chameleon_sequence_wait( chamctxt, sequence );
    free( upperA );

    Li = ((double)1.0 / Ainvnorm) / Anorm;
    Li = normest / 1.1 * Li;

    return Li;
}

/**
 * @brief Perform a QR iteration step to update the polar normal factor
 *
 * This function submit the tasks to compute:
 *      U_{n+1} = alpha * Q1 * Q2' + beta * U_{n}, with:
 *
 * (Q1, Q2) defined such that:
 *                   [ Q1 ] R = [ B1 ] = [ sqrt(c) * U ]
 *                   [ Q2 ]     [ B2 ]   [ Id          ]
 *
 * with Q1, B1 of size m-by-n, and Q2, B2 of size n-by-n.
 *
 * The algorithm actually performs the following steps:
 *    B1 <= qr(B1)
 *    B2 <= qr(B2)
 *    TPQRT( B1, B2 )
 * and then apply the associated operations to generate Q1 and Q2, by first
 * applying TPGQRT( Q1,Q2), and then generating the two separated Qs.
 *
 * @param[in] do_qr
 *        Specify if the matrix Q1 has already been factorized or not.
 *
 * @param[in] last
 *        Specify if this is the last expected iteration or not, to eventually
 *        backup U_{n} in B1.
 *        If last, U is copied into B1 to be able to compute ||U_{n+1} - U_{n}||.
 *        If !last, no need to backup as more iterations are required.
 *
 * @param[in] params
 *        The parameters computed by chameleon_pzgepdf_parameters().
 *
 * @param[in] qrtreeT
 *        Defines the QR tree used to perform the factorization of B1. Must be
 *        identical to the one already used to factorize it if !do_qr.
 *
 * @param[in] qrtreeB
 *        Defines the QR tree used to perform the factorization of B2, and its
 *        combination with B1.
 *
 * @param[in,out] U
 *        The polar normal factor matrix U of size m-by-n to update.
 *        On entry, contains U_{n}. On exit, contains U_{n+1}.
 *
 * @param[in,out] B1
 *        B1 is a workspace of size m-by-n to perform the QR factorization.
 *        On entry, if !do_qr, B1 is expected to contain qr(U).
 *        otherwise if do_qr, B1 is expected to contain U.
 *        On exit, if last, B1 contains U_{n}, otherwise it is destroyed.
 *
 * @param[in,out] TS1
 *        TS1 is the QR workspace used to store the small T during the QR factorization.
 *        On entry, if !do_qr, TS1 is the output of the QR factorization of B1
 *        otherwise if do_qr, the content of TS1 is just destroyed.
 *        On exit, TS1 is destroyed.
 *
 * @param[in,out] TT1
 *        TT1 is the QR workspace used to store the small T during the QR factorization.
 *        On entry, if !do_qr, TT1 is the output of the QR factorization of B1
 *        otherwise if do_qr, the content of TT1 is just destroyed.
 *        On exit, TT1 is destroyed.
 *
 * @param[in,out] D1
 *        D1 is the QR workspace used to store the temporary diagonal tiles to
 *        speed up the computation and increase the level of parallelism.
 *        On entry, if !do_qr, D1 is the output of the QR factorization of B1
 *        otherwise if do_qr, the content of D1 is just destroyed.
 *        On exit, D1 is destroyed.
 *
 * @param[in,out] Q1
 *        Q1 is the a workspace of size m-by-n to generate the Q associated to B1.
 *        On exit, Q1 is destroyed.
 *
 * @param[in,out] B2
 *        B2 is a workspace of size n-by-n to perform the QR factorization.
 *        On exit, B2 is destroyed.
 *
 * @param[in,out] TS2
 *        TS2 is the QR workspace used to store the small T during the QR factorization.
 *        On exit, TS2 is destroyed.
 *
 * @param[in,out] TT2
 *        TT2 is the QR workspace used to store the small T during the QR factorization.
 *        On exit, TT2 is destroyed.
 *
 * @param[in,out] D2
 *        D2 is the QR workspace used to store the temporary diagonal tiles to
 *        speed up the computation and increase the level of parallelism.
 *        On exit, D2 is destroyed.
 *
 * @param[in,out] Q2
 *        Q2 is the a workspace of size n-by-n to generate the Q associated to B2.
 *        On exit, Q2 is destroyed.
 *
 */
static inline void
chameleon_pzgepdf_qdwh_qrstep( int do_qr, int last, double *params,
                               const libhqr_tree_t *qrtreeT, const libhqr_tree_t *qrtreeB,
                               CHAM_desc_t *U, void *gemm_ws,
                               CHAM_desc_t *B1, CHAM_desc_t *TS1, CHAM_desc_t *TT1, CHAM_desc_t *Q1, CHAM_desc_t *D1,
                               CHAM_desc_t *B2, CHAM_desc_t *TS2, CHAM_desc_t *TT2, CHAM_desc_t *Q2, CHAM_desc_t *D2,
                               RUNTIME_sequence_t *sequence, RUNTIME_request_t *request )
{
    double a = params[0];
    double b = params[1];
    double c = params[2];
    double sqrt_c = sqrt( c );
    double alpha, beta;

    /*
     * Generate the matrix B = [ B1 ] = [ sqrt(c) * U ]
     *                         [ B2 ] = [ Id          ]
     */
    if( do_qr ) {
        chameleon_pzlacpy( ChamUpperLower, U, B1,
                           sequence, request );
        chameleon_pzlascal( ChamUpperLower, sqrt_c, B1, sequence, request );
    }
    else {
        /* Scale only the R factor as B1 is already factorized */
        chameleon_pzlascal( ChamUpper, sqrt_c, B1, sequence, request );
    }

    chameleon_pzlaset( ChamUpperLower, 0.0, 1.0, B2, sequence, request );

    /*
     * Factorize B, and generate the associated Q.
     */
    chameleon_pzgepdf_qr( _zgepdf_qdwh_opt_genD, do_qr, _zgepdf_qdwh_opt_id, qrtreeT, qrtreeB,
                          B1, TS1, TT1, D1, Q1,
                          B2, TS2, TT2, D2, Q2,
                          sequence, request );

    /*
     * Backup U in B1 on the last iteration
     */
    if ( last ) {
        chameleon_pzlacpy( ChamUpperLower, U, B1,
                           sequence, request );
    }

    /*
     * Gemm to update U
     *  U = ( (a-b/c)/sqrt(c) ) * Q1 * Q2' + (b/c) * U
     */
    beta  = b / c;
    alpha = ( a - beta ) / sqrt_c;
    chameleon_pzgemm( gemm_ws, ChamNoTrans, ChamConjTrans,
                      alpha, Q1, Q2, beta, U,
                      sequence, request );

    CHAMELEON_Desc_Flush( Q1, sequence );
    CHAMELEON_Desc_Flush( Q2, sequence );
}

/**
 * @brief Perform a Cholesky iteration step to update the polar normal factor
 *
 * This function submit the tasks to compute U_{n+1} = alpha X + beta U_{n}, with:
 *
 * X the solution of the system:
 *     (c * U' * U + I) X = U'
 *
 * @param[in] last
 *        Specify if this is the last expected iteration or not, to eventually backup U_{n}.
 *        If last, U is copied into Un to be able to compute ||U_{n+1} - U_{n}||.
 *        If !last, no need to backup as more iterations are required.
 *
 * @param[in] params
 *        The parameters computed by chameleon_pzgepdf_parameters().
 *
 * @param[in,out] U
 *        The polar normal factor matrix U of size m-by-n to update.
 *        On entry, contains U_{n}. On exit, contains U_{n+1}.
 *
 * @param[out] Un
 *        Un must be of the same size as U (m-by-n)
 *        On exit, Un contains a backup of U_{n} for convolution norm estimation if last.
 *        If !last, Un is not referenced.
 *
 * @param[in,out] A
 *        A is workspace of the same size as H (n-by-n)
 *        It is used internally to store (c * U' * U + I), and is destroyed on exit.
 *
 * @param[in,out] X
 *        X is workspace of the same size as U' (n-by-m).
 *        It is used to store the right hand side of the system to solve, and it is destroyed on exit.
 *
 */
static inline void
chameleon_pzgepdf_qdwh_postep( int last, double *params, CHAM_desc_t *U,
                               CHAM_desc_t *Un, CHAM_desc_t *A, CHAM_desc_t *X,
                               RUNTIME_sequence_t *sequence, RUNTIME_request_t *request )
{
    double a = params[0];
    double b = params[1];
    double c = params[2];
    double alpha, beta;

    /*
     * Compute A = c * U' * U + I
     */
    chameleon_pzlaset( ChamUpperLower, 0., 1., A, sequence, request );
    chameleon_pzherk( ChamUpper, ChamConjTrans,
                      c, U, 1.0, A, sequence, request );

    /*
     * Solve A X = B, with B = U'
     */
    chameleon_pzpotrf( ChamUpper, A, sequence, request );

    chameleon_pztradd( ChamUpperLower, ChamConjTrans, 1., U, 0., X, sequence, request );
    chameleon_pztrsm( ChamLeft, ChamUpper, ChamConjTrans, ChamNonUnit, 1.0, A, X, sequence, request );
    chameleon_pztrsm( ChamLeft, ChamUpper, ChamNoTrans,   ChamNonUnit, 1.0, A, X, sequence, request );

    CHAMELEON_Desc_Flush( A, sequence );

    /*
     * Backup U in Un on the last iteration
     */
    if ( last ) {
        chameleon_pzlacpy( ChamUpperLower, U, Un,
                           sequence, request );
    }

    /*
     * Update U = (a-b/c) * X' + (b/c) * U
     */
    beta  = b / c;
    alpha = a - beta;
    chameleon_pztradd( ChamUpperLower, ChamConjTrans,
                       alpha, X, beta, U,
                       sequence, request );

    CHAMELEON_Desc_Flush( X, sequence );
}

/**
 * @brief Symmetrize/Hermitian a matrix by computing .5 * (A + A')
 *
 * @param[in,out] A
 *        On entry, the N-by-N matrix to symmetrize.
 *        On exit, the symmetrized matrix equal to .5 * (A'+A)
 *
 * @param[in,out] W
 *        A workspace of size N-by-N.
 *
 * @param[in,out] sequence
 *        The sequence to which the tasks should be submitted.
 *
 * @param[in,out] request
 *        The request to which the tasks should be submitted.
 *
 */
static inline void
chameleon_pzgepdf_qdwh_sym( CHAM_desc_t *A, CHAM_desc_t *W,
                            RUNTIME_sequence_t *sequence, RUNTIME_request_t *request )
{
    /*
     * Symmetrize/Hermitian A using W as a workspace.
     *      A = (A+A')/2 { W = A; A = 0.5 * ( A + W' ); }
     *
     * Works only for square matrices.
     */
    chameleon_pzlacpy( ChamUpperLower, A, W,
                       sequence, request );
    chameleon_pztradd( ChamUpperLower, ChamConjTrans,
                       .5, W, .5, A, sequence, request );
    CHAMELEON_Desc_Flush( W, sequence );
}

void
chameleon_pzgepdf_qdwh( cham_mtxtype_t mtxtype, CHAM_desc_t *descU, CHAM_desc_t *descH,
                        gepdf_info_t *info,
                        RUNTIME_sequence_t *sequence, RUNTIME_request_t *request )
{
    CHAM_context_t *chamctxt;

    /*
     * Add additional sequence/request to avoid issue with the large number of
     * task types in Parsec DTD
     */
    RUNTIME_sequence_t *sequence_po = NULL;
    RUNTIME_sequence_t *sequence_qr = NULL;
    RUNTIME_request_t   request_po = RUNTIME_REQUEST_INITIALIZER;
    RUNTIME_request_t   request_qr = RUNTIME_REQUEST_INITIALIZER;

    RUNTIME_sequence_t *sequence_it = sequence;
    RUNTIME_request_t  *request_it  = request;

    CHAM_desc_t descA, descUt;
    CHAM_desc_t descB1, descTS1, descTT1, descD1, descQ1, *D1ptr;
    CHAM_desc_t descB2, descTS2, descTT2, descD2, descQ2, *D2ptr;
    libhqr_tree_t qrtreeT, qrtreeB;
    void *gemm_ws;

    double conv = 100.;
    double Li, params[3];
    double normest, Unorm;
    int it, itconv, facto = -1;

#if !defined(CHAMELEON_SIMULATION)
    double eps  = LAPACKE_dlamch_work('e');
#else
#if defined(PRECISION_z) || defined(PRECISION_d)
    double eps  = 1.e-15;
#else
    double eps  = 1.e-7;
#endif
#endif
    double tol1 = 5. * eps;
    double tol3 = pow( tol1, 1./3. );
    double id_flops_ratio = ( _zgepdf_qdwh_opt_id == 1 ) ? .5 : 1.5;

    chamctxt = chameleon_context_self();
    if (sequence->status != CHAMELEON_SUCCESS) {
        return;
    }
    assert( chamctxt->scheduler != RUNTIME_SCHED_PARSEC );

    if ( info ) {
        info->itQR = 0;
        info->itPO = 0;
        info->flops = 0.;
    }

    /* Generate temporary descriptors */
    chameleon_pzgepdf_qdwh_init( descU, descH,
                                 &qrtreeT, &qrtreeB,
                                 &descA, &descUt,
                                 &descB1, &descTS1, &descTT1, &descQ1, &descD1,
                                 &descB2, &descTS2, &descTT2, &descQ2, &descD2, &gemm_ws );
    if ( _zgepdf_qdwh_opt_genD ) {
        D1ptr = &descD1;
        D2ptr = &descD2;
    }
    else {
        D1ptr = NULL;
        D2ptr = NULL;
    }

    if ( chamctxt->scheduler == RUNTIME_SCHED_PARSEC ) {
        chameleon_sequence_create( chamctxt, &sequence_qr );
        chameleon_sequence_create( chamctxt, &sequence_po );
    }

    /*
     * Backup the original A, to compute H at the end of the algorithm
     */
    chameleon_pzlacpy( ChamUpperLower, descU, &descA,
                       sequence, request );

    /*
     * Two norm estimation
     * Warning: a sequence_wait here is included at the end of the genm2 call
     */
    chameleon_pzgenm2( 1.e-1, descU, &normest, sequence, request );

    if ( _zgepdf_qdwh_verbose && (CHAMELEON_Comm_rank() == 0) ) {
        fprintf( stderr, "pzgepdf_qdwh: ||A||_2 = %e\n", normest );
    }

    /*
     * Scale the original U to form the U0 of the iterative loop, and copy U0 into B1.
     */
    chameleon_pzlascal( ChamUpperLower, 1.0 / normest, descU,
                        sequence, request );
    chameleon_pzlacpy( ChamUpperLower, descU, &descB1,
                       sequence, request );

    switch( mtxtype ) {
#if defined(PRECISION_z) || defined(PRECISION_c)
    case ChamHermitian:
        chameleon_pzlansy_generic( ChamOneNorm, ChamUpper, ChamConjTrans, descU, &Unorm, sequence, request );
        break;
#endif
    case ChamSymmetric:
        chameleon_pzlansy_generic( ChamOneNorm, ChamUpper, ChamTrans, descU, &Unorm, sequence, request );
        break;
    default:
        chameleon_pzlange_generic( ChamOneNorm, ChamUpperLower, ChamNonUnit, descU, &Unorm, sequence, request );
    }
    /* Wait for Unorm */
    chameleon_sequence_wait( chamctxt, sequence );

    /*
     * Estimate of the condition number
     */
    Li = chameleon_pzgeqdwh_condest_qr( chamctxt, Unorm, normest, &qrtreeT,
                                        &descB1, &descTS1, &descTT1, D1ptr,
                                        &descB2, sequence, request );
    if ( info ) {
        info->flops += flops_zgeqrf( descB1.m, descB1.n );
        info->flops += flops_ztrtri( descB1.n );
    }

    /*
     * Estimate the required number of iteration
     */
    {
        double Liconv = Li;
        itconv = 0;
        while( (itconv == 0) || (fabs(1-Liconv) > tol1) ) {
            /* To find the minimum number of iterations to converge.
             * itconv = number of iterations needed until |Li - 1| < tol1
             * This should converge in less than 50 iterations
             */
            if ( itconv > 100 ) {
                chameleon_fatal_error("chameleon_pzgepdf_qdwh", "Failed to converge in less than 100 iterations");
                break;
            }
            itconv++;
            Liconv = chameleon_pzgepdf_parameters( Liconv, params );
        }
    }

    it = 0;
    while( (conv > tol3) || (it < itconv) ) {
        int last;

        Li = chameleon_pzgepdf_parameters( Li, params );
        it++;
        last = ( it >= itconv );

        if ( params[2] > 100 ) {
            int do_qr = (!_zgepdf_qdwh_opt_qr) || (it > 1);

            if ( (chamctxt->scheduler == RUNTIME_SCHED_PARSEC) &&
                 ( sequence_it != sequence_qr ) )
            {
                chameleon_sequence_wait( chamctxt, sequence_it );
                sequence_it = sequence_qr;
                request_it = &request_qr;
            }

            chameleon_pzgepdf_qdwh_qrstep( do_qr, last, params,
                                           &qrtreeT, &qrtreeB, descU, gemm_ws,
                                           &descB1, &descTS1, &descTT1, &descQ1, D1ptr,
                                           &descB2, &descTS2, &descTT2, &descQ2, D2ptr,
                                           sequence_it, request_it );

            facto = 0;
            if ( info ) {
                info->itQR++;
                info->flops += do_qr ? flops_zgeqrf( descB1.m, descB1.n ) : 0;
                info->flops += flops_zungqr( descB1.m, descB1.n, descB1.n );

                info->flops += id_flops_ratio * ( flops_zgeqrf( descB2.m, descB2.n ) +
                                                  flops_zungqr( descB2.m, descB2.n, descB2.n ) );

                info->flops += flops_zgemm( descB1.m, descB1.n, descB1.n );
            }
        }
        else
        {
            if ( (chamctxt->scheduler == RUNTIME_SCHED_PARSEC) &&
                 ( sequence_it != sequence_po ) )
            {
                chameleon_sequence_wait( chamctxt, sequence_it );
                sequence_it = sequence_po;
                request_it = &request_po;
            }

            chameleon_pzgepdf_qdwh_postep( last, params, descU,
                                           &descB1, &descQ2, &descUt,
                                           sequence_it, request_it );

            facto = 1;
            if ( info ) {
                info->itPO++;
                info->flops += flops_zherk( descB1.m, descB1.n );
                info->flops += flops_zpotrf( descB1.n );
                info->flops += 2. * flops_ztrsm( ChamLeft, descUt.m, descUt.n );
            }
        }

        if ( mtxtype != ChamGeneral ) {
            /* In the hermitian case U is square as B2, so B2 can be used as a workspace */
            chameleon_pzgepdf_qdwh_sym( descU, &descB2, sequence_it, request_it );
        }

        /*
         * To find the conv-norm, compute || U_{n+1} - U_{n} ||
         */
        conv = 10.;
        if( last ){
            chameleon_pztradd( ChamUpperLower, ChamNoTrans, 1.0, descU, -1.0, &descB1,
                               sequence, request  );

            chameleon_pzlange_generic( ChamFrobeniusNorm, ChamUpperLower, ChamNonUnit,
                                       &descB1, &conv, sequence_it, request_it );
            CHAMELEON_Desc_Flush( &descB1, sequence );
            chameleon_sequence_wait( chamctxt, sequence_it );
        }

        if ( _zgepdf_qdwh_verbose && (CHAMELEON_Comm_rank() == 0) ) {
            fprintf( stderr, "pzgepdf_qdwh: %02d %-5s %e\n",
                     it, facto == 0 ? "QR" : "PO", conv );
        }
    }

    if ( (chamctxt->scheduler == RUNTIME_SCHED_PARSEC) &&
         ( sequence_it != sequence ) )
    {
        chameleon_sequence_wait( chamctxt, sequence_it );
        chameleon_sequence_destroy( chamctxt, sequence_qr );
        chameleon_sequence_destroy( chamctxt, sequence_po );
    }

    /*
     * Compute H through the QDWH subroutine that
     * H = U'*A; H = .5(H+H')
     */
    switch( mtxtype ) {
#if defined(PRECISION_z) || defined(PRECISION_c)
    case ChamHermitian:
        chameleon_pzhemm( ChamRight, ChamUpper,
                          1., descU, &descA,
                          0., descH, sequence, request );
        if ( info ) {
            info->flops += flops_zhemm( ChamRight, descU->m, descU->n );
        }
        break;
#endif
    case ChamSymmetric:
        chameleon_pzsymm( ChamRight, ChamUpper,
                          1., descU, &descA,
                          0., descH, sequence, request );
        if ( info ) {
            info->flops += flops_zsymm( ChamRight, descU->m, descU->n );
        }
        break;
    default:
        chameleon_pzgemm( gemm_ws, ChamConjTrans, ChamNoTrans,
                          1., descU, &descA,
                          0., descH,
                          sequence, request );
        if ( info ) {
            info->flops += flops_zgemm( descH->m, descH->n, descA.m );
        }
    }
    CHAMELEON_Desc_Flush( &descA, sequence );

    /* Make sure H is hermitian */
    chameleon_pzgepdf_qdwh_sym( descH, &descB2, sequence, request );

    CHAMELEON_Desc_Flush( &descB2, sequence );

    chameleon_sequence_wait( chamctxt, sequence );
    chameleon_pzgepdf_qdwh_fini( &qrtreeT, &qrtreeB, &descA, &descUt,
                                 &descB1, &descTS1, &descTT1, &descQ1, &descD1,
                                 &descB2, &descTS2, &descTT2, &descQ2, &descD2 );

    return;
}
