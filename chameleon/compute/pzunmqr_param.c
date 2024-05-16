/**
 *
 * @file pzunmqr_param.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2021 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zunmqr_param parallel algorithm
 *
 * @version 1.0.0
 * @author Mathieu Faverge
 * @author Raphael Boucherie
 * @date 2020-03-03
 * @precisions normal z -> s d c
 *
 */
#include "control/common.h"
#include <stdlib.h>

#define A(m,n) A,  m,  n
#define C(m,n) C,  m,  n
#define T(m,n) T,  m,  n
#define D(m,n) D,  m,  n

/**
 *  Parallel application of Q using tile V - QR factorization - dynamic scheduling
 */
void chameleon_pzunmqr_param( int genD, const libhqr_tree_t *qrtree,
                              cham_side_t side, cham_trans_t trans,
                              CHAM_desc_t *A, CHAM_desc_t *C,
                              CHAM_desc_t *TS, CHAM_desc_t *TT, CHAM_desc_t *D,
                              RUNTIME_sequence_t *sequence, RUNTIME_request_t *request )
{
    CHAM_context_t *chamctxt;
    RUNTIME_option_t options;
    CHAM_desc_t *T;
    size_t ws_worker = 0;
    size_t ws_host = 0;

    int k, m, n, i, p;
    int temppm, temppn, tempmm, tempnn, tempkn,tempkmin;
    int ib, KT, L;
    int node, nbtiles, *tiles;

    chamctxt = chameleon_context_self();
    if (sequence->status != CHAMELEON_SUCCESS) {
        return;
    }
    RUNTIME_options_init(&options, chamctxt, sequence, request);

    ib = CHAMELEON_IB;

    if ( D == NULL ) {
        D    = A;
        genD = 0;
    }

    /*
     * zunmqr  = A->nb * ib
     * ztpmqrt = A->nb * ib
     */
    ws_worker = A->nb * ib;

#if defined(CHAMELEON_USE_CUDA)
    /* Worker space
     *
     * zunmqr  =     A->nb * ib
     * ztpmqrt = 3 * A->nb * ib
     */
    ws_worker = chameleon_max( ws_worker, ib * A->nb * 3 );
#endif

    ws_worker *= sizeof(CHAMELEON_Complex64_t);
    ws_host   *= sizeof(CHAMELEON_Complex64_t);

    RUNTIME_options_ws_alloc( &options, ws_worker, ws_host );

    /* Initialisation of tiles */
    tiles = (int*)calloc( qrtree->mt, sizeof(int) );

    KT = chameleon_min( A->mt, A->nt );
    if (side == ChamLeft ) {
        if (trans == ChamConjTrans) {
            /*
             *  ChamLeft / ChamConjTrans
             */
            for (k = 0; k < KT; k++) {
                RUNTIME_iteration_push(chamctxt, k);

                tempkn = k == A->nt - 1 ? A->n - k * A->nb : A->nb;

                T = TS;
                for (i = 0; i < qrtree->getnbgeqrf(qrtree, k); i++) {
                    p = qrtree->getm(qrtree, k, i);

                    temppm   = p == C->mt-1 ? C->m - p * C->mb : C->mb;
                    tempkmin = chameleon_min( temppm, tempkn );


                    if ( genD ) {
                        int tempDpm = p == D->mt-1 ? D->m-p*D->mb : D->mb;

                        INSERT_TASK_zlacpy(
                            &options,
                            ChamLower, tempDpm, tempkmin, A->nb,
                            A(p, k),
                            D(p, k) );
#if defined(CHAMELEON_USE_CUDA)
                        INSERT_TASK_zlaset(
                            &options,
                            ChamUpper, tempDpm, tempkmin,
                            0., 1.,
                            D(p, k) );
#endif
                    }
                    for (n = 0; n < C->nt; n++) {
                        tempnn = n == C->nt-1 ? C->n-n*C->nb : C->nb;
                        INSERT_TASK_zunmqr(
                            &options, side, trans,
                            temppm, tempnn, tempkmin, ib, T->nb,
                            D(p, k),
                            T(p, k),
                            C(p, n));
                    }
                    RUNTIME_data_flush( sequence, D(p, k) );
                    RUNTIME_data_flush( sequence, T(p, k) );
                }

                /* Setting the order of the tiles*/
                nbtiles = libhqr_walk_stepk( qrtree, k, tiles );

                for (i = 0; i < nbtiles; i++) {
                    m = tiles[i];
                    p = qrtree->currpiv(qrtree, k, m);

                    tempmm = m == C->mt-1 ? C->m-m*C->mb : C->mb;

                    if( qrtree->gettype(qrtree, k, m) == LIBHQR_KILLED_BY_TS ) {
                        /* TS kernel */
                        L = 0;
                        T = TS;
                    }
                    else {
                        /* TT kernel */
                        L = tempmm;
                        T = TT;
                    }
                    for (n = 0; n < C->nt; n++) {
                        tempnn = n == C->nt-1 ? C->n-n*C->nb : C->nb;

                        node = C->get_rankof( C, m, n );
                        RUNTIME_data_migrate( sequence, C(p, n), node );
                        RUNTIME_data_migrate( sequence, C(m, n), node );

                        INSERT_TASK_ztpmqrt(
                            &options, side, trans,
                            tempmm, tempnn, tempkn, L, ib, T->nb,
                            A(m, k),
                            T(m, k),
                            C(p, n),
                            C(m, n));
                    }
                    RUNTIME_data_flush( sequence, A(m, k) );
                    RUNTIME_data_flush( sequence, T(m, k) );
                }

                /* Restore the original location of the tiles */
                for (n = 0; n < C->nt; n++) {
                    RUNTIME_data_migrate( sequence, C(k, n),
                                          C->get_rankof( C, k, n ) );
                }

                RUNTIME_iteration_pop(chamctxt);
            }
        }
        /*
         *  ChamLeft / ChamNoTrans
         */
        else {
            for (k = KT-1; k >= 0; k--) {
                RUNTIME_iteration_push(chamctxt, k);

                tempkn = k == A->nt-1 ? A->n-k*A->nb : A->nb;

                /* Setting the order of the tiles*/
                nbtiles = libhqr_walk_stepk( qrtree, k, tiles );

                for (i = nbtiles-1; i >= 0; i--) {
                    m = tiles[i];
                    p = qrtree->currpiv(qrtree, k, m);

                    tempmm = m == C->mt-1 ? C->m-m*C->mb : C->mb;

                    if( qrtree->gettype(qrtree, k, m) == LIBHQR_KILLED_BY_TS ) {
                        /* TS kernel */
                        L = 0;
                        T = TS;
                    }
                    else {
                        /* TT kernel */
                        L = tempmm;
                        T = TT;
                    }
                    for (n = 0; n < C->nt; n++) {
                        tempnn = n == C->nt-1 ? C->n-n*C->nb : C->nb;

                        node = C->get_rankof( C, m, n );
                        RUNTIME_data_migrate( sequence, C(p, n), node );
                        RUNTIME_data_migrate( sequence, C(m, n), node );

                        INSERT_TASK_ztpmqrt(
                            &options, side, trans,
                            tempmm, tempnn, tempkn, L, ib, T->nb,
                            A(m, k),
                            T(m, k),
                            C(p, n),
                            C(m, n));
                    }
                    RUNTIME_data_flush( sequence, A(m, k) );
                    RUNTIME_data_flush( sequence, T(m, k) );
                }

                T = TS;
                for (i = 0; i < qrtree->getnbgeqrf(qrtree, k); i++) {
                    p = qrtree->getm(qrtree, k, i);

                    temppm   = p == C->mt-1 ? C->m-p*C->mb : C->mb;
                    tempkmin = chameleon_min( temppm, tempkn );


                    if ( genD ) {
                        int tempDpm = p == D->mt-1 ? D->m-p*D->mb : D->mb;

                        INSERT_TASK_zlacpy(
                            &options,
                            ChamLower, tempDpm, tempkmin, A->nb,
                            A(p, k),
                            D(p, k) );
#if defined(CHAMELEON_USE_CUDA)
                        INSERT_TASK_zlaset(
                            &options,
                            ChamUpper, tempDpm, tempkmin,
                            0., 1.,
                            D(p, k) );
#endif
                    }

                    for (n = 0; n < C->nt; n++) {
                        tempnn = n == C->nt-1 ? C->n-n*C->nb : C->nb;

                        RUNTIME_data_migrate( sequence, C(p, n),
                                              C->get_rankof( C, p, n ) );

                        INSERT_TASK_zunmqr(
                            &options, side, trans,
                            temppm, tempnn, tempkmin, ib, T->nb,
                            D(p, k),
                            T(p, k),
                            C(p, n));
                    }
                    RUNTIME_data_flush( sequence, D(p, k) );
                    RUNTIME_data_flush( sequence, T(p, k) );
                }
                RUNTIME_iteration_pop(chamctxt);
            }
        }
    }
    /*
     *  ChamRight / ChamConjTrans
     */
    else {
        if (trans == ChamConjTrans) {
            for (k = KT-1; k >= 0; k--) {
                RUNTIME_iteration_push(chamctxt, k);

                tempkn = k == A->nt-1 ? A->n - k * A->nb : A->nb;

                /* Setting the order of the tiles*/
                nbtiles = libhqr_walk_stepk( qrtree, k, tiles );

                for (i = nbtiles-1; i >= 0; i--) {
                    n = tiles[i];
                    p = qrtree->currpiv(qrtree, k, n);

                    tempnn = n == C->nt-1 ? C->n-n*C->nb : C->nb;

                    if( qrtree->gettype(qrtree, k, n) == LIBHQR_KILLED_BY_TS ) {
                        /* TS kernel */
                        L = 0;
                        T = TS;
                    }
                    else {
                        /* TT kernel */
                        L = A->mb;
                        T = TT;
                    }

                    for (m = 0; m < C->mt; m++) {
                        tempmm = m == C->mt-1 ? C->m-m*C->mb : C->mb;

                        node = C->get_rankof( C, m, n );
                        RUNTIME_data_migrate( sequence, C(m, p), node );
                        RUNTIME_data_migrate( sequence, C(m, n), node );

                        INSERT_TASK_ztpmqrt(
                            &options, side, trans,
                            tempmm, tempnn, tempkn, chameleon_min( L, tempmm ), ib, T->nb,
                            A(n, k),
                            T(n, k),
                            C(m, p),
                            C(m, n));
                    }
                    RUNTIME_data_flush( sequence, A(n, k) );
                    RUNTIME_data_flush( sequence, T(n, k) );
                }

                T = TS;
                for (i = 0; i < qrtree->getnbgeqrf(qrtree, k); i++) {
                    p = qrtree->getm(qrtree, k, i);

                    temppn   = p == C->nt-1 ? C->n - p * C->nb : C->nb;
                    tempkmin = chameleon_min(temppn, tempkn);

                    if ( genD ) {
                        int tempDpm = p == D->mt-1 ? D->m-p*D->mb : D->mb;

                        INSERT_TASK_zlacpy(
                            &options,
                            ChamLower, tempDpm, tempkmin, A->nb,
                            A(p, k),
                            D(p, k) );
#if defined(CHAMELEON_USE_CUDA)
                        INSERT_TASK_zlaset(
                            &options,
                            ChamUpper, tempDpm, tempkmin,
                            0., 1.,
                            D(p, k) );
#endif
                    }

                    for (m = 0; m < C->mt; m++) {
                        tempmm = m == C->mt-1 ? C->m-m*C->mb : C->mb;

                        RUNTIME_data_migrate( sequence, C(m, p),
                                              C->get_rankof( C, m, p ) );

                        INSERT_TASK_zunmqr(
                            &options, side, trans,
                            tempmm, temppn, tempkmin, ib, T->nb,
                            D(p, k),
                            T(p, k),
                            C(m, p));
                    }
                    RUNTIME_data_flush( sequence, D(p, k) );
                    RUNTIME_data_flush( sequence, T(p, k) );
                }
                RUNTIME_iteration_pop(chamctxt);
            }
        }
        /*
         *  ChamRight / ChamNoTrans
         */
        else {
            for (k = 0; k < KT; k++) {
                RUNTIME_iteration_push(chamctxt, k);

                tempkn = k == A->nt-1 ? A->n-k*A->nb : A->nb;

                T = TS;
                for (i = 0; i < qrtree->getnbgeqrf(qrtree, k); i++) {
                    p = qrtree->getm(qrtree, k, i);

                    temppn   = p == C->nt - 1 ? C->n - p * C->nb : C->nb;
                    tempkmin = chameleon_min( temppn, tempkn );

                    if ( genD ) {
                        int tempDpm = p == D->mt-1 ? D->m-p*D->mb : D->mb;

                        INSERT_TASK_zlacpy(
                            &options,
                            ChamLower, tempDpm, tempkmin, A->nb,
                            A(p, k),
                            D(p, k) );
#if defined(CHAMELEON_USE_CUDA)
                        INSERT_TASK_zlaset(
                            &options,
                            ChamUpper, tempDpm, tempkmin,
                            0., 1.,
                            D(p, k) );
#endif
                    }

                    for (m = 0; m < C->mt; m++) {
                        tempmm = m == C->mt-1 ? C->m-m*C->mb : C->mb;
                        INSERT_TASK_zunmqr(
                            &options, side, trans,
                            tempmm, temppn, tempkmin, ib, T->nb,
                            D(p, k),
                            T(p, k),
                            C(m, p));
                    }
                    RUNTIME_data_flush( sequence, D(p, k) );
                    RUNTIME_data_flush( sequence, T(p, k) );
                }

                /* Setting the order of tiles */
                nbtiles = libhqr_walk_stepk( qrtree, k, tiles );

                for (i = 0; i < nbtiles; i++) {
                    n = tiles[i];
                    p = qrtree->currpiv(qrtree, k, n);

                    tempnn = n == C->nt-1 ? C->n-n*C->nb : C->nb;

                    if( qrtree->gettype(qrtree, k, n) == LIBHQR_KILLED_BY_TS ) {
                        /* TS kernel */
                        L = 0;
                        T = TS;
                    }
                    else {
                        /* TT kernel */
                        L = A->mb;
                        T = TT;
                    }

                    for (m = 0; m < C->mt; m++) {
                        tempmm = m == C->mt-1 ? C->m-m*C->mb : C->mb;

                        node = C->get_rankof( C, m, n );
                        RUNTIME_data_migrate( sequence, C(m, p), node );
                        RUNTIME_data_migrate( sequence, C(m, n), node );

                        INSERT_TASK_ztpmqrt(
                            &options, side, trans,
                            tempmm, tempnn, tempkn, chameleon_min( L, tempmm ), ib, T->nb,
                            A(n, k),
                            T(n, k),
                            C(m, p),
                            C(m, n));
                    }
                    RUNTIME_data_flush( sequence, A(n, k) );
                    RUNTIME_data_flush( sequence, T(n, k) );
                }

                /* Restore the original location of the tiles */
                for (m = 0; m < C->mt; m++) {
                    RUNTIME_data_migrate( sequence, C( m, k ),
                                          C->get_rankof( C, m, k ) );
                }

                RUNTIME_iteration_pop(chamctxt);
            }
        }
    }

    free(tiles);
    RUNTIME_options_ws_free(&options);
    RUNTIME_options_finalize(&options, chamctxt);
}
