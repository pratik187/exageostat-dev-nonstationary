/**
 *
 * @file pzunmlqrh.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2021 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zunmlqrh parallel algorithm
 *
 * @version 1.0.0
 * @comment This file has been automatically generated
 *          from Plasma 2.5.0 for CHAMELEON 0.9.2
 * @author Hatem Ltaief
 * @author Jakub Kurzak
 * @author Dulceneia Becker
 * @author Mathieu Faverge
 * @author Emmanuel Agullo
 * @author Cedric Castagnede
 * @author Florent Pruvost
 * @author Samuel Thibault
 * @date 2020-03-03
 * @precisions normal z -> s d c
 *
 */
#include "control/common.h"

#define A(m,n)  A, (m), (n)
#define C(m,n)  C, (m), (n)
#define T(m,n)  T, (m), (n)
#define T2(m,n) T, (m), ((n)+A->nt)
#define D(m,n)  D, (m), (n)

/**
 *  Parallel application of Q using tile V - LQ factorization (reduction
 *  Householder) - dynamic scheduling
 */
void chameleon_pzunmlqrh( int genD, int BS, cham_side_t side, cham_trans_t trans,
                          CHAM_desc_t *A, CHAM_desc_t *C, CHAM_desc_t *T, CHAM_desc_t *D,
                          RUNTIME_sequence_t *sequence, RUNTIME_request_t *request )
{
    CHAM_context_t *chamctxt;
    RUNTIME_option_t options;
    size_t ws_worker = 0;
    size_t ws_host = 0;

    int k, m, n, p;
    int KT, RD, lastRD;
    int temppm, temppn, tempkm, tempnn, tempmm, tempkmin;
    int ib, node;

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
     * zunmlq  = A->nb * ib
     * ztpmlqt = A->nb * ib
     */
    ws_worker = A->nb * ib;

#if defined(CHAMELEON_USE_CUDA)
    /* Worker space
     *
     * zunmlq  =     A->nb * ib
     * ztpmlqt = 3 * A->nb * ib
     */
    ws_worker = chameleon_max( ws_worker, ib * A->nb * 3 );
#endif

    ws_worker *= sizeof(CHAMELEON_Complex64_t);
    ws_host   *= sizeof(CHAMELEON_Complex64_t);

    RUNTIME_options_ws_alloc( &options, ws_worker, ws_host );

    KT = chameleon_min( A->mt, A->nt );
    if (side == ChamLeft ) {
        if (trans == ChamNoTrans) {
            /*
             *  ChamLeft / ChamNoTrans
             */
            for (k = 0; k < KT; k++) {
                RUNTIME_iteration_push(chamctxt, k);

                tempkm = k == A->mt - 1 ? A->m - k * A->mb : A->mb;


                for (p = k; p < C->mt; p += BS) {

                    temppm   = p == C->mt-1 ? C->m - p * C->mb : C->mb;
                    tempkmin = chameleon_min( temppm, tempkm );


                    if ( genD ) {
                        int tempDpn = p == D->nt-1 ? D->n-p*D->nb : D->nb;

                        INSERT_TASK_zlacpy(
                            &options,
                            ChamUpper, tempkmin, tempDpn, A->nb,
                            A(k, p),
                            D(k, p) );
#if defined(CHAMELEON_USE_CUDA)
                        INSERT_TASK_zlaset(
                            &options,
                            ChamLower, tempkmin, tempDpn,
                            0., 1.,
                            D(k, p) );
#endif
                    }
                    for (n = 0; n < C->nt; n++) {
                        tempnn = n == C->nt-1 ? C->n-n*C->nb : C->nb;
                        INSERT_TASK_zunmlq(
                            &options,
                            side, trans,
                            temppm, tempnn, tempkmin, ib, T->nb,
                            D(k, p),
                            T(k, p),
                            C(p, n));
                    }
                    RUNTIME_data_flush( sequence, D(k, p) );
                    RUNTIME_data_flush( sequence, T(k, p) );

                    for (m = p+1; m < chameleon_min(p+BS, C->mt); m++) {
                        tempmm = m == C->mt-1 ? C->m-m*C->mb : C->mb;

                        for (n = 0; n < C->nt; n++) {
                            tempnn = n == C->nt-1 ? C->n-n*C->nb : C->nb;

                            node = C->get_rankof( C, m, n );
                            RUNTIME_data_migrate( sequence, C(p, n), node );
                            RUNTIME_data_migrate( sequence, C(m, n), node );

                            /* TS kernel */
                            INSERT_TASK_ztpmlqt(
                                &options, side, trans,
                                tempmm, tempnn, tempkm, 0, ib, T->nb,
                                A(k, m),
                                T(k, m),
                                C(p, n),
                                C(m, n));
                        }
                        RUNTIME_data_flush( sequence, A(k, m) );
                        RUNTIME_data_flush( sequence, T(k, m) );
                    }
                }
                for (RD = BS; RD < C->mt-k; RD *= 2) {
                    for (p = k; p+RD < C->mt; p += 2*RD) {
                        m = p+RD;

                        tempmm = m == C->mt-1 ? C->m-m*C->mb : C->mb;

                        for (n = 0; n < C->nt; n++) {
                            tempnn = n == C->nt-1 ? C->n-n*C->nb : C->nb;

                            node = C->get_rankof( C, m, n );
                            RUNTIME_data_migrate( sequence, C(p, n), node );
                            RUNTIME_data_migrate( sequence, C(m, n), node );

                            /* TT kernel */
                            INSERT_TASK_ztpmlqt(
                                &options,
                                side, trans,
                                tempmm, tempnn, tempkm, tempnn, ib, T->nb,
                                A (k, m),
                                T2(k, m),
                                C (p, n),
                                C (m, n));
                        }
                        RUNTIME_data_flush( sequence, A (k, m) );
                        RUNTIME_data_flush( sequence, T2(k, m) );
                    }
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
         *  ChamLeft / ChamConjTrans
         */
        else {
            for (k = KT-1; k >= 0; k--) {
                RUNTIME_iteration_push(chamctxt, k);

                tempkm = k == A->mt-1 ? A->m-k*A->mb : A->mb;

                lastRD = 0;
                for (RD = BS; RD < C->mt-k; RD *= 2)
                    lastRD = RD;
                for (RD = lastRD; RD >= BS; RD /= 2) {
                    for (p = k; p+RD < C->mt; p += 2*RD) {
                        m = p+RD;

                        tempmm = m == C->mt-1 ? C->m-m*C->mb : C->mb;

                        for (n = 0; n < C->nt; n++) {
                            tempnn = n == C->nt-1 ? C->n-n*C->nb : C->nb;

                            node = C->get_rankof( C, m, n );
                            RUNTIME_data_migrate( sequence, C(p, n), node );
                            RUNTIME_data_migrate( sequence, C(m, n), node );

                            /* TT kernel */
                            INSERT_TASK_ztpmlqt(
                                &options, side, trans,
                                tempmm, tempnn, tempkm, tempnn, ib, T->nb,
                                A (k, m),
                                T2(k, m),
                                C (p, n),
                                C (m, n));
                        }
                        RUNTIME_data_flush( sequence, A (k, m) );
                        RUNTIME_data_flush( sequence, T2(k, m) );
                    }
                }
                for (p = k; p < C->mt; p += BS) {

                    for (m = chameleon_min(p+BS, C->mt)-1; m > p; m--) {
                        tempmm = m == C->mt-1 ? C->m-m*C->mb : C->mb;

                        for (n = 0; n < C->nt; n++) {
                            tempnn = n == C->nt-1 ? C->n-n*C->nb : C->nb;

                            node = C->get_rankof( C, m, n );
                            RUNTIME_data_migrate( sequence, C(p, n), node );
                            RUNTIME_data_migrate( sequence, C(m, n), node );

                            /* TS kernel */
                            INSERT_TASK_ztpmlqt(
                                &options, side, trans,
                                tempmm, tempnn, tempkm, 0, ib, T->nb,
                                A(k, m),
                                T(k, m),
                                C(p, n),
                                C(m, n));
                        }
                        RUNTIME_data_flush( sequence, A(k, m) );
                        RUNTIME_data_flush( sequence, T(k, m) );
                    }

                    temppm   = p == C->mt-1 ? C->m-p*C->mb : C->mb;
                    tempkmin = chameleon_min( temppm, tempkm );

                    if ( genD ) {
                        int tempDpn = p == D->nt-1 ? D->n-p*D->nb : D->nb;

                        INSERT_TASK_zlacpy(
                            &options,
                            ChamUpper, tempkmin, tempDpn, A->nb,
                            A(k, p),
                            D(k, p) );
#if defined(CHAMELEON_USE_CUDA)
                        INSERT_TASK_zlaset(
                            &options,
                            ChamLower, tempkmin, tempDpn,
                            0., 1.,
                            D(k, p) );
#endif
                    }

                    for (n = 0; n < C->nt; n++) {
                        tempnn = n == C->nt-1 ? C->n-n*C->nb : C->nb;

                        RUNTIME_data_migrate( sequence, C(p, n),
                                              C->get_rankof( C, p, n ) );

                        INSERT_TASK_zunmlq(
                            &options, side, trans,
                            temppm, tempnn, tempkmin, ib, T->nb,
                            D(k, p),
                            T(k, p),
                            C(p, n));
                    }
                    RUNTIME_data_flush( sequence, D(k, p) );
                    RUNTIME_data_flush( sequence, T(k, p) );
                }
                RUNTIME_iteration_pop(chamctxt);
            }
        }
    }
    /*
     *  ChamRight / ChamNoTrans
     */
    else {
        if (trans == ChamNoTrans) {
            for (k = KT-1; k >= 0; k--) {
                RUNTIME_iteration_push(chamctxt, k);

                tempkm = k == A->mt-1 ? A->m-k*A->mb : A->mb;
                lastRD = 0;
                for (RD = BS; RD < C->nt-k; RD *= 2)
                    lastRD = RD;
                for (RD = lastRD; RD >= BS; RD /= 2) {
                    for (p = k; p+RD < C->nt; p += 2*RD) {
                        n = p+RD;

                        tempnn = n == C->nt-1 ? C->n-n*C->nb : C->nb;

                        for (m = 0; m < C->mt; m++) {
                            tempmm = m == C->mt-1 ? C->m-m*C->mb : C->mb;

                            node = C->get_rankof( C, m, n );
                            RUNTIME_data_migrate( sequence, C(m, p), node );
                            RUNTIME_data_migrate( sequence, C(m, n), node );

                            /* TT kernel */
                            INSERT_TASK_ztpmlqt(
                                &options, side, trans,
                                tempmm, tempnn, tempkm, tempnn, ib, T->nb,
                                A (k, n),
                                T2(k, n),
                                C (m, p),
                                C (m, n));
                        }
                        RUNTIME_data_flush( sequence, A (k, n) );
                        RUNTIME_data_flush( sequence, T2(k, n) );
                    }
                }
                for (p = k; p < C->nt; p += BS) {

                    for (n = chameleon_min(p+BS, C->nt)-1; n > p; n--) {

                        tempnn = n == C->nt-1 ? C->n-n*C->nb : C->nb;

                        for (m = 0; m < C->mt; m++) {
                            tempmm = m == C->mt-1 ? C->m-m*C->mb : C->mb;

                            node = C->get_rankof( C, m, n );
                            RUNTIME_data_migrate( sequence, C(m, p), node );
                            RUNTIME_data_migrate( sequence, C(m, n), node );

                            /* TS kernel */
                            INSERT_TASK_ztpmlqt(
                                &options, side, trans,
                                tempmm, tempnn, tempkm, 0, ib, T->nb,
                                A(k, n),
                                T(k, n),
                                C(m, p),
                                C(m, n));
                        }
                        RUNTIME_data_flush( sequence, A(k, n) );
                        RUNTIME_data_flush( sequence, T(k, n) );
                    }

                    temppn   = p == C->nt-1 ? C->n - p * C->nb : C->nb;
                    tempkmin = chameleon_min( temppn, tempkm );

                    if ( genD ) {
                        int tempDpn = p == D->nt-1 ? D->n-p*D->nb : D->nb;

                        INSERT_TASK_zlacpy(
                            &options,
                            ChamUpper, tempkmin, tempDpn, A->nb,
                            A(k, p),
                            D(k, p) );
#if defined(CHAMELEON_USE_CUDA)
                        INSERT_TASK_zlaset(
                            &options,
                            ChamLower, tempkmin, tempDpn,
                            0., 1.,
                            D(k, p) );
#endif
                    }

                    for (m = 0; m < C->mt; m++) {
                        tempmm = m == C->mt-1 ? C->m-m*C->mb : C->mb;

                        RUNTIME_data_migrate( sequence, C(m, p),
                                              C->get_rankof( C, m, p ) );

                        INSERT_TASK_zunmlq(
                            &options, side, trans,
                            tempmm, temppn, tempkmin, ib, T->nb,
                            D(k, p),
                            T(k, p),
                            C(m, p));
                    }
                    RUNTIME_data_flush( sequence, D(k, p) );
                    RUNTIME_data_flush( sequence, T(k, p) );
                }
                RUNTIME_iteration_pop(chamctxt);
            }
        }
        /*
         *  ChamRight / ChamConjTrans
         */
        else {
            for (k = 0; k < KT; k++) {
                RUNTIME_iteration_push(chamctxt, k);

                tempkm = k == A->mt-1 ? A->m-k*A->mb : A->mb;

                for (p = k; p < C->nt; p += BS) {
                    temppn   = p == C->nt - 1 ? C->n - p * C->nb : C->nb;
                    tempkmin = chameleon_min( temppn, tempkm );

                    if ( genD ) {
                        int tempDpn = p == D->nt-1 ? D->n-p*D->nb : D->nb;

                        INSERT_TASK_zlacpy(
                            &options,
                            ChamUpper, tempkmin, tempDpn, A->nb,
                            A(k, p),
                            D(k, p) );
#if defined(CHAMELEON_USE_CUDA)
                        INSERT_TASK_zlaset(
                            &options,
                            ChamLower, tempkmin, tempDpn,
                            0., 1.,
                            D(k, p) );
#endif
                    }

                    for (m = 0; m < C->mt; m++) {
                        tempmm = m == C->mt-1 ? C->m-m*C->mb : C->mb;
                        INSERT_TASK_zunmlq(
                            &options, side, trans,
                            tempmm, temppn, tempkmin, ib, T->nb,
                            D(k, p),
                            T(k, p),
                            C(m, p));
                    }
                    RUNTIME_data_flush( sequence, D(k, p) );
                    RUNTIME_data_flush( sequence, T(k, p) );

                    for (n = p+1; n < chameleon_min(p+BS, C->nt); n++) {
                        tempnn = n == C->nt-1 ? C->n-n*C->nb : C->nb;
                        for (m = 0; m < C->mt; m++) {
                            tempmm = m == C->mt-1 ? C->m-m*C->mb : C->mb;

                            node = C->get_rankof( C, m, n );
                            RUNTIME_data_migrate( sequence, C(m, p), node );
                            RUNTIME_data_migrate( sequence, C(m, n), node );

                            /* TS kernel */
                            INSERT_TASK_ztpmlqt(
                                &options, side, trans,
                                tempmm, tempnn, tempkm, 0, ib, T->nb,
                                A(k, n),
                                T(k, n),
                                C(m, p),
                                C(m, n));
                        }
                        RUNTIME_data_flush( sequence, A(k, n) );
                        RUNTIME_data_flush( sequence, T(k, n) );
                    }
                }
                for (RD = BS; RD < C->nt-k; RD *= 2) {
                    for (p = k; p+RD < C->nt; p += 2*RD) {
                        n = p + RD;
                        tempnn = n == C->mt-1 ? C->m-n*C->mb : C->mb;

                        for (m = 0; m < C->mt; m++) {
                            tempmm = m == C->mt-1 ? C->m-m*C->mb : C->mb;

                            node = C->get_rankof( C, m, n );
                            RUNTIME_data_migrate( sequence, C(m, p), node );
                            RUNTIME_data_migrate( sequence, C(m, n), node );

                            /* TT kernel */
                            INSERT_TASK_ztpmlqt(
                                &options,
                                side, trans,
                                tempmm, tempnn, tempkm, tempnn, ib, T->nb,
                                A (k, n),
                                T2(k, n),
                                C (m, p),
                                C (m, n));
                        }
                        RUNTIME_data_flush( sequence, A (k, n) );
                        RUNTIME_data_flush( sequence, T2(k, n) );
                    }
                }

                /* Restore the original location of the tiles */
                for (m = 0; m < C->mt; m++) {
                    RUNTIME_data_migrate( sequence, C(m, k),
                                          C->get_rankof( C, m, k ) );
                }

                RUNTIME_iteration_pop(chamctxt);
            }
        }
    }

    RUNTIME_options_ws_free(&options);
    RUNTIME_options_finalize(&options, chamctxt);
    (void)D;
}
