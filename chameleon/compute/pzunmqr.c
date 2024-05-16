/**
 *
 * @file pzunmqr.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2021 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zunmqr parallel algorithm
 *
 * @version 1.0.0
 * @comment This file has been automatically generated
 *          from Plasma 2.5.0 for CHAMELEON 0.9.2
 * @author Hatem Ltaief
 * @author Jakub Kurzak
 * @author Azzam Haidar
 * @author Mathieu Faverge
 * @author Emmanuel Agullo
 * @author Cedric Castagnede
 * @author Florent Pruvost
 * @author Raphael Boucherie
 * @author Samuel Thibault
 * @date 2020-03-03
 * @precisions normal z -> s d c
 *
 */
#include "control/common.h"

#define A(m,n) A,  m,  n
#define C(m,n) C,  m,  n
#define T(m,n) T,  m,  n
#define D(k)   D,  k,  k

/**
 *  Parallel application of Q using tile V - QR factorization - dynamic scheduling
 */
void chameleon_pzunmqr( int genD, cham_side_t side, cham_trans_t trans,
                        CHAM_desc_t *A, CHAM_desc_t *C, CHAM_desc_t *T, CHAM_desc_t *D,
                        RUNTIME_sequence_t *sequence, RUNTIME_request_t *request )
{
    CHAM_context_t *chamctxt;
    RUNTIME_option_t options;
    size_t ws_worker = 0;
    size_t ws_host = 0;

    int k, m, n;
    int tempkm, tempkn, tempkmin, tempmm, tempnn;
    int ib, KT, K;

    chamctxt = chameleon_context_self();
    if (sequence->status != CHAMELEON_SUCCESS) {
        return;
    }
    RUNTIME_options_init(&options, chamctxt, sequence, request);

    ib = CHAMELEON_IB;

    if (A->m > A->n) {
        KT = A->nt;
        K  = A->n;
    } else {
        KT = A->mt;
        K  = A->m;
    }

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
     * ztpmqrt = 2 * A->nb * ib
     */
    ws_worker = chameleon_max( ws_worker, ib * A->nb * 2 );
#endif

    ws_worker *= sizeof(CHAMELEON_Complex64_t);
    ws_host   *= sizeof(CHAMELEON_Complex64_t);

    RUNTIME_options_ws_alloc( &options, ws_worker, ws_host );

    if (side == ChamLeft ) {
        if (trans == ChamConjTrans) {
            /*
             *  ChamLeft / ChamConjTrans
             */
            for (k = 0; k < KT; k++) {
                RUNTIME_iteration_push(chamctxt, k);

                tempkm   = k == C->mt - 1 ? C->m - k * C->mb : C->mb;
                tempkmin = k == KT    - 1 ? K    - k * A->nb : A->nb;


                if ( genD ) {
                    int tempDkm = k == D->mt-1 ? D->m-k*D->mb : D->mb;

                    INSERT_TASK_zlacpy(
                        &options,
                        ChamLower, tempDkm, tempkmin, A->nb,
                        A(k, k),
                        D(k) );
#if defined(CHAMELEON_USE_CUDA)
                    INSERT_TASK_zlaset(
                        &options,
                        ChamUpper, tempDkm, tempkmin,
                        0., 1.,
                        D(k) );
#endif
                }
                for (n = 0; n < C->nt; n++) {
                    tempnn = n == C->nt-1 ? C->n-n*C->nb : C->nb;
                    INSERT_TASK_zunmqr(
                        &options,
                        side, trans,
                        tempkm, tempnn, tempkmin, ib, T->nb,
                        D(k),
                        T(k, k),
                        C(k, n));
                }

                RUNTIME_data_flush( sequence, D(k)    );
                RUNTIME_data_flush( sequence, T(k, k) );

                for (m = k+1; m < C->mt; m++) {
                    tempmm = m == C->mt-1 ? C->m-m*C->mb : C->mb;
                    for (n = 0; n < C->nt; n++) {
                        tempnn = n == C->nt-1 ? C->n-n*C->nb : C->nb;

                        RUNTIME_data_migrate( sequence, C(k, n),
                                              C->get_rankof( C, m, n ) );

                        /* TS kernel */
                        INSERT_TASK_ztpmqrt(
                            &options,
                            side, trans,
                            tempmm, tempnn, tempkmin, 0, ib, T->nb,
                            A(m, k),
                            T(m, k),
                            C(k, n),
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

                tempkm   = k == C->mt - 1 ? C->m - k * C->mb : C->mb;
                tempkmin = k == KT    - 1 ? K    - k * A->nb : A->nb;


                for (m = C->mt-1; m > k; m--) {
                    tempmm = m == C->mt-1 ? C->m-m*C->mb : C->mb;
                    for (n = 0; n < C->nt; n++) {
                        tempnn = n == C->nt-1 ? C->n-n*C->nb : C->nb;

                        RUNTIME_data_migrate( sequence, C(k, n),
                                              C->get_rankof( C, m, n ) );

                        /* TS kernel */
                        INSERT_TASK_ztpmqrt(
                            &options,
                            side, trans,
                            tempmm, tempnn, tempkmin, 0, ib, T->nb,
                            A(m, k),
                            T(m, k),
                            C(k, n),
                            C(m, n));
                    }
                    RUNTIME_data_flush( sequence, A(m, k) );
                    RUNTIME_data_flush( sequence, T(m, k) );
                }

                if ( genD ) {
                    int tempDkm = k == D->mt-1 ? D->m-k*D->mb : D->mb;

                    INSERT_TASK_zlacpy(
                        &options,
                        ChamLower, tempDkm, tempkmin, A->nb,
                        A(k, k),
                        D(k) );
#if defined(CHAMELEON_USE_CUDA)
                    INSERT_TASK_zlaset(
                        &options,
                        ChamUpper, tempDkm, tempkmin,
                        0., 1.,
                        D(k) );
#endif
                }
                for (n = 0; n < C->nt; n++) {
                    tempnn = n == C->nt-1 ? C->n-n*C->nb : C->nb;

                    RUNTIME_data_migrate( sequence, C(k, n),
                                          C->get_rankof( C, k, n ) );

                    INSERT_TASK_zunmqr(
                        &options,
                        side, trans,
                        tempkm, tempnn, tempkmin, ib, T->nb,
                        D(k),
                        T(k, k),
                        C(k, n));
                }
                RUNTIME_data_flush( sequence, D(k)    );
                RUNTIME_data_flush( sequence, T(k, k) );
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

                tempkn   = k == C->nt - 1 ? C->n - k * C->nb : C->nb;
                tempkmin = k == KT    - 1 ? K    - k * A->nb : A->nb;

                for (n = C->nt-1; n > k; n--) {
                    tempnn = n == C->nt-1 ? C->n-n*C->nb : C->nb;
                    for (m = 0; m < C->mt; m++) {
                        tempmm = m == C->mt-1 ? C->m-m*C->mb : C->mb;

                        RUNTIME_data_migrate( sequence, C(m, k),
                                              C->get_rankof( C, m, n ) );

                        /* TS kernel */
                        INSERT_TASK_ztpmqrt(
                            &options,
                            side, trans,
                            tempmm, tempnn, tempkmin, 0, ib, T->nb,
                            A(n, k),
                            T(n, k),
                            C(m, k),
                            C(m, n));
                    }

                    RUNTIME_data_flush( sequence, A(n, k) );
                    RUNTIME_data_flush( sequence, T(n, k) );
                }

                if ( genD ) {
                    int tempDkm = k == D->mt-1 ? D->m-k*D->mb : D->mb;

                    INSERT_TASK_zlacpy(
                        &options,
                        ChamLower, tempDkm, tempkmin, A->nb,
                        A(k, k),
                        D(k) );
#if defined(CHAMELEON_USE_CUDA)
                    INSERT_TASK_zlaset(
                        &options,
                        ChamUpper, tempDkm, tempkmin,
                        0., 1.,
                        D(k) );
#endif
                }
                for (m = 0; m < C->mt; m++) {
                    tempmm = m == C->mt-1 ? C->m-m*C->mb : C->mb;

                    RUNTIME_data_migrate( sequence, C(m, k),
                                          C->get_rankof( C, m, k ) );

                    INSERT_TASK_zunmqr(
                        &options,
                        side, trans,
                        tempmm, tempkn, tempkmin, ib, T->nb,
                        D(k),
                        T(k, k),
                        C(m, k));
                }

                RUNTIME_data_flush( sequence, D(k)    );
                RUNTIME_data_flush( sequence, T(k, k) );

                RUNTIME_iteration_pop(chamctxt);
            }
        }
        /*
         *  ChamRight / ChamNoTrans
         */
        else {
            for (k = 0; k < KT; k++) {
                RUNTIME_iteration_push(chamctxt, k);

                tempkn   = k == C->nt - 1 ? C->n - k * C->nb : C->nb;
                tempkmin = k == KT    - 1 ? K    - k * A->nb : A->nb;

                if ( genD ) {
                    int tempDkm = k == D->mt - 1 ? D->m - k * D->mb : D->mb;

                    INSERT_TASK_zlacpy(
                        &options,
                        ChamLower, tempDkm, tempkmin, A->nb,
                        A(k, k),
                        D(k) );
#if defined(CHAMELEON_USE_CUDA)
                    INSERT_TASK_zlaset(
                        &options,
                        ChamUpper, tempDkm, tempkmin,
                        0., 1.,
                        D(k) );
#endif
                }
                for (m = 0; m < C->mt; m++) {
                    tempmm = m == C->mt-1 ? C->m-m*C->mb : C->mb;
                    INSERT_TASK_zunmqr(
                        &options,
                        side, trans,
                        tempmm, tempkn, tempkmin, ib, T->nb,
                        D(k),
                        T(k, k),
                        C(m, k));
                }

                RUNTIME_data_flush( sequence, D(k)    );
                RUNTIME_data_flush( sequence, T(k, k) );

                for (n = k+1; n < C->nt; n++) {
                    tempnn = n == C->nt-1 ? C->n-n*C->nb : C->nb;
                    for (m = 0; m < C->mt; m++) {
                        tempmm = m == C->mt-1 ? C->m-m*C->mb : C->mb;

                        RUNTIME_data_migrate( sequence, C(m, k),
                                              C->get_rankof( C, m, n ) );

                        /* TS kernel */
                        INSERT_TASK_ztpmqrt(
                            &options,
                            side, trans,
                            tempmm, tempnn, tempkmin, 0, ib, T->nb,
                            A(n, k),
                            T(n, k),
                            C(m, k),
                            C(m, n));
                    }

                    RUNTIME_data_flush( sequence, A(n, k) );
                    RUNTIME_data_flush( sequence, T(n, k) );
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
}
