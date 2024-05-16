/**
 *
 * @file pztradd.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2021 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon ztradd parallel algorithm
 *
 * @version 1.0.0
 * @comment This file has been automatically generated
 *          from Plasma 2.5.0 for CHAMELEON 0.9.2
 * @author Emmanuel Agullo
 * @author Mathieu Faverge
 * @author Florent Pruvost
 * @date 2020-03-03
 * @precisions normal z -> s d c
 *
 */
#include "control/common.h"

#define A(m, n) A,  m,  n
#define B(m, n) B,  m,  n

/**
 *  Parallel tile matrix-matrix multiplication - dynamic scheduling
 */
void chameleon_pztradd(cham_uplo_t uplo, cham_trans_t trans,
                   CHAMELEON_Complex64_t alpha, CHAM_desc_t *A,
                   CHAMELEON_Complex64_t beta,  CHAM_desc_t *B,
                   RUNTIME_sequence_t *sequence, RUNTIME_request_t *request)
{
    CHAM_context_t *chamctxt;
    RUNTIME_option_t options;

    int tempmm, tempnn, tempmn, tempnm;
    int m, n;

    chamctxt = chameleon_context_self();
    if (sequence->status != CHAMELEON_SUCCESS) {
        return;
    }
    RUNTIME_options_init(&options, chamctxt, sequence, request);

    switch(uplo){
    case ChamLower:
        if (trans == ChamNoTrans) {
            for (n = 0; n < chameleon_min(B->mt,B->nt); n++) {
                tempnm = n == B->mt-1 ? B->m-n*B->mb : B->mb;
                tempnn = n == B->nt-1 ? B->n-n*B->nb : B->nb;

                INSERT_TASK_ztradd(
                    &options,
                    uplo, trans, tempnm, tempnn, B->mb,
                    alpha, A(n, n),
                    beta,  B(n, n));

                for (m = n+1; m < B->mt; m++) {
                    tempmm = m == B->mt-1 ? B->m-B->mb*m : B->nb;

                    INSERT_TASK_zgeadd(
                        &options,
                        trans, tempmm, tempnn, B->mb,
                        alpha, A(m, n),
                        beta,  B(m, n));
                }
            }
        }
        else {
            for (n = 0; n < chameleon_min(B->mt,B->nt); n++) {
                tempnm = n == B->mt-1 ? B->m-n*B->mb : B->mb;
                tempnn = n == B->nt-1 ? B->n-n*B->nb : B->nb;

                INSERT_TASK_ztradd(
                    &options,
                    uplo, trans, tempnm, tempnn, B->mb,
                    alpha, A(n, n),
                    beta,  B(n, n));

                for (m = n+1; m < B->mt; m++) {
                    tempmm = m == B->mt-1 ? B->m-B->mb*m : B->nb;

                    INSERT_TASK_zgeadd(
                        &options,
                        trans, tempmm, tempnn, B->mb,
                        alpha, A(n, m),
                        beta,  B(m, n));
                }
            }
        }
        break;
    case ChamUpper:
        if (trans == ChamNoTrans) {
            for (m = 0; m < chameleon_min(B->mt,B->nt); m++) {
                tempmm = m == B->mt-1 ? B->m-B->mb*m : B->nb;
                tempmn = m == B->nt-1 ? B->n-m*B->nb : B->nb;

                INSERT_TASK_ztradd(
                    &options,
                    uplo, trans, tempmm, tempmn, B->mb,
                    alpha, A(m, m),
                    beta,  B(m, m));

                for (n = m+1; n < B->nt; n++) {
                    tempnn = n == B->nt-1 ? B->n-n*B->nb : B->nb;

                    INSERT_TASK_zgeadd(
                        &options,
                        trans, tempmm, tempnn, B->mb,
                        alpha, A(m, n),
                        beta,  B(m, n));
                }
            }
        }
        else {
            for (m = 0; m < chameleon_min(B->mt,B->nt); m++) {
                tempmm = m == B->mt-1 ? B->m-B->mb*m : B->nb;
                tempmn = m == B->nt-1 ? B->n-m*B->nb : B->nb;

                INSERT_TASK_ztradd(
                    &options,
                    uplo, trans, tempmm, tempmn, B->mb,
                    alpha, A(m, m),
                    beta,  B(m, m));

                for (n = m+1; n < B->nt; n++) {
                    tempnn = n == B->nt-1 ? B->n-n*B->nb : B->nb;

                    INSERT_TASK_zgeadd(
                        &options,
                        trans, tempmm, tempnn, B->mb,
                        alpha, A(n, m),
                        beta,  B(m, n));
                }
            }
        }
        break;
    case ChamUpperLower:
    default:
        if (trans == ChamNoTrans) {
            for (m = 0; m < B->mt; m++) {
                tempmm = m == B->mt-1 ? B->m-B->mb*m : B->nb;

                for (n = 0; n < B->nt; n++) {
                    tempnn = n == B->nt-1 ? B->n-n*B->nb : B->nb;

                    INSERT_TASK_zgeadd(
                        &options,
                        trans, tempmm, tempnn, B->mb,
                        alpha, A(m, n),
                        beta,  B(m, n));
                }
            }
        }
        else {
            for (m = 0; m < B->mt; m++) {
                tempmm = m == B->mt-1 ? B->m-B->mb*m : B->nb;

                for (n = 0; n < B->nt; n++) {
                    tempnn = n == B->nt-1 ? B->n-n*B->nb : B->nb;

                    INSERT_TASK_zgeadd(
                        &options,
                        trans, tempmm, tempnn, B->mb,
                        alpha, A(n, m),
                        beta,  B(m, n));
                }
            }
        }
    }

    RUNTIME_options_finalize(&options, chamctxt);
}
