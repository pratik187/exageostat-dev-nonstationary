/**
 *
 * @file quark/codelet_zher2k.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2021 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zher2k Quark codelet
 *
 * @version 1.1.0
 * @author Hatem Ltaief
 * @author Jakub Kurzak
 * @author Mathieu Faverge
 * @author Emmanuel Agullo
 * @author Cedric Castagnede
 * @author Florent Pruvost
 * @date 2020-06-29
 * @precisions normal z -> c
 *
 */
#include "chameleon_quark.h"
#include "chameleon/tasks_z.h"
#include "coreblas/coreblas_ztile.h"

void CORE_zher2k_quark(Quark *quark)
{
    cham_uplo_t uplo;
    cham_trans_t trans;
    int n;
    int k;
    CHAMELEON_Complex64_t alpha;
    CHAM_tile_t *tileA;
    CHAM_tile_t *tileB;
    double beta;
    CHAM_tile_t *tileC;

    quark_unpack_args_9(quark, uplo, trans, n, k, alpha, tileA, tileB, beta, tileC);
    TCORE_zher2k(uplo, trans,
                n, k, alpha, tileA, tileB, beta, tileC);
}

void
INSERT_TASK_zher2k( const RUNTIME_option_t *options,
                    cham_uplo_t uplo, cham_trans_t trans,
                    int n, int k, int nb,
                    CHAMELEON_Complex64_t alpha, const CHAM_desc_t *A, int Am, int An,
                                                 const CHAM_desc_t *B, int Bm, int Bn,
                    double beta,                 const CHAM_desc_t *C, int Cm, int Cn )
{
    if ( alpha == 0. ) {
        return INSERT_TASK_zlascal( options, uplo, n, n, nb,
                                    beta, C, Cm, Cn );
    }

    quark_option_t *opt = (quark_option_t*)(options->schedopt);
    int accessC = ( beta == 0. ) ? OUTPUT : INOUT;

    DAG_CORE_HER2K;
    QUARK_Insert_Task(opt->quark, CORE_zher2k_quark, (Quark_Task_Flags*)opt,
        sizeof(int),                &uplo,      VALUE,
        sizeof(int),                &trans,     VALUE,
        sizeof(int),                        &n,         VALUE,
        sizeof(int),                        &k,         VALUE,
        sizeof(CHAMELEON_Complex64_t),         &alpha,     VALUE,
        sizeof(void*), RTBLKADDR(A, CHAMELEON_Complex64_t, Am, An),                 INPUT,
        sizeof(void*), RTBLKADDR(B, CHAMELEON_Complex64_t, Bm, Bn),                 INPUT,
        sizeof(double),                     &beta,      VALUE,
        sizeof(void*), RTBLKADDR(C, CHAMELEON_Complex64_t, Cm, Cn),                 accessC,
        0);
}
