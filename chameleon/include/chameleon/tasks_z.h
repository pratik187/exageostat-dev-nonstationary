/**
 *
 * @file chameleon_tasks_z.h
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2021 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon CHAMELEON_Complex64_t elementary tasks header
 *
 * @version 1.1.0
 * @comment This file has been automatically generated
 *          from Plasma 2.5.0 for CHAMELEON 0.9.2
 * @author Jakub Kurzak
 * @author Hatem Ltaief
 * @author Azzam Haidar
 * @author Mathieu Faverge
 * @author Emmanuel Agullo
 * @author Cedric Castagnede
 * @author Florent Pruvost
 * @date 2020-10-12
 * @precisions normal z -> c d s
 *
 */
#ifndef _chameleon_tasks_z_h_
#define _chameleon_tasks_z_h_

/**
 *  Declarations of QUARK wrappers (called by CHAMELEON) - alphabetical order
 */
void INSERT_TASK_dlag2z( const RUNTIME_option_t *options,
                         cham_uplo_t uplo, int m, int n,
                         const CHAM_desc_t *A, int Am, int An,
                         const CHAM_desc_t *B, int Bm, int Bn );
void INSERT_TASK_dzasum( const RUNTIME_option_t *options,
                         cham_store_t storev, cham_uplo_t uplo, int M, int N,
                         const CHAM_desc_t *A, int Am, int An,
                         const CHAM_desc_t *B, int Bm, int Bn );
void INSERT_TASK_zaxpy( const RUNTIME_option_t *options,
                        int M, CHAMELEON_Complex64_t alpha,
                        const CHAM_desc_t *A, int Am, int An, int incA,
                        const CHAM_desc_t *B, int Bm, int Bn, int incB );
void INSERT_TASK_zbuild( const RUNTIME_option_t *options,
                         const CHAM_desc_t *A, int Am, int An,
                         void *user_data, void* user_build_callback );
void INSERT_TASK_zgeadd( const RUNTIME_option_t *options,
                         cham_trans_t trans, int m, int n, int nb,
                         CHAMELEON_Complex64_t alpha, const CHAM_desc_t *A, int Am, int An,
                         CHAMELEON_Complex64_t beta,  const CHAM_desc_t *B, int Bm, int Bn );
void INSERT_TASK_zgelqt( const RUNTIME_option_t *options,
                         int m, int n, int ib, int nb,
                         const CHAM_desc_t *A, int Am, int An,
                         const CHAM_desc_t *T, int Tm, int Tn );
void INSERT_TASK_zgemv( const RUNTIME_option_t *options,
                        cham_trans_t trans, int m, int n,
                        CHAMELEON_Complex64_t alpha, const CHAM_desc_t *A, int Am, int An,
                                                     const CHAM_desc_t *X, int Xm, int Xn, int incX,
                        CHAMELEON_Complex64_t beta,  const CHAM_desc_t *Y, int Ym, int Yn, int incY );
void INSERT_TASK_zgemm( const RUNTIME_option_t *options,
                        cham_trans_t transA, cham_trans_t transB,
                        int m, int n, int k, int nb,
                        CHAMELEON_Complex64_t alpha, const CHAM_desc_t *A, int Am, int An,
                        const CHAM_desc_t *B, int Bm, int Bn,
                        CHAMELEON_Complex64_t beta, const CHAM_desc_t *C, int Cm, int Cn );
void INSERT_TASK_zgeqrt( const RUNTIME_option_t *options,
                         int m, int n, int ib, int nb,
                         const CHAM_desc_t *A, int Am, int An,
                         const CHAM_desc_t *T, int Tm, int Tn );
void INSERT_TASK_zgessm( const RUNTIME_option_t *options,
                         int m, int n, int k, int ib, int nb,
                         int *IPIV,
                         const CHAM_desc_t *L, int Lm, int Ln,
                         const CHAM_desc_t *D, int Dm, int Dn,
                         const CHAM_desc_t *A, int Am, int An );
void INSERT_TASK_zgessq( const RUNTIME_option_t *options,
                         cham_store_t storev, int m, int n,
                         const CHAM_desc_t *A, int Am, int An,
                         const CHAM_desc_t *SCALESUMSQ, int SCALESUMSQm, int SCALESUMSQn );
void INSERT_TASK_zgetrf( const RUNTIME_option_t *options,
                         int m, int n, int nb,
                         const CHAM_desc_t *A, int Am, int An,
                         int *IPIV,
                         cham_bool_t check_info, int iinfo );
void INSERT_TASK_zgetrf_incpiv( const RUNTIME_option_t *options,
                                int m, int n, int ib, int nb,
                                const CHAM_desc_t *A, int Am, int An,
                                const CHAM_desc_t *L, int Lm, int Ln,
                                int *IPIV,
                                cham_bool_t check_info, int iinfo );
void INSERT_TASK_zgetrf_nopiv( const RUNTIME_option_t *options,
                               int m, int n, int ib, int nb,
                               const CHAM_desc_t *A, int Am, int An, int iinfo );
void INSERT_TASK_zhe2ge( const RUNTIME_option_t *options,
                         cham_uplo_t uplo,
                         int m, int n, int mb,
                         const CHAM_desc_t *A, int Am, int An,
                         const CHAM_desc_t *B, int Bm, int Bn );
void INSERT_TASK_zhemm( const RUNTIME_option_t *options,
                        cham_side_t side, cham_uplo_t uplo,
                        int m, int n, int nb,
                        CHAMELEON_Complex64_t alpha, const CHAM_desc_t *A, int Am, int An,
                        const CHAM_desc_t *B, int Bm, int Bn,
                        CHAMELEON_Complex64_t beta, const CHAM_desc_t *C, int Cm, int Cn );
void INSERT_TASK_zher2k( const RUNTIME_option_t *options,
                         cham_uplo_t uplo, cham_trans_t trans,
                         int n, int k, int nb,
                         CHAMELEON_Complex64_t alpha, const CHAM_desc_t *A, int Am, int An,
                         const CHAM_desc_t *B, int Bm, int Bn,
                         double beta, const CHAM_desc_t *C, int Cm, int Cn );
void INSERT_TASK_zherfb( const RUNTIME_option_t *options,
                         cham_uplo_t uplo,
                         int n, int k, int ib, int nb,
                         const CHAM_desc_t *A, int Am, int An,
                         const CHAM_desc_t *T, int Tm, int Tn,
                         const CHAM_desc_t *C, int Cm, int Cn );
void INSERT_TASK_zherk( const RUNTIME_option_t *options,
                        cham_uplo_t uplo, cham_trans_t trans,
                        int n, int k, int nb,
                        double alpha, const CHAM_desc_t *A, int Am, int An,
                        double beta, const CHAM_desc_t *C, int Cm, int Cn );
void INSERT_TASK_zhessq( const RUNTIME_option_t *options,
                         cham_store_t storev, cham_uplo_t uplo, int n,
                         const CHAM_desc_t *A, int Am, int An,
                         const CHAM_desc_t *SCALESUMSQ, int SCALESUMSQm, int SCALESUMSQn );
void INSERT_TASK_zlacpy( const RUNTIME_option_t *options,
                         cham_uplo_t uplo, int m, int n, int mb,
                         const CHAM_desc_t *A, int Am, int An,
                         const CHAM_desc_t *B, int Bm, int Bn );
void INSERT_TASK_zlacpyx( const RUNTIME_option_t *options,
                          cham_uplo_t uplo, int m, int n, int mb,
                          int displA, const CHAM_desc_t *A, int Am, int An,
                          int displB, const CHAM_desc_t *B, int Bm, int Bn );
void INSERT_TASK_zlange( const RUNTIME_option_t *options,
                         cham_normtype_t norm, int M, int N, int NB,
                         const CHAM_desc_t *A, int Am, int An,
                         const CHAM_desc_t *B, int Bm, int Bn );
void INSERT_TASK_zlange_max( const RUNTIME_option_t *options,
                             const CHAM_desc_t *A, int Am, int An,
                             const CHAM_desc_t *B, int Bm, int Bn );
void INSERT_TASK_zlanhe( const RUNTIME_option_t *options,
                         cham_normtype_t norm, cham_uplo_t uplo, int N, int NB,
                         const CHAM_desc_t *A, int Am, int An,
                         const CHAM_desc_t *B, int Bm, int Bn );
void INSERT_TASK_zlansy( const RUNTIME_option_t *options,
                         cham_normtype_t norm, cham_uplo_t uplo, int N, int NB,
                         const CHAM_desc_t *A, int Am, int An,
                         const CHAM_desc_t *B, int Bm, int Bn );
void INSERT_TASK_zlantr( const RUNTIME_option_t *options,
                         cham_normtype_t norm, cham_uplo_t uplo, cham_diag_t diag,
                         int M, int N, int NB,
                         const CHAM_desc_t *A, int Am, int An,
                         const CHAM_desc_t *B, int Bm, int Bn );
void INSERT_TASK_zlascal( const RUNTIME_option_t *options,
                          cham_uplo_t uplo,
                          int m, int n, int nb,
                          CHAMELEON_Complex64_t alpha,
                          const CHAM_desc_t *A, int Am, int An );
void INSERT_TASK_zlaset( const RUNTIME_option_t *options,
                         cham_uplo_t uplo, int n1, int n2,
                         CHAMELEON_Complex64_t alpha, CHAMELEON_Complex64_t beta,
                         const CHAM_desc_t *tileA, int tileAm, int tileAn );
void INSERT_TASK_zlaset2( const RUNTIME_option_t *options,
                          cham_uplo_t uplo, int n1, int n2, CHAMELEON_Complex64_t alpha,
                          const CHAM_desc_t *tileA, int tileAm, int tileAn );
void INSERT_TASK_zlatro( const RUNTIME_option_t *options,
                         cham_uplo_t uplo, cham_trans_t trans, int m, int n, int mb,
                         const CHAM_desc_t *A, int Am, int An,
                         const CHAM_desc_t *B, int Bm, int Bn );
void INSERT_TASK_zlauum( const RUNTIME_option_t *options,
                         cham_uplo_t uplo, int n, int nb,
                         const CHAM_desc_t *A, int Am, int An );
void INSERT_TASK_zplghe( const RUNTIME_option_t *options,
                         double bump, int m, int n, const CHAM_desc_t *A, int Am, int An,
                         int bigM, int m0, int n0, unsigned long long int seed );
void INSERT_TASK_zplgsy( const RUNTIME_option_t *options,
                         CHAMELEON_Complex64_t bump, int m, int n, const CHAM_desc_t *A, int Am, int An,
                         int bigM, int m0, int n0, unsigned long long int seed );
void INSERT_TASK_zplrnt( const RUNTIME_option_t *options,
                         int m, int n, const CHAM_desc_t *A, int Am, int An,
                         int bigM, int m0, int n0, unsigned long long int seed );
void INSERT_TASK_zplssq( const RUNTIME_option_t *options,
                         cham_store_t storev, int M, int N,
                         const CHAM_desc_t *SCALESUMSQ, int SCALESUMSQm, int SCALESUMSQn,
                         const CHAM_desc_t *SCLSSQ,     int SCLSSQm,     int SCLSSQn );
void INSERT_TASK_zplssq2( const RUNTIME_option_t *options, int N,
                          const CHAM_desc_t *RESULT, int RESULTm, int RESULTn );
void INSERT_TASK_zpotrf( const RUNTIME_option_t *options,
                         cham_uplo_t uplo, int n, int nb,
                         const CHAM_desc_t *A, int Am, int An,
                         int iinfo );
void INSERT_TASK_zssssm( const RUNTIME_option_t *options,
                         int m1, int n1, int m2, int n2, int k, int ib, int nb,
                         const CHAM_desc_t *A1, int A1m, int A1n,
                         const CHAM_desc_t *A2, int A2m, int A2n,
                         const CHAM_desc_t *L1, int L1m, int L1n,
                         const CHAM_desc_t *L2, int L2m, int L2n,
                         const int *IPIV );
void INSERT_TASK_zsymm( const RUNTIME_option_t *options,
                        cham_side_t side, cham_uplo_t uplo,
                        int m, int n, int nb,
                        CHAMELEON_Complex64_t alpha, const CHAM_desc_t *A, int Am, int An,
                        const CHAM_desc_t *B, int Bm, int Bn,
                        CHAMELEON_Complex64_t beta, const CHAM_desc_t *C, int Cm, int Cn );
void INSERT_TASK_zsyr2k( const RUNTIME_option_t *options,
                         cham_uplo_t uplo, cham_trans_t trans,
                         int n, int k, int nb,
                         CHAMELEON_Complex64_t alpha, const CHAM_desc_t *A, int Am, int An,
                         const CHAM_desc_t *B, int Bm, int Bn,
                         CHAMELEON_Complex64_t beta, const CHAM_desc_t *C, int Cm, int Cn );
void INSERT_TASK_zsyrk( const RUNTIME_option_t *options,
                        cham_uplo_t uplo, cham_trans_t trans,
                        int n, int k, int nb,
                        CHAMELEON_Complex64_t alpha, const CHAM_desc_t *A, int Am, int An,
                        CHAMELEON_Complex64_t beta, const CHAM_desc_t *C, int Cm, int Cn );
void INSERT_TASK_zsyssq( const RUNTIME_option_t *options,
                         cham_store_t storev, cham_uplo_t uplo, int n,
                         const CHAM_desc_t *A, int Am, int An,
                         const CHAM_desc_t *SCALESUMSQ, int SCALESUMSQm, int SCALESUMSQn );
void INSERT_TASK_zsytrf_nopiv( const RUNTIME_option_t *options,
                               cham_uplo_t uplo, int n, int nb,
                               const CHAM_desc_t *A, int Am, int An,
                               int iinfo );
void INSERT_TASK_ztplqt( const RUNTIME_option_t *options,
                         int m, int n, int l, int ib, int nb,
                         const CHAM_desc_t *A, int Am, int An,
                         const CHAM_desc_t *B, int Bm, int Bn,
                         const CHAM_desc_t *T, int Tm, int Tn );
void INSERT_TASK_ztpmlqt( const RUNTIME_option_t *options,
                          cham_side_t side, cham_trans_t trans,
                          int M, int N, int K, int L, int ib, int nb,
                          const CHAM_desc_t *V, int Vm, int Vn,
                          const CHAM_desc_t *T, int Tm, int Tn,
                          const CHAM_desc_t *A, int Am, int An,
                          const CHAM_desc_t *B, int Bm, int Bn );
void INSERT_TASK_ztpmqrt( const RUNTIME_option_t *options,
                          cham_side_t side, cham_trans_t trans,
                          int m, int n, int k, int l, int ib, int nb,
                          const CHAM_desc_t *V, int Vm, int Vn,
                          const CHAM_desc_t *T, int Tm, int Tn,
                          const CHAM_desc_t *A, int Am, int An,
                          const CHAM_desc_t *B, int Bm, int Bn );
void INSERT_TASK_ztpqrt( const RUNTIME_option_t *options,
                         int m, int n, int l, int ib, int nb,
                         const CHAM_desc_t *A, int Am, int An,
                         const CHAM_desc_t *B, int Bm, int Bn,
                         const CHAM_desc_t *T, int Tm, int Tn );
void INSERT_TASK_ztradd( const RUNTIME_option_t *options,
                         cham_uplo_t uplo, cham_trans_t trans, int m, int n, int nb,
                         CHAMELEON_Complex64_t alpha, const CHAM_desc_t *A, int Am, int An,
                         CHAMELEON_Complex64_t beta,  const CHAM_desc_t *B, int Bm, int Bn );
void INSERT_TASK_ztrasm( const RUNTIME_option_t *options,
                         cham_store_t storev, cham_uplo_t uplo, cham_diag_t diag, int M, int N,
                         const CHAM_desc_t *A, int Am, int An,
                         const CHAM_desc_t *B, int Bm, int Bn );
void INSERT_TASK_ztrmm( const RUNTIME_option_t *options,
                        cham_side_t side, cham_uplo_t uplo, cham_trans_t transA, cham_diag_t diag,
                        int m, int n, int nb,
                        CHAMELEON_Complex64_t alpha, const CHAM_desc_t *A, int Am, int An,
                        const CHAM_desc_t *B, int Bm, int Bn );
void INSERT_TASK_ztrsm( const RUNTIME_option_t *options,
                        cham_side_t side, cham_uplo_t uplo, cham_trans_t transA, cham_diag_t diag,
                        int m, int n, int nb,
                        CHAMELEON_Complex64_t alpha, const CHAM_desc_t *A, int Am, int An,
                        const CHAM_desc_t *B, int Bm, int Bn );
void INSERT_TASK_ztrssq( const RUNTIME_option_t *options,
                         cham_uplo_t uplo, cham_diag_t diag,
                         int m, int n,
                         const CHAM_desc_t *A, int Am, int An,
                         const CHAM_desc_t *SCALESUMSQ, int SCALESUMSQm, int SCALESUMSQn );
void INSERT_TASK_ztrtri( const RUNTIME_option_t *options,
                         cham_uplo_t uplo, cham_diag_t diag, int n, int nb,
                         const CHAM_desc_t *A, int Am, int An,

                         int iinfo );
void INSERT_TASK_ztsmlq_hetra1( const RUNTIME_option_t *options,
                                cham_side_t side, cham_trans_t trans,
                                int m1, int n1, int m2, int n2, int k, int ib, int nb,
                                const CHAM_desc_t *A1, int A1m, int A1n,
                                const CHAM_desc_t *A2, int A2m, int A2n,
                                const CHAM_desc_t *V, int Vm, int Vn,
                                const CHAM_desc_t *T, int Tm, int Tn );
void INSERT_TASK_ztsmqr_hetra1( const RUNTIME_option_t *options,
                                cham_side_t side, cham_trans_t trans,
                                int m1, int n1, int m2, int n2, int k, int ib, int nb,
                                const CHAM_desc_t *A1, int A1m, int A1n,
                                const CHAM_desc_t *A2, int A2m, int A2n,
                                const CHAM_desc_t *V, int Vm, int Vn,
                                const CHAM_desc_t *T, int Tm, int Tn );
void INSERT_TASK_ztstrf( const RUNTIME_option_t *options,
                         int m, int n, int ib, int nb,
                         const CHAM_desc_t *U, int Um, int Un,
                         const CHAM_desc_t *A, int Am, int An,
                         const CHAM_desc_t *L, int Lm, int Ln,
                         int *IPIV,
                         cham_bool_t check_info, int iinfo );
void INSERT_TASK_zunmlq( const RUNTIME_option_t *options,
                         cham_side_t side, cham_trans_t trans,
                         int m, int n, int ib,  int nb, int k,
                         const CHAM_desc_t *A, int Am, int An,
                         const CHAM_desc_t *T, int Tm, int Tn,
                         const CHAM_desc_t *C, int Cm, int Cn );
void INSERT_TASK_zunmqr( const RUNTIME_option_t *options,
                         cham_side_t side, cham_trans_t trans,
                         int m, int n, int k, int ib, int nb,
                         const CHAM_desc_t *A, int Am, int An,
                         const CHAM_desc_t *T, int Tm, int Tn,
                         const CHAM_desc_t *C, int Cm, int Cn );

/**
 * Keep these insert_task for retro-compatibility
 */
static inline void
INSERT_TASK_ztslqt( const RUNTIME_option_t *options,
                    int m, int n, int ib, int nb,
                    const CHAM_desc_t *A1, int A1m, int A1n,
                    const CHAM_desc_t *A2, int A2m, int A2n,
                    const CHAM_desc_t *T, int Tm, int Tn )
{
    INSERT_TASK_ztplqt( options, m, n, 0, ib, nb,
                        A1, A1m, A1n,
                        A2, A2m, A2n,
                        T,  Tm,  Tn );
}

static inline void
INSERT_TASK_ztsqrt( const RUNTIME_option_t *options,
                    int m, int n, int ib, int nb,
                    const CHAM_desc_t *A1, int A1m, int A1n,
                    const CHAM_desc_t *A2, int A2m, int A2n,
                    const CHAM_desc_t *T, int Tm, int Tn )
{
    INSERT_TASK_ztpqrt( options, m, n, 0, ib, nb,
                        A1, A1m, A1n,
                        A2, A2m, A2n,
                        T,  Tm,  Tn );
}

static inline void
INSERT_TASK_zttlqt( const RUNTIME_option_t *options,
                    int m, int n, int ib, int nb,
                    const CHAM_desc_t *A1, int A1m, int A1n,
                    const CHAM_desc_t *A2, int A2m, int A2n,
                    const CHAM_desc_t *T, int Tm, int Tn )
{
    INSERT_TASK_ztplqt( options, m, n, n, ib, nb,
                        A1, A1m, A1n,
                        A2, A2m, A2n,
                        T,  Tm,  Tn );
}

static inline void
INSERT_TASK_zttqrt( const RUNTIME_option_t *options,
                    int m, int n, int ib, int nb,
                    const CHAM_desc_t *A1, int A1m, int A1n,
                    const CHAM_desc_t *A2, int A2m, int A2n,
                    const CHAM_desc_t *T, int Tm, int Tn )
{
    INSERT_TASK_ztpqrt( options, m, n, m, ib, nb,
                        A1, A1m, A1n,
                        A2, A2m, A2n,
                        T,  Tm,  Tn );
}

static inline void
INSERT_TASK_ztsmlq( const RUNTIME_option_t *options,
                    cham_side_t side, cham_trans_t trans,
                    int m1, int n1, int m2, int n2, int k, int ib, int nb,
                    const CHAM_desc_t *A1, int A1m, int A1n,
                    const CHAM_desc_t *A2, int A2m, int A2n,
                    const CHAM_desc_t *V, int Vm, int Vn,
                    const CHAM_desc_t *T, int Tm, int Tn )
{
    (void)m1;
    (void)n1;
    INSERT_TASK_ztpmlqt( options, side, trans, m2, n2, k, 0, ib, nb,
                         V, Vm, Vn, T, Tm, Tn,
                         A1, A1m, A1n, A2, A2m, A2n );
}

static inline void
INSERT_TASK_ztsmqr( const RUNTIME_option_t *options,
                    cham_side_t side, cham_trans_t trans,
                    int m1, int n1, int m2, int n2, int k, int ib, int nb,
                    const CHAM_desc_t *A1, int A1m, int A1n,
                    const CHAM_desc_t *A2, int A2m, int A2n,
                    const CHAM_desc_t *V, int Vm, int Vn,
                    const CHAM_desc_t *T, int Tm, int Tn )
{
    (void)m1;
    (void)n1;
    INSERT_TASK_ztpmqrt( options, side, trans, m2, n2, k, 0, ib, nb,
                         V, Vm, Vn, T, Tm, Tn,
                         A1, A1m, A1n, A2, A2m, A2n );
}

static inline void
INSERT_TASK_zttmlq( const RUNTIME_option_t *options,
                    cham_side_t side, cham_trans_t trans,
                    int m1, int n1, int m2, int n2, int k, int ib, int nb,
                    const CHAM_desc_t *A1, int A1m, int A1n,
                    const CHAM_desc_t *A2, int A2m, int A2n,
                    const CHAM_desc_t *V, int Vm, int Vn,
                    const CHAM_desc_t *T, int Tm, int Tn )
{
    (void)m1;
    (void)n1;
    INSERT_TASK_ztpmlqt( options, side, trans, m2, n2, k, n2, ib, nb,
                         V, Vm, Vn, T, Tm, Tn,
                         A1, A1m, A1n, A2, A2m, A2n );
}

static inline void
INSERT_TASK_zttmqr( const RUNTIME_option_t *options,
                    cham_side_t side, cham_trans_t trans,
                    int m1, int n1, int m2, int n2, int k, int ib, int nb,
                    const CHAM_desc_t *A1, int A1m, int A1n,
                    const CHAM_desc_t *A2, int A2m, int A2n,
                    const CHAM_desc_t *V, int Vm, int Vn,
                    const CHAM_desc_t *T, int Tm, int Tn )
{
    (void)m1;
    (void)n1;
    INSERT_TASK_ztpmqrt( options, side, trans, m2, n2, k, m2, ib, nb,
                         V, Vm, Vn, T, Tm, Tn,
                         A1, A1m, A1n, A2, A2m, A2n );
}

/**
 * Gram prototype
 */
void INSERT_TASK_zgram( const RUNTIME_option_t *options,
                        cham_uplo_t uplo,
                        int m, int n, int mt, int nt,
                        const CHAM_desc_t *Di, int Dim, int Din,
                        const CHAM_desc_t *Dj, int Djm, int Djn,
                        const CHAM_desc_t *D, int Dm, int Dn,
                        CHAM_desc_t *A, int Am, int An);

#endif /* _chameleon_tasks_z_h_ */
