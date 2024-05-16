/**
 *
 * @file coreblas_zc.h
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2021 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon CPU complex mixed precision kernels header
 *
 * @version 1.0.0
 * @comment This file has been automatically generated
 *          from Plasma 2.5.0 for CHAMELEON 0.9.2
 * @author Jakub Kurzak
 * @author Hatem Ltaief
 * @author Mathieu Faverge
 * @author Emmanuel Agullo
 * @author Cedric Castagnede
 * @author Florent Pruvost
 * @date 2020-03-03
 * @precisions mixed zc -> ds
 *
 */
#ifndef _coreblas_zc_h_
#define _coreblas_zc_h_

/**
 *  Declarations of serial kernels - alphabetical order
 */
void CORE_clag2z(int m, int n,
                 const CHAMELEON_Complex32_t *A, int lda,
                 CHAMELEON_Complex64_t *B, int ldb);
void CORE_zlag2c(int m, int n,
                 const CHAMELEON_Complex64_t *A, int lda,
                 CHAMELEON_Complex32_t *B, int ldb, int *info);

#endif /* _coreblas_zc_h_ */
