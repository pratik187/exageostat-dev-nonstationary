/**
 *
 * @file starpu/codelet_zcallback.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2021 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zcallback StarPU codelet
 *
 * @version 1.1.0
 *  @author Mathieu Faverge
 *  @author Cedric Augonnet
 *  @author Florent Pruvost
 *  @date 2021-03-16
 *  @precisions normal z -> c d s
 *
 */
#include "chameleon_starpu.h"
#include "runtime_codelet_z.h"

#if defined(PRECISION_z) || defined(PRECISION_c)
CHAMELEON_CL_CB(dlag2z,        cti_handle_get_m(task->handles[1]), cti_handle_get_n(task->handles[1]), 0,                                      M*N)
#endif
CHAMELEON_CL_CB(dzasum,        cti_handle_get_m(task->handles[0]), cti_handle_get_n(task->handles[0]), 0,                                      M*N)
CHAMELEON_CL_CB(zaxpy,         cti_handle_get_m(task->handles[0]), cti_handle_get_m(task->handles[1]), 0,                                      M)
CHAMELEON_CL_CB(zgeadd,        cti_handle_get_m(task->handles[0]), cti_handle_get_n(task->handles[0]), 0,                                      M*N)
CHAMELEON_CL_CB(ztradd,        cti_handle_get_m(task->handles[0]), cti_handle_get_n(task->handles[0]), 0,                                      M*N/2.)
CHAMELEON_CL_CB(zlascal,       cti_handle_get_m(task->handles[0]), cti_handle_get_n(task->handles[0]), 0,                                      M*N)
CHAMELEON_CL_CB(zgelqt,        cti_handle_get_m(task->handles[0]), cti_handle_get_n(task->handles[0]), 0,                                      (4./3.)*M*N*K)
CHAMELEON_CL_CB(zgemv,         cti_handle_get_m(task->handles[0]), cti_handle_get_n(task->handles[0]), 0,                                      2. *M*N  )
CHAMELEON_CL_CB(zgemm,         cti_handle_get_m(task->handles[2]), cti_handle_get_n(task->handles[2]), cti_handle_get_n(task->handles[0]),     2. *M*N*K) /* If A^t, computation is wrong */
CHAMELEON_CL_CB(zgeqrt,        cti_handle_get_m(task->handles[0]), cti_handle_get_n(task->handles[0]), 0,                                      (4./3.)*M*M*N)
CHAMELEON_CL_CB(zgessm,        cti_handle_get_m(task->handles[2]), cti_handle_get_m(task->handles[2]), cti_handle_get_m(task->handles[2]),     2. *M*N*K)
CHAMELEON_CL_CB(zgessq,        cti_handle_get_m(task->handles[0]), cti_handle_get_m(task->handles[0]), 0,                                      4.*M*N)
CHAMELEON_CL_CB(zgetrf,        cti_handle_get_m(task->handles[0]), cti_handle_get_m(task->handles[0]), cti_handle_get_m(task->handles[0]), (2./3.)*M*N*K)
CHAMELEON_CL_CB(zgetrf_incpiv, cti_handle_get_m(task->handles[0]), cti_handle_get_m(task->handles[0]), cti_handle_get_m(task->handles[0]), (2./3.)*M*N*K)
CHAMELEON_CL_CB(zgetrf_nopiv,  cti_handle_get_m(task->handles[0]), cti_handle_get_m(task->handles[0]), cti_handle_get_m(task->handles[0]), (2./3.)*M*N*K)
CHAMELEON_CL_CB(zgram,         cti_handle_get_m(task->handles[3]), cti_handle_get_n(task->handles[3]), 0,                                                M*N)
CHAMELEON_CL_CB(zhe2ge,        cti_handle_get_m(task->handles[0]), cti_handle_get_n(task->handles[0]), 0,                                       (1./2.0)*M*N)
CHAMELEON_CL_CB(zherfb,        cti_handle_get_m(task->handles[0]), cti_handle_get_n(task->handles[0]), 0,                                         2. *M* M*M)
#if defined(PRECISION_z) || defined(PRECISION_c)
CHAMELEON_CL_CB(zhemm,         cti_handle_get_m(task->handles[2]), cti_handle_get_n(task->handles[2]), 0,                                          2.*M*M *N)
CHAMELEON_CL_CB(zher2k,        cti_handle_get_m(task->handles[0]), cti_handle_get_n(task->handles[0]), 0,                                     ( 1.+2.*M*N)*M)
CHAMELEON_CL_CB(zherk,         cti_handle_get_m(task->handles[0]), cti_handle_get_n(task->handles[0]), 0,                                     ( 1.+   M)*M*N)
#endif
CHAMELEON_CL_CB(zlacpy,        cti_handle_get_m(task->handles[0]), cti_handle_get_n(task->handles[0]), 0,                                                M*N)
CHAMELEON_CL_CB(zlange,        cti_handle_get_m(task->handles[0]), cti_handle_get_n(task->handles[0]), 0,                                                M*N)
CHAMELEON_CL_CB(zlaset,        cti_handle_get_m(task->handles[0]), cti_handle_get_n(task->handles[0]), 0,                                                M*N)
CHAMELEON_CL_CB(zlaset2,       cti_handle_get_m(task->handles[0]), cti_handle_get_n(task->handles[0]), 0,                                                M*N)
CHAMELEON_CL_CB(zlatro,        cti_handle_get_m(task->handles[0]), cti_handle_get_n(task->handles[0]), 0,                                                M*N)
CHAMELEON_CL_CB(zlauum,        cti_handle_get_m(task->handles[0]), cti_handle_get_m(task->handles[0]), cti_handle_get_m(task->handles[0]), (1./3.)*M* M*M)
#if defined(PRECISION_z) || defined(PRECISION_c)
CHAMELEON_CL_CB(zplghe,        cti_handle_get_m(task->handles[0]), cti_handle_get_n(task->handles[0]), 0,                                                M*N)
CHAMELEON_CL_CB(zsytrf_nopiv,  cti_handle_get_m(task->handles[0]), 0, 0,                                                                           (1./3.)*M* M*M)
#endif
CHAMELEON_CL_CB(zplgsy,        cti_handle_get_m(task->handles[0]), cti_handle_get_n(task->handles[0]), 0,                                                M*N)
CHAMELEON_CL_CB(zplrnt,        cti_handle_get_m(task->handles[0]), cti_handle_get_n(task->handles[0]), 0,                                                M*N)
CHAMELEON_CL_CB(zbuild,        cti_handle_get_m(task->handles[0]), cti_handle_get_n(task->handles[0]), 0,                                                M*N)
CHAMELEON_CL_CB(zplssq,        cti_handle_get_m(task->handles[0]), cti_handle_get_n(task->handles[0]), 0,                                                M*N)
CHAMELEON_CL_CB(zplssq2,       cti_handle_get_m(task->handles[0]), cti_handle_get_n(task->handles[0]), 0,                                                2*N)
CHAMELEON_CL_CB(zpotrf,        cti_handle_get_m(task->handles[0]), 0, 0,                                                                           (1./3.)*M* M*M)
CHAMELEON_CL_CB(zssssm,        cti_handle_get_m(task->handles[0]), cti_handle_get_m(task->handles[0]), cti_handle_get_m(task->handles[0]), M*M*(2.*M+cti_handle_get_m(task->handles[2])))
CHAMELEON_CL_CB(zsymm,         cti_handle_get_m(task->handles[2]), cti_handle_get_n(task->handles[2]), 0,                                           2.*M*M *N)
CHAMELEON_CL_CB(zsyr2k,        cti_handle_get_m(task->handles[0]), cti_handle_get_n(task->handles[0]), 0,                                      ( 1.+2.*M*N)*M)
CHAMELEON_CL_CB(zsyrk,         cti_handle_get_m(task->handles[0]), cti_handle_get_n(task->handles[0]), 0,                                      ( 1.+   M)*M*N)
CHAMELEON_CL_CB(ztplqt,        cti_handle_get_m(task->handles[1]), cti_handle_get_n(task->handles[1]), cti_handle_get_m(task->handles[0]),       2.*M*N*K)
CHAMELEON_CL_CB(ztpqrt,        cti_handle_get_m(task->handles[1]), cti_handle_get_n(task->handles[1]), cti_handle_get_m(task->handles[0]),       2.*M*N*K)
CHAMELEON_CL_CB(ztpmlqt,       cti_handle_get_m(task->handles[3]), cti_handle_get_n(task->handles[3]), cti_handle_get_m(task->handles[2]),       4.*M*N*K)
CHAMELEON_CL_CB(ztpmqrt,       cti_handle_get_m(task->handles[3]), cti_handle_get_n(task->handles[3]), cti_handle_get_m(task->handles[2]),       4.*M*N*K)
CHAMELEON_CL_CB(ztrasm,        cti_handle_get_m(task->handles[0]), cti_handle_get_n(task->handles[0]), 0,                                         0.5*M*(M+1))
CHAMELEON_CL_CB(ztrmm,         cti_handle_get_m(task->handles[1]), cti_handle_get_n(task->handles[1]), 0,                                               M*M*N)
CHAMELEON_CL_CB(ztrsm,         cti_handle_get_m(task->handles[1]), cti_handle_get_n(task->handles[1]), 0,                                               M*M*N)
CHAMELEON_CL_CB(ztrtri,        cti_handle_get_m(task->handles[0]), cti_handle_get_m(task->handles[0]), cti_handle_get_m(task->handles[0]), (1./3.)*M *M*M)
CHAMELEON_CL_CB(ztsmlq_hetra1, cti_handle_get_m(task->handles[0]), cti_handle_get_m(task->handles[0]), cti_handle_get_m(task->handles[0]), (4.0*M+cti_handle_get_m(task->handles[3]))*M*M)
CHAMELEON_CL_CB(ztsmqr_hetra1, cti_handle_get_m(task->handles[0]), cti_handle_get_m(task->handles[0]), cti_handle_get_m(task->handles[0]), (4.0*M+cti_handle_get_m(task->handles[3]))*M*M)
CHAMELEON_CL_CB(ztstrf,        cti_handle_get_m(task->handles[0]), cti_handle_get_m(task->handles[0]), cti_handle_get_m(task->handles[0]),         M* M*M)
CHAMELEON_CL_CB(zunmlq,        cti_handle_get_m(task->handles[0]), cti_handle_get_m(task->handles[0]), cti_handle_get_m(task->handles[0]),     2. *M* M*M)
CHAMELEON_CL_CB(zunmqr,        cti_handle_get_m(task->handles[0]), cti_handle_get_m(task->handles[0]), cti_handle_get_m(task->handles[0]),     2. *M* M*M)
