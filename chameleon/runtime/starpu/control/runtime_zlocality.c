/**
 *
 * @file starpu/runtime_zlocality.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2021 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon StarPU CHAMELEON_Complex64_t kernel locality management
 *
 * @version 1.0.0
 * @author Cedric Augonnet
 * @author Mathieu Faverge
 * @author Cedric Castagnede
 * @author Florent Pruvost
 * @date 2020-03-03
 * @precisions normal z -> s d c
 *
 */
#include "chameleon_starpu.h"
#include "runtime_codelet_z.h"

#ifdef CHAMELEON_USE_CUDA
/* Only codelets with multiple choices are present here */
void RUNTIME_zlocality_allrestrict( uint32_t where )
{

    /* Blas 3 */
    cl_zgemm_restrict_where( where );
#if defined(PRECISION_z) || defined(PRECISION_c)
    cl_zhemm_restrict_where( where );
    cl_zher2k_restrict_where( where );
    cl_zherk_restrict_where( where );
    cl_zsytrf_nopiv_restrict_where( where );
#endif
    cl_zsymm_restrict_where( where );
    cl_zsyr2k_restrict_where( where );
    cl_zsyrk_restrict_where( where );
    cl_ztrmm_restrict_where( where );
    cl_ztrsm_restrict_where( where );

    /*
     * Lapack
     */
    /* Cholesky */
    cl_zpotrf_restrict_where( where );
    cl_zlauum_restrict_where( where );
    cl_ztrtri_restrict_where( where );

    /* LU */
    cl_zgetrf_incpiv_restrict_where( where );
    cl_zgetrf_nopiv_restrict_where( where );
    cl_zgessm_restrict_where( where );
    cl_zssssm_restrict_where( where );
    cl_ztstrf_restrict_where( where );

    /* QR */
    cl_zgeqrt_restrict_where( where );
    cl_ztpqrt_restrict_where( where );
    cl_zunmqr_restrict_where( where );
    cl_ztpmqrt_restrict_where( where );

    /* LQ */
   cl_zgelqt_restrict_where( where );
   cl_ztplqt_restrict_where( where );
   cl_zunmlq_restrict_where( where );
   cl_ztpmlqt_restrict_where( where );

}

void RUNTIME_zlocality_onerestrict( cham_tasktype_t kernel, uint32_t where )
{
    switch( kernel ) {
    /* Blas 3 */
    case TASK_GEMM:   cl_zgemm_restrict_where( where );  break;
#if defined(PRECISION_z) || defined(PRECISION_c)
    case TASK_HEMM:   cl_zhemm_restrict_where( where );  break;
    case TASK_HER2K:  cl_zher2k_restrict_where( where ); break;
    case TASK_HERK:   cl_zherk_restrict_where( where );  break;
    case TASK_SYTRF_NOPIV: cl_zsytrf_nopiv_restrict_where( where );  break;
#endif
    case TASK_SYMM:   cl_zhemm_restrict_where( where );  break;
    case TASK_SYR2K:  cl_zher2k_restrict_where( where ); break;
    case TASK_SYRK:   cl_zherk_restrict_where( where );  break;
    case TASK_TRMM:   cl_ztrmm_restrict_where( where );  break;
    case TASK_TRSM:   cl_ztrsm_restrict_where( where );  break;

    /*
     * Lapack
     */
    /* Cholesky */
    case TASK_POTRF:  cl_zpotrf_restrict_where( where ); break;
    case TASK_LAUUM:  cl_zlauum_restrict_where( where ); break;
    case TASK_TRTRI:  cl_ztrtri_restrict_where( where ); break;

    /* LU */
    case TASK_GETRF_INCPIV: cl_zgetrf_incpiv_restrict_where( where ); break;
    case TASK_GETRF_NOPIV: cl_zgetrf_nopiv_restrict_where( where ); break;
    case TASK_GESSM:  cl_zgessm_restrict_where( where ); break;
    case TASK_SSSSM:  cl_zssssm_restrict_where( where ); break;
    case TASK_TSTRF:  cl_ztstrf_restrict_where( where ); break;

    /* QR */
    case TASK_GEQRT:  cl_zgeqrt_restrict_where( where ); break;
    case TASK_TPQRT:  cl_ztpqrt_restrict_where( where ); break;
    case TASK_UNMQR:  cl_zunmqr_restrict_where( where ); break;
    case TASK_TPMQRT: cl_ztpmqrt_restrict_where( where ); break;

    /* LQ */
   case TASK_GELQT:  cl_zgelqt_restrict_where( where ); break;
   case TASK_TPLQT:  cl_ztplqt_restrict_where( where ); break;
   case TASK_UNMLQ:  cl_zunmlq_restrict_where( where ); break;
   case TASK_TPMLQT: cl_ztpmlqt_restrict_where( where ); break;

    default:
      return;
    }
}

void RUNTIME_zlocality_allrestore( )
{
    /* Blas 3 */
    cl_zgemm_restore_where();
#if defined(PRECISION_z) || defined(PRECISION_c)
    cl_zhemm_restore_where();
    cl_zher2k_restore_where();
    cl_zherk_restore_where();
    cl_zsytrf_nopiv_restore_where();
#endif
    cl_zsymm_restore_where();
    cl_zsyr2k_restore_where();
    cl_zsyrk_restore_where();
    cl_ztrmm_restore_where();
    cl_ztrsm_restore_where();

    /*
     * Lapack
     */
    /* Cholesky */
    cl_zpotrf_restore_where();
    cl_zlauum_restore_where();
    cl_ztrtri_restore_where();

    /* LU incpiv */
    cl_zgetrf_incpiv_restore_where();
    cl_zgessm_restore_where();
    cl_zssssm_restore_where();
    cl_ztstrf_restore_where();

    /* QR */
    cl_zgeqrt_restore_where();
    cl_ztpqrt_restore_where();
    cl_zunmqr_restore_where();
    cl_ztpmqrt_restore_where();

    /* LQ */
   cl_zgelqt_restore_where();
   cl_ztplqt_restore_where();
   cl_zunmlq_restore_where();
   cl_ztpmlqt_restore_where();

}

void RUNTIME_zlocality_onerestore( cham_tasktype_t kernel )
{
    switch( kernel ) {
    /* Blas 3 */
    case TASK_GEMM:   cl_zgemm_restore_where();  break;
#if defined(PRECISION_z) || defined(PRECISION_c)
    case TASK_HEMM:   cl_zhemm_restore_where();  break;
    case TASK_HER2K:  cl_zher2k_restore_where(); break;
    case TASK_HERK:   cl_zherk_restore_where();  break;
    case TASK_SYTRF_NOPIV: cl_zsytrf_nopiv_restore_where();  break;
#endif
    case TASK_SYMM:   cl_zhemm_restore_where();  break;
    case TASK_SYR2K:  cl_zher2k_restore_where(); break;
    case TASK_SYRK:   cl_zherk_restore_where();  break;
    case TASK_TRMM:   cl_ztrmm_restore_where();  break;
    case TASK_TRSM:   cl_ztrsm_restore_where();  break;

    /*
     * Lapack
     */
    /* Cholesky */
    case TASK_POTRF:  cl_zpotrf_restore_where(); break;
    case TASK_LAUUM:  cl_zlauum_restore_where(); break;
    case TASK_TRTRI:  cl_ztrtri_restore_where(); break;

    /* LU */
    case TASK_GETRF_INCPIV: cl_zgetrf_incpiv_restore_where(); break;
    case TASK_GETRF_NOPIV: cl_zgetrf_nopiv_restore_where(); break;
    case TASK_GESSM:  cl_zgessm_restore_where(); break;
    case TASK_SSSSM:  cl_zssssm_restore_where(); break;
    case TASK_TSTRF:  cl_ztstrf_restore_where(); break;

    /* QR */
    case TASK_GEQRT:  cl_zgeqrt_restore_where(); break;
    case TASK_TPQRT:  cl_ztpqrt_restore_where(); break;
    case TASK_UNMQR:  cl_zunmqr_restore_where(); break;
    case TASK_TPMQRT: cl_ztpmqrt_restore_where(); break;

    /* LQ */
   case TASK_GELQT:  cl_zgelqt_restore_where(); break;
   case TASK_TPLQT:  cl_ztplqt_restore_where(); break;
   case TASK_UNMLQ:  cl_zunmlq_restore_where(); break;
   case TASK_TPMLQT: cl_ztpmlqt_restore_where(); break;

    default:
      return;
    }
}
#else

void RUNTIME_zlocality_allrestrict( uint32_t where )
{
    (void)where;
}

void RUNTIME_zlocality_onerestrict( cham_tasktype_t kernel, uint32_t where )
{
    (void)kernel;
    (void)where;
}

void RUNTIME_zlocality_allrestore( ) {}

void RUNTIME_zlocality_onerestore( cham_tasktype_t kernel )
{
    (void)kernel;
}

#endif
