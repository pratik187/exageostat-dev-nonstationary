/**
 *
 * @file openmp/runtime_zlocality.c
 *
 * @copyright 2012-2017 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2021 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon OpenMP CHAMELEON_Complex64_t kernel locality management
 *
 * @version 1.1.0
 * @author Mathieu Faverge
 * @date 2020-10-15
 * @precisions normal z -> s d c
 *
 */
#include "chameleon_openmp.h"

void RUNTIME_zlocality_allrestrict( uint32_t where )
{
    (void)where;
    chameleon_warning("RUNTIME_zlocality_allrestrict(Openmp)", "Kernel locality cannot be specified with Openmp");
}

void RUNTIME_zlocality_onerestrict( cham_tasktype_t kernel, uint32_t where )
{
    (void)kernel;
    (void)where;
    chameleon_warning("RUNTIME_zlocality_onerestrict(Openmp)", "Kernel locality cannot be specified with Openmp");
}

void RUNTIME_zlocality_allrestore( )
{
    chameleon_warning("RUNTIME_zlocality_allrestore(Openmp)", "Kernel locality cannot be specified with Openmp");
}

void RUNTIME_zlocality_onerestore( cham_tasktype_t kernel )
{
    (void)kernel;
    chameleon_warning("RUNTIME_zlocality_onerestore(Openmp)", "Kernel locality cannot be specified with Openmp");
}
