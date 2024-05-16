/**
 *
 * @file testing_ztrtri.c
 *
 * @copyright 2019-2021 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon ztrtri testing
 *
 * @version 1.1.0
 * @author Lucas Barros de Assis
 * @author Mathieu Faverge
 * @date 2020-11-19
 * @precisions normal z -> c d s
 *
 */
#include <chameleon.h>
#include "testings.h"
#include "testing_zcheck.h"
#include <chameleon/flops.h>

int
testing_ztrtri( run_arg_list_t *args, int check )
{
    int          hres   = 0;
    CHAM_desc_t *descA;

    /* Reads arguments */
    intptr_t    mtxfmt = parameters_getvalue_int( "mtxfmt" );
    int         nb     = run_arg_get_int( args, "nb", 320 );
    int         P      = parameters_getvalue_int( "P" );
    cham_uplo_t uplo   = run_arg_get_uplo( args, "uplo", ChamUpper );
    cham_diag_t diag   = run_arg_get_diag( args, "diag", ChamNonUnit );
    int         N      = run_arg_get_int( args, "N", 1000 );
    int         LDA    = run_arg_get_int( args, "LDA", N );
    int         seedA  = run_arg_get_int( args, "seedA", random() );
    int         Q      = parameters_compute_q( P );
    cham_fixdbl_t t, gflops;
    cham_fixdbl_t flops = flops_ztrtri( N );

    CHAMELEON_Set( CHAMELEON_TILE_SIZE, nb );

    /* Creates the matrices */
    CHAMELEON_Desc_Create(
        &descA, (void*)(-mtxfmt), ChamComplexDouble, nb, nb, nb * nb, LDA, N, 0, 0, N, N, P, Q );

    /* Initialises the matrices with the same values */
    CHAMELEON_zplghe_Tile( (double)N, uplo, descA, seedA );

    /* Calculates the inversed matrices */
    START_TIMING( t );
    hres = CHAMELEON_ztrtri_Tile( uplo, diag, descA );
    STOP_TIMING( t );
    gflops = flops * 1.e-9 / t;
    run_arg_add_fixdbl( args, "time", t );
    run_arg_add_fixdbl( args, "gflops", ( hres == CHAMELEON_SUCCESS ) ? gflops : -1. );

    /* Checks the inverse */
    if ( check ) {
        CHAM_desc_t *descA0 = CHAMELEON_Desc_Copy( descA, NULL );
        CHAMELEON_zplghe_Tile( (double)N, uplo, descA0, seedA );

        hres += check_ztrtri( args, ChamTriangular, uplo, diag, descA0, descA );

        CHAMELEON_Desc_Destroy( &descA0 );
    }

    CHAMELEON_Desc_Destroy( &descA );

    return hres;
}

testing_t   test_ztrtri;
const char *ztrtri_params[] = { "mtxfmt", "nb","uplo", "diag", "n", "lda", "seedA", NULL };
const char *ztrtri_output[] = { NULL };
const char *ztrtri_outchk[] = { "RETURN", NULL };

/**
 * @brief Testing registration function
 */
void testing_ztrtri_init( void ) __attribute__( ( constructor ) );
void
testing_ztrtri_init( void )
{
    test_ztrtri.name        = "ztrtri";
    test_ztrtri.helper      = "Triangular matrix inversion";
    test_ztrtri.params      = ztrtri_params;
    test_ztrtri.output      = ztrtri_output;
    test_ztrtri.outchk      = ztrtri_outchk;
    test_ztrtri.fptr        = testing_ztrtri;
    test_ztrtri.next        = NULL;

    testing_register( &test_ztrtri );
}
