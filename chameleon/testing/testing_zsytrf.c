/**
 *
 * @file testing_zsytrf.c
 *
 * @copyright 2019-2021 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zsytrf testing
 *
 * @version 1.1.0
 * @author Lucas Barros de Assis
 * @author Mathieu Faverge
 * @date 2020-11-19
 * @precisions normal z -> c
 *
 */
#include <chameleon.h>
#include "testings.h"
#include "testing_zcheck.h"
#include <chameleon/flops.h>

int
testing_zsytrf( run_arg_list_t *args, int check )
{
    int          hres   = 0;
    CHAM_desc_t *descA;

    /* Reads arguments */
    intptr_t    mtxfmt = parameters_getvalue_int( "mtxfmt" );
    int         nb     = run_arg_get_int( args, "nb", 320 );
    int         P      = parameters_getvalue_int( "P" );
    cham_uplo_t uplo   = run_arg_get_uplo( args, "uplo", ChamUpper );
    int         N      = run_arg_get_int( args, "N", 1000 );
    int         LDA    = run_arg_get_int( args, "LDA", N );
    int         seedA  = run_arg_get_int( args, "seedA", random() );
    int         Q      = parameters_compute_q( P );
    cham_fixdbl_t t, gflops;
    cham_fixdbl_t flops = flops_zpotrf( N );

    CHAMELEON_Set( CHAMELEON_TILE_SIZE, nb );

    /* Creates the matrices */
    CHAMELEON_Desc_Create(
        &descA, (void*)(-mtxfmt), ChamComplexDouble, nb, nb, nb * nb, LDA, N, 0, 0, N, N, P, Q );

    /* Fills the matrix with random values */
    CHAMELEON_zplgsy_Tile( (double)N, uplo, descA, seedA );

    /* Calculates the solution */
    START_TIMING( t );
    hres = CHAMELEON_zsytrf_Tile( uplo, descA );
    STOP_TIMING( t );
    gflops = flops * 1.e-9 / t;
    run_arg_add_fixdbl( args, "time", t );
    run_arg_add_fixdbl( args, "gflops", ( hres == CHAMELEON_SUCCESS ) ? gflops : -1. );

    /* Checks the factorisation and residue */
    if ( check ) {
        CHAM_desc_t *descA0 = CHAMELEON_Desc_Copy( descA, NULL );
        CHAMELEON_zplgsy_Tile( (double)N, uplo, descA0, seedA );

        hres += check_zxxtrf( args, ChamSymmetric, uplo, descA0, descA );

        CHAMELEON_Desc_Destroy( &descA0 );
    }

    CHAMELEON_Desc_Destroy( &descA );

    return hres;
}

testing_t   test_zsytrf;
const char *zsytrf_params[] = { "mtxfmt", "nb","uplo", "n", "lda", "seedA", NULL };
const char *zsytrf_output[] = { NULL };
const char *zsytrf_outchk[] = { "RETURN", NULL };

/**
 * @brief Testing registration function
 */
void testing_zsytrf_init( void ) __attribute__( ( constructor ) );
void
testing_zsytrf_init( void )
{
    test_zsytrf.name        = "zsytrf";
    test_zsytrf.helper      = "Symmetric trinagular factorization";
    test_zsytrf.params      = zsytrf_params;
    test_zsytrf.output      = zsytrf_output;
    test_zsytrf.outchk      = zsytrf_outchk;
    test_zsytrf.fptr        = testing_zsytrf;
    test_zsytrf.next        = NULL;

    testing_register( &test_zsytrf );
}
