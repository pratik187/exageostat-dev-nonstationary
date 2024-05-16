/**
 *
 * @file testing_zgenm2.c
 *
 * @copyright 2019-2021 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zgenm2 testing
 *
 * @version 1.1.0
 * @author Lucas Barros de Assis
 * @author Mathieu Faverge
 * @date 2020-12-01
 * @precisions normal z -> c d s
 *
 */
#include <chameleon.h>
#if !defined(CHAMELEON_SIMULATION)
#include <coreblas/lapacke.h>
#include <coreblas/cblas.h>
#include <coreblas.h>
#endif
#include "testings.h"
#include "testing_zcheck.h"
#include <chameleon/flops.h>

static cham_fixdbl_t
flops_zgenm2( int M, int N )
{
    double coefabs = 1.;
#if defined( PRECISION_z ) || defined( PRECISION_c )
    coefabs = 3.;
#endif

    return coefabs * M * N;
}

int
testing_zgenm2( run_arg_list_t *args, int check )
{
    int          hres   = 0;
    double       norm;
    CHAM_desc_t *descA;
    double      *D, dmax = 1.;

    /* Reads arguments */
    intptr_t        mtxfmt = parameters_getvalue_int( "mtxfmt" );
    int             nb     = run_arg_get_int( args, "nb", 320 );
    int             P      = parameters_getvalue_int( "P" );
    int             N      = run_arg_get_int( args, "N", 1000 );
    int             M      = run_arg_get_int( args, "M", N );
    int             LDA    = run_arg_get_int( args, "LDA", M );
    int             seedA  = run_arg_get_int( args, "seedA", random() );
    int             Q      = parameters_compute_q( P );
    int             minMN  = chameleon_min( M, N );
    double          cond   = run_arg_get_double( args, "cond", 1.e16 );
    int             mode   = run_arg_get_int( args, "mode", 4 );
    double          tol    = 1.e-1;
    cham_fixdbl_t t, gflops;
    cham_fixdbl_t flops = flops_zgenm2( M, N );

    CHAMELEON_Set( CHAMELEON_TILE_SIZE, nb );

    /* Generate the diagonal of eigen/singular values */
    D = malloc( minMN * sizeof(double) );
#if !defined(CHAMELEON_SIMULATION)
    hres = CORE_dlatm1( mode, cond, 0, ChamDistUniform, seedA, D, minMN );
    if ( hres != 0 ) {
        free( D );
        return hres;
    }

    /* Save the largest absolute value */
    hres = cblas_idamax( minMN, D, 1 );
    dmax = fabs( D[hres] );
#else
    (void)mode;
#endif

    /* Creates the matrix */
    CHAMELEON_Desc_Create(
        &descA, (void*)(-mtxfmt), ChamComplexDouble, nb, nb, nb * nb, LDA, N, 0, 0, M, N, P, Q );

    /* Fills the matrix with random values */
    hres = CHAMELEON_zlatms_Tile(
        ChamDistUniform, seedA, ChamNonsymPosv, D, 0, cond, 0., descA );
    free( D );
    if ( hres != 0 ) {
        return hres;
    }

    /* Calculates the norm */
    START_TIMING( t );
    norm = CHAMELEON_zgenm2_Tile( tol, descA );
    STOP_TIMING( t );
    gflops = flops * 1.e-9 / t;
    run_arg_add_fixdbl( args, "time", t );
    run_arg_add_fixdbl( args, "gflops", ( norm >= 0. ) ? gflops : -1. );

    /* Checks the solution */
    hres = 0;
    if ( check ) {
        double res = fabs(dmax - norm) / (dmax * tol);

        run_arg_add_double( args, "||A||", dmax );
        run_arg_add_double( args, "||B||", norm );
        run_arg_add_double( args, "||R||", res );

        if ( isnan(res) || isinf(res) || (res > 10.0) ) {
            hres = 1;
        }
    }

    CHAMELEON_Desc_Destroy( &descA );

    return hres;
}

testing_t   test_zgenm2;
const char *zgenm2_params[] = { "mtxfmt", "nb", "m", "n", "lda", "seedA", "cond", "mode", NULL };
const char *zgenm2_output[] = { NULL };
const char *zgenm2_outchk[] = { "||A||", "||B||", "||R||", "RETURN", NULL };

/**
 * @brief Testing registration function
 */
void testing_zgenm2_init( void ) __attribute__( ( constructor ) );
void
testing_zgenm2_init( void )
{
    test_zgenm2.name        = "zgenm2";
    test_zgenm2.helper      = "General matrix two-norm estimator";
    test_zgenm2.params      = zgenm2_params;
    test_zgenm2.output      = zgenm2_output;
    test_zgenm2.outchk      = zgenm2_outchk;
    test_zgenm2.fptr        = testing_zgenm2;
    test_zgenm2.next        = NULL;

    testing_register( &test_zgenm2 );
}
