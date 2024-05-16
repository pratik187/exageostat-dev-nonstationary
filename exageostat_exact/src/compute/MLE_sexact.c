/**
 *
 * Copyright (c) 2017-2023, King Abdullah University of Science and Technology
 * All rights reserved.
 *
 * ExaGeoStat is a software package provided by KAUST
 **/
/**
 *
 * @file MLE_exact.c
 *
 * ExaGeoStat exact computation main functions (i.e., generate synthetic dataset, evaluate ML function, and predication).
 *
 * @version 1.2.0
 *
 * @author Sameh Abdulah
 * @date 2020-01-19
 *
 **/
#include "../include/MLE_exact_s.h"

//***************************************************************************************
void EXAGEOSTAT_MLE_szvg_Tile(MLE_data *data, float *Nrand, double*initial_theta, int n, int dts, int log)
//! Generate Observations Vector (Z) for testing Maximum
/*! Likelihood function -- CHAM-sync 
 * Returns Z observation vector
 * @param[in] data: MLE_data struct with different MLE inputs.
 * @param[in] Nrand: A uniform random vector with size n that is used to generate Z .
 * @param[in] initial_theta: Theta vector with three parameter (Variance, Range, Smoothness)
 * 	                     that is used to to generate the Covariance Matrix.
 * @param[in] n: Problem size (number spatial locations).
 * @param[in] dts: tile size (MB) is used only in the case of HiCMA not CHAM.
 * @param[in] test: if 0 -> real data mode, 1 ->test data mode.
 * @param[in] log: equals one if the user needs to generate log files for his problem.
 * */
{
    RUNTIME_sequence_t *msequence = (RUNTIME_sequence_t *) data->sequence;
    RUNTIME_request_t *mrequest = (RUNTIME_request_t *) data->request;

    //Generate the co-variance matrix C
    VERBOSE("Initializing Covariance Matrix (Synthetic Dataset Generation Phase - single precision) .....");
    //EXAGEOSTAT_MLE_scmg_Tile_Async(ChamLower, data->descC, msequence, mrequest, &data->l1, &data->l1, initial_theta, data->dm);
    EXAGEOSTAT_MLE_scmg_Tile_Async(ChamLower, data->descC, &data->l1, &data->l1, &data->lm, initial_theta, data->dm,
                                  data->kernel_fun, msequence, mrequest);
    VERBOSE(" Done.\n");

    //Copy Nrand to Z
    VERBOSE("Generate Normal Random Distribution Vector Z (Synthetic Dataset Generation Phase - single precision) .....");
    EXAGEOSTAT_MLE_szcpy_Tile_Async(data->descZ, Nrand, msequence, mrequest);
    VERBOSE(" Done.\n");

    //Cholesky factorization for the Co-variance matrix C
    VERBOSE("Cholesky factorization of Sigma (Synthetic Dataset Generation Phase - single precision) .....");
    int success = CHAMELEON_spotrf_Tile(ChamLower, data->descC);
    SUCCESS(success, "Factorization cannot be performed..\n The matrix is not positive definite\n\n");
    VERBOSE(" Done.\n");

    //Triangular matrix-matrix multiplication    
    VERBOSE("Triangular matrix-matrix multiplication Z=L.e (Synthetic Dataset Generation Phase - single precision) .....");
    CHAMELEON_strmm_Tile(ChamLeft, ChamLower, ChamNoTrans, ChamNonUnit, 1, data->descC, data->descZ);
    VERBOSE(" Done.\n");

    //if log==1 write vector to disk
    if (log == 1) {
        double*z;
        CHAM_desc_t *CHAM_descZ = (CHAM_desc_t *) (data->descZ);
#if defined(CHAMELEON_USE_MPI)
        z = (double*) malloc(n * sizeof(double));
        CHAMELEON_Tile_to_Lapack( CHAM_descZ, z, n);
#else
        z = CHAM_descZ->mat;
#endif
        write_vectors(z, data, n);
#if defined(CHAMELEON_USE_MPI)
        free(z);
#endif
    }

    CHAMELEON_slaset_Tile(ChamUpperLower, 0, 0, data->descC);
    VERBOSE("Done Z Vector Generation Phase. (Chameleon Synchronous - single precision)\n");
    VERBOSE("************************************************************\n");
}

void EXAGEOSTAT_MLE_szcpy(MLE_data *data, double*streamdata) {
    RUNTIME_sequence_t *msequence = (RUNTIME_sequence_t *) data->sequence;
    RUNTIME_request_t *mrequest = (RUNTIME_request_t *) data->request;
    VERBOSE("Copy Z from vector to decriptor.\n");
    EXAGEOSTAT_MLE_szcpy_Tile_Async(data->descZ, streamdata, msequence, mrequest);
    CHAMELEON_Sequence_Wait(msequence);
    VERBOSE("Done Z copying step.\n");
    VERBOSE("************************************************************\n");
}

void EXAGEOSTAT_MLE_szvg_Tile_Async(MLE_data *data, float *Nrand, double*initial_theta, int n, int dts, int log)
//! Generate Observations Vector (Z) for testing Maximum
/*! Likelihood function -- CHAM-Async
 * Returns Z observation vector
 * @param[in] data: MLE_data struct with different MLE inputs.
 * @param[in] Nrand: A uniform random vector with size n that is used to generate Z .
 * @param[in] initial_theta: Theta vector with three parameter (Variance, Range, Smoothness)
 *                           that is used to to generate the Covariance Matrix.
 * @param[in] n: Problem size (number spatial locations).
 * @param[in] dts: tile size (MB) is used only in the case of HiCMA not CHAM.
 * @param[in] test: if 0 -> real data mode, 1 ->test data mode.
 * @param[in] log: equals one if the user needs to generate log files for his problem.
 * */
{
    RUNTIME_sequence_t *msequence = (RUNTIME_sequence_t *) data->sequence;
    RUNTIME_request_t *mrequest = (RUNTIME_request_t *) data->request;

    //Generate the co-variance matrix C
    VERBOSE("Initializing Covariance Matrix (Synthetic Dataset Generation Phase - single precision).....");
    EXAGEOSTAT_MLE_scmg_Tile_Async(ChamLower, data->descC, &data->l1, &data->l1, &data->lm, initial_theta, data->dm,
                                  data->kernel_fun, msequence, mrequest);
    VERBOSE(" Done.\n");

    //Copy Nrand to Z
    VERBOSE("Generate Normal Random Distribution Vector Z (Synthetic Dataset Generation Phase - single precision) .....");
    EXAGEOSTAT_MLE_szcpy_Tile_Async(data->descZ, Nrand, msequence, mrequest);
    VERBOSE(" Done.\n");

    //Cholesky factorization for the Co-variance matrix C
    VERBOSE("Cholesky factorization of Sigma (Synthetic Dataset Generation Phase - single precision) .....");
    int success = CHAMELEON_spotrf_Tile_Async(ChamLower, data->descC, msequence, mrequest);
    SUCCESS(success, "Factorization cannot be performed..\n The matrix is not positive definite\n\n");
    VERBOSE(" Done.\n");

    //Triangular matrix-matrix multiplication
    VERBOSE("Triangular matrix-matrix multiplication Z=L.e (Synthetic Dataset Generation Phase - single precision) .....");
    CHAMELEON_strmm_Tile_Async(ChamLeft, ChamLower, ChamNoTrans, ChamNonUnit, 1, data->descC, data->descZ, msequence,
                               mrequest);
    VERBOSE(" Done.\n");

    //if log == 1 write vector to disk
    if (log == 1) {
        double*z;
        CHAM_desc_t *CHAM_descZ = (CHAM_desc_t *) (data->descZ);
#if defined(CHAMELEON_USE_MPI)
        z = (double*) malloc(n * sizeof(double));
        CHAMELEON_Tile_to_Lapack( CHAM_descZ, z, n);
#else
        z = CHAM_descZ->mat;
#endif
        write_vectors(z, data, n);

#if defined(CHAMELEON_USE_MPI)
        free(z);
#endif
    }

    VERBOSE("Done Z Vector Generation Phase. (Chameleon Asynchronous - single precision)\n");
    VERBOSE("************************************************************\n");
}

double EXAGEOSTAT_smle_Tile(unsigned n, const double*theta, double*grad, void *CHAM_data) {
    //! Maximum Likelihood Evaluation (MLE)
    /*!  -- CHAM-sync-single precision.
     * Returns the loglikelihhod value for the given theta.
     * @param[in] n: unsigned variable used by NLOPT library.
     * @param[in] theta: theta Vector with three parameter (Variance, Range, Smoothness)
     *                           that is used to to generate the Covariance Matrix.
     * @param[in] grad: double variable used by NLOPT library. 
     * @param[in] CHAM_data: MLE_data struct with different MLE inputs.
     * */
    //Initialization
    double loglik = 0.0, logdet = 0.0, time_facto = 0.0, time_solve = 0.0, logdet_calculate = 0.0, matrix_gen_time = 0.0, zcpy_time = 0.0;
    int N, NRHS, success;
    double flops = 0.0;

    MLE_data *data = ((MLE_data *) CHAM_data);
    data->det = 0;
    data->dotp = 0;

    CHAM_desc_t *CHAM_descC = (CHAM_desc_t *) data->descC;
    CHAM_desc_t *CHAM_descZ = (CHAM_desc_t *) data->descZ;
    CHAM_desc_t *CHAM_descZcpy = (CHAM_desc_t *) data->descZcpy;
    CHAM_desc_t *CHAM_descdet = (CHAM_desc_t *) data->descdet;
    CHAM_desc_t *CHAM_descproduct = (CHAM_desc_t *) data->descproduct;
    RUNTIME_sequence_t *msequence = (RUNTIME_sequence_t *) data->sequence;
    RUNTIME_request_t *mrequest = (RUNTIME_request_t *) data->request;

    N = CHAM_descC->m;
    NRHS = CHAM_descZ->n;
    START_TIMING(zcpy_time);
    if (data->iter_count == 0) {
        //Save a copy of descZ into descZcpy for restoring each iteration (Only for the first iteration)
        CHAMELEON_slacpy_Tile(ChamUpperLower, CHAM_descZ, CHAM_descZcpy);
    } else {
        VERBOSE("Re-store the original Z vector...");
        CHAMELEON_slacpy_Tile(ChamUpperLower, CHAM_descZcpy, CHAM_descZ);
        VERBOSE(" Done.\n");
    }
    STOP_TIMING(zcpy_time);

    //Generate new co-variance matrix C based on new theta	
    VERBOSE("Generate New Covariance Matrix (single precision)...");
    START_TIMING(matrix_gen_time);
    EXAGEOSTAT_MLE_scmg_Tile_Async(ChamLower, CHAM_descC, &data->l1, &data->l1, &data->lm, (double*) theta, data->dm,
                                  data->kernel_fun, msequence, &mrequest[0]);
    CHAMELEON_Sequence_Wait(msequence);
    STOP_TIMING(matrix_gen_time);
    VERBOSE(" Done.\n");

    //Calculate Cholesky Factorization (C=LL-1)
    VERBOSE("Cholesky factorization of Sigma (single precision)...");
    START_TIMING(time_facto);
    success = CHAMELEON_spotrf_Tile(ChamLower, CHAM_descC);
    STOP_TIMING(time_facto);
    SUCCESS(success, "Factorization cannot be performed..\n The matrix is not positive definite\n\n");
    flops = flops + FLOPS_SPOTRF(N);
    VERBOSE(" Done.\n");

    //Calculate log(|C|) --> log(square(|L|))
    VERBOSE("Calculating the log determinant (single precision) ...");
    START_TIMING(logdet_calculate);
    EXAGEOSTAT_MLE_smdet_Tile_Async(CHAM_descC, msequence, &mrequest[0], CHAM_descdet);
    CHAMELEON_Sequence_Wait(msequence);
    logdet = 2 * (float) data->det;
    STOP_TIMING(logdet_calculate);
    VERBOSE(" Done.\n");

    //Solving Linear System (L*X=Z)--->inv(L)*Z
    VERBOSE("Solving the linear system (single precision)...\n");
    START_TIMING(time_solve);
    CHAMELEON_strsm_Tile(ChamLeft, ChamLower, ChamNoTrans, ChamNonUnit, 1, CHAM_descC, CHAM_descZ);
    STOP_TIMING(time_solve);
    flops = flops + FLOPS_STRSM(ChamLeft, N, NRHS);
    VERBOSE(" Done.\n");


    //Calculate MLE likelihood
    VERBOSE("Calculating the MLE likelihood function (single precision) ...");

    CHAMELEON_sgemm_Tile(ChamTrans, ChamNoTrans, 1, CHAM_descZ, CHAM_descZ, 0, CHAM_descproduct);

    loglik = -0.5 * data->dotp - 0.5 * logdet - (double) (N / 2.0) * log(2.0 * PI);
    data->variance = theta[0];

    //Distribute the values in the case of MPI
#if defined(CHAMELEON_USE_MPI)
    MPI_Bcast(&loglik,1, MPI_DOUBLE, 0, MPI_COMM_WORLD );
    MPI_Bcast(theta,3, MPI_DOUBLE, 0, MPI_COMM_WORLD );
    if(CHAMELEON_My_Mpi_Rank() == 0)
    {
#endif
    //Print Iteration Summary
    fprintf(stderr, "\n------ddotproduct: %.17g ", data->sdotp);
    fprintf(stderr, "\n------logdet: %.17g ", logdet);
    //reformat
    printf(" %3d- Model Parameters (variance, range, smoothness): (%2.6f, %2.6f, %2.6f) ----> LogLi: %.17g\n",
           data->iter_count + 1, data->variance, theta[1], theta[2], loglik);

    if (data->log == 1)
        fprintf(data->pFileLog,
                " %3d- Model Parameters (variance, range, smoothness): (%2.6f, %2.6f, %2.6f) ----> LogLi: %.3g\n",
                data->iter_count + 1, data->variance, theta[1], theta[2], loglik);

    printf(" ---- Facto Time: %6.2f\n", time_facto);
    printf(" ---- logdet Time: %6.2f\n", logdet_calculate);
    printf(" ---- dtrsm Time: %6.2f\n", time_solve);
    printf(" ---- Matrix Generation Time: %6.2f\n", matrix_gen_time);
    printf(" ---- Total Time: %6.2f\n", matrix_gen_time + time_facto + logdet_calculate + time_solve);
    printf(" ---- Gflop/s: %6.2f\n", flops / 1e9 / (time_facto + time_solve));
#if defined(CHAMELEON_USE_MPI)
    }
#endif

    data->iter_count++;
    // for experiments
    data->avg_exec_time_per_iter +=/*matrix_gen_time*/+time_facto + logdet_calculate + time_solve;
    data->avg_flops_per_iter += flops / 1e9 / (time_facto + time_solve);
    data->final_loglik = loglik;

    //output
    results.final_loglik = loglik;
    return loglik;
}

double EXAGEOSTAT_smle_Tile_Async(unsigned n, const double*theta, double*grad, void *CHAM_data) {
    //! Maximum Likelihood Evaluation (MLE)
    /*!  -- CHAM-Async
     * Returns the loglikelihhod value for the given theta.
     * @param[in] n: unsigned variable used by NLOPT library.
     * @param[in] theta: theta Vector with three parameter (Variance, Range, Smoothness)
     *                           that is used to to generate the Covariance Matrix.
     * @param[in] grad: double variable used by NLOPT library.
     * @param[in] CHAM_data: MLE_data struct with different MLE inputs.
     * */
    //Initialization
    double loglik = 0.0, logdet = 0.0, time_facto = 0.0, time_solve = 0.0, logdet_calculate = 0.0, matrix_gen_time = 0.0, zcpy_time = 0.0, flops = 0.0;
    int N, NRHS, success;

    MLE_data *data = ((MLE_data *) CHAM_data);
    data->det = 0;
    data->dotp = 0;

    CHAM_desc_t *CHAM_descC = (CHAM_desc_t *) data->descC;
    CHAM_desc_t *CHAM_descZ = (CHAM_desc_t *) data->descZ;
    CHAM_desc_t *CHAM_descZcpy = (CHAM_desc_t *) data->descZcpy;
    CHAM_desc_t *CHAM_descdet = (CHAM_desc_t *) data->descdet;
    CHAM_desc_t *CHAM_descproduct = (CHAM_desc_t *) data->descproduct;
    RUNTIME_sequence_t *msequence = (RUNTIME_sequence_t *) data->sequence;
    RUNTIME_request_t *mrequest = (RUNTIME_request_t *) data->request;

    N = CHAM_descC->m;
    NRHS = CHAM_descZ->n;
    START_TIMING(zcpy_time);
    if (data->iter_count == 0)
        //Save a copy of descZ into descZcpy for restoring each iteration (Only for the first iteration)
        CHAMELEON_slacpy_Tile_Async(ChamUpperLower, CHAM_descZ, CHAM_descZcpy, msequence, mrequest);
    else {
        VERBOSE("re-store the original Z vector (single precision)...");
        CHAMELEON_slacpy_Tile_Async(ChamUpperLower, CHAM_descZcpy, CHAM_descZ, msequence, mrequest);
        VERBOSE(" Done.\n");
    }
    STOP_TIMING(zcpy_time);

    //Generate new co-variance matrix C based on new theta	
    VERBOSE("Generate New Covariance Matrix (single precision)...");
    START_TIMING(matrix_gen_time);
    //EXAGEOSTAT_MLE_scmg_Tile_Async(ChamLower, CHAM_descC, msequence, mrequest, &data->l1, &data->l1, (double*) theta, data->dm);
    EXAGEOSTAT_MLE_scmg_Tile_Async(ChamLower, CHAM_descC, &data->l1, &data->l1, &data->lm, (double*) theta, data->dm,
                                  data->kernel_fun, msequence, &mrequest[0]);
    CHAMELEON_Sequence_Wait(msequence);
    STOP_TIMING(matrix_gen_time);
    VERBOSE(" Done.\n");

    //Calculate Cholesky Factorization (C=LL-1)
    VERBOSE("Cholesky factorization of Sigma (single precision)...");
    START_TIMING(time_facto);
    success = CHAMELEON_spotrf_Tile_Async(ChamLower, CHAM_descC, msequence, mrequest);
    STOP_TIMING(time_facto);
    SUCCESS(success, "Factorization cannot be performed..\n The matrix is not positive definite\n\n");
    flops = flops + FLOPS_SPOTRF(N);
    VERBOSE(" Done.\n");

    //Calculate log(|C|) --> log(square(|L|))
    VERBOSE("Calculating the log determinant (single precision) ...");
    START_TIMING(logdet_calculate);
    EXAGEOSTAT_MLE_smdet_Tile_Async(CHAM_descC, msequence, &mrequest[0], CHAM_descdet);
    CHAMELEON_Sequence_Wait(msequence);
    logdet = 2 * data->det;
    STOP_TIMING(logdet_calculate);
    VERBOSE(" Done.\n");

    //Solving Linear System (L*X=Z)--->inv(L)*Z
    VERBOSE("Solving the linear system (single precision) ...\n");
    START_TIMING(time_solve);
    CHAMELEON_strsm_Tile_Async(ChamLeft, ChamLower, ChamNoTrans, ChamNonUnit, 1, CHAM_descC, CHAM_descZ, msequence,
                               mrequest);
    STOP_TIMING(time_solve);
    flops = flops + FLOPS_STRSM(ChamLeft, N, NRHS);
    VERBOSE(" Done.\n");

    //Claculate MLE likelihood
    VERBOSE("Calculating the MLE likelihood function (single precision) ...");
    void *ws;
    CHAMELEON_sgemm_Tile_Async(ChamTrans, ChamNoTrans, 1, CHAM_descZ, CHAM_descZ, 0, CHAM_descproduct, ws, msequence,
                               mrequest);
    loglik = -0.5 * data->dotp - 0.5 * logdet - (double) (N / 2.0) * log(2.0 * PI);
    VERBOSE(" Done.\n");

    //Distribute the values in the case of MPI
#if defined(CHAMELEON_USE_MPI)
    MPI_Bcast(&loglik,1, MPI_DOUBLE, 0, MPI_COMM_WORLD );
    MPI_Bcast(theta,3, MPI_DOUBLE, 0, MPI_COMM_WORLD );
    if(CHAMELEON_My_Mpi_Rank() == 0)
    {
#endif

    //reformat
    fprintf(stderr, " %3d- Model Parameters (variance, range, smoothness): (%2.6f, %2.6f, %2.6f) ----> LogLi: %2.6f\n",
            data->iter_count + 1, theta[0], theta[1], theta[2], loglik);

    if (data->log == 1)
        fprintf(data->pFileLog,
                " %3d- Model Parameters (variance, range, smoothness): (%2.6f, %2.6f, %2.6f) ----> LogLi: %2.6f\n",
                data->iter_count + 1, theta[0], theta[1], theta[2], loglik);

    fprintf(stderr, " ---- Facto Time: %6.2f\n", time_facto);
    fprintf(stderr, " ---- logdet Time: %6.2f\n", logdet_calculate);
    fprintf(stderr, " ---- dtrsm Time: %6.2f\n", time_solve);
    fprintf(stderr, " ---- Matrix Generation Time: %6.2f\n", matrix_gen_time);
    fprintf(stderr, " ---- Total Time: %6.2f\n", matrix_gen_time + time_facto + logdet_calculate + time_solve);
    fprintf(stderr, " ---- Gflop/s: %6.2f\n", flops / 1e9 / (time_facto + time_solve));
#if defined(CHAMELEON_USE_MPI)
    }
#endif

    data->iter_count++;
    // for experiments
    data->avg_exec_time_per_iter += matrix_gen_time + time_facto + logdet_calculate + time_solve;
    data->avg_flops_per_iter += flops / 1e9 / (time_facto + time_solve);
    data->final_loglik = loglik;

    return loglik;
}

void EXAGEOSTAT_smle_Predict_Allocate(MLE_data *CHAM_data, int nZmiss, int nZobs, int dts, int p_grid, int q_grid,
                                     int mse_flag) {

    CHAM_desc_t *CHAM_descZmiss = NULL;
    CHAM_desc_t *CHAM_descC12 = NULL;
    CHAM_desc_t *CHAM_descC22 = NULL;
    CHAM_desc_t *CHAM_descmse = NULL;
    CHAM_desc_t *CHAM_descZactual = NULL;
    CHAM_desc_t *CHAM_descZobs = NULL;
    MLE_data *data = (MLE_data *) CHAM_data;

    if (nZmiss <= 0) {
        fprintf(stderr, " Number of missing values should be positive value\n");
        return;
    }
    //Descriptors Creation
    EXAGEOSTAT_ALLOCATE_MATRIX_TILE(&CHAM_descZobs, NULL, ChamRealFloat, dts, dts, dts * dts, nZobs, 1, 0, 0, nZobs, 1,
                                    p_grid, q_grid);
    if (mse_flag == 1) {
        EXAGEOSTAT_ALLOCATE_MATRIX_TILE(&CHAM_descZactual, NULL, ChamRealFloat, dts, dts, dts * dts, nZmiss, 1, 0, 0,
                                        nZmiss, 1, p_grid, q_grid);

        EXAGEOSTAT_ALLOCATE_MATRIX_TILE(&CHAM_descmse, &data->mserror, ChamRealFloat, dts, dts, dts * dts, 1, 1, 0, 0,
                                        1, 1, p_grid, q_grid);
    }
    EXAGEOSTAT_ALLOCATE_MATRIX_TILE(&CHAM_descZmiss, NULL, ChamRealFloat, dts, dts, dts * dts, nZmiss, 1, 0, 0, nZmiss,
                                    1, p_grid, q_grid);
    EXAGEOSTAT_ALLOCATE_MATRIX_TILE(&CHAM_descC12, NULL, ChamRealFloat, dts, dts, dts * dts, nZmiss, nZobs, 0, 0,
                                    nZmiss, nZobs, p_grid, q_grid);
    EXAGEOSTAT_ALLOCATE_MATRIX_TILE(&CHAM_descC22, NULL, ChamRealFloat, dts, dts, dts * dts, nZobs, nZobs, 0, 0, nZobs,
                                    nZobs, p_grid, q_grid);

    //Initiate data descriptors
    data->descZmiss = CHAM_descZmiss;
    data->descC12 = CHAM_descC12;
    data->descC22 = CHAM_descC22;
    data->descmse = CHAM_descmse;
    data->descZactual = CHAM_descZactual;
    data->descZobs = CHAM_descZobs;
}


double
EXAGEOSTAT_smle_Predict_Tile(MLE_data *CHAM_data, double*theta, int nZmiss, int nZobs, double*Zobs, double*Zactual,
                            double*Zmiss, int n)
//! //Predict missing values base on a set of given values and covariance matrix
/*!  -- CHAM-sync
 * Returns the prediction Mean Square Error (MSE) as double
 * @param[in] CHAM_data: MLE_data struct with different MLE inputs.
 * @param[in] theta: theta Vector with three parameter (Variance, Range, Smoothness)
 *                           that is used to to generate the Covariance Matrix.
 * @param[in] nZmiss: number of missing values (unknown observations).
 * @param[in] nZobs: number of observed values (known observations).
 * @param[in] n: number of spatial locations.
 * */
{

    //initialization	
    double time_solve = 0.0;
    double mat_gen_time = 0.0;
    double time_gemm = 0.0;
    double time_mse = 0.0;
    double flops = 0.0;

    CHAM_desc_t *CHAM_descZmiss = NULL;
    CHAM_desc_t *CHAM_descC12 = NULL;
    CHAM_desc_t *CHAM_descC22 = NULL;
    CHAM_desc_t *CHAM_descmse = NULL;
    CHAM_desc_t *CHAM_descZactual = NULL;
    CHAM_desc_t *CHAM_descZobs = NULL;
    MLE_data *data = (MLE_data *) CHAM_data;
    RUNTIME_sequence_t *msequence = (RUNTIME_sequence_t *) data->sequence;
    RUNTIME_request_t *mrequest = (RUNTIME_request_t *) data->request;
    data->mserror = 0;

    if (nZmiss <= 0) {
        fprintf(stderr, " Number of missing values should be positive value\n");
        return -1;
    }

    //Initiate data descriptors
    CHAM_descZmiss = data->descZmiss;
    CHAM_descC12 = data->descC12;
    CHAM_descC22 = data->descC22;
    CHAM_descmse = data->descmse;
    CHAM_descZactual = data->descZactual;
    CHAM_descZobs = data->descZobs;

#if defined(CHAMELEON_USE_MPI)
    MPI_Bcast(&data->variance,1, MPI_DOUBLE, 0, MPI_COMM_WORLD );
#endif
    theta[0] = data->variance;

    printf("estimated parameters: %f - %f - %f\n", theta[0], theta[1], theta[2]);
    //Copy data to vectors 
    VERBOSE("Copy measurments vector to descZobs descriptor...");
    CHAMELEON_Lapack_to_Tile(Zobs, nZobs, CHAM_descZobs);
    VERBOSE(" Done.\n");

    if (Zactual != NULL) {
        VERBOSE("Copy actual measurments vector to descZactual descriptor...");
        //EXAGEOSTAT_MLE_szcpy_Tile_Async(CHAM_descZactual, Zactual, msequence, mrequest);
        CHAMELEON_Lapack_to_Tile(Zactual, nZmiss, CHAM_descZactual);
        VERBOSE(" Done.\n");
    }

    CHAMELEON_Sequence_Wait(msequence);

    START_TIMING(mat_gen_time);
    //Generate C22 covariance matrix
    VERBOSE("Generate C22 Covariance Matrix... (Prediction Stage - single precision)");
    EXAGEOSTAT_MLE_scmg_Tile_Async(ChamLower, CHAM_descC22, &data->lobs, &data->lobs, &data->lm, theta, data->dm,
                                  data->kernel_fun, msequence, mrequest);
    CHAMELEON_Sequence_Wait(msequence);
    VERBOSE(" Done.\n");

    //Generate C12 covariance matrix
    VERBOSE("Generate C12 Covariance Matrix... (Prediction Stage - single precision)");
    //EXAGEOSTAT_MLE_scmg_Tile_Async(ChamLower, CHAM_descC12, msequence, mrequest,  &data->lmiss, &data->lobs, theta, data->dm);
    EXAGEOSTAT_MLE_dcmg_Tile_Async(ChamLower, CHAM_descC12, &data->lmiss, &data->lobs, &data->lm, theta, &data->ns_params,data->dm,
                                  data->kernel_fun, data->Nsplit ,msequence, mrequest);
    CHAMELEON_Sequence_Wait(msequence);
    //flops = flops + FLOPS_DPOTRF(nZmiss);
    VERBOSE(" Done.\n");
    STOP_TIMING(mat_gen_time);

    START_TIMING(time_solve);
    //Start prediction
    VERBOSE("Calculate dposv C22 Covariance Matrix... (Prediction Stage - single precision)");
    CHAMELEON_sposv_Tile(ChamLower, CHAM_descC22, CHAM_descZobs);
    flops = flops + FLOPS_SPOTRF(nZobs);
    flops = flops + FLOPS_STRSM(ChamLeft, nZobs, nZobs);
    VERBOSE(" Done.\n");
    STOP_TIMING(time_solve);

    START_TIMING(time_gemm);
    VERBOSE("Calculate dgemm Zmiss= C12 * Zobs Covariance Matrix... (Prediction Stage - single precision)");
    CHAMELEON_sgemm_Tile(ChamNoTrans, ChamNoTrans, 1, CHAM_descC12, CHAM_descZobs, 0, CHAM_descZmiss);
    flops = flops + FLOPS_SGEMM(nZmiss, nZobs, nZobs);
    VERBOSE(" Done.\n");
    STOP_TIMING(time_gemm);

    //return back descZmiss to zmiss vector
    CHAMELEON_Tile_to_Lapack(CHAM_descZmiss, Zmiss, nZmiss);
    //Estimate Mean Square Error
    if (Zactual != NULL) {
        START_TIMING(time_mse);
        VERBOSE("Calculate Mean Square Error (MSE) ... (Prediction Stage) \n");
        EXAGEOSTAT_MLE_smse_Tile_Async(CHAM_descZactual, CHAM_descZmiss, CHAM_descmse, msequence, mrequest);
        CHAMELEON_Sequence_Wait(msequence);
        VERBOSE(" Done.\n");
        STOP_TIMING(time_mse);
        data->mserror /= nZmiss;
    } else
        data->mserror = -1;

#if defined(CHAMELEON_USE_MPI)
    if(CHAMELEON_My_Mpi_Rank() == 0)
    {
#endif
    if (data->log == 1)
        fprintf(data->pFileLog,
                "\n\n# of missing observations :%d\n\nPrediction Execution Time: %2.6f, Flops: %2.6f, Mean Square Error (MSE): %2.6f\n\n",
                nZmiss, (mat_gen_time + time_solve + time_mse), (flops / 1e9 / (time_solve)), data->mserror);

    write_prediction_result("predict_result.dat", n, data->hicma_acc, 0, 0, 0, data->mserror,
                            (mat_gen_time + time_solve + time_gemm), (flops / 1e9 / (time_solve)));

#if defined(CHAMELEON_USE_MPI)
    }
#endif
    return data->mserror;
}

double EXAGEOSTAT_smle_Predict_Tile_Async(MLE_data *CHAM_data, double*theta, int nZmiss, int nZobs, double*Zobs,
                                         double*Zactual, double*Zmiss, int n)
//! //Predict missing values base on a set of given values and covariance matrix
/*!  -- CHAM-Async
 * Returns the prediction Mean Square Error (MSE) as double
 * @param[in] CHAM_data: MLE_data struct with different MLE inputs.
 * @param[in] theta: theta Vector with three parameter (Variance, Range, Smoothness)
 *                           that is used to to generate the Covariance Matrix.
 * @param[in] nZmiss: number of missing values (unknown observations).
 * @param[in] nZobs: number of observed values (known observations).
 * @param[in] n: number of spatial locations.
 * */
{

    //initialization
    double time_solve = 0.0;
    double mat_gen_time = 0.0;
    double time_gemm = 0.0;
    double time_mse = 0.0;
    double flops = 0.0;

    CHAM_desc_t *CHAM_descZmiss = NULL;
    CHAM_desc_t *CHAM_descC12 = NULL;
    CHAM_desc_t *CHAM_descC22 = NULL;
    CHAM_desc_t *CHAM_descmse = NULL;
    CHAM_desc_t *CHAM_descZactual = NULL;
    CHAM_desc_t *CHAM_descZobs = NULL;
    MLE_data *data = (MLE_data *) CHAM_data;
    RUNTIME_sequence_t *msequence = (RUNTIME_sequence_t *) data->sequence;
    RUNTIME_request_t *mrequest = (RUNTIME_request_t *) data->request;
    data->mserror = 0;

    if (nZmiss <= 0) {
        fprintf(stderr, " Number of missing values should be positive value\n");
        return -1;
    }

    //Initiate data descriptors.
    CHAM_descZmiss = data->descZmiss;
    CHAM_descC12 = data->descC12;
    CHAM_descC22 = data->descC22;
    CHAM_descmse = data->descmse;
    CHAM_descZactual = data->descZactual;
    CHAM_descZobs = data->descZobs;

    //Copy data to vectors.
    VERBOSE("Copy measurments vector to descZobs descriptor (single precision)...");
    //CHAMELEON_Lapack_to_Tile_Async( Zobs, nZobs, CHAM_descZobs, msequence, mrequest);
    VERBOSE(" Done.\n");

    if (Zactual != NULL) {
        VERBOSE("Copy actual measurments vector to descZactual descriptor (single precision)...");
        VERBOSE(" Done.\n");
    }


    START_TIMING(mat_gen_time);
    //Generate C22 covariance matrix.
    VERBOSE("Generate C22 Covariance Matrix... (Prediction Stage - single precision)");
    EXAGEOSTAT_MLE_scmg_Tile_Async(ChamLower, CHAM_descC22, &data->lobs, &data->lobs, &data->lm, theta, data->dm,
                                  data->kernel_fun, msequence, mrequest);
    VERBOSE(" Done.\n");

    //Generate C12 covariance matrix.
    VERBOSE("Generate C12 Covariance Matrix... (Prediction Stage - single precision)");
    EXAGEOSTAT_MLE_scmg_Tile_Async(ChamLower, CHAM_descC12, &data->lmiss, &data->lobs, &data->lm, theta, data->dm,
                                  data->kernel_fun, msequence, mrequest);
    CHAMELEON_Sequence_Wait(msequence);
    VERBOSE(" Done.\n");
    STOP_TIMING(mat_gen_time);

    START_TIMING(time_solve);
    //Start prediction.
    VERBOSE("Calculate dposv C22 Covariance Matrix... (Prediction Stage - single precision)");
    CHAMELEON_sposv_Tile_Async(ChamLower, CHAM_descC22, CHAM_descZobs, msequence, mrequest);
    flops = flops + FLOPS_SPOTRF(nZobs);
    flops = flops + FLOPS_STRSM(ChamLeft, nZobs, nZobs);
    VERBOSE(" Done.\n");
    STOP_TIMING(time_solve);

    START_TIMING(time_gemm);
    VERBOSE("Calculate dgemm Zmiss= C12 * Zobs Covariance Matrix... (Prediction Stage - single precision)");
    void *ws;
    CHAMELEON_sgemm_Tile_Async(ChamNoTrans, ChamNoTrans, 1, CHAM_descC12, CHAM_descZobs, 0, CHAM_descZmiss, ws,
                               msequence, mrequest);
    flops = flops + FLOPS_SGEMM(nZmiss, nZobs, nZobs);
    VERBOSE(" Done.\n");
    STOP_TIMING(time_gemm);

    //Estimate Mean Square Error
    if (Zactual != NULL) {
        START_TIMING(time_mse);
        VERBOSE("Calculate Mean Square Error (MSE) ... (Prediction Stage - single precision) \n");
        EXAGEOSTAT_MLE_smse_Tile_Async(CHAM_descZactual, CHAM_descZmiss, CHAM_descmse, msequence, mrequest);
        VERBOSE(" Done.\n");
        STOP_TIMING(time_mse);
        data->mserror /= nZmiss;
    } else
        data->mserror = -1;

#if defined(CHAMELEON_USE_MPI)
    if(CHAMELEON_My_Mpi_Rank() == 0)
    {
#endif
    if (data->log == 1)
        fprintf(data->pFileLog,
                "\n\n# of missing observations :%d\n\nPrediction Execution Time: %2.6f, Flops: %2.6f, Mean Square Error (MSE): %2.6f\n\n",
                nZmiss, (mat_gen_time + time_solve + time_mse), (flops / 1e9 / (time_solve)), data->mserror);

    write_prediction_result("predict_result.dat", n, data->hicma_acc, 0, 0, 0, data->mserror,
                            (mat_gen_time + time_solve + time_gemm), (flops / 1e9 / (time_solve)));
#if defined(CHAMELEON_USE_MPI)
    }
#endif

    return data->mserror;
}

//init Chameleon descriptors
void
EXAGEOSTAT_smle_Call(MLE_data *data, int ncores, int gpus, int dts, int p_grid, int q_grid, int N, int nZobs, int nZmiss)
//! //Initiate CHAM and allocate different descriptors for
/*!  CHAMELEON
 * Returns MLE_data data with initial values and new descriptors locations.
 * @param[in] data: MLE_data struct with different MLE inputs.
 * @param[in] ncores: number of CPU workers.
 * @param[in] gpus: number of GPU workers.
 * @param[in] dts: tile size (MB).
 * @param[in] p_grid: p_grid.
 * @param[in] q_grid: q_grid.
 * @param[in] N: number of spatial locations.	
 * @param[in] nZmiss: number of missing values (unknown observations).
 * @param[in] nZobs: number of observed values (known observations).
 * */
{

    RUNTIME_sequence_t *msequence;
    RUNTIME_request_t mrequest[2] = {CHAMELEON_SUCCESS, CHAMELEON_SUCCESS};
    CHAM_desc_t *CHAM_descC = NULL;
    CHAM_desc_t *CHAM_descZ = NULL;
    CHAM_desc_t *CHAM_descZcpy = NULL;
    CHAM_desc_t *CHAM_descproduct = NULL;
    CHAM_desc_t *CHAM_descdet = NULL;

    // For ditributed system and should be removed
    float *Zcpy = (float *) malloc(N * sizeof(float));

    //Identifies a set of routines sharing common exception handling.
    CHAMELEON_Sequence_Create(&msequence);

    EXAGEOSTAT_ALLOCATE_MATRIX_TILE(&CHAM_descC, NULL, ChamRealFloat, dts, dts, dts * dts, N, N, 0, 0, N, N, p_grid,
                                    q_grid);
    EXAGEOSTAT_ALLOCATE_MATRIX_TILE(&CHAM_descZ, NULL, ChamRealFloat, dts, dts, dts * dts, N, 1, 0, 0, N, 1, p_grid,
                                    q_grid);
    EXAGEOSTAT_ALLOCATE_MATRIX_TILE(&CHAM_descZcpy, Zcpy, ChamRealFloat, dts, dts, dts * dts, N, 1, 0, 0, N, 1, p_grid,
                                    q_grid);
    EXAGEOSTAT_ALLOCATE_MATRIX_TILE(&CHAM_descproduct, &data->sdotp, ChamRealFloat, dts, dts, dts * dts, 1, 1, 0, 0, 1,
                                    1, p_grid, q_grid);
    EXAGEOSTAT_ALLOCATE_MATRIX_TILE(&CHAM_descdet, &data->det, ChamRealFloat, dts, dts, dts * dts, 1, 1, 0, 0, 1, 1,
                                    p_grid, q_grid);

    //Fill data struct
    data->descC = CHAM_descC;
    data->descZ = CHAM_descZ;
    data->descZcpy = CHAM_descZcpy;
    data->descdet = CHAM_descdet;
    data->descproduct = CHAM_descproduct;
    data->sequence = msequence;
    data->request = mrequest;
    //stop gsl error handler
    gsl_set_error_handler_off();
}