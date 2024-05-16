/**
 *
 * Copyright (c) 2017-2023  King Abdullah University of Science and Technology
 * All rights reserved.
 *
 * ExaGeoStat is a software package provided by KAUST
 **/
/**
 *
 * @file MLE_misc.c
 *
 * Auxiliary functions that are needed by ExaGeoStat.
 *
 * @version 1.2.0
 *
 * @author Sameh Abdulah
 * @date 2022-11-09
 *
 **/
#include "../include/MLE_misc.h"
#ifdef __cplusplus
#include "../include/inference.hpp"
#endif
#include <../src/include/MLE.h>
//***************************************************************************************
#define pow_e (2.71828182845904)
void init_data_values(MLE_data *data)
//! initiate data struct with default values
/*!
 * @param[in] data: MLE_data struct with different MLE inputs.
 * */
{
    data->variance = 1.0;                            ///< Variance parameter.
    data->variance1 = 1.0;                           ///< Variance1 parameter.
    data->variance2 = 1.0;                           ///< Variance2 parameter.
    data->variance3 = 1.0;                           ///< Variance3 parameter.
    data->computation = "exact";                     ///< Exact or approx computation.
    data->c_fun = "univariate_matern_stationary";    ///< Matern or pow-exp kernels.
    data->test = 0;                                  ///< Test or real-data computation.
    data->iter_count = 0;                            ///< Number of iterations to converge.
    data->l1.x = NULL;                               ///< 2D locations for the first dataset (X vector).
    data->l1.y = NULL;                               ///< 2D locations for the first dataset (Y vector).
    data->l1.z = NULL;                               ///< 3D locations for the first dataset (Z vector)
    data->lmiss.x = NULL;                            ///< 2D locations for the missing data (X vector) (prediction stage).
    data->lmiss.y = NULL;                            ///< 2D locations for the missing data (Y vector) (prediction stage).
    data->lmiss.z = NULL;                            ///< 2D locations for the missing data (Z vector) (prediction stage).
    data->lobs.x = NULL;                             ///< 2D locations for the observed data (X vector) (prediction stage).
    data->lobs.y = NULL;                             ///< 2D locations for the observed data (Y vector) (prediction stage).
    data->lobs.z = NULL;                             ///< 2D locations for the observed data (Z vector) (prediction stage).
    data->lm.x = NULL;                               ///< 2D locations for the median data point (X vector).
    data->lm.y = NULL;                               ///< 2D locations for the median data point (Y vector).
    data->lm.z = NULL;                               ///< 2D locations for the median data point (Z vector).
    data->ns_params.sigma_arr1      = NULL;          ///< sigma array for loc 1. (nonstat case)
    data->ns_params.lamda_arr1     = NULL;           ///< lamda array for loc 1. (nonstat case)
    data->ns_params.sigma_arr2      = NULL;          ///< sigma array for loc 2. (nonstat case) 
    data->ns_params.lamda_arr2      = NULL;          ///< lamda array for loc 2. (nonstat case)
    data->ns_params.nu_arr1         = NULL;          ///< nu array for loc 1. (nonstat case) 
    data->ns_params.nu_arr2         = NULL;          ///< nu array for loc 2. (nonstat case)
    data->descC = NULL;                              ///< Covariance matrix C descriptor.
    data->descZ = NULL;                              ///< Measurements Z descriptor.
    data->descZ1 = NULL;                             ///< Measurements Z1 descriptor.
    data->descZ2 = NULL;                             ///< Measurements Z2 descriptor.
    data->descZ3 = NULL;                             ///< Measurements Z3 descriptor.
    data->Adense = NULL;                             ///< Dense matrix descriptor in the case of approximation mode - accuracy check.
    data->Adense2 = NULL;                            ///< Dense matrix descriptor2 in the case of approximation mode - accuracy check.
    data->descZcpy = NULL;                           ///< A copy of Measurements Z descriptor.
    data->descdet = NULL;                            ///< Determinant descriptor.
    data->descsum - NULL;                            ///< Sum descriptor (non-Gaussian).
    data->descproduct = NULL;                        ///< Dot product descriptor.
    data->descproduct1 = NULL;                       ///< Dot product descriptor.
    data->descproduct2 = NULL;                       ///< Dot product descriptor.
    data->descproduct3 = NULL;                       ///< Dot product descriptor.
    data->descZmiss = NULL;                          ///< Missing measurements descriptor.
    data->descC11 = NULL;                            ///< Covariance Matrix C11 descriptor.
    data->descC21 = NULL;                            ///< Covariance Matrix C21 descriptor.
    data->descC12 = NULL;                            ///< Covariance Matrix C12 descriptor.
    data->descC22 = NULL;                            ///< Covariance Matrix C22 descriptor.
    data->descZactual = NULL;                        ///< Actual Measurements Z descriptor.
    data->descZobs = NULL;                           ///< observed Measurements Z descriptor.
    data->descmse = NULL;                            ///< Mean Square Error (MSE) descriptor.
    data->descmse1 = NULL;                           ///< Mean Square Error (MSE) descriptor.
    data->descmse2 = NULL;                           ///< Mean Square Error (MSE) descriptor.
    data->descmse3 = NULL;                           ///< Mean Square Error (MSE) descriptor.
    data->descZtrace = NULL;                         ///< Trace descriptor.
    data->sequence = NULL;                           ///< CHAMELEON sequence.
    data->request = NULL;                            ///< CHAMELEON request.
    data->verbose = 0;                               ///< Verbose indicator.
    data->check = 0;                                 ///< Check indicator -- approximation mode.
    data->log = 0;                                   ///< Log files generation indicator, 0-->no, 1-->yes.
    data->avg_exec_time_per_iter = 0;                ///< Average execution time per iteration (only used in verbose mode).
    data->total_exec_time = 0.0;                     ///< Total execution time (only used in verbose mode).
    data->avg_flops_per_iter = 0.0;                  ///< Average flops per iteration (only used in verbose mode).
    data->avg_generation_per_iter = 0.0;
    data->avg_solve_per_iter = 0.0;
    data->avg_cholesky_per_iter = 0.0;
    data->final_loglik = 0.0;                        ///< Final log likelihood value.
    data->locsFPath = "";                            ///< Locations file path -- in the case of real dataset (real mode).
    data->timeFPath = "";                            ///< Time file path -- in the case of real dataset (real mode -- space-time kernel).
    data->obsFPath = "";                             ///< Observations file path --  in the case of real dataset (real mode).
    data->obsFPath2 = "";                            ///< Observations file path --  in the case of real dataset (real mode).
    data->obsFPath3 = "";                            ///< Observations file path --  in the case of real dataset (real mode).
    data->actualZFPath = "";                         ///< Actual observations file path -- in the case of prediction.
    data->actualZFPath2 = "";                        ///< Actual observations file path -- in the case of prediction.
    data->actualZFPath3 = "";                        ///< Actual observations file path -- in the case of prediction.
    data->actualZLocFPath = "";                      ///< Actual locations file path -- in the case of prediction.
    data->det = 0.0;                                 ///< Determinant value.
    data->sum = 0.0;                                 ///< Sum value (non-Gaussian).
    data->dotp = 0.0;                                ///< Dot product value.
    data->dotp1 = 0.0;                               ///< Double dot1 product value.
    data->dotp2 = 0.0;                               ///< Double dot2 product value.
    data->dotp3 = 0.0;                               ///< Double dot3 product value.
    data->sdotp = 0.0;                               ///< Single dot product value.
    data->mserror = 0.0;                             ///< Mean Square Error (MSE) value.
    data->mserror1 = 0.0;                            ///< Mean Square Error (MSE) value.
    data->mserror2 = 0.0;                            ///< Mean Square Error (MSE) value.
    data->mserror3 = 0.0;                            ///< Mean Square Error (MSE) value.
    data->dm = "ed";                                 ///< Distance metric to be used ed->Euclidian Distance -- gcd->Great Circle Distance.
    data->diag_thick = 0;                            ///< The thick of used diagonal in the case of diagonal approximation approach.
    data->nFileLog = NULL;                           ///< Log file name (only used if log -->1).
    data->pFileLog = NULL;                           ///< Log file path (only used if log -->1).
    data->hicma_maxrank = 0;                         ///< Max Rank in the case of LR-HiCMA approx
    data->hicma_data_type = 0;                       ///< To define the problem typr to HiCMA (HICMA_STARSH_PROB_GEOSTAT (Synthetic) or HICMA_STARSH_PROB_GEOSTAT_POINT (real))
    data->hicma_descC = NULL;                        ///< HiCMA descC descriptor (for accuracy check).
    data->hicma_descZ = NULL;                        ///< HiCMA descZ descriptor.
    data->hicma_descCD = NULL;                       ///< HiCMA descCD descriptor.
    data->hicma_descCUV = NULL;                      ///< HiCMA descCUV descriptor.
    data->hicma_descCrk = NULL;                      ///< HiCMA descCrk descriptor.
    data->hicma_descZcpy = NULL;                     ///< A copy of Measurements Z descriptor.
    data->hicma_descdet = NULL;                      ///< Determinant descriptor.
    data->hicma_descproduct = NULL;                  ///< Dot product descriptor.
    data->hicma_descC12D = NULL;                     ///< HiCMA descCD descriptor.
    data->hicma_descC12UV = NULL;                    ///< HiCMA descCUV descriptor.
    data->hicma_descC12rk = NULL;                    ///< HiCMA descCrk descriptor.
    data->hicma_descC22D = NULL;                     ///< HiCMA descCD descriptor.
    data->hicma_descC22UV = NULL;                    ///< HiCMA descCUV descriptor.
    data->hicma_descC22rk = NULL;                    ///< HiCMA descCrk descriptor.
    data->hicma_acc = 9;                             ///< Accuracy in the case of LR-HiCMA approx.
    data->hsequence = NULL;                          ///< HiCMA sequence.
    data->hrequest = NULL;                           ///< HiCMA request.
    data->opt_tol = 5;                               ///< The parameter tol is a tolerance that is used for the purpose of stopping criteria only.
    data->opt_max_iters = -1;                        ///< Maximum number of mle iterations.
    data->ooc = 0;                                   ///< Support Out-Of-Core execution, 0-->no, 1-->yes.
    data->kernel_fun = "univariate_matern_stationary";///< stationary_matern, or non_stationary_matern.
    data->precision = 0;                             ///< Double, single, or mixed.
    //Mixed Precision
    data->desctemp = NULL;                           ///< Temporary descriptor for mixed precision Cholesky factorization.
    data->desctemp22 = NULL;                         ///< Temporary descriptor for mixed precision Cholesky factorization.
    //MLOE and MMOM
    data->desck_t = NULL;
    data->desck_a = NULL;
    data->desck_ttmp = NULL;
    data->desck_atmp = NULL;
    data->descK_t = NULL;
    data->descK_a = NULL;
    data->descexpr1 = NULL;
    data->descexpr2 = NULL;
    data->descexpr3 = NULL;
    data->descexpr4 = NULL;
    data->descestimatedalpha = NULL;
    data->desctruthalpha = NULL;
    data->desc_mloe_mmom = NULL;
    data->descmloe = NULL;
    data->descmmom = NULL;
    data->expr1 = 0.0;
    data->expr2 = 0.0;
    data->expr3 = 0.0;
    data->expr4 = 0.0;
    data->mloe = 0.0;
    data->mmom = 0.0;
    data->mloe_mmom = 0;
    data->mloe_mmom_async = 0;
    data->mspe = 0;
    data->recovery_file = "";
    data->checkpoint_file = "";
    data->recovery_file = "";
    data->time_slots = 0;
    data->idw = 0;
    data->fisher = 0;
    //for non-Gaussian prediction
    data->descr = NULL;
    data->descrcpy = NULL;

    //init result struct
    results.problem_size = -1;
    results.computation = NULL;
    results.kernel = NULL;
    results.ds_type = NULL;
    results.precision = NULL;
    results.z_sample = -1;
    results.dense_ts = -1;
    results.lr_ts = -1;
    results.lr_acc = -1;
    results.lr_maxrank = -1;
    results.ncores = -1;
    results.ngpus = -1;
    results.p = -1;
    results.q = -1;
    results.num_params = -1;
    results.initial_theta = NULL;
    results.starting_theta = NULL;
    results.estimated_theta = NULL;
    results.final_loglik = -1;
    results.time_per_iteration = -1;
    results.cholesky_per_iter = -1;
    results.flops_per_iteration = -1;
    results.generation_per_iter = -1;
    results.solve_per_iter = -1;
    results.total_mle_time = -1;
    results.mse_pred1 = -1;
    results.mse_pred2 = -1;
    results.mse_pred3 = -1;
    results.mse_pred = -1;
    results.trace_pred_sum = -1;
    results.trace_pred_mean = -1;
    results.total_pred_time = -1;
    results.total_pred_flops = -1;
    results.mloe = -1;
    results.mmom = -1;
    results.mloe_exec = NULL;
    results.total_mloe_mmom_time = -1;
    results.matrix_gen_mloe_mmom_time = -1;
    results.cho_fact_mloe_mmom_time = -1;
    results.loop_mloe_mmom_time = -1;
    results.total_mloe_mmom_flops = -1;
    results.fisher_00 = -1;
    results.fisher_11 = -1;
    results.fisher_22 = -1;
    results.fisher_time = -1;
}


static uint32_t Compact1By1(uint32_t x)
//! Collect every second bit into lower part of input
{
    x &= 0x55555555;
    // x = -f-e -d-c -b-a -9-8 -7-6 -5-4 -3-2 -1-0
    x = (x ^ (x >> 1)) & 0x33333333;
    // x = --fe --dc --ba --98 --76 --54 --32 --10
    x = (x ^ (x >> 2)) & 0x0f0f0f0f;
    // x = ---- fedc ---- ba98 ---- 7654 ---- 3210
    x = (x ^ (x >> 4)) & 0x00ff00ff;
    // x = ---- ---- fedc ba98 ---- ---- 7654 3210
    x = (x ^ (x >> 8)) & 0x0000ffff;
    // x = ---- ---- ---- ---- fedc ba98 7654 3210
    return x;
}

static uint32_t DecodeMorton2X(uint32_t code)
//! Decode first input
{
    return Compact1By1(code >> 0);
}

static uint32_t DecodeMorton2Y(uint32_t code)
//! Decode second input
{
    return Compact1By1(code >> 1);
}

static uint32_t Part1By1(uint32_t x)
//! Spread lower bits of input
{
    x &= 0x0000ffff;
    // x = ---- ---- ---- ---- fedc ba98 7654 3210
    x = (x ^ (x << 8)) & 0x00ff00ff;
    // x = ---- ---- fedc ba98 ---- ---- 7654 3210
    x = (x ^ (x << 4)) & 0x0f0f0f0f;
    // x = ---- fedc ---- ba98 ---- 7654 ---- 3210
    x = (x ^ (x << 2)) & 0x33333333;
    // x = --fe --dc --ba --98 --76 --54 --32 --10
    x = (x ^ (x << 1)) & 0x55555555;
    // x = -f-e -d-c -b-a -9-8 -7-6 -5-4 -3-2 -1-0
    return x;
}

static int compare_uint32(const void *a, const void *b)
//! Compare two uint32_t
{
    uint32_t _a = *(uint32_t *) a;
    uint32_t _b = *(uint32_t *) b;
    if (_a < _b) return -1;
    if (_a == _b) return 0;
    return 1;
}


static int compare_sdata(const void *a, const void *b)
//! Compare two sdata
{
    sdata _a = *(sdata *) a;
    sdata _b = *(sdata *) b;
    if (_a.xy < _b.xy) return -1;
    if (_a.xy == _b.xy) return 0;
    return 1;
}


static uint32_t EncodeMorton2(uint32_t x, uint32_t y)
//! Encode two inputs into one
{
    return (Part1By1(y) << 1) + Part1By1(x);
}

static void zsort_locations(int n, location *locations)
//! Sort in Morton order (input points must be in [0;1]x[0;1] square])
{
    // Some sorting, required by spatial statistics code
    int i;
    uint16_t x, y;
    uint32_t z[n];
    // Encode data into vector z
    for (i = 0; i < n; i++) {
        x = (uint16_t) (locations->x[i] * (double) UINT16_MAX + .5);
        y = (uint16_t) (locations->y[i] * (double) UINT16_MAX + .5);
        z[i] = EncodeMorton2(x, y);
    }
    // Sort vector z
    qsort(z, n, sizeof(uint32_t), compare_uint32);
    // Decode data from vector z
    for (i = 0; i < n; i++) {
        x = DecodeMorton2X(z[i]);
        y = DecodeMorton2Y(z[i]);
        locations->x[i] = (double) x / (double) UINT16_MAX;
        locations->y[i] = (double) y / (double) UINT16_MAX;
    }
}

static uint64_t Part1By3(uint64_t x)
// Spread lower bits of input
{
    x &= 0x000000000000ffff;
    // x = ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- fedc ba98 7654 3210
    x = (x ^ (x << 24)) & 0x000000ff000000ff;
    // x = ---- ---- ---- ---- ---- ---- fedc ba98 ---- ---- ---- ---- ---- ---- 7654 3210
    x = (x ^ (x << 12)) & 0x000f000f000f000f;
    // x = ---- ---- ---- fedc ---- ---- ---- ba98 ---- ---- ---- 7654 ---- ---- ---- 3210
    x = (x ^ (x << 6)) & 0x0303030303030303;
    // x = ---- --fe ---- --dc ---- --ba ---- --98 ---- --76 ---- --54 ---- --32 ---- --10
    x = (x ^ (x << 3)) & 0x1111111111111111;
    // x = ---f ---e ---d ---c ---b ---a ---9 ---8 ---7 ---6 ---5 ---4 ---3 ---2 ---1 ---0
    return x;
}

static uint64_t EncodeMorton3(uint64_t x, uint64_t y, uint64_t z)
// Encode 3 inputs into one
{
    return (Part1By3(z) << 2) + (Part1By3(y) << 1) + Part1By3(x);
}

static uint64_t Compact1By3(uint64_t x)
// Collect every 4-th bit into lower part of input
{
    x &= 0x1111111111111111;
    // x = ---f ---e ---d ---c ---b ---a ---9 ---8 ---7 ---6 ---5 ---4 ---3 ---2 ---1 ---0
    x = (x ^ (x >> 3)) & 0x0303030303030303;
    // x = ---- --fe ---- --dc ---- --ba ---- --98 ---- --76 ---- --54 ---- --32 ---- --10
    x = (x ^ (x >> 6)) & 0x000f000f000f000f;
    // x = ---- ---- ---- fedc ---- ---- ---- ba98 ---- ---- ---- 7654 ---- ---- ---- 3210
    x = (x ^ (x >> 12)) & 0x000000ff000000ff;
    // x = ---- ---- ---- ---- ---- ---- fedc ba98 ---- ---- ---- ---- ---- ---- 7654 3210
    x = (x ^ (x >> 24)) & 0x000000000000ffff;
    // x = ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- fedc ba98 7654 3210
    return x;
}

static uint64_t DecodeMorton3X(uint64_t code)
// Decode first input
{
    return Compact1By3(code >> 0);
}

static uint64_t DecodeMorton3Y(uint64_t code)
// Decode second input
{
    return Compact1By3(code >> 1);
}

static uint64_t DecodeMorton3Z(uint64_t code)
// Decode third input
{
    return Compact1By3(code >> 2);
}

static void zsort_locations_3d(int n, location *locations)
//! Sort in Morton order (input points must be in [0;1]x[0;1] square])
{
    // Some sorting, required by spatial statistics code
    int i;
    uint16_t x, y, z;
    uint64_t Z[n];
    // Encode data into vector z
    for (i = 0; i < n; i++) {
        x = (uint16_t) (locations->x[i] * (double) UINT16_MAX + .5);
        y = (uint16_t) (locations->y[i] * (double) UINT16_MAX + .5);
        z = (uint16_t) (locations->z[i] * (double) UINT16_MAX + .5);
        Z[i] = EncodeMorton3(x, y, z);
    }
    // Sort vector z
    qsort(Z, n, sizeof(uint64_t), compare_uint32);
    // Decode data from vector z
    for (i = 0; i < n; i++) {
        x = DecodeMorton3X(Z[i]);
        y = DecodeMorton3Y(Z[i]);
        z = DecodeMorton3Z(Z[i]);
        locations->x[i] = (double) x / (double) UINT16_MAX;
        locations->y[i] = (double) y / (double) UINT16_MAX;
        locations->z[i] = (double) z / (double) UINT16_MAX;
    }
}


static void radix_sort_recursive(uint32_t *data, int count, int ndim,
                                 int *order, int *tmp_order, int sdim, int sbit,
                                 int lo, int hi)
// Hierarchical radix sort to get Z-order of particles.
// This function is static not to be visible outside this module.
{
    int i, lo_last = lo, hi_last = hi;
    uint32_t *sdata = data + sdim * count;
    uint32_t check = 1 << sbit;
    for (i = lo; i <= hi; i++) {
        if ((sdata[order[i]] & check) == 0) {
            tmp_order[lo_last] = order[i];
            lo_last++;
        } else {
            tmp_order[hi_last] = order[i];
            hi_last--;
        }
    }
    for (i = lo; i <= hi; i++)
        order[i] = tmp_order[i];
    if (sdim > 0) {
        if (lo_last - lo > 1)
            radix_sort_recursive(data, count, ndim, order, tmp_order, sdim - 1,
                                 sbit, lo, lo_last - 1);
        if (hi - hi_last > 1)
            radix_sort_recursive(data, count, ndim, order, tmp_order, sdim - 1,
                                 sbit, hi_last + 1, hi);
    } else if (sbit > 0) {
        if (lo_last - lo > 1)
            radix_sort_recursive(data, count, ndim, order, tmp_order, ndim - 1,
                                 sbit - 1, lo, lo_last - 1);
        if (hi - hi_last > 1)
            radix_sort_recursive(data, count, ndim, order, tmp_order, ndim - 1,
                                 sbit - 1, hi_last + 1, hi);
    }
}

static int radix_sort(uint32_t *data, int count, int ndim,
                      int *order)
// Auxiliary sorting function for starsh_particles_zsort_inpace().
// This function is static not to be visible outside this module.
{
    int *tmp_order = (int *) malloc(count * sizeof(int));
    radix_sort_recursive(data, count, ndim, order, tmp_order, ndim - 1, 31, 0, count - 1);
    free(tmp_order);
    return 0;
}

int locations_obs_zsort_inplace(int n, location *locations, double* z)
//! Sort particles in Z-order (Morton order).
/*! This function must be used after initializing @ref STARSH_particles with
 * your own data by starsh_particles_init() or starsh_particles_new().
 *
 * @sa starsh_particles_init(), starsh_particles_new().
 * @ingroup app-particles
 * */
{
    int i;
    int j;//new_j, tmp_j;
    int count = n;
    int ndim = 2;
    int info;
    double* point = (double* ) malloc(ndim * count * sizeof(double));
    double* ptr1;
    double* minmax; // min is stored in lower part, max is stored in upper part
    //double tmp_x;
    for (i = 0; i < count; i++) {
        point[i] = locations->x[i];
        point[i + count] = locations->y[i];
    }

    // I think count should be replaced by ndim
    minmax = (double* ) malloc(2 * count * sizeof(double));

    for (i = 0; i < ndim; i++) {
        ptr1 = point + i * count; // i-th dimension
        minmax[i] = ptr1[0];
        minmax[i + ndim] = minmax[i];
        for (j = 1; j < count; j++) {
            if (minmax[i] > ptr1[j])
                minmax[i] = ptr1[j];
            else if (minmax[i + ndim] < ptr1[j])
                minmax[i + ndim] = ptr1[j];
        }
    }
    // Now minmax[0:ndim] and minmax[ndim:2*ndim] store minimal and maximal
    // values of coordinates
    uint32_t *uint_point = (uint32_t *) malloc(ndim * count * sizeof(uint32_t));
    uint32_t *uint_ptr1;
    double min, range;
    for (i = 0; i < ndim; i++) {
        uint_ptr1 = uint_point + i * count;
        ptr1 = point + i * count;
        min = minmax[i];
        range = minmax[i + ndim] - min;
        for (j = 0; j < count; j++)
            uint_ptr1[j] = (ptr1[j] - min) / range * UINT32_MAX;
    }
    free(minmax);
    // Now uint_ptr1 contains initial coordinates, rescaled to range
    // [0, UINT32_MAX] and converted to uint32_t type to use special radix sort
    // Prepare indexes to store sort order
    int *order = (int *) malloc(count * sizeof(int));
    for (j = 0; j < count; j++)
        order[j] = j;
    info = radix_sort(uint_point, count, ndim, order);
    if (info != 0) {
        free(uint_point);
        free(order);
        return info;
    }
    double* new_point = (double* ) malloc(ndim * count * sizeof(double));
    double* new_z = (double* ) malloc(count * sizeof(double));
    for (j = 0; j < count; j++) {
        for (i = 0; i < ndim; i++)
            new_point[count * i + j] = point[count * i + order[j]];
        new_z[j] = z[order[j]];
    }
    for (i = 0; i < count; i++) {
        locations->x[i] = new_point[i];
        locations->y[i] = new_point[i + count];
        z[i] = new_z[i];
    }
    free(new_point);
    free(point);
    free(new_z);
    free(uint_point);
    free(order);
    return 0;
}


int locations_obs_zsort_inplace_bivariate(int n, location *locations, double* z, double* z2)
//! Sort particles in Z-order (Morton order).
/*! This function must be used after initializing @ref STARSH_particles with
 * your own data by starsh_particles_init() or starsh_particles_new().
 *
 * @sa starsh_particles_init(), starsh_particles_new().
 * @ingroup app-particles
 * */
{
    int i;
    int j;//new_j, tmp_j;
    int count = n;
    int ndim = 2;
    int info;
    double* point = (double* ) malloc(ndim * count * sizeof(double));
    double* ptr1;
    double* minmax; // min is stored in lower part, max is stored in upper part
    //double tmp_x;
    for (i = 0; i < count; i++) {
        point[i] = locations->x[i];
        point[i + count] = locations->y[i];
    }

    // I think count should be replaced by ndim
    minmax = (double* ) malloc(2 * count * sizeof(double));

    for (i = 0; i < ndim; i++) {
        ptr1 = point + i * count; // i-th dimension
        minmax[i] = ptr1[0];
        minmax[i + ndim] = minmax[i];
        for (j = 1; j < count; j++) {
            if (minmax[i] > ptr1[j])
                minmax[i] = ptr1[j];
            else if (minmax[i + ndim] < ptr1[j])
                minmax[i + ndim] = ptr1[j];
        }
    }
    // Now minmax[0:ndim] and minmax[ndim:2*ndim] store minimal and maximal
    // values of coordinates
    uint32_t *uint_point = (uint32_t *) malloc(ndim * count * sizeof(uint32_t));
    uint32_t *uint_ptr1;
    double min, range;
    for (i = 0; i < ndim; i++) {
        uint_ptr1 = uint_point + i * count;
        ptr1 = point + i * count;
        min = minmax[i];
        range = minmax[i + ndim] - min;
        for (j = 0; j < count; j++)
            uint_ptr1[j] = (ptr1[j] - min) / range * UINT32_MAX;
    }
    free(minmax);
    // Now uint_ptr1 contains initial coordinates, rescaled to range
    // [0, UINT32_MAX] and converted to uint32_t type to use special radix sort
    // Prepare indexes to store sort order
    int *order = (int *) malloc(count * sizeof(int));
    for (j = 0; j < count; j++)
        order[j] = j;
    info = radix_sort(uint_point, count, ndim, order);
    if (info != 0) {
        free(uint_point);
        free(order);
        return info;
    }
    double* new_point = (double* ) malloc(ndim * count * sizeof(double));
    double* new_z = (double* ) malloc(count * sizeof(double));
    double* new_z2 = (double* ) malloc(count * sizeof(double));
    for (j = 0; j < count; j++) {
        for (i = 0; i < ndim; i++)
            new_point[count * i + j] = point[count * i + order[j]];
        new_z[j] = z[order[j]];
        new_z2[j] = z2[order[j]];
    }
    for (i = 0; i < count; i++) {
        locations->x[i] = new_point[i];
        locations->y[i] = new_point[i + count];
        z[i] = new_z[i];
        z2[i] = new_z2[i];
    }
    free(new_point);
    free(point);
    free(new_z);
    free(new_z2);
    free(uint_point);
    free(order);
    return 0;
}


int locations_obs_zsort_inplace_trivariate(int n, location *locations, double* z, double* z2, double* z3)
//! Sort particles in Z-order (Morton order).
/*! This function must be used after initializing @ref STARSH_particles with
 * your own data by starsh_particles_init() or starsh_particles_new().
 *
 * @sa starsh_particles_init(), starsh_particles_new().
 * @ingroup app-particles
 * */
{
    int i;
    int j;//new_j, tmp_j;
    int count = n;
    int ndim = 2;
    int info;
    double* point = (double* ) malloc(ndim * count * sizeof(double));
    double* ptr1;
    double* minmax; // min is stored in lower part, max is stored in upper part
    //double tmp_x;
    for (i = 0; i < count; i++) {
        point[i] = locations->x[i];
        point[i + count] = locations->y[i];
    }

    // I think count should be replaced by ndim
    minmax = (double* ) malloc(2 * count * sizeof(double));

    for (i = 0; i < ndim; i++) {
        ptr1 = point + i * count; // i-th dimension
        minmax[i] = ptr1[0];
        minmax[i + ndim] = minmax[i];
        for (j = 1; j < count; j++) {
            if (minmax[i] > ptr1[j])
                minmax[i] = ptr1[j];
            else if (minmax[i + ndim] < ptr1[j])
                minmax[i + ndim] = ptr1[j];
        }
    }
    // Now minmax[0:ndim] and minmax[ndim:2*ndim] store minimal and maximal
    // values of coordinates
    uint32_t *uint_point = (uint32_t *) malloc(ndim * count * sizeof(uint32_t));
    uint32_t *uint_ptr1;
    double min, range;
    for (i = 0; i < ndim; i++) {
        uint_ptr1 = uint_point + i * count;
        ptr1 = point + i * count;
        min = minmax[i];
        range = minmax[i + ndim] - min;
        for (j = 0; j < count; j++)
            uint_ptr1[j] = (ptr1[j] - min) / range * UINT32_MAX;
    }
    free(minmax);
    // Now uint_ptr1 contains initial coordinates, rescaled to range
    // [0, UINT32_MAX] and converted to uint32_t type to use special radix sort
    // Prepare indexes to store sort order
    int *order = (int *) malloc(count * sizeof(int));
    for (j = 0; j < count; j++)
        order[j] = j;
    info = radix_sort(uint_point, count, ndim, order);
    if (info != 0) {
        free(uint_point);
        free(order);
        return info;
    }
    double* new_point = (double* ) malloc(ndim * count * sizeof(double));
    double* new_z = (double* ) malloc(count * sizeof(double));
    double* new_z2 = (double* ) malloc(count * sizeof(double));
    double* new_z3 = (double* ) malloc(count * sizeof(double));
    for (j = 0; j < count; j++) {
        for (i = 0; i < ndim; i++)
            new_point[count * i + j] = point[count * i + order[j]];
        new_z[j] = z[order[j]];
        new_z2[j] = z2[order[j]];
        new_z3[j] = z3[order[j]];
    }
    for (i = 0; i < count; i++) {
        locations->x[i] = new_point[i];
        locations->y[i] = new_point[i + count];
        z[i] = new_z[i];
        z2[i] = new_z2[i];
        z3[i] = new_z3[i];

    }
    free(new_point);
    free(point);
    free(new_z);
    free(new_z2);
    free(new_z3);
    free(uint_point);
    free(order);
    return 0;
}
//! Sort particles in Z-order (Morton order).
/*! This function must be used after initializing @ref STARSH_particles with
 * your own data by starsh_particles_init() or starsh_particles_new().
 *
 * @sa starsh_particles_init(), starsh_particles_new().
 * @ingroup app-particles
 * */

void zsort_locations_obs(int n, location *locations, double* z)
//! Sort in Morton order (input points must be in [0;1]x[0;1] square])
{
    // Some sorting, required by spatial statistics code
    int i;
    uint16_t x, y;
    sdata z_struct[n];

    // Encode data into vector z
    for (i = 0; i < n; i++) {
        x = (uint16_t) (locations->x[i] * (double) UINT16_MAX + .5);
        y = (uint16_t) (locations->y[i] * (double) UINT16_MAX + .5);
        z_struct[i].xy = EncodeMorton2(x, y);
        z_struct[i].z = z[i];
    }

    // Sort vector z
    qsort(z_struct, n, sizeof(sdata), compare_sdata);
    // Decode data from vector z
    for (i = 0; i < n; i++) {
        x = DecodeMorton2X(z_struct[i].xy);
        y = DecodeMorton2Y(z_struct[i].xy);
        locations->x[i] = (double) x / (double) UINT16_MAX;
        locations->y[i] = (double) y / (double) UINT16_MAX;
        z[i] = z_struct[i].z;
    }
}


double uniform_distribution(double rangeLow, double rangeHigh)
//! Generate uniform distribution between rangeLow , rangeHigh
{
    // unsigned int *seed = &exageostat_seed;
    double myRand = (double) rand() / (double) (1.0 + RAND_MAX);
    double range = rangeHigh - rangeLow;
    double myRand_scaled = (myRand * range) + rangeLow;
    return myRand_scaled;
}

location *GenerateXYLoc(int n, int seed)
//! Generate XY location for exact computation (MORSE)
{
    //initalization
    int i = 0, index = 0, j = 0;
    srand(seed);
    location *locations = (location *) malloc(sizeof(location *));
    //Allocate memory
    locations->x = (double* ) malloc(n * sizeof(double));
    locations->y = (double* ) malloc(n * sizeof(double));
    locations->z = NULL;

    int sqrtn = ceil(sqrt(n));

    int *grid = (int *) calloc((int) sqrtn, sizeof(int));

    for (i = 0; i < sqrtn; i++) {
        grid[i] = i + 1;
    }

    for (i = 0; i < sqrtn && index < n; i++)
        for (j = 0; j < sqrtn && index < n; j++) {
            locations->x[index] = (grid[i] - 0.5 + uniform_distribution(-0.4, 0.4)) / sqrtn;
            locations->y[index] = (grid[j] - 0.5 + uniform_distribution(-0.4, 0.4)) / sqrtn;
            index++;
        }
    free(grid);
    zsort_locations(n, locations);
    return locations;
}

location *GenerateXYLoc_R(int n)
//! Generate XY location for exact computation (MORSE)
{
    printf("%d\n", n);
    //initalization
    int i = 0, index = 0, j = 0;
    location *locations = (location *) malloc(sizeof(location *));
    //Allocate memory
    locations->x = (double* ) malloc(n * sizeof(double));
    locations->y = (double* ) malloc(n * sizeof(double));
    locations->z = NULL;


    int sqrtn = ceil(sqrt(n));

    printf("%f", sqrtn);
    for (i = 0; i <= sqrtn && index < n; i++)
        for (j = 0; j <= sqrtn && index < n; j++) {
            locations->x[index] = i / (sqrtn - 1);
            locations->y[index] = j / (sqrtn - 1);
            index++;
        }
    zsort_locations(n, locations);

    return locations;
}

location *GenerateXYZLoc(int n, int seed)
//! Generate XY location for exact computation (MOORSE)
{
    //initalization
    int i = 0, index = 0, j = 0, k = 0;
    srand(seed);
    location *locations = (location *) malloc(sizeof(location));
    //Allocate memory
    locations->x = (double* ) malloc(n * sizeof(double));
    locations->y = (double* ) malloc(n * sizeof(double));
    locations->z = (double* ) malloc(n * sizeof(double));


    int cbrtn = ceil(cbrt(n));

    int *grid = (int *) calloc((int) cbrtn, sizeof(int));
    for (i = 0; i < cbrtn; i++) {
        grid[i] = i + 1;
    }

    for (i = 0; i < cbrtn && index < n; i++)
        for (j = 0; j < cbrtn && index < n; j++)
            for (k = 0; k < cbrtn && index < n; k++) {
                locations->x[index] = (grid[i] - 0.5 + uniform_distribution(-0.4, 0.4)) / cbrtn;
                locations->y[index] = (grid[j] - 0.5 + uniform_distribution(-0.4, 0.4)) / cbrtn;
                locations->z[index] = (grid[k] - 0.5 + uniform_distribution(-0.4, 0.4)) / cbrtn;
                index++;
            }
    free(grid);
    zsort_locations_3d(n, locations);
    return locations;
}

/*****************************************************************************************************************/
location *GenerateXYLoc_reg(int n, int seed)
//! Generate XY location for exact computation (MORSE)
{
    //initalization
    int i = 0, index = 0, j = 0;
    srand(seed);
    location *locations = (location *) malloc(sizeof(location *));
    //Allocate memory
    locations->x = (double* ) malloc(n * sizeof(double));
    locations->y = (double* ) malloc(n * sizeof(double));
    locations->z = NULL;

    int sqrtn = ceil(sqrt(n));

    int *grid = (int *) calloc((int) sqrtn, sizeof(int));

    for (i = 0; i < sqrtn; i++) {
        grid[i] = i + 1;
    }

    for (i = 0; i < sqrtn && index < n; i++)
        for (j = 0; j < sqrtn && index < n; j++) {
            locations->x[index] = ((double) grid[i]) * 3 / sqrtn;
            locations->y[index] = ((double) grid[j]) * 3 / sqrtn;
            printf("%f, %f\n", locations->x[index], locations->y[index]);
            index++;
        }
    free(grid);
    zsort_locations(n, locations);
    return locations;
}

location *GenerateXYLoc_ST(int n, int t_slots, int seed)
//! Generate XY location for exact computation (MOORSE)
{
    //initalization
    int i = 0, index = 0, j = 0;

    srand(seed);
    location *locations = (location *) malloc(sizeof(location));
    //Allocate memory
    locations->x = (double* ) malloc(n * t_slots * sizeof(double));
    locations->y = (double* ) malloc(n * t_slots * sizeof(double));
    locations->z = (double* ) malloc(n * t_slots * sizeof(double));

    int sqrtn = ceil(sqrt(n));

    int *grid = (int *) calloc((int) sqrtn, sizeof(int));

    for (i = 0; i < sqrtn; i++) {
        grid[i] = i + 1;
    }

    for (i = 0; i < sqrtn && index < n; i++)
        for (j = 0; j < sqrtn && index < n; j++) {
            locations->x[index] = (grid[i] - 0.5 + uniform_distribution(-0.4, 0.4)) / sqrtn;
            locations->y[index] = (grid[j] - 0.5 + uniform_distribution(-0.4, 0.4)) / sqrtn;
            locations->z[index] = 1.0;
            index++;
        }


    for (j = 1; j < t_slots; j++) {
        for (i = 0; i < n; i++) {
            locations->x[i + j * n] = locations->x[i];
            locations->y[i + j * n] = locations->y[i];
            locations->z[i + j * n] = (double) (j + 1);
        }
    }
    return locations;
}

void print_dmatrix(char *desc, int m, int n, double* a, int lda)
//! print matrix contents (only for testing accuracy should be removed in release)
{
    int i, j;
    fprintf(stderr, "\n %s\n", desc);
    for (i = 0; i < m; i++) {
        for (j = 0; j < n; j++)
            fprintf(stderr, " %6.4e", a[i + lda * j]);
        fprintf(stderr, "\n");
    }
}


void print_diagonal(char *desc, int m, double* a, int lda)
//! print matrix contents (only for testing accuracy should be removed in release)
{
    int i, j;
    fprintf(stderr, "\n %s\n", desc);
    for (i = 0; i < m; i++) {
        fprintf(stderr, " %6.4e - ", a[i + lda * i]);
    }
    fprintf(stderr, "\n");
}


void print_smatrix(char *desc, int m, int n, float *a, int lda)
//! print matrix contents (only for testing accuracy should be removed in release)
{
    int i, j;
    fprintf(stderr, "\n %s\n", desc);
    for (i = 0; i < m; i++) {
        for (j = 0; j < n; j++)
            fprintf(stderr, " %6.4e", a[i + lda * j]);
        fprintf(stderr, "\n");
    }
}


int countlines(char *filename)
//!count the number of samples in real running mode
{
    FILE *fp = fopen(filename, "r");
    int ch = 0;
    int lines = 0;

    if (fp == NULL) {
        fprintf(stderr, "cannot open locations file\n");
        return 0;
    }

    while (!feof(fp)) {
        ch = fgetc(fp);
        if (ch == '\n')
            lines++;
    }

    fclose(fp);

    //Excluding header line
    return (lines);
}


void write_to_file(char *path, int matrix_size, int ncores, int tile_size, int test, char *computation, int async,
                   char *obsFPath, double total_exec_time, double avg_exec_time_per_iter, double avg_flops_per_iter,
                   int p_grid, int q_grid, double final_loglik, int n)
//! write results in  detail (only for testing accuracy should be removed in release)
{
    FILE *pFile;
    double peakperfomrance;
    double percent;
    pFile = fopen(path, "a");
    if (pFile == NULL) {
        fprintf(stderr, "Cannot access the results path(1)\n");
        exit(0);
    }
    peakperfomrance = q_grid * p_grid * ncores * 16 * 2.3;
    percent = avg_flops_per_iter / peakperfomrance;

    fprintf(pFile, "%d\t", n);
    fprintf(pFile, "%d\t", ncores);
    fprintf(pFile, "%d\t", q_grid * p_grid);
    fprintf(pFile, "%f\t", total_exec_time);
    fprintf(pFile, "%f\t", avg_exec_time_per_iter);
    fprintf(pFile, "%f\t\t", avg_flops_per_iter);
    fprintf(pFile, "%f\t", peakperfomrance);
    fprintf(pFile, "%f\t", percent);
    fprintf(pFile, "%d-", tile_size);
    fprintf(pFile, "%s-", computation);
    if (async == 0)
        fprintf(pFile, "SYNC-");
    else
        fprintf(pFile, "ASYNC-");
    if (test == 1)
        fprintf(pFile, "%d-", matrix_size);
    else
        fprintf(pFile, "%s-", obsFPath);
    fprintf(pFile, "%f \n", final_loglik);
    fclose(pFile);
}


void theta_parser2(double* theta_vec, char *kern, int num_params)
//! parse the theta vector, example: "1:0.5:0.1" -> {1, 0.5, 0.1}
{
    int i = 0;
    if (!strcmp(kern, "")) {
        for (i = 0; i < num_params; i++)
            theta_vec[i] = -1;
    }

    char *token;
    while ((token = strsep(&kern, ":")) != NULL) {
        if (strcmp(token, "?"))
            theta_vec[i] = (double) strtod(token, NULL);
        else
            theta_vec[i] = -1;
        i++;
    }
}

void write_pred_vector(double* zvec, MLE_data *data, int n, double mserror) {
    /*!
     * Returns initial_theta, starting_theta, target_theta.
     * @param[in] zvec: measurements vector.
     * @param[in] data: MLE_data struct with different MLE inputs.
     * @param[in] n: number of spatial locations
     * */

    FILE *pFile;
    pFile = fopen("predict_values.txt", "w+");
    if (pFile == NULL) {
        fprintf(stderr, "Cannot access the results path(1)\n");
        exit(0);
    }

    fprintf(pFile, "MSPE: %2.6f\n", mserror);
    int i = 0;
    for (i = 0; i < n; i++) {
        //	printf("zvec[i]: %f\n", zvec[i]);
        fprintf(pFile, "%f\n", zvec[i]);
    }
    fclose(pFile);
}

void write_vectors(double* zvec, MLE_data *data, int n)
//! store locations, measurements, and log files if log=1
/*!
 * Returns initial_theta, starting_theta, target_theta.
 * @param[in] zvec: measurements vector.
 * @param[in] data: MLE_data struct with different MLE inputs.
 * @param[in] n: number of spatial locations
 * */
{

    int i = 1;
    FILE *pFileZ, *pFileZ2, *pFileZ3, *pFileXY;
    location *l = &data->l1;
    struct stat st = {0};
    char *nFileZ = (char *) malloc(100 * sizeof(char));
    char *nFileZ2 = (char *) malloc(100 * sizeof(char));
    char *nFileZ3 = (char *) malloc(100 * sizeof(char));
    char *temp = (char *) malloc(100 * sizeof(char));
    char *nFileXY = (char *) malloc(100 * sizeof(char));
    data->nFileLog = (char *) malloc(100 * sizeof(char));
    int p = 1;
    //Create New directory if not exist
    if (stat("./synthetic_ds", &st) == -1)
        mkdir("./synthetic_ds", 0700);
    if (strcmp(data->kernel_fun, "bivariate_matern_parsimonious") == 0 ||
        strcmp(data->kernel_fun, "bivariate_matern_parsimonious_profile") == 0)
        p = 2;
    else if (strcmp(data->kernel_fun, "bivariate_matern_parsimonious2") == 0 ||
             strcmp(data->kernel_fun, "bivariate_matern_parsimonious2_profile") == 0)
        p = 2;
    else if (strcmp(data->kernel_fun, "univariate_spacetime_matern_stationary") == 0)
        p = 1;
    else if (strcmp(data->kernel_fun, "bivariate_spacetime_matern_stationary") == 0)
        p = 2;
    else if (strcmp(data->kernel_fun, "trivariate_matern_parsimonious") == 0 ||
             strcmp(data->kernel_fun, "trivariate_matern_parsimonious_profile") == 0)
        p = 3;
    snprintf(nFileZ, 100, "%s%d%s%s%s", "./synthetic_ds/Z1_", n / p, "_", data->kernel_fun, "_");
    snprintf(nFileZ2, 100, "%s%d%s%s%s", "./synthetic_ds/Z2_", n / p, "_", data->kernel_fun, "_");
    snprintf(nFileZ3, 100, "%s%d%s%s%s", "./synthetic_ds/Z3_", n / p, "_", data->kernel_fun, "_");
    snprintf(nFileXY, 100, "%s%d%s%s%s", "./synthetic_ds/LOC_", n / p, "_", data->kernel_fun, "_");
    snprintf(data->nFileLog, 100, "%s%d%s%s%s", "./synthetic_ds/log_", n / p, "_", data->kernel_fun, "_");

    snprintf(temp, 100, "%s%d", data->nFileLog, i);
    while (doesFileExist(temp) == 1) {
        i++;
        snprintf(temp, 100, "%s%d", data->nFileLog, i);
    }

    sprintf(temp, "%d", i);
    strcat(nFileZ, temp);
    strcat(nFileXY, temp);
    strcat(nFileZ2, temp);
    strcat(nFileZ3, temp);
    strcat(data->nFileLog, temp);

    pFileZ = fopen(nFileZ, "w+");
    pFileXY = fopen(nFileXY, "w+");

    if (p == 1)
        for (i = 0; i < n; i++)
            fprintf(pFileZ, "%f\n", zvec[i]);
    else if (p == 2) {
        pFileZ2 = fopen(nFileZ2, "w+");
        for (i = 0; i < n; i += 2) {
            fprintf(pFileZ, "%f\n", zvec[i]);
            fprintf(pFileZ2, "%f\n", zvec[i + 1]);
        }
        fclose(pFileZ2);
    } else if (p == 3) {
        pFileZ2 = fopen(nFileZ2, "w+");
        pFileZ3 = fopen(nFileZ3, "w+");
        for (i = 0; i < n; i += 3) {
            fprintf(pFileZ, "%f\n", zvec[i]);
            fprintf(pFileZ2, "%f\n", zvec[i + 1]);
            fprintf(pFileZ3, "%f\n", zvec[i + 2]);
        }
        fclose(pFileZ2);
        fclose(pFileZ3);
    }
    for (i = 0; i < n / p; i++) {
        if (l->z == NULL)
            fprintf(pFileXY, "%f,%f\n", l->x[i], l->y[i]);
        else
            fprintf(pFileXY, "%f,%f,%f\n", l->x[i], l->y[i], l->z[i]);
    }

    fclose(pFileZ);
    fclose(pFileXY);
}

void write_to_thetafile(char *path, double* theta, int num_params,
                        int n, double time_per_iter,
                        int total_no_iters,
                        double prediction_error,
                        double mloe,
                        double mmom)
//! write results (only for testing accuracy should be removed in release)
{
    FILE *pFile;
    pFile = fopen(path, "a");

    if (pFile == NULL) {
        printf("Cannot access the results path(3)\n");
        exit(0);
    }
    for (int i = 0; i < num_params; i++) {
        fprintf(pFile, "%f,", theta[i]);
    }

    fprintf(pFile, "%d,", n);
    fprintf(pFile, "%d,", total_no_iters);
    fprintf(pFile, "%f,", time_per_iter);
    fprintf(pFile, "%f,", prediction_error);
    fprintf(pFile, "%f,", mloe);
    fprintf(pFile, "%f\n", mmom);
    fclose(pFile);

}


void write_to_estimatedtheta(char *path, double* theta, int num_params,
                             int n, double prediction_time,
                             double mloe_mmom_time,
                             double prediction_error1,
                             double prediction_error2,
                             double prediction_error3,
                             double prediction_error,
                             double mloe,
                             double mmom,
                             int zvecs)
//! write results (only for testing accuracy should be removed in release)
{
    FILE *pFile;
#if defined(CHAMELEON_USE_MPI)
    if ( MORSE_My_Mpi_Rank() == 0 )
    {
#endif
    pFile = fopen(path, "a");

    if (pFile == NULL) {
        printf("Cannot access the results path(3)\n");
        exit(0);
    }

    fprintf(pFile, "%d, ", n);
    for (int i = 0; i < num_params; i++) {
        fprintf(pFile, "%f ", theta[i]);
    }

    fprintf(pFile, "  ---- %f ", prediction_time);
    fprintf(pFile, "%g ", prediction_error1);
    fprintf(pFile, "%g ", prediction_error2);
    fprintf(pFile, "%g ", prediction_error3);
    fprintf(pFile, "%g --- ", prediction_error);
    fprintf(pFile, "%f ", mloe_mmom_time);
    fprintf(pFile, "%f ", mloe);
    fprintf(pFile, "%f ", mmom);
    fprintf(pFile, "%d\n", zvecs);
    fclose(pFile);
#if defined(CHAMELEON_USE_MPI)
    }
#endif

}

void shuffle(double* array, location *locations, size_t n)
//! shuffle an array
{
    if (n > 1) {
        size_t i;

        for (i = 0; i < n - 1; i++) {
            size_t j = i + rand() / (RAND_MAX / (n - i) + 1);

            double t = array[j];
            array[j] = array[i];
            array[i] = t;
            double xtemp = locations->x[j];
            locations->x[j] = locations->x[i];
            locations->x[i] = xtemp;
            double ytemp = locations->y[j];
            locations->y[j] = locations->y[i];
            locations->y[i] = ytemp;

        }
    }
}

void shuffle2(double* array, double* array2, location *locations, size_t n)
//! shuffle an array
{
    if (n > 1) {
        size_t i;

        for (i = 0; i < n - 1; i++) {
            size_t j = i + rand() / (RAND_MAX / (n - i) + 1);

            double t = array[j];
            array[j] = array[i];
            array[i] = t;
            double t2 = array2[j];
            array2[j] = array2[i];
            array2[i] = t2;
            double xtemp = locations->x[j];
            locations->x[j] = locations->x[i];
            locations->x[i] = xtemp;
            double ytemp = locations->y[j];
            locations->y[j] = locations->y[i];
            locations->y[i] = ytemp;

        }
    }
}


void shuffle3(double* array, double* array2, double* array3, location *locations, size_t n)
//! shuffle an array
{
    if (n > 1) {
        size_t i;
        //    unsigned int *seed = &exageostat_seed;
        //    size_t jj = i + rand() / (RAND_MAX / (n - i) + 1);
        for (i = 0; i < n - 1; i++) {
            size_t j = i + rand() / (RAND_MAX / (n - i) + 1);

            double t = array[j];
            array[j] = array[i];
            array[i] = t;
            double t2 = array2[j];
            array2[j] = array2[i];
            array2[i] = t2;
            double t3 = array3[j];
            array3[j] = array3[i];
            array3[i] = t3;
            double xtemp = locations->x[j];
            locations->x[j] = locations->x[i];
            locations->x[i] = xtemp;
            double ytemp = locations->y[j];
            locations->y[j] = locations->y[i];
            locations->y[i] = ytemp;

        }
    }
}

void theta_parser(double* initial_theta, double* target_theta, double* starting_theta,
                  char *ikernel, char *kernel, double* lb,
                  double* up, int test, int num_params)
//! Parse initial_theta, target_theta, starting_theta, inputs
/*!
 * Returns initial_theta, starting_theta, target_theta.
 * @param[out] initial_theta: initial_theta Vector with three parameter (Variance, Range, Smoothness)
 that is used to to generate the Covariance Matrix and initial Z vector.
 * @param[out] starting_theta: theta Vector with three parameter (Variance, Range, Smoothness)
 that is used to to generate the Covariance Matrix of the first MLE iteration.
 * @param[out] target_theta: target theta Vector with three parameter (Variance, Range, Smoothness) unknown theta parameter should be shown as '?'.
 * @param[in] ikernel: initial_theta Vector as string.
 * @param[in] kernel:  target_theta Vector as string.
 * @param[in] lb: optimization lower bounds vector ( lb_1, lb_2, lb_3).
 * @param[in] up: optimization upper bounds vector ( ub_1, ub_2, ub_3).
 * @param[in] test: if test=1 ->test running mode, test=0-> real running mode.
 * */
{

    int i = 0;

    theta_parser2(initial_theta, ikernel, num_params);
    theta_parser2(target_theta, kernel, num_params);

    for (i = 0; i < num_params; i++) {
        if (target_theta[i] != -1) {
            lb[i] = target_theta[i];
            up[i] = target_theta[i];
            starting_theta[i] = target_theta[i];
        }
    }

}

void init_optimizer(nlopt_opt *opt, double* lb, double* up, double tol)
//! Initialize the NLOPT optimizer
/*!
 * Returns nlopt_opt object.
 * @param[in] lb: optimization lower bounds vector ( lb_1, lb_2, lb_3).
 * @param[in] up: optimization upper bounds vector ( ub_1, ub_2, ub_3).
 * @param[in] tol: a tolerance that is used for the purpose of stopping criteria only.
 * @param[out] opt: nlopt_opt object.
 * */
{
    //initalizing opt library
    nlopt_set_lower_bounds(*opt, lb);
    nlopt_set_upper_bounds(*opt, up);
    nlopt_set_ftol_abs(*opt, tol);
}

void print_summary(int test, int N, int ncores, int gpus, int ts, int lts, char *computation, int zvecs, int p_grid,
                   int q_grid, int precision)
//! print the summary of MLE inputs.
{
#if defined(CHAMELEON_USE_MPI)
    if ( MORSE_My_Mpi_Rank() == 0 )
    {
#endif
    fprintf(stderr, "********************SUMMARY**********************\n");
    if (test == 0)
        fprintf(stderr, "#Execution Mode: Real Dataset\n");
    else
        fprintf(stderr, "#Synthetic Dataset\n");
    fprintf(stderr, "Number of Locations: %d\n", N);
    fprintf(stderr, "#Threads per node: %d\n", ncores);
    fprintf(stderr, "#GPUs: %d\n", gpus);
    if (precision == 0)
        fprintf(stderr, "#Double Precision!\n");
    else if (precision == 1)
        fprintf(stderr, "#Single Precision!\n");
    else if (precision == 2)
        fprintf(stderr, "#Single/Double Precision!\n");
    if (strcmp(computation, "exact") == 0 || strcmp(computation, "diag_approx") == 0)
        fprintf(stderr, "#Dense Tile Size: %d\n", ts);
    else
        fprintf(stderr, "LR Tile Size: %d\n", lts);
    fprintf(stderr, "#%s computation\n", computation);
    fprintf(stderr, "#Obervation Vectors (Z): %d\n", zvecs);
    fprintf(stderr, "p=%d, q=%d\n", p_grid, q_grid);
    fprintf(stderr, "***************************************************\n");
#if defined(CHAMELEON_USE_MPI)
    }
#endif
}


int print_result(MLE_data *data, double* starting_theta, int N, int zvecs, int ncores, int ts, int test,
                 double* initial_theta, char *computation, int p_grid, int q_grid, double final_loglik,
                 double prediction_error)
//! print results (only for testing accuracy should be removed in release)
{
    int num_params,split;
#if defined(CHAMELEON_USE_MPI)
    if ( MORSE_My_Mpi_Rank() == 0 )
    {
#endif

    if (strcmp(data->kernel_fun, "univariate_matern_stationary") == 0 ||
        strcmp(data->kernel_fun, "univariate_pow_exp_stationary") == 0)
        num_params = 3;
    else if (strcmp(data->kernel_fun, "univariate_matern_nuggets_stationary") == 0)
        num_params = 4;
    else if (strcmp(data->kernel_fun, "univariate_matern_non_stationary") == 0)
        num_params = 9;
    else if (strcmp(data->kernel_fun, "bivariate_matern_flexible") == 0)
        num_params = 11;
    else if (strcmp(data->kernel_fun, "bivariate_matern_parsimonious") == 0)
        num_params = 6;
    else if (strcmp(data->kernel_fun, "bivariate_matern_parsimonious2") == 0)
        num_params = 6;
    else if (strcmp(data->kernel_fun, "bivariate_matern_parsimonious_profile") == 0)
        num_params = 6;
    else if (strcmp(data->kernel_fun, "bivariate_matern_parsimonious2_profile") == 0)
        num_params = 6;
    else if (strcmp(data->kernel_fun, "univariate_spacetime_matern_stationary") == 0)
        num_params = 7;
    else if (strcmp(data->kernel_fun, "bivariate_spacetime_matern_stationary") == 0)
        num_params = 10;
    else if (strcmp(data->kernel_fun, "trivariate_matern_parsimonious") == 0 ||
             strcmp(data->kernel_fun, "trivariate_matern_parsimonious_profile") == 0)
        num_params = 10;
    else if (strcmp(data->kernel_fun, "univariate_matern_non_stat_estimation") == 0){
        split = data->Nsplit;
        num_params = (int)3*split;
    }
    else if (strcmp(data->kernel_fun, "univariate_matern_non_gaussian") == 0 ||
               strcmp(data->kernel_fun, "univariate_exp_non_gaussian") == 0) 
        num_params = 6;
    else {
        fprintf(stderr, "Choosen kernel is not exist(7)!\n");
        fprintf(stderr, "Called function is: %s\n", __func__);
        exit(0);
    }
    //num_params = strcmp(data->kernel_fun, "stationary_kernel")   == 0? 3 : 9;
    printf("********************SUMMARY*****************************\n");
    printf("#Total Number of Iterations=%d\n", data->iter_count);
    printf("#Total Optimization Time= %6.2f\n", data->total_exec_time);
    printf("#Found Maximum at (");
    int i = 0;
    if (strcmp(data->kernel_fun, "bivariate_matern_parsimonious_profile") == 0 ||
        strcmp(data->kernel_fun, "bivariate_matern_parsimonious2_profile") == 0) {
        printf("%2.6f, %2.6f,", data->variance1, data->variance2);
        i = 2;
    } else if (strcmp(data->kernel_fun, "trivariate_matern_parsimonious_profile") == 0) {
        printf("%2.6f, %2.6f, %2.6f,", data->variance1, data->variance2, data->variance3);
        i = 3;
    } else
        i = 0;

    for (; i < num_params - 1; i++)
        printf("%0.8f, ", starting_theta[i]);
    printf("%0.8f", starting_theta[i]);
    printf(") \n");

    printf("#Number of Locations: %d\n", results.problem_size);
    printf("#Computation: %s\n", results.computation);
    printf("#Kernel: %s\n", results.kernel);
    printf("#Dataset: %s\n", results.ds_type);
    printf("#Precision: %s\n", results.precision);
    printf("#Sample: %d\n", results.z_sample);
    printf("#Dense tile size: %d\n", results.dense_ts);
    printf("#LR tile size (if used): %d\n", results.lr_ts);
    printf("#LR accuracy/max_rank (if used): %d/%d\n", results.lr_acc, results.lr_maxrank);
    printf("#CPU cores per node: %d\n", results.ncores);
    printf("#GPUs per node: %d\n", results.ngpus);
    printf("#p_grid: %d, q_grid: %d\n", results.p, results.q);
    printf("***************************************************\n\n");

    if (data->log == 1) {

        fprintf(data->pFileLog, "Total Number of Iterations=%d\n", data->iter_count);
        fprintf(data->pFileLog, "Total Optimization Time= %6.2f secs\n", data->total_exec_time);
        fprintf(data->pFileLog, "Found Maximum at (");
        if (strcmp(data->kernel_fun, "bivariate_matern_parsimonious_profile") == 0 ||
            strcmp(data->kernel_fun, "bivariate_matern_parsimonious2_profile") == 0) {
            printf("%2.6f, %2.6f,", data->variance1, data->variance2);
            i = 2;
        } else if (strcmp(data->kernel_fun, "trivariate_matern_parsimonious_profile") == 0) {
            printf("%2.6f, %2.6f, %2.6f,", data->variance1, data->variance2, data->variance3);
            i = 3;
        } else
            i = 0;

        for (; i < num_params; i++)
            fprintf(data->pFileLog, "%g, ", starting_theta[i]);
        fprintf(stderr, ")\n");
    }

    results.time_per_iteration = data->avg_exec_time_per_iter / data->iter_count;
    results.flops_per_iteration = data->avg_flops_per_iter / data->iter_count;
    results.total_mle_time = data->avg_exec_time_per_iter;

    FILE *pFile;
    pFile = fopen("results.log", "a");

    if (pFile == NULL) {
        printf("Cannot access the results path(3)\n");
        return -1;
    }

    fprintf(pFile, "%d ", results.problem_size);
    fprintf(pFile, "%s(", results.computation);
    fprintf(pFile, "%d) ", results.lr_acc);

    for (int i = 0; i < results.num_params; i++)
        fprintf(pFile, "%6.6f ", results.initial_theta[i]);

    for (int i = 0; i < results.num_params; i++)
        fprintf(pFile, "%6.6f ", results.estimated_theta[i]);

    fprintf(pFile, "%6.6f ", results.final_loglik);
    fprintf(pFile, "%6.6f ", results.time_per_iteration);
    fprintf(pFile, "%6.6f ", results.flops_per_iteration);
    fprintf(pFile, "%6.6f ", results.total_mle_time);

    fprintf(pFile, "%6.6f ", results.mse_pred1);
    fprintf(pFile, "%6.6f ", results.mse_pred2);
    fprintf(pFile, "%6.6f ", results.mse_pred);
    fprintf(pFile, "%6.4e ", results.trace_pred_sum);
    fprintf(pFile, "%6.4e ", results.trace_pred_mean);
    fprintf(pFile, "%6.6f ", results.total_pred_time);
    fprintf(pFile, "%6.6f ", results.total_pred_flops);

    fprintf(pFile, "%6.6f ", results.mloe);
    fprintf(pFile, "%6.6f ", results.mmom);
    fprintf(pFile, "%s ", results.mloe_exec);
    fprintf(pFile, "%6.2f ", results.total_mloe_mmom_time);
    fprintf(pFile, "%6.2f ", results.matrix_gen_mloe_mmom_time);
    fprintf(pFile, "%6.2f ", results.cho_fact_mloe_mmom_time);
    fprintf(pFile, "%6.2f ", results.loop_mloe_mmom_time);
    fprintf(pFile, "%6.2f ", results.total_mloe_mmom_flops);
    fprintf(pFile, "%6.8f ", results.fisher_00);
    fprintf(pFile, "%6.8f ", results.fisher_11);
    fprintf(pFile, "%6.8f ", results.fisher_22);
    fprintf(pFile, "%6.2f\n\n", results.fisher_time);

    fclose(pFile);
#if defined(CHAMELEON_USE_MPI)
    }
#endif
}


double cWtime(void)
//! get time
{
    struct timeval tp;
    gettimeofday(&tp, NULL);
    return tp.tv_sec + 1e-6 * tp.tv_usec;
}

//! Read real location file
/*!
 * Returns locations l1.
 * @param[in] loc_file: locations file  path.
 * @param[in]  n : number of spatial locations (observations).
 * @param[out] l1 :  location struct with (x,y) locations.
 */

void
write_prediction_result(char *path, int matrix_size, int no_missing, double MSE1, double MSE2, double MSE3, double MSE,
                        double solve_time, double flops)
//! write prediction results (only for testing accuracy should be removed in release).
{
    FILE *pFile;

    pFile = fopen(path, "a");

    if (pFile == NULL) {

        fprintf(stderr, "Cannot access the results path(2)\n");
        exit(0);
    }

    fprintf(pFile, "%d\t", matrix_size);
    fprintf(pFile, "%d\t", no_missing);
    fprintf(pFile, "%f\t", MSE1);
    fprintf(pFile, "%f\t", MSE2);
    fprintf(pFile, "%f\t", MSE3);
    fprintf(pFile, "%f\t", MSE);
    fprintf(pFile, "%6.2f secs\t", solve_time);
    fprintf(pFile, "%6.2f Gflops\n", flops);

    fclose(pFile);
}

// check if the file exist or not
int doesFileExist(const char *filename)
/*! variables from argument
 * Returns 0 if file not exist  and 1 if file exist.
 * @param[in] filename: file path.
 */
{
    struct stat st;
    int result = stat(filename, &st);
    return result == 0;
}

void init_log(MLE_data *data)
//! init and open log files if log==1
/*!
 * Returns MLE_data struct with new log files.
 * @param[out] data: MLE_data struct with different MLE inputs.
 */
{
#if defined(CHAMELEON_USE_MPI)
    if ( MORSE_My_Mpi_Rank() == 0 )
    {
#endif
    data->pFileLog = fopen(data->nFileLog, "w+");
    fprintf(data->pFileLog, "\t\tlog file is generated by ExaGeoStat application\n");
    fprintf(data->pFileLog, "\t\t============================================\n");
#if defined(CHAMELEON_USE_MPI)
    }
#endif
}

void finalize_log(MLE_data *data)
//! finalize and close log files if log==1
/*!
 * Returns MLE_data struct with new log files.
 * @param[out] data: MLE_data struct with different MLE inputs.
 */
{
#if defined(CHAMELEON_USE_MPI)
    if ( MORSE_My_Mpi_Rank() == 0 )
    {
#endif
    fclose(data->pFileLog);
#if defined(CHAMELEON_USE_MPI)
    }
#endif
}

void split_data(MLE_data *data, location *locations, double* Z, double* Zactual, int *N, int nZmiss) {

    int i = 0;
    int nZobs = *N - nZmiss;
    shuffle(Z, locations, *N);
    data->l1.x = (double* ) malloc(nZobs * sizeof(double));
    data->l1.y = (double* ) malloc(nZobs * sizeof(double));
    data->lmiss.x = (double* ) malloc(nZmiss * sizeof(double));
    data->lmiss.y = (double* ) malloc(nZmiss * sizeof(double));

    for (i = 0; i < nZobs; i++) {
        data->l1.x[i] = locations->x[i];
        data->l1.y[i] = locations->y[i];
    }

    data->lobs.x = data->l1.x;
    data->lobs.y = data->l1.y;

    for (i = 0; i < nZmiss; i++) {
        data->lmiss.x[i] = locations->x[i + nZobs];
        data->lmiss.y[i] = locations->y[i + nZobs];
        Zactual[i] = Z[i + nZobs];
    }
    *N = nZobs;

}

void pick_random_points_noshuffle(MLE_data *data, double* Zobs, double* Zactual, int nZmiss, int nZobs, int N) {
    //initialization
    location l;
    location *lmiss;
    location *lobs;
    double* Z;
    int i = 0;

    //memory allocation
    Z = (double* ) malloc(N * sizeof(double));
    l.x = (double* ) malloc(N / 2 * sizeof(double));
    l.y = (double* ) malloc(N / 2 * sizeof(double));

    //copy observed measurments
    MLE_get_zobs(data, Z, N);

    for (i = 0; i < N / 2; i++) {
        l.x[i] = data->l1.x[i];
        l.y[i] = data->l1.y[i];
    }

    for (i = 0; i < nZobs; i++)
        Zobs[i] = Z[nZmiss + i];

    for (i = 0; i < nZmiss; i++)
        Zactual[i] = Z[i];

    lmiss = &(data->lmiss);
    lobs = &(data->lobs);
    lmiss->x = l.x;
    lmiss->y = l.y;
    lobs->x = &l.x[nZmiss];
    lobs->y = &l.y[nZmiss];

    locations_obs_zsort_inplace(nZobs, lobs, Zobs);
    locations_obs_zsort_inplace(nZmiss, lmiss, Zactual);

    free(Z);
}

void pick_random_points(MLE_data *data, double* Zobs, double* Zactual, int nZmiss, int nZobs, int N) {
    //initialization
    location l;
    location *lmiss;
    location *lobs;


    double* Z;
    int i = 0;

    //memory allocation
    Z = (double* ) malloc(N * sizeof(double));
    l.x = (double* ) malloc(N * sizeof(double));
    l.y = (double* ) malloc(N * sizeof(double));

    //copy observed measurments
    MLE_get_zobs(data, Z, N);

    for (i = 0; i < N; i++) {
        l.x[i] = data->l1.x[i];
        l.y[i] = data->l1.y[i];
    }

    shuffle(Z, &l, N);
    for (i = 0; i < nZobs; i++)
        Zobs[i] = Z[nZmiss + i];

    for (i = 0; i < nZmiss; i++)
        Zactual[i] = Z[i];

    lmiss = &(data->lmiss);
    lobs = &(data->lobs);
    lmiss->x = l.x;
    lmiss->y = l.y;
    lobs->x = &l.x[nZmiss];
    lobs->y = &l.y[nZmiss];

    locations_obs_zsort_inplace(nZobs, lobs, Zobs);
    locations_obs_zsort_inplace(nZmiss, lmiss, Zactual);

    free(Z);
}


void pick_random_points2(MLE_data *data, double* Zobs, double* Zactual, int nZmiss, int nZobs, int N) {
    //initialization
    location l;
    location *lmiss;
    location *lobs;
    double* Z;
    double* Z1;
    double* Z2;
    int i = 0;
    int p = 2;

    //memory allocation
    Z = (double* ) malloc(N * sizeof(double));
    l.x = (double* ) malloc(N / p * sizeof(double));
    l.y = (double* ) malloc(N / p * sizeof(double));
    Z1 = (double* ) malloc(N * sizeof(double));
    Z2 = (double* ) malloc(N * sizeof(double));

    //copy observed measurments
    MLE_get_zobs(data, Z, N);

    double sum = 0;
    int j = 0;
    for (i = 0; i < N; i += 2) {
        Z1[j] = Z[i];
        Z2[j++] = Z[i + 1];
    }


    for (i = 0; i < N / p; i++) {
        l.x[i] = data->l1.x[i];
        l.y[i] = data->l1.y[i];
    }

    shuffle2(Z1, Z2, &l, N / p);

    //actual vector
    j = 0;
    for (i = 0; i < nZmiss; i++) {
        Zactual[j] = Z1[i];
        Zactual[j + 1] = Z2[i];
        printf("%f, %f, ---%f, %f\n", l.x[i], l.y[i], Z1[i], Z2[i]);
        j += 2;
    }

    //observation vector
    j = 0;
    for (i = 0; i < nZobs; i++) {
        Zobs[j] = Z1[nZmiss + i];
        Zobs[j + 1] = Z2[nZmiss + i];
        j += 2;
    }

    lmiss = &(data->lmiss);
    lobs = &(data->lobs);
    lmiss->x = l.x;
    lmiss->y = l.y;
    lobs->x = &l.x[nZmiss];
    lobs->y = &l.y[nZmiss];

    //TODO: check why using them give wrong answer
    //	locations_obs_zsort_inplace(nZobs, lobs, Zobs);
    //	locations_obs_zsort_inplace(nZmiss, lmiss, Zactual);

    free(Z1);
    free(Z2);
}


void pick_random_points3(MLE_data *data, double* Zobs, double* Zactual, int nZmiss, int nZobs, int N) {
    //initialization
    location l;
    location *lmiss;
    location *lobs;
    double* Z;
    double* Z1;
    double* Z2;
    double* Z3;
    int i = 0;
    int p = 3;

    //memory allocation
    Z = (double* ) malloc(N * sizeof(double));
    l.x = (double* ) malloc(N / p * sizeof(double));
    l.y = (double* ) malloc(N / p * sizeof(double));
    Z1 = (double* ) malloc(N * sizeof(double));
    Z2 = (double* ) malloc(N * sizeof(double));
    Z3 = (double* ) malloc(N * sizeof(double));
    //copy observed measurments
    MLE_get_zobs(data, Z, N);

    double sum = 0;
    int j = 0;
    for (i = 0; i < N; i += 3) {
        Z1[j] = Z[i];
        Z2[j++] = Z[i + 1];
        Z3[j++] = Z[i + 2];
    }

    for (i = 0; i < N / p; i++) {
        l.x[i] = data->l1.x[i];
        l.y[i] = data->l1.y[i];
    }

    shuffle3(Z1, Z2, Z3, &l, N / p);

    //actual vector
    j = 0;
    for (i = 0; i < nZmiss; i++) {
        Zactual[j] = Z1[i];
        Zactual[j + 1] = Z2[i];
        Zactual[j + 2] = Z3[i];
        j += 3;
    }

    //observation vector
    j = 0;
    for (i = 0; i < nZobs; i++) {
        Zobs[j] = Z1[nZmiss + i];
        Zobs[j + 1] = Z2[nZmiss + i];
        Zobs[j + 2] = Z3[nZmiss + i];
        j += 3;
    }

    lmiss = &(data->lmiss);
    lobs = &(data->lobs);
    lmiss->x = l.x;
    lmiss->y = l.y;
    lobs->x = &l.x[nZmiss];
    lobs->y = &l.y[nZmiss];

    //TODO: check why using them give wrong answer
    //  locations_obs_zsort_inplace(nZobs, lobs, Zobs);
    //  locations_obs_zsort_inplace(nZmiss, lmiss, Zactual);

    free(Z1);
    free(Z2);
    free(Z3);
}

void generate_interior_points(MLE_data *data, double* Zobs, double* Zactual, int nZmiss, int nZobs, int N) {
    //initialization
    location l;
    location *lmiss;
    location *lobs;
    double* Z;
    int i = 0;

    if (nZmiss >= nZobs) {
        fprintf(stderr, "Cannot generate missing locations larger than or equal the observed ones\n");
        return;
    }

    //memory allocation
    Z = (double* ) malloc(N * sizeof(double));
    l.x = (double* ) malloc(N * sizeof(double));
    l.y = (double* ) malloc(N * sizeof(double));

    //copy observed measurments
    MLE_get_zobs(data, Z, N);

    for (i = 0; i < N; i++) {
        l.x[i] = data->l1.x[i];
        l.y[i] = data->l1.y[i];
    }

    shuffle(Z, &l, N);

    for (i = 0; i < nZobs; i++)
        Zobs[i] = Z[nZmiss + i];

    lmiss = &(data->lmiss);
    lobs = &(data->lobs);

    for (i = 0; i < nZmiss; i++) {
        l.x[i] = (l.x[i] + l.x[i + 1]) / 2.0;
        l.y[i] = (l.y[i] + l.y[i + 1]) / 2.0;
    }

    lmiss->x = l.x;
    lmiss->y = l.y;
    lobs->x = &l.x[nZmiss];
    lobs->y = &l.y[nZmiss];
    free(Z);
}

// This function converts decimal degrees to radians
static double deg2rad(double deg) {
    return (deg * PI / 180);
}

//  This function converts radians to decimal degrees
static double rad2deg(double rad) {
    return (rad * 180 / PI);
}

/**
 * Returns the distance between two points on the Earth.
 * Direct translation from http://en.wikipedia.org/wiki/Haversine_formula
 * @param lat1d Latitude of the first point in degrees
 * @param lon1d Longitude of the first point in degrees
 * @param lat2d Latitude of the second point in degrees
 * @param lon2d Longitude of the second point in degrees
 * @return The distance between the two points in kilometers
 */
static double distanceEarth(double lat1d, double lon1d, double lat2d, double lon2d) {
    double lat1r, lon1r, lat2r, lon2r, u, v;
    lat1r = deg2rad(lat1d);
    lon1r = deg2rad(lon1d);
    lat2r = deg2rad(lat2d);
    lon2r = deg2rad(lon2d);
    u = sin((lat2r - lat1r) / 2);
    v = sin((lon2r - lon1r) / 2);
    return 2.0 * earthRadiusKm * asin(sqrt(u * u + cos(lat1r) * cos(lat2r) * v * v));
}

static double calculateDistance(double x1, double y1, double x2, double y2, int distance_metric) {

    if (distance_metric == 1)
        return distanceEarth(x1, y1, x2, y2);
    return sqrt(pow((x2 - x1), 2) + pow((y2 - y1), 2));
}


double core_matern_vector(double x0, double y0, double x1, double y1, double* localtheta, int distance_metric) {

    int i, j;
    double expr = 0.0;
    double con = 0.0;
    double sigma_square = localtheta[0];// * localtheta[0];

    con = pow(2, (localtheta[2] - 1)) * tgamma(localtheta[2]);
    con = 1.0 / con;
    con = sigma_square * con;

    expr = calculateDistance(x0, y0, x1, y1, distance_metric) / localtheta[1];
    if (expr == 0)
        return sigma_square /*+ 1e-4*/;
    else
        return con * pow(expr, localtheta[2]) * gsl_sf_bessel_Knu(localtheta[2], expr); // Matern Function
}

void fwrite_array(int m, int n, int ld, double* arr, char *file) {
    FILE *fp = fopen(file, "w");
    if (fp == NULL) {
        fprintf(stderr, "File %s cannot be opened to write\n", file);
        exit(1);
    }
    int i, j;
    fprintf(fp, "%d %d\n", m, n);
    for (i = 0; i < m; i++) {
        for (j = 0; j < n; j++) {
            fprintf(fp, "%d\t", (int) arr[ld * j + i]);
        }
        fprintf(fp, "\n");
    }
    fclose(fp);
}

void print_array(int m, int n, int ld, double* arr, FILE *fp) {
    int i, j;
    fprintf(fp, "%d %d\n", m, n);
    for (i = 0; i < m; i++) {
        for (j = 0; j < n; j++) {
            fprintf(fp, "%d\t", (int) arr[ld * j + i]);
        }
        fprintf(fp, "\n");
    }
}

void checkpointing(char *path, int iter_count, double* theta, double loglik, int num_params)
//! write prediction results (only for testing accuracy should be removed in release).
{
    FILE *pFile;

    pFile = fopen(path, "a");

    if (pFile == NULL) {

        fprintf(stderr, "Cannot access the results path\n");
        exit(0);
    }

    fprintf(pFile, "%d ", iter_count);
    for (int i = 0; i < num_params; i++) {
        fprintf(pFile, "%.17g ", theta[i]);
    }
    fprintf(pFile, "%.17g ", loglik);

    fprintf(pFile, "\n");
    fclose(pFile);
}

bool recover(char *path, int iter_count, double* theta, double* loglik, int num_params) {

    FILE *fp;
    char *line = NULL;
    size_t len = 0;
    ssize_t read;
    int count = 0;
    int i = 0;
    char *pch;
    fp = fopen(path, "r");
    if (fp == NULL) {
        printf("cannot open observations file\n");
        exit(EXIT_FAILURE);
    }

    while ((read = getline(&line, &len, fp)) != -1) {
        pch = strtok(line, " ");
        count = atoi(pch);
        if (count == iter_count) {
            pch = strtok(NULL, " ");
            for (i = 0; i < num_params; i++) {
                theta[i] = atof(pch);
                pch = strtok(NULL, " ");
            }
            *loglik = atof(pch);
            fclose(fp);
            free(line);

            return true;
        }
        count++;
    }

    fclose(fp);
    free(line);

    return false;
}

double* pred_idw(MLE_data *data, double* z_miss, double* z_actual, double* z_obs, int nZmiss, int nZobs) {
    int i = 0, j = 0;
    double sigma_1 = 0;
    double sigma_2 = 0;
    double error = 0, dij = 0;
    double error1 = 0, error2 = 0;
    double x1, y1, x2, y2;
    int n = 2;
    location *lmiss = &data->lmiss;
    location *lobs = &data->lobs;
    double* mspe = (double* ) malloc(3 * sizeof(double));
    for (j = 0; j < nZmiss; j++) {
        x2 = lmiss->x[j];
        y2 = lmiss->y[j];
        for (i = 0; i < nZobs; i++) {
            x1 = lobs->x[i];
            y1 = lobs->y[i];
            dij = calculateDistance(x1, y1, x2, y2, 0);
            if (dij != 0) {
                sigma_1 += z_obs[i] / pow(dij, n);
                sigma_2 += 1.0 / pow(dij, n);
            }
        }
        z_miss[j] = sigma_1 / sigma_2;
        if (j % 2 == 0)
            error1 += pow((z_actual[j] - z_miss[j]), 2);
        else
            error2 += pow((z_actual[j] - z_miss[j]), 2);
        error += pow((z_actual[j] - z_miss[j]), 2);
        sigma_1 = 0;
        sigma_2 = 0;
    }

    mspe[0] = error / nZmiss;
    mspe[1] = error1 / (nZmiss / 2);
    mspe[2] = error2 / (nZmiss / 2);
    return mspe;
}

int print_predicted_values(location *lmiss, double* Zactual, double* Zmiss, int nZmiss, int p)
//! print results (only for testing accuracy should be removed in release)

{
    int i;
    FILE *pFile;
#if defined(CHAMELEON_USE_MPI)
    if ( MORSE_My_Mpi_Rank() == 0 )
    {
#endif

    char buffer[100];

    time_t t = time(NULL);
    struct tm *pp = localtime(&t);

    strftime(buffer, 100, "predicted_values_%B_%d_%Y_%T", pp);

    printf("%s\n", buffer);

    pFile = fopen(buffer, "w");

    if (pFile == NULL) {
        printf("Cannot access the path\n");
        return -1;
    }

    fprintf(pFile, "(    x   ,     y   ) ----> var:i (Z_actual, Z_predict)\n");
    int j = 0, k = 0;
    for (int index = 0; index < nZmiss; index++) {
        fprintf(pFile, "(%3.6f, %3.6f) ----> ", lmiss->x[index], lmiss->y[index]);

        for (k = 0; k < p - 1; k++) {
            fprintf(pFile, "var:%d (%3.6f, %3.6f) - ", k, Zactual[j], Zmiss[j]);
            j++;
        }

        fprintf(pFile, "var:%d (%3.6f, %3.6f)\n", k, Zactual[j], Zmiss[j]);
        j++;
    }
    fclose(pFile);
#if defined(CHAMELEON_USE_MPI)
    }
#endif
    return 0;
}

int data_MinMax_nonstat(double *array, double *a, int N){
    int i;
    a[0]=a[1]=array[0];
    // size_t n = sizeof(array) / sizeof(array[0]);
    int n = N;
    // printf("%d \t %d\n",sizeof(array),sizeof(array[0]));
    // printf("x vector size is %f", array[3120]);
    for(i=1; i<n; i++)
    {
        if(a[0]>array[i])
		  a[0]=array[i];   
		if(a[1]<array[i])
		    a[1]=array[i];       
    }
    // printf("minimum is %f \n",a[0]);
    // printf("maximum is %f \n",a[1]);
    return(0);
}

int check_stationarity(double *x,double *y, double *z, double startx, double endx, double starty, double endy, double startz, double endz, bool z_user_range, int N){
    double probability;
    // printf("%f\t%f\t%f\t%f\t%f\t%f\t%f\n",startx,endx,starty,endy,x[1],y[1],z[1]);
    // if(sizeof(predict(x,y,z,startx,endx,starty,endy,N)) == sizeof(int) ) {
    //     printf("%f",predict(x,y,z,startx,endx,starty,endy,N));
    // }else{
    //     printf("no its not\n");
    // }

    // printf("------%f\n", predict(x,y,z,startx,endx,starty,endy,N));
    
    probability = predict(x,y,z,startx,endx,starty,endy,startz,endz, z_user_range,N)/10000.0;
    printf("\n\nprobability value in C : %f\n",probability);
    // exit(0);
    if(probability > 0.72){
        return 1;
    }
    else return 0;
}

void split_4(full_data *l, double *x, double *y, double *z,  int N, double startx, double endx, double starty, double endy){
    int n = 2;
    int splitX;
    int splitY;
    double Mx[2];
    double My[2];
    double Mz[2];
    data_MinMax_nonstat(x,Mx,N);
    data_MinMax_nonstat(y,My,N);
    startx = Mx[0];
    endx = Mx[1];
    starty = My[0];
    endy = My[1];
    double rangex = Mx[1] - Mx[0];
    double rangey = My[1] - My[0];
    if(rangex > rangey){
        splitX = 2;
        splitY = 1;
    }else{
        splitX = 1;
        splitY = 2;
    }
    int i, j, k, a[splitX][splitY];
    printf("counting.....\n");
    int count = 0;
    for (j = 0; j < splitX; j++){
        for (k = 0;k < splitY; k++){
            a[j][k] = count;
            count = count +1;  
        }
    }
    
    for (i = 0; i < n; i++){
        l[i].x = (double *) malloc(N * sizeof(double));
        l[i].y = (double *) malloc(N * sizeof(double));
        l[i].z = (double *) malloc(N * sizeof(double));
        l[i].index = 0;
    }
    
    
    
    
    // printf("%f\t%f\n",rangex,rangey);

    for (j = 0; j < splitX; j++){
        for (k = 0; k < splitY; k++){
            for (i = 0; i < N; i++){
                // printf("%d\t%d\t%d\tj_values: %f\t%f\tk_values: %f\t%f\n",j,k,i,Mx[0]+(j*rangex/(double)nsplit),Mx[0]+((j+1)*rangex/(double)nsplit),My[0]+(k*rangey/(double)nsplit),My[0]+((k+1)*rangey/(double)nsplit));
                // printf("%f\t%f\n",x[i],y[i]);
                if ((x[i] >= Mx[0]+(j*rangex/(double)splitX)) && (x[i] <= Mx[0]+((j+1)*rangex/(double)splitX))){
                    
                    if ((y[i] >= My[0]+(k*rangey/(double)splitY)) && (y[i] <= My[0]+((k+1)*rangey/(double)splitY))){
                        // printf("done till now........ \n");
                        l[a[j][k]].x[l[a[j][k]].index] = x[i];
                        l[a[j][k]].y[l[a[j][k]].index] = y[i];
                        l[a[j][k]].z[l[a[j][k]].index] = z[i];
                        l[a[j][k]].index = l[a[j][k]].index + 1;
                        // printf("%f\t%f\n",l[a[j][k]].x[l[a[j][k]].index],l[a[j][k]].y[l[a[j][k]].index]);

                    }
                    
                }
            }
        }
    }
printf("splitting done!!!!\n");
}

void data_split_nonstat_dynamic(double *x, double *y, double *z, MLE_data *data, int N){
    // full_data* l[n];
    int max_split = 20;
    double Mx[2];
    double My[2];
    double Mz[2];
    int n = 2;
    int i,j;
    printf("inside dynamic splitting function : \n");
    data_MinMax_nonstat(x,Mx,N);
    data_MinMax_nonstat(y,My,N);
    data_MinMax_nonstat(z,Mz,N);
    double startx = Mx[0];
    double endx = Mx[1];
    double starty = My[0];
    double endy = My[1];
    double startz = Mz[0];
    double endz = Mz[1];
    bool z_user_range = false;
    int val = check_stationarity(x,y,z,startx,endx,starty,endy,startz, endz, z_user_range,N);
    printf("\n \n data points in first split : %d\n", N);
    if(val == 1){
        full_data* l = malloc(max_split*sizeof *l);
        split_4(l,x,y,z,N,startx,endx,starty,endy);
        z_user_range = false;
        for(i = 0;i<n;i++){
            // int layer = 1;
            // printf("%f\t%f\t%f\t%f\n",startx,endx,starty,endy);
            // exit(0);
            // for(int t=0;t<15;t++){
            //     printf("%f___%f\n",l[i].x[t],l[i].y[t])
            // }
            data_MinMax_nonstat(l[i].x,Mx,l[i].index);
            data_MinMax_nonstat(l[i].y,My,l[i].index);
            startx = Mx[0];
            endx = Mx[1];
            starty = My[0];
            endy = My[1];
            printf("%f,%f,%f,%f,%d\n",startx,endx,starty,endy,l[i].index);
            val = check_stationarity(l[i].x,l[i].y,l[i].z,startx,endx,starty,endy,startz, endz, z_user_range,l[i].index);
            printf("\n \n data points in %d split : %d___ %d\n", i, l[i].index, val);
            // if(i == 0 || i==1)
            //     val = 0;
            if((val == 1) && (n+1<max_split) && (l[i].index > 10000)){

                for(j=n-1; j>i ;j--)
                    l[j+1] = l[j];

                full_data* l1 = malloc(2*sizeof *l1);
                split_4(l1,l[i].x,l[i].y,l[i].z,l[i].index,startx,endx,starty,endy);
                int x = 0;
                for(j=i;j<i+2;j++)
                   { l[j]=l1[x];
                    x = x+1;}

                n =n+1;
                i = i-1 ;
            }
            else{}
        }
        for(i = 0;i<n;i++){
            
            data_MinMax_nonstat(l[i].x,Mx,l[i].index);
            data_MinMax_nonstat(l[i].y,My,l[i].index);
            startx = Mx[0];
            endx = Mx[1];
            starty = My[0];
            endy = My[1];
            val = check_stationarity(l[i].x,l[i].y,l[i].z,startx,endx,starty,endy,startz, endz, z_user_range,l[i].index);
            printf("val of %dth splitted data_______ %d\n",i,val);
        }
        printf("%d\n",n);
        data->Nsplit = n;
        data->lm.x = (double *) malloc(n * sizeof(double));
        data->lm.y = (double *) malloc(n * sizeof(double));
        for(i=0; i<n; i++){
            data_MinMax_nonstat(l[i].x,Mx,l[i].index);
            data_MinMax_nonstat(l[i].y,My,l[i].index);
            double med_x = Mx[0] + (Mx[1] - Mx[0])/2;
            double med_y = My[0] + (My[1] - My[0])/2;
            data->lm.x[i] = med_x;
            data->lm.y[i] = med_y;
            printf("%f,%f\n",data->lm.x[i],data->lm.y[i]);
        }
        

        free(l);

    }
    else{
        printf("stationary data. Use the stationary kernel\n");
        exit(0);
    }

    
}



 void mergearray(double a[], double b[], int arr1size, int arr2size){
  
    // resultant Array Size Declaration
    int arr_resultsize = arr1size + arr2size;
    int i, j;
  
    // copying array 2 elements in to c array
    for (i = 0, j = arr1size;
         j < arr_resultsize && i < arr2size; i++, j++) {
        a[j] = b[i];
    }
}
double array_average(double a[], int n){
 int i;
 double sum = 0;
 double average = 0.0;
 /* calculate the sum of grades using for loop*/
    for(i = 0; i < n; i++){
       sum = sum + a[i];
    }
 average = sum/(double)n;

 return average;
}

void copy_array(double a[], double b[], int arr2size){
    for(int i=0;i<arr2size;i++){
        a[i] = b[i];
    }
}

double measure_distance(double x1, double x2, double y1, double y2){
    double distance = pow((x1-x2),2) + pow((y1-y2),2);
    return distance;
}

int data_Min_index(double *array, int N){
    int i;
    int a=0;
    double b = array[0];
    // size_t n = sizeof(array) / sizeof(array[0]);
    // printf("%d \t %d\n",sizeof(array),sizeof(array[0]));
    // printf("x vector size is %f", array[3120]);
    for(i=1; i<N; i++)
    {
        if(b > array[i])
		  a=i;    
    }
    // printf("minimum is %f \n",a[0]);
    // printf("maximum is %f \n",a[1]);
    return(a);
}

int data_Max_index(double *array, int N){
    int i;
    int a=0;
    double b = array[0];
    // size_t n = sizeof(array) / sizeof(array[0]);
    // printf("%d \t %d\n",sizeof(array),sizeof(array[0]));
    // printf("x vector size is %f", array[3120]);
    for(i=1; i<N; i++)
    {
        if(b>=array[i])
		  a=i;    
    }
    // printf("minimum is %f \n",a[0]);
    // printf("maximum is %f \n",a[1]);
    return(a);
}


void data_merge_nonstat_clustering(double *x, double *y, double *z, MLE_data *data, int N, int splitN){
    clock_t start, end;
    double cpu_time_used;
    struct stat st = {0};
    start = clock();
    // srand(0);
    int split = splitN;
    int n = 100;
    int i,j,k,a,random;
    double Mx[2];
    double My[2];
    double proba[n];
    double distance[split];

    double startx ;
    double endx ;
    double starty ;
    double endy ;
    bool z_user_range = false;
    data->lm.x = (double *) malloc(split * sizeof(double));
    data->lm.y = (double *) malloc(split * sizeof(double));
    full_data *l_nodes = malloc(n*sizeof *l_nodes);
    for (i = 0; i < n; i++){
        l_nodes[i].x = (double *) malloc(split * sizeof(double));
        l_nodes[i].y = (double *) malloc(split * sizeof(double));
        l_nodes[i].z = (double *) malloc(split * sizeof(double));
        l_nodes[i].index = 0;
    }
    full_data *l = malloc(split*sizeof *l);
    for (i = 0; i < n; i++){
        
        for(j=0;j<split;j++){
            if(j == 0){
                random = rand() % N;
                l_nodes[i].x[j] = x[random];
                l_nodes[i].y[j] = y[random];
                l_nodes[i].z[j] = z[random];
            }
            else{
                double distance1[j];
                a = 0;
                distance1[0] = 0;
                while(sqrt(distance1[a])< 0.5){
                    random = rand() % N;
                    for(k=j-1;k>=0;k--){
                        distance1[k] = measure_distance(l_nodes[i].x[k],x[random],l_nodes[i].y[k],y[random]);
                        // printf("%f\n",distance1[k]);
                    }
                    // printf("###############\n");
                    a = data_Min_index(distance1,j); 
                     
                }
                // printf("____%f____\n",distance1[a]);
                l_nodes[i].x[j] = x[random];
                l_nodes[i].y[j] = y[random];
                l_nodes[i].z[j] = z[random];
                
            }
        }
        // printf("%f__%f__%f__%f__%f__%f__%f__%f\n",l_nodes[i].x[0],l_nodes[i].y[0],l_nodes[i].x[1],l_nodes[i].y[1],
                                                // l_nodes[i].x[2],l_nodes[i].y[2],l_nodes[i].x[3],l_nodes[i].y[3]);
    }

    // for(i=0;i<n;i++){
    //     l_nodes[i].x[j] = x[random];
    //     l_nodes[i].y[j] = y[random];
    //     l_nodes[i].z[j] = z[random];
    // }
    l_nodes[0].x[0] = 0.25;
    l_nodes[0].y[0] = 0.25;
    l_nodes[0].z[0] = 0.25;

    l_nodes[0].x[1] = 0.75;
    l_nodes[0].y[1] = 0.75;
    l_nodes[0].z[1] = 0.25;

    // l_nodes[0].x[2] = 0.75;
    // l_nodes[0].y[2] = 0.25;
    // l_nodes[0].z[2] = 0.25;

    // l_nodes[0].x[3] = 0.75;
    // l_nodes[0].y[3] = 0.75;
    // l_nodes[0].z[3] = 0.25;

    // l_nodes[1].x[0] = 0.25;
    // l_nodes[1].y[0] = 0.30;
    // l_nodes[1].z[0] = 0.25;

    // l_nodes[1].x[1] = 0.65;
    // l_nodes[1].y[1] = 0.85;
    // l_nodes[1].z[1] = 0.25;

    // l_nodes[1].x[2] = 0.55;
    // l_nodes[1].y[2] = 0.50;
    // l_nodes[1].z[2] = 0.25;

    // l_nodes[1].x[3] = 0.35;
    // l_nodes[1].y[3] = 0.85;
    // l_nodes[1].z[3] = 0.25;

    for(i=0;i<n;i++){

        for (j = 0; j < split; j++){
            l[j].x = (double *) malloc(N * sizeof(double));
            l[j].y = (double *) malloc(N * sizeof(double));
            l[j].z = (double *) malloc(N * sizeof(double));
            l[j].index = 0;
        }
        
        for(k=0;k<N;k++){
            for(j=0;j<split;j++){
                distance[j] = measure_distance(l_nodes[i].x[j],x[k],l_nodes[i].y[j],y[k]);
            }
            a = data_Min_index(distance,split);
            l[a].x[l[a].index] = x[k];
            l[a].y[l[a].index] = y[k];
            l[a].z[l[a].index] = z[k];
            l[a].index = l[a].index +1;
        }
        
        proba[i] = 0;
        for(j=0;j<split;j++){

            data_MinMax_nonstat(x,Mx,N);
            printf("%f__%f\n",Mx[0],Mx[1]);
            Mx[0]=Mx[1]=My[0]=My[1] = 0;
            data_MinMax_nonstat(l[j].x,Mx,l[j].index);
            data_MinMax_nonstat(l[j].y,My,l[j].index);
            if(l[j].index < 1000){
                proba[i] = proba[i] + 100;
            }
            else{
                proba[i] = proba[i] + predict(l[j].x,l[j].y,l[j].z,Mx[0],Mx[1],My[0],My[1],Mx[0],Mx[1], z_user_range,l[j].index)/10000.0;
            }
            // printf("inside running loop 2 !!!!!!!!!!\n");

            // printf("%f__%f__%f__%f__%d\n",Mx[0],Mx[1],My[0],My[1],l[j].index);
            l_nodes[i].x[j] = array_average(l[j].x, l[j].index);
            l_nodes[i].y[j] = array_average(l[j].y, l[j].index);



            printf("%f\n",proba[i]);
            // l[j].x = NULL;
            // l[j].y = NULL;
            // l[j].z = NULL;
            // l[j].index = 0;
            // Mx[0]=Mx[1]=My[0]=My[1] = 0;

        }
    }
    a = data_Min_index(proba,n);
    printf("probability : %f\n",proba[a]);
    // char * nFileLoc  = (char *) malloc(100 * sizeof(char));
    // snprintf(nFileLoc, 100, "%s", "./probabilities_2sub.txt");
    // FILE *fptr1;
    
    // fptr1 = fopen(nFileLoc,"w");
    
    // fprintf(fptr1,"%f\n",proba[a]);
    // fclose(fptr1);

    for(i=0;i<split;i++){
        data->lm.x[i] = l_nodes[a].x[i];
        data->lm.y[i] = l_nodes[a].y[i];
        printf("%f_____%f\n",data->lm.x[i],data->lm.y[i]);
    }
    if (stat("./subregion_stationary", &st) == -1)
        mkdir("./subregion_stationary", 0700);
    for(j=0;j<split;j++){
                
        char * nFileLoc  = (char *) malloc(100 * sizeof(char));
        char * nFileZ  = (char *) malloc(100 * sizeof(char));
        // Put "file" then k then ".txt" in to filename.
        
        // snprintf(nFileZ, sizeof(char) * 32, "./subregion_stationary/",zi,j);
        snprintf(nFileLoc, 100, "%s%d", "./subregion_stationary/LOC_",j);
        snprintf(nFileZ, 100, "%s%d", "./subregion_stationary/Z_",j);
        printf("%s\n",nFileLoc);
        printf("%s\n",nFileZ);
        FILE *fptr1;
        FILE *fptr2;
        fptr1 = fopen(nFileLoc,"w");
        fptr2 = fopen(nFileZ,"w");
        for(i=0; i<l[j].index; i++){
            fprintf(fptr1,"%f,%f\n",l[j].x[i],l[j].y[i]);
            fprintf(fptr2,"%f\n",l[j].z[i]);
        }
        fclose(fptr1);
        fclose(fptr2);
        printf("\n######################################\n");
    }
end = clock();
cpu_time_used = ((double) (end - start))/CLOCKS_PER_SEC;
printf("time taken : %f\n", cpu_time_used);
}


// void data_merge_nonstat_dynamic(double *x, double *y, double *z, MLE_data *data, int N){
//     // full_data* l[n];
//     int split = 3;

//     int n = (int) split*split;
//     int merge_to = 3;
//     double Mx[2];
//     double My[2];
//     double Mz[2];
//     data->lm.x = (double *) malloc(merge_to * sizeof(double));
//     data->lm.y = (double *) malloc(merge_to * sizeof(double));
//     full_data_merge* l = malloc(n*sizeof *l);
//     full_data* l_merge = malloc(merge_to*sizeof *l_merge);
//     int i, j, k, i0,j0, a[split][split], b[split][split];
//     printf("counting.....\n");
//     int count = 0;
//     for (j = 0; j < split; j++){
//         for (k = 0;k < split; k++){
//             a[j][k] = count;
//             count = count +1;  
//         }
//     }
//     for (j = 0; j < split; j++){
//         for (k = 0;k < split; k++){
//             b[j][k] = a[j][split-k-1];
//         }
//     }

// //     for (int r=0; r<split; r++)
// // {
// //     for(int c=0; c<split; c++)
// //     {
// //          printf("%d     ", a[r][c]);
// //     }
// //     printf("\n");
// //     }


//     for (i = 0; i < n; i++){
//         l[i].x = (double *) malloc((int)(N/5) * sizeof(double));
//         l[i].y = (double *) malloc((int)(N/5) * sizeof(double));
//         l[i].z = (double *) malloc((int)(N/5) * sizeof(double));
//         l[i].index = 0;
//         l[i].subregion = (int)0;
//     }
    
//     data_MinMax_nonstat(x,Mx,N);
//     data_MinMax_nonstat(y,My,N);
//     data_MinMax_nonstat(z,Mz,N);
//     double rangex = Mx[1] - Mx[0];
//     double rangey = My[1] - My[0];

//     for (j = 0; j < split; j++){
//         for (k = 0; k < split; k++){
//             for (i = 0; i < N; i++){
//                 // printf("%d\t%d\t%d\tj_values: %f\t%f\tk_values: %f\t%f\n",j,k,i,Mx[0]+(j*rangex/(double)nsplit),Mx[0]+((j+1)*rangex/(double)nsplit),My[0]+(k*rangey/(double)nsplit),My[0]+((k+1)*rangey/(double)nsplit));
//                 // printf("%f\t%f\n",x[i],y[i]);
//                 if ((x[i] >= Mx[0]+(j*rangex/(double)split)) && (x[i] <= Mx[0]+((j+1)*rangex/(double)split))){
                    
//                     if ((y[i] >= My[0]+(k*rangey/(double)split)) && (y[i] <= My[0]+((k+1)*rangey/(double)split))){
//                         // printf("done till now........ \n");
//                         l[a[j][k]].x[l[a[j][k]].index] = x[i];
//                         l[a[j][k]].y[l[a[j][k]].index] = y[i];
//                         l[a[j][k]].z[l[a[j][k]].index] = z[i];
//                         l[a[j][k]].index = l[a[j][k]].index + 1;
//                         // printf("done!!!!!!!!!!!!!!!\n");

//                     }
                    
//                 }
//             }
//         }
//     }
//     for (i = 0; i < merge_to; i++){
//         l_merge[i].x = (double *) malloc(N * sizeof(double));
//         l_merge[i].y = (double *) malloc(N * sizeof(double));
//         l_merge[i].z = (double *) malloc(N * sizeof(double));
//         l_merge[i].index = 0;
//     }

//     int K = 1;
    
    
//     while(K < merge_to){
//         int val = 0;
//         i = 0;
//         j = 0;
        
//         while(l[a[i][j]].subregion != 0){
            
//             if(j == split-1){
//                 i = i + 1;
//                 i0 = i; 
//                 j = 0;
//             }else{
//                 j = j+1;
//                 j0 = j;
//             }
//         }
//         i0 = i;
//         j0 = j;
//         copy_array(l_merge[K-1].x,l[a[i][j]].x,l[a[i][j]].index);
//         copy_array(l_merge[K-1].y,l[a[i][j]].y,l[a[i][j]].index);
//         copy_array(l_merge[K-1].z,l[a[i][j]].z,l[a[i][j]].index);
//         l_merge[K-1].index = l[a[i][j]].index;
//         l[a[i][j]].subregion = K;

//         // printf("%f___%f\n",l_merge[K-1].x[0],l[a[i][j]].x[0]);
//         printf("___%d___%d\n",i,j);
//         // exit(0);
//         while(val == 0){
//             if(j == split -1){
//                 if(j0 == 0 || i0 == split-1 || l[a[i0+1][j0-1]].subregion != 0){
//                     // printf("%d___%d___%d___%d\n",i,j,l[a[i][j]].subregion,l[a[i][j]].index);
//                     full_data* l_test = malloc(1*sizeof *l_test);
//                     l_test[0].x = (double *) malloc((int)(N) * sizeof(double));
//                     l_test[0].y = (double *) malloc((int)(N) * sizeof(double));
//                     l_test[0].z = (double *) malloc((int)(N) * sizeof(double));
//                     l_test[0].index = 0;

//                     copy_array(l_test[0].x,l_merge[K-1].x,l_merge[K-1].index);
//                     copy_array(l_test[0].y,l_merge[K-1].y,l_merge[K-1].index);
//                     copy_array(l_test[0].z,l_merge[K-1].z,l_merge[K-1].index);
                    
//                     l_test[0].index =  l_test[0].index + l_merge[K-1].index;
//                     mergearray(l_test[0].x,l[a[i][j]].x,l_test[0].index, l[a[i][j]].index); 
//                     mergearray(l_test[0].y,l[a[i][j]].y,l_test[0].index, l[a[i][j]].index); 
//                     mergearray(l_test[0].z,l[a[i][j]].z,l_test[0].index, l[a[i][j]].index);  
//                     l_test[0].index =  l_test[0].index + l[a[i][j]].index;

//                     // printf("%f___%f___%d___%d\n",l_test[0].y[0],l[a[i][j-1]].y[0],i,j-1);
//                     // printf("%f___%f___%d___%d\n",l_test[0].y[l_test[0].index-1],l[a[i][j]].y[l[a[i][j]].index-1],i,j);
//                     // exit(0);

//                     data_MinMax_nonstat(l_test[0].x,Mx,l_test[0].index);
//                     data_MinMax_nonstat(l_test[0].y,My,l_test[0].index);
//                     data_MinMax_nonstat(l_test[0].z,Mz,l_test[0].index);
//                     double startx = Mx[0];
//                     double endx = Mx[1];
//                     double starty = My[0];
//                     double endy = My[1];
//                     double startz = Mz[0];
//                     double endz = Mz[1];
//                     bool z_user_range = false;
//                     // printf("%f,%f,%f,%f,%d\n",startx,endx,starty,endy,l[i].index);
//                     if(l_test[0].index < 10000){
//                         printf("___%d___%d\n",i,j);
//                         val = 0;
//                     }else{
//                         printf("___%d___%d\n",i,j);
//                         val = check_stationarity(l_test[0].x,l_test[0].y,l_test[0].z,
//                                     startx,endx,starty,endy,startz, endz, z_user_range,l_test[0].index);
//                     }

//                     if(val == 0){
                        
//                         copy_array(l_merge[K-1].x,l_test[0].x,l_test[0].index);
//                         copy_array(l_merge[K-1].y,l_test[0].y,l_test[0].index);
//                         copy_array(l_merge[K-1].z,l_test[0].z,l_test[0].index);
//                         l_merge[K-1].index = l_test[0].index;
//                         l[a[i][j]].subregion = K;
//                         // printf("___%d___%d\n",i,j);
//                         // printf("___%d\n",l_merge[K-1].index);
//                         j = 0 ;
//                         i = i+1;
//                         i0=i;
//                         j0=j;
//                         printf("I am here 1 loop !!!!\n");
//                     }  
//                     else if(val == 1 && j >= 1 && (i < split -1)){
//                         i = i+1;
//                         j = 0;
//                         val = 0;
//                         i0=i;
//                         j0=j;
//                         printf("I am here 3 loop !!!!!\n");
//                     }
//                     else if(i == split -1){
//                         printf("inside break level....\n");
//                         break;
//                     } 
//                     else{
//                         K = K+1;
//                     }
//                     free(l_test);  
//                 }
//                 else if(l[a[i0+1][j0-1]].subregion == 0){
//                     k = 0;
//                     // printf("%d___%d___%d___%d\n",i,j,l[a[i][j]].subregion,l[a[i][j]].index);
//                     full_data* l_test = malloc(1*sizeof *l_test);
//                     l_test[0].x = (double *) malloc((int)(N) * sizeof(double));
//                     l_test[0].y = (double *) malloc((int)(N) * sizeof(double));
//                     l_test[0].z = (double *) malloc((int)(N) * sizeof(double));
//                     l_test[0].index = 0;

//                     copy_array(l_test[0].x,l_merge[K-1].x,l_merge[K-1].index);
//                     copy_array(l_test[0].y,l_merge[K-1].y,l_merge[K-1].index);
//                     copy_array(l_test[0].z,l_merge[K-1].z,l_merge[K-1].index);
                    
//                     l_test[0].index =  l_test[0].index + l_merge[K-1].index;
//                     mergearray(l_test[0].x,l[a[i][j]].x,l_test[0].index, l[a[i][j]].index); 
//                     mergearray(l_test[0].y,l[a[i][j]].y,l_test[0].index, l[a[i][j]].index); 
//                     mergearray(l_test[0].z,l[a[i][j]].z,l_test[0].index, l[a[i][j]].index);  
//                     l_test[0].index =  l_test[0].index + l[a[i][j]].index;

//                     // printf("%f___%f___%d___%d\n",l_test[0].y[0],l[a[i][j-1]].y[0],i,j-1);
//                     // printf("%f___%f___%d___%d\n",l_test[0].y[l_test[0].index-1],l[a[i][j]].y[l[a[i][j]].index-1],i,j);
//                     // exit(0);

//                     data_MinMax_nonstat(l_test[0].x,Mx,l_test[0].index);
//                     data_MinMax_nonstat(l_test[0].y,My,l_test[0].index);
//                     data_MinMax_nonstat(l_test[0].z,Mz,l_test[0].index);
//                     double startx = Mx[0];
//                     double endx = Mx[1];
//                     double starty = My[0];
//                     double endy = My[1];
//                     double startz = Mz[0];
//                     double endz = Mz[1];
//                     bool z_user_range = false;
//                     // printf("%f,%f,%f,%f,%d\n",startx,endx,starty,endy,l[i].index);
//                     if(l_test[0].index < 10000){
//                         printf("___%d___%d\n",i,j);
//                         val = 0;
//                     }else{
//                         printf("___%d___%d\n",i,j);
//                         val = check_stationarity(l_test[0].x,l_test[0].y,l_test[0].z,
//                                     startx,endx,starty,endy,startz, endz, z_user_range,l_test[0].index);
//                     }

//                     if(val == 0){
                        
//                         copy_array(l_merge[K-1].x,l_test[0].x,l_test[0].index);
//                         copy_array(l_merge[K-1].y,l_test[0].y,l_test[0].index);
//                         copy_array(l_merge[K-1].z,l_test[0].z,l_test[0].index);
//                         l_merge[K-1].index = l_test[0].index;
//                         l[a[i][j]].subregion = K;
//                         // printf("___%d___%d\n",i,j);
//                         // printf("___%d\n",l_merge[K-1].index);
//                         j = 0 ;
//                         i = i+1;
//                         i0=i;
//                         j0=j;
//                         printf("I am here 1 loop !!!!\n");
//                     }  
//                     else if(val == 1 && j >= 1 && (i < split -1)){
//                         i = i+1;
//                         j = 0;
//                         val = 0;
//                         i0=i;
//                         j0=j;
//                         printf("I am here 3 loop !!!!!\n");
//                     }
//                     else if(i == split -1){
//                         printf("inside break level....\n");
//                         break;
//                     } 
//                     else{
//                         K = K+1;
//                     }
//                     free(l_test);  

//                 }

//             }
//             while(l[a[i][j]].subregion != 0){
//                 if(j == split-1){
//                     i = i + 1; 
//                     j = 0;
//                 }else{
//                     j = j+1;
//                 }
//             }
//             // printf("%d___%d___%d___%d\n",i,j,l[a[i][j]].subregion,l[a[i][j]].index);
//             full_data* l_test = malloc(1*sizeof *l_test);
//             l_test[0].x = (double *) malloc((int)(N) * sizeof(double));
//             l_test[0].y = (double *) malloc((int)(N) * sizeof(double));
//             l_test[0].z = (double *) malloc((int)(N) * sizeof(double));
//             l_test[0].index = 0;

//             copy_array(l_test[0].x,l_merge[K-1].x,l_merge[K-1].index);
//             copy_array(l_test[0].y,l_merge[K-1].y,l_merge[K-1].index);
//             copy_array(l_test[0].z,l_merge[K-1].z,l_merge[K-1].index);
            
//             l_test[0].index =  l_test[0].index + l_merge[K-1].index;
//             mergearray(l_test[0].x,l[a[i][j]].x,l_test[0].index, l[a[i][j]].index); 
//             mergearray(l_test[0].y,l[a[i][j]].y,l_test[0].index, l[a[i][j]].index); 
//             mergearray(l_test[0].z,l[a[i][j]].z,l_test[0].index, l[a[i][j]].index);  
//             l_test[0].index =  l_test[0].index + l[a[i][j]].index;

//             // printf("%f___%f___%d___%d\n",l_test[0].y[0],l[a[i][j-1]].y[0],i,j-1);
//             // printf("%f___%f___%d___%d\n",l_test[0].y[l_test[0].index-1],l[a[i][j]].y[l[a[i][j]].index-1],i,j);
//             // exit(0);

//             data_MinMax_nonstat(l_test[0].x,Mx,l_test[0].index);
//             data_MinMax_nonstat(l_test[0].y,My,l_test[0].index);
//             data_MinMax_nonstat(l_test[0].z,Mz,l_test[0].index);
//             double startx = Mx[0];
//             double endx = Mx[1];
//             double starty = My[0];
//             double endy = My[1];
//             double startz = Mz[0];
//             double endz = Mz[1];
//             bool z_user_range = false;
//             // printf("%f,%f,%f,%f,%d\n",startx,endx,starty,endy,l[i].index);
//             if(l_test[0].index < 10000){
//                 printf("___%d___%d\n",i,j);
//                 val = 0;
//             }else{
//                 printf("___%d___%d\n",i,j);
//                 val = check_stationarity(l_test[0].x,l_test[0].y,l_test[0].z,
//                             startx,endx,starty,endy,startz, endz, z_user_range,l_test[0].index);
//             }

//             if(val == 0 && j <= split -2){
                
//                 copy_array(l_merge[K-1].x,l_test[0].x,l_test[0].index);
//                 copy_array(l_merge[K-1].y,l_test[0].y,l_test[0].index);
//                 copy_array(l_merge[K-1].z,l_test[0].z,l_test[0].index);
//                 l_merge[K-1].index = l_test[0].index;
//                 l[a[i][j]].subregion = K;
//                 // printf("___%d___%d\n",i,j);
//                 // printf("___%d\n",l_merge[K-1].index);
//                 j = j+1 ;
//                 printf("I am here 1 loop !!!!\n");
//             }  
//             else if(val == 0 && (j) == split - 1 && (i < split -1)){
//                 copy_array(l_merge[K-1].x,l_test[0].x,l_test[0].index);
//                 copy_array(l_merge[K-1].y,l_test[0].y,l_test[0].index);
//                 copy_array(l_merge[K-1].z,l_test[0].z,l_test[0].index);
//                 l_merge[K-1].index = l_test[0].index;
//                 l[a[i][j]].subregion = K;
//                 i = i +1;
//                 j = 0;
                
//                 printf("I am here 2 loop !!!!\n");
//             } 
//             else if(val == 1 && j >= 1 && (i < split -1)){
//                 i = i+1;
//                 j = 0;
//                 val = 0;
//                 printf("I am here 3 loop !!!!!\n");
//             }
//             else if(i == split -1 && j == split -1){
//                 printf("inside break level....\n");
//                 break;
//             } 
//             else{
//                 K = K+1;
//             }
//             free(l_test);  

//         }
//     }

//     int i0 = 0;
//     int j0 = 0;
//     while(l[a[i0][j0]].subregion != 0){
//         if(j0 == split-1){
//             i0 = i0 + 1; 
//             j0 = 0;
//         }else{
//             j0 = j0+1;
//         }
        
        
//     }

//     copy_array(l_merge[K-1].x,l[a[i0][j0]].x,l[a[i0][j0]].index);
//     copy_array(l_merge[K-1].y,l[a[i0][j0]].y,l[a[i0][j0]].index);
//     copy_array(l_merge[K-1].z,l[a[i0][j0]].z,l[a[i0][j0]].index);
//     l_merge[K-1].index = l[a[i0][j0]].index;
//     if(j0 == split - 1){
//         i0 = i0 + 1;
//         j0 = 0;
//     }else{
//         j0 = j0 + 1;
//     }
//     for (i = i0; i < split; i++){
//         for (j = j0; j < split; j++){
//             mergearray(l_merge[K-1].x,l[a[i][j]].x,l_merge[K-1].index, l[a[i][j]].index); 
//             mergearray(l_merge[K-1].y,l[a[i][j]].y,l_merge[K-1].index, l[a[i][j]].index); 
//             mergearray(l_merge[K-1].z,l[a[i][j]].z,l_merge[K-1].index, l[a[i][j]].index);  
//             l_merge[K-1].index =  l_merge[K-1].index + l[a[i][j]].index;
//         }
//     }
//     printf("%f___%f\n",l_merge[0].x[0],l_merge[0].x[l_merge[0].index - 1]);
//     double u = array_average(l_merge[i].x,l_merge[i].index);
//     printf("%f\n",u);
//     // printf("%f___%f\n",l_merge[0].x[0],l_merge[0].x[l_merge[0].index - 1]);
//     printf("%f___%f\n",l_merge[1].x[0],l_merge[1].x[l_merge[1].index - 1]);
//     for(i=0; i<merge_to; i++){
        
//         data->lm.x[i] = array_average(l_merge[i].x,l_merge[i].index);
//         data->lm.y[i] = array_average(l_merge[i].y,l_merge[i].index);
//         printf("%f,%f\n",data->lm.x[i],data->lm.y[i]);
//     }

//     free(l);
//     free(l_merge);

// }



void data_split_nonstat_userDefined(double *x, double *y, double *z, MLE_data *data, int N, int splitX, int splitY, int zi){
    // full_data* l[n];
    int n = (int) splitX*splitY;
    double Mx[2];
    double My[2];
    double Mz[2];
    data->lm.x = (double *) malloc(n * sizeof(double));
    data->lm.y = (double *) malloc(n * sizeof(double));
    full_data* l = malloc(n*sizeof *l);
    int i, j, k, a[splitX][splitY];
    printf("counting.....\n");
    int count = 0;
    for (j = 0; j < splitX; j++){
        for (k = 0;k < splitY; k++){
            a[j][k] = count;
            count = count +1;  
        }
    }

    if(strcmp(data->kernel_fun, "univariate_matern_non_stat")   == 0){
        for (j = 0; j < splitX; j++){
            for (k = 0; k < splitY; k++){
                data->lm.x[a[j][k]] = ((j)/(double)splitX) + (((j+1)/(double)splitX) - ((j)/(double)splitX))/2.0;
                data->lm.y[a[j][k]] = ((k)/(double)splitY) + (((k+1)/(double)splitY) - ((k)/(double)splitY))/2.0;       
            }
        }
        printf("generation stage....\n");
        // data->lm.x = mx;
        // data->lm.y = my;
    }
    else if(strcmp(data->kernel_fun, "univariate_matern_non_stat_estimation")   == 0){
        for (i = 0; i < n; i++){
            l[i].x = (double *) malloc(N * sizeof(double));
            l[i].y = (double *) malloc(N * sizeof(double));
            l[i].z = (double *) malloc(N * sizeof(double));
            l[i].index = 0;
        }
        
        data_MinMax_nonstat(x,Mx,N);
        data_MinMax_nonstat(y,My,N);
        data_MinMax_nonstat(z,Mz,N);
        double rangex = Mx[1] - Mx[0];
        double rangey = My[1] - My[0];

        for (j = 0; j < splitX; j++){
            for (k = 0; k < splitY; k++){
                for (i = 0; i < N; i++){
                    // printf("%d\t%d\t%d\tj_values: %f\t%f\tk_values: %f\t%f\n",j,k,i,Mx[0]+(j*rangex/(double)nsplit),Mx[0]+((j+1)*rangex/(double)nsplit),My[0]+(k*rangey/(double)nsplit),My[0]+((k+1)*rangey/(double)nsplit));
                    // printf("%f\t%f\n",x[i],y[i]);
                    if ((x[i] >= Mx[0]+(j*rangex/(double)splitX)) && (x[i] <= Mx[0]+((j+1)*rangex/(double)splitX))){
                        
                        if ((y[i] >= My[0]+(k*rangey/(double)splitY)) && (y[i] <= My[0]+((k+1)*rangey/(double)splitY))){
                            // printf("done till now........ \n");
                            l[a[j][k]].x[l[a[j][k]].index] = x[i];
                            l[a[j][k]].y[l[a[j][k]].index] = y[i];
                            l[a[j][k]].z[l[a[j][k]].index] = z[i];
                            l[a[j][k]].index = l[a[j][k]].index + 1;
                            // printf("done!!!!!!!!!!!!!!!\n");

                        }
                        
                    }
                }
            }
        }
        if (zi != 1000000){
            for(j=0;j<n;j++){
                
                char * nFileLoc  = (char *) malloc(100 * sizeof(char));
                char * nFileZ  = (char *) malloc(100 * sizeof(char));
                // Put "file" then k then ".txt" in to filename.
                
                // snprintf(nFileZ, sizeof(char) * 32, "./subregion_stationary/",zi,j);
                snprintf(nFileLoc, 100, "%s%d%d", "./subregion_stationary/LOC_",zi,j);
                snprintf(nFileZ, 100, "%s%d%d", "./subregion_stationary/Z_",zi,j);
                printf("%s\n",nFileLoc);
                printf("%s\n",nFileZ);
                FILE *fptr1;
                FILE *fptr2;
                fptr1 = fopen(nFileLoc,"w");
                fptr2 = fopen(nFileZ,"w");
                for(i=0; i<l[j].index; i++){
                    fprintf(fptr1,"%f,%f\n",l[j].x[i],l[j].y[i]);
                    fprintf(fptr2,"%f\n",l[j].z[i]);
                }
                fclose(fptr1);
                fclose(fptr2);
                printf("\n######################################\n");
            }
        }
        

        


        for(i=0; i<n; i++)
        {
            data_MinMax_nonstat(l[i].x,Mx,l[i].index);
            data_MinMax_nonstat(l[i].y,My,l[i].index);
            double med_x = Mx[0] + (Mx[1] - Mx[0])/2;
            double med_y = My[0] + (My[1] - My[0])/2;
            data->lm.x[i] = med_x;
            data->lm.y[i] = med_y;
        }
        // data->lm.x = mx;
        // data->lm.y = my;
        for(i=0; i<n; i++){
            printf("%f,%f\n",data->lm.x[i],data->lm.y[i]);
        }

    }
    free(l);

}

static double gaussian_kernel(double x1, double y1, double x2, double y2, double h){
    double distance = pow((x2 - x1), 2) + pow((y2 - y1), 2);
    double k = pow(pow_e, (-distance)/(2*h));
    return k;
}

void parameter_array(MLE_data *data, int N1, int N2, location *l1, location *l2, location *lm, double *localtheta){
    double l1x, l1y, l2x, l2y, lmx, lmy;
    int i,j;
    // double expr = 0.0;
    // double con, sigma_square, beta, nu;
    double h = 0.09;

    data->ns_params.sigma_arr1 = NULL;
    data->ns_params.lamda_arr1 = NULL;
    data->ns_params.sigma_arr2 = NULL;
    data->ns_params.lamda_arr2 = NULL;
    data->ns_params.nu_arr1 = NULL;
    data->ns_params.nu_arr2 = NULL;

    data->ns_params.sigma_arr1 = (double *) malloc(N1 * sizeof(double));
    data->ns_params.lamda_arr1 = (double *) malloc(N1 * sizeof(double));
    data->ns_params.nu_arr1 = (double *) malloc(N1 * sizeof(double));
    data->ns_params.sigma_arr2 = (double *) malloc(N2 * sizeof(double));
    data->ns_params.lamda_arr2 = (double *) malloc(N2 * sizeof(double));
    data->ns_params.nu_arr2 = (double *) malloc(N2 * sizeof(double));

    // double sum1 =0;
    // double sum2 = 0;
    // for(i=0; i < 4 ; i ++){
    //     sum1+=lm->x[i];
    //     sum2+=lm->y[i];
    // }
    // printf("%f__%f\n",sum1,sum2);

    // if(N1 == N2){
    //     if(l1->x[0] == l2->x[0] && l1->x[1] == l2->x[1] && l1->x[2] == l2->x[2]){
    //         for(i=0; i<N1; i++){
    //             l1x = l1->x[i];
    //             l1y = l1->y[i];
    //             double weight = 0.0;
    //             double gaussian_kernel_sum = 0;
    //             double sigma1 = 0.0;
    //             double lambda1 = 0.0;
    //             double nu1 = 0.0;
    //             double *gaussian_kernel_list = calloc(data->Nsplit,sizeof(double));
    //             for(j=0; j<data->Nsplit; j++){
    //                 lmx = lm->x[j];
    //                 lmy = lm->y[j];
    //                 gaussian_kernel_list[j] = gaussian_kernel(l1x,l1y,lmx,lmy,h);
    //                 gaussian_kernel_sum = gaussian_kernel_sum + gaussian_kernel_list[j];
    //             }
    //             // printf("%f\n",gaussian_kernel_sum);
    //             for(j=0; j<data->Nsplit; j++){
    //                 weight = gaussian_kernel_list[j]/gaussian_kernel_sum;
    //                 // sigma1 = sigma1 + (weight*localtheta[j*2]);
    //                 // lambda1 = lambda1 + (weight*localtheta[j*2 + 1]);
    //                 lambda1 = lambda1 + (weight*localtheta[j*2]);
    //                 sigma1 = sigma1 + (weight*localtheta[j*2 + 1]);
    //                 nu1 = nu1 + (weight*localtheta[j*2 + 2]);
    //                 // printf("%f\n",weight);
    //                 // printf("%f,%f,%f\n",weight,localtheta[j*3 + 1],localtheta[j*3 + 2]);
    //                 // printf("%f________\n",nu1);
    //             }
    //             // printf("####3#####\n");
    //             data->ns_params.sigma_arr1[i] = sigma1;
    //             data->ns_params.lamda_arr1[i] = lambda1;
    //             data->ns_params.nu_arr1[i]    = nu1;
    //             data->ns_params.sigma_arr2[i] = sigma1;
    //             data->ns_params.lamda_arr2[i] = lambda1;
    //             data->ns_params.nu_arr2[i]    = nu1;
                
    //             sigma1=0.0;
    //             lambda1=0.0;
    //             nu1 = 0.0;
    //             // printf("%f final sigma\n",nu_arr_1[i]);
    //         }

    //     }
    // }
    // else{
        for(i=0; i<N1; i++){
                l1x = l1->x[i];
                l1y = l1->y[i];
                double weight = 0.0;
                double gaussian_kernel_sum = 0;
                double sigma1 = 0.0;
                double lambda1 = 0.0;
                double nu1 = 0.0;
                double *gaussian_kernel_list = calloc(data->Nsplit,sizeof(double));
                for(j=0; j<data->Nsplit; j++){
                    lmx = lm->x[j];
                    lmy = lm->y[j];
                    gaussian_kernel_list[j] = gaussian_kernel(l1x,l1y,lmx,lmy,h);
                    gaussian_kernel_sum = gaussian_kernel_sum + gaussian_kernel_list[j];
                }
                // printf("%f\n",gaussian_kernel_sum);
                for(j=0; j<data->Nsplit; j++){
                    weight = gaussian_kernel_list[j]/gaussian_kernel_sum;
                    sigma1 = sigma1 + (weight*localtheta[j*3 + 1]);
                    lambda1 = lambda1 + (weight*localtheta[j*3]);
                    nu1 = nu1 + (weight*localtheta[j*3 + 2]);
                    // printf("%f\n",weight);
                    // printf("%f,%f,%f\n",weight,localtheta[j*3 + 1],localtheta[j*3 + 2]);
                    // printf("%f________\n",nu1);
                }
                // printf("####3#####\n");
                data->ns_params.sigma_arr1[i] = sigma1;
                data->ns_params.lamda_arr1[i] = lambda1;
                data->ns_params.nu_arr1[i]    = nu1;
                sigma1=0.0;
                lambda1=0.0;
                nu1=0.0;
        
        }
        for(i=0; i<N2; i++){
                l2x = l2->x[i];
                l2y = l2->y[i];
                double weight = 0.0;
                double gaussian_kernel_sum = 0;
                double sigma1 = 0.0;
                double lambda1 = 0.0;
                double nu1 = 0.0;
                double *gaussian_kernel_list = calloc(data->Nsplit,sizeof(double));
                for(j=0; j<data->Nsplit; j++){
                    lmx = lm->x[j];
                    lmy = lm->y[j];
                    gaussian_kernel_list[j] = gaussian_kernel(l2x,l2y,lmx,lmy,h);
                    gaussian_kernel_sum = gaussian_kernel_sum + gaussian_kernel_list[j];
                }
                // printf("%f\n",gaussian_kernel_sum);
                for(j=0; j<data->Nsplit; j++){
                    weight = gaussian_kernel_list[j]/gaussian_kernel_sum;
                    sigma1 = sigma1 + (weight*localtheta[j*3 + 1]);
                    lambda1 = lambda1 + (weight*localtheta[j*3]);
                    nu1 = nu1 + (weight*localtheta[j*3 + 2]);
                    // printf("%f\n",weight);
                    // printf("%f,%f,%f\n",weight,localtheta[j*3 + 1],localtheta[j*3 + 2]);
                    // printf("%f________\n",nu1);
                }
                // printf("####3#####\n");
                data->ns_params.sigma_arr2[i] = sigma1;
                data->ns_params.lamda_arr2[i] = lambda1;
                data->ns_params.nu_arr2[i]    = nu1;
                sigma1=0.0;
                lambda1=0.0;
                nu1=0.0;
        }
// }
    
}
