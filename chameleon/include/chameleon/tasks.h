/**
 *
 * @file chameleon_tasks.h
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2021 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon elementary tasks main header
 *
 * @version 1.0.0
 * @author Mathieu Faverge
 * @author Cedric Augonnet
 * @author Florent Pruvost
 * @date 2020-03-03
 *
 */
#ifndef _chameleon_tasks_h_
#define _chameleon_tasks_h_

#include "chameleon/config.h"

BEGIN_C_DECLS

/**
 * @brief Kernel enum
 *
 * Those enums are Used to apply operations on specific kernels, and or for
 * tracing/profiling.
 */
typedef enum chameleon_tasktype_e {

  TASK_GEMM,
  TASK_HEMM,
  TASK_HER2K,
  TASK_HERK,
  TASK_SYTRF_NOPIV,
  TASK_SYMM,
  TASK_SYR2K,
  TASK_SYRK,
  TASK_TRSM,
  TASK_TRMM,

  TASK_GELQT,
  TASK_GEQRT,
  TASK_GESSM,
  TASK_GETRF,
  TASK_GETRF_INCPIV,
  TASK_GETRF_NOPIV,
  TASK_LAUUM,
  TASK_ORMLQ,
  TASK_ORMQR,
  TASK_POTRF,
  TASK_SSSSM,
  TASK_TPLQT,
  TASK_TPMLQT,
  TASK_TPMQRT,
  TASK_TPQRT,
  TASK_TRTRI,
  TASK_TSTRF,
  TASK_UNMLQ,
  TASK_UNMQR,

  TASK_GEADD,
  TASK_LASCAL,
  TASK_LACPY,
  TASK_LAG2C,
  TASK_LAG2Z,
  TASK_LANGE,
  TASK_LANHE,
  TASK_LANSY,
  TASK_LASET,
  TASK_LASET2,
  TASK_PEMV,
  TASK_PLGHE,
  TASK_PLGSY,
  TASK_PLRNT,
  TASK_TILE_ZERO,

  TASK_NBKERNELS
} cham_tasktype_t;

#define TASK_TSLQT TASK_TPLQT
#define TASK_TSMLQ TASK_TPMLQT
#define TASK_TSMQR TASK_TPMQRT
#define TASK_TSQRT TASK_TPQRT
#define TASK_TTLQT TASK_TPLQT
#define TASK_TTMLQ TASK_TPMLQT
#define TASK_TTMQR TASK_TPMQRT
#define TASK_TTQRT TASK_TPQRT

typedef int (*cham_unary_operator_t)( const CHAM_desc_t *desc,
                                      cham_uplo_t uplo, int m, int n,
                                      CHAM_tile_t *data, void *op_args );

void INSERT_TASK_map( const RUNTIME_option_t *options,
                      cham_uplo_t uplo, const CHAM_desc_t *A, int Am, int An,
                      cham_unary_operator_t op_fct, void *op_args );

#include "chameleon/tasks_z.h"
#include "chameleon/tasks_d.h"
#include "chameleon/tasks_c.h"
#include "chameleon/tasks_s.h"
#include "chameleon/tasks_zc.h"
#include "chameleon/tasks_ds.h"

END_C_DECLS

#endif /* _chameleon_tasks_h_ */
