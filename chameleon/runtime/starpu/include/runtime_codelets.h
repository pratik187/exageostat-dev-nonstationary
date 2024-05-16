 /**
 *
 * @file starpu/runtime_codelets.h
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2021 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon StarPU codelets main header
 *
 * @version 1.1.0
 * @author Cedric Augonnet
 * @author Mathieu Faverge
 * @author Cedric Castagnede
 * @author Florent Pruvost
 * @date 2020-06-29
 *
 */
#ifndef _runtime_codelets_h_
#define _runtime_codelets_h_

#include "chameleon/config.h"
#include "runtime_codelet_profile.h"

//#undef STARPU_CUDA_ASYNC
#ifdef STARPU_CUDA_ASYNC
#define CODELET_CUDA_FLAGS(flags) .cuda_flags = {(flags)},
#else
#define CODELET_CUDA_FLAGS(flags)
#endif

#define CODELETS_ALL(cl_name, cpu_func_name, cuda_func_name, _original_location_, cuda_flags) \
    struct starpu_perfmodel cl_##cl_name##_fake = {                     \
        .type   = STARPU_HISTORY_BASED,                                 \
        .symbol = "fake_"#cl_name                                       \
    };                                                                  \
                                                                        \
    struct starpu_perfmodel cl_##cl_name##_model = {                    \
        .type   = STARPU_HISTORY_BASED,                                 \
        .symbol = ""#cl_name                                            \
    };                                                                  \
                                                                        \
    struct starpu_codelet cl_##cl_name = {                              \
        .where     = (_original_location_),                             \
        .cpu_func  = ((cpu_func_name)),                                 \
        CODELET_CUDA_FLAGS(cuda_flags)                                  \
        .cuda_func = ((cuda_func_name)),                                \
        .nbuffers  = STARPU_VARIABLE_NBUFFERS,                          \
        .model     = &cl_##cl_name##_model,                             \
        .name      = #cl_name                                           \
    };                                                                  \
                                                                        \
    void cl_##cl_name##_restrict_where(uint32_t where)                  \
    {                                                                   \
        if ( cl_##cl_name.where & where )                               \
            cl_##cl_name.where = (cl_##cl_name.where & where);          \
    }                                                                   \
                                                                        \
    void cl_##cl_name##_restore_where(void)                             \
    {                                                                   \
        cl_##cl_name.where = (_original_location_);                     \
    }                                                                   \
                                                                        \
    void cl_##cl_name##_restore_model(void)                             \
    {                                                                   \
        cl_##cl_name.model = &cl_##cl_name##_model;                     \
    }

#if defined(CHAMELEON_SIMULATION)
#define CODELETS_CPU(name, cpu_func_name)                    \
    CODELETS_ALL( name, (starpu_cpu_func_t) 1, NULL, STARPU_CPU, 0 )
#else
#define CODELETS_CPU(name, cpu_func_name)                    \
    CODELETS_ALL( name, cpu_func_name, NULL, STARPU_CPU, 0 )
#endif

#define CODELETS_GPU(name, cpu_func_name, cuda_func_name, cuda_flags) \
    CODELETS_ALL( name, cpu_func_name, cuda_func_name, STARPU_CPU  | STARPU_CUDA, cuda_flags )

#define CODELETS_ALL_HEADER(name)                            \
     CHAMELEON_CL_CB_HEADER(name);                           \
     void cl_##name##_load_fake_model(void);                 \
     void cl_##name##_restore_model(void);                   \
     extern struct starpu_codelet cl_##name;                 \
     void cl_##name##_restrict_where(uint32_t where);        \
     void cl_##name##_restore_where(void)

#if defined(CHAMELEON_SIMULATION)
#if defined(CHAMELEON_USE_CUDA)
#define CODELETS(name, cpu_func_name, cuda_func_name, cuda_flags) \
    CODELETS_GPU(name, (starpu_cpu_func_t) 1, (starpu_cuda_func_t) 1, cuda_flags)

#define CODELETS_HEADER(name)  CODELETS_ALL_HEADER(name)
#else
#define CODELETS(name, cpu_func_name, cuda_func_name, cuda_flags) \
    CODELETS_CPU(name, (starpu_cpu_func_t) 1)

#define CODELETS_HEADER(name)  CODELETS_ALL_HEADER(name)
#endif
#elif defined(CHAMELEON_USE_CUDA)
#define CODELETS(name, cpu_func_name, cuda_func_name, cuda_flags) \
    CODELETS_GPU(name, cpu_func_name, cuda_func_name, cuda_flags)

#define CODELETS_HEADER(name)  CODELETS_ALL_HEADER(name)
#else
#define CODELETS(name, cpu_func_name, cuda_func_name, cuda_flags) \
    CODELETS_CPU(name, cpu_func_name)

#define CODELETS_HEADER(name)  CODELETS_ALL_HEADER(name)
#endif

CODELETS_HEADER(map);

#endif /* _runtime_codelets_h_ */
