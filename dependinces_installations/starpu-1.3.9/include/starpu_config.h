/* include/starpu_config.h.  Generated from starpu_config.h.in by configure.  */
/* StarPU --- Runtime system for heterogeneous multicore architectures.
 *
 * Copyright (C) 2009-2021  Université de Bordeaux, CNRS (LaBRI UMR 5800), Inria
 * Copyright (C) 2020       Federal University of Rio Grande do Sul (UFRGS)
 *
 * StarPU is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation; either version 2.1 of the License, or (at
 * your option) any later version.
 *
 * StarPU is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 *
 * See the GNU Lesser General Public License in COPYING.LGPL for more details.
 */

/*
 * This is the public config.h file, installed along libstarpu.
 *
 * It should only contain the build-time #defines which have an effect on the
 * API & ABI.
 */

#ifndef __STARPU_CONFIG_PUBLIC_H__
#define __STARPU_CONFIG_PUBLIC_H__

/**
   Define the major version of StarPU. This is the version used when
   compiling the application.
   @ingroup API_Versioning
*/
#define STARPU_MAJOR_VERSION 1

/**
   @ingroup API_Versioning
   Define the minor version of StarPU. This is the version used when
   compiling the application.
*/
#define STARPU_MINOR_VERSION 3

/**
   Define the release version of StarPU. This is the version used when
   compiling the application.
   @ingroup API_Versioning
*/
#define STARPU_RELEASE_VERSION 9

#define STARPU_USE_CPU 1

/**
   Defined when StarPU has been installed with
   CUDA support. It should be used in your code to detect the
   availability of CUDA.
   @ingroup API_CUDA_Extensions
*/
/* #undef STARPU_USE_CUDA */

/**
   Defined when StarPU has been installed with
   NVidia-ML support. It should be used in your code to detect the
   availability of NVML-related functions.
   @ingroup API_CUDA_Extensions
*/
/* #undef STARPU_HAVE_LIBNVIDIA_ML */

/**
   Defined when StarPU has been installed with OpenCL support. It
   should be used in your code to detect the availability of OpenCL as
   shown in \ref FullSourceCodeVectorScal.
   @ingroup API_OpenCL_Extensions
*/
/* #undef STARPU_USE_OPENCL */

/**
   Defined when StarPU has been installed with MIC support. It should
   be used in your code to detect the availability of MIC.
   @ingroup API_MIC_Extensions
*/
/* #undef STARPU_USE_MIC */

/**
   Defined when StarPU has been installed with MPI Master Slave
   support. It should be used in your code to detect the availability
   of MPI Master Slave.
   @ingroup API_MPI_Support
*/
/* #undef STARPU_USE_MPI_MASTER_SLAVE */

/**
   Defined when StarPU has been installed with OpenMP Runtime support.
   It should be used in your code to detect the availability of the
   runtime support for OpenMP.
   @ingroup API_OpenMP_Runtime_Support
*/
/* #undef STARPU_OPENMP */

/* #undef STARPU_CLUSTER */

/* #undef STARPU_SIMGRID */
/* #undef STARPU_SIMGRID_MC */
/* #undef STARPU_SIMGRID_HAVE_XBT_BARRIER_INIT */
/* #undef STARPU_HAVE_SIMGRID_MSG_H */
/* #undef STARPU_HAVE_MSG_MSG_H */
/* #undef STARPU_HAVE_SIMGRID_ACTOR_H */
/* #undef STARPU_HAVE_SIMGRID_SEMAPHORE_H */
/* #undef STARPU_HAVE_SIMGRID_MUTEX_H */
/* #undef STARPU_HAVE_SIMGRID_COND_H */
/* #undef STARPU_HAVE_SIMGRID_BARRIER_H */
/* #undef STARPU_HAVE_XBT_SYNCHRO_H */
/* #undef STARPU_HAVE_VALGRIND_H */
/* #undef STARPU_HAVE_MEMCHECK_H */
/* #undef STARPU_VALGRIND_FULL */
/* #undef STARPU_SANITIZE_LEAK */
#define STARPU_NON_BLOCKING_DRIVERS 1
/* workers must call callbacks on sleep/wake-up */
/* #undef STARPU_WORKER_CALLBACKS */

/* #undef STARPU_HAVE_ICC */

/**
   Defined when StarPU has been installed with MPI support. It should
   be used in your code to detect the availability of MPI.
   @ingroup API_MPI_Support
*/
/* #undef STARPU_USE_MPI */
/* #undef STARPU_USE_MPI_MPI */
/* #undef STARPU_USE_MPI_NMAD */

/* #undef STARPU_ATLAS */
/* #undef STARPU_GOTO */
/* #undef STARPU_OPENBLAS */
#define STARPU_MKL 1
/* #undef STARPU_ARMPL */
/* #undef STARPU_SYSTEM_BLAS */
#define STARPU_HAVE_BLAS 1

/**
   Define the directory in which the OpenCL codelets of the
   applications provided with StarPU have been installed.
   @ingroup API_OpenCL_Extensions
*/
/* #undef STARPU_OPENCL_DATADIR */
/* #undef STARPU_HAVE_MAGMA */

/* #undef STARPU_OPENGL_RENDER */
/* #undef STARPU_USE_GTK */
#define STARPU_HAVE_X11 1

#define STARPU_HAVE_POSIX_MEMALIGN 1

#define STARPU_HAVE_MEMALIGN 1

#define STARPU_HAVE_MALLOC_H 1

#define STARPU_HAVE_SYNC_BOOL_COMPARE_AND_SWAP 1
#define STARPU_HAVE_SYNC_VAL_COMPARE_AND_SWAP 1
#define STARPU_HAVE_SYNC_FETCH_AND_ADD 1
#define STARPU_HAVE_SYNC_FETCH_AND_OR 1
#define STARPU_HAVE_SYNC_LOCK_TEST_AND_SET 1
#define STARPU_HAVE_SYNC_SYNCHRONIZE 1
#define STARPU_HAVE_ATOMIC_EXCHANGE_N 1

/* #undef STARPU_DEVEL */
/* #undef STARPU_MODEL_DEBUG */
/* #undef STARPU_NO_ASSERT */
/* #undef STARPU_DEBUG */
/* #undef STARPU_VERBOSE */
#define STARPU_GDB_PATH "/usr/bin/gdb"

/* #undef STARPU_HAVE_FFTW */
/* #undef STARPU_HAVE_FFTWF */
/* #undef STARPU_HAVE_FFTWL */
/* #undef STARPU_HAVE_CUFFTDOUBLECOMPLEX */

/* #undef STARPU_HAVE_CURAND */

/**
   Define the maximum number of memory nodes managed by StarPU. The
   default value can be modified at configure by using the option \ref
   enable-maxnodes "--enable-maxnodes". Reducing it allows to
   considerably reduce memory used by StarPU data structures.
   @ingroup API_Workers_Properties
*/
#define STARPU_MAXNODES 4

/**
   Define the maximum number of buffers that tasks will be able to
   take as parameters. The default value is 8, it can be changed by
   using the configure option \ref enable-maxbuffers
   "--enable-maxbuffers".
   @ingroup API_Codelet_And_Tasks
*/
#define STARPU_NMAXBUFS 8

/**
   Define the maximum number of CPU workers managed by StarPU. The
   default value can be modified at configure by using the option \ref
   enable-maxcpus "--enable-maxcpus".
   @ingroup API_Workers_Properties
*/
#define STARPU_MAXCPUS 64

/**
   Define the maximum number of NUMA nodes managed by StarPU. The
   default value can be modified at configure by using the option \ref
   enable-maxnumanodes "--enable-maxnumanodes".
   @ingroup API_Workers_Properties
*/
#define STARPU_MAXNUMANODES 2

/**
   Define the maximum number of CUDA devices that are supported by StarPU.
   @ingroup API_CUDA_Extensions
*/
#define STARPU_MAXCUDADEVS 4

/**
   Define the maximum number of OpenCL devices that are supported by
   StarPU.
   @ingroup API_OpenCL_Extensions
*/
#define STARPU_MAXOPENCLDEVS 8

/**
   Define the maximum number of MIC devices that are supported by
   StarPU.
   @ingroup API_MIC_Extensions
*/
#define STARPU_MAXMICDEVS 0

/**
   Define the maximum number of workers managed by StarPU.
   @ingroup API_Workers_Properties
*/
#define STARPU_NMAXWORKERS 64

/**
   Define the maximum number of scheduling contexts managed by StarPU.
   The default value can be modified at configure by using the option
   \ref enable-max-sched-ctxs "--enable-max-sched-ctxs".
   @ingroup API_Scheduling_Policy
*/
#define STARPU_NMAX_SCHED_CTXS 10

/**
   Define the maximum number of implementations per architecture. The
   default value can be modified at configure by using the option \ref
   enable-maximplementations "--enable-maximplementations".
   @ingroup API_Scheduling_Policy
*/
#define STARPU_MAXIMPLEMENTATIONS 4

#define STARPU_MAXMPKERNELS 10
/* #undef STARPU_USE_SC_HYPERVISOR */
/* #undef STARPU_SC_HYPERVISOR_DEBUG */
/* #undef STARPU_HAVE_GLPK_H */

/* #undef STARPU_HAVE_CUDA_MEMCPY_PEER */
/* #undef STARPU_HAVE_LIBNUMA */

/* #undef STARPU_HAVE_WINDOWS */
#define STARPU_LINUX_SYS 1
#define STARPU_HAVE_SETENV 1
#define STARPU_HAVE_UNSETENV 1
#define STARPU_HAVE_UNISTD_H 1
/* #undef STARPU_HAVE_HDF5 */

/* #undef STARPU_USE_FXT */
/* #undef STARPU_FXT_LOCK_TRACES */

#ifdef _MSC_VER
typedef long starpu_ssize_t;
#define __starpu_func__ __FUNCTION__
#else
#  include <sys/types.h>
typedef ssize_t starpu_ssize_t;
#define __starpu_func__ __func__
#endif

#if defined(c_plusplus) || defined(__cplusplus)
/* inline is part of C++ */
#  define __starpu_inline inline
#elif defined(_MSC_VER) || defined(__HP_cc)
#  define __starpu_inline __inline
#else
#  define __starpu_inline __inline__
#endif

/* #undef STARPU_QUICK_CHECK */
/* #undef STARPU_LONG_CHECK */
#define STARPU_USE_DRAND48 1
#define STARPU_USE_ERAND48_R 1
#define STARPU_HAVE_NEARBYINTF 1
#define STARPU_HAVE_RINTF 1

#define STARPU_HAVE_HWLOC 1
#define STARPU_HAVE_PTHREAD_SPIN_LOCK 1
#define STARPU_HAVE_PTHREAD_BARRIER 1
#define STARPU_HAVE_PTHREAD_SETNAME_NP 1
#define STARPU_HAVE_STRUCT_TIMESPEC 1

/** This is only for building examples */
/* #undef STARPU_HAVE_HELGRIND_H */

/** Enable Fortran to C MPI interface */
/* #undef HAVE_MPI_COMM_F2C */

/* #undef STARPU_HAVE_DARWIN */

#define STARPU_HAVE_CXX11 1
#define STARPU_HAVE_STRERROR_R 1
#define STARPU_HAVE_STATEMENT_EXPRESSIONS 1
#define STARPU_PERF_MODEL_DIR "/home/nagp/exageostat-dev/dependinces_installations"

#endif
