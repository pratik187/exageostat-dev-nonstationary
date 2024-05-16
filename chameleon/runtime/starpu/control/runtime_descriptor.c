/**
 *
 * @file starpu/runtime_descriptor.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2021 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon StarPU descriptor routines
 *
 * @version 1.1.0
 * @author Cedric Augonnet
 * @author Mathieu Faverge
 * @author Cedric Castagnede
 * @author Florent Pruvost
 * @author Guillaume Sylvand
 * @author Raphael Boucherie
 * @author Samuel Thibault
 * @date 2020-12-01
 *
 */
#include <stdlib.h>
#include <unistd.h>
#include "chameleon_starpu.h"

/**
 *  Set the tag sizes
 */
#if defined(CHAMELEON_USE_MPI)

/* Take 24 bits for the tile id, and 7 bits for descriptor id.
 These values can be changed through the call CHAMELEON_user_tag_size(int tag_width, int tag_sep) */
#define TAG_WIDTH_MIN 20
static int tag_width = 64;
static int tag_sep   = 40;
static int _tag_mpi_initialized_ = 0;

static inline int
chameleon_starpu_tag_init( int user_tag_width,
                           int user_tag_sep )
{
    if (!_tag_mpi_initialized_) {
        int ok = 0;
        uintptr_t tag_ub;

        tag_width = user_tag_width;
        tag_sep   = user_tag_sep;

        void *tag_ub_p = NULL;
#if defined(HAVE_STARPU_MPI_COMM_GET_ATTR)
        starpu_mpi_comm_get_attr(MPI_COMM_WORLD, STARPU_MPI_TAG_UB, &tag_ub_p, &ok);
#else
        MPI_Comm_get_attr(MPI_COMM_WORLD, MPI_TAG_UB, &tag_ub_p, &ok);
#endif
        tag_ub = (uintptr_t)tag_ub_p;

        if ( !ok ) {
            chameleon_error("RUNTIME_desc_create", "MPI_TAG_UB not known by StarPU");
        }

        while ( ((uintptr_t)((1UL<<tag_width) - 1) > tag_ub ) &&
                (tag_width >= TAG_WIDTH_MIN) )
        {
            tag_width--;
            tag_sep--;
        }

        if ( tag_width < TAG_WIDTH_MIN ) {
            chameleon_error("RUNTIME_desc_create", "MPI_TAG_UB too small to identify all the data");
            return CHAMELEON_ERR_OUT_OF_RESOURCES;
        }

        _tag_mpi_initialized_ = 1;
        return CHAMELEON_SUCCESS;
    }
    else {
        return CHAMELEON_ERR_REINITIALIZED;
    }
}


#ifndef HAVE_STARPU_MPI_DATA_REGISTER
#define starpu_mpi_data_register( handle_, tag_, owner_ )       \
    do {                                                        \
        starpu_data_set_rank( (handle_), (owner_) );            \
        starpu_data_set_tag( (handle_), (tag_) );               \
    } while(0)
#endif

#endif

void RUNTIME_comm_set_tag_sizes( int user_tag_width,
                                 int user_tag_sep )
{
#if defined(CHAMELEON_USE_MPI)
    int rc;
    rc = chameleon_starpu_tag_init( user_tag_width, user_tag_sep );
    if ( rc != CHAMELEON_SUCCESS ) {
        chameleon_error("RUNTIME_user_tag_size",
                    "must be called before creating any Cham descriptor with CHAMELEON_Desc_create(). The tag sizes will not be modified.");
    }
#endif
    (void)user_tag_width; (void)user_tag_sep;
}

/**
 *  Malloc/Free of the data
 */
#ifdef STARPU_MALLOC_SIMULATION_FOLDED
#define FOLDED STARPU_MALLOC_SIMULATION_FOLDED
#else
#define FOLDED 0
#endif

void *RUNTIME_malloc( size_t size )
{
#if defined(CHAMELEON_SIMULATION) && !defined(STARPU_MALLOC_SIMULATION_FOLDED) && !defined(CHAMELEON_USE_MPI)
    return (void*) 1;
#else
    void *ptr;

    if (starpu_malloc_flags(&ptr, size, STARPU_MALLOC_PINNED|FOLDED|STARPU_MALLOC_COUNT) != 0) {
        return NULL;
    }
    return ptr;
#endif
}

void RUNTIME_free( void  *ptr,
                   size_t size )
{
#if defined(CHAMELEON_SIMULATION) && !defined(STARPU_MALLOC_SIMULATION_FOLDED) && !defined(CHAMELEON_USE_MPI)
    (void)ptr; (void)size;
    return;
#else
    starpu_free_flags(ptr, size, STARPU_MALLOC_PINNED|FOLDED|STARPU_MALLOC_COUNT);
#endif
}

/**
 *  Create data descriptor
 */
void RUNTIME_desc_create( CHAM_desc_t *desc )
{
    int64_t lmt = desc->lmt;
    int64_t lnt = desc->lnt;

    desc->occurences = 1;

    /*
     * Allocate starpu_handle_t array (handlers are initialized on the fly when
     * discovered by any algorithm to save space)
     */
    desc->schedopt = (void*)calloc(lnt*lmt,sizeof(starpu_data_handle_t));
    assert(desc->schedopt);

#if defined(CHAMELEON_USE_CUDA) && !defined(CHAMELEON_SIMULATION)
    /*
     * Register allocated memory as CUDA pinned memory
     */
    if ( (desc->use_mat == 1) && (desc->register_mat == 1) )
    {
        int64_t eltsze = CHAMELEON_Element_Size(desc->dtyp);
        size_t size = (size_t)(desc->llm) * (size_t)(desc->lln) * eltsze;
        cudaError_t rc;

        /* Register the matrix as pinned memory */
        rc = cudaHostRegister( desc->mat, size, cudaHostRegisterPortable );
        if ( rc != cudaSuccess )
        {
            /* Disable the unregister as register failed */
            desc->register_mat = 0;
            chameleon_warning("RUNTIME_desc_create(StarPU): cudaHostRegister - ", cudaGetErrorString( rc ));
        }
    }
#endif

    if (desc->ooc) {
        char   *backend = getenv("STARPU_DISK_SWAP_BACKEND");

        if (backend && strcmp(backend, "unistd_o_direct") == 0) {
            int     lastmm   = desc->lm - (desc->lmt-1) * desc->mb;
            int     lastnn   = desc->ln - (desc->lnt-1) * desc->nb;
            int64_t eltsze   = CHAMELEON_Element_Size(desc->dtyp);
            int     pagesize = getpagesize();

            if ( ((desc->mb * desc->nb * eltsze) % pagesize != 0) ||
                 ((lastmm   * desc->nb * eltsze) % pagesize != 0) ||
                 ((desc->mb * lastnn   * eltsze) % pagesize != 0) ||
                 ((lastmm   * lastnn   * eltsze) % pagesize != 0) )
            {
                chameleon_error("RUNTIME_desc_create", "Matrix and tile size not suitable for out-of-core: all tiles have to be multiples of the system page size. Tip : choose 'n' and 'nb' as both multiples of 32.");
                return;
            }
        }
    }

#if defined(CHAMELEON_USE_MPI)
    /*
     * Check that we are not going over MPI tag limitations
     */
    {
        chameleon_starpu_tag_init( tag_width, tag_sep );

        /* Check that we won't create overflow in tags used */
        if ( ((uintptr_t)(lnt*lmt)) > ((uintptr_t)(1UL<<tag_sep)) ) {
            chameleon_fatal_error("RUNTIME_desc_create", "Too many tiles in the descriptor for MPI tags");
            return;
        }
        assert( (lmt*lnt) <= (1UL<<tag_sep) );

        if ( ((uintptr_t)desc->id) >= (uintptr_t)(1UL<<(tag_width-tag_sep)) ) {
            chameleon_fatal_error("RUNTIME_desc_create", "Number of descriptor available in MPI mode out of stock");
            return;
        }
        assert( ((uintptr_t)desc->id) < (uintptr_t)(1UL<<(tag_width-tag_sep)) );
    }
#endif
}

/**
 *  Destroy data descriptor
 */
void RUNTIME_desc_destroy( CHAM_desc_t *desc )
{
    desc->occurences--;

    /*
     * If this is the last descriptor using the matrix, we release the handle
     * and unregister the GPU data
     */
    if ( desc->occurences == 0 ) {
        starpu_data_handle_t *handle = (starpu_data_handle_t*)(desc->schedopt);
        int lmt = desc->lmt;
        int lnt = desc->lnt;
        int m, n;

        for (n = 0; n < lnt; n++) {
            for (m = 0; m < lmt; m++)
            {
                if (*handle != NULL) {
                    starpu_data_unregister(*handle);
                    *handle = NULL;
                }
                handle++;
            }
        }

#if defined(CHAMELEON_USE_CUDA) && !defined(CHAMELEON_SIMULATION)
        if ( (desc->use_mat == 1) && (desc->register_mat == 1) )
        {
            /* Unmap the pinned memory associated to the matrix */
            if (cudaHostUnregister(desc->mat) != cudaSuccess)
            {
                chameleon_warning("RUNTIME_desc_destroy(StarPU)",
                                  "cudaHostUnregister failed to unregister the "
                                  "pinned memory associated to the matrix");
            }
        }
#endif /* defined(CHAMELEON_USE_CUDA) */

        free(desc->schedopt);
    }
}

/**
 *  Acquire data
 */
int RUNTIME_desc_acquire( const CHAM_desc_t *desc )
{
    starpu_data_handle_t *handle = (starpu_data_handle_t*)(desc->schedopt);
    int lmt = desc->lmt;
    int lnt = desc->lnt;
    int m, n;

    for (n = 0; n < lnt; n++) {
        for (m = 0; m < lmt; m++)
        {
            if ( (*handle == NULL) ||
                 !chameleon_desc_islocal( desc, m, n ) )
            {
                handle++;
                continue;
            }
            starpu_data_acquire(*handle, STARPU_R);
            handle++;
        }
    }
    return CHAMELEON_SUCCESS;
}

/**
 *  Release data
 */
int RUNTIME_desc_release( const CHAM_desc_t *desc )
{
    starpu_data_handle_t *handle = (starpu_data_handle_t*)(desc->schedopt);
    int lmt = desc->lmt;
    int lnt = desc->lnt;
    int m, n;

    for (n = 0; n < lnt; n++) {
        for (m = 0; m < lmt; m++)
        {
            if ( (*handle == NULL) ||
                 !chameleon_desc_islocal( desc, m, n ) )
            {
                handle++;
                continue;
            }
            starpu_data_release(*handle);
            handle++;
        }
    }
    return CHAMELEON_SUCCESS;
}

/**
 *  Flush cached data
 */
void RUNTIME_flush()
{
#if defined(CHAMELEON_USE_MPI)
    starpu_mpi_cache_flush_all_data(MPI_COMM_WORLD);
#endif
}

/**
 * Different implementations of the flush call based on StarPU version
 */
#ifdef HAVE_STARPU_DATA_WONT_USE

static inline void
chameleon_starpu_data_wont_use( starpu_data_handle_t handle ) {
    starpu_data_wont_use( handle );
}

#elif defined HAVE_STARPU_IDLE_PREFETCH

static inline void
chameleon_starpu_data_flush( starpu_data_handle_t handle)
{
    starpu_data_idle_prefetch_on_node(handle, STARPU_MAIN_RAM, 1);
    starpu_data_release_on_node(handle, -1);
}

static inline void
chameleon_starpu_data_wont_use( starpu_data_handle_t handle ) {
    starpu_data_acquire_on_node_cb( handle, -1, STARPU_R,
                                    chameleon_starpu_data_flush, handle );
}

#else

static inline void
chameleon_starpu_data_wont_use( starpu_data_handle_t handle ) {
    starpu_data_acquire_cb( handle, STARPU_R,
                            (void (*)(void*))&starpu_data_release, handle );
}

#endif

void RUNTIME_desc_flush( const CHAM_desc_t        *desc,
                         const RUNTIME_sequence_t *sequence )
{
    int mt = desc->mt;
    int nt = desc->nt;
    int m, n;

    for (n = 0; n < nt; n++)
    {
        for (m = 0; m < mt; m++)
        {
            RUNTIME_data_flush( sequence, desc, m, n );
        }
    }
}

void RUNTIME_data_flush( const RUNTIME_sequence_t *sequence,
                         const CHAM_desc_t *A, int m, int n )
{
    int64_t mm = m + (A->i / A->mb);
    int64_t nn = n + (A->j / A->nb);

    starpu_data_handle_t *handle = A->schedopt;
    handle += ((int64_t)A->lmt) * nn + mm;

    if (*handle == NULL) {
        return;
    }

#if defined(CHAMELEON_USE_MPI)
    starpu_mpi_cache_flush( MPI_COMM_WORLD, *handle );
#endif

    if ( chameleon_desc_islocal( A, m, n ) ) {
        chameleon_starpu_data_wont_use( *handle );
    }

    (void)sequence;
}

#if defined(CHAMELEON_USE_MIGRATE)
void RUNTIME_data_migrate( const RUNTIME_sequence_t *sequence,
                           const CHAM_desc_t *A, int Am, int An, int new_rank )
{
#if defined(HAVE_STARPU_MPI_DATA_MIGRATE)
    int old_rank;
    starpu_data_handle_t *handle = (starpu_data_handle_t*)(A->schedopt);
    starpu_data_handle_t lhandle;
    handle += ((int64_t)(A->lmt) * (int64_t)An + (int64_t)Am);

    lhandle = *handle;
    if ( lhandle == NULL ) {
        /* Register the data */
        lhandle = RUNTIME_data_getaddr( A, Am, An );
    }
    old_rank = starpu_mpi_data_get_rank( lhandle );

    if ( old_rank != new_rank ) {
        starpu_mpi_data_migrate( MPI_COMM_WORLD, lhandle, new_rank );
    }

    (void)sequence;
#else
    (void)sequence; (void)A; (void)Am; (void)An; (void)new_rank;
#endif
}
#endif

/**
 *  Get data addr
 */
/* For older revision of StarPU, STARPU_MAIN_RAM is not defined */
#ifndef STARPU_MAIN_RAM
#define STARPU_MAIN_RAM 0
#endif

void *RUNTIME_data_getaddr( const CHAM_desc_t *A, int m, int n )
{
    int64_t mm = m + (A->i / A->mb);
    int64_t nn = n + (A->j / A->nb);

    starpu_data_handle_t *ptrtile = A->schedopt;
    ptrtile += ((int64_t)A->lmt) * nn + mm;

    if (*ptrtile == NULL) {
        int home_node = -1;
        int myrank = A->myrank;
        int owner  = A->get_rankof( A, m, n );
        CHAM_tile_t *tile = A->get_blktile( A, m, n );

        if ( myrank == owner ) {
            if ( tile->mat != NULL )
            {
                home_node = STARPU_MAIN_RAM;
            }
        }

        starpu_cham_tile_register( ptrtile, home_node, tile, A->dtyp );

#if defined(HAVE_STARPU_DATA_SET_OOC_FLAG)
        if ( A->ooc == 0 ) {
            starpu_data_set_ooc_flag( *ptrtile, 0 );
        }
#endif

#if defined(HAVE_STARPU_DATA_SET_COORDINATES)
        starpu_data_set_coordinates( *ptrtile, 2, m, n );
#endif

#if defined(CHAMELEON_USE_MPI)
        {
            int64_t block_ind = A->lmt * nn + mm;
            starpu_mpi_data_register(*ptrtile, (((int64_t)A->id) << tag_sep) | block_ind, owner);
        }
#endif /* defined(CHAMELEON_USE_MPI) */
    }

    return *ptrtile;
}
