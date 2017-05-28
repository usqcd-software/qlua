#ifndef QLM_H__DAB3D8AA_24DB_4C80_8268_62DEB1688C8E
#define QLM_H__DAB3D8AA_24DB_4C80_8268_62DEB1688C8E

#include "modules.h"                                                 /* DEPS */
#include "qlua.h"                                                    /* DEPS */
#include "qcomplex.h"                                                /* DEPS */
#include "qvector.h"                                                 /* DEPS */
#include "lattice.h"                                                 /* DEPS */
#include "qlayout.h"                                                 /* DEPS */
#include "latreal.h"                                                 /* DEPS */
#include "latcomplex.h"                                              /* DEPS */
#include "latcolmat.h"                                               /* DEPS */
#include "latcolvec.h"                                               /* DEPS */
#include "latdirferm.h"                                              /* DEPS */
#include "latdirprop.h"                                              /* DEPS */
#include "crc32.h"                                                   /* DEPS */
#include "qend.h"                                                    /* DEPS */
#include "qlanczos.h"                                                /* DEPS */
#include "qmp.h"

#include <string.h>
#include <math.h>
#include <sys/stat.h>
#include <unistd.h>
#include <errno.h>
#include <time.h>


/* XXX NOTES
    * some enum values are defined explicitly to be able to change 
      aux functions (len/size/etc) into define's and have them inlined 
    * to apply eo-prec matvec, double conversion is required: block<->qdp,qdp<->internal
 */

/* internal storage */
typedef enum {
    QLM_ERR_SUCCESS         = 0,
    QLM_ERR_ENOMEM          = 1,        /* out of memory (malloc) */
    QLM_ERR_INVALID_VEC     = 2,        /* vec.index out of range */
    QLM_ERR_INVALID_FIELD   = 3,        /* on export in C functions */
    QLM_ERR_INVALID_QTYPE   = 4,        /* on import in C, Qlua functions */
    QLM_ERR_TYPE_MISMATCH   = 5,        /* mismatch of metadata in func with >=2 qlmData's */
    QLM_ERR_GEOM_MISMATCH   = 6         /* mismatch of metadata in func with >=2 qlmData's */
} qlmError;
const char *qlm_strerror(qlmError e);
extern qlmError qlm_error;  /* TODO init to QLM_ERR_SUCCESS */

typedef enum {
    QLM_REAL    = 1,
    QLM_CPLX    = 2
} qlmDomain;
/* FIXME use existing enum? */
typedef enum {
    QLM_LATREAL         =  1,
    QLM_LATCOMPLEX      =  3,
    QLM_LATCOLVEC       =  5,
    QLM_LATCOLMAT       =  7,
    QLM_LATDIRFERM      = 11,
    QLM_LATDIRPROP      = 13
} qlmLatField;
/* internal layout */
typedef enum {
    QLM_LAYOUT_QDP      = 1,
    QLM_LAYOUT_BLOCK    = 2
} qlmLayout;
typedef enum qlmSubset {
    QLM_SUBSET_FULL     = 1,
    QLM_SUBSET_EOPREC   = 2
} qlmSubset;
typedef enum {
    QLM_PREC_FLOAT   = 4,
    QLM_PREC_DOUBLE  = 8
} qlmPrec;

/* storage of sites: site data is contiguous; complex numbers are contiguous
    * QLM_LAYOUT_QDP: QDP site order, only for QLM_SUBSET_FULL
    * QLM_LAYOUT_BLOCK: blocked layout (mg, evec compression ; blocks must be contiguous for BLAS)
      + local dimensions must be divisible by block
      + if not QLM_SUBSET_FULL, then `blocksize[0]' must be even; otherwise, blocks have different lengths
      + storage of blocks is lexicographic 
      + sites within block are lexicographic
   site data

   overall index for dfa [ivec, iblock, isite, i5, is, ic]
 */

typedef struct {
    mLattice    *S;         /* associated lattice ; note that operations are not subject to mask */
    qlmLatField lftype;
    qlmSubset   sset;
    qlmLayout   layout;
    qlmPrec     prec;
    qlmDomain   domain;     /* deduced from lftype */
    
    /* vector space */
    int         nvec;                               /* number of stored vectors */
    int         nvec_basis;                         /* number of basis vectors [QLM_LAYOUT_BLOCK] 
                                                       QLM_LAYOUT_BLOCK ? <=nvec : ==nvec */
    /* float|double data (~ Fortran matrix) 
       [ivec(m), iblock(b), {isiteb,iarr,is^?,ic^?}(q)]
        ivec    (basis) vector index
        iblock  lex-linear block index                      [QLM_LAYOUT_BLOCK]
        isiteb  lex-linear (prec)site index within block    [QLM_LAYOUT_BLOCK]
        iarr    array index(ie 5th dim for DWF)
        is, ic  spin/color indices
     */
    void        *vec_data;                          /* size=nvec_basis * vec_size */
    /* vector coeffs. in block space 
        [jvec, ]
       NULL if not blocking 
     */
    void        *blk_coeff;                         /* block proj.coeff. data [QLM_LAYOUT_BLOCK] */

    /* latfield and layout-specific parameters */
    int         is_array;
    int         nf, ns, nc;                         /* color, spin, number of fields [QLM_LAT*_A] 
                                                       for import/export */
    int         x0[QLUA_MAX_LATTICE_RANK];          /* initial coordinate of local vec */
    int         parity,                             /* parity [[!QLM_SUBSET_FULL] */ 
                x0_parity;                          /* parity of initial coordinate [!QLM_SUBSET_FULL] */


    /* "private" geometry vars */
    /* local vol dimensions and length (sites) */
    int         vec_site_dim[QLUA_MAX_LATTICE_RANK],
                vec_site_len;       /* == QLM_SUBSET_FULL ? prod(vec_site_dim) 
                                                          : prod(vec_site_dim)/2 */
    /* local vol dimensions and length (blocks) */
    int         vec_blk_dim[QLUA_MAX_LATTICE_RANK],
                vec_blk_len;        /* vec_block_len = prod(vec_block_dim) */
    /* block dimensions and length(site data records) [QLM_LAYOUT_BLOCK] */
    int         blk_site_dim[QLUA_MAX_LATTICE_RANK],
                blk_site_len;       /* == QLM_SUBSET_FULL ? prod(blk_site_dim) 
                                                          : prod(blk_site_dim)/2 */
    int         vec_num_len,        /* == site_num_len * vec_site_len */
                blk_num_len,        /* == site_num_len * blk_site_len */
                site_num_len;       /* == nf * ns^? * nc^? */
    int         vec_size,           /* == num_size * vec_num_len */
                blk_size,           /* == num_size * blk_num_len */
                site_size,          /* == num_size * site_num_len */
                num_size;           /* == (QLM_REAL ? 1 : 2) * (QLM_PREC_FLOAT:4,QLM_PREC_DOUBLE:8) */
} qlmData;

qlmData *qlm_data_alloc(
        mLattice *S, const int blk_site_dim[],
        qlmLatField lftype, qlmSubset sset, qlmLayout layout, qlmPrec p, 
        int parity, int is_array, int nf, int ns, int nc,
        int nvec, int nvec_basis);
void qlm_data_free(qlmData *qlm);

/* QDP type primitives */
QDP_D_Real **           qlm_vector_alloc_R(int n);
QDP_D_Complex **        qlm_vector_alloc_C(int n);
QDP_D_ColorVector **    qlm_vector_alloc_V(int n);
QDP_D_ColorMatrix **    qlm_vector_alloc_M(int n);
QDP_D_DiracFermion **   qlm_vector_alloc_D(int n);
QDP_D_DiracPropagator **qlm_vector_alloc_P(int n);
void qlm_vector_free_R(QDP_D_Real **,            int n);
void qlm_vector_free_C(QDP_D_Complex **,         int n);
void qlm_vector_free_V(QDP_D_ColorVector **,     int n);
void qlm_vector_free_M(QDP_D_ColorMatrix **,     int n);
void qlm_vector_free_D(QDP_D_DiracFermion **,    int n);
void qlm_vector_free_P(QDP_D_DiracPropagator **, int n);
QLA_D_Real **           qlm_vector_expose_R(QDP_D_Real **,            int n);
QLA_D_Complex **        qlm_vector_expose_C(QDP_D_Complex **,         int n);
QLA_D_ColorVector **    qlm_vector_expose_V(QDP_D_ColorVector **,     int n);
QLA_D_ColorMatrix **    qlm_vector_expose_M(QDP_D_ColorMatrix **,     int n);
QLA_D_DiracFermion **   qlm_vector_expose_D(QDP_D_DiracFermion **,    int n);
QLA_D_DiracPropagator **qlm_vector_expose_P(QDP_D_DiracPropagator **, int n);
void qlm_vector_reset_R(QDP_D_Real **,            int n);
void qlm_vector_reset_C(QDP_D_Complex **,         int n);
void qlm_vector_reset_V(QDP_D_ColorVector **,     int n);
void qlm_vector_reset_M(QDP_D_ColorMatrix **,     int n);
void qlm_vector_reset_D(QDP_D_DiracFermion **,    int n);
void qlm_vector_reset_P(QDP_D_DiracPropagator **, int n);


/* C-export functions : export from memory into pre-allocated lattice objects */
int qlm_export_vector_R(QDP_D_Real **dst,            const qlmData *d, const void *x);
int qlm_export_vector_C(QDP_D_Complex **dst,         const qlmData *d, const void *x);
int qlm_export_vector_V(QDP_D_ColorVector **dst,     const qlmData *d, const void *x);
int qlm_export_vector_M(QDP_D_ColorMatrix **dst,     const qlmData *d, const void *x);
int qlm_export_vector_D(QDP_D_DiracFermion **dst,    const qlmData *d, const void *x);
int qlm_export_vector_P(QDP_D_DiracPropagator **dst, const qlmData *d, const void *x);
/* C-extract functions : extract from vector space */
int qlm_extract_vector_R(QDP_D_Real **dst,            const qlmData *d, int ivec);
int qlm_extract_vector_C(QDP_D_Complex **dst,         const qlmData *d, int ivec);
int qlm_extract_vector_V(QDP_D_ColorVector **dst,     const qlmData *d, int ivec);
int qlm_extract_vector_M(QDP_D_ColorMatrix **dst,     const qlmData *d, int ivec);
int qlm_extract_vector_D(QDP_D_DiracFermion **dst,    const qlmData *d, int ivec);
int qlm_extract_vector_P(QDP_D_DiracPropagator **dst, const qlmData *d, int ivec);
/* Qlua-export function: push table onto stack */
int qlm_extract_vector_qlua(lua_State *L, qlmData *d, int ivec);
/* C-import functions */
int qlm_import_vector_R(const qlmData *d, void *x, QDP_D_Real **src);
int qlm_import_vector_C(const qlmData *d, void *x, QDP_D_Complex **src);
int qlm_import_vector_V(const qlmData *d, void *x, QDP_D_ColorVector **src);
int qlm_import_vector_M(const qlmData *d, void *x, QDP_D_ColorMatrix **src);
int qlm_import_vector_D(const qlmData *d, void *x, QDP_D_DiracFermion **src);
int qlm_import_vector_P(const qlmData *d, void *x, QDP_D_DiracPropagator **src);
int qlm_insert_vector_R(qlmData *d, int ivec, QDP_D_Real **src);
int qlm_insert_vector_C(qlmData *d, int ivec, QDP_D_Complex **src);
int qlm_insert_vector_V(qlmData *d, int ivec, QDP_D_ColorVector **src);
int qlm_insert_vector_M(qlmData *d, int ivec, QDP_D_ColorMatrix **src);
int qlm_insert_vector_D(qlmData *d, int ivec, QDP_D_DiracFermion **src);
int qlm_insert_vector_P(qlmData *d, int ivec, QDP_D_DiracPropagator **src);
/* Qlua-import function: get table from stack */
int qlm_insert_vector_qlua(lua_State *L, qlmData *d, int ivec);




#endif/*QLM_H__DAB3D8AA_24DB_4C80_8268_62DEB1688C8E*/

