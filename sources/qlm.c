/* support for "latmat" - matrices (Lat*internal)*ivec
   - eigenspaces
   - multigrid bases
 */
#include "qlm.h"
#include <assert.h>

static qlmDomain
qlm_domain(qlmLatField t) 
{
    if (QLM_LATREAL == t) return QLM_REAL; 
    else return QLM_CPLX;
}
static size_t 
qlm_real_size(qlmPrec p)
{
    switch (p) {
    case QLM_PREC_FLOAT  :   return sizeof(float);
    case QLM_PREC_DOUBLE :   return sizeof(double);
    }
    return 0;
}
static size_t 
qlm_num_reals(qlmDomain d)
{
    switch (d) {
    case QLM_REAL   :   return 1;
    case QLM_CPLX   :   return 2;
    }
    return 0;
}
static size_t 
qlm_site_num_len(const qlmData *d)
{
    int k = 0;
    switch (d->lftype) {
    case QLM_LATREAL        : 
    case QLM_LATCOMPLEX     :
        k = 1; break;
    case QLM_LATCOLVEC      :
        k = d->nc; break;
    case QLM_LATCOLMAT      :
        k = d->nc * d->nc ; break;
    case QLM_LATDIRFERM     :
        k = d->nc * d->ns ; break;
    case QLM_LATDIRPROP     :
        k = d->nc * d->ns ; k = k * k ; break;
    }
    return (d->is_array ? k * d->nf : k);
}

qlmData *
qlm_data_alloc(
        mLattice *S, const int blk_site_dim[],
        qlmLatField lftype, qlmSubset sset, qlmLayout layout, qlmPrec prec, 
        int parity, int is_array, int nf, int ns, int nc,
        int nvec, int nvec_basis)
{
    int prod, i;
    qlmData *d  = malloc(sizeof(qlmData));
    if (NULL == d) 
        return NULL;
    d->S        = S;
    d->lftype   = lftype;
    d->sset     = sset;
    d->layout   = layout;
    d->prec     = prec;
    d->domain   = qlm_domain(lftype);
    d->parity   = parity;

    d->is_array = is_array;
    d->nf       = (is_array ? nf : 1);
    d->ns       = ns;
    d->nc       = nc;

    d->num_size     = qlm_num_reals(d->domain) * qlm_real_size(prec);
    d->site_num_len = qlm_site_num_len(d);
    d->site_size    = d->num_size * d->site_num_len;

    d->vec_site_len = 1;
    for (i = 0; i < S->rank ; i++) {
        d->vec_site_dim[i] = S->dim[i];
        d->vec_site_len *= d->vec_site_dim[i];
    }
    if (QLM_SUBSET_EOPREC == sset) {
        if (0 != d->vec_site_len % 2) {
            free(d);
            qlm_error = QLM_ERR_GEOM_MISMATCH;
            return NULL;
        }
        d->vec_site_len /= 2;
    }
    d->vec_num_len  = d->site_num_len * d->vec_site_len;
    d->vec_size     = d->num_size * d->vec_num_len;

    d->nvec         = nvec;
    d->nvec_basis   = nvec_basis;

    if (QLM_LAYOUT_QDP == layout) {
        assert(nvec_basis == nvec);
        for (i = 0; i < S->rank ; i++) {
            d->blk_site_dim[i]  = 0;
            d->vec_blk_dim[i]   = 0;
        }
        d->blk_site_len = 0;
        d->vec_blk_len  = 0;
    }
    else if (QLM_LAYOUT_BLOCK == layout) {
        if (NULL == blk_site_dim) {
            free(d);
            qlm_error = QLM_ERR_GEOM_MISMATCH;
            return NULL;
        }
        assert(nvec_basis <= nvec);
        d->blk_site_len = 1;
        d->vec_blk_len  = 1;
        for (i = 0; i < S->rank ; i++) {
            d->blk_site_dim[i] = blk_site_dim[i];
            d->blk_site_len *= d->blk_site_dim[i];
            /* subdivide vecdim into blkdim */
            if (0 != d->vec_site_dim[i] % d->blk_site_dim[i]) {
                free(d);
                qlm_error = QLM_ERR_GEOM_MISMATCH;
                return NULL;
            }
            d->vec_blk_dim[i] = d->vec_site_dim[i] / d->blk_site_dim[i];
            d->vec_blk_len  *= d->vec_blk_dim[i];
        }
        if (QLM_SUBSET_EOPREC == sset) {
            if (0 != d->blk_site_len % 2) {
                free(d);
                qlm_error = QLM_ERR_GEOM_MISMATCH;
                return NULL;
            }
            d->blk_site_len /= 2;
        }
    } 
    d->blk_num_len  = d->site_num_len * d->blk_site_len;
    d->blk_size     = d->num_size * d->blk_num_len;

    d->vec_data     = malloc(d->nvec_basis * d->vec_size);
    if (NULL == d->vec_data) {
        free(d); 
        qlm_error = QLM_ERR_ENOMEM;
        return NULL;
    }
    if (QLM_LAYOUT_BLOCK == layout) {
        d->blk_coeff    = malloc(d->nvec * d->nvec_basis * d->vec_blk_len * d->num_size);
        if (NULL == d->blk_coeff) {
            free(d->vec_data);
            free(d);
            qlm_error = QLM_ERR_ENOMEM;
            return NULL;
        }
    }
    return d;
}


void qlm_data_free(qlmData *qlm)
{
    if (NULL == qlm) return;
    if (NULL != qlm->vec_data)  free(qlm->vec_data);
    if (NULL != qlm->blk_coeff) free(qlm->blk_coeff);
}



#define QT(a)    a##Real
#define QA(a)    a##R
#define QAx(a,b) a##R##b
#include "qlm-x.c"
#define QT(a)    a##Complex
#define QA(a)    a##C
#define QAx(a,b) a##C##b
#include "qlm-x.c"
#define QT(a)    a##ColorVector
#define QA(a)    a##V
#define QAx(a,b) a##V##b
#include "qlm-x.c"
#define QT(a)    a##ColorMatrix
#define QA(a)    a##M
#define QAx(a,b) a##M##b
#include "qlm-x.c"
#define QT(a)    a##DiracFermion
#define QA(a)    a##D
#define QAx(a,b) a##D##b
#include "qlm-x.c"
#define QT(a)    a##DiracPropagator
#define QA(a)    a##P
#define QAx(a,b) a##P##b
#include "qlm-x.c"
