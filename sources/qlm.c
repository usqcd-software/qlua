/* support for "latmat" - matrices (Lat*internal)*ivec
   - eigenspaces
   - multigrid bases
 */
#include "qlm.h"
#include <assert.h>
#include <string.h>
#include <stdlib.h>
#include "qdp_d2.h"
#include "qdp_d3.h"
#include "qdp_dn.h"

#include "modules.h"                                                 /* DEPS */
#include "qlua.h"                                                    /* DEPS */
#include "qcomplex.h"                                                /* DEPS */
#include "qvector.h"                                                 /* DEPS */
#include "lattice.h"                                                 /* DEPS */
#include "qlayout.h"                                                 /* DEPS */

#include "latreal.h"                                                 /* DEPS */
#include "latcomplex.h"                                              /* DEPS */
#include "latcolvec.h"                                               /* DEPS */
#include "latcolmat.h"                                               /* DEPS */
#include "latdirferm.h"                                              /* DEPS */
#include "latdirprop.h"                                              /* DEPS */

#define NORMALIZE_INDEX(L, idx) do if (idx < 0) idx = lua_gettop(L) + 1 + idx; while (0)

qlmError qlm_error;
const char *
qlm_strerror(qlmError e) 
{
    switch (e) {
    case QLM_ERR_SUCCESS        : return "success" ;
    case QLM_ERR_ENOMEM         : return "not enough memory" ;
    case QLM_ERR_INVALID_VEC    : return "invalid vector index" ;
    case QLM_ERR_INVALID_FIELD  : return "incorrect lattice field" ;
    case QLM_ERR_TYPE_MISMATCH  : return "mismatched metadata" ;
    case QLM_ERR_GEOM_MISMATCH  : return "mismatched geometry" ;
    case QLM_ERR_ENUM_ERROR     : return "unsupported constant value (internal error)";
    case QLM_ERR_GEOM_SITE      : return "operation on an incorrect site";
    }
    /* should not get here : handle all cases above */
    return NULL;
}

static size_t 
qlm_field_num_len(const qlmData *d)
{
    switch (d->lftype) {
    case QLM_LATREAL        : 
    case QLM_LATCOMPLEX     :   return 1;
    case QLM_LATCOLVEC      :   return d->nc;
    case QLM_LATCOLMAT      :   return d->nc * d->nc ;
    case QLM_LATDIRFERM     :   return d->nc * d->ns ;
    case QLM_LATDIRPROP     :   return (d->nc * d->nc) * (d->ns * d->ns); 
    case QLM_LAT_NONE       :   return 0;
    }
    return 0;
}

static size_t 
qlm_site_num_len(const qlmData *d)
{
    return (d->is_array ? d->arr_len : 1) * qlm_field_num_len(d);
}
const char *
qlm_lftype_str(qlmLatField lftype) {
    switch(lftype) {
    case QLM_LATREAL    : return "Real"           ; break ;
    case QLM_LATCOMPLEX : return "Complex"        ; break ;
    case QLM_LATCOLVEC  : return "ColorVector"    ; break ;
    case QLM_LATCOLMAT  : return "ColorMatrix"    ; break ;
    case QLM_LATDIRFERM : return "DiracFermion"   ; break ;
    case QLM_LATDIRPROP : return "DiracPropagator"; break ;
    default : return "<none>";
    }
}
qlmLatField 
qlm_str2lftype(const char *str) 
{
    if      (!strcmp("Real"           , str)) return QLM_LATREAL     ;
    else if (!strcmp("Complex"        , str)) return QLM_LATCOMPLEX  ;
    else if (!strcmp("ColorVector"    , str)) return QLM_LATCOLVEC   ;
    else if (!strcmp("ColorMatrix"    , str)) return QLM_LATCOLMAT   ;
    else if (!strcmp("DiracFermion"   , str)) return QLM_LATDIRFERM  ;
    else if (!strcmp("DiracPropagator", str)) return QLM_LATDIRPROP  ;
    else {
        qlm_error = QLM_ERR_ENUM_ERROR;
        return QLM_LAT_NONE; 
    }
}
qlmPrec
qlm_str2prec(const char *str)
{
    if      (!strcmp("float",  str)) return QLM_PREC_FLOAT;
    else if (!strcmp("double", str)) return QLM_PREC_DOUBLE;
    else {
        qlm_error = QLM_ERR_ENUM_ERROR;
        return QLM_PREC_NONE;
    }
}
qlmSublat 
qlm_str2sublat(const char *str)
{
    if      (!strcmp("full", str)) return QLM_SUBLAT_FULL;
    else if (!strcmp("even", str)) return QLM_SUBLAT_EVEN;
    else if (!strcmp("odd" , str)) return QLM_SUBLAT_ODD;
    else {
        qlm_error = QLM_SUBLAT_NONE;
        return QLM_PREC_NONE;
    }
}

qlmData *
qlm_data_alloc(
        mLattice *S, 
        qlmSublat sublat, 
        int nvec, 
        qlmLatField lftype, 
        int is_array, 
        int arr_len, 
        int ns, 
        int nc,
        qlmPrec prec, 
        int is_blocked, 
        const int blk_site_dim[], 
        int nvec_basis)
{
    qlmData *d  = NULL;
    d           = malloc(sizeof(qlmData));
    if (NULL == d) 
        return NULL;
    d->S        = S;
    d->lftype   = lftype;
    d->sublat   = sublat;
    d->prec     = prec;
    d->domain   = qlm_domain(lftype);

    d->is_array = is_array;
    d->arr_len  = (is_array ? arr_len : 1);
    d->ns       = ns;
    d->nc       = nc;

    d->num_size     = qlm_num_reals(d->domain) * qlm_real_size(prec);
    d->site_num_len = qlm_site_num_len(d);
    d->site_size    = d->num_size * d->site_num_len;

    d->ndim     = S->rank;
    int xlo[QLUA_MAX_LATTICE_RANK], 
        xhi[QLUA_MAX_LATTICE_RANK];
    qlua_sublattice(xlo, xhi, QDP_this_node, (void *)S);
    d->x0_parity = 0;
    d->vec_site_len = 1;
    for (int i = 0; i < S->rank ; i++) {
        d->x0[i] = xlo[i];
        d->x0_parity ^= (d->x0[i] & 1);
        d->vec_site_dim[i] = xhi[i] - xlo[i];
        d->vec_site_len *= d->vec_site_dim[i];
    }
    if (IS_SUBLAT_EOPC(sublat)) {
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
    d->is_blocked   = is_blocked;
    d->nvec_basis   = nvec_basis;

    if (!is_blocked) {
        assert(nvec_basis == nvec);
        /* safe values: exactly 1 block */
        for (int i = 0; i < S->rank ; i++) {
            d->blk_site_dim[i]  = d->vec_site_dim[i];
            d->vec_blk_dim[i]   = 1;
        }
        d->blk_site_len = d->vec_site_len;
        d->vec_blk_len  = 1;
    }
    else {
        if (NULL == blk_site_dim) {
            free(d);
            qlm_error = QLM_ERR_GEOM_MISMATCH;
            return NULL;
        }
        assert(nvec_basis <= nvec);
        d->blk_site_len = 1;
        d->vec_blk_len  = 1;
        for (int i = 0; i < S->rank ; i++) {
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
        if (IS_SUBLAT_EOPC(d->sublat)) {
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

    int vec_data_size = d->nvec_basis * d->vec_size;
    d->vec_data     = malloc(vec_data_size);
    if (NULL == d->vec_data) {
        free(d); 
        qlm_error = QLM_ERR_ENOMEM;
        return NULL;
    }
    memset(d->vec_data, 0, vec_data_size);

    if (is_blocked) {
        int blk_coeff_size = d->nvec * d->nvec_basis * d->vec_blk_len * d->num_size;
        d->blk_coeff    = malloc(blk_coeff_size);
        if (NULL == d->blk_coeff) {
            free(d->vec_data);
            free(d);
            qlm_error = QLM_ERR_ENOMEM;
            return NULL;
        }
        memset(d->blk_coeff, 0, blk_coeff_size);
    } else
        d->blk_coeff = NULL;


    return d;
}


void qlm_data_free(qlmData *qlm)
{
    if (NULL == qlm) return;
    if (NULL != qlm->vec_data)  free(qlm->vec_data);
    if (NULL != qlm->blk_coeff) free(qlm->blk_coeff);
    free(qlm);
}


/* TODO move to library */
/* print array of "type"s separated by delimeters
   XXX USE WITH CARE - MOSTLY FOR DEBUGGING */
#define len_arr(a)   (sizeof(a)/sizeof(a[0]))
#define snprintf_arr(str_, size_, format, delim, type, ptr_, n_) \
do { \
    int size = size_, dlen, delim_len = strlen(delim); \
    type *ptr = ptr_; \
    char *str = str_; \
    for (int n = n_ ; 0 < n && 0 < size ; n--, ptr++) { \
        dlen = snprintf(str, size, format, *ptr); \
        size -= dlen;  str += dlen; \
        if (size <= 0) break; \
        if (1 < n) { \
            if (delim_len < size) { \
                strncpy(str, delim, size); size -= delim_len;  str += delim_len; \
            } else break; \
        }\
    }\
} while(0)


/****************************************************************************/
/* standardize handling of (arrays of) qdp fields */
typedef struct {
    qlmLatField lftype;
    int nc;
    int is_array, arr_len;
    int is_exposed;
    union {
        void *x;
        void **a;
    } p;    /* for lattice fields (QDP ptrs) */
    union {
        void *x;
        void **a;
    } q;    /* for exposed fields (QLA ptrs) */
} qdp_vtype;

static qdp_vtype *
qdp_vtype_alloc(qlmLatField lftype, int is_array, int arr_len, int nc)
{
    qdp_vtype *res = malloc(sizeof(qdp_vtype));
    if (NULL == res) return NULL;
    res->lftype     = lftype;
    res->nc         = nc;
    res->is_array   = is_array;
    res->arr_len    = arr_len;
    res->is_exposed = 0;
    if (is_array) {
        res->p.a = malloc(arr_len * sizeof(void *));
        res->q.a = malloc(arr_len * sizeof(void *));
        if (NULL == res->p.a || NULL == res->q.a) {
            free(res);
            return NULL;
        }
    } else
        res->p.a = res->q.a = NULL;
    return res;
}
/****** expose/reset ******/
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


static void *
qdptype_create(mLattice *S, qlmLatField lftype, int nc)
{ 
    switch(lftype) {
    case QLM_LATREAL:       return (void *)QDP_D_create_R_L(S->lat);
    case QLM_LATCOMPLEX:    return (void *)QDP_D_create_C_L(S->lat);
    case QLM_LATCOLVEC:     
        if     (2 == nc)    return (void *)QDP_D2_create_V_L(S->lat);
        else if(3 == nc)    return (void *)QDP_D3_create_V_L(S->lat);
        else                return (void *)QDP_DN_create_V_L(nc, S->lat);
    case QLM_LATCOLMAT:     
        if     (2 == nc)    return (void *)QDP_D2_create_M_L(S->lat);
        else if(3 == nc)    return (void *)QDP_D3_create_M_L(S->lat);
        else                return (void *)QDP_DN_create_M_L(nc, S->lat);
    case QLM_LATDIRFERM:
        if     (2 == nc)    return (void *)QDP_D2_create_D_L(S->lat);
        else if(3 == nc)    return (void *)QDP_D3_create_D_L(S->lat);
        else                return (void *)QDP_DN_create_D_L(nc, S->lat);
    case QLM_LATDIRPROP:
        if     (2 == nc)    return (void *)QDP_D2_create_P_L(S->lat);
        else if(3 == nc)    return (void *)QDP_D3_create_P_L(S->lat);
        else                return (void *)QDP_DN_create_P_L(nc, S->lat);
    default:
        qlm_error = QLM_ERR_INVALID_FIELD;
        return NULL;
    }
    return NULL; 
}
/* push lattice field to the stack top and return (void *)(QDP_Field *) pointer 
   luastack:-0+1 */
static void *
qdptype_push(lua_State *L, int Sidx, qlmLatField lftype, int nc)
{ 
    switch(lftype) {
    case QLM_LATREAL:       return (void *)(qlua_newZeroLatReal(L, Sidx)->ptr);
    case QLM_LATCOMPLEX:    return (void *)(qlua_newZeroLatComplex(L, Sidx)->ptr);
    case QLM_LATCOLVEC:     
        if     (2 == nc)    return (void *)(qlua_newZeroLatColVec2(L, Sidx, nc)->ptr);
        else if(3 == nc)    return (void *)(qlua_newZeroLatColVec3(L, Sidx, nc)->ptr);
        else                return (void *)(qlua_newZeroLatColVecN(L, Sidx, nc)->ptr);
    case QLM_LATCOLMAT:     
        if     (2 == nc)    return (void *)(qlua_newZeroLatColMat2(L, Sidx, nc)->ptr);
        else if(3 == nc)    return (void *)(qlua_newZeroLatColMat3(L, Sidx, nc)->ptr);
        else                return (void *)(qlua_newZeroLatColMatN(L, Sidx, nc)->ptr);
    case QLM_LATDIRFERM:
        if     (2 == nc)    return (void *)(qlua_newZeroLatDirFerm2(L, Sidx, nc)->ptr);
        else if(3 == nc)    return (void *)(qlua_newZeroLatDirFerm3(L, Sidx, nc)->ptr);
        else                return (void *)(qlua_newZeroLatDirFermN(L, Sidx, nc)->ptr);
    case QLM_LATDIRPROP:
        if     (2 == nc)    return (void *)(qlua_newZeroLatDirProp2(L, Sidx, nc)->ptr);
        else if(3 == nc)    return (void *)(qlua_newZeroLatDirProp3(L, Sidx, nc)->ptr);
        else                return (void *)(qlua_newZeroLatDirPropN(L, Sidx, nc)->ptr);
    default:
        qlm_error = QLM_ERR_INVALID_FIELD;
        return NULL;
    }
    return NULL; 
}
static void
qdptype_destroy(qlmLatField lftype, int nc, void *x)
{
    switch(lftype) {
    case QLM_LATREAL:       QDP_D_destroy_R((QDP_D_Real *)x);               break;
    case QLM_LATCOMPLEX:    QDP_D_destroy_C((QDP_D_Complex *)x);            break;
    case QLM_LATCOLVEC:     
        if     (2 == nc)    QDP_D2_destroy_V((QDP_D2_ColorVector *)x);
        else if(3 == nc)    QDP_D3_destroy_V((QDP_D3_ColorVector *)x);
        else                QDP_DN_destroy_V((QDP_DN_ColorVector *)x);      
        break;
    case QLM_LATCOLMAT:     
        if     (2 == nc)    QDP_D2_destroy_M((QDP_D2_ColorMatrix *)x);
        else if(3 == nc)    QDP_D3_destroy_M((QDP_D3_ColorMatrix *)x);
        else                QDP_DN_destroy_M((QDP_DN_ColorMatrix *)x);      
        break;
    case QLM_LATDIRFERM:
        if     (2 == nc)    QDP_D2_destroy_D((QDP_D2_DiracFermion *)x);
        else if(3 == nc)    QDP_D3_destroy_D((QDP_D3_DiracFermion *)x);
        else                QDP_DN_destroy_D((QDP_DN_DiracFermion *)x);     
        break;
    case QLM_LATDIRPROP:
        if     (2 == nc)    QDP_D2_destroy_P((QDP_D2_DiracPropagator *)x);
        else if(3 == nc)    QDP_D3_destroy_P((QDP_D3_DiracPropagator *)x);
        else                QDP_DN_destroy_P((QDP_DN_DiracPropagator *)x);  
        break;
    default:
        qlm_error = QLM_ERR_INVALID_FIELD;
    }
}
static void *
qdptype_expose(qlmLatField lftype, int nc, void *x)
{
    switch(lftype) {
    case QLM_LATREAL:       return (void *)QDP_D_expose_R((QDP_D_Real *)x);
    case QLM_LATCOMPLEX:    return (void *)QDP_D_expose_C((QDP_D_Complex *)x);
    case QLM_LATCOLVEC:     
        if     (2 == nc)    return (void *)QDP_D2_expose_V((QDP_D2_ColorVector *)x);
        else if(3 == nc)    return (void *)QDP_D3_expose_V((QDP_D3_ColorVector *)x);
        else                return (void *)QDP_DN_expose_V((QDP_DN_ColorVector *)x);
    case QLM_LATCOLMAT:     
        if     (2 == nc)    return (void *)QDP_D2_expose_M((QDP_D2_ColorMatrix *)x);
        else if(3 == nc)    return (void *)QDP_D3_expose_M((QDP_D3_ColorMatrix *)x);
        else                return (void *)QDP_DN_expose_M((QDP_DN_ColorMatrix *)x);
    case QLM_LATDIRFERM:
        if     (2 == nc)    return (void *)QDP_D2_expose_D((QDP_D2_DiracFermion *)x);
        else if(3 == nc)    return (void *)QDP_D3_expose_D((QDP_D3_DiracFermion *)x);
        else                return (void *)QDP_DN_expose_D((QDP_DN_DiracFermion *)x);
    case QLM_LATDIRPROP:
        if     (2 == nc)    return (void *)QDP_D2_expose_P((QDP_D2_DiracPropagator *)x);
        else if(3 == nc)    return (void *)QDP_D3_expose_P((QDP_D3_DiracPropagator *)x);
        else                return (void *)QDP_DN_expose_P((QDP_DN_DiracPropagator *)x);
    default:
        qlm_error = QLM_ERR_INVALID_FIELD;
        return NULL;
    }
    return NULL;
}
/* return value 0:OK, 1:Fail */
static int
qdptype_reset(qlmLatField lftype, int nc, void *x)
{
    switch(lftype) {
    case QLM_LATREAL:       QDP_D_reset_R((QDP_D_Real *)x);                 break;
    case QLM_LATCOMPLEX:    QDP_D_reset_C((QDP_D_Complex *)x);              break;
    case QLM_LATCOLVEC:     
        if     (2 == nc)    QDP_D2_reset_V((QDP_D2_ColorVector *)x);     
        else if(3 == nc)    QDP_D3_reset_V((QDP_D3_ColorVector *)x);
        else                QDP_DN_reset_V((QDP_DN_ColorVector *)x);        
        break;
    case QLM_LATCOLMAT:
        if     (2 == nc)    QDP_D2_reset_M((QDP_D2_ColorMatrix *)x);     
        else if(3 == nc)    QDP_D3_reset_M((QDP_D3_ColorMatrix *)x);
        else                QDP_DN_reset_M((QDP_DN_ColorMatrix *)x);        
        break;
    case QLM_LATDIRFERM:
        if     (2 == nc)    QDP_D2_reset_D((QDP_D2_DiracFermion *)x);
        else if(3 == nc)    QDP_D3_reset_D((QDP_D3_DiracFermion *)x);
        else                QDP_DN_reset_D((QDP_DN_DiracFermion *)x);       
        break;
    case QLM_LATDIRPROP:
        if     (2 == nc)    QDP_D2_reset_P((QDP_D2_DiracPropagator *)x);
        else if(3 == nc)    QDP_D3_reset_P((QDP_D3_DiracPropagator *)x);
        else                QDP_DN_reset_P((QDP_DN_DiracPropagator *)x);    
        break;
    default:
        qlm_error = QLM_ERR_INVALID_FIELD;
        return 1;
    }
    return 0;
}
#if 0 /* commented out because `qdp_vtype' is only a presentation of 
         lattice objects stored on stack; should not be created except by 
         Lua-push operation or destroyed except by garbate collector */
/* create fields 
   return value 0:OK 1:fail (set qlm_error) */
static int
qdp_vtype_create(mLattice *S, qdp_vtype *x)
{
    assert(NULL != x);
    if (x->is_array) {
        for (int i = 0 ; i < x->arr_len ; i++) {
            x->p.a[i] = qdptype_create(S, x->lftype, x->nc);
            if (NULL == x->p.a[i]) {
                for (int j = 0 ; j < i ; j++)
                    qdptype_destroy(x->lftype, x->nc, x->p.a[j]);
                qlm_error = QLM_ERR_ENOMEM; 
                return 1;
            }
        }
    }
    else {
        x->p.x = qdptype_create(S, x->lftype, x->nc);
        if (NULL == x->p.x) {
            qlm_error = QLM_ERR_ENOMEM; 
            return 1;
        }
    }
    x->is_exposed = 0;
    return 1;
}
/* destroy fields 
   return value 0:OK 1:fail (set qlm_error) */
static void
qdp_vtype_destroy(qdp_vtype *x)
{
    assert(NULL != x);
    assert(0 == x->is_exposed);
    if (x->is_array) {
        for (int i = 0 ; i < x->arr_len ; i++)
            if (NULL != x->p.a[i]) 
                qdptype_destroy(x->lftype, x->nc, x->p.a[i]);
    }
    else
        if (NULL != x->p.x)
            qdptype_destroy(x->lftype, x->nc, x->p.x);
}
#endif
/* push qdp_vtype (field or array) on top of stack;
    luastack:-0+1 */
static void
qdp_vtype_push(lua_State *L, int Sidx, qdp_vtype *x)
{
    assert(NULL != x);
    NORMALIZE_INDEX(L, Sidx);
    if (x->is_array) {
        lua_createtable(L, x->arr_len, 0);
        int tabidx = lua_gettop(L);
        for (int i = 0 ; i < x->arr_len ; i++) {
            x->p.a[i] = qdptype_push(L, Sidx, x->lftype, x->nc); /* luastack:-0+1 */
            lua_rawseti(L, tabidx, 1 + i);
        }
    } else
        x->p.x = qdptype_push(L, Sidx, x->lftype, x->nc);
    x->is_exposed = 0;
}

static void
qdp_vtype_expose(qdp_vtype *x)
{
    if (x->is_array) {
        for (int i = 0 ; i < x->arr_len ; i++)
            x->q.a[i] = qdptype_expose(x->lftype, x->nc, x->p.a[i]);
    }
    else
        x->q.x = qdptype_expose(x->lftype, x->nc, x->p.x);
    x->is_exposed = 1;
}

static void
qdp_vtype_reset(qdp_vtype *x)
{
    if (x->is_array) {
        for (int i = 0 ; i < x->arr_len ; i++)
            qdptype_reset(x->lftype, x->nc, x->p.a[i]);
    }
    else
        qdptype_reset(x->lftype, x->nc, x->p.x);
    x->is_exposed = 0;
}

static void
qdp_vtype_free(qdp_vtype *x)
{
    if (NULL == x) return;
    if (x->is_exposed) {
        qdp_vtype_reset(x);
        x->is_exposed = 0;
    }
    if (x->is_array) {
        if (NULL != x->p.a) free(x->p.a);
        x->p.a = NULL;
        if (NULL != x->q.a) free(x->q.a);
        x->q.a = NULL;
    }
    free(x);
}


typedef struct {void *ptr;} qlua_check_nocolor_s;
qlua_check_nocolor_s qlua_check_nocolor_s_NULL = {NULL};
static qlua_check_nocolor_s*
qlua_check_nocolor(lua_State *L, int idx, mLattice *S, int nc)
{
    luaL_error(L, "compiled without nc=%d support", nc);
    return &qlua_check_nocolor_s_NULL;
}
#if USE_Nc2
#define Q_NC2(x) x##2
#else 
#define Q_NC2(x) qlua_check_nocolor
#endif

#if USE_Nc3
#define Q_NC3(x) x##3
#else 
#define Q_NC3(x) qlua_check_nocolor
#endif

#if USE_NcN
#define Q_NCN(x) x##N
#else 
#define Q_NCN(x) qlua_check_nocolor
#endif

#define Qs(x,c) (2==(c) ? Q_NC2(x) : (3==(c) ? Q_NC3(x) : Q_NCN(x)))

/****** check that QDP type the object on stack is compatible with latmat */
static void *
qlm_check_qdptype(lua_State *L, int idx, const qlmData *d)
{
    switch(d->lftype) {
    case QLM_LATREAL:       return (void *)(qlua_checkLatReal(L, idx, d->S)->ptr);
    case QLM_LATCOMPLEX:    return (void *)(qlua_checkLatComplex(L, idx, d->S)->ptr);
    case QLM_LATCOLVEC:     
        if (2 == d->nc)     return (void *)(Q_NC2(qlua_checkLatColVec)(L, idx, d->S, d->nc)->ptr);
        else if(3 == d->nc) return (void *)(Q_NC3(qlua_checkLatColVec)(L, idx, d->S, d->nc)->ptr);
        else                return (void *)(Q_NCN(qlua_checkLatColVec)(L, idx, d->S, d->nc)->ptr);
    case QLM_LATCOLMAT:     
        if (2 == d->nc)     return (void *)(Q_NC2(qlua_checkLatColMat)(L, idx, d->S, d->nc)->ptr);
        else if (3 == d->nc)return (void *)(Q_NC3(qlua_checkLatColMat)(L, idx, d->S, d->nc)->ptr);
        else                return (void *)(Q_NCN(qlua_checkLatColMat)(L, idx, d->S, d->nc)->ptr);
    case QLM_LATDIRFERM:
        if (2 == d->nc)     return (void *)(Q_NC2(qlua_checkLatDirFerm)(L, idx, d->S, d->nc)->ptr);
        else if (3 == d->nc)return (void *)(Q_NC3(qlua_checkLatDirFerm)(L, idx, d->S, d->nc)->ptr);
        else                return (void *)(Q_NCN(qlua_checkLatDirFerm)(L, idx, d->S, d->nc)->ptr);
    case QLM_LATDIRPROP:
        if (2 == d->nc)     return (void *)(Q_NC2(qlua_checkLatDirProp)(L, idx, d->S, d->nc)->ptr);
        else if (3 == d->nc)return (void *)(Q_NC3(qlua_checkLatDirProp)(L, idx, d->S, d->nc)->ptr);
        else                return (void *)(Q_NCN(qlua_checkLatDirProp)(L, idx, d->S, d->nc)->ptr);
    default:
       luaL_error(L, "bad latfield type %d", d->lftype);
       return NULL;
    }
    return NULL;
}
/* check that a table contains QDP objects compatible with latmat 
   qlmData `is_array', `arr_len' are ignored - `len' is used for checks if 0 <= len
   returns an array of pointers in res_ if not NULL; must be able to accept 
   the pointers: either call with (len>0, res_=<ptr>) for existing storage
   or (len<=0, res_=NULL) to allocate storage, which must be deallocated externally
   returns a pointer to the storage
   FIXME make a struct with alloc/free to specifically handle (arrays of) QDP object ptrs? 
 */
void **
qlm_check_qdptype_array(lua_State *L, int idx, int len, void **res_, const qlmData *d)
{
    luaL_checktype(L, idx, LUA_TTABLE);
    if (0 <= len && lua_objlen(L, idx) != len) {
        luaL_error(L, "expect array[%d] of %s", len, qlm_lftype_str(d->lftype));
        return NULL;
    }
    void **res = (NULL == res_ ? malloc(len * sizeof(void *)) : res_);
    if (NULL == res) {
        luaL_error(L, "not enough memory");
        return NULL;
    }
    lua_pushvalue(L, idx);          /* in case idx is < 0 or special */
    for (int i = 0 ; i < len ; i++) {
        lua_pushinteger(L, 1 + i);  /* sic! Lua indexing */
        lua_gettable(L, -2);
        res[i] = qlm_check_qdptype(L, -1, d);
        lua_pop(L, 1);
    }
    lua_pop(L, 1);                  /* pop table ref copy */
    return res;
}
/* reset qdp lattice object */

static qdp_vtype *
qlm_check_qdp_vtype(lua_State *L, int idx, const qlmData *d)
{
    qdp_vtype *res = qdp_vtype_alloc(d->lftype, d->is_array, d->arr_len, d->nc);
    if (NULL == res) {
        luaL_error(L, "not enough memory");
        return NULL;
    }
    if (d->is_array)
        qlm_check_qdptype_array(L, idx, d->arr_len, res->p.a, d);
    else
        res->p.x = qlm_check_qdptype(L, idx, d);
    return res;
}
/* 0:OK 1:fail (set qlm_error) */
static qdp_vtype *
qlm_new_qdp_vtype(lua_State *L, int Sidx, const qlmData *d)
{
    qdp_vtype *res = qdp_vtype_alloc(d->lftype, d->is_array, d->arr_len, d->nc);
    if (NULL == res) {
        luaL_error(L, "(qlm_new_qdp_vtype) not enough memory\n");
        return NULL;
    }
    qdp_vtype_push(L, Sidx, res);
    return res;
}


/****************************************************************************/
/* site layout functions: conversion between coord and linear indices 
   GLOSSARY
    x2i     coordinate to index
    i2x     index to coordinate
    eo      even/odd (can take parity information from qlmData)

    dim     array of dimensions
    coord   array of coordinates, [0; dim)
    len     length of storage indexed by `idx'
    idx     linear index, [0; len)
    XY
 */
/* aux functions: coord <-> lexicographic index */
static void
index_i2x(int ndim, const int *restrict dim, int *restrict x, int ilex)
{
    for (int i = 0 ; i < ndim ; i++) {
        x[i]    = ilex % dim[i];
        ilex   /= dim[i];
    }
}
static void 
index_x2i(int ndim, const int *restrict dim, int *restrict ilex, const int *restrict x)
{
    int res = 0;
    for (int i = ndim ; i-- ;)
        res = res * dim[i] + x[i];
    *ilex = res;
}
static void
/* aux functions: coord <-> lexicographic index, even/odd */
index_i2x_eo(int ndim, const int *restrict dim, int *restrict x, int ilex_eo, int par)
{
    int ilex = 2 * ilex_eo,
        p = 0;
    for (int i = 0 ; i < ndim ; i++) {
        x[i]    = ilex % dim[i];
        p      += x[i];
        ilex   /= dim[i];
    }
    if ((p & 1) ^ (par & 1)) {  /* wrong parity */
        ilex = 2 * ilex_eo + 1;
        for (int i = 0 ; i < ndim ; i++) {
            x[i]    = ilex % dim[i];
            ilex   /= dim[i];
        }
    }
}
static void 
index_x2i_eo(int ndim, const int *restrict dim, int *restrict ilex, int *restrict par, const int *restrict x)
{
    index_x2i(ndim, dim, ilex, x);
    *par    = *ilex % 2;
    *ilex  /= 2;
}


/* prototypes for coord <-> linear index conversion */
typedef void (*site_layout_i2x)(int *restrict x, const int site_idx, void *arg_);
typedef void (*site_layout_x2i)(int *restrict site_idx, const int *restrict x, void *arg_);

/* qdp lattice site layout */
static void 
site_layout_x2i_qdp(int *restrict vs_idx, const int *restrict vs_coord, void *arg_)
{
    qlmData *d = (qlmData *)arg_;
    int ls_coord[QLUA_MAX_LATTICE_RANK];
    for (int i = 0 ; i < d->ndim; i++)
        ls_coord[i] = d->x0[i] + vs_coord[i];
    /* TODO check coordinate range? */
    if (QDP_node_number_L(d->S->lat, ls_coord) != QDP_this_node)
        *vs_idx = -1;
    else
        *vs_idx = QDP_index_L(d->S->lat, ls_coord);
}
static void 
site_layout_i2x_qdp(int *restrict vs_coord, const int vs_idx, void *arg_)
{
    qlmData *d = (qlmData *)arg_;
    int ls_coord[QLUA_MAX_LATTICE_RANK];
    QDP_get_coords_L(d->S->lat, ls_coord, QDP_this_node, vs_idx);
    for (int i = 0 ; i < d->ndim; i++)
        vs_coord[i] = ls_coord[i] - d->x0[i];
}

/* lexicographic site layout - default for FULL sublat */
static void 
site_layout_x2i_lex(int *restrict vs_idx, const int *restrict vs_coord, void *arg_)
{
    qlmData *d = (qlmData *)arg_;
    assert(QLM_SUBLAT_FULL == d->sublat);
    index_x2i(d->ndim, d->vec_site_dim, vs_idx, vs_coord);
}
static void 
site_layout_i2x_lex(int *restrict vs_coord, const int vs_idx, void *arg_)
{
    qlmData *d = (qlmData *)arg_;
    assert(QLM_SUBLAT_FULL == d->sublat);
    index_i2x(d->ndim, d->vec_site_dim, vs_coord, vs_idx);
}


/* blocked site layout */
static void
site_layout_x2i_blk(int *restrict vs_idx, const int *restrict vs_coord, void *arg_)
{
    qlmData *d = (qlmData *)arg_;
    assert(QLM_SUBLAT_FULL == d->sublat);
    int vb_coord[QLUA_MAX_LATTICE_RANK],    /* block in node */
        bs_coord[QLUA_MAX_LATTICE_RANK];    /* site in block */
    for (int i = 0 ; i < d->ndim ; i++) {
        bs_coord[i] = vs_coord[i] % d->blk_site_dim[i];
        vb_coord[i] = vs_coord[i] / d->blk_site_dim[i];
    }
    int vb_idx, bs_idx;
    index_x2i(d->ndim, d->vec_blk_dim, &vb_idx, vb_coord);
    index_x2i(d->ndim, d->blk_site_dim, &bs_idx, bs_coord);
    *vs_idx = bs_idx + d->blk_site_len * vb_idx;
}
static void
site_layout_i2x_blk(int *restrict vs_coord, const int vs_idx, void *arg_)
{
    qlmData *d = (qlmData *)arg_;
    assert(QLM_SUBLAT_FULL == d->sublat);
    int vb_coord[QLUA_MAX_LATTICE_RANK],    /* block in node */
        bs_coord[QLUA_MAX_LATTICE_RANK];    /* site in block */
    int vb_idx, bs_idx;
    bs_idx  = vs_idx % d->blk_site_len;
    vb_idx  = vs_idx / d->blk_site_len;
    index_i2x(d->ndim, d->vec_blk_dim, vb_coord, vb_idx);
    index_i2x(d->ndim, d->blk_site_dim, bs_coord, bs_idx);
    for (int i = 0 ; i < d->ndim ; i++)
        vs_coord[i] = bs_coord[i] + d->blk_site_dim[i] * vb_coord[i];
}


/* lexicographic site layout, even/odd */
static void 
site_layout_x2i_lex_eo(int *restrict vs_idx, const int *restrict vs_coord, void *arg_)
{
    qlmData *d = (qlmData *)arg_;
    assert(IS_SUBLAT_EOPC(d->sublat));
    int vs_par;
    index_x2i_eo(d->ndim, d->vec_site_dim, vs_idx, &vs_par, vs_coord);
    if (SUBLAT_PARITY(d->sublat) != (vs_par ^ d->x0_parity))
        *vs_idx = 0;
}
static void 
site_layout_i2x_lex_eo(int *restrict vs_coord, const int vs_idx, void *arg_)
{
    qlmData *d = (qlmData *)arg_;
    assert(IS_SUBLAT_EOPC(d->sublat));
    index_i2x_eo(d->ndim, d->vec_site_dim, vs_coord, vs_idx, 
            SUBLAT_PARITY(d->sublat) ^ d->x0_parity);
}


/* blocked site layout, even/odd */
static void
site_layout_x2i_blk_eo(int *restrict vs_idx, const int *restrict vs_coord, void *arg_)
{
    qlmData *d = (qlmData *)arg_;
    assert(IS_SUBLAT_EOPC(d->sublat));
    int vb_coord[QLUA_MAX_LATTICE_RANK],    /* block in node */
        bs_coord[QLUA_MAX_LATTICE_RANK];    /* site in block */
    int b0_par = d->x0_parity;              /* parity of LO blk corner */
    for (int i = 0 ; i < d->ndim ; i++) {
        bs_coord[i] = vs_coord[i] % d->blk_site_dim[i];
        vb_coord[i] = vs_coord[i] / d->blk_site_dim[i];
        b0_par     ^= (bs_coord[i] * d->blk_site_dim[i]) & 1;
    }
    int vb_idx, bs_idx, bs_par;
    index_x2i(d->ndim, d->vec_blk_dim, &vb_idx, vb_coord);
    index_x2i_eo(d->ndim, d->blk_site_dim, &bs_idx, &bs_par, bs_coord);
    if (SUBLAT_PARITY(d->sublat) != (b0_par ^ bs_par))
        *vs_idx = -1;   /* not in sublat */
    else
        *vs_idx = bs_idx + d->blk_site_len * vb_idx;
}
static void
site_layout_i2x_blk_eo(int *restrict vs_coord, const int vs_idx, void *arg_)
{
    qlmData *d = (qlmData *)arg_;
    assert(IS_SUBLAT_EOPC(d->sublat));
    int vb_coord[QLUA_MAX_LATTICE_RANK],    /* block in node */
        bs_coord[QLUA_MAX_LATTICE_RANK];    /* site in block */
    int vb_idx, bs_idx;
    int bs_par = d->x0_parity ^ SUBLAT_PARITY(d->sublat);
    bs_idx  = vs_idx % d->blk_site_len;
    vb_idx  = vs_idx / d->blk_site_len;
    index_i2x(d->ndim, d->vec_blk_dim, vb_coord, vb_idx);
    for (int i = 0 ; i < d->ndim ; i++)
        bs_par ^= (vb_coord[i] * d->vec_blk_dim[i]) & 1;
    index_i2x_eo(d->ndim, d->blk_site_dim, bs_coord, bs_idx, bs_par);
    for (int i = 0 ; i < d->ndim ; i++)
        vs_coord[i] = bs_coord[i] + d->blk_site_dim[i] * vb_coord[i];
}


/****** import / export */
/* <QDPt>[<color>]<int_prec><qdp_prec>[<color>] */
typedef void (*field_copy)(void *restrict dst_, int dst_idx, 
                           const void *restrict src_, int src_idx, int iarr, void *arg_);
#define Qpp(func) func ## fd
#define pint    float
#define QT(x)   QLA_D##x
#include <qlm_layout-x.c>

#define Qpp(func) func ## dd
#define pint    double
#define QT(x)   QLA_D##x
#include <qlm_layout-x.c>


#define QLM_P(a)    (QLM_PREC_FLOAT==d->prec ? a##fd : (QLM_PREC_DOUBLE==d->prec ? a##dd : NULL))
#define QLM_PC(a)   (2==d->nc ? QLM_P(a##2) : (3==d->nc ? QLM_P(a##3) : NULL ))

/* RETURN 0 : ok, 1 : fail (set qlm_error) */
int
qlm_import_qdp_vtype(qlmData *d, int ivec, const qdp_vtype *qx)
{
    site_layout_x2i src_x2i = site_layout_x2i_qdp;
    site_layout_i2x dst_i2x;
    if (QLM_SUBLAT_FULL == d->sublat)
        dst_i2x = d->is_blocked ? site_layout_i2x_blk : site_layout_i2x_lex ;
    else if (IS_SUBLAT_EOPC(d->sublat))
        dst_i2x = d->is_blocked ? site_layout_i2x_blk_eo : site_layout_i2x_lex_eo ;
    else {
        qlm_error = QLM_ERR_ENUM_ERROR;
        return 1;
    }
    field_copy fc = NULL;
    switch(d->lftype) {
    case QLM_LATREAL:       fc = QLM_P (fc_qla2int_R); break;
    case QLM_LATCOMPLEX:    fc = QLM_P (fc_qla2int_C); break;
    case QLM_LATCOLVEC:     fc = QLM_PC(fc_qla2int_V); break;
    case QLM_LATCOLMAT:     fc = QLM_PC(fc_qla2int_M); break;
    case QLM_LATDIRFERM:    fc = QLM_PC(fc_qla2int_D); break;
    case QLM_LATDIRPROP:    fc = QLM_PC(fc_qla2int_P); break;
    default:
        qlm_error = QLM_ERR_ENUM_ERROR;
       return 1;
    }
    if (NULL == fc) {
        qlm_error = QLM_ERR_ENUM_ERROR;
        return 1;
    }

    int vs_coord[QLUA_MAX_LATTICE_RANK];
    int vs_idx_src;
    void *ivec_data = d->vec_data + ivec * d->vec_size;
    for (int vs_idx_dst = 0 ; vs_idx_dst < d->vec_site_len ; vs_idx_dst++) {
        dst_i2x(vs_coord, vs_idx_dst, (void *)d);
        src_x2i(&vs_idx_src, vs_coord, (void *)d);
        if (vs_idx_src < 0) {
            qlm_error = QLM_ERR_GEOM_SITE;
            return 1;
        }
        if (d->is_array) {
            for (int a = 0 ; a < d->arr_len ; a++)
                fc(ivec_data, vs_idx_dst, qx->q.a, vs_idx_src, a, d);
        } else
            fc(ivec_data, vs_idx_dst, &(qx->q.x), vs_idx_src, 0, d);
    }
    return 0;
}

/* return value 0:OK, 1:Fail (set qlm_error) */
int
qlm_export_qdp_vtype(qdp_vtype *qx, qlmData *d, int ivec)
{
    site_layout_x2i dst_x2i = site_layout_x2i_qdp;
    site_layout_i2x src_i2x = NULL;
    if (QLM_SUBLAT_FULL == d->sublat)
        src_i2x = d->is_blocked ? site_layout_i2x_blk : site_layout_i2x_lex ;
    else if (IS_SUBLAT_EOPC(d->sublat))
        src_i2x = d->is_blocked ? site_layout_i2x_blk_eo : site_layout_i2x_lex_eo ;
    else {
        qlm_error = QLM_ERR_ENUM_ERROR;
        return 1;
    }
    field_copy fc = NULL;
    switch(d->lftype) {
    case QLM_LATREAL:       fc = QLM_P (fc_int2qla_R); break;
    case QLM_LATCOMPLEX:    fc = QLM_P (fc_int2qla_C); break;
    case QLM_LATCOLVEC:     fc = QLM_PC(fc_int2qla_V); break;
    case QLM_LATCOLMAT:     fc = QLM_PC(fc_int2qla_M); break;
    case QLM_LATDIRFERM:    fc = QLM_PC(fc_int2qla_D); break;
    case QLM_LATDIRPROP:    fc = QLM_PC(fc_int2qla_P); break;
    default:
        qlm_error = QLM_ERR_ENUM_ERROR;
        return 1;
    }
    if (NULL == fc) {
        qlm_error = QLM_ERR_ENUM_ERROR;
        return 1;
    }

    int vs_coord[QLUA_MAX_LATTICE_RANK];
    int vs_idx_dst;
    void *ivec_data = d->vec_data + ivec * d->vec_size;
    for (int vs_idx_src = 0 ; vs_idx_src < d->vec_site_len ; vs_idx_src++) {
        src_i2x(vs_coord, vs_idx_src, (void *)d);
        dst_x2i(&vs_idx_dst, vs_coord, (void *)d); 
        if (vs_idx_dst < 0) {
            qlm_error = QLM_ERR_GEOM_SITE;
            return 1;
        }
        if (d->is_array) {
            for (int a = 0 ; a < d->arr_len ; a++)
                fc(qx->q.a, vs_idx_dst, ivec_data, vs_idx_src, a, d);
        } else
            fc(&(qx->q.x), vs_idx_dst, ivec_data, vs_idx_src, 0, d);
    }
    return 0;
}
#undef QLM_P
#undef QLM_PC



/***********************************************************************
 * interface
 ***********************************************************************/
typedef struct { /* XXX why do I wrap it? */
    qlmData     *d;
} mQlm;
static const char qlm_str[] = "latmat";
static const char qlm_name[] = "latmat";
static const char qlm_meta_name[] = "qcd.latmat";
static mQlm *qlua_checkQlm(lua_State *L, int idx, mLattice *S);
mQlm *qlua_newQlm(lua_State *L, int Sidx, qlmData *d);

/******** latmat methods ********/    
static int
q_qlm_fmt(lua_State *L)
{
    char fmt[256];
    char fmt_array[16]  = "";
    char fmt_sublat[16] = "";
    char fmt_block[64]  = "";
    char fmt_layout[64] = "";
    char fmt_field[64]  = "";
    const char *lftype_str = "";
    mQlm *c = qlua_checkQlm(L, 1, NULL);
    qlmData *d = c->d;
    lftype_str = qlm_lftype_str(d->lftype);
    if (d->is_array)
        snprintf(fmt_array, sizeof(fmt_array),
                " array[%d]", d->arr_len);
    if (IS_SUBLAT_EOPC(d->sublat)) 
        snprintf(fmt_sublat, sizeof(fmt_sublat), ",eo(%d)", SUBLAT_PARITY(d->sublat));
    if (d->is_blocked) {
        snprintf_arr(fmt_block, sizeof(fmt_block), "%d", ",", 
                int, d->blk_site_dim, d->S->rank);
        snprintf(fmt_layout, sizeof(fmt_layout), ",BLOCK(%s,nbasis=%d)", 
                fmt_block, d->nvec_basis);
    }
    snprintf(fmt_field, sizeof(fmt_field), " nc=%d ns=%d", d->nc, d->ns);

    snprintf(fmt, sizeof(fmt), 
            "%s(%dx %s%s%s%s%s)", qlm_str, d->nvec, lftype_str, 
            fmt_array, fmt_field, fmt_sublat, fmt_layout);
    lua_pushstring(L, fmt);

    return 1;
}
static int
q_qlm_gc(lua_State *L)
{
    mQlm *c = qlua_checkQlm(L, 1, NULL);

    if (NULL != c->d) {
        qlm_data_free(c->d);
        c->d = NULL;
    }

    return 0;
}
static int 
q_qlm_index(lua_State *L)
{
    mQlm *c; 
    c = qlua_checkQlm(L, 1, NULL);
    const char *what_str = luaL_checkstring(L, 2);

    if (!strcmp("lattice", what_str)) {
        lua_getfenv(L, 1);
        lua_getfield(L, -1, "lattice");
        return 1;
    } 
    /* add other read-only fields */
    else {
        return luaL_error(L, "invalid field: %s", what_str);
    }
}
static int 
q_qlm_is_array(lua_State *L)
{
    mQlm *c = qlua_checkQlm(L, 1, NULL);
    lua_pushboolean(L, c->d->is_array);
    return 1;
}
static int q_qlm_array_len(lua_State *L)
{
    mQlm *c = qlua_checkQlm(L, 1, NULL);
    lua_pushinteger(L, c->d->arr_len);
    return 1;
}
static int q_qlm_get_lattice(lua_State *L) 
{
    mQlm *c = NULL;
    c = qlua_checkQlm(L, 1, NULL);
    lua_getfenv(L, 1);
    lua_getfield(L, -1, "lattice");
    return 1;
}
static int q_qlm_get_nvec(lua_State *L)
{
    mQlm *c = qlua_checkQlm(L, 1, NULL);
    lua_pushinteger(L, c->d->nvec);
    return 1;
}
/* [lua] latmat:insert(ivec, VEC) */
static int 
q_qlm_insert(lua_State *L)
{
    mQlm *c = qlua_checkQlm(L, 1, NULL);
    qlmData *d = c->d;

    int ivec = luaL_checkint(L, 2);
    if (ivec < 0 || c->d->nvec <= ivec) 
        luaL_error(L, "index %d is outside range [0;%d)", ivec, c->d->nvec);

    qdp_vtype *vec = qlm_check_qdp_vtype(L, 3, d);
    if (NULL == vec)
        return luaL_error(L, "expect vtype at argument 3; %s", qlm_strerror(qlm_error));
    qdp_vtype_expose(vec);
    if (qlm_import_qdp_vtype(d, ivec, vec)) {
        qdp_vtype_reset(vec);
        qdp_vtype_free(vec);
        return luaL_error(L, "import error: %s", qlm_strerror(qlm_error));
    }

    qdp_vtype_reset(vec);
    qdp_vtype_free(vec);
    return 0;   /* no return */
}
/* [lua] latmat:extract(ivec) ; return a (table of) latfields */
static int 
q_qlm_extract(lua_State *L)
{
    mQlm *c     = qlua_checkQlm(L, 1, NULL);
    qlmData *d  = c->d;
    mLattice *S = NULL;
    S = qlua_ObjLattice(L, 1);
    int Sidx    = lua_gettop(L);

    int ivec = luaL_checkint(L, 2);
    if (ivec < 0 || c->d->nvec <= ivec) 
        luaL_error(L, "index %d is outside range [0;%d)", ivec, c->d->nvec);

    qdp_vtype *vec = qlm_new_qdp_vtype(L, Sidx, d);
    if (NULL == vec)
        return luaL_error(L, "failed to create lattice object (qlm_new_qdp_vtype); %s", 
                qlm_strerror(qlm_error));
    qdp_vtype_expose(vec);

    if (qlm_export_qdp_vtype(vec, d, ivec)) {
        qdp_vtype_reset(vec);
        qdp_vtype_free(vec);
        return luaL_error(L, "export error: %s", qlm_strerror(qlm_error));
    }

    qdp_vtype_reset(vec);
    qdp_vtype_free(vec);
    return 1;
}


static struct luaL_Reg mtQlm[] = {
    { "__tostring",                 q_qlm_fmt                       },
    { "__gc",                       q_qlm_gc                        },
//    { "__index",                    q_qlm_index     },
    { "extract",                    q_qlm_extract                   },
    { "insert",                     q_qlm_insert                  },
//    { "copy",                       q_qlm_copy                      },
//    { "read_compressed",            q_qlm_read_compressed           },
    { NULL,                         NULL                            }
};

static mQlm *
qlua_checkQlm(lua_State *L, int idx, mLattice *S) {
    mQlm *c = qlua_checkLatticeType(L, idx, qLatmat, qlm_meta_name);
    if (NULL != S) {
        mLattice *S1 = qlua_ObjLattice(L, idx);
        if (S1->id != S->id)
            luaL_error(L, "%s on a wrong lattice", qlm_meta_name);
        lua_pop(L, 1);
    }
        
    if (NULL == c->d)
        luaL_error(L, "NULL qlm object");
    
    return c;
}

/* lattice object is stored in the userdata's env */
mQlm *
qlua_newQlm(lua_State *L, int Sidx, qlmData *d)
{
    mQlm *c = lua_newuserdata(L, sizeof(mQlm));
    c->d    = d;
    qlua_createLatticeTable(L, Sidx, mtQlm, qLatmat, qlm_meta_name);
    if ( !lua_setmetatable(L, -2) )
        luaL_error(L, "cannot set meta");

    /* initialize environment */
    lua_newtable(L);
    lua_pushvalue(L, Sidx);
    lua_setfield(L, -2, "lattice");
    /* more values to store in env */
    if (! lua_setfenv(L, -2))
        luaL_error(L, "cannot set env");

    /* TODO */
    return c;
}
/* qcd.latmat.create(lat, type, nvec, 
        {array=<arr_len>, ns=<ns>, nc=<nc>, 
        prec[ision]="float"|"double", 
        sublat="even"|"odd"|"full"},
        block={bx,by,bz,bt},
        nvec_basis=<nvec_basis>})
 */
#define QLM_NC_DEFAULT  3
#define QLM_NS_DEFAULT  4
static int
q_qlm_create(lua_State *L)
{
    if (lua_gettop(L) < 3)
        return qlua_badconstr(L, qlm_str);
    mLattice *S = qlua_checkLattice(L, 1);
    int nvec = 0; 
    qlmLatField lftype;
    const char *lftype_str = luaL_checkstring(L, 2);
    if (QLM_LAT_NONE == (lftype = qlm_str2lftype(lftype_str)))
        luaL_error(L, "bad lftype='%s'", lftype_str);
   
    if ((nvec = luaL_checkinteger(L, 3))<= 0)
        luaL_error(L, "bad nvec=%d", nvec);

    int is_array = 0, 
        arr_len = 1;
    int is_blocked = 0,
        nvec_basis = nvec;
    int *bs_dim = NULL;
    int ns = QLM_NS_DEFAULT, 
        nc = QLM_NC_DEFAULT;
    qlmPrec prec = QLM_PREC_FLOAT;
    qlmSublat sublat = QLM_SUBLAT_FULL;

    const int oi = 4;
    if (qlua_checkopt_paramtable(L, oi)) {
        const char *opt_str;
        opt_str = qlua_tabkey_stringopt(L, oi, "sublat", "full");
        if (QLM_SUBLAT_NONE == (sublat = qlm_str2sublat(opt_str)))
            luaL_error(L, "bad sublat='%s'", opt_str);

        opt_str = qlua_tabkey_stringopt(L, oi, "prec", "double");
        if (QLM_PREC_NONE == (prec = qlm_str2prec(opt_str)))
            luaL_error(L, "bad prec='%s'", opt_str);

        is_array = qlua_tabkey_present(L, oi, "array");
        if (is_array) {
            if ((arr_len = qlua_tabkey_intopt(L, oi, "array", -1)) <= 0)
                luaL_error(L, "bad arr_len=%d", arr_len);
        } else arr_len = 1;

        is_blocked = qlua_tabkey_present(L, oi, "block");
        if (is_blocked) {
            qlua_tabkey(L, oi, "block");
            bs_dim = qlua_checkintarray(L, -1, S->rank, NULL);
            lua_pop(L, 1);
            if (NULL == bs_dim) 
                luaL_error(L, "bad blocksize");
            for (int i = 0 ; i < S->rank ; i++) 
                if (bs_dim[i] <= 0)
                    luaL_error(L, "bad blocksize[%d]=%d\n", i, bs_dim[i]);

            nvec_basis = qlua_tabkey_int(L, oi, "nvec_basis");
            if (nvec_basis <= 0 || nvec < nvec_basis)
                luaL_error(L, "bad nvec_basis=%d", nvec_basis);
        }

        if ((ns = qlua_tabkey_intopt(L, oi, "ns", QLM_NS_DEFAULT)) <= 0)
            luaL_error(L, "bad Ns=%d", ns);
        if ((nc = qlua_tabkey_intopt(L, oi, "nc", QLM_NC_DEFAULT)) <= 0)
            luaL_error(L, "bad Nc=%d", nc);
    }

    qlmData *d  = qlm_data_alloc(S, sublat, nvec, lftype, 
            is_array, arr_len, ns, nc, prec, 
            is_blocked, bs_dim, nvec_basis);
    if (NULL == d) 
        luaL_error(L, "qlm_data_alloc: return NULL: %s", qlm_strerror(qlm_error));

    mQlm *c = NULL;
    c  = qlua_newQlm(L, 1, d);

    if (NULL != bs_dim) 
        qlua_free(L, bs_dim);
    return 1;
}

void node2coord(int *x, int n, mLattice *S);
int coord2node(const int *x, mLattice *S);

/* qlm_sandbox(what_str, ...) */
int q_qlm_sandbox(lua_State *L)
{
    const char *what_str = luaL_checkstring(L, 1);
    if (!strcmp("replace_func_env", what_str)) {
        if (lua_gettop(L) < 3)
            return luaL_error(L, "%s: expect 3 args", what_str);
        lua_pushvalue(L, 2);
        lua_gettable(L, LUA_ENVIRONINDEX);
        lua_pushvalue(L, 2);
        lua_pushvalue(L, 3);
        lua_settable(L, LUA_ENVIRONINDEX);
        return 1;
    } else if (!strcmp("print_geom", what_str)) {
        int lo[QLUA_MAX_LATTICE_RANK],
            hi[QLUA_MAX_LATTICE_RANK],
            nD[QLUA_MAX_LATTICE_RANK];

        mLattice *S = qlua_checkLattice(L, 2);
        node2coord(nD, QDP_this_node, S);
        qlua_sublattice(lo, hi, QDP_this_node, S);
        
        char buf[1024], fmt_nD[1024], fmt_lo[128], fmt_hi[128];
        snprintf_arr(fmt_nD, sizeof(fmt_nD), "%2d", ",", int, nD, S->rank);
        snprintf_arr(fmt_lo, sizeof(fmt_lo), "%2d", ",", int, lo, S->rank);
        snprintf_arr(fmt_hi, sizeof(fmt_hi), "%2d", ",", int, hi, S->rank);
        snprintf(buf, sizeof(buf), "[%4d] {%s} : {%s}..{%s}", 
                 QDP_this_node, fmt_lo, fmt_hi, fmt_nD);
        printf("%s\n", buf);
        return 0;
    } else {
        return luaL_error(L, "unknown what_str='%s'", what_str);
    }
    /* should not get here */
    return 0;
}



static int
qlm_selftest()
{
//    return 1;   /* not implemented */
    return 0;
}

static struct luaL_Reg fQlm[] = {
    { "create",             q_qlm_create     },
    { "qlm_sandbox",        q_qlm_sandbox    },
    { NULL,            NULL             }
};

int init_qlm(lua_State *L)
{
    int status;
    if (0 != (status = qlm_selftest()))
        return status;     /* internal checks for layouts, etc */

    lua_getglobal(L, qcdlib);
    lua_newtable(L);
    luaL_register(L, NULL, fQlm);
    lua_setfield(L, -2, qlm_name);
    lua_pop(L, 1);
    return 0;
}

void fini_qlm(void)
{
}
