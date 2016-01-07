#include "qlopt.h"                                                   /* DEPS */
#include "qcomplex.h"                                                /* DEPS */
#include "modules.h"                                                 /* DEPS */
#include "lattice.h"                                                 /* DEPS */
#include "qlayout.h"                                                 /* DEPS */
#include "latint.h"                                                  /* DEPS */
#include "latcomplex.h"                                              /* DEPS */
#include "latreal.h"                                                 /* DEPS */
#include "latcolvec.h"                                               /* DEPS */
#include "latcolmat.h"                                               /* DEPS */
#include "latdirferm.h"                                              /* DEPS */
#include "latdirprop.h"                                              /* DEPS */
#include "qgamma.h"                                                  /* DEPS */

#include <assert.h>
#include <string.h>

/* aux: calculate a product of integer array */
static int 
prod_int_(int n, int *a) 
{
    int res;
    for (res = 1; n--; res *= *(a++));
    return res;
}



/* qlopt_read_type_internal - read generic element from stack
   base == NULL : use t->p as the destination
   base != NULL : use base + t->off as the destination
 */
static const char * 
qlopt_type_name(qloptType argtype)
{
    switch(argtype) {
    case QLOPT_NONE         : return "nil";
    case QLOPT_INT          : return "int";
    case QLOPT_BOOL         : return "bool";
    case QLOPT_REAL         : return "real";
    case QLOPT_COMPLEX      : return "complex";
    case QLOPT_STRING       : return "string";
    case QLOPT_LATINT       : return "LatInt";
    case QLOPT_LATREAL      : return "LatReal";
    case QLOPT_LATCOMPLEX   : return "LatComplex";
    case QLOPT_LATCOLVEC    : return "LatColVec";
    case QLOPT_LATCOLMAT    : return "LatColMat";
    case QLOPT_LATDIRFERM   : return "LatDirFerm";
    case QLOPT_LATDIRPROP   : return "LatDirProp";
    case QLOPT_GAMMA        : return "Gamma";
    case QLOPT_SUBSET       : return "Subset";
    case QLOPT_TABLE        : return "table";
    case QLOPT_ARRAY        : return "array";
    default : return "(unknown)";
    }
}

static size_t
qlopt_type_size(qloptType argtype)
{
    switch(argtype) {
    case QLOPT_INT          : 
    case QLOPT_BOOL         : 
        return sizeof(qlopt_Int);

    case QLOPT_REAL         : 
        return sizeof(qlopt_Real);

    case QLOPT_COMPLEX      : 
        return sizeof(qlopt_Complex);

    case QLOPT_STRING       : 
    case QLOPT_LATINT       : 
    case QLOPT_LATREAL      : 
    case QLOPT_LATCOMPLEX   : 
    case QLOPT_LATCOLVEC    : 
    case QLOPT_LATCOLMAT    : 
    case QLOPT_LATDIRFERM   : 
    case QLOPT_LATDIRPROP   : 
    case QLOPT_GAMMA        : 
    case QLOPT_SUBSET       : 
        return sizeof(void *);  /* sic! all pointers of the same size */

    case QLOPT_NONE         : 
    case QLOPT_TABLE        : 
    case QLOPT_ARRAY        : 
    default : 
        return SIZE_MAX ;       /* size is meaningless, must be handled separately */
    }
}

/* invoke on type mismatch to print warning or terminate, depending on flags */
static int
qlopt_arg_mismatch(lua_State *L, int pos, int flags, 
        qloptType argtype, const char *msg_head)
{
    if (flags & QLOPT_BAD_SKIP) return 0;
    if (flags & QLOPT_BAD_WARN) {
        if (QDP_this_node == qlua_master_node)
            printf("%s mismatch: expect %s", msg_head, 
                   qlopt_type_name(argtype));
        return 0;
    }
    return luaL_error(L, "%s mismatch: expect %s", msg_head, 
                      qlopt_type_name(argtype));
}
/* invoke on missing arg to print warning or terminate, depending on flags */
static int
qlopt_arg_missing(lua_State *L, int pos, int flags, 
        qloptType argtype, const char *msg_head)
{
    if (flags & QLOPT_OPT) return 0;
    if (flags & QLOPT_OPT_WARN) {
        if (QDP_this_node == qlua_master_node) 
            printf("%s missing: expect %s", msg_head, 
                   qlopt_type_name(argtype));
        return 0;
    }
    return luaL_error(L, "%s missing: expect %s", msg_head, 
                      qlopt_type_name(argtype));
}

static int 
qlopt_read_arg_internal(lua_State *L, int pos, int flags, 
        qloptType argtype, qloptData t, void *base,
        const char *msg_head);

static int
qlopt_read_table_internal(lua_State *L, int pos, 
            qloptElem *elem_list, void *base, const char *msg_head);

static int 
qlopt_read_array_checkdim_alloc(lua_State *L, int pos, 
        qloptNdarray *a, const char *msg_head);
static int 
qlopt_read_array_internal(lua_State *L, int i, int pos_i, qloptNdarray *a, 
        void *buf_i, const char *msg_head);


/* parse an argument from a Lua stack
   XXX strict typing is imposed in some cases, e.g. parsing will fail 
   when string <-> number conversion is necessary; this conversion is 
   typically performed by lua_tonumber and lua_tostring that have side effects 
 */
static int 
qlopt_read_arg_internal(lua_State *L, int pos, int flags, 
        qloptType argtype, qloptData t, void *base,
        const char *msg_head)
{
    int status = 0;
    char strbuf[1024];
    if (NULL == msg_head) {
        snprintf(strbuf, sizeof(strbuf), "(arg %d)", pos);
        msg_head = strbuf;
    }
    if (lua_gettop(L) < pos 
        || ( (LUA_TNONE == lua_type(L, pos) || LUA_TNIL == lua_type(L, pos))
             && QLOPT_NONE != argtype)) 
        return qlopt_arg_missing(L, pos, flags, argtype, msg_head);

    void *p = (NULL != base ? base + t.off : t.ptr);
    switch(argtype) {
    case QLOPT_NONE         : return 0;
    case QLOPT_BOOL         : {
        if (lua_type(L, pos) != LUA_TBOOLEAN)
            return qlopt_arg_mismatch(L, pos, flags, argtype, msg_head);
        *(qlopt_Int *)p = lua_toboolean(L, pos);
        return 0; 
    }
    case QLOPT_INT          : {
        if (lua_type(L, pos) != LUA_TNUMBER) /* using lua_isnumber may be slower */
            return qlopt_arg_mismatch(L, pos, flags, argtype, msg_head);
        *(qlopt_Int *)p = luaL_checkinteger(L, pos);
        return 0; 
    }
    case QLOPT_REAL         : {
        if (lua_type(L, pos) != LUA_TNUMBER) /* using lua_isnumber may be slower */
            return qlopt_arg_mismatch(L, pos, flags, argtype, msg_head);
        *(qlopt_Real *)p    = luaL_checknumber(L, pos);
        return 0; 
    }
    case QLOPT_COMPLEX      : {
        if (qlua_qtype(L, pos) != qComplex)
            return qlopt_arg_mismatch(L, pos, flags, argtype, msg_head);
        QLA_D_Complex *c = qlua_checkComplex(L, pos);
        *(double complex *)p = QLA_real(*c) + I*QLA_imag(*c);
        return 0; 
    }
    case QLOPT_STRING       : {
        if (lua_type(L, pos) != LUA_TSTRING)
            return qlopt_arg_mismatch(L, pos, flags, argtype, msg_head);
        *(const char **)p = luaL_checkstring(L, pos);
        return 0; 
    }
    case QLOPT_LATINT       : {
        if (qlua_qtype(L, pos) != qLatInt)
            return qlopt_arg_mismatch(L, pos, flags, argtype, msg_head);
        *(QDP_Int **)p = qlua_checkLatInt(L, pos, NULL)->ptr;
        return 0; 
    }
    case QLOPT_LATREAL      : {
        if (qlua_qtype(L, pos) != qLatReal)
            return qlopt_arg_mismatch(L, pos, flags, argtype, msg_head);
        *(QDP_Real **)p = qlua_checkLatReal(L, pos, NULL)->ptr;
        return 0; 
    }
    case QLOPT_LATCOMPLEX   : { 
        if (qlua_qtype(L, pos) != qLatComplex)
            return qlopt_arg_mismatch(L, pos, flags, argtype, msg_head);
        *(QDP_Complex **)p = qlua_checkLatComplex(L, pos, NULL)->ptr;
        return 0; 
    }
    /* TODO iterate over Nc */
    case QLOPT_LATCOLVEC    : {
        if (qlua_qtype(L, pos) != qLatColVec3)
            return qlopt_arg_mismatch(L, pos, flags, argtype, msg_head);
        *(QDP_D3_ColorVector **)p = qlua_checkLatColVec3(L, pos, NULL, 3)->ptr;
        return 0; 
    }
    case QLOPT_LATCOLMAT    : {
        if (qlua_qtype(L, pos) != qLatColMat3)
            return qlopt_arg_mismatch(L, pos, flags, argtype, msg_head);
        *(QDP_D3_ColorMatrix **)p = qlua_checkLatColMat3(L, pos, NULL, 3)->ptr;
        return 0; 
    }
    case QLOPT_LATDIRFERM   : {
        if (qlua_qtype(L, pos) != qLatDirFerm3)
            return qlopt_arg_mismatch(L, pos, flags, argtype, msg_head);
        *(QDP_D3_DiracFermion **)p = qlua_checkLatDirFerm3(L, pos, NULL, 3)->ptr;
        return 0;
    }
    case QLOPT_LATDIRPROP   : {
        if (qlua_qtype(L, pos) != qLatDirProp3)
            return qlopt_arg_mismatch(L, pos, flags, argtype, msg_head);
        *(QDP_D3_DiracPropagator **)p = qlua_checkLatDirProp3(L, pos, NULL, 3)->ptr;
        return 0;
    }
    case QLOPT_GAMMA        : {
        if (qlua_qtype(L, pos) != qGamma)
            return qlopt_arg_mismatch(L, pos, flags, argtype, msg_head);
        *(mGamma **)p = qlua_checkClifford(L, pos)->g;
        return 0; 
    }
    case QLOPT_TABLE        : {
        if (lua_type(L, pos) != LUA_TTABLE) 
            qlopt_arg_mismatch(L, pos, flags, argtype, msg_head);
        return qlopt_read_table_internal(L, pos, t.tab, base, msg_head);
    }
    case QLOPT_ARRAY        : {
        if (NULL != base)
            return luaL_error(L, "%s cannot parse array within array", msg_head); 
        if (lua_type(L, pos) != LUA_TTABLE) 
            qlopt_arg_mismatch(L, pos, flags, argtype, msg_head);
        if (0 != (status = qlopt_read_array_checkdim_alloc(L, pos, t.arr, msg_head)))
            return status;
        if (NULL == t.arr->buf)
            luaL_error(L, "%s allocation failed", msg_head);
        /* call top-level read_array */
        return qlopt_read_array_internal(L, 0, pos, t.arr, t.arr->buf, msg_head); 
    }
    case QLOPT_SUBSET       :
        return luaL_error(L, "%s unsupported opt.type %s", msg_head, 
                          qlopt_type_name(argtype));
    default : return luaL_error(L, "%s unsupported arg.type %d", msg_head, argtype);
    }
}

/* check dimensions of an array on stack
   * fill those that are undef (<=0)
   * check that dim_i does not exceed maxdim_i
   * if buf == NULL, allocate space (prod(a->dim) * a->argsize), 
     save pointer to a->buf 
 */
static int 
qlopt_read_array_checkdim_alloc(lua_State *L, int pos, qloptNdarray *a, const char *msg_head)
{
    char strbuf[1024];
    if (NULL == msg_head) {
        snprintf(strbuf, sizeof(strbuf), "(arg %d)", pos);
        msg_head = strbuf;
    }
    /* FIXME argument # when reporting errors must use optional info from 
       callers if parsing a nested argument list */
    if (lua_gettop(L) < pos) {
        luaL_error(L, "%s expect %d-d array at %d", msg_head, a->ndim, pos);
        return 1;
    }
    int i = 0;
    int i_pos = pos;
    do {
        if (lua_type(L, i_pos) != LUA_TTABLE) 
            return luaL_error(L, "%s expect %d-d array", 
                              msg_head, a->ndim);
        int dim_i = lua_objlen(L, i_pos);
        if (dim_i <= 0) 
            return luaL_error(L, "%s bad %d-nested array len = %d", 
                              msg_head, i, dim_i);
        if (0 < a->dim[i] && dim_i != a->dim[i])
            return luaL_error(L, "%s %d-nested array len mismatch: %d != %d(expected)", 
                              msg_head, i, dim_i, a->dim[i]);
        if(NULL != a->maxdim && 0 < a->maxdim[i] && a->maxdim[i] < dim_i)
            return luaL_error(L, "%s %d-nested array len exceeds limit: %d > %d(max)",
                              msg_head, i, dim_i, a->maxdim[i]);
        /* next nested array */
        a->dim[i] = dim_i;
        i++;
        lua_pushinteger(L, 1);
        lua_rawget(L, i_pos);
        i_pos = lua_gettop(L);
    } while (i < a->ndim);

    /* clean up exploratory pushes */
    lua_pop(L, i);

    /* allocate array if necessary */
    int nmemb = prod_int_(a->ndim, a->dim);
    if (! (QLOPT_ARRAY == a->argtype 
            || QLOPT_TABLE == a->argtype
            || QLOPT_NONE == a->argtype)) {
        if (qlopt_type_size(a->argtype) != a->argsize)
            return luaL_error(L, "%s internal type size mismatch ; submit bug report",
                        msg_head);
    }
    int bufsize = nmemb * a->argsize;
    if (NULL == a->buf) {
        a->buf = malloc(bufsize);
        if (NULL == a->buf)
            return luaL_error(L, "%s not enough memory (try alloc %llu)", 
                              msg_head, (size_t)(bufsize));
        memset(a->buf, 0, bufsize); 
    } 
    /* if fill value specified, fill the buffer */
    if (NULL != a->fill) {
        for (int i = 0 ; i < nmemb ; i++)
            memcpy(a->buf + a->argsize * i, a->fill, a->argsize);
    }
    return 0;
}

/* read array from stack into a->buf 
   * pos_i is the location of the current element on stack 
     (either a value or a nested array)
   * a->buf must be allocated and have the right size (by checkdim_alloc)
   * a->dim must all be defined (by checkdim_alloc)
   * buf_i points to the current location in a->buf 
     (a->buf is not used directly)
*/
static int 
qlopt_read_array_internal(lua_State *L, int i, int pos_i, qloptNdarray *a, 
        void *buf_i, const char *msg_head)
{
    int status = 0;
    char strbuf0[128], strbuf[1024];
    if (NULL == msg_head) {
        snprintf(strbuf0, sizeof(strbuf0), "(arg %d)", pos_i);
        msg_head = strbuf0;
    }
    if (NULL == buf_i)
        return luaL_error(L, "%s buffer has not been allocated", msg_head);
    if (i < a->ndim) {
        if (a->dim[i] <= 0) 
            return luaL_error(L, "%s array dim[%d] is undefinded", msg_head, i);
        if (LUA_TTABLE != lua_type(L, pos_i))
            return luaL_error(L, "%s array expected", msg_head);
        if (lua_objlen(L, pos_i) != a->dim[i])
            return luaL_error(L, "%s bad array dimension", msg_head);
        /* the remaining fastest-changing dimensions */
        int stride_memb = prod_int_(a->ndim - i - 1, a->dim + i + 1);
        for (int k = 0 ; k < a->dim[i] ; k++) {
            /* message prefix for nested parsing */
            if (NULL == msg_head) 
                snprintf(strbuf, sizeof(strbuf), "%d", k);
            else 
                snprintf(strbuf, sizeof(strbuf), "%s:%d", msg_head, k);

            lua_pushinteger(L, k + 1); /* sic! Lua numbering */
            lua_rawget(L, pos_i);
            if (0 != (status = qlopt_read_array_internal(L, i + 1, lua_gettop(L), a,
                        buf_i + k * stride_memb * a->argsize, strbuf))) {
                lua_pop(L, 1);
                return status;
            }
            lua_pop(L, 1);
        }
        return status;
    } else if (i == a->ndim) {
        /* TODO special handling for TABLE : give elem_list in #5 */
        qloptData val;
        if (QLOPT_TABLE == a->argtype)
            val.tab = a->elem_list ;
        else
            val.off = 0 ;
        return qlopt_read_arg_internal(L, pos_i, 0, a->argtype, val, buf_i, msg_head);
    } else /* should not be here */
        return luaL_error(L, "%s internal error: submit bug report", msg_head);
}
/* parse a table from Lua stack */
static int
qlopt_read_table_internal(lua_State *L, int pos, 
            qloptElem *elem_list, void *base, const char *msg_head)
{
    char strbuf0[128], strbuf[1024];
    if (NULL == msg_head) {
        snprintf(strbuf0, sizeof(strbuf0), "(arg %d)", pos);
        msg_head = strbuf0;
    }
    if (lua_gettop(L) < pos 
            || LUA_TTABLE != lua_type(L, pos))
        return luaL_error(L, "%s table expected", msg_head);


    qloptElem *elem = elem_list;
    int pos_next = 1;
    while (QLOPT_END != elem->pos) {
        if (QLOPT_KEY == elem->pos) {
            snprintf(strbuf, sizeof(strbuf), "%s:%s", msg_head, elem->key);
            lua_pushstring(L, elem->key);
            lua_rawget(L, pos);
        } else if (QLOPT_NEXT == elem->pos) {
            if (QLOPT_KEY == pos_next)
                return luaL_error(L, "%s no arg next to key", msg_head);
            if (pos_next <= 0)
                return luaL_error(L, "%s bad pos_next=%d", msg_head, pos_next);
            snprintf(strbuf, sizeof(strbuf), "%s:%d", msg_head, pos_next);
            lua_pushinteger(L, pos_next);
            lua_rawget(L, pos);
            pos_next ++;
        } else {
            if (elem->pos <= 0)
                return luaL_error(L, "%s bad pos=%d", msg_head, elem->pos);
            snprintf(strbuf, sizeof(strbuf), "%s:%d", msg_head, elem->pos);
            lua_pushinteger(L, elem->pos);
            lua_rawget(L, pos);
            pos_next = elem->pos + 1;
        }
        int arg_pos = lua_gettop(L);
        int status = qlopt_read_arg_internal(L, arg_pos, elem->flags, 
                elem->argtype, elem->val, base, strbuf);
        if (status)
            return luaL_error(L, "%s : cannot parse argument", strbuf);
        lua_pop(L, 1);
        elem++;
    }
    return 0;
}

int 
qlopt_read_array(lua_State *L, int pos, qloptNdarray *a)
{
    int status = 0;
    char strbuf[128];
    snprintf(strbuf, sizeof(strbuf), "(array at %d)", pos);
    if (0 != (status = qlopt_read_array_checkdim_alloc(L, pos, a, strbuf))) {
        return status;
    }
    return qlopt_read_array_internal(L, 0, pos, a, a->buf, strbuf);
}

int 
qlopt_read_table(lua_State *L, int pos, qloptElem *elem_list)
{
    int status = 0;
    char strbuf[128];
    snprintf(strbuf, sizeof(strbuf), "(table at %d)", pos);
    return qlopt_read_table_internal(L, pos, elem_list, NULL, strbuf);
}

int 
qlopt_read_stack(lua_State *L, qloptElem *elem_list)
{
    int status = 0;
    char strbuf[128];
    int pos_next = 1,
        pos;
    int i_elem = 1;
    qloptElem *elem = elem_list;
    while (QLOPT_END != elem->pos) {
        if (QLOPT_KEY == elem->pos)
            return luaL_error(L, "(elem %d) cannot use keys to ref.to args", i_elem);
        if (QLOPT_NEXT == elem->pos) {
            if (QLOPT_KEY == pos_next)
                return luaL_error(L, "(elem %d) no arg next to key", i_elem);
            if (pos_next <= 0)
                return luaL_error(L, "(elem %d) bad pos_next=%d", i_elem, pos_next);
            pos = pos_next;
            pos_next ++;
        } else {
            if (elem->pos <= 0)
                return luaL_error(L, "(elem %d) bad pos=%d", i_elem, elem->pos);
            pos = elem->pos;
            pos_next = pos + 1;
        }
        snprintf(strbuf, sizeof(strbuf), "(arg %d)", pos);
        int status = qlopt_read_arg_internal(L, pos, elem->flags, 
                elem->argtype, elem->val, NULL, strbuf);
        if (status)
            return luaL_error(L, "%s : cannot parse argument", strbuf);

        elem++;
        i_elem++;
    }
    return 0;

}

#define printf0(...) do { if (QDP_this_node == qlua_master_node) printf(__VA_ARGS__); } while(0)
static void 
print_arr_qloptint(int ndim, int *dim, qlopt_Int *arr, const char *msg_head) 
{
    char strbuf[1024];
    if (0 == ndim)
        printf0("%s = %lld\n", msg_head, (long long int)(*arr));
    else if (0 < ndim) {
        int memb_stride = prod_int_(ndim - 1, dim + 1);
        for (int i = 0 ; i < *dim ; i++) {
            snprintf(strbuf, sizeof(strbuf), "%s[%d]", msg_head, i);
            print_arr_qloptint(ndim - 1, dim + 1, arr + memb_stride * i, strbuf);
        }
    }
    /* if ndim < 0, do nothing */
}
static void 
print_arr_qloptreal(int ndim, int *dim, qlopt_Real *arr, const char *msg_head) 
{
    char strbuf[1024];
    if (0 == ndim)
        printf0("%s = %e\n", msg_head, (double)(*arr));
    else if (0 < ndim) {
        int memb_stride = prod_int_(ndim - 1, dim + 1);
        for (int i = 0 ; i < *dim ; i++) {
            snprintf(strbuf, sizeof(strbuf), "%s[%d]", msg_head, i);
            print_arr_qloptreal(ndim - 1, dim + 1, arr + memb_stride * i, strbuf);
        }
    }
    /* if ndim < 0, do nothing */
}
static void 
print_arr_qloptcomplex(int ndim, int *dim, qlopt_Complex *arr, const char *msg_head) 
{
    char strbuf[1024];
    if (0 == ndim)
        printf0("%s = (%e, %e)\n", msg_head, creal(*arr), cimag(*arr));
    else if (0 < ndim) {
        int memb_stride = prod_int_(ndim - 1, dim + 1);
        for (int i = 0 ; i < *dim ; i++) {
            snprintf(strbuf, sizeof(strbuf), "%s[%d]", msg_head, i);
            print_arr_qloptcomplex(ndim - 1, dim + 1, arr + memb_stride * i, strbuf);
        }
    }
    /* if ndim < 0, do nothing */
}
static void
print_arr_general(int ndim, int *dim, void *buf, size_t argsize, 
        void (*print_rec)(void *buf, const char *msg_head),  /* print a record in 1 line */ 
        const char *msg_head)
{
    char strbuf[1024];
    if (0 == ndim)
        print_rec(buf, msg_head);
    else if (0 < ndim) {
        int memb_stride = prod_int_(ndim - 1, dim + 1);
        for (int i = 0 ; i < *dim ; i++) {
            snprintf(strbuf, sizeof(strbuf), "%s[%d]", msg_head, i);
            print_arr_general(ndim - 1, dim + 1, buf + memb_stride * argsize * i, 
                    argsize, print_rec, strbuf);
        }
    }
    /* if ndim < 0, do nothing */
}



static void
print_dim(int n, const int dim[], const char *msg_head)
{
    printf0("%s@{", msg_head);
    for (int i = 0 ; i < n ; i++)
        printf0("%d ", dim[i]);
    printf0("}\n");
}

/* example for reading an array of records */
typedef struct recx_s {
    QDP_D_Complex   *cplx;
    char            *name;
    QLA_Real        bc_t;
} recx;
qloptElem qlopt_tab_recx[] = {
    { 1, NULL, QLOPT_STRING,    0, {.off = offsetof(recx, name) } },
    { 2, NULL, QLOPT_LATCOMPLEX,0, {.off = offsetof(recx, cplx) } },
    { 3, NULL, QLOPT_REAL,      0, {.off = offsetof(recx, bc_t) } },
    { QLOPT_END }
};

static void
print_recx(void *a_, const char *msg_head) 
{
    recx *a = (recx *)a_;
    QLA_D_Real n2 = 0.;
    QDP_D_r_eq_norm2_C(&n2, a->cplx, QDP_all_L(QDP_D_get_lattice_C(a->cplx)));
    printf0("%s = {name='%s'  f=%p(|f|2=%e)  bc_t=%f}\n", 
            msg_head, a->name, 
            (void*)a->cplx, (double)n2,
            (double)(a->bc_t));
}
static void
print_str(void *a_, const char *msg_head)
{
    printf0("%s = '%s'\n", msg_head, *(char **)a_);
}

int 
q_test_qlopt(lua_State *L)
{
    int i_arr[3];
    double d_arr[5];
    const char *s_arr[7];
    void *xx = NULL;

    /*
    test_parse(
      l_nev, l_ncv, l_maxiter, l_tol,
      { ["eigcg"] = {eigcg_vmax, eigcg_nev, eigcg_tol, eigcg_umax},
        ["cheb_accel"] = { l_poly_n, l_poly_a, l_poly_b, l_poly_x0 },
        ["coeff_list"]  = coeffs_list_Nx2,
        ["mom_list"]    = mom_list_Mx4,
        ["xy"] = xy_2x3x4,
        ["which"] = "LR", -- need the largest ev of polynomial T_n(A)
        ["arpack_logfile"] = arpack_logfile(cfg_key),
        ["inplace"] = true
        }
        abc_list)
    */
    qlopt_Int l_nev = -1,
              l_ncv = -1,
              l_maxiter = -1;
    qlopt_Real l_tol = -1.;


    /* variables to fill, with default values ; 
       while default values are not necessary for mandatory options, 
       the table itself will be mandatory, so better initialize */

    qlopt_Int eigcg_vmax  = -1, 
              eigcg_nev   = -1, 
              eigcg_umax  = -1;
    qlopt_Real eigcg_tol = -1.;
    qloptElem qlopt_tab_eigcg[] = {
        { 1, NULL, QLOPT_INT,   0, { .i = &eigcg_vmax } },
        { 2, NULL, QLOPT_INT,   0, { .i = &eigcg_nev  } },
        { 3, NULL, QLOPT_REAL,  0, { .r = &eigcg_tol  } },
        { 4, NULL, QLOPT_INT,   0, { .i = &eigcg_umax } },
        { QLOPT_END }
    };
    
    qlopt_Int l_n = -1;
    qlopt_Real l_a  = DBL_MAX,
           l_b  = DBL_MAX,
           l_x0 = DBL_MAX;
    qloptElem qlopt_tab_cheb_accel[] = {
        { 1, NULL, QLOPT_INT,   0,   { .i = &l_n } },
        { 2, NULL, QLOPT_REAL,  0,   { .r = &l_a } },
        { 3, NULL, QLOPT_REAL,  0,   { .r = &l_b } },
        { 4, NULL, QLOPT_REAL,  QLOPT_OPT,   { .r = &l_x0 } },
        { QLOPT_END }
    };
    
    int coeffs_dim[]    = { -1, 2 },
        coeffs_dim_max[] = { 10, -1 }; 
    qloptNdarray qlopt_arr_coeffs = { 2, coeffs_dim, coeffs_dim_max, 
                QLOPT_REAL, NULL, sizeof(qlopt_Real), NULL };

    int mom_dim[] = {-1, 4},
        mom_dim_max[] = { 10, -1};
    qloptNdarray qlopt_arr_mom = qlopt_array_scalar(2, mom_dim, mom_dim_max, 
                QLOPT_INT, NULL, NULL);

    int xy_dim[] = { 2, 3, 4 };
    qlopt_Int xy[2][3][4];
    qloptNdarray qlopt_arr_xy = qlopt_array_scalar(3, xy_dim, NULL, 
                QLOPT_INT, xy, NULL);

    const char *which = NULL;
    const char *arpack_logfile = NULL;
    qlopt_Int inplace = 0;
    
    int abc_dim[]    = { -1, -1, -1 },
        abc_dim_max[] = { 10, 10, 10 };
    qloptNdarray qlopt_arr_abc = qlopt_array_scalar(3, abc_dim, abc_dim_max, 
                QLOPT_COMPLEX, NULL, NULL);

    int recx_dim[] = { -1, -1 };
    int recx_dim_max[] = { -1, -1 };
    qloptNdarray qlopt_arr_recx = qlopt_array_struct(2, recx_dim, recx_dim_max, 
                qlopt_tab_recx, sizeof(recx), NULL, NULL);

    int str_dim[]     = { -1, -1 },
        str_dim_max[] = {  4,  4 };
    const char *str_arr_buf[16];
    qloptNdarray qlopt_arr_str = qlopt_array_scalar(2, str_dim, str_dim_max, 
                QLOPT_STRING, str_arr_buf, NULL);


    qloptElem qlopt_tab_opt[] = {
        { QLOPT_KEY, "eigcg",       QLOPT_TABLE,  0, { .tab = qlopt_tab_eigcg } },
        { QLOPT_KEY, "cheb_accel",  QLOPT_TABLE,  0, { .tab = qlopt_tab_cheb_accel } },
        { QLOPT_KEY, "xy",          QLOPT_ARRAY,  0, { .arr = &qlopt_arr_xy } },
        { QLOPT_KEY, "mom_list",    QLOPT_ARRAY,  0, { .arr = &qlopt_arr_mom } },
        { QLOPT_KEY, "coeff_list",  QLOPT_ARRAY,  0, { .arr = &qlopt_arr_coeffs } },
        { QLOPT_KEY, "abc",         QLOPT_ARRAY,  0, { .arr = &qlopt_arr_abc } },
        { QLOPT_KEY, "which",       QLOPT_STRING, 0, { .s = &which } },
        { QLOPT_KEY, "arpack_logfile", QLOPT_STRING, QLOPT_OPT,
                    { .s = &arpack_logfile } },
        { QLOPT_KEY, "inplace",     QLOPT_BOOL, QLOPT_OPT,
                    { .i = &inplace } },
        { QLOPT_END }
    } ;
    
    qloptElem stack_elem[] = {
        { QLOPT_NEXT, NULL, QLOPT_ARRAY, 0, { .arr = &qlopt_arr_recx } },
        { QLOPT_NEXT, NULL, QLOPT_INT,   0, { .i = &l_nev } },
        { QLOPT_NEXT, NULL, QLOPT_INT,   0, { .i = &l_ncv } },
        { QLOPT_NEXT, NULL, QLOPT_INT,   0, { .i = &l_maxiter } },
        { QLOPT_NEXT, NULL, QLOPT_REAL,  0, { .r = &l_tol } },
        { QLOPT_NEXT, NULL, QLOPT_ARRAY, 0, { .arr = &qlopt_arr_str } },
        { QLOPT_NEXT, NULL, QLOPT_TABLE, QLOPT_OPT, { .tab = qlopt_tab_opt } },
        { QLOPT_END }
    };

    int status = 0;
    
    printf0("START lua_gettop(L) = %d\n", lua_gettop(L));
    status = qlopt_read_stack(L, stack_elem);
    printf0("qlopt read status = %d\n", status);
    printf0("STOP lua_gettop(L) = %d\n", lua_gettop(L));
    /* print parsed values */
    printf0("l_nev=%d  l_ncv=%d  l_maxiter=%d  l_tol=%e\n", 
            l_nev, l_ncv, l_maxiter, l_tol);
    printf0("eigcg = {%d, %d, %e, %d}\n", 
            eigcg_vmax, eigcg_nev, eigcg_tol, eigcg_umax);
    printf0("cheb_accel = {%d, %e, %e, %e}\n",
            l_n, l_a, l_b, l_x0);
    
    printf0("which='%s'  arpack_logfile='%s'  inplace=%d\n", which, arpack_logfile, (int)inplace);

    print_dim(2, coeffs_dim, "coeffs");
    print_arr_qloptreal(2, coeffs_dim, qlopt_arr_coeffs.buf, "coeffs");

    print_dim(2, mom_dim, "mom");
    print_arr_qloptint(2, mom_dim, qlopt_arr_mom.buf, "mom");
    
    print_dim(3, xy_dim, "xy");
    print_arr_qloptint(3, xy_dim, (qlopt_Int *)xy, "xy");

    print_dim(3, abc_dim, "abc");
    print_arr_qloptcomplex(3, abc_dim, qlopt_arr_abc.buf, "abc");
    
    print_dim(2, recx_dim, "recx");
    print_arr_general(2, recx_dim, qlopt_arr_recx.buf, sizeof(recx), print_recx, "recx");
    
    print_dim(2, str_dim, "str");
    print_arr_general(2, str_dim, str_arr_buf, sizeof(char *), print_str, "str");
  
#if 0
    /* 1st arg : integer */
    status = qlua_parse_scalar(&l_nev, QLOPT_INT, 0, 1); 
    
    /* 2nd arg : something */
    /* just a wrapper : override position in the structure, call appropriate function depending on the argtype */
    status = qlua_parse_stackitem(&parse_tab_arg[1], 2); 

    /* all args: arg_first, arg_last ?*/
    status = qlua_parse_arglist(parse_tab_arg/*, arg_first, arg_last ?*/);

    /* 5th arg : table */
    qlua_parse_tab(parse_tab_opt, 5);

    /* 6th arg : 3d array */
    abc = qlua_parse_array(3, abc_dim, abc_maxdim, abc, 6);

#endif
    return 0;
}

qloptNdarray 
qlopt_array_scalar(int ndim, int *dim, int *maxdim, 
                qloptType argtype, void *buf, void *fill)
{
    assert(QLOPT_ARRAY != argtype);
    assert(QLOPT_TABLE != argtype);
    assert(QLOPT_NONE != argtype);
    assert(NULL != dim);
    assert(0 < ndim);
    size_t argsize = qlopt_type_size(argtype);
    assert(argsize < SIZE_MAX);
    qloptNdarray res;
    res.ndim        = ndim;
    res.dim         = dim;
    res.maxdim      = maxdim; 
    res.argtype     = argtype;
    res.buf         = buf;
    res.argsize     = argsize;
    res.fill        = fill;
    res.elem_list   = NULL;
    return res;
}

qloptNdarray 
qlopt_array_struct(int ndim, int *dim, int *maxdim, 
                qloptElem *elem_list, size_t argsize, void *buf, void *fill)
{
    assert(NULL != dim);
    assert(0 < ndim);
    qloptNdarray res;
    res.ndim        = ndim;
    res.dim         = dim;
    res.maxdim      = maxdim; 
    res.argtype     = QLOPT_TABLE;
    res.buf         = buf;
    res.argsize     = argsize;
    res.fill        = fill;
    res.elem_list   = elem_list;
    return res;
}
