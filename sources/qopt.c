#include "qopt.h"                                                    /* DEPS */
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



/* qopt_read_type_internal - read generic element from stack
   base == NULL : use t->p as the destination
   base != NULL : use base + t->off as the destination
 */
static const char * 
qopt_type_name(QLUA_Type argtype)
{
    switch(argtype) {
    case qOptNone         : return "nil";
    case qOptInt          : return "int";
    case qOptBool         : return "bool";
    case qReal         : return "real";
    case qComplex      : return "complex";
    case qString       : return "string";
    case qLatInt       : return "LatInt";
    case qLatReal      : return "LatReal";
    case qLatComplex   : return "LatComplex";
    case qLatColVec3    : return "LatColVec";
    case qLatColMat3    : return "LatColMat";
    case qLatDirFerm3   : return "LatDirFerm";
    case qLatDirProp3   : return "LatDirProp";
    case qGamma        : return "Gamma";
    case qLatSubset       : return "Subset";
    case qOptStruct        : return "table";
    case qOptArray        : return "array";
    default : return "(unknown)";
    }
}

static size_t
qopt_type_size(QLUA_Type argtype)
{
    switch(argtype) {
    case qOptInt          : 
    case qOptBool         : 
        return sizeof(qopt_Int);

    case qReal         : 
        return sizeof(qopt_Real);

    case qComplex      : 
        return sizeof(qopt_Complex);

    case qString       : 
    case qLatInt       : 
    case qLatReal      : 
    case qLatComplex   : 
    case qLatColVec3    : 
    case qLatColMat3    : 
    case qLatDirFerm3   : 
    case qLatDirProp3   : 
    case qGamma        : 
    case qLatSubset       : 
        return sizeof(void *);  /* sic! all pointers of the same size */

    case qOptNone         : 
    case qOptStruct        : 
    case qOptArray        : 
    default : 
        return SIZE_MAX ;       /* size is meaningless, must be handled separately */
    }
}

/* invoke on type mismatch to print warning or terminate, depending on flags */
static int
qopt_arg_mismatch(lua_State *L, int pos, int flags, 
        QLUA_Type argtype, const char *msg_head)
{
    if (flags & QOPT_BAD_SKIP) return 0;
    if (flags & QOPT_BAD_WARN) {
        if (QDP_this_node == qlua_master_node)
            printf("%s mismatch: expect %s", msg_head, 
                   qopt_type_name(argtype));
        return 0;
    }
    return luaL_error(L, "%s mismatch: expect %s", msg_head, 
                      qopt_type_name(argtype));
}
/* invoke on missing arg to print warning or terminate, depending on flags */
static int
qopt_arg_missing(lua_State *L, int pos, int flags, 
        QLUA_Type argtype, const char *msg_head)
{
    if (flags & QOPT_OPT) return 0;
    if (flags & QOPT_OPT_WARN) {
        if (QDP_this_node == qlua_master_node) 
            printf("%s missing: expect %s", msg_head, 
                   qopt_type_name(argtype));
        return 0;
    }
    return luaL_error(L, "%s missing: expect %s", msg_head, 
                      qopt_type_name(argtype));
}

static int 
qopt_read_arg_internal(lua_State *L, int pos, int flags, 
        QLUA_Type argtype, qoptData t, void *base,
        const char *msg_head);

static int
qopt_read_table_internal(lua_State *L, int pos, 
            qoptElem *elem_list, void *base, const char *msg_head);

static int 
qopt_read_array_checkdim_alloc(lua_State *L, int pos, 
        qoptArray *a, const char *msg_head);
static int 
qopt_read_array_internal(lua_State *L, int i, int pos_i, qoptArray *a, 
        void *buf_i, const char *msg_head);


/* parse an argument from a Lua stack
   XXX strict typing is imposed in some cases, e.g. parsing will fail 
   when string <-> number conversion is necessary; this conversion is 
   typically performed by lua_tonumber and lua_tostring that have side effects 
 */
static int 
qopt_read_arg_internal(lua_State *L, int pos, int flags, 
        QLUA_Type argtype, qoptData t, void *base,
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
             && qOptNone != argtype)) 
        return qopt_arg_missing(L, pos, flags, argtype, msg_head);

    void *p = (NULL != base ? base + t.off : t.ptr);
    switch(argtype) {
    case qOptNone         : return 0;
    case qOptBool         : {
        if (lua_type(L, pos) != LUA_TBOOLEAN)
            return qopt_arg_mismatch(L, pos, flags, argtype, msg_head);
        *(qopt_Int *)p = lua_toboolean(L, pos);
        return 0; 
    }
    case qOptInt          : {
        if (lua_type(L, pos) != LUA_TNUMBER) /* using lua_isnumber may be slower */
            return qopt_arg_mismatch(L, pos, flags, argtype, msg_head);
        *(qopt_Int *)p = luaL_checkinteger(L, pos);
        return 0; 
    }
    case qReal         : {
        if (lua_type(L, pos) != LUA_TNUMBER) /* using lua_isnumber may be slower */
            return qopt_arg_mismatch(L, pos, flags, argtype, msg_head);
        *(qopt_Real *)p    = luaL_checknumber(L, pos);
        return 0; 
    }
    case qComplex      : {
        if (qlua_qtype(L, pos) != qComplex)
            return qopt_arg_mismatch(L, pos, flags, argtype, msg_head);
        QLA_D_Complex *c = qlua_checkComplex(L, pos);
        *(double complex *)p = QLA_real(*c) + I*QLA_imag(*c);
        return 0; 
    }
    case qString       : {
        if (lua_type(L, pos) != LUA_TSTRING)
            return qopt_arg_mismatch(L, pos, flags, argtype, msg_head);
        *(const char **)p = luaL_checkstring(L, pos);
        return 0; 
    }
    case qLatInt       : {
        if (qlua_qtype(L, pos) != qLatInt)
            return qopt_arg_mismatch(L, pos, flags, argtype, msg_head);
        *(QDP_Int **)p = qlua_checkLatInt(L, pos, NULL)->ptr;
        return 0; 
    }
    case qLatReal      : {
        if (qlua_qtype(L, pos) != qLatReal)
            return qopt_arg_mismatch(L, pos, flags, argtype, msg_head);
        *(QDP_Real **)p = qlua_checkLatReal(L, pos, NULL)->ptr;
        return 0; 
    }
    case qLatComplex   : { 
        if (qlua_qtype(L, pos) != qLatComplex)
            return qopt_arg_mismatch(L, pos, flags, argtype, msg_head);
        *(QDP_Complex **)p = qlua_checkLatComplex(L, pos, NULL)->ptr;
        return 0; 
    }
    /* TODO iterate over Nc */
    case qLatColVec3    : {
        if (qlua_qtype(L, pos) != qLatColVec3)
            return qopt_arg_mismatch(L, pos, flags, argtype, msg_head);
        *(QDP_D3_ColorVector **)p = qlua_checkLatColVec3(L, pos, NULL, 3)->ptr;
        return 0; 
    }
    case qLatColMat3    : {
        if (qlua_qtype(L, pos) != qLatColMat3)
            return qopt_arg_mismatch(L, pos, flags, argtype, msg_head);
        *(QDP_D3_ColorMatrix **)p = qlua_checkLatColMat3(L, pos, NULL, 3)->ptr;
        return 0; 
    }
    case qLatDirFerm3   : {
        if (qlua_qtype(L, pos) != qLatDirFerm3)
            return qopt_arg_mismatch(L, pos, flags, argtype, msg_head);
        *(QDP_D3_DiracFermion **)p = qlua_checkLatDirFerm3(L, pos, NULL, 3)->ptr;
        return 0;
    }
    case qLatDirProp3   : {
        if (qlua_qtype(L, pos) != qLatDirProp3)
            return qopt_arg_mismatch(L, pos, flags, argtype, msg_head);
        *(QDP_D3_DiracPropagator **)p = qlua_checkLatDirProp3(L, pos, NULL, 3)->ptr;
        return 0;
    }
    case qGamma        : {
        if (qlua_qtype(L, pos) != qGamma)
            return qopt_arg_mismatch(L, pos, flags, argtype, msg_head);
        *(mGamma **)p = qlua_checkClifford(L, pos)->g;
        return 0; 
    }
    case qOptStruct        : {
        if (lua_type(L, pos) != LUA_TTABLE) 
            qopt_arg_mismatch(L, pos, flags, argtype, msg_head);
        return qopt_read_table_internal(L, pos, t.tab, base, msg_head);
    }
    case qOptArray        : {
        if (NULL != base)
            return luaL_error(L, "%s cannot parse array within array", msg_head); 
        if (lua_type(L, pos) != LUA_TTABLE) 
            qopt_arg_mismatch(L, pos, flags, argtype, msg_head);
        if (0 != (status = qopt_read_array_checkdim_alloc(L, pos, t.arr, msg_head)))
            return status;
        if (NULL == t.arr->buf)
            luaL_error(L, "%s allocation failed", msg_head);
        /* call top-level read_array */
        return qopt_read_array_internal(L, 0, pos, t.arr, t.arr->buf, msg_head); 
    }
    case qLatSubset       :
        return luaL_error(L, "%s unsupported opt.type %s", msg_head, 
                          qopt_type_name(argtype));
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
qopt_read_array_checkdim_alloc(lua_State *L, int pos, qoptArray *a, const char *msg_head)
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
    if (! (qOptArray == a->argtype 
            || qOptStruct == a->argtype
            || qOptNone == a->argtype)) {
        if (qopt_type_size(a->argtype) != a->argsize)
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
qopt_read_array_internal(lua_State *L, int i, int pos_i, qoptArray *a, 
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
            if (0 != (status = qopt_read_array_internal(L, i + 1, lua_gettop(L), a,
                        buf_i + k * stride_memb * a->argsize, strbuf))) {
                lua_pop(L, 1);
                return status;
            }
            lua_pop(L, 1);
        }
        return status;
    } else if (i == a->ndim) {
        /* TODO special handling for TABLE : give elem_list in #5 */
        qoptData val;
        if (qOptStruct == a->argtype)
            val.tab = a->elem_list ;
        else
            val.off = 0 ;
        return qopt_read_arg_internal(L, pos_i, 0, a->argtype, val, buf_i, msg_head);
    } else /* should not be here */
        return luaL_error(L, "%s internal error: submit bug report", msg_head);
}
/* parse a table from Lua stack */
static int
qopt_read_table_internal(lua_State *L, int pos, 
            qoptElem *elem_list, void *base, const char *msg_head)
{
    char strbuf0[128], strbuf[1024];
    if (NULL == msg_head) {
        snprintf(strbuf0, sizeof(strbuf0), "(arg %d)", pos);
        msg_head = strbuf0;
    }
    if (lua_gettop(L) < pos 
            || LUA_TTABLE != lua_type(L, pos))
        return luaL_error(L, "%s table expected", msg_head);


    qoptElem *elem = elem_list;
    int pos_next = 1;
    while (QOPT_END != elem->pos) {
        if (QOPT_KEY == elem->pos) {
            snprintf(strbuf, sizeof(strbuf), "%s:%s", msg_head, elem->key);
            lua_pushstring(L, elem->key);
            lua_rawget(L, pos);
        } else if (QOPT_NEXT == elem->pos) {
            if (QOPT_KEY == pos_next)
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
        int status = qopt_read_arg_internal(L, arg_pos, elem->flags, 
                elem->argtype, elem->val, base, strbuf);
        if (status)
            return luaL_error(L, "%s : cannot parse argument", strbuf);
        lua_pop(L, 1);
        elem++;
    }
    return 0;
}

int 
qopt_read_array(lua_State *L, int pos, qoptArray *a)
{
    int status = 0;
    char strbuf[128];
    snprintf(strbuf, sizeof(strbuf), "(array at %d)", pos);
    if (0 != (status = qopt_read_array_checkdim_alloc(L, pos, a, strbuf))) {
        return status;
    }
    return qopt_read_array_internal(L, 0, pos, a, a->buf, strbuf);
}

int 
qopt_read_table(lua_State *L, int pos, qoptElem *elem_list)
{
    char strbuf[128];
    snprintf(strbuf, sizeof(strbuf), "(table at %d)", pos);
    return qopt_read_table_internal(L, pos, elem_list, NULL, strbuf);
}

int 
qopt_read_stack(lua_State *L, qoptElem *elem_list)
{
    char strbuf[128];
    int pos_next = 1,
        pos;
    int i_elem = 1;
    qoptElem *elem = elem_list;
    while (QOPT_END != elem->pos) {
        if (QOPT_KEY == elem->pos)
            return luaL_error(L, "(elem %d) cannot use keys to ref.to args", i_elem);
        if (QOPT_NEXT == elem->pos) {
            if (QOPT_KEY == pos_next)
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
        int status = qopt_read_arg_internal(L, pos, elem->flags, 
                elem->argtype, elem->val, NULL, strbuf);
        if (status)
            return luaL_error(L, "%s : cannot parse argument", strbuf);

        elem++;
        i_elem++;
    }
    return 0;

}



qoptArray 
qopt_array_scalar(int ndim, int *dim, int *maxdim, 
                QLUA_Type argtype, void *buf, void *fill)
{
    assert(qOptArray != argtype);
    assert(qOptStruct != argtype);
    assert(qOptNone != argtype);
    assert(NULL != dim);
    assert(0 < ndim);
    size_t argsize = qopt_type_size(argtype);
    assert(argsize < SIZE_MAX);
    qoptArray res;
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

qoptArray 
qopt_array_struct(int ndim, int *dim, int *maxdim, 
                qoptElem *elem_list, size_t argsize, void *buf, void *fill)
{
    assert(NULL != dim);
    assert(0 < ndim);
    qoptArray res;
    res.ndim        = ndim;
    res.dim         = dim;
    res.maxdim      = maxdim; 
    res.argtype     = qOptStruct;
    res.buf         = buf;
    res.argsize     = argsize;
    res.fill        = fill;
    res.elem_list   = elem_list;
    return res;
}


