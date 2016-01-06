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
qlopt_read_array_internal(lua_State *L, int i, int pos_i, void *buf_i, 
        qloptNdarray *a, const char *msg_head);


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
        return qlopt_read_array_internal(L, 0, pos, t.arr, base, msg_head); 
    }
    case QLOPT_SUBSET       :
        return luaL_error(L, "%s unsupported opt.type %s", msg_head, 
                          qlopt_type_name(argtype));
    default : return luaL_error(L, "%s unsupported arg.type %d", msg_head, argtype);
    }
}

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
        lua_pushinteger(L, 0);
        lua_rawget(L, i_pos);
        i_pos = lua_gettop(L);
    } while (i < a->ndim);

    /* clean up exploratory pushes */
    lua_pop(L, i);

    /* allocate array if necessary */
    int nmemb = prod_int_(a->ndim, a->dim);
    int bufsize = nmemb * a->argsize;
    if (NULL == a->buf) {
        a->buf = malloc(bufsize);
        if (NULL == a->buf)
            return luaL_error(L, "%s not enough memory (try alloc %llu)", 
                              msg_head, (size_t)(bufsize));
        a->bufsize = bufsize;
    } else {
        if (a->bufsize < bufsize)
            return luaL_error(L, "%s allocated space is not sufficient (%d < %d)", 
                    msg_head, a->bufsize, bufsize);
    }
    return 0;
}

static int 
qlopt_read_array_internal(lua_State *L, int i, int pos_i, void *buf_i, 
        qloptNdarray *a, const char *msg_head)
{
    int status = 0;
    char strbuf0[128], strbuf[1024];
    if (NULL == msg_head) {
        snprintf(strbuf0, sizeof(strbuf0), "(arg %d)", pos_i);
        msg_head = strbuf0;
    }
    if (i < a->ndim) {
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
            if (0 != (status = qlopt_read_array_internal(L, i + 1, lua_gettop(L), 
                        buf_i + k * stride_memb * a->argsize, a, strbuf))) {
                lua_pop(L, 1);
                return status;
            }
            lua_pop(L, 1);
        }
        return status;
    } else {
        qloptData zero_off = {.off = 0};
        return qlopt_read_arg_internal(L, pos_i, 0, a->argtype, zero_off, buf_i, strbuf);
    }
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
    return qlopt_read_array_internal(L, 0, pos, a->buf, a, strbuf);
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
            snprintf(strbuf, sizeof(strbuf), "(arg %d)", pos_next);
            pos = pos_next;
            pos_next ++;
        } else {
            if (elem->pos <= 0)
                return luaL_error(L, "(elem %d) bad pos=%d", i_elem, elem->pos);
            pos = elem->pos;
            pos_next = pos + 1;
        }
        int status = qlopt_read_arg_internal(L, pos, elem->flags, 
                elem->argtype, elem->val, NULL, strbuf);
        if (status)
            return luaL_error(L, "%s : cannot parse argument", strbuf);

        elem++;
        i_elem++;
    }
    return 0;

}

#if 0
static int 
test_compile_struc_init()
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
    int l_nev = -1,
        l_ncv = -1,
        l_maxiter = -1;
    double l_tol = -1.;

    int l_n = -1;
    double l_a  = DBL_MAX,
           l_b  = DBL_MAX,
           l_x0 = DBL_MAX;
    qloptElem_s parse_tab_cheb_accel[] = {
        { 1, NULL, QLOPT_INT,   0,   { .i = &l_n } },
        { 2, NULL, QLOPT_REAL,  0,   { .r = &l_a } },
        { 3, NULL, QLOPT_REAL,  0,   { .r = &l_a } },
        { 4, NULL, QLOPT_REAL,  QLOPT_OPT,   { .r = &l_x0 } },
        { QLOPT_NONE }
    };

    /* variables to fill, with default values ; 
       while default values are not necessary for mandatory options, 
       the table itself will be mandatory, so better initialize */
    int eigcg_vmax  = -1, 
        eigcg_nev   = -1, 
        eigcg_umax  = -1;
    double eigcg_tol = -1.;
    qloptElem_s parse_tab_eigcg[] = {
        { 1, NULL, QLOPT_INT,   0, { .i = &eigcg_vmax } },
        { 2, NULL, QLOPT_INT,   0, { .i = &eigcg_nev  } },
        { 3, NULL, QLOPT_REAL,  0, { .r = &eigcg_tol  } },
        { 4, NULL, QLOPT_INT,   0, { .i = &eigcg_umax } },
        { QLOPT_END }
    };
    
    int coeffs_dim[]    = { -1, 2 },
        coeffs_maxdim[] = { 10, 3 }; /* '3' is ignored because the coeffs_dim[1] is already known */
    double *coeffs = NULL;
    qloptNdarray_s parse_arr_coeffs = { 
        QLOPT_REAL, 2, coeffs_dim, coeffs_maxdim, &coeffs };

    int mom_dim[] = {10, 4};
    int mom[10][4],  
        *p_mom = &mom[0][0];
    qloptNdarray_s parse_arr_mom   = { 
            QLOPT_INT, 2, mom_dim, NULL, &p_mom };

    int xy_dim[] = {2,3,4};
    int xy[2][3][4], 
        *p_xy = (int*)&xy;

    const char *which = NULL;
    const char *arpack_logfile = NULL;
    int inplace = 0;

    int abc_dim[]    = { -1, -1, -1 },
        abc_maxdim[] = { 10, 10, 10 };
    complex double *abc = NULL;


    qloptElem_s parse_tab_opt[] = {
        { QLOPT_KEY, "eigcg",       QLOPT_TABLE,  0, { .t = parse_tab_eigcg } },
        { QLOPT_KEY, "cheb_accel",  QLOPT_TABLE,  0, { .t = parse_tab_cheb_accel } },
        { QLOPT_KEY, "coeff_list",  QLOPT_ARRAY,  0, { .a = parse_arr_coeffs } },
        { QLOPT_KEY, "mom_list",    QLOPT_ARRAY,  0, { .a = parse_arr_mom } },
        { QLOPT_KEY, "xy",          QLOPT_ARRAY,  0,    
                    { .a = { QLOPT_INT, 3, xy_dim, NULL, &p_xy } } },
        { QLOPT_KEY, "which",       QLOPT_STRING, 0, { .s = &which } },
        { QLOPT_KEY, "arpack_logfile", QLOPT_STRING, QLOPT_OPT,
                    { .s = &arpack_logfile } },
        { QLOPT_KEY, "inplace",     QLOPT_BOOL, QLOPT_OPT,
                    { .i = &inplace } },
        { QLOPT_END }
    } ;
    qloptElem_s parse_tab_arg[] = {
        { 1, NULL, QLOPT_INT,   0, { .i = &l_nev } },
        { 2, NULL, QLOPT_INT,   0, { .i = &l_ncv } },
        { 3, NULL, QLOPT_INT,   0, { .i = &l_maxiter } },
        { 4, NULL, QLOPT_REAL,  0, { .r = &l_tol } },
        { 5, NULL, QLOPT_TABLE, QLOPT_OPT, { .t = parse_tab_opt } },
        { 6, NULL, QLOPT_ARRAY, QLOPT_OPT, 
                    { .a = { QLOPT_COMPLEX, 3, abc_dim, abc_maxdim, &abc } } },
        { QLOPT_END }
    };

    int status;
  
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

    return 0;
}
#endif
