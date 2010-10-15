#include "qlua.h"                                                    /* DEPS */
#include "modules.h"
#include "fix.h"                                                     /* DEPS */
#include "qcomplex.h"                                                /* DEPS */
#include "seqrandom.h"                                               /* DEPS */
#include "seqcolvec.h"                                               /* DEPS */
#include "seqcolmat.h"                                               /* DEPS */
#include "seqdirferm.h"                                              /* DEPS */
#include "seqdirprop.h"                                              /* DEPS */
#include "qvector.h"                                                 /* DEPS */
#include "qxml.h"                                                    /* DEPS */
#include "lattice.h"                                                 /* DEPS */
#include "latsubset.h"                                               /* DEPS */
#include "latmulti.h"                                                /* DEPS */
#include "latint.h"                                                  /* DEPS */
#include "latreal.h"                                                 /* DEPS */
#include "latrandom.h"                                               /* DEPS */
#include "latcomplex.h"                                              /* DEPS */
#include "latcolvec.h"                                               /* DEPS */
#include "latcolmat.h"                                               /* DEPS */
#include "latdirferm.h"                                              /* DEPS */
#include "latdirprop.h"                                              /* DEPS */
#include "qgamma.h"                                                  /* DEPS */
#include "nersc_io.h"                                                /* DEPS */
#include "qdpc_io.h"                                                 /* DEPS */
#include "ddpairs_io.h"                                              /* DEPS */
#include "qdpcc_io.h"                                                /* DEPS */
#include "qmp.h"
#ifdef HAS_AFF
#include "lhpc-aff.h"
#include "aff_io.h"                                                  /* DEPS */
#endif
#ifdef HAS_CLOVER
#include "qclover.h"                                                 /* DEPS */
#endif
#ifdef HAS_MDWF
#include "qmdwf.h"                                                   /* DEPS */
#endif
#ifdef HAS_EXTRAS
#include "extras.h"                                                  /* DEPS */
#endif
#ifdef HAS_GSL
#include "qmatrix.h"                                                 /* DEPS */
#endif

#include <string.h>
#include <stdarg.h>
#include <qmp.h>
#include <assert.h>
#include <libgen.h>

/* ZZZ include other package headers here */

const static char *a_type_key = "a-type";
const static char *lattice_key = "lattice";
const char *progname = "qlua";
const char *qcdlib = "qcd";
int qlua_master_node = -1;

static struct {
    char *name;
    char *value;
} versions[] = {
    {"qlua",  "QLUA version 0.20.00-rc1 $Id$"},
    {"lua",    LUA_VERSION },
    {"qdp",    QDP_VERSION },
#ifdef HAS_AFF
    {"aff",    LHPC_AFF_VERSION },
#endif
#ifdef HAS_CLOVER
    {"clover", CLOVER_VERSION },
#endif
#ifdef HAS_MDWF
    {"mdwf", MDWF_VERSION },
#endif
#ifdef HAS_EXTRAS
    {"extras", "included" },
#endif
#ifdef HAS_GSL
    {"gsl",     GSL_VERSION },
#endif
#ifdef HAS_CBLAS
    {"cblas",     CBLAS_VERSION },
#endif
    {"colors",   (
#if USE_Nc2
                   " 2"
#endif
#if USE_Nc3
                   " 3"
#endif
#if USE_NcN
                   " N"
#endif
                 ) + 1 },
    {NULL,     NULL}
};

/* reporting */
void
message(const char *fmt, ...)
{
    if ((qlua_master_node < 0) || (QDP_this_node == qlua_master_node)) {
        va_list va;

        va_start(va, fmt);
        vfprintf(stderr, fmt, va);
        va_end(va);
    }
}

void
report(lua_State *L, const char *fname, int status)
{
    if ((qlua_master_node < 0) || (QDP_this_node == qlua_master_node)) {
        if (status && !lua_isnil(L, -1)) {
            const char *msg = lua_tostring(L, -1);
            if (msg == NULL) msg = "(error object is not a string)";
            message("%s ERROR:: %s\n", progname, msg);
            lua_pop(L, 1);
        }
    }
}

void
XMP_dist_int_array(int src_node, int count, int *data)
{
    int i;

    if (src_node != QDP_this_node)
        memset(data, 0, count * sizeof (int));
    for (i = 0; i < count; i++)
        QMP_sum_int(&data[i]);
}

void
XMP_dist_double_array(int src_node, int count, double *data)
{
    if (src_node != QDP_this_node)
        memset(data, 0, count * sizeof (double));
    QMP_sum_double_array(data, count);
}

/* memory allocation */
void *
qlua_malloc(lua_State *L, int size)
{
    void *p = malloc(size);
    if (p == 0) {
        lua_gc(L, LUA_GCCOLLECT, 0);
        p = malloc(size);
        if (p == 0)
            luaL_error(L, "not enough memory");
    }
    return p;
}

void
qlua_free(lua_State *L, void *ptr)
{
    if (ptr)
        free(ptr);
}

static void
qlua_fillmeta(lua_State *L, const luaL_Reg *table, QLUA_Type t_id)
{
    int i;

    lua_pushvalue(L, -1);
    lua_setfield(L, -2, "__index");
    
    if (t_id != qNoType) {
        lua_pushstring(L, a_type_key);
        lua_pushnumber(L, t_id);
        lua_settable(L, -3);
    }
    for (i = 0; table[i].func; i++) {
        lua_pushstring(L, table[i].name);
        lua_pushcfunction(L, table[i].func);
        lua_settable(L, -3);
    }
}

void
qlua_metatable(lua_State *L, const char *name, const luaL_Reg *table,
               QLUA_Type t_id)
{
    luaL_newmetatable(L, name);
    luaL_getmetatable(L, name);
    qlua_fillmeta(L, table, t_id);
    lua_pop(L, 1);
}

void
qlua_selftable(lua_State *L, const luaL_Reg *table, QLUA_Type t_id)
{
    lua_createtable(L, 0, 0);
    qlua_fillmeta(L, table, t_id);
}

void
qlua_latticetable(lua_State *L, const luaL_Reg *table, QLUA_Type t_id, int lidx)
{
    lua_createtable(L, 0, 0);
    qlua_fillmeta(L, table, t_id);
    lua_pushstring(L, lattice_key);
    lua_pushvalue(L, lidx);
    lua_settable(L, -3);
}

void
qlua_createLatticeTable(lua_State *L,
                        int Sidx,
                        const struct luaL_Reg *ft,
                        QLUA_Type t_id,
                        const char *name)
{
    if (luaL_getmetafield(L, Sidx, name) != 0)
        goto end;
    lua_getmetatable(L, Sidx);
    qlua_latticetable(L, ft, t_id, Sidx);
    lua_setfield(L, -2, name);
    lua_pop(L, 1);
    luaL_getmetafield(L, Sidx, name);
end:
    return;
}

int
qlua_badconstr(lua_State *L, const char *name)
{
    return luaL_error(L, "bad %s constructor", name);
}

int
qlua_badindex(lua_State *L, const char *type)
{
    return luaL_error(L, "bad index for %s", type);
}

int
qlua_lookup(lua_State *L, int idx, const char *table)
{
    const char *key = lua_tostring(L, idx);

    luaL_getmetatable(L, table);
    lua_getfield(L, -1, key);
    if (lua_isnil(L, -1))
        return luaL_error(L, "bad index `%s' for ::%s", key, table);
    return 1;
}

int
qlua_selflookup(lua_State *L, int idx, const char *key)
{
    if (luaL_getmetafield(L, idx, key) == 0)
        return luaL_error(L, "bad index `%s'", key);

    return 1;
}

int
qlua_getlattice(lua_State *L, int idx)
{
    if (luaL_getmetafield(L, idx, lattice_key) == 0)
        return luaL_error(L, "lattice object expected");
    return 1;
}

int
qlua_checkint(lua_State *L, int idx, const char *fmt, ...)
{
    lua_Integer d = lua_tointeger(L, idx);
    
    if (d == 0 && !lua_isnumber(L, idx)) {
        char buf[72];
        va_list va;
        va_start(va, fmt);
        vsnprintf(buf, sizeof (buf) - 1, fmt, va);
        va_end(va);
        luaL_error(L, "expected int, %s", buf);
    }
    return (int)d;
}

const char *
qlua_checkstring(lua_State *L, int idx, const char *fmt, ...)
{
    const char *d = lua_tostring(L, idx);
    
    if (d == 0) {
        char buf[72];
        va_list va;
        va_start(va, fmt);
        vsnprintf(buf, sizeof (buf) - 1, fmt, va);
        va_end(va);
        luaL_error(L, "expected string, %s", buf);
    }
    return d;
}

void
qlua_checktable(lua_State *L, int idx, const char *fmt, ...)
{
    if (lua_type(L, idx) != LUA_TTABLE) {
        char buf[72];
        va_list va;
        va_start(va, fmt);
        vsnprintf(buf, sizeof (buf) - 1, fmt, va);
        va_end(va);
        luaL_error(L, "expected table, %s", buf);
    }
}

int
qlua_index(lua_State *L, int n, const char *name, int max_value)
{
    int v = -1;

    luaL_checktype(L, n, LUA_TTABLE);
    lua_getfield(L, n, name);
    if (lua_isnumber(L, -1)) {
        v = qlua_checkint(L, -1, "integer in [0,%d) expected", max_value);
        if ((v < 0) || (v >= max_value))
            v = -1;
    }
    lua_pop(L, 1);
    
    return v;
}

void
qlua_checkindex2(lua_State *L, int n, const char *name, int *sl, int *sr)
{
    luaL_checktype(L, n, LUA_TTABLE);
    lua_pushnumber(L, 1);
    lua_gettable(L, n);
    *sl = qlua_checkint(L, -1, "integer expected for index0 in %s", name);
    lua_pop(L, 1);

    lua_pushnumber(L, 2);
    lua_gettable(L, n);
    *sr = qlua_checkint(L, -1, "integer expected for index1 in %s", name);
    lua_pop(L, 1);
}

int
qlua_checkindex(lua_State *L, int n, const char *name, int max_value)
{
    int v = qlua_index(L, n, name, max_value);

    if (v == -1)
        luaL_error(L, "bad index");

    return v;
}

int
qlua_diracindex(lua_State *L, int n)
{
    return qlua_index(L, n, "d", QDP_Ns);
}

int
qlua_checkdiracindex(lua_State *L, int n)
{
    return qlua_checkindex(L, n, "d", QDP_Ns);
}

int
qlua_colorindex(lua_State *L, int n, int nc)
{
    return qlua_index(L, n, "c", nc);
}

int
qlua_checkcolorindex(lua_State *L, int n, int nc)
{
    return qlua_checkindex(L, n, "c", nc);
}

int
qlua_leftindex(lua_State *L, int n, int nc)
{
    return qlua_index(L, n, "a", nc);
}

int
qlua_checkleftindex(lua_State *L, int n, int nc)
{
    return qlua_checkindex(L, n, "a", nc);
}

int
qlua_rightindex(lua_State *L, int n, int nc)
{
    return qlua_index(L, n, "b", nc);
}

int
qlua_checkrightindex(lua_State *L, int n, int nc)
{
    return qlua_checkindex(L, n, "b", nc);
}

int
qlua_gammaindex(lua_State *L, int n)
{
    int d = qlua_index(L, n, "mu", 6);

    if (d == 4)
        return -1;
    return d;
}

int
qlua_checkgammaindex(lua_State *L, int n)
{
    int d = qlua_gammaindex(L, n);

    if (d == -1)
        return luaL_error(L, "bad index");

    return d;
}

int
qlua_gammabinary(lua_State *L, int n)
{
    return qlua_index(L, n, "n", 16);
}

int
qlua_checkgammabinary(lua_State *L, int n)
{
    return qlua_checkindex(L, n, "n", 16);
}

void *
qlua_checkLatticeType(lua_State *L, int idx, QLUA_Type t_id,
                      const char *name)
{
    void *p = lua_touserdata(L, idx);
    int p_id;

    if (p == NULL)
        goto bad_value;
    if (luaL_getmetafield(L, idx, a_type_key) == 0)
        goto bad_value;

    p_id = luaL_checkint(L, -1);
    lua_pop(L, 1);

    if (p_id != t_id)
        goto bad_value;

    return p;
bad_value:
    luaL_error(L, "expecting %s", name);
    return NULL;
}

QLUA_Type
qlua_qtype(lua_State *L, int idx)
{
    luaL_checkany(L, idx);
    switch (lua_type(L, idx)) {
    case LUA_TNUMBER:
        return qReal;
    case LUA_TSTRING:
        return qString;
    case LUA_TTABLE:
        return qTable;
    case LUA_TUSERDATA: {
        QLUA_Type tv = qOther;

        lua_getmetatable(L, idx);
        lua_pushstring(L, a_type_key);
        lua_rawget(L, -2);
        if (lua_isnumber(L, -1))
            tv = luaL_checkint(L, -1);
        lua_pop(L, 2);
        return tv;
    }
    default:
        return qOther;
    }
}

QLUA_Type
qlua_atype(lua_State *L, int idx)
{
    QLUA_Type tv = qlua_qtype(L, idx);

    if (tv > qOther)
        tv = qOther;

    return tv;
}

QLUA_Ztype
qlua_ztype(lua_State *L, int idx)
{
    switch (qlua_qtype(L, idx)) {
    case qReal:     return zReal;
    case qLatInt:   return zLatInt;
    case qLatReal:  return zLatReal;
    default:        return zNoType;
    }
}
/* generic operations dispatchers */
#define Op2Idx(a,b) ((a)*qArithTypeCount + (b))
#define zOp2Idx(a,b) ((a)*zArithTypeCount + (b))

q_op qlua_add_table[qArithTypeCount * qArithTypeCount];
q_op qlua_sub_table[qArithTypeCount * qArithTypeCount];
q_op qlua_mul_table[qArithTypeCount * qArithTypeCount];
q_op qlua_div_table[qArithTypeCount * qArithTypeCount];
q_op qlua_mod_table[qArithTypeCount * qArithTypeCount];
q_op qlua_min_table[zArithTypeCount * zArithTypeCount];
q_op qlua_max_table[zArithTypeCount * zArithTypeCount];
q_op qlua_eq_table[zArithTypeCount * zArithTypeCount];
q_op qlua_ne_table[zArithTypeCount * zArithTypeCount];
q_op qlua_lt_table[zArithTypeCount * zArithTypeCount];
q_op qlua_le_table[zArithTypeCount * zArithTypeCount];
q_op qlua_gt_table[zArithTypeCount * zArithTypeCount];
q_op qlua_ge_table[zArithTypeCount * zArithTypeCount];

void
qlua_reg_op2(const QLUA_Op2 *ops)
{
    int i;

    for (i = 0; ops[i].table; i++) {
        assert(ops[i].ta < qOther);
        assert(ops[i].tb < qOther);
        ops[i].table[Op2Idx(ops[i].ta, ops[i].tb)] = ops[i].op;
    }
}

void
qlua_reg_zop2(const QLUA_ZOp2 *ops)
{
    int i;

    for (i = 0; ops[i].table; i++) {
        ops[i].table[zOp2Idx(ops[i].ta, ops[i].tb)] = ops[i].op;
    }
}

int 
qlua_add(lua_State *L)
{
    int idx = Op2Idx(qlua_atype(L, 1), qlua_atype(L, 2));

    if (qlua_add_table[idx])
        return qlua_add_table[idx](L);
    else
        return luaL_error(L, "bad argument for addition");
}

int 
qlua_sub(lua_State *L)
{
    int idx = Op2Idx(qlua_atype(L, 1), qlua_atype(L, 2));

    if (qlua_sub_table[idx])
        return qlua_sub_table[idx](L);
    else
        return luaL_error(L, "bad argument for subtraction");
}

int 
qlua_mul(lua_State *L)
{
    int idx = Op2Idx(qlua_atype(L, 1), qlua_atype(L, 2));
    
    if (qlua_mul_table[idx])
        return qlua_mul_table[idx](L);
    else
        return luaL_error(L, "bad argument for multiplication");
}

int 
qlua_div(lua_State *L)
{
    int idx = Op2Idx(qlua_atype(L, 1), qlua_atype(L, 2));

    if (qlua_div_table[idx])
        return qlua_div_table[idx](L);
    else
        return luaL_error(L, "bad argument for division");
}

int 
qlua_mod(lua_State *L)
{
    int idx = Op2Idx(qlua_atype(L, 1), qlua_atype(L, 2));

    if (qlua_mod_table[idx])
        return qlua_mod_table[idx](L);
    else
        return luaL_error(L, "bad argument for modulo");
}

#define Op1Idx(a)   (a)
static q_op qt_dot[(qOther + 1)];

void
qlua_reg_dot(QLUA_Type ta, q_op op)
{
    assert(ta < qOther);

    qt_dot[Op1Idx(ta)] = op;
}

static int 
q_dot(lua_State *L)
{
    int idx = Op1Idx(qlua_atype(L, 1));
    if (qt_dot[idx]) 
        return qt_dot[idx](L);
    else
        return luaL_error(L, "bad argument for inner dot");
}

static int 
q_min(lua_State *L)
{
    int idx = zOp2Idx(qlua_ztype(L, 1), qlua_ztype(L, 2));
    
    if (qlua_min_table[idx])
        return qlua_min_table[idx](L);
    else
        return luaL_error(L, "bad argument for <?");
}

static int 
q_max(lua_State *L)
{
    int idx = zOp2Idx(qlua_ztype(L, 1), qlua_ztype(L, 2));
    
    if (qlua_max_table[idx])
        return qlua_max_table[idx](L);
    else
        return luaL_error(L, "bad argument for >?");
}

static int 
q_eq(lua_State *L)
{
    int idx = zOp2Idx(qlua_ztype(L, 1), qlua_ztype(L, 2));
    
    if (qlua_eq_table[idx])
        return qlua_eq_table[idx](L);
    else
        return luaL_error(L, "bad argument for ==");
}

static int 
q_ne(lua_State *L)
{
    int idx = zOp2Idx(qlua_ztype(L, 1), qlua_ztype(L, 2));
    
    if (qlua_ne_table[idx])
        return qlua_ne_table[idx](L);
    else
        return luaL_error(L, "bad argument for !=");
}

static int 
q_lt(lua_State *L)
{
    int idx = zOp2Idx(qlua_ztype(L, 1), qlua_ztype(L, 2));
    
    if (qlua_lt_table[idx])
        return qlua_lt_table[idx](L);
    else
        return luaL_error(L, "bad argument for <");
}

static int 
q_le(lua_State *L)
{
    int idx = zOp2Idx(qlua_ztype(L, 1), qlua_ztype(L, 2));
    
    if (qlua_le_table[idx])
        return qlua_le_table[idx](L);
    else
        return luaL_error(L, "bad argument for <=");
}

static int 
q_gt(lua_State *L)
{
    int idx = zOp2Idx(qlua_ztype(L, 1), qlua_ztype(L, 2));
    
    if (qlua_gt_table[idx])
        return qlua_gt_table[idx](L);
    else
        return luaL_error(L, "bad argument for >");
}

static int 
q_ge(lua_State *L)
{
    int idx = zOp2Idx(qlua_ztype(L, 1), qlua_ztype(L, 2));
    
    if (qlua_ge_table[idx])
        return qlua_ge_table[idx](L);
    else
        return luaL_error(L, "bad argument for >=");
}

static struct luaL_Reg fQCD[] = {
    { "dot",    q_dot }, /* local inner dot */
    { "min",    q_min }, /* local min */
    { "max",    q_max }, /* local min */
    { "ne",     q_ne  }, /* local  != */
    { "eq",     q_eq  }, /* local  == */
    { "lt",     q_lt  }, /* local  <  */
    { "le",     q_le  }, /* local  <= */
    { "gt",     q_gt  }, /* local  >  */
    { "ge",     q_ge  }, /* local  >= */
    { NULL,     NULL}
};

/* gemeric error report on illegal write operations */
int
qlua_nowrite(lua_State *L)
{
    return luaL_error(L, "assignment is not permitted");
}

/* environment setup */
static void
qlua_openlibs (lua_State *L) {
    static const luaL_Reg lualibs[] = {
        {"", luaopen_base},
        {LUA_LOADLIBNAME, luaopen_package},
        {LUA_TABLIBNAME, luaopen_table},
        {LUA_OSLIBNAME, luaopen_os},
        {LUA_STRLIBNAME, luaopen_string},
        {LUA_MATHLIBNAME, luaopen_math},
        {LUA_DBLIBNAME, luaopen_debug},
        {NULL, NULL}
    };
    
    const luaL_Reg *lib = lualibs;
    for (; lib->func; lib++) {
        lua_pushcfunction(L, lib->func);
        lua_pushstring(L, lib->name);
        lua_call(L, 1, 0);
    }
}

void
qlua_init(lua_State *L, int argc, char *argv[])
{
    static const lua_CFunction qcd_inits[] = {
        init_qlua_io,
        init_complex,
        init_seqrandom,
        init_seqcolvec,
        init_seqcolmat,
        init_seqdirferm,
        init_seqdirprop,
        init_xml,
        init_vector,
#ifdef HAS_GSL
        init_matrix,
#endif
        init_lattice,
        init_latint,
        init_latreal,
        init_latrandom,
        init_latcomplex,
        init_latsubset,
        init_latmulti,
        init_latcolvec,
        init_latcolmat,
        init_latdirferm,
        init_latdirprop,
        init_gamma,
#ifdef HAS_AFF
        init_aff_io,
#endif
        init_nersc_io,
        init_qdpc_io,
        init_ddpairs_io,
        init_qdpcc_io,
#ifdef HAS_CLOVER
        init_clover,
#endif
#ifdef HAS_MDWF
        init_mdwf,
#endif
#ifdef HAS_EXTRAS
        init_extras,
#endif
        /* ZZZ add other packages here */
        NULL };

    int i;

    lua_gc(L, LUA_GCSTOP, 0);  /* stop collector during initialization */
    qlua_openlibs(L);  /* open libraries */

    // set exename, exepath and exerealpath
    char ra[PATH_MAX];
    realpath(argv[0], ra);
    char *bn = basename(argv[0]);
    char *dn = dirname(argv[0]);
    char *rdn = dirname(ra);
    lua_pushstring(L, bn);
    lua_setglobal(L, "exename");
    lua_pushstring(L, dn);
    lua_setglobal(L, "exepath");
    lua_pushstring(L, rdn);
    lua_setglobal(L, "exerealpath");

    luaL_register(L, qcdlib, fQCD);
    for (i = 0; qcd_inits[i]; i++) {
        lua_pushcfunction(L, qcd_inits[i]);
        lua_call(L, 0, 0);
    }
    lua_getglobal(L, qcdlib);
    lua_pushnumber(L, QDP_Ns);
    lua_setfield(L, -2, "Ns");
    lua_newtable(L);
    for (i = 0; versions[i].name; i++) {
        lua_pushstring(L, versions[i].value);
        lua_setfield(L, -2, versions[i].name);
    }
    lua_setfield(L, -2, "version");
    lua_pop(L, 1);

    /* collect all script names */
    {
        int i;
        lua_createtable(L, argc - 1, 0);
        for (i = 1; i < argc; i++) {
            lua_pushstring(L, argv[i]);
            lua_rawseti(L, -2, i);
        }
        lua_setglobal(L, "scripts");
    }
    lua_gc(L, LUA_GCRESTART, 0);
}

/* cleanup (mostly housekeeping to make various tools happy */
void
qlua_fini(lua_State *L)
{
    const static lua_CFunction qcd_finis[] = {
        /* ZZZ add other packages here */
#ifdef HAS_EXTRAS
        fini_extras,
#endif
#ifdef HAS_MDWF
        fini_mdwf,
#endif
#ifdef HAS_CLOVER
        fini_clover,
#endif
        fini_qdpcc_io,
        fini_ddpairs_io,
        fini_qdpc_io,
        fini_nersc_io,
#ifdef HAS_AFF
        fini_aff_io,
#endif
        fini_gamma,
        fini_latdirprop,
        fini_latdirferm,
        fini_latcolmat,
        fini_latcolvec,
        fini_latmulti,
        fini_latsubset,
        fini_latcomplex,
        fini_latrandom,
        fini_latreal,
        fini_latint,
        fini_lattice,
#ifdef HAS_GSL
        fini_matrix,
#endif
        fini_vector,
        fini_xml,
        fini_seqdirprop,
        fini_seqdirferm,
        fini_seqcolmat,
        fini_seqcolvec,
        fini_seqrandom,
        fini_complex,
        fini_qlua_io,
        NULL };
    int i;

    for (i = 0; qcd_finis[i]; i++) {
        lua_pushcfunction(L, qcd_finis[i]);
        lua_call(L, 0, 0);
    }
}

/* the driver */
int
main(int argc, char *argv[])
{
    int status = 1;
    int i;
    lua_State *L = NULL;

    if (QDP_initialize(&argc, &argv)) {
        fprintf(stderr, "QDP initialization failed\n");
        return 1;
    }
    QDP_profcontrol(0);
    double node = QDP_this_node;
    QMP_min_double(&node);
    qlua_master_node = node;

    L = lua_open();
    if (L == NULL) {
        message("can not create Lua state");
        goto end;
    }
    qlua_init(L, argc, argv);  /* open libraries */

    if (argc < 2) {
        message("QLUA component versions:\n");
        for (i = 0; versions[i].name; i++)
            message(" %10s: %s\n", versions[i].name, versions[i].value);
    } else {

      for (i = 1; i < argc; i++) {
	if(strcmp(argv[i],"-e")==0) { // process command
	  const char *chunk = argv[i] + 2;
	  if (*chunk == '\0') {
	    if(++i>=argc) {
	      message("missing argument to -e");
	      goto end;
	    }
	    chunk = argv[i];
	  }
	  lua_assert(chunk != NULL);
	  status = luaL_dostring(L, chunk);
	  report(L, "=(command line)", status);
	  if (status) {
	    QDP_abort(1);
	    break;
	  }
	} else { // process file
	  status = luaL_dofile(L, argv[i]);
	  report(L, argv[i], status);
	  if (status) {
	    QDP_abort(1);
	    break;
	  }
	}
      }

    }
    qlua_fini(L);
    lua_close(L);
end:
    QDP_finalize();
    return status;
}
