#include <qlua.h>                                                    /* DEPS */
#include <modules.h>
#include <fix.h>                                                     /* DEPS */
#include <qcomplex.h>                                                /* DEPS */
#include <qgamma.h>                                                  /* DEPS */
#include <qvector.h>                                                 /* DEPS */
#include <lattice.h>                                                 /* DEPS */
#include <latint.h>                                                  /* DEPS */
#include <latrandom.h>                                               /* DEPS */
#include <latreal.h>                                                 /* DEPS */
#include <latcomplex.h>                                              /* DEPS */
#include <latcolvec.h>                                               /* DEPS */
#include <latcolmat.h>                                               /* DEPS */
#include <latdirferm.h>                                              /* DEPS */
#include <latdirprop.h>                                              /* DEPS */
#include <latsubset.h>                                               /* DEPS */
#include <latmulti.h>                                                /* DEPS */
#include <qdpc_io.h>                                                 /* DEPS */
#include <qdpcc_io.h>                                                /* DEPS */
#include <ddpairs_io.h>                                              /* DEPS */
#include <nersc_io.h>                                                /* DEPS */
#include <qxml.h>                                                    /* DEPS */
#ifdef HAS_AFF
#include <lhpc-aff.h>
#include <aff_io.h>                                                  /* DEPS */
#endif
#ifdef HAS_CLOVER
#include <qclover.h>                                                 /* DEPS */
#endif
#ifdef HAS_MDWF
#include <qmdwf.h>                                                   /* DEPS */
#endif

#include <string.h>
#include <qmp.h>

/* ZZZ include other package headers here */

const char *progname = "qlua";
const char *qcdlib = "qcd";
int qlua_primary_node = 1;

static struct {
    char *name;
    char *value;
} versions[] = {
    {"qlua",  "QLUA version 0.9.5+XXX $Id$"},
    {"lua",    LUA_VERSION },
    {"qdp",    QDP_VERSION },
#ifdef HAS_AFF
    {"aff",    AFF_VERSION },
#endif
#ifdef HAS_CLOVER
    {"clover", CLOVER_VERSION },
#endif
#ifdef HAS_MDWF
    {"mdwf",   MDWF_VERSION },
#endif
    {NULL,     NULL}
};

/* reporting */
void
message(const char *fmt, ...)
{
    if (qlua_primary_node) {
        va_list va;

        va_start(va, fmt);
        vfprintf(stderr, fmt, va);
        va_end(va);
    }
}

void
report(lua_State *L, const char *fname, int status)
{
    if (qlua_primary_node) {
        if (status && !lua_isnil(L, -1)) {
            const char *msg = lua_tostring(L, -1);
            if (msg == NULL) msg = "(error object is not a string)";
            message("%s ERROR:: %s\n", progname, msg);
            lua_pop(L, 1);
        }
    }
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

void
qlua_metatable(lua_State *L, const char *name, const luaL_Reg *table)
{
    int i;
    int base = lua_gettop(L);

    luaL_newmetatable(L, name);
    luaL_getmetatable(L, name);
    lua_pushstring(L, "__index");
    luaL_getmetatable(L, name);
    lua_settable(L, -3);
    for (i = 0; table[i].func; i++) {
        lua_pushstring(L, table[i].name);
        lua_pushcfunction(L, table[i].func);
        lua_settable(L, -3);
    }
    lua_settop(L, base);
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
qlua_index(lua_State *L, int n, const char *name, int max_value)
{
    int v = -1;

    luaL_checktype(L, n, LUA_TTABLE);
    lua_getfield(L, n, name);
    if (lua_isnumber(L, -1)) {
        v = luaL_checkint(L, -1);
        if ((v < 0) || (v >= max_value))
            v = -1;
    }
    lua_pop(L, 1);
    
    return v;
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
qlua_colorindex(lua_State *L, int n)
{
    return qlua_index(L, n, "c", QDP_Nc);
}

int
qlua_checkcolorindex(lua_State *L, int n)
{
    return qlua_checkindex(L, n, "c", QDP_Nc);
}

int
qlua_leftindex(lua_State *L, int n)
{
    return qlua_index(L, n, "a", QDP_Nc);
}

int
qlua_checkleftindex(lua_State *L, int n)
{
    return qlua_checkindex(L, n, "a", QDP_Nc);
}

int
qlua_rightindex(lua_State *L, int n)
{
    return qlua_index(L, n, "b", QDP_Nc);
}

int
qlua_checkrightindex(lua_State *L, int n)
{
    return qlua_checkindex(L, n, "b", QDP_Nc);
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

static int
qlua_type(lua_State *L, int idx, const char *mt)
{
    int base = lua_gettop(L);
    int v;

    lua_getmetatable(L, idx);
    luaL_getmetatable(L, mt);
    v = lua_equal(L, -1, -2);
    lua_settop(L, base);

    return v;
}

int
qlua_gettype(lua_State *L, int idx)
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
        static const struct {
            const char *name;
            int ty;
        } t[] = {
            { mtnComplex,       qComplex },
            { mtnGamma,         qGamma },
            { mtnVecInt,        qVecInt },
            { mtnVecReal,       qVecReal },
            { mtnVecComplex,    qVecComplex },
            { mtnLatInt,        qLatInt },
            { mtnLatReal,       qLatReal },
            { mtnLatRandom,     qLatRandom },
            { mtnLatComplex,    qLatComplex },
            { mtnLatColVec,     qLatColVec },
            { mtnLatColMat,     qLatColMat },
            { mtnLatDirFerm,    qLatDirFerm },
            { mtnLatDirProp,    qLatDirProp },
            /* ZZZ other types */
            { NULL,              qOther }
        };
        int i;

        for (i = 0; t[i].name; i++)
            if (qlua_type(L, idx, t[i].name)) return t[i].ty;
        return qOther;
    }
    default:
        return qOther;
    }
}

/* generic operations dispatchers */
#define Op1Idx(a)   (a)
#define Op2Idx(a,b) ((a)*(qOther + 1) + (b))

static q_op qt_add[(qOther + 1) * (qOther + 1)];

void
qlua_reg_add(int ta, int tb, q_op op)
{
    qt_add[Op2Idx(ta, tb)] = op;
}

int 
qlua_add(lua_State *L)
{
    q_op op = qt_add[Op2Idx(qlua_gettype(L, 1), qlua_gettype(L, 2))];

    if (op)
        return op(L);
    else
        return luaL_error(L, "bad argument for addition");
}

static q_op qt_sub[(qOther + 1) * (qOther + 1)];

void
qlua_reg_sub(int ta, int tb, q_op op)
{
    qt_sub[Op2Idx(ta, tb)] = op;
}

int 
qlua_sub(lua_State *L)
{
    q_op op = qt_sub[Op2Idx(qlua_gettype(L, 1), qlua_gettype(L, 2))];

    if (op)
        return op(L);
    else
        return luaL_error(L, "bad argument for subtraction");
}

static q_op qt_mul[(qOther + 1) * (qOther + 1)];

void
qlua_reg_mul(int ta, int tb, q_op op)
{
    qt_mul[Op2Idx(ta, tb)] = op;
}

int 
qlua_mul(lua_State *L)
{
    q_op op = qt_mul[Op2Idx(qlua_gettype(L, 1), qlua_gettype(L, 2))];

    if (op)
        return op(L);
    else
        return luaL_error(L, "bad argument for multiplication");
}

static q_op qt_div[(qOther + 1) * (qOther + 1)];

void
qlua_reg_div(int ta, int tb, q_op op)
{
    qt_div[Op2Idx(ta, tb)] = op;
}

int 
qlua_div(lua_State *L)
{
    q_op op = qt_div[Op2Idx(qlua_gettype(L, 1), qlua_gettype(L, 2))];

    if (op)
        return op(L);
    else
        return luaL_error(L, "bad argument for division");
}

static q_op qt_dot[(qOther + 1)];

void
qlua_reg_dot(int ta, q_op op)
{
    qt_dot[Op1Idx(ta)] = op;
}

static int 
q_dot(lua_State *L)
{
    q_op op = qt_dot[Op1Idx(qlua_gettype(L, 1))];

    if (op)
        return op(L);
    else
        return luaL_error(L, "bad argument for inner dot");
}

static struct luaL_Reg fQCD[] = {
    { "dot",    q_dot}, /* local inner dot */
    { NULL,     NULL}
};

/* environment setup */
static void
qlua_openlibs (lua_State *L) {
    static const luaL_Reg lualibs[] = {
        {"", luaopen_base},
        {LUA_LOADLIBNAME, luaopen_package},
        {LUA_TABLIBNAME, luaopen_table},
        /* {LUA_IOLIBNAME, luaopen_io}, */
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
        init_gamma,
        init_vector,
        init_lattice,
        init_latint,
        init_latrandom,
        init_latreal,
        init_latcomplex,
        init_latcolvec,
        init_latcolmat,
        init_latdirferm,
        init_latdirprop,
        init_latsubset,
        init_latmulti,
        init_qdpc_io,
        init_qdpcc_io,
        init_ddpairs_io,
        init_nersc_io,
        init_xml,
#ifdef HAS_AFF
        init_aff_io,
#endif
#ifdef HAS_CLOVER
        init_clover,
#endif
#ifdef HAS_MDWF
        init_mdwf,
#endif
        /* ZZZ add other packages here */
        NULL };

    int i;

    lua_gc(L, LUA_GCSTOP, 0);  /* stop collector during initialization */
    qlua_openlibs(L);  /* open libraries */
    luaL_register(L, qcdlib, fQCD);
    for (i = 0; qcd_inits[i]; i++) {
        lua_pushcfunction(L, qcd_inits[i]);
        lua_call(L, 0, 0);
    }
    lua_getglobal(L, qcdlib);
    lua_pushnumber(L, QDP_Nc);
    lua_setfield(L, -2, "Nc");
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
#ifdef HAS_MDWF
        fini_mdwf,
#endif
#ifdef HAS_CLOVER
        fini_clover,
#endif
#ifdef HAS_AFF
        fini_aff_io,
#endif
        fini_xml,
        fini_nersc_io,
        fini_ddpairs_io,
        fini_qdpcc_io,
        fini_qdpc_io,
        fini_latmulti,
        fini_latsubset,
        fini_latdirprop,
        fini_latdirferm,
        fini_latcolmat,
        fini_latcolvec,
        fini_latcomplex,
        fini_latreal,
        fini_latrandom,
        fini_latint,
        fini_lattice,
        fini_vector,
        fini_gamma,
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
    qlua_primary_node = QMP_is_primary_node();
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
            status = luaL_dofile(L, argv[i]);
            report(L, argv[i], status);
            if (status) {
                QDP_abort(1);
                break;
            }
        }
    }
    qlua_fini(L);
    lua_close(L);
end:
    QDP_finalize();
    return status;
}
