#include <qlua.h>
#include <qcomplex.h>
#include <qvector.h>
#include <qlint.h>
/* include other package headers here */

const char *progname = "qlua";
const char *qcdlib = "qcd";

/* reporting */
void
message(const char *fmt, ...)
{
    if (QDP_this_node == 0) {
        va_list va;

        va_start(va, fmt);
        vfprintf(stderr, fmt, va);
        va_end(va);
    }
}

void
report(lua_State *L, const char *fname, int status)
{
    if (QDP_this_node == 0) {
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
        if (qlua_type(L, idx, mtnComplex)) return qComplex;
        if (qlua_type(L, idx, mtnVecInt)) return qVecInt;
        if (qlua_type(L, idx, mtnVecDouble)) return qVecDouble;
        if (qlua_type(L, idx, mtnVecComplex)) return qVecComplex;
        if (qlua_type(L, idx, mtnLatInt)) return qLatInt;
    }
    default:
        return qOther;
    }
}

/* generic operations dispatchers */
typedef struct {
    int ta, tb;
    int (*op)(lua_State *L);
} q_Op2Table;

int
qlua_dispatch(lua_State *L, const q_Op2Table *t, const char *name)
{
    int i;
    int ta = qlua_gettype(L, 1);
    int tb = qlua_gettype(L, 2);

    for (i = 0; t[i].op; i++) {
        if ((t[i].ta == ta) && (t[i].tb == tb))
            return t[i].op(L);
    }

    return luaL_error(L, "bad arguments for %s", name);
}

int
qlua_add(lua_State *L)
{
    static const q_Op2Table qadd_table[] = {
        { qReal,    qComplex, q_r_add_c },
        { qComplex, qReal,    q_c_add_r },
        { qComplex, qComplex, q_c_add_c },
        { qLatInt,  qLatInt,  q_I_add_I },
        /* add other packages here */
        { qOther,   qOther,   NULL}
    };
    return qlua_dispatch(L, qadd_table, "addition");
}

int
qlua_sub(lua_State *L)
{
    static const q_Op2Table qsub_table[] = {
        { qReal,    qComplex, q_r_sub_c },
        { qComplex, qReal,    q_c_sub_r },
        { qComplex, qComplex, q_c_sub_c },
        { qLatInt,  qLatInt,  q_I_sub_I },
        /* add other packages here */
        { qOther,   qOther,   NULL}
    };

    return qlua_dispatch(L, qsub_table, "subtraction");
}

int
qlua_mul(lua_State *L)
{
    static const q_Op2Table qmul_table[] = {
        { qReal,    qComplex, q_r_mul_c },
        { qComplex, qReal,    q_c_mul_r },
        { qComplex, qComplex, q_c_mul_c },
        { qReal,    qLatInt,  q_i_mul_I },
        { qLatInt,  qReal,    q_I_mul_i },
        { qLatInt,  qLatInt,  q_I_mul_I },
        /* add other packages here */
        { qOther,   qOther,   NULL}
    };

    return qlua_dispatch(L, qmul_table, "multiplication");
}

int
qlua_div(lua_State *L)
{
    static const q_Op2Table qdiv_table[] = {
        { qReal,    qComplex, q_r_div_c },
        { qComplex, qReal,    q_c_div_r },
        { qComplex, qComplex, q_c_div_c },
        { qLatInt,  qLatInt,  q_I_div_I },
        /* add other packages here */
        { qOther,   qOther,   NULL}
    };

    return qlua_dispatch(L, qdiv_table, "division");
}

/* environment setup */
void
qlua_init(lua_State *L)
{
    static const struct {
        int (*init)(lua_State *L);
    } qcd_inits[] = {
        { init_complex },
        { init_vector },
        { init_latint },
        /* add other packages here */
        { NULL }
    };

    int i;

    lua_gc(L, LUA_GCSTOP, 0);  /* stop collector during initialization */
    luaL_openlibs(L);  /* open libraries */
    for (i = 0; qcd_inits[i].init; i++) {
        lua_pushcfunction(L, qcd_inits[i].init);
        lua_call(L, 0, 0);
    }
    lua_gc(L, LUA_GCRESTART, 0);
}

/* cleanup (mostly housekeeping to make various tools happy */
void
qlua_fini(lua_State *L)
{
    static struct {
        int (*fini)(lua_State *L);
    } qcd_finis[] = { /* keep it in the reverse order with respect to init */
        /* add other packages here */
        { fini_latint },
        { fini_vector },
        { fini_complex },
        { NULL }
    };
    int i;

    for (i = 0; qcd_finis[i].fini; i++) {
        lua_pushcfunction(L, qcd_finis[i].fini);
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
    L = lua_open();
    if (L == NULL) {
        message("can not create Lua state");
        goto end;
    }
    qlua_init(L);  /* open libraries */

    for (i = 1; i < argc; i++) {
        status = luaL_dofile(L, argv[i]);
        report(L, argv[i], status);
        if (status)
            break;
    }
    qlua_fini(L);
    lua_close(L);
end:
    QDP_finalize();
    return status;
}
