#include <stdarg.h>
#include <stdio.h>

#define lua_c
#include <lua.h>
#include <lauxlib.h>
#include <lualib.h>

#include <qlua.h>
#include <qcomplex.h>
#include <qvector.h>

static struct {
    int (*init)(lua_State *L);
} qcd_libs[] = {
    { init_complex },
    { init_vector },
    { NULL }
};

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
    case LUA_TUSERDATA: {
        if (qlua_type(L, idx, mtnComplex)) return qComplex;
        if (qlua_type(L, idx, mtnVecInt)) return qVecInt;
        if (qlua_type(L, idx, mtnVecDouble)) return qVecDouble;
        if (qlua_type(L, idx, mtnVecComplex)) return qVecComplex;
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
        /* add other packages here */
        { qOther,   qOther,   NULL}
    };

    return qlua_dispatch(L, qdiv_table, "division");
}

void
qlua_init(lua_State *L)
{
    int i;

    lua_gc(L, LUA_GCSTOP, 0);  /* stop collector during initialization */
    luaL_openlibs(L);  /* open libraries */
    for (i = 0; qcd_libs[i].init; i++) {
        lua_pushcfunction(L, qcd_libs[i].init);
        lua_call(L, 0, 0);
    }
    lua_gc(L, LUA_GCRESTART, 0);
}
