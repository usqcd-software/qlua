#include <qlua.h>                                                    /* DEPS */
#include <lattice.h>                                                 /* DEPS */
#include <latint.h>                                                  /* DEPS */
#include <latmulti.h>                                                /* DEPS */

/* QLUA multisets are not QDP multisets! */

const char mtnLatMulti[] = "lattice.multi";

void
qlua_checkLatMulti(lua_State *L, int idx)
{
    luaL_checktype(L, idx, LUA_TTABLE);
    if (lua_getmetatable(L, idx) == 0)
        goto no_multi;
    luaL_getmetatable(L, mtnLatMulti);
    if (lua_equal(L, -1, -2) == 0)
        goto no_multi;
    lua_pop(L, 2);

    return;
no_multi:
    luaL_error(L, "lattice:MultiSet expected");
}

int
qlua_LatMultiSize(lua_State *L, int idx)
{
    int v;

    qlua_checkLatMulti(L, idx);
    lua_pushnumber(L, 1);
    lua_gettable(L, idx);

    v = luaL_checkint(L, -1);
    lua_pop(L, 1);

    return v;
}

mLatInt *
qlua_LatMultiIndex(lua_State *L, int idx)
{
    mLatInt *v;

    qlua_checkLatMulti(L, idx);
    lua_pushnumber(L, 2);
    lua_gettable(L, idx);
    v = qlua_checkLatInt(L, -1);
    lua_pop(L, 1); /* if MultiSet is collected, v may be collected as well! */
    
    return v;
}

static int
q_multi_fmt(lua_State *L)
{
    char fmt[72];
    int size = qlua_LatMultiSize(L, 1);

    sprintf(fmt, "MultiSet(%d,...)", size);
    lua_pushstring(L, fmt);
    
    return 1;
}

static int
q_latmulti(lua_State *L)
{
    int size = luaL_checkint(L, 2);
    
    qlua_checkLatInt(L, 3); /* index */
    lua_createtable(L, 2, 0);
    lua_pushnumber(L, size);
    lua_rawseti(L, -2, 1);
    lua_pushvalue(L, 3);
    lua_rawseti(L, -2, 2);
    luaL_getmetatable(L, mtnLatMulti);
    lua_setmetatable(L, -2);

    return 1;
}

static struct luaL_Reg mtLatMulti[] = {
    { "__tostring",     q_multi_fmt },
    { NULL,             NULL        }
};

static struct luaL_Reg fLatMulti[] = {
    { "MultiSet",    q_latmulti },
    { NULL,          NULL       }
};

int
init_latmulti(lua_State *L)
{
    luaL_getmetatable(L, opLattice);
    luaL_register(L, NULL, fLatMulti);
    lua_pop(L, 1);
    qlua_metatable(L, mtnLatMulti, mtLatMulti);

    return 0;
}

int
fini_latmulti(lua_State *L)
{
    return 0;
}
