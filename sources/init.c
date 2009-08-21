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
