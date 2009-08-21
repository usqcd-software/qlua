#ifndef MARK_2DCAC914_635D_4D58_AA60_DC75CD13961F
#define MARK_2DCAC914_635D_4D58_AA60_DC75CD13961F

#include <stdarg.h>
#include <stdio.h>

#define lua_c
#include <lua.h>
#include <lauxlib.h>
#include <lualib.h>

extern const char *progname;
extern const char *qcdlib;
void qlua_metatable(lua_State *L, const char *name, const luaL_Reg *table);
void qlua_init(lua_State *L);

#endif /* !defined(MARK_2DCAC914_635D_4D58_AA60_DC75CD13961F) */
