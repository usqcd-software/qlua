#ifndef MARK_2DCAC914_635D_4D58_AA60_DC75CD13961F
#define MARK_2DCAC914_635D_4D58_AA60_DC75CD13961F

#include <stdarg.h>
#include <stdio.h>

#define lua_c
#include <lua.h>
#include <lauxlib.h>
#include <lualib.h>

#include <qdp.h>

enum {
    qReal,
    qString,
    qTable,
    qComplex,
    qVecInt,
    qVecDouble,
    qVecComplex,
    qLatInt,
    qLatReal,
    qLatComplex,
    qLatColorVector,
    qLatColorMatrix,
    qLatDiracFermion,
    qLatHalfFermion,
    qLatDiracPropagator,
    /* add types for other packages here */
    qOther
};

extern const char *progname;
extern const char *qcdlib;
#if 0 /* XXX */
extern const char *mtnLatReal;
extern const char *mtnLatComplex;
extern const char *mtnLatColorVector;
extern const char *mtnLatColorMatrix;
extern const char *mtnLatDiracFermion;
extern const char *mtnLatHalfFermion;
extern const char *mtnLatDiracPropagator;
#endif
/* add metatable names for other types here */

void qlua_metatable(lua_State *L, const char *name, const luaL_Reg *table);
void qlua_init(lua_State *L);
void qlua_fini(lua_State *L);
int qlua_gettype(lua_State *L, int idx);
void *qlua_malloc(lua_State *L, int size);
void qlua_free(lua_State *L, void *ptr);

int qlua_add(lua_State *L);
int qlua_sub(lua_State *L);
int qlua_mul(lua_State *L);
int qlua_div(lua_State *L);

#endif /* !defined(MARK_2DCAC914_635D_4D58_AA60_DC75CD13961F) */
