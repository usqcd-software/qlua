#ifndef MARK_2DCAC914_635D_4D58_AA60_DC75CD13961F
#define MARK_2DCAC914_635D_4D58_AA60_DC75CD13961F

#define QDP_Precision 'D'
#define QDP_Nc 3
#define QDP_Nf 4 /* fermion dimension, seems undefined in QDP */
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
    qLatRandom,
    qLatReal,
    qLatComplex,
    qLatColVec,
    qLatColMat,
    qLatDirFerm,
    qLatDirProp,
    /* ZZZ add types for other packages here */
    qOther
};

extern const char *progname;
extern const char *qcdlib;
/* add metatable names for other types here */

void qlua_metatable(lua_State *L, const char *name, const luaL_Reg *table);
int qlua_lookup(lua_State *L, int idx, const char *table);
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
