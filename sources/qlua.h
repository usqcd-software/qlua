#ifndef MARK_2DCAC914_635D_4D58_AA60_DC75CD13961F
#define MARK_2DCAC914_635D_4D58_AA60_DC75CD13961F

#define QDP_Precision 'D'
#define QDP_Nc 3
#define QDP_Ns 4 /* fermion dimension, seems undefined in QDP */
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
    qGamma,
    qVecInt,
    qVecReal,
    qVecComplex,
    qLatInt,
    qLatRandom,
    qLatReal,
    qLatComplex,
    qLatColVec,
    qLatColMat,
    qLatDirFerm,
    qLatDirProp,
    qReader,
    qWriter,
    qAffReader,
    qAffWriter,
    qClover,
    qMoebius,
    /* ZZZ add types for other packages here */
    qOther
};

extern const char *progname;
extern const char *qcdlib;
extern int qlua_primary_node;

void qlua_metatable(lua_State *L, const char *name, const luaL_Reg *table);
int qlua_lookup(lua_State *L, int idx, const char *table);
int qlua_gettype(lua_State *L, int idx);
void *qlua_malloc(lua_State *L, int size);
void qlua_free(lua_State *L, void *ptr);

int qlua_index(lua_State *L, int idx, const char *name, int mv);  /* k or -1 */
int qlua_checkindex(lua_State *L, int idx, const char *name, int mv);

int qlua_diracindex(lua_State *L, int idx);               /* d = 0 .. Nf - 1 */
int qlua_checkdiracindex(lua_State *L, int idx);
int qlua_colorindex(lua_State *L, int idx);               /* c = 0 .. Nc - 1 */
int qlua_checkcolorindex(lua_State *L, int idx);
int qlua_checkleftindex(lua_State *L, int idx);           /* a = 0 .. Nc - 1 */
int qlua_leftindex(lua_State *L, int idx);
int qlua_checkrightindex(lua_State *L, int idx);          /* b = 0 .. Nc - 1 */
int qlua_rightindex(lua_State *L, int idx);
int qlua_gammaindex(lua_State *L, int idx);                /* mu = 0..3 or 5 */
int qlua_checkgammaindex(lua_State *L, int idx);
int qlua_gammabinary(lua_State *L, int idx);                  /* n = 0 .. 15 */
int qlua_checkgammabinary(lua_State *L, int idx);

const char *qlua_checkstring(lua_State *L, int idx, const char *fmt, ...);
int qlua_checkint(lua_State *L, int idx, const char *fmt, ...);
double qlua_checknumber(lua_State *L, int idx, const char *fmt, ...);
void qlua_checktable(lua_State *L, int idx, const char *fmt, ...);

typedef int (*q_op)(lua_State *L);

void qlua_reg_add(int ta, int tb, q_op op);
int qlua_add(lua_State *L);
void qlua_reg_sub(int ta, int tb, q_op op);
int qlua_sub(lua_State *L);
void qlua_reg_mul(int ta, int tb, q_op op);
int qlua_mul(lua_State *L);
void qlua_reg_div(int ta, int tb, q_op op);
int qlua_div(lua_State *L);
void qlua_reg_mod(int ta, int tb, q_op op);
int qlua_mod(lua_State *L);

void qlua_reg_dot(int ta, q_op op);

int qlua_badconstr(lua_State *L, const char *name);
int qlua_badindex(lua_State *L, const char *type);

/* strict memory management: collect garbage befor any QDP operation */
#define CALL_QDP(L) do lua_gc(L, LUA_GCCOLLECT, 0); while (0)

#endif /* !defined(MARK_2DCAC914_635D_4D58_AA60_DC75CD13961F) */
