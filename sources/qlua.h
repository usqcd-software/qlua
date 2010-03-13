#ifndef MARK_2DCAC914_635D_4D58_AA60_DC75CD13961F
#define MARK_2DCAC914_635D_4D58_AA60_DC75CD13961F

#define QDP_Precision 'D' /* DO NOT CHANGE THESE DEFAULTS */
#define QDP_Nc 3          /* DO NOT CHANGE THESE DEFAULTS */
#define QDP_Ns 4 /* fermion dimension, seems undefined in QDP */
#include <stdarg.h>
#include <stdio.h>

#define lua_c
#include <lua.h>
#include <lauxlib.h>
#include <lualib.h>

#include <qdp.h>

typedef enum {
    /* start with all type that have any of arithmetic operations defined */
    qReal,
    qComplex,
    qGamma,
    qMatReal,
    qMatComplex,
    qLatInt,
    qLatRealF,
    qLatRealD,
    qLatComplexF,
    qLatComplexD,
    qLatColVecF2,
    qLatColVecD2,
    qLatColVecF3,
    qLatColVecD3,
    qLatColVecFn,
    qLatColVecDn,
    qLatColMatF2,
    qLatColMatD2,
    qLatColMatF3,
    qLatColMatD3,
    qLatColMatFn,
    qLatColMatDn,
    qLatDirFermF2,
    qLatDirFermD2,
    qLatDirFermF3,
    qLatDirFermD3,
    qLatDirFermFn,
    qLatDirFermDn,
    qLatDirPropF2,
    qLatDirPropD2,
    qLatDirPropF3,
    qLatDirPropD3,
    qLatDirPropFn,
    qLatDirPropDn,
    qOther,          /* no operations for this type */
    qArithTypeCount, /* number of types in arith dispatch tables */
    qLattice,
    qString,
    qTable,
    qVecInt,
    qVecReal,
    qVecComplex,
    qLatRandom,
    qReader,
    qWriter,
    qAffReader,
    qAffWriter,
    qClover,
    /* ZZZ add types for other packages here */
    qNoType
} QLUA_Type;

extern const char *progname;
extern const char *qcdlib;
extern int qlua_primary_node;

void qlua_metatable(lua_State *L, const char *name, const luaL_Reg *table,
                    QLUA_Type t_id);
int qlua_lookup(lua_State *L, int idx, const char *table);
QLUA_Type qlua_qtype(lua_State *L, int idx);
QLUA_Type qlua_atype(lua_State *L, int idx); /* non-arith types => qOther */
const char *qlua_ptype(lua_State *L, int idx);
void *qlua_malloc(lua_State *L, int size);
void qlua_free(lua_State *L, void *ptr);

int qlua_index(lua_State *L, int idx, const char *name, int mv);  /* k or -1 */
int qlua_checkindex(lua_State *L, int idx, const char *name, int mv);
void qlua_checkindex2(lua_State *L, int idx, const char *nm, int *i0, int *i1);

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

typedef int (*q_op)(lua_State *L, void *env);
typedef struct {
    q_op  op;
    void *env;
} q_op_table;

void qlua_reg_add(QLUA_Type ta, QLUA_Type tb, q_op op, void *env);
int qlua_add(lua_State *L);
void qlua_reg_sub(QLUA_Type ta, QLUA_Type tb, q_op op, void *env);
int qlua_sub(lua_State *L);
void qlua_reg_mul(QLUA_Type ta, QLUA_Type tb, q_op op, void *env);
int qlua_mul(lua_State *L);
void qlua_reg_div(QLUA_Type ta, QLUA_Type tb, q_op op, void *env);
int qlua_div(lua_State *L);
void qlua_reg_mod(QLUA_Type ta, QLUA_Type tb, q_op op, void *env);
int qlua_mod(lua_State *L);

void qlua_reg_dot(QLUA_Type ta, q_op op, void *env);

int qlua_badconstr(lua_State *L, const char *name);
int qlua_badindex(lua_State *L, const char *type);

/* generic error reporter for __newindex in metatables */
int qlua_nowrite(lua_State *L);

/* strict memory management: collect garbage befor any QDP operation */
#define CALL_QDP(L) do lua_gc(L, LUA_GCCOLLECT, 0); while (0)

#endif /* !defined(MARK_2DCAC914_635D_4D58_AA60_DC75CD13961F) */
