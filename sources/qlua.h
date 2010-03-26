#ifndef MARK_2DCAC914_635D_4D58_AA60_DC75CD13961F
#define MARK_2DCAC914_635D_4D58_AA60_DC75CD13961F

#define QDP_Nc 'Z' /* Try to confuse QDP & QLA to warn about undefined Nc */
#define QDP_Precision 'D' /* DO NOT CHANGE THESE DEFAULTS */
#define QDP_Ns 4 /* fermion dimension, seems undefined in QDP */
#include <stdarg.h>
#include <stdio.h>

#define lua_c
#include "lua.h"
#include "lauxlib.h"
#include "lualib.h"

#include "qdp.h"
#include "qdp_d2.h"
#include "qdp_d3.h"
#include "qdp_dn.h"
#include "qla_types.h"
#include "qla.h"
#include "qla_d2.h"
#include "qla_d3.h"
#include "qla_dn.h"

typedef enum {
    /* start with all type that have any of arithmetic operations defined */
    qReal,
    qComplex,
    qGamma,
    qMatReal,
    qMatComplex,
    qLatInt,
    qLatReal,
    qLatComplex,
    qLatColVec2,
    qLatColVec3,
    qLatColVecN,
    qLatColMat2,
    qLatColMat3,
    qLatColMatN,
    qLatDirFerm2,
    qLatDirFerm3,
    qLatDirFermN,
    qLatDirProp2,
    qLatDirProp3,
    qLatDirPropN,
    qOther,          /* no operations for this type */
    qArithTypeCount, /* number of types in arith dispatch tables */
    qLattice,
    qLatMulti,
    qLatSubset,
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
    qDeflator,
    qDeflatorState,
    /* ZZZ add types for other packages here */
    qNoType
} QLUA_Type;

extern const char *progname;
extern const char *qcdlib;
extern int qlua_primary_node;

void qlua_metatable(lua_State *L,
                    const char *name,
                    const luaL_Reg *table,
                    QLUA_Type t_id);
void qlua_selftable(lua_State *L,
                    const luaL_Reg *table,
                    QLUA_Type t_id);
void qlua_latticetable(lua_State *L,
                       const luaL_Reg *table,
                       QLUA_Type t_id,
                       int Sidx);
int qlua_getlattice(lua_State *L, int idx); /* => [idx]'s lattice object */
int qlua_lookup(lua_State *L, int idx, const char *table);
int qlua_selflookup(lua_State *L, int idx, const char *key);
void *qlua_checkLatticeType(lua_State *L,
                            int idx,
                            QLUA_Type t_id,
                            const char *name);
void qlua_createLatticeTable(lua_State *L,
                             int Sidx,
                             const struct luaL_Reg *ft,
                             QLUA_Type t_id,
                             const char *name);
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
int qlua_colorindex(lua_State *L, int idx, int nc);       /* c = 0 .. Nc - 1 */
int qlua_checkcolorindex(lua_State *L, int idx, int nc);
int qlua_checkleftindex(lua_State *L, int idx, int nc);   /* a = 0 .. Nc - 1 */
int qlua_leftindex(lua_State *L, int idx, int nc);
int qlua_checkrightindex(lua_State *L, int idx, int nc);  /* b = 0 .. Nc - 1 */
int qlua_rightindex(lua_State *L, int idx, int nc);
int qlua_gammaindex(lua_State *L, int idx);                /* mu = 0..3 or 5 */
int qlua_checkgammaindex(lua_State *L, int idx);
int qlua_gammabinary(lua_State *L, int idx);                  /* n = 0 .. 15 */
int qlua_checkgammabinary(lua_State *L, int idx);

const char *qlua_checkstring(lua_State *L, int idx, const char *fmt, ...);
int qlua_checkint(lua_State *L, int idx, const char *fmt, ...);
double qlua_checknumber(lua_State *L, int idx, const char *fmt, ...);
void qlua_checktable(lua_State *L, int idx, const char *fmt, ...);

typedef int (*q_op)(lua_State *L);

typedef struct {
    q_op     *table;
    QLUA_Type ta;
    QLUA_Type tb;
    q_op      op;
} QLUA_Op2;

extern q_op qlua_add_table[];
extern q_op qlua_sub_table[];
extern q_op qlua_mul_table[];
extern q_op qlua_div_table[];
extern q_op qlua_mod_table[];

void qlua_reg_op2(const QLUA_Op2 *ops);

int qlua_add(lua_State *L);
int qlua_sub(lua_State *L);
int qlua_mul(lua_State *L);
int qlua_div(lua_State *L);
int qlua_mod(lua_State *L);

void qlua_reg_dot(QLUA_Type ta, q_op op);

int qlua_badconstr(lua_State *L, const char *name);
int qlua_badindex(lua_State *L, const char *type);

/* generic error reporter for __newindex in metatables */
int qlua_nowrite(lua_State *L);

/* strict memory management: collect garbage befor any QDP operation */
#define CALL_QDP(L) do lua_gc(L, LUA_GCCOLLECT, 0); while (0)

#endif /* !defined(MARK_2DCAC914_635D_4D58_AA60_DC75CD13961F) */
