#include "qlua.h"                                                    /* DEPS */
#include "qcomplex.h"                                                /* DEPS */
#include "lattice.h"                                                 /* DEPS */
#include "latint.h"                                                  /* DEPS */
#include "latrandom.h"                                               /* DEPS */
#include "latreal.h"                                                 /* DEPS */
#include "latcomplex.h"                                              /* DEPS */
#include "latcolvec.h"                                               /* DEPS */
#include "latcolmat.h"                                               /* DEPS */
#include "latdirferm.h"                                              /* DEPS */
#include "qmp.h"

#define QNc  '2'
#define Qcolors "2"
#define Qs(a)   a ## 2
#define Qx(a,b)  a ## 2 ## b
#define QC(x)    2
#include "latdirferm-x.c"                                           /* DEPS */

#define QNc  '3'
#define Qcolors "3"
#define Qs(a)   a ## 3
#define Qx(a,b)  a ## 3 ## b
#define QC(x)    3
#include "latdirferm-x.c"                                           /* DEPS */

#define QNc  'N'
#define Qcolors "N"
#define Qs(a)   a ## N
#define Qx(a,b)  a ## N ## b
#define QC(x)    (x)->nc
#include "latdirferm-x.c"                                           /* DEPS */

static int
do_gaussian(lua_State *L, mLatRandom *a, mLattice *S, int nc)
{
    switch (nc) {
    case 2: {
        mLatDirFerm2 *r = qlua_newLatDirFerm2(L, lua_gettop(L), 2);
        CALL_QDP(L);
        QDP_D2_D_eq_gaussian_S(r->ptr, a->ptr, *S->qss);
        return 1;
    }
    case 3: {
        mLatDirFerm3 *r = qlua_newLatDirFerm3(L, lua_gettop(L), 3);
        CALL_QDP(L);
        QDP_D3_D_eq_gaussian_S(r->ptr, a->ptr, *S->qss);
        return 1;
    }
    default: {
        mLatDirFermN *r = qlua_newLatDirFermN(L, lua_gettop(L), nc);
        CALL_QDP(L);
        QDP_DN_D_eq_gaussian_S(r->ptr, a->ptr, *S->qss);
        return 1;
    }
    }
}

/* Random vectors of default colors
 * S:gaussian_DiracFermion()
 */
int
q_D_gaussian(lua_State *L)
{
    mLatRandom *a = qlua_checkLatRandom(L, 1, NULL);
    mLattice *S = qlua_ObjLattice(L, 1);

    return do_gaussian(L, a, S, S->nc);
}

/* Random vectors of NC colors
 * S:gaussian_DiracFermionN(nc)
 */
int
q_D_gaussian_N(lua_State *L)
{
    mLatRandom *a = qlua_checkLatRandom(L, 1, NULL);
    mLattice *S = qlua_ObjLattice(L, 1);
    int nc = luaL_checkint(L, 2);

    if (nc < 2)
        return luaL_error(L, "bad number of colors");

    return do_gaussian(L, a, S, nc);
}

/* Lattice Dirac Fermions:
 *  L:DiracFermion()             -- zero fermion in default colors
 *  L:DiracFermion(C,{c=i,d=j})  -- default colors, D[i,j] = C zero otherwise
 *  L:DiracFermion(V,{d=j})      -- default colors, D[.,j] = V zero otherwise
 */
static int
q_latdirferm(lua_State *L)
{
    mLattice *S = qlua_checkLattice(L, 1);

    switch (S->nc) {
    case 2: return q_latdirferm_2(L, S, 2, 0);
    case 3: return q_latdirferm_3(L, S, 3, 0);
    default: return q_latdirferm_N(L, S, S->nc, 0);
    }
}

/* Lattice Dirac Fermions:
 *  L:DiracFermionN(n)              -- zero fermion in n colors
 *  L:DiracFermionN(n,C,{c=i,d=j})  -- n colors, D[i,j] = C zero otherwise
 *  L:DiracFermionN(n,V,{d=j})      -- n colors, D[.,j] = V zero otherwise
 */
static int
q_latdirfermN(lua_State *L)
{
    mLattice *S = qlua_checkLattice(L, 1);
    int nc = luaL_checkint(L, 2);

    switch (nc) {
    case 2: return q_latdirferm_2(L, S, 2, 1);
    case 3: return q_latdirferm_3(L, S, 3, 1);
    default: return q_latdirferm_N(L, S, nc, 1);
    }
}

static struct luaL_Reg fLatDirFerm[] = {
    { "DiracFermion",       q_latdirferm  },
    { "DiracFermionN",      q_latdirfermN },
    { NULL,                 NULL          }
};

int
init_latdirferm(lua_State *L)
{
    static const QLUA_Op2 ops[] = {
        { qlua_add_table, qLatDirFerm2,  qLatDirFerm2,  q_D_add_D_2 },
        { qlua_add_table, qLatDirFerm3,  qLatDirFerm3,  q_D_add_D_3 },
        { qlua_add_table, qLatDirFermN,  qLatDirFermN,  q_D_add_D_N },
        { qlua_sub_table, qLatDirFerm2,  qLatDirFerm2,  q_D_sub_D_2 },
        { qlua_sub_table, qLatDirFerm3,  qLatDirFerm3,  q_D_sub_D_3 },
        { qlua_sub_table, qLatDirFermN,  qLatDirFermN,  q_D_sub_D_N },
        { qlua_mul_table, qReal,         qLatDirFerm2,  q_r_mul_D_2 },
        { qlua_mul_table, qReal,         qLatDirFerm3,  q_r_mul_D_3 },
        { qlua_mul_table, qReal,         qLatDirFermN,  q_r_mul_D_N },
        { qlua_mul_table, qLatDirFerm2,  qReal,         q_D_mul_r_2 },
        { qlua_mul_table, qLatDirFerm3,  qReal,         q_D_mul_r_3 },
        { qlua_mul_table, qLatDirFermN,  qReal,         q_D_mul_r_N },
        { qlua_mul_table, qComplex,      qLatDirFerm2,  q_c_mul_D_2 },
        { qlua_mul_table, qComplex,      qLatDirFerm3,  q_c_mul_D_3 },
        { qlua_mul_table, qComplex,      qLatDirFermN,  q_c_mul_D_N },
        { qlua_mul_table, qLatDirFerm2,  qComplex,      q_D_mul_c_2 },
        { qlua_mul_table, qLatDirFerm3,  qComplex,      q_D_mul_c_3 },
        { qlua_mul_table, qLatDirFermN,  qComplex,      q_D_mul_c_N },
        { qlua_mul_table, qLatReal,      qLatDirFerm2,  q_R_mul_D_2 },
        { qlua_mul_table, qLatReal,      qLatDirFerm3,  q_R_mul_D_3 },
        { qlua_mul_table, qLatReal,      qLatDirFermN,  q_R_mul_D_N },
        { qlua_mul_table, qLatDirFerm2,  qLatReal,      q_D_mul_R_2 },
        { qlua_mul_table, qLatDirFerm3,  qLatReal,      q_D_mul_R_3 },
        { qlua_mul_table, qLatDirFermN,  qLatReal,      q_D_mul_R_N },
        { qlua_mul_table, qLatComplex,   qLatDirFerm2,  q_C_mul_D_2 },
        { qlua_mul_table, qLatComplex,   qLatDirFerm3,  q_C_mul_D_3 },
        { qlua_mul_table, qLatComplex,   qLatDirFermN,  q_C_mul_D_N },
        { qlua_mul_table, qLatDirFerm2,  qLatComplex,   q_D_mul_C_2 },
        { qlua_mul_table, qLatDirFerm3,  qLatComplex,   q_D_mul_C_3 },
        { qlua_mul_table, qLatDirFermN,  qLatComplex,   q_D_mul_C_N },
        { qlua_mul_table, qLatColMat2,   qLatDirFerm2,  q_M_mul_D_2 },
        { qlua_mul_table, qLatColMat3,   qLatDirFerm3,  q_M_mul_D_3 },
        { qlua_mul_table, qLatColMatN,   qLatDirFermN,  q_M_mul_D_N },
        { qlua_div_table, qLatDirFerm2,  qReal,         q_D_div_r_2 },
        { qlua_div_table, qLatDirFerm3,  qReal,         q_D_div_r_3 },
        { qlua_div_table, qLatDirFermN,  qReal,         q_D_div_r_N },
        { qlua_div_table, qLatDirFerm2,  qComplex,      q_D_div_c_2 },
        { qlua_div_table, qLatDirFerm3,  qComplex,      q_D_div_c_3 },
        { qlua_div_table, qLatDirFermN,  qComplex,      q_D_div_c_N },
        { NULL,           qNoType,       qNoType,       NULL        }
    };
    luaL_getmetatable(L, opLattice);
    luaL_register(L, NULL, fLatDirFerm);
    lua_pop(L, 1);
    qlua_reg_op2(ops);
    qlua_reg_dot(qLatDirFerm2,  q_D_dot_2);
    qlua_reg_dot(qLatDirFerm3,  q_D_dot_3);
    qlua_reg_dot(qLatDirFermN,  q_D_dot_N);

    return 0;
}

int
fini_latdirferm(lua_State *L)
{
    return 0;
}

