#include "modules.h"                                                 /* DEPS */
#include "qlua.h"                                                    /* DEPS */
#include "qcomplex.h"                                                /* DEPS */
#include "seqrandom.h"                                               /* DEPS */
#include "seqcolvec.h"                                               /* DEPS */
#include "seqcolmat.h"                                               /* DEPS */
#include "seqdirferm.h"                                              /* DEPS */
#include "qla.h"

#if USE_Nc2
#define QNc  '2'
#define Qcolors "2"
#define Qs(a)   a ## 2
#define Qx(a,b)  a ## 2 ## b
#define QC(x)    2
#define QNC(x)
#include "seqdirferm-x.c"                                           /* DEPS */
#endif

#if USE_Nc3
#define QNc  '3'
#define Qcolors "3"
#define Qs(a)   a ## 3
#define Qx(a,b)  a ## 3 ## b
#define QC(x)    3
#define QNC(x)
#include "seqdirferm-x.c"                                           /* DEPS */
#endif

#if USE_NcN
#define QNc  'N'
#define Qcolors "N"
#define Qs(a)   a ## N
#define Qx(a,b)  a ## N ## b
#define QC(x)    (x)->nc
#define QNC(x)   (x), 
#include "seqdirferm-x.c"                                           /* DEPS */
#endif

/* Random vectors of NC colors
 * S:gaussian_DiracFermionN(nc)
 */
int
q_d_gaussian_N(lua_State *L)
{
    mSeqRandom *a = qlua_checkSeqRandom(L, 1);
    int nc = luaL_checkint(L, 2);

    switch (nc) {
#if USE_Nc2
    case 2: {
        mSeqDirFerm2 *r = qlua_newSeqDirFerm2(L, 2);
        QLA_D2_D_eq_gaussian_S(r->ptr, &a->state);
        return 1;
    }
#endif
#if USE_Nc3
    case 3: {
        mSeqDirFerm3 *r = qlua_newSeqDirFerm3(L, 3);
        QLA_D3_D_eq_gaussian_S(r->ptr, &a->state);
        return 1;
    }
#endif
#if USE_NcN
    default: {
        mSeqDirFermN *r = qlua_newSeqDirFermN(L, nc);
        QLA_DN_D_eq_gaussian_S(nc, r->ptr, &a->state);
        return 1;
    }
#else
    default: return luaL_error(L, "bad number of colors");
#endif
    }
}

/* Sequential Dirac Fermions:
 *  qcd.DiracFermionN(n)              -- zero fermion in n colors
 *  qcd.DiracFermionN(n,C,{c=i,d=j})  -- n colors, D[i,j] = C zero otherwise
 *  qcd.DiracFermionN(n,V,{d=j})      -- n colors, D[.,j] = V zero otherwise
 */
static int
q_seqdirfermN(lua_State *L)
{
    int nc = luaL_checkint(L, 1);
    lua_remove(L, 1);

    switch (nc) {
#if USE_Nc2
    case 2: return q_seqdirferm_2(L, 2);
#endif
#if USE_Nc3
    case 3: return q_seqdirferm_3(L, 3);
#endif
#if USE_NcN
    default: return q_seqdirferm_N(L, nc);
#else
    default: return luaL_error(L, "bad number of colors");
#endif
    }
}

static struct luaL_Reg fSeqDirFerm[] = {
    { "DiracFermionN",      q_seqdirfermN },
    { NULL,                 NULL          }
};

int
init_seqdirferm(lua_State *L)
{
    luaL_register(L, qcdlib, fSeqDirFerm);
#if USE_Nc2
    qlua_metatable(L, mtnSeqDirFerm2, mtSeqDirFerm2, qSeqDirFerm2);
    qlua_reg_op2(ops2);
    qlua_reg_dot(qSeqDirFerm2,  q_d_dot_2);
#endif
#if USE_Nc3
    qlua_metatable(L, mtnSeqDirFerm3, mtSeqDirFerm3, qSeqDirFerm3);
    qlua_reg_op2(ops3);
    qlua_reg_dot(qSeqDirFerm3,  q_d_dot_3);
#endif
#if USE_NcN
    qlua_metatable(L, mtnSeqDirFermN, mtSeqDirFermN, qSeqDirFermN);
    qlua_reg_op2(opsN);
    qlua_reg_dot(qSeqDirFermN,  q_d_dot_N);
#endif

    return 0;
}

int
fini_seqdirferm(lua_State *L)
{
    return 0;
}
