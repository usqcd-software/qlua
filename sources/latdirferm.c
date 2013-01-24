#include "modules.h"                                                 /* DEPS */
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
#include "seqdirferm.h"                                              /* DEPS */
#include "latmulti.h"                                                /* DEPS */
#include "qmp.h"
#include "qla.h"

#if USE_Nc2
#define QNc  '2'
#define Qcolors "2"
#define Qs(a)   a ## 2
#define Qx(a,b)  a ## 2 ## b
#define QC(x)    2
#define QNC(x)
#include "latdirferm-x.c"                                           /* DEPS */
#endif

#if USE_Nc3
#define QNc  '3'
#define Qcolors "3"
#define Qs(a)   a ## 3
#define Qx(a,b)  a ## 3 ## b
#define QC(x)    3
#define QNC(x)
#include "latdirferm-x.c"                                           /* DEPS */
#endif

#if USE_NcN
#define QNc  'N'
#define Qcolors "N"
#define Qs(a)   a ## N
#define Qx(a,b)  a ## N ## b
#define QC(x)    (x)->nc
#define QNC(x)   (x),
#include "latdirferm-x.c"                                           /* DEPS */
#endif

static int
do_gaussian(lua_State *L, mLatRandom *a, mLattice *S, int nc)
{
    switch (nc) {
#if USE_Nc2
    case 2: {
        mLatDirFerm2 *r = qlua_newLatDirFerm2(L, lua_gettop(L), 2);
        CALL_QDP(L);
        QDP_D2_D_eq_gaussian_S(r->ptr, a->ptr, *S->qss);
        return 1;
    }
#endif
#if USE_Nc3
    case 3: {
        mLatDirFerm3 *r = qlua_newLatDirFerm3(L, lua_gettop(L), 3);
        CALL_QDP(L);
        QDP_D3_D_eq_gaussian_S(r->ptr, a->ptr, *S->qss);
        return 1;
    }
#endif
#if USE_NcN
    default: {
        mLatDirFermN *r = qlua_newLatDirFermN(L, lua_gettop(L), nc);
        CALL_QDP(L);
        QDP_DN_D_eq_gaussian_S(r->ptr, a->ptr, *S->qss);
        return 1;
    }
#else
    default: return luaL_error(L, "bad number of colors");
#endif
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

    if (nc < 1)
        return luaL_error(L, "bad number of colors");

    return do_gaussian(L, a, S, nc);
}

/* Lattice Dirac Fermions:
 *  L:DiracFermion()             -- zero fermion in default colors
 *  L:DiracFermion(d)            -- constant fermion fill
 *  L:DiracFermion(D)            -- lattice fermion copy
 *  L:DiracFermion(C,{c=i,d=j})  -- default colors, D[i,j] = C zero otherwise
 *  L:DiracFermion(V,{d=j})      -- default colors, D[.,j] = V zero otherwise
 */
static int
q_latdirferm(lua_State *L)
{
    mLattice *S = qlua_checkLattice(L, 1);

    if (lua_gettop(L) == 2) {
        switch (qlua_qtype(L, 2)) {
#if USE_Nc2
        case qSeqDirFerm2: return q_latdirferm_seq_2(L, S, 2);
        case qLatDirFerm2: return q_latdirferm_lat_2(L, S, 2);
#endif
#if USE_Nc3
        case qSeqDirFerm3: return q_latdirferm_seq_3(L, S, 3);
        case qLatDirFerm3: return q_latdirferm_lat_3(L, S, 3);
#endif
#if USE_NcN
        case qSeqDirFermN: {
            mSeqDirFermN *v = qlua_checkSeqDirFermN(L, 2, -1);
            return q_latdirferm_seq_N(L, S, v->nc);
        }
        case qLatDirFermN: {
            mLatDirFermN *v = qlua_checkLatDirFermN(L, 2, S, -1);
            return q_latdirferm_lat_N(L, S, v->nc);
        }
#endif
        default:
            break;
        }
    }

    switch (S->nc) {
#if USE_Nc2
    case 2: return q_latdirferm_2(L, S, 2, 0);
#endif
#if USE_Nc3
    case 3: return q_latdirferm_3(L, S, 3, 0);
#endif
#if USE_NcN
    default: return q_latdirferm_N(L, S, S->nc, 0);
#else
    default: return luaL_error(L, "bad number of colors");
#endif
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
#if USE_Nc2 || USE_Nc3 || USE_NcN
    mLattice *S = qlua_checkLattice(L, 1);
#endif
    int nc = luaL_checkint(L, 2);

    switch (nc) {
#if USE_Nc2
    case 2: return q_latdirferm_2(L, S, 2, 1);
#endif
#if USE_Nc3
    case 3: return q_latdirferm_3(L, S, 3, 1);
#endif
#if USE_NcN
    default: return q_latdirferm_N(L, S, nc, 1);
#else
    default: return luaL_error(L, "bad number of colors");
#endif
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
    luaL_getmetatable(L, opLattice);
    luaL_register(L, NULL, fLatDirFerm);
    lua_pop(L, 1);
#if USE_Nc2
    qlua_reg_op2(ops2);
    qlua_reg_dot(qLatDirFerm2,  q_D_dot_2);
#endif
#if USE_Nc3
    qlua_reg_op2(ops3);
    qlua_reg_dot(qLatDirFerm3,  q_D_dot_3);
#endif
#if USE_NcN
    qlua_reg_op2(opsN);
    qlua_reg_dot(qLatDirFermN,  q_D_dot_N);
#endif

    return 0;
}

int
fini_latdirferm(lua_State *L)
{
    return 0;
}

