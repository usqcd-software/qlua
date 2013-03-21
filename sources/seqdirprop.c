#include "modules.h"                                                 /* DEPS */
#include "qlua.h"                                                    /* DEPS */
#include "qcomplex.h"                                                /* DEPS */
#include "seqrandom.h"                                               /* DEPS */
#include "seqcolmat.h"                                               /* DEPS */
#include "seqdirferm.h"                                              /* DEPS */
#include "seqdirprop.h"                                              /* DEPS */

#if USE_Nc2
#define QNc  '2'
#define Qcolors "2"
#define Qs(a)   a ## 2
#define Qx(a,b)  a ## 2 ## b
#define QC(x)    2
#define QNC(x)
#include "seqdirprop-x.c"                                           /* DEPS */
#endif

#if USE_Nc3
#define QNc  '3'
#define Qcolors "3"
#define Qs(a)   a ## 3
#define Qx(a,b)  a ## 3 ## b
#define QC(x)    3
#define QNC(x)
#include "seqdirprop-x.c"                                           /* DEPS */
#endif

#if USE_NcN
#define QNc  'N'
#define Qcolors "N"
#define Qs(a)   a ## N
#define Qx(a,b)  a ## N ## b
#define QC(x)    (x)->nc
#define QNC(x)   (x),
#include "seqdirprop-x.c"                                           /* DEPS */
#endif

/* Random vectors of NC colors
 * S:gaussian_DiracPropagatorN(nc)
 */
int
q_p_gaussian_N(lua_State *L)
{
    mSeqRandom *a = qlua_checkSeqRandom(L, 1);
    int nc = luaL_checkint(L, 2);

    switch (nc) {
#if USE_Nc2
    case 2: {
        mSeqDirProp2 *r = qlua_newSeqDirProp2(L, 2);
        QLA_D2_P_eq_gaussian_S(r->ptr, &a->state);
        return 1;
    }
#endif
#if USE_Nc3
    case 3: {
        mSeqDirProp3 *r = qlua_newSeqDirProp3(L, 3);
        QLA_D3_P_eq_gaussian_S(r->ptr, &a->state);
        return 1;
    }
#endif
#if USE_NcN
    default: {
        mSeqDirPropN *r = qlua_newSeqDirPropN(L, nc);
        QLA_DN_P_eq_gaussian_S(nc, r->ptr, &a->state);
        return 1;
    }
#else
    default:
        return luaL_error(L, "bad number of colors");
#endif
    }
}

/* Sequential Dirac Propagators:
 *  qcd.DiracPropagatorN(n)              -- zero propagator in n colors
 *  qcd.DiracPropagatorN(n,M)            -- n colors, P[.,k,.,k] = M for all k
 *  qcd.DiracPropagatorN(n,D,{c=i,d=j})  -- n colors, P[i,j,.,.] = D
 */
static int
q_seqdirpropN(lua_State *L)
{
    int nc = luaL_checkint(L, 2);

    switch (nc) {
#if USE_Nc2
    case 2: return q_seqdirprop_2(L, 2);
#endif
#if USE_Nc3
    case 3: return q_seqdirprop_3(L, 3);
#endif
#if USE_NcN
    default: return q_seqdirprop_N(L, nc);
#else
    default: return luaL_error(L, "bad number of colors");
#endif
    }
}

static struct luaL_Reg fSeqDirProp[] = {
    { "DiracPropagatorN",  q_seqdirpropN },
    { NULL,                NULL          }
};

int
init_seqdirprop(lua_State *L)
{
    luaL_register(L, qcdlib, fSeqDirProp);
#if USE_Nc2
    qlua_metatable(L, mtnSeqDirProp2, mtSeqDirProp2, qSeqDirProp2);
    qlua_reg_op2(ops2);
    qlua_reg_dot(qSeqDirProp2,  q_p_dot_2);
#endif
#if USE_Nc3
    qlua_metatable(L, mtnSeqDirProp3, mtSeqDirProp3, qSeqDirProp3);
    qlua_reg_op2(ops3);
    qlua_reg_dot(qSeqDirProp3,  q_p_dot_3);
#endif
#if USE_NcN
    qlua_metatable(L, mtnSeqDirPropN, mtSeqDirPropN, qSeqDirPropN);
    qlua_reg_op2(opsN);
    qlua_reg_dot(qSeqDirPropN,  q_p_dot_N);
#endif

    return 0;
}

int
fini_seqdirprop(lua_State *L)
{
    return 0;
}
