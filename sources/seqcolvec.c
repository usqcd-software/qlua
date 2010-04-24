#include "modules.h"                                                 /* DEPS */
#include "qlua.h"                                                    /* DEPS */
#include "seqcolvec.h"                                               /* DEPS */
#include "qcomplex.h"                                                /* DEPS */
#include "seqrandom.h"                                               /* DEPS */
#include "qla.h"

#if USE_Nc2
#define QNc  '2'
#define Qcolors "2"
#define Qs(a)   a ## 2
#define Qx(a,b)  a ## 2 ## b
#define QC(x)    2
#define QNC(x)
#include "seqcolvec-x.c"                                            /* DEPS */
#endif

#if USE_Nc3
#define QNc  '3'
#define Qcolors "3"
#define Qs(a)   a ## 3
#define Qx(a,b)  a ## 3 ## b
#define QC(x)    3
#define QNC(x)
#include "seqcolvec-x.c"                                            /* DEPS */
#endif

#if USE_NcN
#define QNc  'N'
#define Qcolors "N"
#define Qs(a)   a ## N
#define Qx(a,b)  a ## N ## b
#define QC(x)    (x)->nc
#define QNC(x)   (x),
#include "seqcolvec-x.c"                                            /* DEPS */
#endif

/* Random vectors of NC colors
 * s:gaussian_ColorVectorN(nc)
 */
int
q_v_gaussian_N(lua_State *L)
{
    mSeqRandom *a = qlua_checkSeqRandom(L, 1);
    int nc = luaL_checkint(L, 2);

    if (nc < 2)
        return luaL_error(L, "bad number of colors");

    switch (nc) {
#if USE_Nc2
    case 2: {
        mSeqColVec2 *r = qlua_newSeqColVec2(L, 2);
        QLA_D2_V_eq_gaussian_S(r->ptr, &a->state);
        return 1;
    }
#endif
#if USE_Nc3
    case 3: {
        mSeqColVec3 *r = qlua_newSeqColVec3(L, 3);
        QLA_D3_V_eq_gaussian_S(r->ptr, &a->state);
        return 1;
    }
#endif
#if USE_NcN
    default: {
        mSeqColVecN *r = qlua_newSeqColVecN(L, nc);
        QLA_DN_V_eq_gaussian_S(nc, r->ptr, &a->state);
        return 1;
    }
#else
    default:
        return luaL_error(L, "Bad number of colors in gaussian");
#endif
    }
}

/* Sequential color vectors:
 *  qcd.ColorVectorN(n) -- zero vector in n colors
 *  qcd.ColorVectorN(n,C,m) -- n colors, X[m] = C
 */
static int
q_seqcolvecN(lua_State *L)
{
    int nc = luaL_checkint(L, 1);

    switch (nc) {
#if USE_Nc2
    case 2: return q_seqcolvec_2(L, 2);
#endif
#if USE_Nc3
    case 3: return q_seqcolvec_3(L, 3);
#endif
#if USE_NcN
    default: return q_seqcolvec_N(L, nc);
#else
    default: return luaL_error(L, "bad number of colors");
#endif
    }
}

static struct luaL_Reg fSeqColVec[] = {
    { "ColorVectorN",      q_seqcolvecN  },
    { NULL,                NULL          }
};

int
init_seqcolvec(lua_State *L)
{
    luaL_register(L, qcdlib, fSeqColVec);
#if USE_Nc2
    qlua_metatable(L, mtnSeqColVec2, mtSeqColVec2, qSeqColVec2);
    qlua_reg_op2(ops2);
    qlua_reg_dot(qSeqColVec2,  q_v_dot_2);
#endif
#if USE_Nc3
    qlua_metatable(L, mtnSeqColVec3, mtSeqColVec3, qSeqColVec3);
    qlua_reg_op2(ops3);
    qlua_reg_dot(qSeqColVec3,  q_v_dot_3);
#endif
#if USE_NcN
    qlua_metatable(L, mtnSeqColVecN, mtSeqColVecN, qSeqColVecN);
    qlua_reg_op2(opsN);
    qlua_reg_dot(qSeqColVecN,  q_v_dot_N);
#endif

    return 0;
}

int
fini_seqcolvec(lua_State *L)
{
    return 0;
}
