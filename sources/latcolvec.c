#include "modules.h"                                                 /* DEPS */
#include "qlua.h"                                                    /* DEPS */
#include "lattice.h"                                                 /* DEPS */
#include "latreal.h"                                                 /* DEPS */
#include "latcolvec.h"                                               /* DEPS */
#include "seqcolvec.h"                                               /* DEPS */
#include "qcomplex.h"                                                /* DEPS */
#include "latint.h"                                                  /* DEPS */
#include "latcomplex.h"                                              /* DEPS */
#include "latrandom.h"                                               /* DEPS */
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
#include "latcolvec-x.c"                                            /* DEPS */
#endif

#if USE_Nc3
#define QNc  '3'
#define Qcolors "3"
#define Qs(a)   a ## 3
#define Qx(a,b)  a ## 3 ## b
#define QC(x)    3
#define QNC(x)
#include "latcolvec-x.c"                                            /* DEPS */
#endif

#if USE_NcN
#define QNc  'N'
#define Qcolors "N"
#define Qs(a)   a ## N
#define Qx(a,b)  a ## N ## b
#define QC(x)    (x)->nc
#define QNC(x)   (x), 
#include "latcolvec-x.c"                                            /* DEPS */
#endif

static int
do_gaussian(lua_State *L, mLatRandom *a, mLattice *S, int nc)
{
    switch (nc) {
#if USE_Nc2
    case 2: {
        mLatColVec2 *r = qlua_newLatColVec2(L, lua_gettop(L), 2);
        CALL_QDP(L);
        QDP_D2_V_eq_gaussian_S(r->ptr, a->ptr, *S->qss);
        return 1;
    }
#endif
#if USE_Nc3
    case 3: {
        mLatColVec3 *r = qlua_newLatColVec3(L, lua_gettop(L), 3);
        CALL_QDP(L);
        QDP_D3_V_eq_gaussian_S(r->ptr, a->ptr, *S->qss);
        return 1;
    }
#endif
#if USE_NcN
    default: {
        mLatColVecN *r = qlua_newLatColVecN(L, lua_gettop(L), nc);
        CALL_QDP(L);
        QDP_DN_V_eq_gaussian_S(r->ptr, a->ptr, *S->qss);
        return 1;
    }
#else
    default:
        return luaL_error(L, "Bad number of colors in gaussian");
#endif
    }
}

/* Random vectors of default colors
 * S:gaussian_ColorVector()
 */
int
q_V_gaussian(lua_State *L)
{
    mLatRandom *a = qlua_checkLatRandom(L, 1, NULL);
    mLattice *S = qlua_ObjLattice(L, 1);

    return do_gaussian(L, a, S, S->nc);
}

/* Random vectors of NC colors
 * S:gaussian_ColorVectorN(nc)
 */
int
q_V_gaussian_N(lua_State *L)
{
    mLatRandom *a = qlua_checkLatRandom(L, 1, NULL);
    mLattice *S = qlua_ObjLattice(L, 1);
    int nc = luaL_checkint(L, 2);

    if (nc < 2)
        return luaL_error(L, "bad number of colors");

    return do_gaussian(L, a, S, nc);
}

/* Lattice color vectors:
 *  L:ColorVector()    -- zero vector in default colors
 *  L:ColorVector(v)   -- constant vector fill
 *  L:ColorVector(V)   -- lattice vector copy
 *  L:ColorVector(C,m) -- default colors, X[m] = C
 */
static int
q_latcolvec(lua_State *L)
{
    mLattice *S = qlua_checkLattice(L, 1);

    if (lua_gettop(L) == 2) {
        switch (qlua_qtype(L, 2)) {
#if USE_Nc2
        case qSeqColVec2: return q_latcolvec_seq_2(L, S, 2);
#endif      
#if USE_Nc3
        case qSeqColVec3: return q_latcolvec_seq_3(L, S, 3);
#endif      
#if USE_NcN
        case qSeqColVecN: {
            mSeqColVecN *v = qlua_checkSeqColVecN(L, 2, -1);
            return q_latcolvec_seq_N(L, S, v->nc);
        }
#endif   
        default:
            break;
        }
    }

    switch (S->nc) {
#if USE_Nc2
    case 2: return q_latcolvec_2(L, S, 2, 0);
#endif
#if USE_Nc3
    case 3: return q_latcolvec_3(L, S, 3, 0);
#endif
#if USE_NcN
    default: return q_latcolvec_N(L, S, S->nc, 0);
#else
    default: return luaL_error(L, "bad number of colors");
#endif
    }
}

/* Lattice color vectors:
 *  L:ColorVectorN(n) -- zero vector in n colors
 *  L:ColorVectorN(n,C,m) -- n colors, X[m] = C
 */
static int
q_latcolvecN(lua_State *L)
{
#if USE_Nc2 || USE_Nc3 || USE_NcN
    mLattice *S = qlua_checkLattice(L, 1);
#else
    qlua_checkLattice(L, 1);
#endif
    int nc = luaL_checkint(L, 2);

    switch (nc) {
#if USE_Nc2
    case 2: return q_latcolvec_2(L, S, 2, 1);
#endif
#if USE_Nc3
    case 3: return q_latcolvec_3(L, S, 3, 1);
#endif
#if USE_NcN
    default: return q_latcolvec_N(L, S, nc, 1);
#else
    default: return luaL_error(L, "bad number of colors");
#endif
    }
}

static struct luaL_Reg fLatColVec[] = {
    { "ColorVector",       q_latcolvec   },
    { "ColorVectorN",      q_latcolvecN  },
    { NULL,                NULL          }
};

int
init_latcolvec(lua_State *L)
{
    luaL_getmetatable(L, opLattice);
    luaL_register(L, NULL, fLatColVec);
    lua_pop(L, 1);
#if USE_Nc2
    qlua_reg_op2(ops2);
    qlua_reg_dot(qLatColVec2,  q_V_dot_2);
#endif
#if USE_Nc3
    qlua_reg_op2(ops3);
    qlua_reg_dot(qLatColVec3,  q_V_dot_3);
#endif
#if USE_NcN
    qlua_reg_op2(opsN);
    qlua_reg_dot(qLatColVecN,  q_V_dot_N);
#endif

    return 0;
}

int
fini_latcolvec(lua_State *L)
{
    return 0;
}
