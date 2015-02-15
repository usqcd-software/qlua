#include "modules.h"                                                 /* DEPS */
#include "qlua.h"                                                    /* DEPS */
#include "qcomplex.h"                                                /* DEPS */
#include "seqcolmat.h"                                               /* DEPS */
#include "seqrandom.h"                                               /* DEPS */
#include "seqcolvec.h"                                               /* DEPS */
#include "qla.h"

#define COLOR_EPS 1e-10

#if USE_Nc2
#define QNc  '2'
#define Qcolors "2"
#define Qs(a)   a ## 2
#define Qx(a,b)  a ## 2 ## b
#define QC(x)    2
#define QNC(x)
#include "seqcolmat-x.c"                                            /* DEPS */
#endif

#if USE_Nc3
#define QNc  '3'
#define Qcolors "3"
#define Qs(a)   a ## 3
#define Qx(a,b)  a ## 3 ## b
#define QC(x)    3
#define QNC(x)
#include "seqcolmat-x.c"                                            /* DEPS */
#endif

#if USE_NcN
#define QNc  'N'
#define Qcolors "N"
#define Qs(a)   a ## N
#define Qx(a,b)  a ## N ## b
#define QC(x)    (x)->nc
#define QNC(x)   (x),
#include "seqcolmat-x.c"                                            /* DEPS */
#endif

/* Random vectors of NC colors
 * S:gaussian_ColorVectorN(nc)
 */
int
q_m_gaussian_N(lua_State *L)
{
    mSeqRandom *a = qlua_checkSeqRandom(L, 1);
    int nc = luaL_checkint(L, 2);

    switch (nc) {
#if USE_Nc2
    case 2: {
        mSeqColMat2 *r = qlua_newSeqColMat2(L, 2);
        QLA_D2_M_eq_gaussian_S(r->ptr, &a->state);
        return 1;
    }
#endif
#if USE_Nc3
    case 3: {
        mSeqColMat3 *r = qlua_newSeqColMat3(L, 3);
        QLA_D3_M_eq_gaussian_S(r->ptr, &a->state);
        return 1;
    }
#endif
#if USE_NcN
    default: {
        mSeqColMatN *r = qlua_newSeqColMatN(L, nc);
        QLA_DN_M_eq_gaussian_S(nc, r->ptr, &a->state);
        return 1;
    }
#else
    default:
        return luaL_error(L, "bad number of colors");
#endif
    }
}

/* Sequential color matrices
 * qcd.ColorMatrixN(n)            -- zero matrix in n colors
 * qcd.ColorMatrixN(n,r)           -- r * unit matrix in n colors
 * qcd.ColorMatrixN(n,c)           -- c * unit matrix in n colors
 * qcd.ColorMatrixN(n,c,{a=i,b=j}) -- n colors M[i,j] = C, zero otherwise
 * qcd.ColorMatrixN(n,v,{b=j})     -- n colors M[*,j] = V, zero otherwise
 * qcd.ColorMatrixN(n,v1,v2)       -- n colors M[i,j] = V1[i] * V2[j]:conj
 */
static int
q_seqcolmatN(lua_State *L)
{
    int nc = luaL_checkint(L, 1);
    
    switch (nc) {
#if USE_Nc2
    case 2: return q_seqcolmat_2(L, 2);
#endif
#if USE_Nc3
    case 3: return q_seqcolmat_3(L, 3);
#endif
#if USE_NcN
    default: return q_seqcolmat_N(L, nc);
#else
    default: return luaL_error(L, "bad number of colors");
#endif
    }
}

static struct luaL_Reg fSeqColMat[] = {
    { "ColorMatrixN",      q_seqcolmatN },
    { NULL,                NULL         }
};

int
init_seqcolmat(lua_State *L)
{
    luaL_register(L, qcdlib, fSeqColMat);
#if USE_Nc2
    qlua_metatable(L, mtnSeqColMat2, mtSeqColMat2, qSeqColMat2);
    qlua_reg_op2(ops2);
    qlua_reg_dot(qSeqColMat2,  q_m_dot_2);
#endif
#if USE_Nc3
    qlua_metatable(L, mtnSeqColMat3, mtSeqColMat3, qSeqColMat3);
    qlua_reg_op2(ops3);
    qlua_reg_dot(qSeqColMat3,  q_m_dot_3);
#endif
#if USE_NcN
    qlua_metatable(L, mtnSeqColMatN, mtSeqColMatN, qSeqColMatN);
    qlua_reg_op2(opsN);
    qlua_reg_dot(qSeqColMatN,  q_m_dot_N);
#endif

    return 0;
}

void
fini_seqcolmat(void)
{
}
