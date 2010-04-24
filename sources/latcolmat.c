#include "modules.h"                                                 /* DEPS */
#include "qlua.h"                                                    /* DEPS */
#include "qcomplex.h"                                                /* DEPS */
#include "lattice.h"                                                 /* DEPS */
#include "latint.h"                                                  /* DEPS */
#include "latreal.h"                                                 /* DEPS */
#include "latcomplex.h"                                              /* DEPS */
#include "latcolmat.h"                                               /* DEPS */
#include "latrandom.h"                                               /* DEPS */
#include "latcolvec.h"                                               /* DEPS */
#include "qmp.h"

#include <math.h>

#define COLOR_EPS 1e-10

#if USE_Nc2
#define QNc  '2'
#define Qcolors "2"
#define Qs(a)   a ## 2
#define Qx(a,b)  a ## 2 ## b
#define QC(x)    2
#include "latcolmat-x.c"                                            /* DEPS */
#endif

#if USE_Nc3
#define QNc  '3'
#define Qcolors "3"
#define Qs(a)   a ## 3
#define Qx(a,b)  a ## 3 ## b
#define QC(x)    3
#include "latcolmat-x.c"                                            /* DEPS */
#endif

#if USE_NcN
#define QNc  'N'
#define Qcolors "N"
#define Qs(a)   a ## N
#define Qx(a,b)  a ## N ## b
#define QC(x)    (x)->nc
#include "latcolmat-x.c"                                            /* DEPS */
#endif

static int
do_gaussian(lua_State *L, mLatRandom *a, mLattice *S, int nc)
{
    switch (nc) {
#if USE_Nc2
    case 2: {
        mLatColMat2 *r = qlua_newLatColMat2(L, lua_gettop(L), 2);
        CALL_QDP(L);
        QDP_D2_M_eq_gaussian_S(r->ptr, a->ptr, *S->qss);
        return 1;
    }
#endif
#if USE_Nc3
    case 3: {
        mLatColMat3 *r = qlua_newLatColMat3(L, lua_gettop(L), 3);
        CALL_QDP(L);
        QDP_D3_M_eq_gaussian_S(r->ptr, a->ptr, *S->qss);
        return 1;
    }
#endif
#if USE_NcN
    default: {
        mLatColMatN *r = qlua_newLatColMatN(L, lua_gettop(L), nc);
        CALL_QDP(L);
        QDP_DN_M_eq_gaussian_S(r->ptr, a->ptr, *S->qss);
        return 1;
    }
#else
    default:
        return luaL_error(L, "bad number of colors");
#endif
    }
}

/* Random vectors of default colors
 * S:gaussian_ColorVector()
 */
int
q_M_gaussian(lua_State *L)
{
    mLatRandom *a = qlua_checkLatRandom(L, 1, NULL);
    mLattice *S = qlua_ObjLattice(L, 1);

    return do_gaussian(L, a, S, S->nc);
}

/* Random vectors of NC colors
 * S:gaussian_ColorVectorN(nc)
 */
int
q_M_gaussian_N(lua_State *L)
{
    mLatRandom *a = qlua_checkLatRandom(L, 1, NULL);
    mLattice *S = qlua_ObjLattice(L, 1);
    int nc = luaL_checkint(L, 2);

    if (nc < 2)
        return luaL_error(L, "bad number of colors");

    return do_gaussian(L, a, S, nc);
}

/* Lattice color matrices
 * L:ColorMatrix()            -- zero matrix in default colors
 * L:ColorMatrix(r)           -- r * unit matrix in default colors
 * L:ColorMatrix(c)           -- c * unit matrix in default colors
 * L:ColorMatrix(C,{a=i,b=j}) -- default colors M[i,j] = C, zero otherwise
 * L:ColorMatrix(V,{b=j})     -- default colors M[*,j] = V, zero otherwise
 * L:ColorMatrix(V1,V2)       -- default colors M[i,j] = V1[i] * V2[j]:conj
 */
static int
q_latcolmat(lua_State *L)
{
    mLattice *S = qlua_checkLattice(L, 1);
    
    switch (S->nc) {
#if USE_Nc2
    case 2: return q_latcolmat_2(L, S, 2, 0);
#endif
#if USE_Nc3
    case 3: return q_latcolmat_3(L, S, 3, 0);
#endif
#if USE_NcN
    default: return q_latcolmat_N(L, S, S->nc, 0);
#else
    default: return luaL_error(L, "bad number of colors");
#endif
    }
}

/* Lattice color matrices
 * L:ColorMatrixN(n)            -- zero matrix in n colors
 * L:ColorMatrixN(n,r)           -- r * unit matrix in n colors
 * L:ColorMatrixN(n,c)           -- c * unit matrix in n colors
 * L:ColorMatrixN(n,C,{a=i,b=j}) -- n colors M[i,j] = C, zero otherwise
 * L:ColorMatrixN(n,V,{b=j})     -- n colors M[*,j] = V, zero otherwise
 * L:ColorMatrixN(n,V1,V2)       -- n colors M[i,j] = V1[i] * V2[j]:conj
 */
static int
q_latcolmatN(lua_State *L)
{
#if USE_Nc2 || USE_Nc3 || USE_NcN
    mLattice *S = qlua_checkLattice(L, 1);
#endif
    int nc = luaL_checkint(L, 2);
    
    switch (nc) {
#if USE_Nc2
    case 2: return q_latcolmat_2(L, S, 2, 1);
#endif
#if USE_Nc3
    case 3: return q_latcolmat_3(L, S, 3, 1);
#endif
#if USE_NcN
    default: return q_latcolmat_N(L, S, nc, 1);
#else
    default: return luaL_error(L, "bad number of colors");
#endif
    }
}

static struct luaL_Reg fLatColMat[] = {
    { "ColorMatrix",       q_latcolmat  },
    { "ColorMatrixN",      q_latcolmatN },
    { NULL,                NULL         }
};

int
init_latcolmat(lua_State *L)
{
    luaL_getmetatable(L, opLattice);
    luaL_register(L, NULL, fLatColMat);
    lua_pop(L, 1);
#if USE_Nc2
    qlua_reg_op2(ops2);
    qlua_reg_dot(qLatColMat2,  q_M_dot_2);
#endif
#if USE_Nc3
    qlua_reg_op2(ops3);
    qlua_reg_dot(qLatColMat3,  q_M_dot_3);
#endif
#if USE_NcN
    qlua_reg_op2(opsN);
    qlua_reg_dot(qLatColMatN,  q_M_dot_N);
#endif

    return 0;
}

int
fini_latcolmat(lua_State *L)
{
    return 0;
}
