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
#include "qdp_d2.h"
#include "qdp_d3.h"
#include "qdp_dn.h"
#include "qla_types.h"
#include "qla_d2.h"
#include "qla_d3.h"
#include "qla_dn.h"
#include <math.h>

#define COLOR_EPS 1e-10

#define QNc  '2'
#define Qcolors "2"
#define Qs(a)   a ## 2
#define Qx(a,b)  a ## 2 ## b
#define QC(x)    2
#include "latcolmat-x.c"                                            /* DEPS */

#define QNc  '3'
#define Qcolors "3"
#define Qs(a)   a ## 3
#define Qx(a,b)  a ## 3 ## b
#define QC(x)    3
#include "latcolmat-x.c"                                            /* DEPS */

#define QNc  'N'
#define Qcolors "N"
#define Qs(a)   a ## N
#define Qx(a,b)  a ## N ## b
#define QC(x)    (x)->nc
#include "latcolmat-x.c"                                            /* DEPS */

static int
do_gaussian(lua_State *L, mLatRandom *a, mLattice *S, int nc)
{
    switch (nc) {
    case 2: {
        mLatColMat2 *r = qlua_newLatColMat2(L, lua_gettop(L), 2);
        CALL_QDP(L);
        QDP_D2_M_eq_gaussian_S(r->ptr, a->ptr, *S->qss);
        return 1;
    }
    case 3: {
        mLatColMat3 *r = qlua_newLatColMat3(L, lua_gettop(L), 3);
        CALL_QDP(L);
        QDP_D3_M_eq_gaussian_S(r->ptr, a->ptr, *S->qss);
        return 1;
    }
    default: {
        mLatColMatN *r = qlua_newLatColMatN(L, lua_gettop(L), nc);
        CALL_QDP(L);
        QDP_DN_M_eq_gaussian_S(nc, r->ptr, a->ptr, *S->qss);
        return 1;
    }
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

    if ((nc < 2) || (nc > QLA_MAX_Nc))
        return luaL_error(L, "bad number of colors %d (MAX=%d)",
                          nc, QLA_MAX_Nc);

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
    case 2: return q_latcolmat_2(L, S, 2, 0);
    case 3: return q_latcolmat_3(L, S, 3, 0);
    default: return q_latcolmat_N(L, S, S->nc, 0);
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
    mLattice *S = qlua_checkLattice(L, 1);
    int nc = luaL_checkint(L, 2);
    
    switch (nc) {
    case 2: return q_latcolmat_2(L, S, 2, 1);
    case 3: return q_latcolmat_3(L, S, 3, 1);
    default: return q_latcolmat_N(L, S, nc, 1);
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
    static const QLUA_Op2 ops[] = {
        { qlua_add_table, qLatColMat2,  qLatColMat2,  q_M_add_M_2 },
        { qlua_add_table, qLatColMat3,  qLatColMat3,  q_M_add_M_3 },
        { qlua_add_table, qLatColMatN,  qLatColMatN,  q_M_add_M_N },
        { qlua_sub_table, qLatColMat2,  qLatColMat2,  q_M_sub_M_2 },
        { qlua_sub_table, qLatColMat3,  qLatColMat3,  q_M_sub_M_3 },
        { qlua_sub_table, qLatColMatN,  qLatColMatN,  q_M_sub_M_N },
        { qlua_mul_table, qReal,        qLatColMat2,  q_r_mul_M_2 },
        { qlua_mul_table, qReal,        qLatColMat3,  q_r_mul_M_3 },
        { qlua_mul_table, qReal,        qLatColMatN,  q_r_mul_M_N },
        { qlua_mul_table, qLatColMat2,  qReal,        q_M_mul_r_2 },
        { qlua_mul_table, qLatColMat3,  qReal,        q_M_mul_r_3 },
        { qlua_mul_table, qLatColMatN,  qReal,        q_M_mul_r_N },
        { qlua_mul_table, qComplex,     qLatColMat2,  q_c_mul_M_2 },
        { qlua_mul_table, qComplex,     qLatColMat3,  q_c_mul_M_3 },
        { qlua_mul_table, qComplex,     qLatColMatN,  q_c_mul_M_N },
        { qlua_mul_table, qLatColMat2,  qComplex,     q_M_mul_c_2 },
        { qlua_mul_table, qLatColMat3,  qComplex,     q_M_mul_c_3 },
        { qlua_mul_table, qLatColMatN,  qComplex,     q_M_mul_c_N },
        { qlua_mul_table, qLatReal,     qLatColMat2,  q_R_mul_M_2 },
        { qlua_mul_table, qLatReal,     qLatColMat3,  q_R_mul_M_3 },
        { qlua_mul_table, qLatReal,     qLatColMatN,  q_R_mul_M_N },
        { qlua_mul_table, qLatColMat2,  qLatReal,     q_M_mul_R_2 },
        { qlua_mul_table, qLatColMat3,  qLatReal,     q_M_mul_R_3 },
        { qlua_mul_table, qLatColMatN,  qLatReal,     q_M_mul_R_N },
        { qlua_mul_table, qLatComplex,  qLatColMat2,  q_C_mul_M_2 },
        { qlua_mul_table, qLatComplex,  qLatColMat3,  q_C_mul_M_3 },
        { qlua_mul_table, qLatComplex,  qLatColMatN,  q_C_mul_M_N },
        { qlua_mul_table, qLatColMat2,  qLatComplex,  q_M_mul_C_2 },
        { qlua_mul_table, qLatColMat3,  qLatComplex,  q_M_mul_C_3 },
        { qlua_mul_table, qLatColMatN,  qLatComplex,  q_M_mul_C_N },
        { qlua_mul_table, qLatColMat2,  qLatColVec2,  q_M_mul_V_2 },
        { qlua_mul_table, qLatColMat3,  qLatColVec3,  q_M_mul_V_3 },
        { qlua_mul_table, qLatColMatN,  qLatColVecN,  q_M_mul_V_N },
        { qlua_mul_table, qLatColMat2,  qLatColMat2,  q_M_mul_M_2 },
        { qlua_mul_table, qLatColMat3,  qLatColMat3,  q_M_mul_M_3 },
        { qlua_mul_table, qLatColMatN,  qLatColMatN,  q_M_mul_M_N },
        { qlua_div_table, qLatColMat2,  qReal,        q_M_div_r_2 },
        { qlua_div_table, qLatColMat3,  qReal,        q_M_div_r_3 },
        { qlua_div_table, qLatColMatN,  qReal,        q_M_div_r_N },
        { qlua_div_table, qLatColMat2,  qComplex,     q_M_div_c_2 },
        { qlua_div_table, qLatColMat3,  qComplex,     q_M_div_c_3 },
        { qlua_div_table, qLatColMatN,  qComplex,     q_M_div_c_N },
        { NULL,           qNoType,      qNoType,      NULL        }
    };
    luaL_getmetatable(L, opLattice);
    luaL_register(L, NULL, fLatColMat);
    lua_pop(L, 1);
    qlua_reg_op2(ops);
    qlua_reg_dot(qLatColMat2,  q_M_dot_2);
    qlua_reg_dot(qLatColMat3,  q_M_dot_3);
    qlua_reg_dot(qLatColMatN,  q_M_dot_N);

    return 0;
}

int
fini_latcolmat(lua_State *L)
{
    return 0;
}
