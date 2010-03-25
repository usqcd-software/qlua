#include "qlua.h"                                                    /* DEPS */
#include "lattice.h"                                                 /* DEPS */
#include "latreal.h"                                                 /* DEPS */
#include "latcolvec.h"                                               /* DEPS */
#include "qcomplex.h"                                                /* DEPS */
#include "latint.h"                                                  /* DEPS */
#include "latcomplex.h"                                              /* DEPS */
#include "latrandom.h"                                               /* DEPS */
#include "qmp.h"

#define QNc  '2'
#define Qcolors "2"
#define Qs(a)   a ## 2
#define Qx(a,b)  a ## 2 ## b
#define QC(x)    2
#include "latcolvec-x.c"                                            /* DEPS */

#define QNc  '3'
#define Qcolors "3"
#define Qs(a)   a ## 3
#define Qx(a,b)  a ## 3 ## b
#define QC(x)    3
#include "latcolvec-x.c"                                            /* DEPS */

#define QNc  'N'
#define Qcolors "N"
#define Qs(a)   a ## N
#define Qx(a,b)  a ## N ## b
#define QC(x)    (x)->nc
#include "latcolvec-x.c"                                            /* DEPS */

static int
do_gaussian(lua_State *L, mLatRandom *a, mLattice *S, int nc)
{
    switch (nc) {
    case 2: {
        mLatColVec2 *r = qlua_newLatColVec2(L, lua_gettop(L), 2);
        CALL_QDP(L);
        QDP_D2_V_eq_gaussian_S(r->ptr, a->ptr, *S->qss);
        return 1;
    }
    case 3: {
        mLatColVec3 *r = qlua_newLatColVec3(L, lua_gettop(L), 3);
        CALL_QDP(L);
        QDP_D3_V_eq_gaussian_S(r->ptr, a->ptr, *S->qss);
        return 1;
    }
    default: {
        mLatColVecN *r = qlua_newLatColVecN(L, lua_gettop(L), nc);
        CALL_QDP(L);
        QDP_DN_V_eq_gaussian_S(r->ptr, a->ptr, *S->qss);
        return 1;
    }
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
 *  L:ColorVector() -- zero vector in default colors
 *  L:ColorVector(C,m) -- default colors, X[m] = C
 */
static int
q_latcolvec(lua_State *L)
{
    mLattice *S = qlua_checkLattice(L, 1);

    switch (S->nc) {
    case 2: return q_latcolvec_2(L, S, 2, 0);
    case 3: return q_latcolvec_3(L, S, 3, 0);
    default: return q_latcolvec_N(L, S, S->nc, 0);
    }
}

/* Lattice color vectors:
 *  L:ColorVectorN(n) -- zero vector in n colors
 *  L:ColorVectorN(n,C,m) -- n colors, X[m] = C
 */
static int
q_latcolvecN(lua_State *L)
{
    mLattice *S = qlua_checkLattice(L, 1);
    int nc = luaL_checkint(L, 2);

    switch (nc) {
    case 2: return q_latcolvec_2(L, S, 2, 1);
    case 3: return q_latcolvec_3(L, S, 3, 1);
    default: return q_latcolvec_N(L, S, nc, 1);
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
    static const QLUA_Op2 ops[] = {
        { qlua_add_table, qLatColVec2,  qLatColVec2,  q_V_add_V_2 },
        { qlua_add_table, qLatColVec3,  qLatColVec3,  q_V_add_V_3 },
        { qlua_add_table, qLatColVecN,  qLatColVecN,  q_V_add_V_N },
        { qlua_sub_table, qLatColVec2,  qLatColVec2,  q_V_sub_V_2 },
        { qlua_sub_table, qLatColVec3,  qLatColVec3,  q_V_sub_V_3 },
        { qlua_sub_table, qLatColVecN,  qLatColVecN,  q_V_sub_V_N },
        { qlua_mul_table, qReal,        qLatColVec2,  q_r_mul_V_2 },
        { qlua_mul_table, qReal,        qLatColVec3,  q_r_mul_V_3 },
        { qlua_mul_table, qReal,        qLatColVecN,  q_r_mul_V_N },
        { qlua_mul_table, qLatColVec2,  qReal,        q_V_mul_r_2 },
        { qlua_mul_table, qLatColVec3,  qReal,        q_V_mul_r_3 },
        { qlua_mul_table, qLatColVecN,  qReal,        q_V_mul_r_N },
        { qlua_mul_table, qComplex,     qLatColVec2,  q_c_mul_V_2 },
        { qlua_mul_table, qComplex,     qLatColVec3,  q_c_mul_V_3 },
        { qlua_mul_table, qComplex,     qLatColVecN,  q_c_mul_V_N },
        { qlua_mul_table, qLatColVec2,  qComplex,     q_V_mul_c_2 },
        { qlua_mul_table, qLatColVec3,  qComplex,     q_V_mul_c_3 },
        { qlua_mul_table, qLatColVecN,  qComplex,     q_V_mul_c_N },
        { qlua_mul_table, qLatReal,     qLatColVec2,  q_R_mul_V_2 },
        { qlua_mul_table, qLatReal,     qLatColVec3,  q_R_mul_V_3 },
        { qlua_mul_table, qLatReal,     qLatColVecN,  q_R_mul_V_N },
        { qlua_mul_table, qLatColVec2,  qLatReal,     q_V_mul_R_2 },
        { qlua_mul_table, qLatColVec3,  qLatReal,     q_V_mul_R_3 },
        { qlua_mul_table, qLatColVecN,  qLatReal,     q_V_mul_R_N },
        { qlua_mul_table, qLatComplex,  qLatColVec2,  q_C_mul_V_2 },
        { qlua_mul_table, qLatComplex,  qLatColVec3,  q_C_mul_V_3 },
        { qlua_mul_table, qLatComplex,  qLatColVecN,  q_C_mul_V_N },
        { qlua_mul_table, qLatColVec2,  qLatComplex,  q_V_mul_C_2 },
        { qlua_mul_table, qLatColVec3,  qLatComplex,  q_V_mul_C_3 },
        { qlua_mul_table, qLatColVecN,  qLatComplex,  q_V_mul_C_N },
        { qlua_div_table, qLatColVec2,  qReal,        q_V_div_r_2 },
        { qlua_div_table, qLatColVec3,  qReal,        q_V_div_r_3 },
        { qlua_div_table, qLatColVecN,  qReal,        q_V_div_r_N },
        { qlua_div_table, qLatColVec2,  qComplex,     q_V_div_c_2 },
        { qlua_div_table, qLatColVec3,  qComplex,     q_V_div_c_3 },
        { qlua_div_table, qLatColVecN,  qComplex,     q_V_div_c_N },
        { NULL,           qNoType,      qNoType,      NULL        }
    };
    luaL_getmetatable(L, opLattice);
    luaL_register(L, NULL, fLatColVec);
    lua_pop(L, 1);
    qlua_reg_op2(ops);
    qlua_reg_dot(qLatColVec2,  q_V_dot_2);
    qlua_reg_dot(qLatColVec3,  q_V_dot_3);
    qlua_reg_dot(qLatColVecN,  q_V_dot_N);

    return 0;
}

int
fini_latcolvec(lua_State *L)
{
    return 0;
}
