#include "qlua.h"                                                    /* DEPS */
#include "qcomplex.h"                                                /* DEPS */
#include "lattice.h"                                                 /* DEPS */
#include "latrandom.h"                                               /* DEPS */
#include "latreal.h"                                                 /* DEPS */
#include "latcomplex.h"                                              /* DEPS */
#include "latcolmat.h"                                               /* DEPS */
#include "latdirferm.h"                                              /* DEPS */
#include "latdirprop.h"                                              /* DEPS */
#include "latint.h"                                                  /* DEPS */
#include "qmp.h"
#include "qdp_d2.h"
#include "qdp_d3.h"
#include "qdp_dn.h"
#include "qla_types.h"
#include "qla_d2.h"
#include "qla_d3.h"
#include "qla_dn.h"

#define QNc  '2'
#define Qcolors "2"
#define Qs(a)   a ## 2
#define Qx(a,b)  a ## 2 ## b
#define QC(x)    2
#include "latdirprop-x.c"                                           /* DEPS */

#define QNc  '3'
#define Qcolors "3"
#define Qs(a)   a ## 3
#define Qx(a,b)  a ## 3 ## b
#define QC(x)    3
#include "latdirprop-x.c"                                           /* DEPS */

#define QNc  'N'
#define Qcolors "N"
#define Qs(a)   a ## N
#define Qx(a,b)  a ## N ## b
#define QC(x)    (x)->nc
#include "latdirprop-x.c"                                           /* DEPS */

static int
do_gaussian(lua_State *L, mLatRandom *a, mLattice *S, int nc)
{
    switch (nc) {
    case 2: {
        mLatDirProp2 *r = qlua_newLatDirProp2(L, lua_gettop(L), 2);
        CALL_QDP(L);
        QDP_D2_P_eq_gaussian_S(r->ptr, a->ptr, *S->qss);
        return 1;
    }
    case 3: {
        mLatDirProp3 *r = qlua_newLatDirProp3(L, lua_gettop(L), 3);
        CALL_QDP(L);
        QDP_D3_P_eq_gaussian_S(r->ptr, a->ptr, *S->qss);
        return 1;
    }
    default: {
        mLatDirPropN *r = qlua_newLatDirPropN(L, lua_gettop(L), nc);
        CALL_QDP(L);
        QDP_DN_P_eq_gaussian_S(r->ptr, a->ptr, *S->qss);
        return 1;
    }
    }
}

/* Random vectors of default colors
 * S:gaussian_DiracPropagator()
 */
int
q_P_gaussian(lua_State *L)
{
    mLatRandom *a = qlua_checkLatRandom(L, 1, NULL);
    mLattice *S = qlua_ObjLattice(L, 1);

    return do_gaussian(L, a, S, S->nc);
}

/* Random vectors of NC colors
 * S:gaussian_DiracPropagatorN(nc)
 */
int
q_P_gaussian_N(lua_State *L)
{
    mLatRandom *a = qlua_checkLatRandom(L, 1, NULL);
    mLattice *S = qlua_ObjLattice(L, 1);
    int nc = luaL_checkint(L, 2);

    if (nc < 2)
        return luaL_error(L, "bad number of colors");

    return do_gaussian(L, a, S, nc);
}

/* Lattice Dirac Propagators:
 *  L:DiracPropagator()             -- zero propagator in default colors
 *  L:DiracPropagator(M)            -- default colors, P[.,k,.,k] = M for all k
 *  L:DiracPropagator(D,{c=i,d=j})  -- default colors, P[i,j,.,.] = D
 */
static int
q_latdirprop(lua_State *L)
{
    mLattice *S = qlua_checkLattice(L, 1);

    switch (S->nc) {
    case 2: return q_latdirprop_2(L, S, 2, 0);
    case 3: return q_latdirprop_3(L, S, 3, 0);
    default: return q_latdirprop_N(L, S, S->nc, 0);
    }
}

/* Lattice Dirac Propagators:
 *  L:DiracPropagatorN(n)              -- zero propagator in n colors
 *  L:DiracPropagatorN(n,M)            -- n colors, P[.,k,.,k] = M for all k
 *  L:DiracPropagatorN(n,D,{c=i,d=j})  -- n colors, P[i,j,.,.] = D
 */
static int
q_latdirpropN(lua_State *L)
{
    mLattice *S = qlua_checkLattice(L, 1);
    int nc = luaL_checkint(L, 2);

    switch (nc) {
    case 2: return q_latdirprop_2(L, S, 2, 1);
    case 3: return q_latdirprop_3(L, S, 3, 1);
    default: return q_latdirprop_N(L, S, nc, 1);
    }
}

static struct {
    QLA_D3_DiracPropagator *a;
    QLA_D3_DiracPropagator *b;
} QQ_args;

/* TODO force to unroll loops over color&spin indices? */
#define do_QQ_contract_func(contract_idx,A,B,C,D)\
 static void do_QQc ## contract_idx(QLA_D3_DiracPropagator *q3, int idx) \
 {                                                                      \
     QLA_D3_DiracPropagator *q1 = &QQ_args.a[idx];                      \
     QLA_D3_DiracPropagator *q2 = &QQ_args.b[idx];                      \
     static const int eps[3][3] = { { 0, 1, 2}, { 1, 2, 0}, { 2, 0, 1} }; \
     int p_a, p_b;                                                      \
     for (p_a = 0; p_a < 3; p_a++) { /* Nc */                           \
         int i1 = eps[p_a][0], j1 = eps[p_a][1], k1 = eps[p_a][2];      \
         for (p_b = 0; p_b < 3; p_b++) { /* Nc */                       \
             int i2 = eps[p_b][0], j2 = eps[p_b][1], k2 = eps[p_b][2];  \
             int a, b, c;                                               \
             for (a = 0; a < QDP_Ns; a++) {                             \
                 for (b = 0; b < QDP_Ns; b++) {                         \
                     QLA_Complex s3;                                    \
                     QLA_c_eq_r(s3, 0.0);                               \
                     for (c = 0; c < QDP_Ns; c++) {                     \
                         QLA_c_peq_c_times_c(s3,                        \
                                             QLA_elem_P(*q1,i1,(A),i2,(B)), \
                                             QLA_elem_P(*q2,j1,(C),j2,(D))); \
                         QLA_c_meq_c_times_c(s3,                        \
                                             QLA_elem_P(*q1,i1,(A),j2,(B)), \
                                             QLA_elem_P(*q2,j1,(C),i2,(D))); \
                         QLA_c_meq_c_times_c(s3,                        \
                                             QLA_elem_P(*q1,j1,(A),i2,(B)), \
                                             QLA_elem_P(*q2,i1,(C),j2,(D))); \
                         QLA_c_peq_c_times_c(s3,                        \
                                             QLA_elem_P(*q1,j1,(A),j2,(B)), \
                                             QLA_elem_P(*q2,i1,(C),i2,(D))); \
                     }                                                  \
                     QLA_c_eq_c(QLA_elem_P(*q3,k2,a,k1,b), s3);         \
                 }                                                      \
             }                                                          \
         }                                                              \
     }                                                                  \
 }
do_QQ_contract_func(12, c,c,a,b);
do_QQ_contract_func(13, c,a,c,b);
do_QQ_contract_func(14, c,a,b,c);
do_QQ_contract_func(23, a,c,c,b);
do_QQ_contract_func(24, a,c,b,c);
do_QQ_contract_func(34, a,b,c,c);
#undef do_QQ_contract_func

static int
q_su3contract_general(lua_State *L, 
        void (*do_contract)(QLA_D3_DiracPropagator *, int))
{
    mLatDirProp3 *a = qlua_checkLatDirProp3(L, 1, NULL, 3); /* the 1st arg */
    mLattice *S = qlua_ObjLattice(L, 1);
    mLatDirProp3 *b = qlua_checkLatDirProp3(L, 2, S, 3); /* the 2nd arg */
    mLatDirProp3 *r = qlua_newLatDirProp3(L, lua_gettop(L), 3); 

    CALL_QDP(L);
    QQ_args.a = QDP_D3_expose_P(a->ptr);
    QQ_args.b = QDP_D3_expose_P(b->ptr);
    QDP_D3_P_eq_funci(r->ptr, do_contract, *S->qss);
    QDP_D3_reset_P(a->ptr);
    QDP_D3_reset_P(b->ptr);
    QQ_args.a = 0;
    QQ_args.b = 0;

    return 1;
}
#define q_su3contract_func(contract_idx) \
static int q_su3contract ## contract_idx(lua_State *L) \
{ \
    return q_su3contract_general(L, do_QQc ## contract_idx); \
}
q_su3contract_func(12);
q_su3contract_func(13);
q_su3contract_func(14);
q_su3contract_func(23);
q_su3contract_func(24);
q_su3contract_func(34);
#undef q_su3contract_func

static struct luaL_Reg fLatDirProp[] = {
    { "DiracPropagator",   q_latdirprop  },
    { "DiracPropagatorN",  q_latdirpropN },
    { NULL,                NULL          }
};

static struct luaL_Reg fQCDDirProp[] = {
    { "quarkContract12",   q_su3contract12 },
    { "quarkContract13",   q_su3contract13 },
    { "quarkContract14",   q_su3contract14 },
    { "quarkContract23",   q_su3contract23 },
    { "quarkContract24",   q_su3contract24 },
    { "quarkContract34",   q_su3contract34 },
    { NULL,                NULL            }
};

int
init_latdirprop(lua_State *L)
{
    static const QLUA_Op2 ops[] = {
        { qlua_add_table, qLatDirProp2,  qLatDirProp2,  q_P_add_P_2 },
        { qlua_add_table, qLatDirProp3,  qLatDirProp3,  q_P_add_P_3 },
        { qlua_add_table, qLatDirPropN,  qLatDirPropN,  q_P_add_P_N },
        { qlua_sub_table, qLatDirProp2,  qLatDirProp2,  q_P_sub_P_2 },
        { qlua_sub_table, qLatDirProp3,  qLatDirProp3,  q_P_sub_P_3 },
        { qlua_sub_table, qLatDirPropN,  qLatDirPropN,  q_P_sub_P_N },
        { qlua_mul_table, qReal,         qLatDirProp2,  q_r_mul_P_2 },
        { qlua_mul_table, qReal,         qLatDirProp3,  q_r_mul_P_3 },
        { qlua_mul_table, qReal,         qLatDirPropN,  q_r_mul_P_N },
        { qlua_mul_table, qLatDirProp2,  qReal,         q_P_mul_r_2 },
        { qlua_mul_table, qLatDirProp3,  qReal,         q_P_mul_r_3 },
        { qlua_mul_table, qLatDirPropN,  qReal,         q_P_mul_r_N },
        { qlua_mul_table, qComplex,      qLatDirProp2,  q_c_mul_P_2 },
        { qlua_mul_table, qComplex,      qLatDirProp3,  q_c_mul_P_3 },
        { qlua_mul_table, qComplex,      qLatDirPropN,  q_c_mul_P_N },
        { qlua_mul_table, qLatDirProp2,  qComplex,      q_P_mul_c_2 },
        { qlua_mul_table, qLatDirProp3,  qComplex,      q_P_mul_c_3 },
        { qlua_mul_table, qLatDirPropN,  qComplex,      q_P_mul_c_N },
        { qlua_mul_table, qLatReal,      qLatDirProp2,  q_R_mul_P_2 },
        { qlua_mul_table, qLatReal,      qLatDirProp3,  q_R_mul_P_3 },
        { qlua_mul_table, qLatReal,      qLatDirPropN,  q_R_mul_P_N },
        { qlua_mul_table, qLatDirProp2,  qLatReal,      q_P_mul_R_2 },
        { qlua_mul_table, qLatDirProp3,  qLatReal,      q_P_mul_R_3 },
        { qlua_mul_table, qLatDirPropN,  qLatReal,      q_P_mul_R_N },
        { qlua_mul_table, qLatComplex,   qLatDirProp2,  q_C_mul_P_2 },
        { qlua_mul_table, qLatComplex,   qLatDirProp3,  q_C_mul_P_3 },
        { qlua_mul_table, qLatComplex,   qLatDirPropN,  q_C_mul_P_N },
        { qlua_mul_table, qLatDirProp2,  qLatComplex,   q_P_mul_C_2 },
        { qlua_mul_table, qLatDirProp3,  qLatComplex,   q_P_mul_C_3 },
        { qlua_mul_table, qLatDirPropN,  qLatComplex,   q_P_mul_C_N },
        { qlua_mul_table, qLatDirProp2,  qLatDirProp2,  q_P_mul_P_2 },
        { qlua_mul_table, qLatDirProp3,  qLatDirProp3,  q_P_mul_P_3 },
        { qlua_mul_table, qLatDirPropN,  qLatDirPropN,  q_P_mul_P_N },
        { qlua_mul_table, qLatDirProp2,  qLatColMat2,   q_P_mul_M_2 },
        { qlua_mul_table, qLatDirProp3,  qLatColMat3,   q_P_mul_M_3 },
        { qlua_mul_table, qLatDirPropN,  qLatColMatN,   q_P_mul_M_N },
        { qlua_mul_table, qLatColMat2,   qLatDirProp2,  q_M_mul_P_2 },
        { qlua_mul_table, qLatColMat3,   qLatDirProp3,  q_M_mul_P_3 },
        { qlua_mul_table, qLatColMatN,   qLatDirPropN,  q_M_mul_P_N },
        { qlua_div_table, qLatDirProp2,  qReal,         q_P_div_r_2 },
        { qlua_div_table, qLatDirProp3,  qReal,         q_P_div_r_3 },
        { qlua_div_table, qLatDirPropN,  qReal,         q_P_div_r_N },
        { qlua_div_table, qLatDirProp2,  qComplex,      q_P_div_c_2 },
        { qlua_div_table, qLatDirProp3,  qComplex,      q_P_div_c_3 },
        { qlua_div_table, qLatDirPropN,  qComplex,      q_P_div_c_N },
        { NULL,           qNoType,       qNoType,       NULL        }
    };
    luaL_register(L, qcdlib, fQCDDirProp);
    luaL_getmetatable(L, opLattice);
    luaL_register(L, NULL, fLatDirProp);
    lua_pop(L, 1);
    qlua_reg_op2(ops);
    qlua_reg_dot(qLatDirProp2,  q_P_dot_2);
    qlua_reg_dot(qLatDirProp3,  q_P_dot_3);
    qlua_reg_dot(qLatDirPropN,  q_P_dot_N);

    return 0;
}

int
fini_latdirprop(lua_State *L)
{
    return 0;
}
