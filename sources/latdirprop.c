#include "modules.h"                                                 /* DEPS */
#include "qlua.h"                                                    /* DEPS */
#include "qcomplex.h"                                                /* DEPS */
#include "lattice.h"                                                 /* DEPS */
#include "latrandom.h"                                               /* DEPS */
#include "latreal.h"                                                 /* DEPS */
#include "latcomplex.h"                                              /* DEPS */
#include "latcolmat.h"                                               /* DEPS */
#include "latdirferm.h"                                              /* DEPS */
#include "latdirprop.h"                                              /* DEPS */
#include "seqdirprop.h"                                              /* DEPS */
#include "latint.h"                                                  /* DEPS */
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
#include "latdirprop-x.c"                                           /* DEPS */
#endif

#if USE_Nc3
#define QNc  '3'
#define Qcolors "3"
#define Qs(a)   a ## 3
#define Qx(a,b)  a ## 3 ## b
#define QC(x)    3
#define QNC(x)
#include "latdirprop-x.c"                                           /* DEPS */
#endif

#if USE_NcN
#define QNc  'N'
#define Qcolors "N"
#define Qs(a)   a ## N
#define Qx(a,b)  a ## N ## b
#define QC(x)    (x)->nc
#define QNC(x)   (x),
#include "latdirprop-x.c"                                           /* DEPS */
#endif

static int
do_gaussian(lua_State *L, mLatRandom *a, mLattice *S, int nc)
{
    switch (nc) {
#if USE_Nc2
    case 2: {
        mLatDirProp2 *r = qlua_newLatDirProp2(L, lua_gettop(L), 2);
        CALL_QDP(L);
        QDP_D2_P_eq_gaussian_S(r->ptr, a->ptr, *S->qss);
        return 1;
    }
#endif
#if USE_Nc3
    case 3: {
        mLatDirProp3 *r = qlua_newLatDirProp3(L, lua_gettop(L), 3);
        CALL_QDP(L);
        QDP_D3_P_eq_gaussian_S(r->ptr, a->ptr, *S->qss);
        return 1;
    }
#endif
#if USE_NcN
    default: {
        mLatDirPropN *r = qlua_newLatDirPropN(L, lua_gettop(L), nc);
        CALL_QDP(L);
        QDP_DN_P_eq_gaussian_S(r->ptr, a->ptr, *S->qss);
        return 1;
    }
#else
    default:
        return luaL_error(L, "bad number of colors");
#endif
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
 *  L:DiracPropagator(p)            -- constant propagator fill
 *  L:DiracPropagator(P)            -- lattice propagator copy
 *  L:DiracPropagator(M)            -- default colors, P[.,k,.,k] = M for all k
 *  L:DiracPropagator(D,{c=i,d=j})  -- default colors, P[i,j,.,.] = D
 */
static int
q_latdirprop(lua_State *L)
{
    mLattice *S = qlua_checkLattice(L, 1);

    if (lua_gettop(L) == 2) {
        switch (qlua_qtype(L, 2)) {
#if USE_Nc2
        case qSeqDirProp2: return q_latdirprop_seq_2(L, S, 2);
        case qLatDirProp2: return q_latdirprop_lat_2(L, S, 2);
        case qLatColMat2:  return q_latdirprop_mat_2(L, S, 2);
#endif
#if USE_Nc3
        case qSeqDirProp3: return q_latdirprop_seq_3(L, S, 3);
        case qLatDirProp3: return q_latdirprop_lat_3(L, S, 3);
        case qLatColMat3:  return q_latdirprop_mat_3(L, S, 3);
#endif
#if USE_NcN
        case qSeqDirPropN: {
            mSeqDirPropN *v = qlua_checkSeqDirPropN(L, 2, -1);
            return q_latdirprop_seq_N(L, S, v->nc);
        }
        case qLatDirPropN: {
            mLatDirPropN *v = qlua_checkLatDirPropN(L, 2, S, -1);
            return q_latdirprop_lat_N(L, S, v->nc);
        }
        case qLatColMatN: {
            mLatColMatN *v = qlua_checkLatColMatN(L, 2, S, -1);
            return q_latdirprop_mat_N(L, S, v->nc);
        }
#endif
        default:
            break;
        }
    }
    switch (S->nc) {
#if USE_Nc2
    case 2: return q_latdirprop_2(L, S, 2, 0);
#endif
#if USE_Nc3
    case 3: return q_latdirprop_3(L, S, 3, 0);
#endif
#if USE_NcN
    default: return q_latdirprop_N(L, S, S->nc, 0);
#else
    default: return luaL_error(L, "bad number of colors");
#endif
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
#if USE_Nc2 || USE_Nc3 || USE_NcN
    mLattice *S = qlua_checkLattice(L, 1);
#endif
    int nc = luaL_checkint(L, 2);

    switch (nc) {
#if USE_Nc2
    case 2: return q_latdirprop_2(L, S, 2, 1);
#endif
#if USE_Nc3
    case 3: return q_latdirprop_3(L, S, 3, 1);
#endif
#if USE_NcN
    default: return q_latdirprop_N(L, S, nc, 1);
#else
    default: return luaL_error(L, "bad number of colors");
#endif
    }
}

#if USE_Nc3
typedef struct {
    QLA_D3_DiracPropagator *a;
    QLA_D3_DiracPropagator *b;
} QQ_arg;

/* TODO force to unroll loops over color&spin indices? */
#define do_QQ_contract_func(contract_idx,A,B,C,D)\
    static void do_QQc ## contract_idx(QLA_D3_DiracPropagator *q3, \
                                       int idx,                         \
                                       void *env)                       \
 {                                                                      \
     QQ_arg *arg = env;                                                 \
     QLA_D3_DiracPropagator *q1 = &arg->a[idx];                         \
     QLA_D3_DiracPropagator *q2 = &arg->b[idx];                         \
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
                      void (*do_contract)(QLA_D3_DiracPropagator *,int,void *))
{
    mLatDirProp3 *a = qlua_checkLatDirProp3(L, 1, NULL, 3); /* the 1st arg */
    mLattice *S = qlua_ObjLattice(L, 1);
    mLatDirProp3 *b = qlua_checkLatDirProp3(L, 2, S, 3); /* the 2nd arg */
    mLatDirProp3 *r = qlua_newLatDirProp3(L, lua_gettop(L), 3); 
    QQ_arg arg;

    CALL_QDP(L);
    arg.a = QDP_D3_expose_P(a->ptr);
    arg.b = QDP_D3_expose_P(b->ptr);
    QDP_D3_P_eq_funcia(r->ptr, do_contract, &arg, *S->qss);
    QDP_D3_reset_P(a->ptr);
    QDP_D3_reset_P(b->ptr);

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
#endif /* USE_Nc3 */

static struct luaL_Reg fLatDirProp[] = {
    { "DiracPropagator",   q_latdirprop  },
    { "DiracPropagatorN",  q_latdirpropN },
    { NULL,                NULL          }
};

static struct luaL_Reg fQCDDirProp[] = {
#if USE_Nc3
    { "quarkContract12",   q_su3contract12 },
    { "quarkContract13",   q_su3contract13 },
    { "quarkContract14",   q_su3contract14 },
    { "quarkContract23",   q_su3contract23 },
    { "quarkContract24",   q_su3contract24 },
    { "quarkContract34",   q_su3contract34 },
#endif
    { NULL,                NULL            }
};

int
init_latdirprop(lua_State *L)
{
    luaL_register(L, qcdlib, fQCDDirProp);
    luaL_getmetatable(L, opLattice);
    luaL_register(L, NULL, fLatDirProp);
    lua_pop(L, 1);
#if USE_Nc2
    qlua_reg_op2(ops2);
    qlua_reg_dot(qLatDirProp2,  q_P_dot_2);
#endif
#if USE_Nc3
    qlua_reg_op2(ops3);
    qlua_reg_dot(qLatDirProp3,  q_P_dot_3);
#endif
#if USE_NcN
    qlua_reg_op2(opsN);
    qlua_reg_dot(qLatDirPropN,  q_P_dot_N);
#endif

    return 0;
}

int
fini_latdirprop(lua_State *L)
{
    return 0;
}
