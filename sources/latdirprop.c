#include <qlua.h>                                                    /* DEPS */
#include <qcomplex.h>                                                /* DEPS */
#include <lattice.h>                                                 /* DEPS */
#include <latrandom.h>                                               /* DEPS */
#include <latreal.h>                                                 /* DEPS */
#include <latcomplex.h>                                              /* DEPS */
#include <latcolmat.h>                                               /* DEPS */
#include <latdirferm.h>                                              /* DEPS */
#include <latdirprop.h>                                              /* DEPS */
#include <latint.h>                                                  /* DEPS */
#include <qmp.h>

const char mtnLatDirProp[] = "lattice.DiracPropagator";
static const char opLatDirProp[] = "lattice.DiracPropagator.ops";

mLatDirProp *
qlua_newLatDirProp(lua_State *L)
{
    QDP_DiracPropagator *v = QDP_create_P();
    mLatDirProp *hdr;

    if (v == 0) {
        lua_gc(L, LUA_GCCOLLECT, 0);
        v = QDP_create_P();
        if (v == 0)
            luaL_error(L, "not enough memory (QDP_DiracPropagator)");
    }
    hdr = lua_newuserdata(L, sizeof (mLatDirProp));
    hdr->ptr = v;
    luaL_getmetatable(L, mtnLatDirProp);
    lua_setmetatable(L, -2);

    return hdr;
}

mLatDirProp *
qlua_checkLatDirProp(lua_State *L, int idx)
{
    void *v = luaL_checkudata(L, idx, mtnLatDirProp);
    
    luaL_argcheck(L, v != 0, idx, "lattice.DiracPropagator expected");
    
    return v;
}

static int
q_P_fmt(lua_State *L)
{
    char fmt[72];
    mLatDirProp *b = qlua_checkLatDirProp(L, 1);

    sprintf(fmt, "QDP:DiracPropagator(%p)", b->ptr);
    lua_pushstring(L, fmt);

    return 1;
}

static int
q_P_gc(lua_State *L)
{
    mLatDirProp *b = qlua_checkLatDirProp(L, 1);

    QDP_destroy_P(b->ptr);
    b->ptr = 0;

    return 0;
}

static int
q_P_get(lua_State *L)
{
    switch (qlua_gettype(L, 2)) {
    case qTable: {
        mLatDirProp *V = qlua_checkLatDirProp(L, 1);
        int d = qlua_checkdiracindex(L, 2);
        int c = qlua_checkcolorindex(L, 2);
        mLatDirFerm *r = qlua_newLatDirFerm(L);
                
        CALL_QDP(L);
        QDP_D_eq_diracvec_P(r->ptr, V->ptr, c, d, *qCurrent);

        return 1;
    }
    case qString:
        return qlua_lookup(L, 2, opLatDirProp);
    }
    return qlua_badindex(L, "DiracPropagator");
}

static int
q_P_put(lua_State *L)
{
    mLatDirProp *V = qlua_checkLatDirProp(L, 1);
    int d = qlua_checkdiracindex(L, 2);
    int c = qlua_checkcolorindex(L, 2);
    mLatDirFerm *r = qlua_checkLatDirFerm(L, 3);
            
    CALL_QDP(L);
    QDP_P_eq_diracvec_D(V->ptr, r->ptr, c, d, *qCurrent);

    return 0;
}

int
q_P_gaussian(lua_State *L)
{
    mLatRandom *a = qlua_checkLatRandom(L, 1);
    mLatDirProp *r = qlua_newLatDirProp(L);

    CALL_QDP(L);
    QDP_P_eq_gaussian_S(r->ptr, a->ptr, *qCurrent);

    return 1;
}

static int
q_P_norm2(lua_State *L)
{
    mLatDirProp *a = qlua_checkLatDirProp(L, 1);
    QLA_Real n;

    CALL_QDP(L);
    QDP_r_eq_norm2_P(&n, a->ptr, *qCurrent);
    lua_pushnumber(L, n);
    
    return 1;
}

static int
q_P_shift(lua_State *L)
{
    mLatDirProp *a = qlua_checkLatDirProp(L, 1);
    QDP_Shift shift = qlua_checkShift(L, 2);
    QDP_ShiftDir dir = qlua_checkShiftDir(L, 3);
    mLatDirProp *r = qlua_newLatDirProp(L);

    CALL_QDP(L);
    QDP_P_eq_sP(r->ptr, a->ptr, shift, dir, *qCurrent);

    return 1;
}

static int
q_P_conj(lua_State *L)
{
    mLatDirProp *a = qlua_checkLatDirProp(L, 1);
    mLatDirProp *r = qlua_newLatDirProp(L);

    CALL_QDP(L);
    QDP_P_eq_conj_P(r->ptr, a->ptr, *qCurrent);

    return 1;
}

static int
q_P_trans(lua_State *L)
{
    mLatDirProp *a = qlua_checkLatDirProp(L, 1);
    mLatDirProp *r = qlua_newLatDirProp(L);

    CALL_QDP(L);
    QDP_P_eq_transpose_P(r->ptr, a->ptr, *qCurrent);

    return 1;
}

static int
q_P_adjoin(lua_State *L)
{
    mLatDirProp *a = qlua_checkLatDirProp(L, 1);
    mLatDirProp *r = qlua_newLatDirProp(L);

    CALL_QDP(L);
    QDP_P_eq_Pa(r->ptr, a->ptr, *qCurrent);

    return 1;
}

static int
q_P_spintrace(lua_State *L)
{
    mLatDirProp *a = qlua_checkLatDirProp(L, 1);
    mLatColMat *r = qlua_newLatColMat(L);

    CALL_QDP(L);
    QDP_M_eq_spintrace_P(r->ptr, a->ptr, *qCurrent);

    return 1;
}

static struct {
    QLA_DiracPropagator *a;
} Pst_args; /* YYY global state */

static void
do_Pst(QLA_DiracPropagator *r, int idx)
{
    QLA_DiracPropagator *a = &Pst_args.a[idx];
    int is, ic, js, jc;

    for (is = 0; is < QDP_Ns; is++) {
        for (js = 0; js < QDP_Ns; js++) {
            for (ic = 0; ic < QDP_Nc; ic++)  {
                for (jc = 0; jc < QDP_Nc; jc++) {
                    QLA_c_eq_c(QLA_D3_elem_P(*r,ic,is,jc,js),
                               QLA_D3_elem_P(*a,ic,js,jc,is));
                }
            }
        }
    }
}

static int
q_P_spintranspose(lua_State *L)
{
    mLatDirProp *a = qlua_checkLatDirProp(L, 1);
    mLatDirProp *r = qlua_newLatDirProp(L);

    CALL_QDP(L);
    Pst_args.a = QDP_expose_P(a->ptr);
    QDP_P_eq_funci(r->ptr, do_Pst, *qCurrent);
    QDP_reset_P(a->ptr);
    Pst_args.a = 0;

    return 1;
}

static int
q_P_set(lua_State *L)
{
    mLatDirProp *r = qlua_checkLatDirProp(L, 1);
    mLatDirProp *a = qlua_checkLatDirProp(L, 2);

    CALL_QDP(L);
    QDP_P_eq_P(r->ptr, a->ptr, *qCurrent);
    lua_pop(L, 1);

    return 1;
}

static int
q_P_add_P(lua_State *L)
{
    mLatDirProp *a = qlua_checkLatDirProp(L, 1);
    mLatDirProp *b = qlua_checkLatDirProp(L, 2);
    mLatDirProp *c = qlua_newLatDirProp(L);

    CALL_QDP(L);
    QDP_P_eq_P_plus_P(c->ptr, a->ptr, b->ptr, *qCurrent);

    return 1;
}

static int
q_P_sub_P(lua_State *L)
{
    mLatDirProp *a = qlua_checkLatDirProp(L, 1);
    mLatDirProp *b = qlua_checkLatDirProp(L, 2);
    mLatDirProp *c = qlua_newLatDirProp(L);

    CALL_QDP(L);
    QDP_P_eq_P_minus_P(c->ptr, a->ptr, b->ptr, *qCurrent);

    return 1;
}


static int
q_P_mul_r(lua_State *L)
{
    mLatDirProp *a = qlua_checkLatDirProp(L, 1);
    QLA_Real b = luaL_checknumber(L, 2);
    mLatDirProp *c = qlua_newLatDirProp(L);

    CALL_QDP(L);
    QDP_P_eq_r_times_P(c->ptr, &b, a->ptr, *qCurrent);

    return 1;
}

static int
q_r_mul_P(lua_State *L)
{
    QLA_Real a = luaL_checknumber(L, 1);
    mLatDirProp *b = qlua_checkLatDirProp(L, 2);
    mLatDirProp *c = qlua_newLatDirProp(L);

    CALL_QDP(L);
    QDP_P_eq_r_times_P(c->ptr, &a, b->ptr, *qCurrent);

    return 1;
}

static int
q_P_mul_c(lua_State *L)
{
    mLatDirProp *a = qlua_checkLatDirProp(L, 1);
    QLA_Complex *b = qlua_checkComplex(L, 2);
    mLatDirProp *c = qlua_newLatDirProp(L);

    CALL_QDP(L);
    QDP_P_eq_c_times_P(c->ptr, b, a->ptr, *qCurrent);

    return 1;
}

static int
q_c_mul_P(lua_State *L)
{
    QLA_Complex *a = qlua_checkComplex(L, 1);
    mLatDirProp *b = qlua_checkLatDirProp(L, 2);
    mLatDirProp *c = qlua_newLatDirProp(L);

    CALL_QDP(L);
    QDP_P_eq_c_times_P(c->ptr, a, b->ptr, *qCurrent);

    return 1;
}

static struct {
    QLA_Real *a;
    QLA_DiracPropagator *b;
} RPmul_args; /* YYY global state */

static void
do_RPmul(QLA_DiracPropagator *r, int idx)
{
    QLA_P_eq_r_times_P(r, &RPmul_args.a[idx], &RPmul_args.b[idx]);
}

static void
X_P_eq_R_times_P(QDP_DiracPropagator *r, QDP_Real *a, QDP_DiracPropagator *b,
                 QDP_Subset s)
{
    RPmul_args.a = QDP_expose_R(a);
    RPmul_args.b = QDP_expose_P(b);
    QDP_P_eq_funci(r, do_RPmul, s);
    QDP_reset_R(a);
    QDP_reset_P(b);
    RPmul_args.a = 0;
    RPmul_args.b = 0;
}

static int
q_R_mul_P(lua_State *L)
{
    mLatReal *a = qlua_checkLatReal(L, 1);
    mLatDirProp *b = qlua_checkLatDirProp(L, 2);
    mLatDirProp *r = qlua_newLatDirProp(L);

    CALL_QDP(L);
    X_P_eq_R_times_P(r->ptr, a->ptr, b->ptr, *qCurrent);

    return 1;
}

static int
q_P_mul_R(lua_State *L)
{
    mLatDirProp *a = qlua_checkLatDirProp(L, 1);
    mLatReal *b = qlua_checkLatReal(L, 2);
    mLatDirProp *r = qlua_newLatDirProp(L);

    CALL_QDP(L);
    X_P_eq_R_times_P(r->ptr, b->ptr, a->ptr, *qCurrent);

    return 1;
}

static struct {
    QLA_Complex *a;
    QLA_DiracPropagator *b;
} CPmul_args; /* YYY global state */

static void
do_CPmul(QLA_DiracPropagator *r, int idx)
{
    QLA_P_eq_c_times_P(r, &CPmul_args.a[idx], &CPmul_args.b[idx]);
}

static void
X_P_eq_C_times_P(QDP_DiracPropagator *r, QDP_Complex *a, QDP_DiracPropagator *b,
                 QDP_Subset s)
{
    CPmul_args.a = QDP_expose_C(a);
    CPmul_args.b = QDP_expose_P(b);
    QDP_P_eq_funci(r, do_CPmul, s);
    QDP_reset_C(a);
    QDP_reset_P(b);
    CPmul_args.a = 0;
    CPmul_args.b = 0;
}

static int
q_C_mul_P(lua_State *L)
{
    mLatComplex *a = qlua_checkLatComplex(L, 1);
    mLatDirProp *b = qlua_checkLatDirProp(L, 2);
    mLatDirProp *r = qlua_newLatDirProp(L);

    CALL_QDP(L);
    X_P_eq_C_times_P(r->ptr, a->ptr, b->ptr, *qCurrent);

    return 1;
}

static int
q_P_mul_C(lua_State *L)
{
    mLatDirProp *a = qlua_checkLatDirProp(L, 1);
    mLatComplex *b = qlua_checkLatComplex(L, 2);
    mLatDirProp *r = qlua_newLatDirProp(L);

    CALL_QDP(L);
    X_P_eq_C_times_P(r->ptr, b->ptr, a->ptr, *qCurrent);

    return 1;
}

static int
q_P_neg(lua_State *L)
{
    mLatDirProp *a = qlua_checkLatDirProp(L, 1);
    mLatDirProp *r = qlua_newLatDirProp(L);
    QLA_Real m1 = -1;

    CALL_QDP(L);
    QDP_P_eq_r_times_P(r->ptr, &m1, a->ptr, *qCurrent);

    return 1;
}

static int
q_P_mul_P(lua_State *L)
{
    mLatDirProp *a = qlua_checkLatDirProp(L, 1);
    mLatDirProp *b = qlua_checkLatDirProp(L, 2);
    mLatDirProp *c = qlua_newLatDirProp(L);

    CALL_QDP(L);
    QDP_P_eq_P_times_P(c->ptr, a->ptr, b->ptr, *qCurrent);

    return 1;
}

static int
q_P_mul_M(lua_State *L)
{
    mLatDirProp *a = qlua_checkLatDirProp(L, 1);
    mLatColMat *b = qlua_checkLatColMat(L, 2);
    mLatDirProp *c = qlua_newLatDirProp(L);

    CALL_QDP(L);
    QDP_P_eq_P_times_M(c->ptr, a->ptr, b->ptr, *qCurrent);

    return 1;
}

static int
q_M_mul_P(lua_State *L)
{
    mLatColMat *a = qlua_checkLatColMat(L, 1);
    mLatDirProp *b = qlua_checkLatDirProp(L, 2);
    mLatDirProp *c = qlua_newLatDirProp(L);

    CALL_QDP(L);
    QDP_P_eq_M_times_P(c->ptr, a->ptr, b->ptr, *qCurrent);

    return 1;
}

static int
q_P_div_r(lua_State *L)
{
    mLatDirProp *a = qlua_checkLatDirProp(L, 1);
    QLA_Real b = 1 / luaL_checknumber(L, 2);
    mLatDirProp *c = qlua_newLatDirProp(L);

    CALL_QDP(L);
    QDP_P_eq_r_times_P(c->ptr, &b, a->ptr, *qCurrent);

    return 1;
}

static int
q_P_div_c(lua_State *L)
{
    mLatDirProp *a = qlua_checkLatDirProp(L, 1);
    QLA_Complex *b = qlua_checkComplex(L, 2);
    mLatDirProp *c = qlua_newLatDirProp(L);
    double n = 1 / (QLA_real(*b) * QLA_real(*b) + QLA_imag(*b) * QLA_imag(*b));
    QLA_Complex s;

    CALL_QDP(L);
    QLA_real(s) = n * QLA_real(*b);
    QLA_imag(s) = -n * QLA_imag(*b);
    QDP_P_eq_c_times_P(c->ptr, &s, a->ptr, *qCurrent);

    return 1;
}

static int
q_P_dot(lua_State *L)
{
    mLatDirProp *a = qlua_checkLatDirProp(L, 1);
    mLatDirProp *b = qlua_checkLatDirProp(L, 2);
    mLatComplex *s = qlua_newLatComplex(L);

    CALL_QDP(L);
    QDP_C_eq_P_dot_P(s->ptr, a->ptr, b->ptr, *qCurrent);

    return 1;
}

static int
q_latdirprop(lua_State *L)
{
    switch (lua_gettop(L)) {
    case 1: {
        mLatDirProp *v = qlua_newLatDirProp(L);

        CALL_QDP(L);
        QDP_P_eq_zero(v->ptr, *qCurrent);
        
        return 1;
    }
    case 2: {
        mLatDirProp *v = qlua_newLatDirProp(L);
        
        switch (qlua_gettype(L, 2)) {
        case qLatDirProp: {
            mLatDirProp *w = qlua_checkLatDirProp(L, 2);

            CALL_QDP(L);
            QDP_P_eq_P(v->ptr, w->ptr, *qCurrent);

            return 1;
        }
        case qLatColMat: {
            mLatColMat *w = qlua_checkLatColMat(L, 2);
            mLatComplex *c = qlua_newLatComplex(L);
            int ic, jc, ks;

            CALL_QDP(L);
            QDP_P_eq_zero(v->ptr, *qCurrent);
            for (ic = 0; ic < QDP_Nc; ic++) {
                for (jc = 0; jc < QDP_Nc; jc++) {
                    QDP_C_eq_elem_M(c->ptr, w->ptr, ic, jc, *qCurrent);
                    for (ks = 0; ks < QDP_Ns; ks++)
                        QDP_P_eq_elem_C(v->ptr, c->ptr, ic, ks, jc, ks, *qCurrent);
                }
            }
            lua_pop(L, 1);

            return 1;
        }
        default:
            break;
        }
    }
    case 3: {
        mLatDirFerm *z = qlua_checkLatDirFerm(L, 2);
        int d = qlua_checkdiracindex(L, 3);
        int c = qlua_checkcolorindex(L, 3);
        mLatDirProp *v = qlua_newLatDirProp(L);

        CALL_QDP(L);
        QDP_P_eq_zero(v->ptr, *qCurrent);
        QDP_P_eq_diracvec_D(v->ptr, z->ptr, c, d, *qCurrent);

        return 1;
    }
    }
    return qlua_badconstr(L, "DiracPropagator");
}

static struct {
    QLA_DiracPropagator *a;
    QLA_DiracPropagator *b;
} QQ_args;

/* TODO force to unroll loops over color&spin indices? */
#define do_QQ_contract_func(contract_idx,A,B,C,D)\
static void do_QQc ## contract_idx(QLA_DiracPropagator *q3, int idx) \
{ \
    QLA_DiracPropagator *q1 = &QQ_args.a[idx]; \
    QLA_DiracPropagator *q2 = &QQ_args.b[idx]; \
    static const int eps[3][3] = { { 0, 1, 2}, { 1, 2, 0}, { 2, 0, 1} }; \
    int p_a, p_b; \
    for (p_a = 0; p_a < QDP_Nc; p_a++) { \
        int i1 = eps[p_a][0], j1 = eps[p_a][1], k1 = eps[p_a][2]; \
        for (p_b = 0; p_b < QDP_Nc; p_b++) { \
            int i2 = eps[p_b][0], j2 = eps[p_b][1], k2 = eps[p_b][2]; \
            int a, b, c; \
            for (a = 0; a < QDP_Ns; a++) { \
                for (b = 0; b < QDP_Ns; b++) { \
                    QLA_Complex s3; \
                    QLA_c_eq_r(s3, 0.0); \
                    for (c = 0; c < QDP_Ns; c++) { \
                        QLA_c_peq_c_times_c(s3,QLA_elem_P(*q1,i1,(A),i2,(B)), \
                                               QLA_elem_P(*q2,j1,(C),j2,(D)));\
                        QLA_c_meq_c_times_c(s3,QLA_elem_P(*q1,i1,(A),j2,(B)), \
                                               QLA_elem_P(*q2,j1,(C),i2,(D)));\
                        QLA_c_meq_c_times_c(s3,QLA_elem_P(*q1,j1,(A),i2,(B)), \
                                               QLA_elem_P(*q2,i1,(C),j2,(D)));\
                        QLA_c_peq_c_times_c(s3,QLA_elem_P(*q1,j1,(A),j2,(B)), \
                                               QLA_elem_P(*q2,i1,(C),i2,(D)));\
                    } \
                    QLA_c_eq_c(QLA_elem_P(*q3,k2,a,k1,b), s3); \
                } \
            } \
        } \
    } \
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
        void (*do_contract)(QLA_DiracPropagator *, int))
{
    mLatDirProp *a = qlua_checkLatDirProp(L, 1); /* the first argument */
    mLatDirProp *b = qlua_checkLatDirProp(L, 2); /* the second argument */
    mLatDirProp *r = qlua_newLatDirProp(L); 

    if (QDP_Nc != 3)
        return luaL_error(L, "Bad value of qcd.Nc");

    CALL_QDP(L);
    QQ_args.a = QDP_expose_P(a->ptr);
    QQ_args.b = QDP_expose_P(b->ptr);
    QDP_P_eq_funci(r->ptr, do_contract, *qCurrent);
    QDP_reset_P(a->ptr);
    QDP_reset_P(b->ptr);
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

static struct luaL_Reg LatDirPropMethods[] = {
    { "norm2",           q_P_norm2 },
    { "shift",           q_P_shift },
    { "conj",            q_P_conj },
    { "transpose",       q_P_trans },
    { "adjoin",          q_P_adjoin },
    { "spintrace",       q_P_spintrace },
    { "spintranspose",   q_P_spintranspose },
    { "set",             q_P_set },
    { NULL,              NULL }
};

static struct luaL_Reg mtLatDirProp[] = {
    { "__tostring",        q_P_fmt },
    { "__gc",              q_P_gc },
    { "__index",           q_P_get },
    { "__newindex",        q_P_put },
    { "__unm",             q_P_neg },
    { "__add",             qlua_add },
    { "__sub",             qlua_sub },
    { "__mul",             qlua_mul },
    { "__div",             qlua_div },
    { NULL,                NULL }
};

static struct luaL_Reg fLatDirProp[] = {
    { "DiracPropagator",   q_latdirprop },
    { NULL,                NULL }
};

static struct luaL_Reg fQCDDirProp[] = {
    { "quarkContract12",   q_su3contract12 },
    { "quarkContract13",   q_su3contract13 },
    { "quarkContract14",   q_su3contract14 },
    { "quarkContract23",   q_su3contract23 },
    { "quarkContract24",   q_su3contract24 },
    { "quarkContract34",   q_su3contract34 },
    { NULL,                NULL }
};

int
init_latdirprop(lua_State *L)
{
    luaL_register(L, qcdlib, fQCDDirProp);
    luaL_getmetatable(L, opLattice);
    luaL_register(L, NULL, fLatDirProp);
    lua_pop(L, 1);
    qlua_metatable(L, mtnLatDirProp, mtLatDirProp);
    qlua_metatable(L, opLatDirProp, LatDirPropMethods);
    qlua_reg_add(qLatDirProp,  qLatDirProp,  q_P_add_P);
    qlua_reg_sub(qLatDirProp,  qLatDirProp,  q_P_sub_P);
    qlua_reg_mul(qReal,        qLatDirProp,  q_r_mul_P);
    qlua_reg_mul(qLatDirProp,  qReal,        q_P_mul_r);
    qlua_reg_mul(qComplex,     qLatDirProp,  q_c_mul_P);
    qlua_reg_mul(qLatDirProp,  qComplex,     q_P_mul_c);
    qlua_reg_mul(qLatReal,     qLatDirProp,  q_R_mul_P);
    qlua_reg_mul(qLatDirProp,  qLatReal,     q_P_mul_R);
    qlua_reg_mul(qLatComplex,  qLatDirProp,  q_C_mul_P);
    qlua_reg_mul(qLatDirProp,  qLatComplex,  q_P_mul_C);
    qlua_reg_mul(qLatDirProp,  qLatDirProp,  q_P_mul_P);
    qlua_reg_mul(qLatDirProp,  qLatColMat,   q_P_mul_M);
    qlua_reg_mul(qLatColMat,   qLatDirProp,  q_M_mul_P);
    qlua_reg_div(qLatDirProp,  qReal,        q_P_div_r);
    qlua_reg_div(qLatDirProp,  qComplex,     q_P_div_c);
    qlua_reg_dot(qLatDirProp,  q_P_dot);

    return 0;
}

int
fini_latdirprop(lua_State *L)
{
    return 0;
}
