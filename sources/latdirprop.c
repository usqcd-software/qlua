#include <qlua.h>                                                    /* DEPS */
#include <qcomplex.h>                                                /* DEPS */
#include <lattice.h>                                                 /* DEPS */
#include <latcolmat.h>                                               /* DEPS */
#include <latdirferm.h>                                              /* DEPS */
#include <latdirprop.h>                                              /* DEPS */
#include <latrandom.h>                                               /* DEPS */
#include <latcomplex.h>                                              /* DEPS */
#include <latint.h>                                                  /* DEPS */
#include <qmp.h>

const char mtnLatDirProp[] = "qcd.lattice.DiracPropagator";
static const char opLatDirProp[] = "qcd.lattice.DiracPropagator.op";

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
    
    luaL_argcheck(L, v != 0, idx, "qcd.DiracPropagator expected");
    
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
                
        QDP_D_eq_diracvec_P(r->ptr, V->ptr, c, d, qCurrent);

        return 1;
    }
    case qString:
        return qlua_lookup(L, 2, opLatDirProp);
    }
    return luaL_error(L, "bad index");
}

static int
q_P_put(lua_State *L)
{
    mLatDirProp *V = qlua_checkLatDirProp(L, 1);
    int d = qlua_checkdiracindex(L, 2);
    int c = qlua_checkcolorindex(L, 2);
    mLatDirFerm *r = qlua_checkLatDirFerm(L, 3);
            
    QDP_D_eq_diracvec_P(r->ptr, V->ptr, c, d, qCurrent);

    return 0;
}

int
q_P_gaussian(lua_State *L)
{
    mLatRandom *a = qlua_checkLatRandom(L, 1);
    mLatDirProp *r = qlua_newLatDirProp(L);

    QDP_P_eq_gaussian_S(r->ptr, a->ptr, qCurrent);

    return 1;
}

static int
q_P_norm2(lua_State *L)
{
    mLatDirProp *a = qlua_checkLatDirProp(L, 1);
    QLA_Real n;

    QDP_r_eq_norm2_P(&n, a->ptr, qCurrent);
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

    QDP_P_eq_sP(r->ptr, a->ptr, shift, dir, qCurrent);

    return 1;
}

static int
q_P_conj(lua_State *L)
{
    mLatDirProp *a = qlua_checkLatDirProp(L, 1);
    mLatDirProp *r = qlua_newLatDirProp(L);

    QDP_P_eq_conj_P(r->ptr, a->ptr, qCurrent);

    return 1;
}

static int
q_P_trans(lua_State *L)
{
    mLatDirProp *a = qlua_checkLatDirProp(L, 1);
    mLatDirProp *r = qlua_newLatDirProp(L);

    QDP_P_eq_transpose_P(r->ptr, a->ptr, qCurrent);

    return 1;
}

static int
q_P_adjoin(lua_State *L)
{
    mLatDirProp *a = qlua_checkLatDirProp(L, 1);
    mLatDirProp *r = qlua_newLatDirProp(L);

    QDP_P_eq_Pa(r->ptr, a->ptr, qCurrent);

    return 1;
}

static int
q_P_spintrace(lua_State *L)
{
    mLatDirProp *a = qlua_checkLatDirProp(L, 1);
    mLatColMat *r = qlua_newLatColMat(L);

    QDP_M_eq_spintrace_P(r->ptr, a->ptr, qCurrent);

    return 1;
}

static int
q_P_set(lua_State *L)
{
    mLatDirProp *r = qlua_checkLatDirProp(L, 1);
    mLatDirProp *a = qlua_checkLatDirProp(L, 2);

    QDP_P_eq_P(r->ptr, a->ptr, qCurrent);
    lua_pop(L, 1);

    return 1;
}

static int
q_P_add_P(lua_State *L)
{
    mLatDirProp *a = qlua_checkLatDirProp(L, 1);
    mLatDirProp *b = qlua_checkLatDirProp(L, 2);
    mLatDirProp *c = qlua_newLatDirProp(L);

    QDP_P_eq_P_plus_P(c->ptr, a->ptr, b->ptr, qCurrent);

    return 1;
}

static int
q_P_sub_P(lua_State *L)
{
    mLatDirProp *a = qlua_checkLatDirProp(L, 1);
    mLatDirProp *b = qlua_checkLatDirProp(L, 2);
    mLatDirProp *c = qlua_newLatDirProp(L);

    QDP_P_eq_P_minus_P(c->ptr, a->ptr, b->ptr, qCurrent);

    return 1;
}


static int
q_P_mul_r(lua_State *L)
{
    mLatDirProp *a = qlua_checkLatDirProp(L, 1);
    QLA_Real b = luaL_checknumber(L, 2);
    mLatDirProp *c = qlua_newLatDirProp(L);

    QDP_P_eq_r_times_P(c->ptr, &b, a->ptr, qCurrent);

    return 1;
}

static int
q_r_mul_P(lua_State *L)
{
    QLA_Real a = luaL_checknumber(L, 1);
    mLatDirProp *b = qlua_checkLatDirProp(L, 2);
    mLatDirProp *c = qlua_newLatDirProp(L);

    QDP_P_eq_r_times_P(c->ptr, &a, b->ptr, qCurrent);

    return 1;
}

static int
q_P_mul_c(lua_State *L)
{
    mLatDirProp *a = qlua_checkLatDirProp(L, 1);
    QLA_Complex *b = qlua_checkComplex(L, 2);
    mLatDirProp *c = qlua_newLatDirProp(L);

    QDP_P_eq_c_times_P(c->ptr, b, a->ptr, qCurrent);

    return 1;
}

static int
q_c_mul_P(lua_State *L)
{
    QLA_Complex *a = qlua_checkComplex(L, 1);
    mLatDirProp *b = qlua_checkLatDirProp(L, 2);
    mLatDirProp *c = qlua_newLatDirProp(L);

    QDP_P_eq_c_times_P(c->ptr, a, b->ptr, qCurrent);

    return 1;
}

static int
q_P_neg(lua_State *L)
{
    mLatDirProp *a = qlua_checkLatDirProp(L, 1);
    mLatDirProp *r = qlua_newLatDirProp(L);
    QLA_Real m1 = -1;

    QDP_P_eq_r_times_P(r->ptr, &m1, a->ptr, qCurrent);

    return 1;
}

static int
q_P_mul_P(lua_State *L)
{
    mLatDirProp *a = qlua_checkLatDirProp(L, 1);
    mLatDirProp *b = qlua_checkLatDirProp(L, 2);
    mLatDirProp *c = qlua_newLatDirProp(L);

    QDP_P_eq_P_times_P(c->ptr, a->ptr, b->ptr, qCurrent);

    return 1;
}

static int
q_P_mul_M(lua_State *L)
{
    mLatDirProp *a = qlua_checkLatDirProp(L, 1);
    mLatColMat *b = qlua_checkLatColMat(L, 2);
    mLatDirProp *c = qlua_newLatDirProp(L);

    QDP_P_eq_P_times_M(c->ptr, a->ptr, b->ptr, qCurrent);

    return 1;
}

static int
q_M_mul_P(lua_State *L)
{
    mLatColMat *a = qlua_checkLatColMat(L, 1);
    mLatDirProp *b = qlua_checkLatDirProp(L, 2);
    mLatDirProp *c = qlua_newLatDirProp(L);

    QDP_P_eq_M_times_P(c->ptr, a->ptr, b->ptr, qCurrent);

    return 1;
}

static int
q_P_div_r(lua_State *L)
{
    mLatDirProp *a = qlua_checkLatDirProp(L, 1);
    QLA_Real b = 1 / luaL_checknumber(L, 2);
    mLatDirProp *c = qlua_newLatDirProp(L);

    QDP_P_eq_r_times_P(c->ptr, &b, a->ptr, qCurrent);

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

    QLA_real(s) = n * QLA_real(*b);
    QLA_imag(s) = -n * QLA_imag(*b);
    QDP_P_eq_c_times_P(c->ptr, &s, a->ptr, qCurrent);

    return 1;
}

static int
q_P_dot(lua_State *L)
{
    mLatDirProp *a = qlua_checkLatDirProp(L, 1);
    mLatDirProp *b = qlua_checkLatDirProp(L, 2);
    mLatComplex *s = qlua_newLatComplex(L);

    QDP_C_eq_P_dot_P(s->ptr, a->ptr, b->ptr, qCurrent);

    return 1;
}

static int
q_latdirprop(lua_State *L)
{
    switch (lua_gettop(L)) {
    case 1: {
        mLatDirProp *v = qlua_newLatDirProp(L);

        QDP_P_eq_zero(v->ptr, qCurrent);
        
        return 1;
    }
    case 2: {
        mLatDirProp *w = qlua_checkLatDirProp(L, 2);
        mLatDirProp *v = qlua_newLatDirProp(L);
        
        QDP_P_eq_P(v->ptr, w->ptr, qCurrent);
        
        return 1;
    }
    case 3: {
        mLatDirFerm *z = qlua_checkLatDirFerm(L, 2);
        int d = qlua_checkdiracindex(L, 3);
        int c = qlua_checkcolorindex(L, 3);
        mLatDirProp *v = qlua_newLatDirProp(L);

        QDP_P_eq_zero(v->ptr, qCurrent);
        QDP_P_eq_diracvec_D(v->ptr, z->ptr, c, d, qCurrent);

        return 1;
    }
    }
    return qlua_badconstr(L, "DiracPropagator");
}

static struct luaL_Reg LatDirPropMethods[] = {
    { "norm2",      q_P_norm2 },
    { "shift",      q_P_shift },
    { "conj",       q_P_conj },
    { "transpose",  q_P_trans },
    { "adjoin",     q_P_adjoin },
    { "spintrace",  q_P_spintrace },
    { "set",        q_P_set },
    { NULL,         NULL }
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

int
init_latdirprop(lua_State *L)
{
    luaL_getmetatable(L, opLattice);
    luaL_register(L, NULL, fLatDirProp);
    lua_pop(L, 1);
    qlua_metatable(L, mtnLatDirProp, mtLatDirProp);
    qlua_metatable(L, opLatDirProp, LatDirPropMethods);
    qlua_reg_add(qLatDirProp, qLatDirProp, q_P_add_P);
    qlua_reg_sub(qLatDirProp, qLatDirProp, q_P_sub_P);
    qlua_reg_mul(qReal,       qLatDirProp, q_r_mul_P);
    qlua_reg_mul(qLatDirProp, qReal,       q_P_mul_r);
    qlua_reg_mul(qComplex,    qLatDirProp, q_c_mul_P);
    qlua_reg_mul(qLatDirProp, qComplex,    q_P_mul_c);
    qlua_reg_mul(qLatDirProp, qLatDirProp, q_P_mul_P);
    qlua_reg_mul(qLatDirProp, qLatColMat,  q_P_mul_M);
    qlua_reg_mul(qLatColMat,  qLatDirProp, q_M_mul_P);
    qlua_reg_div(qLatDirProp, qReal,       q_P_div_r);
    qlua_reg_div(qLatDirProp, qComplex,    q_P_div_c);
    qlua_reg_dot(qLatDirProp, q_P_dot);

    return 0;
}

int
fini_latdirprop(lua_State *L)
{
    return 0;
}
