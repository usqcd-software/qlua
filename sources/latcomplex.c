#include <qlua.h>                                                    /* DEPS */
#include <qcomplex.h>                                                /* DEPS */
#include <qvector.h>                                                 /* DEPS */
#include <lattice.h>                                                 /* DEPS */
#include <latcomplex.h>                                              /* DEPS */
#include <latreal.h>                                                 /* DEPS */
#include <latint.h>                                                  /* DEPS */
#include <latrandom.h>                                               /* DEPS */
#include <latmulti.h>                                                /* DEPS */
#include <qmp.h>

const char mtnLatComplex[] = "lattice.Complex";
static const char opLatComplex[] = "lattice.Complex.ops";

mLatComplex *
qlua_newLatComplex(lua_State *L)
{
    QDP_Complex *v = QDP_create_C();
    mLatComplex *hdr;

    if (v == 0) {
        lua_gc(L, LUA_GCCOLLECT, 0);
        v = QDP_create_C();
        if (v == 0)
            luaL_error(L, "not enough memory (QDP_Complex)");
    }
    hdr = lua_newuserdata(L, sizeof (mLatComplex));
    hdr->ptr = v;
    luaL_getmetatable(L, mtnLatComplex);
    lua_setmetatable(L, -2);

    return hdr;
}

mLatComplex *
qlua_checkLatComplex(lua_State *L, int idx)
{
    void *v = luaL_checkudata(L, idx, mtnLatComplex);

    luaL_argcheck(L, v != 0, idx, "lattice.Complex expected");
    
    return v;
}

static int
q_C_fmt(lua_State *L)
{
    char fmt[72];
    mLatComplex *b = qlua_checkLatComplex(L, 1);

    sprintf(fmt, "LatComplex(%p)", b->ptr);
    lua_pushstring(L, fmt);

    return 1;
}

static int
q_C_gc(lua_State *L)
{
    mLatComplex *b = qlua_checkLatComplex(L, 1);

    QDP_destroy_C(b->ptr);
    b->ptr = 0;

    return 0;
}

static int
q_C_neg(lua_State *L)
{
    mLatComplex *a = qlua_checkLatComplex(L, 1);
    mLatComplex *r = qlua_newLatComplex(L);
    QLA_Real m1 = -1;

    QDP_C_eq_r_times_C(r->ptr, &m1, a->ptr, qCurrent);
    return 1;
}

static int
q_C_get(lua_State *L)
{
    switch (qlua_gettype(L, 2)) {
    case qTable: {
        mLatComplex *V = qlua_checkLatComplex(L, 1);
        QLA_Complex *W;
        QLA_Complex *locked;
        double z_re;
        double z_im;
        int *idx = 0;
        
        idx = qlua_checklatcoord(L, 2);
        locked = QDP_expose_C(V->ptr);
        if (QDP_node_number(idx) == QDP_this_node) {
            QLA_Complex *zz = &QLA_elem_C(locked[QDP_index(idx)]);

            z_re = QLA_real(*zz);
            z_im = QLA_imag(*zz);
        } else {
            z_re = 0;
            z_im = 0;
        }
        QDP_reset_C(V->ptr);
        qlua_free(L, idx);
        QMP_sum_double(&z_re);
        QMP_sum_double(&z_im);
        W = qlua_newComplex(L);
        QLA_real(*W) = z_re;
        QLA_imag(*W) = z_im;

        return 1;
    }
    case qString:
        return qlua_lookup(L, 2, opLatComplex);
    }
    return qlua_badindex(L, "Complex");
}

static int
q_C_put(lua_State *L)
{
    mLatComplex *V = qlua_checkLatComplex(L, 1);
    QLA_Complex *locked;
    int *idx = 0;
    double z_re, z_im;

    switch (qlua_gettype(L, 3)) {
    case qReal:
        z_re = luaL_checknumber(L, 3);
        break;
    case qComplex: {
        QLA_Complex *z = qlua_checkComplex(L, 3);
        z_re = QLA_real(*z);
        z_im = QLA_imag(*z);
        break;
    }
    default:
        return luaL_error(L, "bad argument for complex put");
    }
    idx = qlua_checklatcoord(L, 2);
    locked = QDP_expose_C(V->ptr);
    if (QDP_node_number(idx) == QDP_this_node) {
        QLA_Complex *zz = &QLA_elem_C(locked[QDP_index(idx)]);

        QLA_real(*zz) = z_re;
        QLA_imag(*zz) = z_im;
    }
    QDP_reset_C(V->ptr);
    qlua_free(L, idx);

    return 0;
}

static int
q_C_real(lua_State *L)
{
    mLatComplex *a = qlua_checkLatComplex(L, 1);
    mLatReal *c = qlua_newLatReal(L);

    QDP_R_eq_re_C(c->ptr, a->ptr, qCurrent);

    return 1;
}

static int
q_C_imag(lua_State *L)
{
    mLatComplex *a = qlua_checkLatComplex(L, 1);
    mLatReal *c = qlua_newLatReal(L);

    QDP_R_eq_im_C(c->ptr, a->ptr, qCurrent);

    return 1;
}

static int
q_C_add_C(lua_State *L)
{
    mLatComplex *a = qlua_checkLatComplex(L, 1);
    mLatComplex *b = qlua_checkLatComplex(L, 2);
    mLatComplex *c = qlua_newLatComplex(L);

    QDP_C_eq_C_plus_C(c->ptr, a->ptr, b->ptr, qCurrent);

    return 1;
}

static int
q_C_sub_C(lua_State *L)
{
    mLatComplex *a = qlua_checkLatComplex(L, 1);
    mLatComplex *b = qlua_checkLatComplex(L, 2);
    mLatComplex *c = qlua_newLatComplex(L);

    QDP_C_eq_C_minus_C(c->ptr, a->ptr, b->ptr, qCurrent);

    return 1;
}

static int
q_C_mul_C(lua_State *L)
{
    mLatComplex *a = qlua_checkLatComplex(L, 1);
    mLatComplex *b = qlua_checkLatComplex(L, 2);
    mLatComplex *c = qlua_newLatComplex(L);

    QDP_C_eq_C_times_C(c->ptr, a->ptr, b->ptr, qCurrent);

    return 1;
}

static int
q_C_mul_c(lua_State *L)
{
    mLatComplex *a = qlua_checkLatComplex(L, 1);
    QLA_Complex *b = qlua_checkComplex(L, 2);
    mLatComplex *c = qlua_newLatComplex(L);

    QDP_C_eq_c_times_C(c->ptr, b, a->ptr, qCurrent);

    return 1;
}

static int
q_c_mul_C(lua_State *L)
{
    QLA_Complex *a = qlua_checkComplex(L, 1);
    mLatComplex *b = qlua_checkLatComplex(L, 2);
    mLatComplex *c = qlua_newLatComplex(L);

    QDP_C_eq_c_times_C(c->ptr, a, b->ptr, qCurrent);

    return 1;
}


static int
q_C_mul_r(lua_State *L)
{
    mLatComplex *a = qlua_checkLatComplex(L, 1);
    QLA_Real b = luaL_checknumber(L, 2);
    mLatComplex *c = qlua_newLatComplex(L);

    QDP_C_eq_r_times_C(c->ptr, &b, a->ptr, qCurrent);

    return 1;
}

static int
q_r_mul_C(lua_State *L)
{
    QLA_Real a = luaL_checknumber(L, 1);
    mLatComplex *b = qlua_checkLatComplex(L, 2);
    mLatComplex *c = qlua_newLatComplex(L);

    QDP_C_eq_r_times_C(c->ptr, &a, b->ptr, qCurrent);

    return 1;
}

static int
q_C_div_C(lua_State *L)
{
    mLatComplex *a = qlua_checkLatComplex(L, 1);
    mLatComplex *b = qlua_checkLatComplex(L, 2);
    mLatComplex *c = qlua_newLatComplex(L);

    QDP_C_eq_C_divide_C(c->ptr, a->ptr, b->ptr, qCurrent);

    return 1;
}

static int
q_C_sum(lua_State *L)
{
    mLatComplex *a = qlua_checkLatComplex(L, 1);

    switch (lua_gettop(L)) {
    case 1: {
        QLA_Complex *s = qlua_newComplex(L);

        QDP_c_eq_sum_C(s, a->ptr, qCurrent);

        return 1;
    }
    case 2: {
        mLatMulti *m = qlua_checkLatMulti(L, 2);
        mVecComplex *r = qlua_newVecComplex(L, m->count);

        r->size = m->count;
        QDP_c_eq_sum_C_multi(r->val, a->ptr, m->subset, m->count);

        return 1;
    }
    }
    return luaL_error(L, "bad arguments for Complex:sum()");
}

static int
q_C_norm2(lua_State *L)
{
    mLatComplex *a = qlua_checkLatComplex(L, 1);
    QLA_Real n;

    QDP_r_eq_norm2_C(&n, a->ptr, qCurrent);
    lua_pushnumber(L, n);
    
    return 1;
}

static int
q_C_shift(lua_State *L)
{
    mLatComplex *a = qlua_checkLatComplex(L, 1);
    QDP_Shift shift = qlua_checkShift(L, 2);
    QDP_ShiftDir dir = qlua_checkShiftDir(L, 3);
    mLatComplex *r = qlua_newLatComplex(L);

    QDP_C_eq_sC(r->ptr, a->ptr, shift, dir, qCurrent);

    return 1;
}

static int
q_C_conj(lua_State *L)
{
    mLatComplex *a = qlua_checkLatComplex(L, 1);
    mLatComplex *r = qlua_newLatComplex(L);

    QDP_C_eq_conj_C(r->ptr, a->ptr, qCurrent);

    return 1;
}

static int
q_C_abs(lua_State *L)
{
    mLatComplex *a = qlua_checkLatComplex(L, 1);
    mLatReal *r = qlua_newLatReal(L);

    QDP_R_eq_norm_C(r->ptr, a->ptr, qCurrent);

    return 1;
}

static int
q_C_arg(lua_State *L)
{
    mLatComplex *a = qlua_checkLatComplex(L, 1);
    mLatReal *r = qlua_newLatReal(L);

    QDP_R_eq_arg_C(r->ptr, a->ptr, qCurrent);

    return 1;
}

static int
q_C_sqrt(lua_State *L)
{
    mLatComplex *a = qlua_checkLatComplex(L, 1);
    mLatComplex *r = qlua_newLatComplex(L);

    QDP_C_eq_csqrt_C(r->ptr, a->ptr, qCurrent);

    return 1;
}

static int
q_C_exp(lua_State *L)
{
    mLatComplex *a = qlua_checkLatComplex(L, 1);
    mLatComplex *r = qlua_newLatComplex(L);

    QDP_C_eq_cexp_C(r->ptr, a->ptr, qCurrent);

    return 1;
}

static int
q_C_log(lua_State *L)
{
    mLatComplex *a = qlua_checkLatComplex(L, 1);
    mLatComplex *r = qlua_newLatComplex(L);

    QDP_C_eq_clog_C(r->ptr, a->ptr, qCurrent);

    return 1;
}

static int
q_C_set(lua_State *L)
{
    mLatComplex *r = qlua_checkLatComplex(L, 1);
    mLatComplex *a = qlua_checkLatComplex(L, 2);

    QDP_C_eq_C(r->ptr, a->ptr, qCurrent);
    lua_pop(L, 1);

    return 1;
}

int
q_C_gaussian(lua_State *L)
{
    mLatRandom *a = qlua_checkLatRandom(L, 1);
    mLatComplex *r = qlua_newLatComplex(L);

    QDP_C_eq_gaussian_S(r->ptr, a->ptr, qCurrent);

    return 1;
}

static int
q_C_dot(lua_State *L)
{
    mLatComplex *a = qlua_checkLatComplex(L, 1);
    mLatComplex *b = qlua_checkLatComplex(L, 2);
    mLatComplex *s = qlua_newLatComplex(L);

    QDP_C_eq_C_dot_C(s->ptr, a->ptr, b->ptr, qCurrent);

    return 1;
}

static int
q_latcomplex(lua_State *L)
{
    switch (lua_gettop(L)) {
    case 2:
        switch (qlua_gettype(L, 2)) {
        case qReal: {
            QLA_Real d = luaL_checknumber(L, 2);
            QLA_Complex z;
            mLatComplex *v = qlua_newLatComplex(L);
            
            QLA_real(z) = d;
            QLA_imag(z) = 0.0;
            QDP_C_eq_c(v->ptr, &z, qCurrent);

            return 1;
        }
        case qComplex: {
            QLA_Complex *z = qlua_checkComplex(L, 2);
            mLatComplex *v = qlua_newLatComplex(L);

            QDP_C_eq_c(v->ptr, z, qCurrent);
            
            return 1;
        }
        case qLatReal: {
            mLatReal *d = qlua_checkLatReal(L, 2);
            mLatComplex *v = qlua_newLatComplex(L);
            
            QDP_C_eq_R(v->ptr, d->ptr, qCurrent);

            return 1;
        }
        case qLatComplex: {
            mLatComplex *d = qlua_checkLatComplex(L, 2);
            mLatComplex *v = qlua_newLatComplex(L);
            
            QDP_C_eq_C(v->ptr, d->ptr, qCurrent);

            return 1;
        }
        default:
            break;
        }
    case 3: {
        mLatReal *a = qlua_checkLatReal(L, 2);
        mLatReal *b = qlua_checkLatReal(L, 3);
        mLatComplex *c = qlua_newLatComplex(L);

        QDP_C_eq_R_plus_i_R(c->ptr, a->ptr, b->ptr, qCurrent);

        return 1;
    }
    }
    return qlua_badconstr(L, "Complex");
}

static struct luaL_Reg LatComplexMethods[] = {
    { "real",   q_C_real },
    { "imag",   q_C_imag },
    { "sum",    q_C_sum },
    { "norm2",  q_C_norm2 },
    { "shift",  q_C_shift },
    { "conj",   q_C_conj },
    { "abs",    q_C_abs },
    { "arg",    q_C_arg },
    { "sqrt",   q_C_sqrt },
    { "exp",    q_C_exp },
    { "log",    q_C_log },
    { "set",    q_C_set },
    { NULL,     NULL}
};

static struct luaL_Reg mtLatComplex[] = {
    { "__tostring",        q_C_fmt },
    { "__gc",              q_C_gc },
    { "__index",           q_C_get },
    { "__newindex",        q_C_put },
    { "__unm",             q_C_neg },
    { "__add",             qlua_add },
    { "__sub",             qlua_sub },
    { "__mul",             qlua_mul },
    { "__div",             qlua_div },
    { NULL,                NULL }
};

static struct luaL_Reg fLatComplex[] = {
    { "Complex",     q_latcomplex },
    { NULL,          NULL }
};

int
init_latcomplex(lua_State *L)
{
    luaL_getmetatable(L, opLattice);
    luaL_register(L, NULL, fLatComplex);
    lua_pop(L, 1);
    qlua_metatable(L, mtnLatComplex, mtLatComplex);
    qlua_metatable(L, opLatComplex, LatComplexMethods);
    qlua_reg_add(qLatComplex, qLatComplex, q_C_add_C);
    qlua_reg_sub(qLatComplex, qLatComplex, q_C_sub_C);
    qlua_reg_mul(qLatComplex, qLatComplex, q_C_mul_C);
    qlua_reg_mul(qLatComplex, qLatComplex, q_c_mul_C);
    qlua_reg_mul(qLatComplex, qLatComplex, q_C_mul_c);
    qlua_reg_mul(qReal,       qLatComplex, q_r_mul_C);
    qlua_reg_mul(qLatComplex, qReal,       q_C_mul_r);
    qlua_reg_div(qLatComplex, qLatComplex, q_C_div_C);
    qlua_reg_dot(qLatComplex, q_C_dot);

    return 0;
}

int
fini_latcomplex(lua_State *L)
{
    return 0;
}
