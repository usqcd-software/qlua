#include <qlua.h>
#include <latcomplex.h>
#include <latreal.h>
#include <latint.h>
#include <latrandom.h>
#include <qcomplex.h>
#include <qmp.h>

const char *mtnLatComplex = "qcd.lattice.complex";
static const char *opLatComplex = "qcd.lattice.complex.op";

mLatComplex *
qlua_newLatComplex(lua_State *L)
{
    QDP_D_Complex *v = QDP_D_create_C();
    mLatComplex *hdr;

    if (v == 0) {
        lua_gc(L, LUA_GCCOLLECT, 0);
        v = QDP_D_create_C();
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

    luaL_argcheck(L, v != 0, idx, "qcd.LatComplex expected");
    
    return v;
}

static int
qLatComplex_fmt(lua_State *L)
{
    char fmt[72];
    mLatComplex *b = qlua_checkLatComplex(L, 1);

    sprintf(fmt, "LatComplex(%p)", b->ptr);
    lua_pushstring(L, fmt);

    return 1;
}

static int
qLatComplex_gc(lua_State *L)
{
    mLatComplex *b = qlua_checkLatComplex(L, 1);

    QDP_D_destroy_C(b->ptr);
    b->ptr = 0;

    return 0;
}

static int
q_neg_C(lua_State *L)
{
    mLatComplex *a = qlua_checkLatComplex(L, 1);
    mLatComplex *r = qlua_newLatComplex(L);
    QLA_D_Real m1 = -1;

    QDP_D_C_eq_r_times_C(r->ptr, &m1, a->ptr, QDP_all);
    return 1;
}

static int
qLatComplex_get(lua_State *L)
{
    switch (qlua_gettype(L, 2)) {
    case qTable: {
        mLatComplex *V = qlua_checkLatComplex(L, 1);
        QLA_D_Complex *W;
        QLA_D_Complex *locked;
        double z_re;
        double z_im;
        int *idx = 0;
        
        idx = qlua_lattice_coord(L, 2);
        locked = QDP_D_expose_C(V->ptr);
        if (QDP_node_number(idx) == QDP_this_node) {
            QLA_D_Complex *zz = &QLA_D_elem_C(locked[QDP_index(idx)]);

            z_re = QLA_real(*zz);
            z_im = QLA_imag(*zz);
        } else {
            z_re = 0;
            z_im = 0;
        }
        QDP_D_reset_C(V->ptr);
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
    return luaL_error(L, "bad index");
}

static int
qLatComplex_put(lua_State *L)
{
    mLatComplex *V = qlua_checkLatComplex(L, 1);
    QLA_D_Complex *locked;
    int *idx = 0;
    double z_re, z_im;

    switch (qlua_gettype(L, 3)) {
    case qReal:
        z_re = luaL_checknumber(L, 3);
        break;
    case qComplex: {
        QLA_D_Complex *z = qlua_checkComplex(L, 3);
        z_re = QLA_real(*z);
        z_im = QLA_imag(*z);
        break;
    }
    default:
        return luaL_error(L, "bad argument");
    }
    idx = qlua_lattice_coord(L, 2);
    locked = QDP_D_expose_C(V->ptr);
    if (QDP_node_number(idx) == QDP_this_node) {
        QLA_D_Complex *zz = &QLA_D_elem_C(locked[QDP_index(idx)]);

        QLA_real(*zz) = z_re;
        QLA_imag(*zz) = z_im;
    }
    QDP_D_reset_C(V->ptr);
    qlua_free(L, idx);

    return 0;
}

static int
q_latcomplex(lua_State *L)
{
    switch (lua_gettop(L)) {
    case 1:
        switch (qlua_gettype(L, 1)) {
        case qReal: {
            QLA_D_Real d = luaL_checknumber(L, 1);
            QLA_D_Complex z;
            mLatComplex *v = qlua_newLatComplex(L);
            
            QLA_real(z) = d;
            QLA_imag(z) = 0.0;
            QDP_D_C_eq_c(v->ptr, &z, QDP_all);

            return 1;
        }
        case qLatReal: {
            mLatReal *d = qlua_checkLatReal(L, 1);
            mLatComplex *v = qlua_newLatComplex(L);
            
            QDP_D_C_eq_R(v->ptr, d->ptr, QDP_all);

            return 1;
        }
        case qLatComplex: {
            mLatComplex *d = qlua_checkLatComplex(L, 1);
            mLatComplex *v = qlua_newLatComplex(L);
            
            QDP_D_C_eq_C(v->ptr, d->ptr, QDP_all);

            return 1;
        }
        default:
            break;
        }
    case 2: {
        mLatReal *a = qlua_checkLatReal(L, 1);
        mLatReal *b = qlua_checkLatReal(L, 2);
        mLatComplex *c = qlua_newLatComplex(L);

        QDP_D_C_eq_R_plus_i_R(c->ptr, a->ptr, b->ptr, QDP_all);
        return 1;
    }
    }
    return luaL_error(L, "bad argument");
}

int
q_C_add_C(lua_State *L)
{
    mLatComplex *a = qlua_checkLatComplex(L, 1);
    mLatComplex *b = qlua_checkLatComplex(L, 2);
    mLatComplex *c = qlua_newLatComplex(L);

    QDP_D_C_eq_C_plus_C(c->ptr, a->ptr, b->ptr, QDP_all);

    return 1;
}

int
q_C_sub_C(lua_State *L)
{
    mLatComplex *a = qlua_checkLatComplex(L, 1);
    mLatComplex *b = qlua_checkLatComplex(L, 2);
    mLatComplex *c = qlua_newLatComplex(L);

    QDP_D_C_eq_C_minus_C(c->ptr, a->ptr, b->ptr, QDP_all);

    return 1;
}

int
q_C_mul_C(lua_State *L)
{
    mLatComplex *a = qlua_checkLatComplex(L, 1);
    mLatComplex *b = qlua_checkLatComplex(L, 2);
    mLatComplex *c = qlua_newLatComplex(L);

    QDP_D_C_eq_C_times_C(c->ptr, a->ptr, b->ptr, QDP_all);

    return 1;
}

int
q_C_mul_c(lua_State *L)
{
    mLatComplex *a = qlua_checkLatComplex(L, 1);
    QLA_D_Complex *b = qlua_checkComplex(L, 2);
    mLatComplex *c = qlua_newLatComplex(L);

    QDP_D_C_eq_c_times_C(c->ptr, b, a->ptr, QDP_all);

    return 1;
}

int
q_c_mul_C(lua_State *L)
{
    QLA_D_Complex *a = qlua_checkComplex(L, 1);
    mLatComplex *b = qlua_checkLatComplex(L, 2);
    mLatComplex *c = qlua_newLatComplex(L);

    QDP_D_C_eq_c_times_C(c->ptr, a, b->ptr, QDP_all);

    return 1;
}


int
q_C_mul_r(lua_State *L)
{
    mLatComplex *a = qlua_checkLatComplex(L, 1);
    QLA_D_Real b = luaL_checknumber(L, 2);
    mLatComplex *c = qlua_newLatComplex(L);

    QDP_D_C_eq_r_times_C(c->ptr, &b, a->ptr, QDP_all);

    return 1;
}

int
q_r_mul_C(lua_State *L)
{
    QLA_D_Real a = luaL_checknumber(L, 1);
    mLatComplex *b = qlua_checkLatComplex(L, 2);
    mLatComplex *c = qlua_newLatComplex(L);

    QDP_D_C_eq_r_times_C(c->ptr, &a, b->ptr, QDP_all);

    return 1;
}

int
q_C_div_C(lua_State *L)
{
    mLatComplex *a = qlua_checkLatComplex(L, 1);
    mLatComplex *b = qlua_checkLatComplex(L, 2);
    mLatComplex *c = qlua_newLatComplex(L);

    QDP_D_C_eq_C_divide_C(c->ptr, a->ptr, b->ptr, QDP_all);

    return 1;
}

static int
q_C_sum(lua_State *L)
{
    mLatComplex *a = qlua_checkLatComplex(L, 1);
    QLA_D_Complex *s = qlua_newComplex(L);

    QDP_D_c_eq_sum_C(s, a->ptr, QDP_all);

    return 1;
}

static int
q_C_norm2(lua_State *L)
{
    mLatComplex *a = qlua_checkLatComplex(L, 1);
    QLA_D_Real n;

    QDP_D_r_eq_norm2_C(&n, a->ptr, QDP_all);
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

    QDP_D_C_eq_sC(r->ptr, a->ptr, shift, dir, QDP_all);

    return 1;
}

static int
q_C_conj(lua_State *L)
{
    mLatComplex *a = qlua_checkLatComplex(L, 1);
    mLatComplex *r = qlua_newLatComplex(L);

    QDP_D_C_eq_conj_C(r->ptr, a->ptr, QDP_all);

    return 1;
}

static int
q_C_abs(lua_State *L)
{
    mLatComplex *a = qlua_checkLatComplex(L, 1);
    mLatReal *r = qlua_newLatReal(L);

    QDP_D_R_eq_norm_C(r->ptr, a->ptr, QDP_all);

    return 1;
}

static int
q_C_arg(lua_State *L)
{
    mLatComplex *a = qlua_checkLatComplex(L, 1);
    mLatReal *r = qlua_newLatReal(L);

    QDP_D_R_eq_arg_C(r->ptr, a->ptr, QDP_all);

    return 1;
}

static int
q_C_sqrt(lua_State *L)
{
    mLatComplex *a = qlua_checkLatComplex(L, 1);
    mLatComplex *r = qlua_newLatComplex(L);

    QDP_D_C_eq_csqrt_C(r->ptr, a->ptr, QDP_all);

    return 1;
}

static int
q_C_exp(lua_State *L)
{
    mLatComplex *a = qlua_checkLatComplex(L, 1);
    mLatComplex *r = qlua_newLatComplex(L);

    QDP_D_C_eq_cexp_C(r->ptr, a->ptr, QDP_all);

    return 1;
}

static int
q_C_log(lua_State *L)
{
    mLatComplex *a = qlua_checkLatComplex(L, 1);
    mLatComplex *r = qlua_newLatComplex(L);

    QDP_D_C_eq_clog_C(r->ptr, a->ptr, QDP_all);

    return 1;
}

int
q_C_dot(lua_State *L)
{
    mLatComplex *a = qlua_checkLatComplex(L, 1);
    mLatComplex *b = qlua_checkLatComplex(L, 2);
    QLA_D_Complex *s = qlua_newComplex(L);

    QDP_D_c_eq_C_dot_C(s, a->ptr, b->ptr, QDP_all);

    return 1;
}

static struct luaL_Reg LatComplexMethods[] = {
    { "sum",    q_C_sum },
    { "norm2",  q_C_norm2 },
    { "shift",  q_C_shift },
    { "conj",   q_C_conj },
    { "abs",    q_C_abs },
    { "arg",    q_C_arg },
    { "sqrt",   q_C_sqrt },
    { "exp",    q_C_exp },
    { "log",    q_C_log },
    { NULL,     NULL}
};

static struct luaL_Reg mtLatComplex[] = {
    { "__tostring",        qLatComplex_fmt },
    { "__gc",              qLatComplex_gc },
    { "__index",           qLatComplex_get },
    { "__newindex",        qLatComplex_put },
    { "__umn",             q_neg_C },
    { "__add",             qlua_add },
    { "__sub",             qlua_sub },
    { "__mul",             qlua_mul },
    { "__div",             qlua_div },
    { NULL,                NULL }
};

static struct luaL_Reg fLatComplex[] = {
    { "lat_complex", q_latcomplex },
    { NULL,          NULL }
};

int
init_latcomplex(lua_State *L)
{
    luaL_register(L, qcdlib, fLatComplex);
    qlua_metatable(L, mtnLatComplex, mtLatComplex);
    qlua_metatable(L, opLatComplex, LatComplexMethods);

    return 0;
}

int
fini_latcomplex(lua_State *L)
{
    return 0;
}
