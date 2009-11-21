#include <qlua.h>                                                    /* DEPS */
#include <qvector.h>                                                 /* DEPS */
#include <lattice.h>                                                 /* DEPS */
#include <latreal.h>                                                 /* DEPS */
#include <latint.h>                                                  /* DEPS */
#include <latrandom.h>                                               /* DEPS */
#include <latcomplex.h>                                              /* DEPS */
#include <latmulti.h>                                                /* DEPS */
#include <qmp.h>

const char mtnLatReal[] = "lattice.Real";
static const char opLatReal[] = "lattice.Real.ops";

mLatReal *
qlua_newLatReal(lua_State *L)
{
    QDP_Real *v = QDP_create_R();
    mLatReal *hdr;

    if (v == 0) {
        lua_gc(L, LUA_GCCOLLECT, 0);
        v = QDP_create_R();
        if (v == 0)
            luaL_error(L, "not enough memory (QDP_Real)");
    }
    hdr = lua_newuserdata(L, sizeof (mLatReal));
    hdr->ptr = v;
    luaL_getmetatable(L, mtnLatReal);
    lua_setmetatable(L, -2);

    return hdr;
}

mLatReal *
qlua_checkLatReal(lua_State *L, int idx)
{
    void *v = luaL_checkudata(L, idx, mtnLatReal);

    luaL_argcheck(L, v != 0, idx, "lattice.Real expected");

    return v;
}

static int
q_R_fmt(lua_State *L)
{
    char fmt[72];
    mLatReal *b = qlua_checkLatReal(L, 1);

    sprintf(fmt, "LatReal(%p)", b->ptr);
    lua_pushstring(L, fmt);

    return 1;
}

static int
q_R_gc(lua_State *L)
{
    mLatReal *b = qlua_checkLatReal(L, 1);

    QDP_destroy_R(b->ptr);
    b->ptr = 0;

    return 0;
}

static int
q_R_neg(lua_State *L)
{
    mLatReal *a = qlua_checkLatReal(L, 1);
    mLatReal *res = qlua_newLatReal(L);
    QLA_Real m1 = -1;

    CALL_QDP(L);
    QDP_R_eq_r_times_R(res->ptr, &m1, a->ptr, *qCurrent);
    return 1;
}

static int
q_R_sum(lua_State *L)
{
    mLatReal *a = qlua_checkLatReal(L, 1);
    
    switch (lua_gettop(L)) {
    case 1: {
        QLA_Real sum;

        CALL_QDP(L);
        QDP_r_eq_sum_R(&sum, a->ptr, *qCurrent);
        lua_pushnumber(L, sum);

        return 1;
    }
    case 2: {
        mLatMulti *m = qlua_checkLatMulti(L, 2);
        int size = m->size;
        QLA_Int *ii = m->idx;
        mVecReal *r = qlua_newVecReal(L, size);
        int k;
        QLA_Real *xx;
        
        r->size = size;
        for (k = 0; k < size; k++)
            r->val[k] = 0;

        CALL_QDP(L);
        xx = QDP_expose_R(a->ptr);
        for (k = 0; k < QDP_sites_on_node; k++, xx++, ii++) {
            int t = *ii;
            if ((t < 0) || (t >= size))
                continue;
            r->val[t] += *xx;
        }
        QDP_reset_R(a->ptr);
        QMP_sum_double_array(r->val, size);

        return 1;
    }
    }
    return luaL_error(L, "bad arguments for Real:sum()");
}

static int
q_R_norm2(lua_State *L)
{
    mLatReal *a = qlua_checkLatReal(L, 1);
    QLA_Real n;

    CALL_QDP(L);
    QDP_r_eq_norm2_R(&n, a->ptr, *qCurrent);
    lua_pushnumber(L, n);

    return 1;
}

static int
q_R_shift(lua_State *L)
{
    mLatReal *a = qlua_checkLatReal(L, 1);
    QDP_Shift shift = qlua_checkShift(L, 2);
    QDP_ShiftDir dir = qlua_checkShiftDir(L, 3);
    mLatReal *r = qlua_newLatReal(L);

    CALL_QDP(L);
    QDP_R_eq_sR(r->ptr, a->ptr, shift, dir, *qCurrent);

    return 1;
}

static int
q_R_sin(lua_State *L)
{
    mLatReal *a = qlua_checkLatReal(L, 1);
    mLatReal *r = qlua_newLatReal(L);

    CALL_QDP(L);
    QDP_R_eq_sin_R(r->ptr, a->ptr, *qCurrent);

    return 1;
}

static int
q_R_cos(lua_State *L)
{
    mLatReal *a = qlua_checkLatReal(L, 1);
    mLatReal *r = qlua_newLatReal(L);

    CALL_QDP(L);
    QDP_R_eq_cos_R(r->ptr, a->ptr, *qCurrent);

    return 1;
}

static int
q_R_tan(lua_State *L)
{
    mLatReal *a = qlua_checkLatReal(L, 1);
    mLatReal *r = qlua_newLatReal(L);

    CALL_QDP(L);
    QDP_R_eq_tan_R(r->ptr, a->ptr, *qCurrent);

    return 1;
}

static int
q_R_asin(lua_State *L)
{
    mLatReal *a = qlua_checkLatReal(L, 1);
    mLatReal *r = qlua_newLatReal(L);

    CALL_QDP(L);
    QDP_R_eq_asin_R(r->ptr, a->ptr, *qCurrent);

    return 1;
}

static int
q_R_acos(lua_State *L)
{
    mLatReal *a = qlua_checkLatReal(L, 1);
    mLatReal *r = qlua_newLatReal(L);

    CALL_QDP(L);
    QDP_R_eq_acos_R(r->ptr, a->ptr, *qCurrent);

    return 1;
}

static int
q_R_atan(lua_State *L)
{
    mLatReal *a = qlua_checkLatReal(L, 1);
    mLatReal *r = qlua_newLatReal(L);

    CALL_QDP(L);
    QDP_R_eq_atan_R(r->ptr, a->ptr, *qCurrent);

    return 1;
}

static int
q_R_sqrt(lua_State *L)
{
    mLatReal *a = qlua_checkLatReal(L, 1);
    mLatReal *r = qlua_newLatReal(L);

    CALL_QDP(L);
    QDP_R_eq_sqrt_R(r->ptr, a->ptr, *qCurrent);

    return 1;
}

static int
q_R_abs(lua_State *L)
{
    mLatReal *a = qlua_checkLatReal(L, 1);
    mLatReal *r = qlua_newLatReal(L);

    CALL_QDP(L);
    QDP_R_eq_fabs_R(r->ptr, a->ptr, *qCurrent);

    return 1;
}

static int
q_R_exp(lua_State *L)
{
    mLatReal *a = qlua_checkLatReal(L, 1);
    mLatReal *r = qlua_newLatReal(L);

    CALL_QDP(L);
    QDP_R_eq_exp_R(r->ptr, a->ptr, *qCurrent);

    return 1;
}

static int
q_R_expi(lua_State *L)
{
    mLatReal *a = qlua_checkLatReal(L, 1);
    mLatComplex *r = qlua_newLatComplex(L);

    CALL_QDP(L);
    QDP_C_eq_cexpi_R(r->ptr, a->ptr, *qCurrent);

    return 1;
}

static int
q_R_log(lua_State *L)
{
    mLatReal *a = qlua_checkLatReal(L, 1);
    mLatReal *r = qlua_newLatReal(L);

    CALL_QDP(L);
    QDP_R_eq_log_R(r->ptr, a->ptr, *qCurrent);

    return 1;
}

static int
q_R_sign(lua_State *L)
{
    mLatReal *a = qlua_checkLatReal(L, 1);
    mLatReal *r = qlua_newLatReal(L);

    CALL_QDP(L);
    QDP_R_eq_sign_R(r->ptr, a->ptr, *qCurrent);

    return 1;
}

static int
q_R_floor(lua_State *L)
{
    mLatReal *a = qlua_checkLatReal(L, 1);
    mLatReal *r = qlua_newLatReal(L);

    CALL_QDP(L);
    QDP_R_eq_floor_R(r->ptr, a->ptr, *qCurrent);

    return 1;
}

static int
q_R_ceil(lua_State *L)
{
    mLatReal *a = qlua_checkLatReal(L, 1);
    mLatReal *r = qlua_newLatReal(L);

    CALL_QDP(L);
    QDP_R_eq_ceil_R(r->ptr, a->ptr, *qCurrent);

    return 1;
}

static int
q_R_sinh(lua_State *L)
{
    mLatReal *a = qlua_checkLatReal(L, 1);
    mLatReal *r = qlua_newLatReal(L);

    CALL_QDP(L);
    QDP_R_eq_sinh_R(r->ptr, a->ptr, *qCurrent);

    return 1;
}

static int
q_R_cosh(lua_State *L)
{
    mLatReal *a = qlua_checkLatReal(L, 1);
    mLatReal *r = qlua_newLatReal(L);

    CALL_QDP(L);
    QDP_R_eq_cosh_R(r->ptr, a->ptr, *qCurrent);

    return 1;
}

static int
q_R_tanh(lua_State *L)
{
    mLatReal *a = qlua_checkLatReal(L, 1);
    mLatReal *r = qlua_newLatReal(L);

    CALL_QDP(L);
    QDP_R_eq_tanh_R(r->ptr, a->ptr, *qCurrent);

    return 1;
}

static int
q_R_log10(lua_State *L)
{
    mLatReal *a = qlua_checkLatReal(L, 1);
    mLatReal *r = qlua_newLatReal(L);

    CALL_QDP(L);
    QDP_R_eq_log10_R(r->ptr, a->ptr, *qCurrent);

    return 1;
}

int
q_R_random(lua_State *L)
{
    mLatRandom *a = qlua_checkLatRandom(L, 1);
    mLatReal *r = qlua_newLatReal(L);

    CALL_QDP(L);
    QDP_R_eq_random_S(r->ptr, a->ptr, *qCurrent);

    return 1;
}

int
q_R_gaussian(lua_State *L)
{
    mLatRandom *a = qlua_checkLatRandom(L, 1);
    mLatReal *r = qlua_newLatReal(L);

    CALL_QDP(L);
    QDP_R_eq_gaussian_S(r->ptr, a->ptr, *qCurrent);

    return 1;
}

static int
q_R_trunc(lua_State *L)
{
    mLatReal *a = qlua_checkLatReal(L, 1);
    mLatInt *r = qlua_newLatInt(L);

    CALL_QDP(L);
    QDP_I_eq_trunc_R(r->ptr, a->ptr, *qCurrent);

    return 1;
}

static int
q_R_round(lua_State *L)
{
    mLatReal *a = qlua_checkLatReal(L, 1);
    mLatInt *r = qlua_newLatInt(L);

    CALL_QDP(L);
    QDP_I_eq_round_R(r->ptr, a->ptr, *qCurrent);

    return 1;
}

static int
q_R_set(lua_State *L)
{
    mLatReal *r = qlua_checkLatReal(L, 1);
    mLatReal *a = qlua_checkLatReal(L, 2);

    CALL_QDP(L);
    QDP_R_eq_R(r->ptr, a->ptr, *qCurrent);
    lua_pop(L, 1);

    return 1;
}

static int
q_R_get(lua_State *L)
{
    switch (qlua_gettype(L, 2)) {
    case qTable: {
        mLatReal *V = qlua_checkLatReal(L, 1);
        QLA_Real *locked;
        int *idx = 0;
        double z;

        idx = qlua_checklatcoord(L, 2);
        CALL_QDP(L);
        locked = QDP_expose_R(V->ptr);
        if (QDP_node_number(idx) == QDP_this_node) {
            z = QLA_elem_R(locked[QDP_index(idx)]);
        } else {
            z = 0;
        }
        QDP_reset_R(V->ptr);
        qlua_free(L, idx);
        QMP_sum_double(&z);
        lua_pushnumber(L, z);

        return 1;
    }
    case qString:
        return qlua_lookup(L, 2, opLatReal);
    }
    return qlua_badindex(L, "Real");
}


static int
q_R_put(lua_State *L)
{
    mLatReal *V = qlua_checkLatReal(L, 1);
    QLA_Real *locked;
    int *idx = 0;
    double z = luaL_checknumber(L, 3);

    idx = qlua_checklatcoord(L, 2);
    CALL_QDP(L);
    locked = QDP_expose_R(V->ptr);
    if (QDP_node_number(idx) == QDP_this_node) {
        QLA_elem_R(locked[QDP_index(idx)]) = z;
    }
    QDP_reset_R(V->ptr);
    qlua_free(L, idx);

    return 0;
}

static int
q_R_add_R(lua_State *L)
{
    mLatReal *a = qlua_checkLatReal(L, 1);
    mLatReal *b = qlua_checkLatReal(L, 2);
    mLatReal *c = qlua_newLatReal(L);

    CALL_QDP(L);
    QDP_R_eq_R_plus_R(c->ptr, a->ptr, b->ptr, *qCurrent);

    return 1;
}

static struct {
    QLA_Real a;
    QLA_Real *b;
} Rop_args; /* YYY global state */

static void
do_rRadd(QLA_Real *r, int idx)
{
    *r = Rop_args.a + Rop_args.b[idx];
}

static void
X_R_eq_r_plus_R(QDP_Real *r, double a, QDP_Real *b, QDP_Subset s)
{
    Rop_args.a = a;
    Rop_args.b = QDP_expose_R(b);
    QDP_R_eq_funci(r, do_rRadd, s);
    QDP_reset_R(b);
    Rop_args.a = 0;
    Rop_args.b = 0;
}

static int
q_r_add_R(lua_State *L)
{
    double a = luaL_checknumber(L, 1);
    mLatReal *b = qlua_checkLatReal(L, 2);
    mLatReal *c = qlua_newLatReal(L);

    CALL_QDP(L);
    X_R_eq_r_plus_R(c->ptr, a, b->ptr, *qCurrent);

    return 1;
}

static int
q_R_add_r(lua_State *L)
{
    mLatReal *a = qlua_checkLatReal(L, 1);
    double b = luaL_checknumber(L, 2);
    mLatReal *c = qlua_newLatReal(L);

    CALL_QDP(L);
    X_R_eq_r_plus_R(c->ptr, b, a->ptr, *qCurrent);

    return 1;
}

static int
q_R_sub_R(lua_State *L)
{
    mLatReal *a = qlua_checkLatReal(L, 1);
    mLatReal *b = qlua_checkLatReal(L, 2);
    mLatReal *c = qlua_newLatReal(L);

    CALL_QDP(L);
    QDP_R_eq_R_minus_R(c->ptr, a->ptr, b->ptr, *qCurrent);

    return 1;
}

static void
do_rRsub(QLA_Real *r, int idx)
{
    *r = Rop_args.a - Rop_args.b[idx];
}

static void
X_R_eq_r_minus_R(QDP_Real *r, double a, QDP_Real *b, QDP_Subset s)
{
    Rop_args.a = a;
    Rop_args.b = QDP_expose_R(b);
    QDP_R_eq_funci(r, do_rRsub, s);
    QDP_reset_R(b);
    Rop_args.a = 0;
    Rop_args.b = 0;
}

static int
q_r_sub_R(lua_State *L)
{
    double a = luaL_checknumber(L, 1);
    mLatReal *b = qlua_checkLatReal(L, 2);
    mLatReal *c = qlua_newLatReal(L);

    CALL_QDP(L);
    X_R_eq_r_minus_R(c->ptr, a, b->ptr, *qCurrent);

    return 1;
}

static int
q_R_sub_r(lua_State *L)
{
    mLatReal *a = qlua_checkLatReal(L, 1);
    double b = luaL_checknumber(L, 2);
    mLatReal *c = qlua_newLatReal(L);

    CALL_QDP(L);
    X_R_eq_r_plus_R(c->ptr, -b, a->ptr, *qCurrent);

    return 1;
}

static int
q_R_mul_R(lua_State *L)
{
    mLatReal *a = qlua_checkLatReal(L, 1);
    mLatReal *b = qlua_checkLatReal(L, 2);
    mLatReal *c = qlua_newLatReal(L);

    CALL_QDP(L);
    QDP_R_eq_R_times_R(c->ptr, a->ptr, b->ptr, *qCurrent);

    return 1;
}

static int
q_r_mul_R(lua_State *L)
{
    QLA_Real b = luaL_checknumber(L, 1);
    mLatReal *a = qlua_checkLatReal(L, 2);
    mLatReal *c = qlua_newLatReal(L);

    CALL_QDP(L);
    QDP_R_eq_r_times_R(c->ptr, &b, a->ptr, *qCurrent);

    return 1;
}

static int
q_R_mul_r(lua_State *L)
{
    mLatReal *a = qlua_checkLatReal(L, 1);
    QLA_Real b = luaL_checknumber(L, 2);
    mLatReal *c = qlua_newLatReal(L);

    CALL_QDP(L);
    QDP_R_eq_r_times_R(c->ptr, &b, a->ptr, *qCurrent);

    return 1;
}

static int
q_R_div_R(lua_State *L)
{
    mLatReal *a = qlua_checkLatReal(L, 1);
    mLatReal *b = qlua_checkLatReal(L, 2);
    mLatReal *c = qlua_newLatReal(L);

    CALL_QDP(L);
    QDP_R_eq_R_divide_R(c->ptr, a->ptr, b->ptr, *qCurrent);

    return 1;
}

static int
q_R_div_r(lua_State *L)
{
    mLatReal *a = qlua_checkLatReal(L, 1);
    double b = 1/luaL_checknumber(L, 2);
    mLatReal *c = qlua_newLatReal(L);

    CALL_QDP(L);
    QDP_R_eq_r_times_R(c->ptr, &b, a->ptr, *qCurrent);

    return 1;
}

static void
do_rRdiv(QLA_Real *r, int idx)
{
    *r = Rop_args.a / Rop_args.b[idx];
}

static void
X_R_eq_r_divide_R(QDP_Real *r, double a, QDP_Real *b, QDP_Subset s)
{
    Rop_args.a = a;
    Rop_args.b = QDP_expose_R(b);
    QDP_R_eq_funci(r, do_rRdiv, s);
    QDP_reset_R(b);
    Rop_args.a = 0;
    Rop_args.b = 0;
}

static int
q_r_div_R(lua_State *L)
{
    double a = luaL_checknumber(L, 1);
    mLatReal *b = qlua_checkLatReal(L, 2);
    mLatReal *c = qlua_newLatReal(L);

    CALL_QDP(L);
    X_R_eq_r_divide_R(c->ptr, a, b->ptr, *qCurrent);

    return 1;
}

static int
q_latreal(lua_State *L)
{
    switch (lua_gettop(L)) {
    case 1: {
        mLatReal *v = qlua_newLatReal(L);

        CALL_QDP(L);
        QDP_R_eq_zero(v->ptr, *qCurrent);

        return 1;
    }
    case 2:
        switch (qlua_gettype(L, 2)) {
        case qReal: {
            QLA_Real d = luaL_checknumber(L, 2);
            mLatReal *v = qlua_newLatReal(L);
            
            CALL_QDP(L);
            QDP_R_eq_r(v->ptr, &d, *qCurrent);
            
            return 1;
        }
        case qLatInt: {
            mLatInt *d = qlua_checkLatInt(L, 2);
            mLatReal *v = qlua_newLatReal(L);
            
            CALL_QDP(L);
            QDP_R_eq_I(v->ptr, d->ptr, *qCurrent);
            
            return 1;
        }
        case qLatReal: {
            mLatReal *d = qlua_checkLatReal(L, 2);
            mLatReal *v = qlua_newLatReal(L);
            
            CALL_QDP(L);
            QDP_R_eq_R(v->ptr, d->ptr, *qCurrent);
            
            return 1;
        }
        }
    }
    return qlua_badconstr(L, "Real");
}

static const struct luaL_Reg LatRealMethods[] = {
    { "sum",       q_R_sum      },
    { "norm2",     q_R_norm2    },
    { "shift",     q_R_shift    },
    { "sin",       q_R_sin      },
    { "cos",       q_R_cos      },
    { "tan",       q_R_tan      },
    { "asin",      q_R_asin     },
    { "acos",      q_R_acos     },
    { "atan",      q_R_atan     },
    { "sqrt",      q_R_sqrt     },
    { "abs",       q_R_abs      },
    { "exp",       q_R_exp      },
    { "log",       q_R_log      },
    { "sign",      q_R_sign     },
    { "ceil",      q_R_ceil     },
    { "floor",     q_R_floor    },
    { "sinh",      q_R_sinh     },
    { "cosh",      q_R_cosh     },
    { "tanh",      q_R_tanh     },
    { "log10",     q_R_log10    },
    { "expi",      q_R_expi     },
    { "trunc",     q_R_trunc    },
    { "round",     q_R_round    },
    { "set",       q_R_set      },
    { NULL,        NULL         }
};

static struct luaL_Reg mtLatReal[] = {
    { "__tostring",   q_R_fmt },
    { "__gc",         q_R_gc },
    { "__index",      q_R_get },
    { "__newindex",   q_R_put },
    { "__unm",        q_R_neg },
    { "__add",        qlua_add },
    { "__sub",        qlua_sub },
    { "__mul",        qlua_mul },
    { "__div",        qlua_div },
    { NULL,           NULL }
};

static struct luaL_Reg fLatReal[] = {
    { "Real",   q_latreal },
    { NULL,         NULL }
};

int
init_latreal(lua_State *L)
{
    luaL_getmetatable(L, opLattice);
    luaL_register(L, NULL, fLatReal);
    lua_pop(L, 1);
    qlua_metatable(L, mtnLatReal, mtLatReal);
    qlua_metatable(L, opLatReal, LatRealMethods);

    qlua_reg_add(qLatReal, qLatReal, q_R_add_R);
    qlua_reg_add(qReal,    qLatReal, q_r_add_R);
    qlua_reg_add(qLatReal, qReal,    q_R_add_r);
    qlua_reg_sub(qLatReal, qLatReal, q_R_sub_R);
    qlua_reg_sub(qReal,    qLatReal, q_r_sub_R);
    qlua_reg_sub(qLatReal, qReal,    q_R_sub_r);
    qlua_reg_mul(qLatReal, qLatReal, q_R_mul_R);
    qlua_reg_mul(qReal,    qLatReal, q_r_mul_R);
    qlua_reg_mul(qLatReal, qReal,    q_R_mul_r);
    qlua_reg_div(qLatReal, qLatReal, q_R_div_R);
    qlua_reg_div(qReal,    qLatReal, q_r_div_R);
    qlua_reg_div(qLatReal, qReal,    q_R_div_r);
    qlua_reg_dot(qLatReal, q_R_mul_R);

    return 0;
}

int
fini_latreal(lua_State *L)
{
    return 0;
}
