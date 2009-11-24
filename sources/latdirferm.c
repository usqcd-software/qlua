#include <qlua.h>                                                    /* DEPS */
#include <qcomplex.h>                                                /* DEPS */
#include <lattice.h>                                                 /* DEPS */
#include <latint.h>                                                  /* DEPS */
#include <latrandom.h>                                               /* DEPS */
#include <latreal.h>                                                 /* DEPS */
#include <latcomplex.h>                                              /* DEPS */
#include <latcolvec.h>                                               /* DEPS */
#include <latcolmat.h>                                               /* DEPS */
#include <latdirferm.h>                                              /* DEPS */
#include <qmp.h>

const char mtnLatDirFerm[] = "lattice.DiracFermion";
static char opLatDirFerm[] = "lattice.DiracFermion.ops";

mLatDirFerm *
qlua_newLatDirFerm(lua_State *L)
{
    QDP_DiracFermion *v = QDP_create_D();
    mLatDirFerm *hdr;

    if (v == 0) {
        lua_gc(L, LUA_GCCOLLECT, 0);
        v = QDP_create_D();
        if (v == 0)
            luaL_error(L, "not enough memory (QDP_ColorMatrix)");
    }
    hdr = lua_newuserdata(L, sizeof (mLatDirFerm));
    hdr->ptr = v;
    luaL_getmetatable(L, mtnLatDirFerm);
    lua_setmetatable(L, -2);

    return hdr;
}

mLatDirFerm *
qlua_checkLatDirFerm(lua_State *L, int idx)
{
    void *v = luaL_checkudata(L, idx, mtnLatDirFerm);
    
    luaL_argcheck(L, v != 0, idx, "lattice.DiracFermion expected");
    
    return v;
}

static int
q_D_fmt(lua_State *L)
{
    char fmt[72];
    mLatDirFerm *b = qlua_checkLatDirFerm(L, 1);

    sprintf(fmt, "QDP:DiracFermion(%p)", b->ptr);
    lua_pushstring(L, fmt);

    return 1;
}

static int
q_D_gc(lua_State *L)
{
    mLatDirFerm *b = qlua_checkLatDirFerm(L, 1);

    QDP_destroy_D(b->ptr);
    b->ptr = 0;

    return 0;
}

static int
q_D_get(lua_State *L)
{
    switch (qlua_gettype(L, 2)) {
    case qTable: {
        mLatDirFerm *V = qlua_checkLatDirFerm(L, 1);
        int d = qlua_checkdiracindex(L, 2);
        int c = qlua_colorindex(L, 2);
        int *idx = qlua_latcoord(L, 2);
        if (idx == NULL) {
            if (c == -1) {
                mLatColVec *r = qlua_newLatColVec(L);

                CALL_QDP(L);
                QDP_V_eq_colorvec_D(r->ptr, V->ptr, d, *qCurrent);
            } else {
                mLatComplex *r = qlua_newLatComplex(L);

                CALL_QDP(L);
                QDP_C_eq_elem_D(r->ptr, V->ptr, c, d, *qCurrent);
            }
        } else {
            if (c == -1) {
                qlua_free(L, idx);
                return qlua_badindex(L, "DiracFermion");
            } else {
                QLA_Complex *W = qlua_newComplex(L);
                QLA_DiracFermion *locked;
                double zri[2];

                qlua_verifylatcoord(L, idx);
                CALL_QDP(L);
                locked = QDP_expose_D(V->ptr);
                if (QDP_node_number(idx) == QDP_this_node) {
                    QLA_Complex *zz = &QLA_elem_D(locked[QDP_index(idx)], c, d);
                    zri[0] = QLA_real(*zz);
                    zri[1] = QLA_imag(*zz);
                } else {
                    zri[0] = 0;
                    zri[1] = 0;
                }
                QDP_reset_D(V->ptr);
                QMP_sum_double_array(zri, 2);
                QLA_c_eq_r_plus_ir(*W, zri[0], zri[1]);

            }
        }
        qlua_free(L, idx);
        return 1;
    }
    case qString:
        return qlua_lookup(L, 2, opLatDirFerm);
    }
    return qlua_badindex(L, "DiracFermion");
}

static int
q_D_put(lua_State *L)
{
    mLatDirFerm *V = qlua_checkLatDirFerm(L, 1);
    int d = qlua_checkdiracindex(L, 2);
    int c = qlua_colorindex(L, 2);
    int *idx = qlua_latcoord(L, 2);

    if (idx == NULL) {
        if (c == -1) {
            mLatColVec *z = qlua_checkLatColVec(L, 3);

            CALL_QDP(L);
            QDP_D_eq_colorvec_V(V->ptr, z->ptr, d, *qCurrent);
        } else {
            mLatComplex *z = qlua_checkLatComplex(L, 3);

            CALL_QDP(L);
            QDP_D_eq_elem_C(V->ptr, z->ptr, c, d, *qCurrent);
        }
    } else {
        if (c == -1) {
            qlua_free(L, idx);
            return qlua_badindex(L, "DiracFermion");
        } else {
            QLA_Complex *z = qlua_checkComplex(L, 3);
            qlua_verifylatcoord(L, idx);
            CALL_QDP(L);
            if (QDP_node_number(idx) == QDP_this_node) {
                QLA_DiracFermion *locked = QDP_expose_D(V->ptr);
                QLA_Complex *zz = &QLA_elem_D(locked[QDP_index(idx)], c, d);

                QLA_c_eq_c(*zz, *z);
                QDP_reset_D(V->ptr);
            }
        }
    }
    qlua_free(L, idx);

    return 0;
}

static int
q_D_neg(lua_State *L)
{
    mLatDirFerm *a = qlua_checkLatDirFerm(L, 1);
    mLatDirFerm *r = qlua_newLatDirFerm(L);
    QLA_Real m1 = -1;

    CALL_QDP(L);
    QDP_D_eq_r_times_D(r->ptr, &m1, a->ptr, *qCurrent);

    return 1;
}

static int
q_D_norm2(lua_State *L)
{
    mLatDirFerm *a = qlua_checkLatDirFerm(L, 1);
    QLA_Real n;

    CALL_QDP(L);
    QDP_r_eq_norm2_D(&n, a->ptr, *qCurrent);
    lua_pushnumber(L, n);
    
    return 1;
}

static int
q_D_shift(lua_State *L)
{
    mLatDirFerm *a = qlua_checkLatDirFerm(L, 1);
    QDP_Shift shift = qlua_checkShift(L, 2);
    QDP_ShiftDir dir = qlua_checkShiftDir(L, 3);
    mLatDirFerm *r = qlua_newLatDirFerm(L);

    CALL_QDP(L);
    QDP_D_eq_sD(r->ptr, a->ptr, shift, dir, *qCurrent);

    return 1;
}

static int
q_D_conj(lua_State *L)
{
    mLatDirFerm *a = qlua_checkLatDirFerm(L, 1);
    mLatDirFerm *r = qlua_newLatDirFerm(L);

    CALL_QDP(L);
    QDP_D_eq_conj_D(r->ptr, a->ptr, *qCurrent);

    return 1;
}

static int
q_D_gamma(lua_State *L)
{
    mLatDirFerm *f = qlua_checkLatDirFerm(L, 1);
    int mu = qlua_gammaindex(L, 2);
    int n = qlua_gammabinary(L, 2);
    
    if (((n == -1) && (mu == -1)) || ((n != -1) && (mu != -1)))
        return qlua_badindex(L, "DiracFermion");
    if (n == -1) {
        mLatDirFerm *r = qlua_newLatDirFerm(L);

        CALL_QDP(L);
        if (mu < 5) {
            QDP_D_eq_gamma_times_D(r->ptr, f->ptr, 1 << mu, *qCurrent);
        } else {
            QDP_D_eq_gamma_times_D(r->ptr, f->ptr, 15, *qCurrent);
        }

        return 1;
    }
    if (mu == -1) {
        mLatDirFerm *r = qlua_newLatDirFerm(L);

        CALL_QDP(L);
        QDP_D_eq_gamma_times_D(r->ptr, f->ptr, n, *qCurrent);
        
        return 1;
    }

    return luaL_error(L, "this could never happen");
}

static int
q_D_set(lua_State *L)
{
    mLatDirFerm *r = qlua_checkLatDirFerm(L, 1);
    mLatDirFerm *a = qlua_checkLatDirFerm(L, 2);

    CALL_QDP(L);
    QDP_D_eq_D(r->ptr, a->ptr, *qCurrent);
    lua_pop(L, 1);

    return 1;
}

static int
q_D_dot(lua_State *L)
{
    mLatDirFerm *a = qlua_checkLatDirFerm(L, 1);
    mLatDirFerm *b = qlua_checkLatDirFerm(L, 2);
    mLatComplex *s = qlua_newLatComplex(L);

    CALL_QDP(L);
    QDP_C_eq_D_dot_D(s->ptr, a->ptr, b->ptr, *qCurrent);

    return 1;
}

int
q_D_gaussian(lua_State *L)
{
    mLatRandom *a = qlua_checkLatRandom(L, 1);
    mLatDirFerm *r = qlua_newLatDirFerm(L);

    CALL_QDP(L);
    QDP_D_eq_gaussian_S(r->ptr, a->ptr, *qCurrent);

    return 1;
}


static int
q_D_add_D(lua_State *L)
{
    mLatDirFerm *a = qlua_checkLatDirFerm(L, 1);
    mLatDirFerm *b = qlua_checkLatDirFerm(L, 2);
    mLatDirFerm *c = qlua_newLatDirFerm(L);

    CALL_QDP(L);
    QDP_D_eq_D_plus_D(c->ptr, a->ptr, b->ptr, *qCurrent);

    return 1;
}

static int
q_D_sub_D(lua_State *L)
{
    mLatDirFerm *a = qlua_checkLatDirFerm(L, 1);
    mLatDirFerm *b = qlua_checkLatDirFerm(L, 2);
    mLatDirFerm *c = qlua_newLatDirFerm(L);

    CALL_QDP(L);
    QDP_D_eq_D_minus_D(c->ptr, a->ptr, b->ptr, *qCurrent);

    return 1;
}

static int
q_D_mul_r(lua_State *L)
{
    mLatDirFerm *a = qlua_checkLatDirFerm(L, 1);
    QLA_Real b = luaL_checknumber(L, 2);
    mLatDirFerm *c = qlua_newLatDirFerm(L);

    CALL_QDP(L);
    QDP_D_eq_r_times_D(c->ptr, &b, a->ptr, *qCurrent);

    return 1;
}

static int
q_r_mul_D(lua_State *L)
{
    QLA_Real a = luaL_checknumber(L, 1);
    mLatDirFerm *b = qlua_checkLatDirFerm(L, 2);
    mLatDirFerm *c = qlua_newLatDirFerm(L);

    CALL_QDP(L);
    QDP_D_eq_r_times_D(c->ptr, &a, b->ptr, *qCurrent);

    return 1;
}

static int
q_D_mul_c(lua_State *L)
{
    mLatDirFerm *a = qlua_checkLatDirFerm(L, 1);
    QLA_Complex *b = qlua_checkComplex(L, 2);
    mLatDirFerm *c = qlua_newLatDirFerm(L);

    CALL_QDP(L);
    QDP_D_eq_c_times_D(c->ptr, b, a->ptr, *qCurrent);

    return 1;
}

static int
q_c_mul_D(lua_State *L)
{
    QLA_Complex *a = qlua_checkComplex(L, 1);
    mLatDirFerm *b = qlua_checkLatDirFerm(L, 2);
    mLatDirFerm *c = qlua_newLatDirFerm(L);

    CALL_QDP(L);
    QDP_D_eq_c_times_D(c->ptr, a, b->ptr, *qCurrent);

    return 1;
}

static struct {
    QLA_Real *a;
    QLA_DiracFermion *b;
} RDmul_args; /* YYY global state */

static void
do_RDmul(QLA_DiracFermion *r, int idx)
{
    QLA_D_eq_r_times_D(r, &RDmul_args.a[idx], &RDmul_args.b[idx]);
}

static void
X_D_eq_R_times_D(QDP_DiracFermion *r, QDP_Real *a, QDP_DiracFermion *b,
                 QDP_Subset s)
{
    RDmul_args.a = QDP_expose_R(a);
    RDmul_args.b = QDP_expose_D(b);
    QDP_D_eq_funci(r, do_RDmul, s);
    QDP_reset_R(a);
    QDP_reset_D(b);
    RDmul_args.a = 0;
    RDmul_args.b = 0;
}

static int
q_R_mul_D(lua_State *L)
{
    mLatReal *a = qlua_checkLatReal(L, 1);
    mLatDirFerm *b = qlua_checkLatDirFerm(L, 2);
    mLatDirFerm *c = qlua_newLatDirFerm(L);

    CALL_QDP(L);
    X_D_eq_R_times_D(c->ptr, a->ptr, b->ptr, *qCurrent);

    return 1;
}

static int
q_D_mul_R(lua_State *L)
{
    mLatDirFerm *a = qlua_checkLatDirFerm(L, 1);
    mLatReal *b = qlua_checkLatReal(L, 2);
    mLatDirFerm *c = qlua_newLatDirFerm(L);

    CALL_QDP(L);
    X_D_eq_R_times_D(c->ptr, b->ptr, a->ptr, *qCurrent);

    return 1;
}

static struct {
    QLA_Complex *a;
    QLA_DiracFermion *b;
} CDmul_args; /* YYY global state */

static void
do_CDmul(QLA_DiracFermion *r, int idx)
{
    QLA_D_eq_c_times_D(r, &CDmul_args.a[idx], &CDmul_args.b[idx]);
}

static void
X_D_eq_C_times_D(QDP_DiracFermion *r, QDP_Complex *a, QDP_DiracFermion *b,
                 QDP_Subset s)
{
    CDmul_args.a = QDP_expose_C(a);
    CDmul_args.b = QDP_expose_D(b);
    QDP_D_eq_funci(r, do_CDmul, s);
    QDP_reset_C(a);
    QDP_reset_D(b);
    CDmul_args.a = 0;
    CDmul_args.b = 0;
}

static int
q_C_mul_D(lua_State *L)
{
    mLatComplex *a = qlua_checkLatComplex(L, 1);
    mLatDirFerm *b = qlua_checkLatDirFerm(L, 2);
    mLatDirFerm *c = qlua_newLatDirFerm(L);

    CALL_QDP(L);
    X_D_eq_C_times_D(c->ptr, a->ptr, b->ptr, *qCurrent);

    return 1;
}

static int
q_D_mul_C(lua_State *L)
{
    mLatDirFerm *a = qlua_checkLatDirFerm(L, 1);
    mLatComplex *b = qlua_checkLatComplex(L, 2);
    mLatDirFerm *c = qlua_newLatDirFerm(L);

    CALL_QDP(L);
    X_D_eq_C_times_D(c->ptr, b->ptr, a->ptr, *qCurrent);

    return 1;
}

static int
q_M_mul_D(lua_State *L)
{
    mLatColMat *a = qlua_checkLatColMat(L, 1);
    mLatDirFerm *b = qlua_checkLatDirFerm(L, 2);
    mLatDirFerm *c = qlua_newLatDirFerm(L);

    CALL_QDP(L);
    QDP_D_eq_M_times_D(c->ptr, a->ptr, b->ptr, *qCurrent);

    return 1;
}

static int
q_D_div_r(lua_State *L)
{
    mLatDirFerm *a = qlua_checkLatDirFerm(L, 1);
    QLA_Real b = 1 / luaL_checknumber(L, 2);
    mLatDirFerm *c = qlua_newLatDirFerm(L);

    CALL_QDP(L);
    QDP_D_eq_r_times_D(c->ptr, &b, a->ptr, *qCurrent);

    return 1;
}

static int
q_D_div_c(lua_State *L)
{
    mLatDirFerm *a = qlua_checkLatDirFerm(L, 1);
    QLA_Complex *b = qlua_checkComplex(L, 2);
    mLatDirFerm *c = qlua_newLatDirFerm(L);
    double n = 1 / (QLA_real(*b) * QLA_real(*b) + QLA_imag(*b) * QLA_imag(*b));
    QLA_Complex s;

    CALL_QDP(L);
    QLA_real(s) = n * QLA_real(*b);
    QLA_imag(s) = -n * QLA_imag(*b);
    QDP_D_eq_c_times_D(c->ptr, &s, a->ptr, *qCurrent);

    return 1;
}

static int
q_latdirferm(lua_State *L)
{
    switch (lua_gettop(L)) {
    case 1: {
        mLatDirFerm *v = qlua_newLatDirFerm(L);

        CALL_QDP(L);
        QDP_D_eq_zero(v->ptr, *qCurrent);

        return 1;
    }
    case 2: {
        mLatDirFerm *f = qlua_checkLatDirFerm(L, 2);
        mLatDirFerm *v = qlua_newLatDirFerm(L);

        CALL_QDP(L);
        QDP_D_eq_D(v->ptr, f->ptr, *qCurrent);
        
        return 1;
    }
    case 3: {
        switch (qlua_gettype(L, 2)) {
        case qLatComplex: {
            mLatComplex *z = qlua_checkLatComplex(L, 2);
            int c = qlua_checkcolorindex(L, 3);
            int d = qlua_checkdiracindex(L, 3);
            mLatDirFerm *v = qlua_newLatDirFerm(L);

            CALL_QDP(L);
            QDP_D_eq_zero(v->ptr, *qCurrent);
            QDP_D_eq_elem_C(v->ptr, z->ptr, c, d, *qCurrent);

            return 1;
        }
        case qLatColVec: {
            mLatColVec *w = qlua_checkLatColVec(L, 2);
            int d = qlua_checkdiracindex(L, 3);
            mLatDirFerm *v = qlua_newLatDirFerm(L);

            CALL_QDP(L);
            QDP_D_eq_zero(v->ptr, *qCurrent);
            QDP_D_eq_colorvec_V(v->ptr, w->ptr, d, *qCurrent);

            return 1;
        }
        }
        break;
    }
    }
    return qlua_badconstr(L, "DiracFermion");
}

static struct luaL_Reg LatDirFermMethods[] = {
    { "norm2",      q_D_norm2 },
    { "shift",      q_D_shift },
    { "conj",       q_D_conj },
    { "gamma",      q_D_gamma },
    { "set",        q_D_set },
    { NULL,         NULL }
};

static struct luaL_Reg mtLatDirFerm[] = {
    { "__tostring",        q_D_fmt },
    { "__gc",              q_D_gc },
    { "__index",           q_D_get },
    { "__newindex",        q_D_put },
    { "__unm",             q_D_neg },
    { "__add",             qlua_add },
    { "__sub",             qlua_sub },
    { "__mul",             qlua_mul },
    { "__div",             qlua_div },
    { NULL,                NULL }
};

static struct luaL_Reg fLatDirFerm[] = {
    { "DiracFermion",      q_latdirferm },
    { NULL,                NULL }
};

int
init_latdirferm(lua_State *L)
{
    luaL_getmetatable(L, opLattice);
    luaL_register(L, NULL, fLatDirFerm);
    lua_pop(L, 1);
    qlua_metatable(L, mtnLatDirFerm, mtLatDirFerm);
    qlua_metatable(L, opLatDirFerm, LatDirFermMethods);
    qlua_reg_add(qLatDirFerm,  qLatDirFerm,  q_D_add_D);
    qlua_reg_sub(qLatDirFerm,  qLatDirFerm,  q_D_sub_D);
    qlua_reg_mul(qReal,        qLatDirFerm,  q_r_mul_D);
    qlua_reg_mul(qLatDirFerm,  qReal,        q_D_mul_r);
    qlua_reg_mul(qComplex,     qLatDirFerm,  q_c_mul_D);
    qlua_reg_mul(qLatDirFerm,  qComplex,     q_D_mul_c);
    qlua_reg_mul(qLatReal,     qLatDirFerm,  q_R_mul_D);
    qlua_reg_mul(qLatDirFerm,  qLatReal,     q_D_mul_R);
    qlua_reg_mul(qLatComplex,  qLatDirFerm,  q_C_mul_D);
    qlua_reg_mul(qLatDirFerm,  qLatComplex,  q_D_mul_C);
    qlua_reg_mul(qLatColMat,   qLatDirFerm,  q_M_mul_D);
    qlua_reg_div(qLatDirFerm,  qReal,        q_D_div_r);
    qlua_reg_div(qLatDirFerm,  qComplex,     q_D_div_c);
    qlua_reg_dot(qLatDirFerm,  q_D_dot);

    return 0;
}

int
fini_latdirferm(lua_State *L)
{
    return 0;
}

