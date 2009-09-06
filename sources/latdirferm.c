#include <qlua.h>                                                    /* DEPS */
#include <qcomplex.h>                                                /* DEPS */
#include <lattice.h>                                                 /* DEPS */
#include <latint.h>                                                  /* DEPS */
#include <latrandom.h>                                               /* DEPS */
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

                QDP_V_eq_colorvec_D(r->ptr, V->ptr, d, *qCurrent);
            } else {
                mLatComplex *r = qlua_newLatComplex(L);

                QDP_C_eq_elem_D(r->ptr, V->ptr, c, d, *qCurrent);
            }
        } else {
            if (c == -1) {
                qlua_free(L, idx);
                return qlua_badindex(L, "DiracFermion");
            } else {
                QLA_Complex *W = qlua_newComplex(L);
                QLA_DiracFermion *locked;
                double z_re, z_im;

                locked = QDP_expose_D(V->ptr);
                if (QDP_node_number(idx) == QDP_this_node) {
                    QLA_Complex *zz = &QLA_elem_D(locked[QDP_index(idx)], c, d);
                    z_re = QLA_real(*zz);
                    z_im = QLA_imag(*zz);
                } else {
                    z_re = 0;
                    z_im = 0;
                }
                QDP_reset_D(V->ptr);
                QMP_sum_double(&z_re);
                QMP_sum_double(&z_im);
                QLA_real(*W) = z_re;
                QLA_imag(*W) = z_im;
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

            QDP_D_eq_colorvec_V(V->ptr, z->ptr, d, *qCurrent);
        } else {
            mLatComplex *z = qlua_checkLatComplex(L, 3);

            QDP_D_eq_elem_C(V->ptr, z->ptr, c, d, *qCurrent);
        }
    } else {
        if (c == -1) {
            qlua_free(L, idx);
            return qlua_badindex(L, "DiracFermion");
        } else {
            QLA_Complex *z = qlua_checkComplex(L, 3);
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

    QDP_D_eq_r_times_D(r->ptr, &m1, a->ptr, *qCurrent);

    return 1;
}

static int
q_D_norm2(lua_State *L)
{
    mLatDirFerm *a = qlua_checkLatDirFerm(L, 1);
    QLA_Real n;

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

    QDP_D_eq_sD(r->ptr, a->ptr, shift, dir, *qCurrent);

    return 1;
}

static int
q_D_conj(lua_State *L)
{
    mLatDirFerm *a = qlua_checkLatDirFerm(L, 1);
    mLatDirFerm *r = qlua_newLatDirFerm(L);

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

        if (mu < 5) {
            QDP_D_eq_gamma_times_D(r->ptr, f->ptr, 1 << mu, *qCurrent);
        } else {
            QDP_D_eq_gamma_times_D(r->ptr, f->ptr, 15, *qCurrent);
        }

        return 1;
    }
    if (mu == -1) {
        mLatDirFerm *r = qlua_newLatDirFerm(L);

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

    QDP_C_eq_D_dot_D(s->ptr, a->ptr, b->ptr, *qCurrent);

    return 1;
}

int
q_D_gaussian(lua_State *L)
{
    mLatRandom *a = qlua_checkLatRandom(L, 1);
    mLatDirFerm *r = qlua_newLatDirFerm(L);

    QDP_D_eq_gaussian_S(r->ptr, a->ptr, *qCurrent);

    return 1;
}


static int
q_D_add_D(lua_State *L)
{
    mLatDirFerm *a = qlua_checkLatDirFerm(L, 1);
    mLatDirFerm *b = qlua_checkLatDirFerm(L, 2);
    mLatDirFerm *c = qlua_newLatDirFerm(L);

    QDP_D_eq_D_plus_D(c->ptr, a->ptr, b->ptr, *qCurrent);

    return 1;
}

static int
q_D_sub_D(lua_State *L)
{
    mLatDirFerm *a = qlua_checkLatDirFerm(L, 1);
    mLatDirFerm *b = qlua_checkLatDirFerm(L, 2);
    mLatDirFerm *c = qlua_newLatDirFerm(L);

    QDP_D_eq_D_minus_D(c->ptr, a->ptr, b->ptr, *qCurrent);

    return 1;
}

static int
q_D_mul_r(lua_State *L)
{
    mLatDirFerm *a = qlua_checkLatDirFerm(L, 1);
    QLA_Real b = luaL_checknumber(L, 2);
    mLatDirFerm *c = qlua_newLatDirFerm(L);

    QDP_D_eq_r_times_D(c->ptr, &b, a->ptr, *qCurrent);

    return 1;
}

static int
q_r_mul_D(lua_State *L)
{
    QLA_Real a = luaL_checknumber(L, 1);
    mLatDirFerm *b = qlua_checkLatDirFerm(L, 2);
    mLatDirFerm *c = qlua_newLatDirFerm(L);

    QDP_D_eq_r_times_D(c->ptr, &a, b->ptr, *qCurrent);

    return 1;
}

static int
q_D_mul_c(lua_State *L)
{
    mLatDirFerm *a = qlua_checkLatDirFerm(L, 1);
    QLA_Complex *b = qlua_checkComplex(L, 2);
    mLatDirFerm *c = qlua_newLatDirFerm(L);

    QDP_D_eq_c_times_D(c->ptr, b, a->ptr, *qCurrent);

    return 1;
}

static int
q_c_mul_D(lua_State *L)
{
    QLA_Complex *a = qlua_checkComplex(L, 1);
    mLatDirFerm *b = qlua_checkLatDirFerm(L, 2);
    mLatDirFerm *c = qlua_newLatDirFerm(L);

    QDP_D_eq_c_times_D(c->ptr, a, b->ptr, *qCurrent);

    return 1;
}

static int
q_M_mul_D(lua_State *L)
{
    mLatColMat *a = qlua_checkLatColMat(L, 1);
    mLatDirFerm *b = qlua_checkLatDirFerm(L, 2);
    mLatDirFerm *c = qlua_newLatDirFerm(L);

    QDP_D_eq_M_times_D(c->ptr, a->ptr, b->ptr, *qCurrent);

    return 1;
}

static int
q_D_div_r(lua_State *L)
{
    mLatDirFerm *a = qlua_checkLatDirFerm(L, 1);
    QLA_Real b = 1 / luaL_checknumber(L, 2);
    mLatDirFerm *c = qlua_newLatDirFerm(L);

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

        QDP_D_eq_zero(v->ptr, *qCurrent);

        return 1;
    }
    case 2: {
        mLatDirFerm *f = qlua_checkLatDirFerm(L, 2);
        mLatDirFerm *v = qlua_newLatDirFerm(L);

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

            QDP_D_eq_zero(v->ptr, *qCurrent);
            QDP_D_eq_elem_C(v->ptr, z->ptr, c, d, *qCurrent);

            return 1;
        }
        case qLatColVec: {
            mLatColVec *w = qlua_checkLatColVec(L, 2);
            int d = qlua_checkdiracindex(L, 3);
            mLatDirFerm *v = qlua_newLatDirFerm(L);

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
    qlua_reg_add(qLatDirFerm, qLatDirFerm, q_D_add_D);
    qlua_reg_sub(qLatDirFerm, qLatDirFerm, q_D_sub_D);
    qlua_reg_mul(qReal,       qLatDirFerm, q_r_mul_D);
    qlua_reg_mul(qLatDirFerm, qReal,       q_D_mul_r);
    qlua_reg_mul(qComplex,    qLatDirFerm, q_c_mul_D);
    qlua_reg_mul(qLatDirFerm, qComplex,    q_D_mul_c);
    qlua_reg_mul(qLatColMat,  qLatDirFerm, q_M_mul_D);
    qlua_reg_div(qLatDirFerm, qReal,       q_D_div_r);
    qlua_reg_div(qLatDirFerm, qComplex,    q_D_div_c);
    qlua_reg_dot(qLatDirFerm, q_D_dot);

    return 0;
}

int
fini_latdirferm(lua_State *L)
{
    return 0;
}

