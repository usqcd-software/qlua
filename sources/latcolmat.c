#include <qlua.h>                                                    /* DEPS */
#include <qcomplex.h>                                                /* DEPS */
#include <lattice.h>                                                 /* DEPS */
#include <latint.h>                                                  /* DEPS */
#include <latreal.h>                                                 /* DEPS */
#include <latcomplex.h>                                              /* DEPS */
#include <latcolmat.h>                                               /* DEPS */
#include <latrandom.h>                                               /* DEPS */
#include <latcolvec.h>                                               /* DEPS */
#include <qmp.h>
#include <math.h>

const char mtnLatColMat[] = "lattice.ColorMatrix";
static char opLatColMat[] = "lattice.ColorMatrix.op";

mLatColMat *
qlua_newLatColMat(lua_State *L)
{
    QDP_ColorMatrix *v = QDP_create_M();
    mLatColMat *hdr;

    if (v == 0) {
        lua_gc(L, LUA_GCCOLLECT, 0);
        v = QDP_create_M();
        if (v == 0)
            luaL_error(L, "not enough memory (QDP_ColorMatrix)");
    }
    hdr = lua_newuserdata(L, sizeof (mLatColMat));
    hdr->ptr = v;
    luaL_getmetatable(L, mtnLatColMat);
    lua_setmetatable(L, -2);

    return hdr;
}

mLatColMat *
qlua_checkLatColMat(lua_State *L, int idx)
{
    void *v = luaL_checkudata(L, idx, mtnLatColMat);
    
    luaL_argcheck(L, v != 0, idx, "lattice.ColorMatrix expected");
    
    return v;
}

static int
q_M_fmt(lua_State *L)
{
    char fmt[72];
    mLatColMat *b = qlua_checkLatColMat(L, 1);

    sprintf(fmt, "QDP:ColorMatrix(%p)", b->ptr);
    lua_pushstring(L, fmt);

    return 1;
}

static int
q_M_gc(lua_State *L)
{
    mLatColMat *b = qlua_checkLatColMat(L, 1);

    QDP_destroy_M(b->ptr);
    b->ptr = 0;

    return 0;
}

static int
q_M_get(lua_State *L)
{
    switch (qlua_gettype(L, 2)) {
    case qTable: {
        mLatColMat *V = qlua_checkLatColMat(L, 1);
        int b = qlua_checkrightindex(L, 2);
        int a = qlua_leftindex(L, 2);
        int *idx = qlua_latcoord(L, 2);
        if (idx == NULL) {
            if (a == -1) {
                mLatColVec *r = qlua_newLatColVec(L);
                
                QDP_V_eq_colorvec_M(r->ptr, V->ptr, b, *qCurrent);
            } else {
                mLatComplex *r = qlua_newLatComplex(L);
                
                QDP_C_eq_elem_M(r->ptr, V->ptr, a, b, *qCurrent);
            }
        } else {
            if (a == -1) {
                qlua_free(L, idx);
                return qlua_badindex(L, "ColorMatrix");
            } else {
                QLA_Complex *W = qlua_newComplex(L);
                QLA_ColorMatrix *locked;
                double z_re, z_im;
                
                locked = QDP_expose_M(V->ptr);
                if (QDP_node_number(idx) == QDP_this_node) {
                    QLA_Complex *zz = &QLA_elem_M(locked[QDP_index(idx)], a, b);
                    z_re = QLA_real(*zz);
                    z_im = QLA_imag(*zz);
                } else {
                    z_re = 0;
                    z_im = 0;
                }
                QDP_reset_M(V->ptr);
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
        return qlua_lookup(L, 2, opLatColMat);
    }
    return qlua_badindex(L, "ColorMatrix");
}

static int
q_M_put(lua_State *L)
{
    mLatColMat *V = qlua_checkLatColMat(L, 1);
    int b = qlua_checkrightindex(L, 2);
    int a = qlua_leftindex(L, 2);
    int *idx = qlua_latcoord(L, 2);

    if (idx == NULL) {
        if (a == -1) {
            mLatColVec *z = qlua_checkLatColVec(L, 3);

            QDP_M_eq_colorvec_V(V->ptr, z->ptr, b, *qCurrent);
        } else {
            mLatComplex *z = qlua_checkLatComplex(L, 3);
            
            QDP_M_eq_elem_C(V->ptr, z->ptr, a, b, *qCurrent);
        }
    } else {
        if (a == -1) {
            qlua_free(L, idx);
            return qlua_badindex(L, "ColorMatrix");
        } else {
            QLA_Complex *z = qlua_checkComplex(L, 3);
            if (QDP_node_number(idx) == QDP_this_node) {
                QLA_ColorMatrix *locked = QDP_expose_M(V->ptr);
                QLA_Complex *zz = &QLA_elem_M(locked[QDP_index(idx)], a, b);
                QLA_c_eq_c(*zz, *z);
                QDP_reset_M(V->ptr);
            }
        }
    }
    qlua_free(L, idx);
    return 0;
}

static int
q_M_norm2(lua_State *L)
{
    mLatColMat *a = qlua_checkLatColMat(L, 1);
    QLA_Real n;

    QDP_r_eq_norm2_M(&n, a->ptr, *qCurrent);
    lua_pushnumber(L, n);
    
    return 1;
}

static int
q_M_shift(lua_State *L)
{
    mLatColMat *a = qlua_checkLatColMat(L, 1);
    QDP_Shift shift = qlua_checkShift(L, 2);
    QDP_ShiftDir dir = qlua_checkShiftDir(L, 3);
    mLatColMat *r = qlua_newLatColMat(L);

    QDP_M_eq_sM(r->ptr, a->ptr, shift, dir, *qCurrent);

    return 1;
}

static int
q_M_conj(lua_State *L)
{
    mLatColMat *a = qlua_checkLatColMat(L, 1);
    mLatColMat *r = qlua_newLatColMat(L);

    QDP_M_eq_conj_M(r->ptr, a->ptr, *qCurrent);

    return 1;
}

static int
q_M_trans(lua_State *L)
{
    mLatColMat *a = qlua_checkLatColMat(L, 1);
    mLatColMat *r = qlua_newLatColMat(L);

    QDP_M_eq_transpose_M(r->ptr, a->ptr, *qCurrent);

    return 1;
}

static int
q_M_adjoin(lua_State *L)
{
    mLatColMat *a = qlua_checkLatColMat(L, 1);
    mLatColMat *r = qlua_newLatColMat(L);

    QDP_M_eq_Ma(r->ptr, a->ptr, *qCurrent);

    return 1;
}

static int
q_M_trace(lua_State *L)
{
    mLatColMat *a = qlua_checkLatColMat(L, 1);
    mLatComplex *r = qlua_newLatComplex(L);

    QDP_C_eq_trace_M(r->ptr, a->ptr, *qCurrent);

    return 1;
}

struct {
    QLA_ColorMatrix *a;
} Mexp_args; /* YYY global state */

static void
do_Mexp(QLA_ColorMatrix *r, int idx)
{
    static const int n = 15; /* terms in the Taylor series */
    QLA_Complex cone;
    QLA_ColorMatrix mone;
    QLA_ColorMatrix axs;
    QLA_ColorMatrix v;
    QLA_ColorMatrix *a = &Mexp_args.a[idx];
    int i, j, k;
    double max_el;
    double s;
    int pw;

    for (max_el = 0, i = 0; i < QDP_Nc; i++) {
        for (j = 0; j < QDP_Nc; j++) {
            QLA_Complex z;
            double v;
            QLA_C_eq_elem_M(&z, a, i, j);
            v = fabs(QLA_real(z));
            if (max_el < v)
                max_el = v;
            v = fabs(QLA_imag(z));
            if (max_el < v)
                max_el = v;
        }
    }
    pw = (int)(ceil(max_el * 2 * QDP_Nc));
    QLA_c_eq_r_plus_ir(cone, 1.0, 0.0);
    QLA_M_eq_c(&mone, &cone);

    s = 1.0 / (n * pw);
    QLA_M_eq_r_times_M_plus_M(&v, &s, a, &mone);
    for (k = n; --k;) {
        QLA_M_eq_M_times_M(&axs, a, &v);
        s = 1.0 / (k * pw);
        QLA_M_eq_r_times_M_plus_M(&v, &s, &axs, &mone);
    }

    QLA_M_eq_M(r, &mone);
    while (pw) {
        if (pw & 1) {
            QLA_M_eq_M_times_M(&axs, r, &v);
            QLA_M_eq_M(r, &axs);
        }
        pw >>= 1;
        if (pw > 0) {
            QLA_M_eq_M_times_M(&axs, &v, &v);
            QLA_M_eq_M(&v, &axs);
        }
    }
}

static void
X_M_eq_exp_M(QDP_ColorMatrix *r, QDP_ColorMatrix *a, QDP_Subset s)
{
    Mexp_args.a = QDP_expose_M(a);
    QDP_M_eq_funci(r, do_Mexp, s);
    QDP_reset_M(a);
    Mexp_args.a = 0;
}

static int
q_M_exp(lua_State *L)
{
    mLatColMat *a = qlua_checkLatColMat(L, 1);
    mLatColMat *r = qlua_newLatColMat(L);

    X_M_eq_exp_M(r->ptr, a->ptr, *qCurrent);

    return 1;
}

static int
q_M_set(lua_State *L)
{
    mLatColMat *r = qlua_checkLatColMat(L, 1);
    mLatColMat *a = qlua_checkLatColMat(L, 2);

    QDP_M_eq_M(r->ptr, a->ptr, *qCurrent);
    lua_pop(L, 1);

    return 1;
}

static int
q_M_dot(lua_State *L)
{
    mLatColMat *a = qlua_checkLatColMat(L, 1);
    mLatColMat *b = qlua_checkLatColMat(L, 2);
    mLatComplex *s = qlua_newLatComplex(L);

    QDP_C_eq_M_dot_M(s->ptr, a->ptr, b->ptr, *qCurrent);

    return 1;
}

int
q_M_gaussian(lua_State *L)
{
    mLatRandom *a = qlua_checkLatRandom(L, 1);
    mLatColMat *r = qlua_newLatColMat(L);

    QDP_M_eq_gaussian_S(r->ptr, a->ptr, *qCurrent);

    return 1;
}


static int
q_M_add_M(lua_State *L)
{
    mLatColMat *a = qlua_checkLatColMat(L, 1);
    mLatColMat *b = qlua_checkLatColMat(L, 2);
    mLatColMat *c = qlua_newLatColMat(L);

    QDP_M_eq_M_plus_M(c->ptr, a->ptr, b->ptr, *qCurrent);

    return 1;
}

static int
q_M_sub_M(lua_State *L)
{
    mLatColMat *a = qlua_checkLatColMat(L, 1);
    mLatColMat *b = qlua_checkLatColMat(L, 2);
    mLatColMat *c = qlua_newLatColMat(L);

    QDP_M_eq_M_minus_M(c->ptr, a->ptr, b->ptr, *qCurrent);

    return 1;
}

static int
q_M_neg(lua_State *L)
{
    mLatColMat *a = qlua_checkLatColMat(L, 1);
    mLatColMat *r = qlua_newLatColMat(L);
    QLA_Real m1 = -1;

    QDP_M_eq_r_times_M(r->ptr, &m1, a->ptr, *qCurrent);

    return 1;
}

static int
q_M_mul_M(lua_State *L)
{
    mLatColMat *a = qlua_checkLatColMat(L, 1);
    mLatColMat *b = qlua_checkLatColMat(L, 2);
    mLatColMat *c = qlua_newLatColMat(L);

    QDP_M_eq_M_times_M(c->ptr, a->ptr, b->ptr, *qCurrent);

    return 1;
}

static int
q_M_mul_V(lua_State *L)
{
    mLatColMat *a = qlua_checkLatColMat(L, 1);
    mLatColVec *b = qlua_checkLatColVec(L, 2);
    mLatColVec *c = qlua_newLatColVec(L);

    QDP_V_eq_M_times_V(c->ptr, a->ptr, b->ptr, *qCurrent);

    return 1;
}

static int
q_M_mul_r(lua_State *L)
{
    mLatColMat *a = qlua_checkLatColMat(L, 1);
    QLA_Real b = luaL_checknumber(L, 2);
    mLatColMat *c = qlua_newLatColMat(L);

    QDP_M_eq_r_times_M(c->ptr, &b, a->ptr, *qCurrent);

    return 1;
}

static int
q_r_mul_M(lua_State *L)
{
    QLA_Real a = luaL_checknumber(L, 1);
    mLatColMat *b = qlua_checkLatColMat(L, 2);
    mLatColMat *c = qlua_newLatColMat(L);

    QDP_M_eq_r_times_M(c->ptr, &a, b->ptr, *qCurrent);

    return 1;
}

static int
q_M_mul_c(lua_State *L)
{
    mLatColMat *a = qlua_checkLatColMat(L, 1);
    QLA_Complex *b = qlua_checkComplex(L, 2);
    mLatColMat *c = qlua_newLatColMat(L);

    QDP_M_eq_c_times_M(c->ptr, b, a->ptr, *qCurrent);

    return 1;
}

static int
q_c_mul_M(lua_State *L)
{
    QLA_Complex *a = qlua_checkComplex(L, 1);
    mLatColMat *b = qlua_checkLatColMat(L, 2);
    mLatColMat *c = qlua_newLatColMat(L);

    QDP_M_eq_c_times_M(c->ptr, a, b->ptr, *qCurrent);

    return 1;
}

static struct {
    QLA_Real *a;
    QLA_ColorMatrix *b;
} RMmul_args; /* YYY global state */

static void
do_RMmul(QLA_ColorMatrix *r, int idx)
{
    QLA_M_eq_r_times_M(r, &RMmul_args.a[idx], &RMmul_args.b[idx]);
}

static void
X_M_eq_R_times_M(QDP_ColorMatrix *r, QDP_Real *a, QDP_ColorMatrix *b,
                 QDP_Subset s)
{
    RMmul_args.a = QDP_expose_R(a);
    RMmul_args.b = QDP_expose_M(b);
    QDP_M_eq_funci(r, do_RMmul, s);
    QDP_reset_R(a);
    QDP_reset_M(b);
    RMmul_args.a = 0;
    RMmul_args.b = 0;
}

static int
q_R_mul_M(lua_State *L)
{
    mLatReal *a = qlua_checkLatReal(L, 1);
    mLatColMat *b = qlua_checkLatColMat(L, 2);
    mLatColMat *r = qlua_newLatColMat(L);

    X_M_eq_R_times_M(r->ptr, a->ptr, b->ptr, *qCurrent);

    return 1;
}

static int
q_M_mul_R(lua_State *L)
{
    mLatColMat *a = qlua_checkLatColMat(L, 1);
    mLatReal *b = qlua_checkLatReal(L, 2);
    mLatColMat *r = qlua_newLatColMat(L);

    X_M_eq_R_times_M(r->ptr, b->ptr, a->ptr, *qCurrent);

    return 1;
}

static struct {
    QLA_Complex *a;
    QLA_ColorMatrix *b;
} CMmul_args; /* YYY global state */

static void
do_CMmul(QLA_ColorMatrix *r, int idx)
{
    QLA_M_eq_c_times_M(r, &CMmul_args.a[idx], &CMmul_args.b[idx]);
}

static void
X_M_eq_C_times_M(QDP_ColorMatrix *r, QDP_Complex *a, QDP_ColorMatrix *b,
                 QDP_Subset s)
{
    CMmul_args.a = QDP_expose_C(a);
    CMmul_args.b = QDP_expose_M(b);
    QDP_M_eq_funci(r, do_CMmul, s);
    QDP_reset_C(a);
    QDP_reset_M(b);
    CMmul_args.a = 0;
    CMmul_args.b = 0;
}

static int
q_C_mul_M(lua_State *L)
{
    mLatComplex *a = qlua_checkLatComplex(L, 1);
    mLatColMat *b = qlua_checkLatColMat(L, 2);
    mLatColMat *r = qlua_newLatColMat(L);

    X_M_eq_C_times_M(r->ptr, a->ptr, b->ptr, *qCurrent);

    return 1;
}

static int
q_M_mul_C(lua_State *L)
{
    mLatColMat *a = qlua_checkLatColMat(L, 1);
    mLatComplex *b = qlua_checkLatComplex(L, 2);
    mLatColMat *r = qlua_newLatColMat(L);

    X_M_eq_C_times_M(r->ptr, b->ptr, a->ptr, *qCurrent);

    return 1;
}

static int
q_M_div_r(lua_State *L)
{
    mLatColMat *a = qlua_checkLatColMat(L, 1);
    QLA_Real b = 1 / luaL_checknumber(L, 2);
    mLatColMat *c = qlua_newLatColMat(L);

    QDP_M_eq_r_times_M(c->ptr, &b, a->ptr, *qCurrent);

    return 1;
}

static int
q_M_div_c(lua_State *L)
{
    mLatColMat *a = qlua_checkLatColMat(L, 1);
    QLA_Complex *b = qlua_checkComplex(L, 2);
    mLatColMat *c = qlua_newLatColMat(L);
    double n = 1 / (QLA_real(*b) * QLA_real(*b) + QLA_imag(*b) * QLA_imag(*b));
    QLA_Complex s;

    QLA_real(s) = n * QLA_real(*b);
    QLA_imag(s) = -n * QLA_imag(*b);
    QDP_M_eq_c_times_M(c->ptr, &s, a->ptr, *qCurrent);

    return 1;
}

static int
q_latcolmat(lua_State *L)
{
    switch (lua_gettop(L)) {
    case 1: {
        mLatColMat *v = qlua_newLatColMat(L);

        QDP_M_eq_zero(v->ptr, *qCurrent);
        
        return 1;
    }
    case 2: {
        switch (qlua_gettype(L, 2)) {
        case qReal: {
            QLA_Real x = luaL_checknumber(L, 2);
            mLatColMat *v = qlua_newLatColMat(L);
            QLA_Complex z;

            QLA_real(z) = x;
            QLA_imag(z) = 0;
            QDP_M_eq_c(v->ptr, &z, *qCurrent);

            return 1;
        }
        case qComplex: {
            QLA_Complex *z = qlua_checkComplex(L, 2);
            mLatColMat *v = qlua_newLatColMat(L);

            QDP_M_eq_c(v->ptr, z, *qCurrent);

            return 1;
        }
        case qLatColMat: {
            mLatColMat *w = qlua_checkLatColMat(L, 2);
            mLatColMat *v = qlua_newLatColMat(L);

            QDP_M_eq_M(v->ptr, w->ptr, *qCurrent);

            return 1;
        }
        }
        break;
    }
    case 3: {
        switch (qlua_gettype(L, 2)) {
        case qLatComplex: {
            mLatComplex *z = qlua_checkLatComplex(L, 2);
            int a = qlua_checkleftindex(L, 3);
            int b = qlua_checkrightindex(L, 3);
            mLatColMat *v = qlua_newLatColMat(L);

            QDP_M_eq_zero(v->ptr, *qCurrent);
            QDP_M_eq_elem_C(v->ptr, z->ptr, a, b, *qCurrent);

            return 1;
        }
        case qLatColVec: {
            mLatColVec *f = qlua_checkLatColVec(L, 2);
            switch (qlua_gettype(L, 3)) {
            case qTable: {
                int b = qlua_checkrightindex(L, 3);
                mLatColMat *v = qlua_newLatColMat(L);
                
                QDP_M_eq_zero(v->ptr, *qCurrent);
                QDP_M_eq_colorvec_V(v->ptr, f->ptr, b, *qCurrent);
                
                return 1;
            }
            case qLatColVec: {
                mLatColVec *g = qlua_checkLatColVec(L, 3);
                mLatColMat *v = qlua_newLatColMat(L);

                QDP_M_eq_V_times_Va(v->ptr, f->ptr, g->ptr, *qCurrent);

                return 1;
            }
            }
            break;
        }
        }
        break;
    }
    }
    return qlua_badconstr(L, "ColorMatrix");
}

static struct luaL_Reg LatColMatMethods[] = {
    { "norm2",      q_M_norm2 },
    { "shift",      q_M_shift },
    { "conj",       q_M_conj },
    { "transpose",  q_M_trans },
    { "adjoin",     q_M_adjoin },
    { "trace",      q_M_trace },
    { "set",        q_M_set },
    { "exp",        q_M_exp },
    { NULL,         NULL }
};

static struct luaL_Reg mtLatColMat[] = {
    { "__tostring",        q_M_fmt },
    { "__gc",              q_M_gc },
    { "__index",           q_M_get },
    { "__newindex",        q_M_put },
    { "__unm",             q_M_neg },
    { "__add",             qlua_add },
    { "__sub",             qlua_sub },
    { "__mul",             qlua_mul },
    { "__div",             qlua_div },
    { NULL,                NULL }
};

static struct luaL_Reg fLatColMat[] = {
    { "ColorMatrix",       q_latcolmat },
    { NULL,                NULL }
};

int
init_latcolmat(lua_State *L)
{
    luaL_getmetatable(L, opLattice);
    luaL_register(L, NULL, fLatColMat);
    lua_pop(L, 1);
    qlua_metatable(L, mtnLatColMat, mtLatColMat);
    qlua_metatable(L, opLatColMat, LatColMatMethods);
    qlua_reg_add(qLatColMat,  qLatColMat,  q_M_add_M);
    qlua_reg_sub(qLatColMat,  qLatColMat,  q_M_sub_M);
    qlua_reg_mul(qReal,       qLatColMat,  q_r_mul_M);
    qlua_reg_mul(qLatColMat,  qReal,       q_M_mul_r);
    qlua_reg_mul(qComplex,    qLatColMat,  q_c_mul_M);
    qlua_reg_mul(qLatColMat,  qComplex,    q_M_mul_c);
    qlua_reg_mul(qLatReal,    qLatColMat,  q_R_mul_M);
    qlua_reg_mul(qLatColMat,  qLatReal,    q_M_mul_R);
    qlua_reg_mul(qLatComplex, qLatColMat,  q_C_mul_M);
    qlua_reg_mul(qLatColMat,  qLatComplex, q_M_mul_C);
    qlua_reg_mul(qLatColMat,  qLatColVec,  q_M_mul_V);
    qlua_reg_mul(qLatColMat,  qLatColMat,  q_M_mul_M);
    qlua_reg_div(qLatColMat,  qReal,       q_M_div_r);
    qlua_reg_div(qLatColMat,  qComplex,    q_M_div_c);
    qlua_reg_dot(qLatColMat,  q_M_dot);

    return 0;
}

int
fini_latcolmat(lua_State *L)
{
    return 0;
}
