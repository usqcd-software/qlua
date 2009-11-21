#include <qlua.h>                                                    /* DEPS */
#include <qcomplex.h>                                                /* DEPS */
#include <qvector.h>                                                 /* DEPS */
#include <lattice.h>                                                 /* DEPS */
#include <latcomplex.h>                                              /* DEPS */
#include <latreal.h>                                                 /* DEPS */
#include <latint.h>                                                  /* DEPS */
#include <latrandom.h>                                               /* DEPS */
#include <latmulti.h>                                                /* DEPS */
#include <math.h>
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

    CALL_QDP(L);
    QDP_C_eq_r_times_C(r->ptr, &m1, a->ptr, *qCurrent);
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
        CALL_QDP(L);
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
        z_im = 0;
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
    CALL_QDP(L);
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

    CALL_QDP(L);
    QDP_R_eq_re_C(c->ptr, a->ptr, *qCurrent);

    return 1;
}

static int
q_C_imag(lua_State *L)
{
    mLatComplex *a = qlua_checkLatComplex(L, 1);
    mLatReal *c = qlua_newLatReal(L);

    CALL_QDP(L);
    QDP_R_eq_im_C(c->ptr, a->ptr, *qCurrent);

    return 1;
}

static int
q_C_add_C(lua_State *L)
{
    mLatComplex *a = qlua_checkLatComplex(L, 1);
    mLatComplex *b = qlua_checkLatComplex(L, 2);
    mLatComplex *c = qlua_newLatComplex(L);

    CALL_QDP(L);
    QDP_C_eq_C_plus_C(c->ptr, a->ptr, b->ptr, *qCurrent);

    return 1;
}

static struct {
    QLA_Real *r;
    QLA_Complex *c;
} CRarg; /* YYY global state */

static void
X_C_eq_R_op_C(QDP_Complex *r, QDP_Real *a, QDP_Complex *b,
              void (*op)(QLA_Complex *a, int idx),
              QDP_Subset s)
{
    CRarg.r = QDP_expose_R(a);
    CRarg.c = QDP_expose_C(b);
    QDP_C_eq_funci(r, op, s);
    CRarg.r = 0;
    CRarg.c = 0;
    QDP_reset_R(a);
    QDP_reset_C(b);
}

static void
doRCadd(QLA_Complex *r, int idx)
{
    QLA_Complex *v = &CRarg.c[idx];

    QLA_real(*r) = QLA_real(*v) + CRarg.r[idx];
    QLA_imag(*r) = QLA_imag(*v);
}

static int
q_C_add_R(lua_State *L)
{
    mLatComplex *a = qlua_checkLatComplex(L, 1);
    mLatReal *b = qlua_checkLatReal(L, 2);
    mLatComplex *r = qlua_newLatComplex(L);

    CALL_QDP(L);
    X_C_eq_R_op_C(r->ptr, b->ptr, a->ptr, doRCadd, *qCurrent);

    return 1;
}

static int
q_R_add_C(lua_State *L)
{
    mLatReal *a = qlua_checkLatReal(L, 1);
    mLatComplex *b = qlua_checkLatComplex(L, 2);
    mLatComplex *r = qlua_newLatComplex(L);

    CALL_QDP(L);
    X_C_eq_R_op_C(r->ptr, a->ptr, b->ptr, doRCadd, *qCurrent);

    return 1;
}

static void
X_C_eq_r_op_C(QDP_Complex *r, QLA_Real a, QDP_Complex *b,
              void (*op)(QLA_Complex *a, int idx),
              QDP_Subset s)
{
    CRarg.r = &a;
    CRarg.c = QDP_expose_C(b);
    QDP_C_eq_funci(r, op, s);
    CRarg.r = 0;
    CRarg.c = 0;
    QDP_reset_C(b);
}

static void
dorCadd(QLA_Complex *r, int idx)
{
    QLA_Complex *v = &CRarg.c[idx];

    QLA_real(*r) = QLA_real(*v) + *CRarg.r;
    QLA_imag(*r) = QLA_imag(*v);
}

static int
q_r_add_C(lua_State *L)
{
    QLA_Real a = luaL_checknumber(L, 1);
    mLatComplex *b = qlua_checkLatComplex(L, 2);
    mLatComplex *r = qlua_newLatComplex(L);

    CALL_QDP(L);
    X_C_eq_r_op_C(r->ptr, a, b->ptr, dorCadd, *qCurrent);

    return 1;
}

static int
q_C_add_r(lua_State *L)
{
    mLatComplex *a = qlua_checkLatComplex(L, 1);
    QLA_Real b = luaL_checknumber(L, 2);
    mLatComplex *r = qlua_newLatComplex(L);

    CALL_QDP(L);
    X_C_eq_r_op_C(r->ptr, b, a->ptr, dorCadd, *qCurrent);

    return 1;
}

static struct {
    QLA_Complex *a;
    QLA_Complex *b;
} cCarg; /* YYY global state */

static void
X_C_eq_c_op_C(QDP_Complex *r, QLA_Complex *a, QDP_Complex *b,
              void (*op)(QLA_Complex *a, int idx),
              QDP_Subset s)
{
    cCarg.a = a;
    cCarg.b = QDP_expose_C(b);
    QDP_C_eq_funci(r, op, s);
    cCarg.a = 0;
    cCarg.b = 0;
    QDP_reset_C(b);
}

static void
docCadd(QLA_Complex *r, int idx)
{
    QLA_C_eq_C_plus_C(r, cCarg.a, &cCarg.b[idx]);
}

static int
q_c_add_C(lua_State *L)
{
    QLA_Complex *a = qlua_checkComplex(L, 1);
    mLatComplex *b = qlua_checkLatComplex(L, 2);
    mLatComplex *r = qlua_newLatComplex(L);

    CALL_QDP(L);
    X_C_eq_c_op_C(r->ptr, a, b->ptr, docCadd, *qCurrent);

    return 1;
}

static int
q_C_add_c(lua_State *L)
{
    mLatComplex *a = qlua_checkLatComplex(L, 1);
    QLA_Complex *b = qlua_checkComplex(L, 2);
    mLatComplex *r = qlua_newLatComplex(L);

    CALL_QDP(L);
    X_C_eq_c_op_C(r->ptr, b, a->ptr, docCadd, *qCurrent);

    return 1;
}

static struct {
    QLA_Complex *a;
    QLA_Real *b;
} cRarg; /* YYY global state */

static void
X_C_eq_c_op_R(QDP_Complex *r, QLA_Complex *a, QDP_Real *b,
              void (*op)(QLA_Complex *a, int idx),
              QDP_Subset s)
{
    cRarg.a = a;
    cRarg.b = QDP_expose_R(b);
    QDP_C_eq_funci(r, op, s);
    cRarg.a = 0;
    cRarg.b = 0;
    QDP_reset_R(b);
}

static void
docRadd(QLA_Complex *r, int idx)
{
    QLA_real(*r) = QLA_real(*cRarg.a) + cRarg.b[idx];
    QLA_imag(*r) = QLA_imag(*cRarg.a);
}

static int
q_c_add_R(lua_State *L)
{
    QLA_Complex *a = qlua_checkComplex(L, 1);
    mLatReal *b = qlua_checkLatReal(L, 2);
    mLatComplex *r = qlua_newLatComplex(L);

    CALL_QDP(L);
    X_C_eq_c_op_R(r->ptr, a, b->ptr, docRadd, *qCurrent);

    return 1;
}

static int
q_R_add_c(lua_State *L)
{
    mLatReal *a = qlua_checkLatReal(L, 1);
    QLA_Complex *b = qlua_checkComplex(L, 2);
    mLatComplex *r = qlua_newLatComplex(L);

    CALL_QDP(L);
    X_C_eq_c_op_R(r->ptr, b, a->ptr, docRadd, *qCurrent);

    return 1;
}

static int
q_C_sub_C(lua_State *L)
{
    mLatComplex *a = qlua_checkLatComplex(L, 1);
    mLatComplex *b = qlua_checkLatComplex(L, 2);
    mLatComplex *c = qlua_newLatComplex(L);

    CALL_QDP(L);
    QDP_C_eq_C_minus_C(c->ptr, a->ptr, b->ptr, *qCurrent);

    return 1;
}

static void
doRCsub(QLA_Complex *r, int idx)
{
    QLA_Complex *v = &CRarg.c[idx];

    QLA_real(*r) = -QLA_real(*v) + CRarg.r[idx];
    QLA_imag(*r) = -QLA_imag(*v);
}

static void
doCRsub(QLA_Complex *r, int idx)
{
    QLA_Complex *v = &CRarg.c[idx];

    QLA_real(*r) = QLA_real(*v) - CRarg.r[idx];
    QLA_imag(*r) = QLA_imag(*v);
}

static int
q_C_sub_R(lua_State *L)
{
    mLatComplex *a = qlua_checkLatComplex(L, 1);
    mLatReal *b = qlua_checkLatReal(L, 2);
    mLatComplex *r = qlua_newLatComplex(L);

    CALL_QDP(L);
    X_C_eq_R_op_C(r->ptr, b->ptr, a->ptr, doCRsub, *qCurrent);

    return 1;
}

static int
q_R_sub_C(lua_State *L)
{
    mLatReal *a = qlua_checkLatReal(L, 1);
    mLatComplex *b = qlua_checkLatComplex(L, 2);
    mLatComplex *r = qlua_newLatComplex(L);

    CALL_QDP(L);
    X_C_eq_R_op_C(r->ptr, a->ptr, b->ptr, doRCsub, *qCurrent);

    return 1;
}

static void
dorCsub(QLA_Complex *r, int idx)
{
    QLA_Complex *v = &CRarg.c[idx];

    QLA_real(*r) = -QLA_real(*v) + *CRarg.r;
    QLA_imag(*r) = -QLA_imag(*v);
}

static int
q_r_sub_C(lua_State *L)
{
    QLA_Real a = luaL_checknumber(L, 1);
    mLatComplex *b = qlua_checkLatComplex(L, 2);
    mLatComplex *r = qlua_newLatComplex(L);

    CALL_QDP(L);
    X_C_eq_r_op_C(r->ptr, a, b->ptr, dorCsub, *qCurrent);

    return 1;
}

static int
q_C_sub_r(lua_State *L)
{
    mLatComplex *a = qlua_checkLatComplex(L, 1);
    QLA_Real b = -luaL_checknumber(L, 2);
    mLatComplex *r = qlua_newLatComplex(L);

    CALL_QDP(L);
    X_C_eq_r_op_C(r->ptr, b, a->ptr, dorCadd, *qCurrent);

    return 1;
}

static void
docCsub(QLA_Complex *r, int idx)
{
    QLA_C_eq_C_minus_C(r, cCarg.a, &cCarg.b[idx]);
}

static int
q_c_sub_C(lua_State *L)
{
    QLA_Complex *a = qlua_checkComplex(L, 1);
    mLatComplex *b = qlua_checkLatComplex(L, 2);
    mLatComplex *r = qlua_newLatComplex(L);

    X_C_eq_c_op_C(r->ptr, a, b->ptr, docCsub, *qCurrent);

    return 1;
}

static int
q_C_sub_c(lua_State *L)
{
    mLatComplex *a = qlua_checkLatComplex(L, 1);
    QLA_Complex *b = qlua_checkComplex(L, 2);
    mLatComplex *r = qlua_newLatComplex(L);
    QLA_Complex v;

    CALL_QDP(L);
    QLA_C_eqm_C(&v, b);
    X_C_eq_c_op_C(r->ptr, &v, a->ptr, docCadd, *qCurrent);

    return 1;
}

static void
docRsub(QLA_Complex *r, int idx)
{
    QLA_real(*r) = QLA_real(*cRarg.a) - cRarg.b[idx];
    QLA_imag(*r) = QLA_imag(*cRarg.a);
}

static int
q_c_sub_R(lua_State *L)
{
    QLA_Complex *a = qlua_checkComplex(L, 1);
    mLatReal *b = qlua_checkLatReal(L, 2);
    mLatComplex *r = qlua_newLatComplex(L);

    CALL_QDP(L);
    X_C_eq_c_op_R(r->ptr, a, b->ptr, docRsub, *qCurrent);

    return 1;
}

static int
q_R_sub_c(lua_State *L)
{
    mLatReal *a = qlua_checkLatReal(L, 1);
    QLA_Complex *b = qlua_checkComplex(L, 2);
    mLatComplex *r = qlua_newLatComplex(L);
    QLA_Complex v;

    CALL_QDP(L);
    QLA_C_eqm_C(&v, b);
    X_C_eq_c_op_R(r->ptr, &v, a->ptr, docRadd, *qCurrent);

    return 1;
}

static int
q_C_mul_C(lua_State *L)
{
    mLatComplex *a = qlua_checkLatComplex(L, 1);
    mLatComplex *b = qlua_checkLatComplex(L, 2);
    mLatComplex *c = qlua_newLatComplex(L);

    CALL_QDP(L);
    QDP_C_eq_C_times_C(c->ptr, a->ptr, b->ptr, *qCurrent);

    return 1;
}

static int
q_C_mul_c(lua_State *L)
{
    mLatComplex *a = qlua_checkLatComplex(L, 1);
    QLA_Complex *b = qlua_checkComplex(L, 2);
    mLatComplex *c = qlua_newLatComplex(L);

    CALL_QDP(L);
    QDP_C_eq_c_times_C(c->ptr, b, a->ptr, *qCurrent);

    return 1;
}

static int
q_c_mul_C(lua_State *L)
{
    QLA_Complex *a = qlua_checkComplex(L, 1);
    mLatComplex *b = qlua_checkLatComplex(L, 2);
    mLatComplex *c = qlua_newLatComplex(L);

    CALL_QDP(L);
    QDP_C_eq_c_times_C(c->ptr, a, b->ptr, *qCurrent);

    return 1;
}

static int
q_C_mul_r(lua_State *L)
{
    mLatComplex *a = qlua_checkLatComplex(L, 1);
    QLA_Real b = luaL_checknumber(L, 2);
    mLatComplex *c = qlua_newLatComplex(L);

    CALL_QDP(L);
    QDP_C_eq_r_times_C(c->ptr, &b, a->ptr, *qCurrent);

    return 1;
}

static int
q_r_mul_C(lua_State *L)
{
    QLA_Real a = luaL_checknumber(L, 1);
    mLatComplex *b = qlua_checkLatComplex(L, 2);
    mLatComplex *c = qlua_newLatComplex(L);

    CALL_QDP(L);
    QDP_C_eq_r_times_C(c->ptr, &a, b->ptr, *qCurrent);

    return 1;
}

static void
doRCmul(QLA_Complex *r, int idx)
{
    QLA_Complex *v = &CRarg.c[idx];

    QLA_real(*r) = QLA_real(*v) * CRarg.r[idx];
    QLA_imag(*r) = QLA_imag(*v) * CRarg.r[idx];
}

static int
q_C_mul_R(lua_State *L)
{
    mLatComplex *a = qlua_checkLatComplex(L, 1);
    mLatReal *b = qlua_checkLatReal(L, 2);
    mLatComplex *r = qlua_newLatComplex(L);

    CALL_QDP(L);
    X_C_eq_R_op_C(r->ptr, b->ptr, a->ptr, doRCmul, *qCurrent);

    return 1;
}

static int
q_R_mul_C(lua_State *L)
{
    mLatReal *a = qlua_checkLatReal(L, 1);
    mLatComplex *b = qlua_checkLatComplex(L, 2);
    mLatComplex *r = qlua_newLatComplex(L);

    CALL_QDP(L);
    X_C_eq_R_op_C(r->ptr, a->ptr, b->ptr, doRCmul, *qCurrent);

    return 1;
}

static void
docRmul(QLA_Complex *r, int idx)
{
    QLA_real(*r) = QLA_real(*cRarg.a) * cRarg.b[idx];
    QLA_imag(*r) = QLA_imag(*cRarg.a) * cRarg.b[idx];
}

static int
q_c_mul_R(lua_State *L)
{
    QLA_Complex *a = qlua_checkComplex(L, 1);
    mLatReal *b = qlua_checkLatReal(L, 2);
    mLatComplex *r = qlua_newLatComplex(L);

    CALL_QDP(L);
    X_C_eq_c_op_R(r->ptr, a, b->ptr, docRmul, *qCurrent);

    return 1;
}

static int
q_R_mul_c(lua_State *L)
{
    mLatReal *a = qlua_checkLatReal(L, 1);
    QLA_Complex *b = qlua_checkComplex(L, 2);
    mLatComplex *r = qlua_newLatComplex(L);

    CALL_QDP(L);
    X_C_eq_c_op_R(r->ptr, b, a->ptr, docRmul, *qCurrent);

    return 1;
}

static int
q_C_div_C(lua_State *L)
{
    mLatComplex *a = qlua_checkLatComplex(L, 1);
    mLatComplex *b = qlua_checkLatComplex(L, 2);
    mLatComplex *c = qlua_newLatComplex(L);

    CALL_QDP(L);
    QDP_C_eq_C_divide_C(c->ptr, a->ptr, b->ptr, *qCurrent);

    return 1;
}

static void
doRCdiv(QLA_Complex *r, int idx)
{
    QLA_Complex *v = &CRarg.c[idx];
    double n = 1/hypot(QLA_real(*v), QLA_imag(*v));

    QLA_Complex w;

    QLA_real(w) = (QLA_real(*v)*n) * n;
    QLA_imag(w) = (QLA_imag(*v)*n) * n;

    QLA_real(*r) = CRarg.r[idx] * QLA_real(w);
    QLA_imag(*r) = -CRarg.r[idx] * QLA_imag(*v);
}

static void
doCRdiv(QLA_Complex *r, int idx)
{
    QLA_Complex *v = &CRarg.c[idx];
    double n = 1 / CRarg.r[idx];
    
    QLA_real(*r) = QLA_real(*v) * n;
    QLA_imag(*r) = QLA_imag(*v) * n;
}

static int
q_C_div_R(lua_State *L)
{
    mLatComplex *a = qlua_checkLatComplex(L, 1);
    mLatReal *b = qlua_checkLatReal(L, 2);
    mLatComplex *r = qlua_newLatComplex(L);
    
    CALL_QDP(L);
    X_C_eq_R_op_C(r->ptr, b->ptr, a->ptr, doCRdiv, *qCurrent);

    return 1;
}

static int
q_R_div_C(lua_State *L)
{
    mLatReal *a = qlua_checkLatReal(L, 1);
    mLatComplex *b = qlua_checkLatComplex(L, 2);
    mLatComplex *r = qlua_newLatComplex(L);

    CALL_QDP(L);
    X_C_eq_R_op_C(r->ptr, a->ptr, b->ptr, doRCdiv, *qCurrent);

    return 1;
}

static void
dorCdiv(QLA_Complex *r, int idx)
{
    QLA_Complex *v = &CRarg.c[idx];
    double n = 1/hypot(QLA_real(*v), QLA_imag(*v));

    QLA_Complex w;

    QLA_real(w) = (QLA_real(*v)*n) * n;
    QLA_imag(w) = (QLA_imag(*v)*n) * n;

    QLA_real(*r) = *CRarg.r * QLA_real(w);
    QLA_imag(*r) = -*CRarg.r * QLA_imag(*v);
}

static int
q_r_div_C(lua_State *L)
{
    QLA_Real a = luaL_checknumber(L, 1);
    mLatComplex *b = qlua_checkLatComplex(L, 2);
    mLatComplex *r = qlua_newLatComplex(L);

    CALL_QDP(L);
    X_C_eq_r_op_C(r->ptr, a, b->ptr, dorCdiv, *qCurrent);

    return 1;
}

static int
q_C_div_r(lua_State *L)
{
    mLatComplex *a = qlua_checkLatComplex(L, 1);
    QLA_Real b = 1 / luaL_checknumber(L, 2);
    mLatComplex *r = qlua_newLatComplex(L);

    CALL_QDP(L);
    QDP_C_eq_r_times_C(r->ptr, &b, a->ptr, *qCurrent);

    return 1;
}

static void
docCdiv(QLA_Complex *r, int idx)
{
    QLA_C_eq_C_divide_C(r, cCarg.a, &cCarg.b[idx]);
}

static int
q_c_div_C(lua_State *L)
{
    QLA_Complex *a = qlua_checkComplex(L, 1);
    mLatComplex *b = qlua_checkLatComplex(L, 2);
    mLatComplex *r = qlua_newLatComplex(L);

    CALL_QDP(L);
    X_C_eq_c_op_C(r->ptr, a, b->ptr, docCdiv, *qCurrent);

    return 1;
}

static int
q_C_div_c(lua_State *L)
{
    mLatComplex *a = qlua_checkLatComplex(L, 1);
    QLA_Complex *b = qlua_checkComplex(L, 2);
    mLatComplex *r = qlua_newLatComplex(L);
    double h = 1/hypot(QLA_real(*b), QLA_imag(*b));
    QLA_Complex v;

    QLA_real(v) = (QLA_real(*b) * h) * h;
    QLA_imag(v) = -(QLA_imag(*b) * h) * h;

    QDP_C_eq_c_times_C(r->ptr, &v, a->ptr, *qCurrent);

    return 1;
}

static void
docRdiv(QLA_Complex *r, int idx)
{
    QLA_real(*r) = QLA_real(*cRarg.a) / cRarg.b[idx];
    QLA_imag(*r) = QLA_imag(*cRarg.a) / cRarg.b[idx];
}

static int
q_c_div_R(lua_State *L)
{
    QLA_Complex *a = qlua_checkComplex(L, 1);
    mLatReal *b = qlua_checkLatReal(L, 2);
    mLatComplex *r = qlua_newLatComplex(L);

    CALL_QDP(L);
    X_C_eq_c_op_R(r->ptr, a, b->ptr, docRdiv, *qCurrent);

    return 1;
}

static int
q_R_div_c(lua_State *L)
{
    mLatReal *a = qlua_checkLatReal(L, 1);
    QLA_Complex *b = qlua_checkComplex(L, 2);
    mLatComplex *r = qlua_newLatComplex(L);
    double h = 1/hypot(QLA_real(*b), QLA_imag(*b));
    QLA_Complex v;

    CALL_QDP(L);
    QLA_real(v) = (QLA_real(*b) * h) * h;
    QLA_imag(v) = -(QLA_imag(*b) * h) * h;
    X_C_eq_c_op_R(r->ptr, &v, a->ptr, docRmul, *qCurrent);

    return 1;
}

static int
q_C_project(lua_State *L)
{
    mLatComplex *a = qlua_checkLatComplex(L, 1);
    mLatComplex *b = qlua_checkLatComplex(L, 2);
    mLatMulti *m = qlua_checkLatMulti(L, 3);
    int size = m->size;
    mVecComplex *r = qlua_newVecComplex(L, size);
    QLA_Real *rr = qlua_malloc(L, 2 * size * sizeof (QLA_Real));
    QLA_Int *ii = m->idx;
    QLA_Complex *aa;
    QLA_Complex *bb;
    int k;
    
    r->size = size;
    for (k = 0; k < 2 * size; k++)
        rr[k] = 0;
    CALL_QDP(L);
    aa = QDP_expose_C(a->ptr);
    bb = QDP_expose_C(b->ptr);
    for (k = 0; k < QDP_sites_on_node; k++, ii++, aa++, bb++) {
        int idx = *ii;
        double var, vai, vbr, vbi;
        if ((idx < 0) || (idx >= size))
            continue;
        var = QLA_real(*aa);
        vai = QLA_imag(*aa);
        vbr = QLA_real(*bb);
        vbi = QLA_imag(*bb);
        rr[2 * idx    ] += var * vbr - vai * vbi;
        rr[2 * idx + 1] += var * vbi + vai * vbr;
    }
    QDP_reset_C(a->ptr);
    QDP_reset_C(b->ptr);
    QMP_sum_double_array(rr, 2 * size);
    for (k = 0; k < size; k++) {
        QLA_c_eq_r_plus_ir(r->val[k], rr[2 * k], rr[2 * k + 1]);
    }

    return 1;
}

static int
q_C_sum(lua_State *L)
{
    mLatComplex *a = qlua_checkLatComplex(L, 1);

    switch (lua_gettop(L)) {
    case 1: {
        QLA_Complex *s = qlua_newComplex(L);

        CALL_QDP(L);
        QDP_c_eq_sum_C(s, a->ptr, *qCurrent);

        return 1;
    }
    case 2: {
        mLatMulti *m = qlua_checkLatMulti(L, 2);
        int size = m->size;
        QLA_Int *ii = m->idx;
        mVecComplex *r = qlua_newVecComplex(L, size);
        int k;
        QLA_Complex *xx;
        QLA_Real *rr;
        
        r->size = size;
        rr = qlua_malloc(L, size * 2 * sizeof (QLA_Real));
        for (k = 0; k < 2 * size; k++)
            rr[k] = 0;

        CALL_QDP(L);
        xx = QDP_expose_C(a->ptr);
        for (k = 0; k < QDP_sites_on_node; k++, xx++, ii++) {
            int t = *ii;
            if ((t < 0) || (t >= size))
                continue;
            rr[2 * t] += QLA_real(*xx);
            rr[2 * t + 1] += QLA_imag(*xx);
        }
        QDP_reset_C(a->ptr);
        QMP_sum_double_array(rr, 2 * size);
        for (k = 0; k < size; k++) {
            QLA_c_eq_r_plus_ir(r->val[k], rr[2 * k], rr[2 * k + 1]);
        }
        qlua_free(L, rr);

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

    CALL_QDP(L);
    QDP_r_eq_norm2_C(&n, a->ptr, *qCurrent);
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

    CALL_QDP(L);
    QDP_C_eq_sC(r->ptr, a->ptr, shift, dir, *qCurrent);

    return 1;
}

static int
q_C_conj(lua_State *L)
{
    mLatComplex *a = qlua_checkLatComplex(L, 1);
    mLatComplex *r = qlua_newLatComplex(L);

    CALL_QDP(L);
    QDP_C_eq_conj_C(r->ptr, a->ptr, *qCurrent);

    return 1;
}

static int
q_C_abs(lua_State *L)
{
    mLatComplex *a = qlua_checkLatComplex(L, 1);
    mLatReal *r = qlua_newLatReal(L);

    CALL_QDP(L);
    QDP_R_eq_norm_C(r->ptr, a->ptr, *qCurrent);

    return 1;
}

static int
q_C_arg(lua_State *L)
{
    mLatComplex *a = qlua_checkLatComplex(L, 1);
    mLatReal *r = qlua_newLatReal(L);

    CALL_QDP(L);
    QDP_R_eq_arg_C(r->ptr, a->ptr, *qCurrent);

    return 1;
}

static int
q_C_sqrt(lua_State *L)
{
    mLatComplex *a = qlua_checkLatComplex(L, 1);
    mLatComplex *r = qlua_newLatComplex(L);

    CALL_QDP(L);
    QDP_C_eq_csqrt_C(r->ptr, a->ptr, *qCurrent);

    return 1;
}

static int
q_C_exp(lua_State *L)
{
    mLatComplex *a = qlua_checkLatComplex(L, 1);
    mLatComplex *r = qlua_newLatComplex(L);

    CALL_QDP(L);
    QDP_C_eq_cexp_C(r->ptr, a->ptr, *qCurrent);

    return 1;
}

static int
q_C_log(lua_State *L)
{
    mLatComplex *a = qlua_checkLatComplex(L, 1);
    mLatComplex *r = qlua_newLatComplex(L);

    CALL_QDP(L);
    QDP_C_eq_clog_C(r->ptr, a->ptr, *qCurrent);

    return 1;
}

static int
q_C_set(lua_State *L)
{
    mLatComplex *r = qlua_checkLatComplex(L, 1);
    mLatComplex *a = qlua_checkLatComplex(L, 2);

    CALL_QDP(L);
    QDP_C_eq_C(r->ptr, a->ptr, *qCurrent);
    lua_pop(L, 1);

    return 1;
}

int
q_C_gaussian(lua_State *L)
{
    mLatRandom *a = qlua_checkLatRandom(L, 1);
    mLatComplex *r = qlua_newLatComplex(L);

    CALL_QDP(L);
    QDP_C_eq_gaussian_S(r->ptr, a->ptr, *qCurrent);

    return 1;
}

static int
q_C_dot(lua_State *L)
{
    mLatComplex *a = qlua_checkLatComplex(L, 1);
    mLatComplex *b = qlua_checkLatComplex(L, 2);
    mLatComplex *s = qlua_newLatComplex(L);

    CALL_QDP(L);
    QDP_C_eq_C_dot_C(s->ptr, a->ptr, b->ptr, *qCurrent);

    return 1;
}

static int
q_latcomplex(lua_State *L)
{
    switch (lua_gettop(L)) {
    case 1: {
        mLatComplex *v = qlua_newLatComplex(L);
        
        CALL_QDP(L);
        QDP_C_eq_zero(v->ptr, *qCurrent);

        return 1;
    }
    case 2:
        switch (qlua_gettype(L, 2)) {
        case qReal: {
            QLA_Real d = luaL_checknumber(L, 2);
            QLA_Complex z;
            mLatComplex *v = qlua_newLatComplex(L);
            
            CALL_QDP(L);
            QLA_real(z) = d;
            QLA_imag(z) = 0.0;
            QDP_C_eq_c(v->ptr, &z, *qCurrent);

            return 1;
        }
        case qComplex: {
            QLA_Complex *z = qlua_checkComplex(L, 2);
            mLatComplex *v = qlua_newLatComplex(L);

            CALL_QDP(L);
            QDP_C_eq_c(v->ptr, z, *qCurrent);
            
            return 1;
        }
        case qLatReal: {
            mLatReal *d = qlua_checkLatReal(L, 2);
            mLatComplex *v = qlua_newLatComplex(L);
            
            CALL_QDP(L);
            QDP_C_eq_R(v->ptr, d->ptr, *qCurrent);

            return 1;
        }
        case qLatComplex: {
            mLatComplex *d = qlua_checkLatComplex(L, 2);
            mLatComplex *v = qlua_newLatComplex(L);
            
            CALL_QDP(L);
            QDP_C_eq_C(v->ptr, d->ptr, *qCurrent);

            return 1;
        }
        default:
            break;
        }
    case 3: {
        mLatReal *a = qlua_checkLatReal(L, 2);
        mLatReal *b = qlua_checkLatReal(L, 3);
        mLatComplex *c = qlua_newLatComplex(L);

        CALL_QDP(L);
        QDP_C_eq_R_plus_i_R(c->ptr, a->ptr, b->ptr, *qCurrent);

        return 1;
    }
    }
    return qlua_badconstr(L, "Complex");
}

static struct luaL_Reg LatComplexMethods[] = {
    { "real",     q_C_real },
    { "imag",     q_C_imag },
    { "sum",      q_C_sum },
    { "norm2",    q_C_norm2 },
    { "shift",    q_C_shift },
    { "conj",     q_C_conj },
    { "abs",      q_C_abs },
    { "arg",      q_C_arg },
    { "sqrt",     q_C_sqrt },
    { "exp",      q_C_exp },
    { "log",      q_C_log },
    { "set",      q_C_set },
    { "project",  q_C_project },
    { NULL,       NULL}
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
    qlua_reg_add(qLatComplex, qLatReal,    q_C_add_R);
    qlua_reg_add(qLatReal,    qLatComplex, q_R_add_C);
    qlua_reg_add(qLatComplex, qReal,       q_C_add_r);
    qlua_reg_add(qReal,       qLatComplex, q_r_add_C);
    qlua_reg_add(qLatComplex, qComplex,    q_C_add_c);
    qlua_reg_add(qComplex,    qLatComplex, q_c_add_C);
    qlua_reg_add(qComplex,    qLatReal,    q_c_add_R);
    qlua_reg_add(qLatReal,    qComplex,    q_R_add_c);

    qlua_reg_sub(qLatComplex, qLatComplex, q_C_sub_C);
    qlua_reg_sub(qLatComplex, qLatReal,    q_C_sub_R);
    qlua_reg_sub(qLatReal,    qLatComplex, q_R_sub_C);
    qlua_reg_sub(qLatComplex, qReal,       q_C_sub_r);
    qlua_reg_sub(qReal,       qLatComplex, q_r_sub_C);
    qlua_reg_sub(qLatComplex, qComplex,    q_C_sub_c);
    qlua_reg_sub(qComplex,    qLatComplex, q_c_sub_C);
    qlua_reg_sub(qComplex,    qLatReal,    q_c_sub_R);
    qlua_reg_sub(qLatReal,    qComplex,    q_R_sub_c);

    qlua_reg_mul(qLatComplex, qLatComplex, q_C_mul_C);
    qlua_reg_mul(qLatComplex, qLatReal,    q_C_mul_R);
    qlua_reg_mul(qLatReal,    qLatComplex, q_R_mul_C);
    qlua_reg_mul(qLatComplex, qReal,       q_C_mul_r);
    qlua_reg_mul(qReal,       qLatComplex, q_r_mul_C);
    qlua_reg_mul(qLatComplex, qComplex,    q_C_mul_c);
    qlua_reg_mul(qComplex,    qLatComplex, q_c_mul_C);
    qlua_reg_mul(qComplex,    qLatReal,    q_c_mul_R);
    qlua_reg_mul(qLatReal,    qComplex,    q_R_mul_c);

    qlua_reg_div(qLatComplex, qLatComplex, q_C_div_C);
    qlua_reg_div(qLatComplex, qLatReal,    q_C_div_R);
    qlua_reg_div(qLatReal,    qLatComplex, q_R_div_C);
    qlua_reg_div(qLatComplex, qReal,       q_C_div_r);
    qlua_reg_div(qReal,       qLatComplex, q_r_div_C);
    qlua_reg_div(qLatComplex, qComplex,    q_C_div_c);
    qlua_reg_div(qComplex,    qLatComplex, q_c_div_C);
    qlua_reg_div(qComplex,    qLatReal,    q_c_div_R);
    qlua_reg_div(qLatReal,    qComplex,    q_R_div_c);

    qlua_reg_dot(qLatComplex, q_C_dot);

    return 0;
}

int
fini_latcomplex(lua_State *L)
{
    return 0;
}
