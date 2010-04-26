#include "qlua.h"                                                    /* DEPS */
#include "qcomplex.h"                                                /* DEPS */
#include "qvector.h"                                                 /* DEPS */
#include "lattice.h"                                                 /* DEPS */
#include "latsubset.h"                                               /* DEPS */
#include "latcomplex.h"                                              /* DEPS */
#include "latreal.h"                                                 /* DEPS */
#include "latint.h"                                                  /* DEPS */
#include "latrandom.h"                                               /* DEPS */
#include "latmulti.h"                                                /* DEPS */
#include <math.h>
#include "qmp.h"

static const char LatComplexName[] = "lattice.Complex";

static int
q_C_fmt(lua_State *L)
{
    char fmt[72];
    mLatComplex *b = qlua_checkLatComplex(L, 1, NULL);

    sprintf(fmt, "LatComplex(%p)", b->ptr);
    lua_pushstring(L, fmt);

    return 1;
}

static int
q_C_gc(lua_State *L)
{
    mLatComplex *b = qlua_checkLatComplex(L, 1, NULL);

    QDP_destroy_C(b->ptr);
    b->ptr = 0;

    return 0;
}

static int
q_C_neg(lua_State *L)
{
    mLatComplex *a = qlua_checkLatComplex(L, 1, NULL);
    mLattice *S = qlua_ObjLattice(L, 1);
    mLatComplex *r = qlua_newLatComplex(L, lua_gettop(L));
    QLA_Real m1 = -1;

    CALL_QDP(L);
    QDP_C_eq_r_times_C(r->ptr, &m1, a->ptr, *S->qss);
    return 1;
}

static int
q_C_get(lua_State *L)
{
    switch (qlua_qtype(L, 2)) {
    case qTable: {
        mLatComplex *V = qlua_checkLatComplex(L, 1, NULL);
        mLattice *S = qlua_ObjLattice(L, 1);
        QLA_Complex *W;
        QLA_Complex *locked;
        double zri[2];
        int *idx = 0;
        
        idx = qlua_checklatcoord(L, 2, S);
        CALL_QDP(L);
        locked = QDP_expose_C(V->ptr);
        int site_node = QDP_node_number_L(S->lat, idx);
        if (site_node == QDP_this_node) {
            QLA_Complex *zz = &QLA_elem_C(locked[QDP_index_L(S->lat, idx)]);

            zri[0] = QLA_real(*zz);
            zri[1] = QLA_imag(*zz);
        }
        QDP_reset_C(V->ptr);
        qlua_free(L, idx);
        XMP_dist_double_array(site_node, 2, zri);
        W = qlua_newComplex(L);
        QLA_c_eq_r_plus_ir(*W, zri[0], zri[1]);

        return 1;
    }
    case qString:
        return qlua_selflookup(L, 1, luaL_checkstring(L, 2));
    default:
        break;
    }
    return qlua_badindex(L, "Complex");
}

static int
q_C_put(lua_State *L)
{
    mLatComplex *V = qlua_checkLatComplex(L, 1, NULL);
    mLattice *S = qlua_ObjLattice(L, 1);
    QLA_Complex *locked;
    int *idx = 0;
    double z_re, z_im;

    switch (qlua_qtype(L, 3)) {
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
    idx = qlua_checklatcoord(L, 2, S);
    CALL_QDP(L);
    locked = QDP_expose_C(V->ptr);
    if (QDP_node_number_L(S->lat, idx) == QDP_this_node) {
        QLA_Complex *zz = &QLA_elem_C(locked[QDP_index_L(S->lat, idx)]);

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
    mLatComplex *a = qlua_checkLatComplex(L, 1, NULL);
    mLattice *S = qlua_ObjLattice(L, 1);
    mLatReal *c = qlua_newLatReal(L, lua_gettop(L));

    CALL_QDP(L);
    QDP_R_eq_re_C(c->ptr, a->ptr, *S->qss);

    return 1;
}

static int
q_C_imag(lua_State *L)
{
    mLatComplex *a = qlua_checkLatComplex(L, 1, NULL);
    mLattice *S = qlua_ObjLattice(L, 1);
    mLatReal *c = qlua_newLatReal(L, lua_gettop(L));

    CALL_QDP(L);
    QDP_R_eq_im_C(c->ptr, a->ptr, *S->qss);

    return 1;
}

static int
q_C_add_C(lua_State *L)
{
    mLatComplex *a = qlua_checkLatComplex(L, 1, NULL);
    mLattice *S = qlua_ObjLattice(L, 1);
    mLatComplex *b = qlua_checkLatComplex(L, 2, S);
    mLatComplex *c = qlua_newLatComplex(L, lua_gettop(L));

    CALL_QDP(L);
    QDP_C_eq_C_plus_C(c->ptr, a->ptr, b->ptr, *S->qss);

    return 1;
}

static struct {
    QLA_Real *r;
    QLA_Complex *c;
} CRarg; /* YYY global state */

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
    mLatComplex *a = qlua_checkLatComplex(L, 1, NULL);
    mLattice *S = qlua_ObjLattice(L, 1);
    mLatReal *b = qlua_checkLatReal(L, 2, S);
    mLatComplex *r = qlua_newLatComplex(L, lua_gettop(L));

    CALL_QDP(L);

    CRarg.r = QDP_expose_R(b->ptr);
    CRarg.c = QDP_expose_C(a->ptr);
    QDP_C_eq_funci(r->ptr, doRCadd, *S->qss);
    CRarg.r = 0;
    CRarg.c = 0;
    QDP_reset_R(b->ptr);
    QDP_reset_C(a->ptr);

    return 1;
}

static int
q_R_add_C(lua_State *L)
{
    mLatReal *a = qlua_checkLatReal(L, 1, NULL);
    mLattice *S = qlua_ObjLattice(L, 1);
    mLatComplex *b = qlua_checkLatComplex(L, 2, S);
    mLatComplex *r = qlua_newLatComplex(L, lua_gettop(L));

    CALL_QDP(L);
    CRarg.r = QDP_expose_R(a->ptr);
    CRarg.c = QDP_expose_C(b->ptr);
    QDP_C_eq_funci(r->ptr, doRCadd, *S->qss);
    CRarg.r = 0;
    CRarg.c = 0;
    QDP_reset_R(a->ptr);
    QDP_reset_C(b->ptr);

    return 1;
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
    mLatComplex *b = qlua_checkLatComplex(L, 2, NULL);
    mLattice *S = qlua_ObjLattice(L, 2);
    mLatComplex *r = qlua_newLatComplex(L, lua_gettop(L));

    CALL_QDP(L);
    CRarg.r = &a;
    CRarg.c = QDP_expose_C(b->ptr);
    QDP_C_eq_funci(r->ptr, dorCadd, *S->qss);
    CRarg.r = 0;
    CRarg.c = 0;
    QDP_reset_C(b->ptr);

    return 1;
}

static int
q_C_add_r(lua_State *L)
{
    mLatComplex *a = qlua_checkLatComplex(L, 1, NULL);
    mLattice *S = qlua_ObjLattice(L, 1);
    QLA_Real b = luaL_checknumber(L, 2);
    mLatComplex *r = qlua_newLatComplex(L, lua_gettop(L));

    CALL_QDP(L);
    CRarg.r = &b;
    CRarg.c = QDP_expose_C(a->ptr);
    QDP_C_eq_funci(r->ptr, dorCadd, *S->qss);
    CRarg.r = 0;
    CRarg.c = 0;
    QDP_reset_C(a->ptr);

    return 1;
}

static struct {
    QLA_Complex *a;
    QLA_Complex *b;
} cCarg; /* YYY global state */

static void
docCadd(QLA_Complex *r, int idx)
{
    QLA_C_eq_C_plus_C(r, cCarg.a, &cCarg.b[idx]);
}

static int
q_c_add_C(lua_State *L)
{
    QLA_Complex *a = qlua_checkComplex(L, 1);
    mLatComplex *b = qlua_checkLatComplex(L, 2, NULL);
    mLattice *S = qlua_ObjLattice(L, 2);
    mLatComplex *r = qlua_newLatComplex(L, lua_gettop(L));

    CALL_QDP(L);
    cCarg.a = a;
    cCarg.b = QDP_expose_C(b->ptr);
    QDP_C_eq_funci(r->ptr, docCadd, *S->qss);
    cCarg.a = 0;
    cCarg.b = 0;
    QDP_reset_C(b->ptr);

    return 1;
}

static int
q_C_add_c(lua_State *L)
{
    mLatComplex *a = qlua_checkLatComplex(L, 1, NULL);
    mLattice *S = qlua_ObjLattice(L, 1);
    QLA_Complex *b = qlua_checkComplex(L, 2);
    mLatComplex *r = qlua_newLatComplex(L, lua_gettop(L));

    CALL_QDP(L);
    cCarg.a = b;
    cCarg.b = QDP_expose_C(a->ptr);
    QDP_C_eq_funci(r->ptr, docCadd, *S->qss);
    cCarg.a = 0;
    cCarg.b = 0;
    QDP_reset_C(a->ptr);

    return 1;
}

static struct {
    QLA_Complex *a;
    QLA_Real *b;
} cRarg; /* YYY global state */

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
    mLatReal *b = qlua_checkLatReal(L, 2, NULL);
    mLattice *S = qlua_ObjLattice(L, 2);
    mLatComplex *r = qlua_newLatComplex(L, lua_gettop(L));

    CALL_QDP(L);
    cRarg.a = a;
    cRarg.b = QDP_expose_R(b->ptr);
    QDP_C_eq_funci(r->ptr, docRadd, *S->qss);
    cRarg.a = 0;
    cRarg.b = 0;
    QDP_reset_R(b->ptr);

    return 1;
}

static int
q_R_add_c(lua_State *L)
{
    mLatReal *a = qlua_checkLatReal(L, 1, NULL);
    mLattice *S = qlua_ObjLattice(L, 1);
    QLA_Complex *b = qlua_checkComplex(L, 2);
    mLatComplex *r = qlua_newLatComplex(L, lua_gettop(L));

    CALL_QDP(L);
    cRarg.a = b;
    cRarg.b = QDP_expose_R(a->ptr);
    QDP_C_eq_funci(r->ptr, docRadd, *S->qss);
    cRarg.a = 0;
    cRarg.b = 0;
    QDP_reset_R(a->ptr);

    return 1;
}

static int
q_C_sub_C(lua_State *L)
{
    mLatComplex *a = qlua_checkLatComplex(L, 1, NULL);
    mLattice *S = qlua_ObjLattice(L, 1);
    mLatComplex *b = qlua_checkLatComplex(L, 2, S);
    mLatComplex *c = qlua_newLatComplex(L, lua_gettop(L));

    CALL_QDP(L);
    QDP_C_eq_C_minus_C(c->ptr, a->ptr, b->ptr, *S->qss);

    return 1;
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
    mLatComplex *a = qlua_checkLatComplex(L, 1, NULL);
    mLattice *S = qlua_ObjLattice(L, 1);
    mLatReal *b = qlua_checkLatReal(L, 2, S);
    mLatComplex *r = qlua_newLatComplex(L, lua_gettop(L));

    CALL_QDP(L);
    CRarg.r = QDP_expose_R(b->ptr);
    CRarg.c = QDP_expose_C(a->ptr);
    QDP_C_eq_funci(r->ptr, doCRsub, *S->qss);
    CRarg.r = 0;
    CRarg.c = 0;
    QDP_reset_C(a->ptr);
    QDP_reset_R(b->ptr);

    return 1;
}

static void
doRCsub(QLA_Complex *r, int idx)
{
    QLA_Complex *v = &CRarg.c[idx];

    QLA_real(*r) = -QLA_real(*v) + CRarg.r[idx];
    QLA_imag(*r) = -QLA_imag(*v);
}

static int
q_R_sub_C(lua_State *L)
{
    mLatReal *a = qlua_checkLatReal(L, 1, NULL);
    mLattice *S = qlua_ObjLattice(L, 1);
    mLatComplex *b = qlua_checkLatComplex(L, 2, S);
    mLatComplex *r = qlua_newLatComplex(L, lua_gettop(L));

    CALL_QDP(L);
    CRarg.r = QDP_expose_R(a->ptr);
    CRarg.c = QDP_expose_C(b->ptr);
    QDP_C_eq_funci(r->ptr, doRCsub, *S->qss);
    CRarg.r = 0;
    CRarg.c = 0;
    QDP_reset_R(a->ptr);
    QDP_reset_C(b->ptr);

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
    mLatComplex *b = qlua_checkLatComplex(L, 2, NULL);
    mLattice *S = qlua_ObjLattice(L, 2);
    mLatComplex *r = qlua_newLatComplex(L, lua_gettop(L));

    CALL_QDP(L);
    CRarg.r = &a;
    CRarg.c = QDP_expose_C(b->ptr);
    QDP_C_eq_funci(r->ptr, dorCsub, *S->qss);
    CRarg.r = 0;
    CRarg.c = 0;
    QDP_reset_C(b->ptr);

    return 1;
}

static int
q_C_sub_r(lua_State *L)
{
    mLatComplex *a = qlua_checkLatComplex(L, 1, NULL);
    mLattice *S = qlua_ObjLattice(L, 1);
    QLA_Real b = -luaL_checknumber(L, 2);
    mLatComplex *r = qlua_newLatComplex(L, lua_gettop(L));

    CALL_QDP(L);
    CRarg.r = &b;
    CRarg.c = QDP_expose_C(a->ptr);
    QDP_C_eq_funci(r->ptr, dorCsub, *S->qss);
    CRarg.r = 0;
    CRarg.c = 0;
    QDP_reset_C(a->ptr);

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
    mLatComplex *b = qlua_checkLatComplex(L, 2, NULL);
    mLattice *S = qlua_ObjLattice(L, 2);
    mLatComplex *r = qlua_newLatComplex(L, lua_gettop(L));

    CALL_QDP(L);
    cCarg.a = a;
    cCarg.b = QDP_expose_C(b->ptr);
    QDP_C_eq_funci(r->ptr, docCsub, *S->qss);
    cCarg.a = 0;
    cCarg.b = 0;
    QDP_reset_C(b->ptr);

    return 1;
}

static int
q_C_sub_c(lua_State *L)
{
    mLatComplex *a = qlua_checkLatComplex(L, 1, NULL);
    mLattice *S = qlua_ObjLattice(L, 1);
    QLA_Complex *b = qlua_checkComplex(L, 2);
    mLatComplex *r = qlua_newLatComplex(L, lua_gettop(L));
    QLA_Complex v;

    CALL_QDP(L);
    QLA_C_eqm_C(&v, b);
    cCarg.a = &v;
    cCarg.b = QDP_expose_C(a->ptr);
    QDP_C_eq_funci(r->ptr, docCadd, *S->qss);
    cCarg.a = 0;
    cCarg.b = 0;
    QDP_reset_C(a->ptr);

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
    mLatReal *b = qlua_checkLatReal(L, 2, NULL);
    mLattice *S = qlua_ObjLattice(L, 2);
    mLatComplex *r = qlua_newLatComplex(L, lua_gettop(L));

    CALL_QDP(L);
    cRarg.a = a;
    cRarg.b = QDP_expose_R(b->ptr);
    QDP_C_eq_funci(r->ptr, docRsub, *S->qss);
    cRarg.a = 0;
    cRarg.b = 0;
    QDP_reset_R(b->ptr);

    return 1;
}

static int
q_R_sub_c(lua_State *L)
{
    mLatReal *a = qlua_checkLatReal(L, 1, NULL);
    mLattice *S = qlua_ObjLattice(L, 1);
    QLA_Complex *b = qlua_checkComplex(L, 2);
    mLatComplex *r = qlua_newLatComplex(L, lua_gettop(L));
    QLA_Complex v;

    CALL_QDP(L);
    QLA_C_eqm_C(&v, b);
    cRarg.a = &v;
    cRarg.b = QDP_expose_R(a->ptr);
    QDP_C_eq_funci(r->ptr, docRsub, *S->qss);
    cRarg.a = 0;
    cRarg.b = 0;
    QDP_reset_R(a->ptr);

    return 1;
}

static int
q_C_mul_C(lua_State *L)
{
    mLatComplex *a = qlua_checkLatComplex(L, 1, NULL);
    mLattice *S = qlua_ObjLattice(L, 1);
    mLatComplex *b = qlua_checkLatComplex(L, 2, S);
    mLatComplex *c = qlua_newLatComplex(L, lua_gettop(L));

    CALL_QDP(L);
    QDP_C_eq_C_times_C(c->ptr, a->ptr, b->ptr, *S->qss);

    return 1;
}

static int
q_C_mul_c(lua_State *L)
{
    mLatComplex *a = qlua_checkLatComplex(L, 1, NULL);
    mLattice *S = qlua_ObjLattice(L, 1);
    QLA_Complex *b = qlua_checkComplex(L, 2);
    mLatComplex *c = qlua_newLatComplex(L, lua_gettop(L));

    CALL_QDP(L);
    QDP_C_eq_c_times_C(c->ptr, b, a->ptr, *S->qss);

    return 1;
}

static int
q_c_mul_C(lua_State *L)
{
    QLA_Complex *a = qlua_checkComplex(L, 1);
    mLatComplex *b = qlua_checkLatComplex(L, 2, NULL);
    mLattice *S = qlua_ObjLattice(L, 2);
    mLatComplex *c = qlua_newLatComplex(L, lua_gettop(L));

    CALL_QDP(L);
    QDP_C_eq_c_times_C(c->ptr, a, b->ptr, *S->qss);

    return 1;
}

static int
q_C_mul_r(lua_State *L)
{
    mLatComplex *a = qlua_checkLatComplex(L, 1, NULL);
    mLattice *S = qlua_ObjLattice(L, 1);
    QLA_Real b = luaL_checknumber(L, 2);
    mLatComplex *c = qlua_newLatComplex(L, lua_gettop(L));

    CALL_QDP(L);
    QDP_C_eq_r_times_C(c->ptr, &b, a->ptr, *S->qss);

    return 1;
}

static int
q_r_mul_C(lua_State *L)
{
    QLA_Real a = luaL_checknumber(L, 1);
    mLatComplex *b = qlua_checkLatComplex(L, 2, NULL);
    mLattice *S = qlua_ObjLattice(L, 2);
    mLatComplex *c = qlua_newLatComplex(L, lua_gettop(L));

    CALL_QDP(L);
    QDP_C_eq_r_times_C(c->ptr, &a, b->ptr, *S->qss);

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
    mLatComplex *a = qlua_checkLatComplex(L, 1, NULL);
    mLattice *S = qlua_ObjLattice(L, 1);
    mLatReal *b = qlua_checkLatReal(L, 2, S);
    mLatComplex *r = qlua_newLatComplex(L, lua_gettop(L));

    CALL_QDP(L);
    CRarg.r = QDP_expose_R(b->ptr);
    CRarg.c = QDP_expose_C(a->ptr);
    QDP_C_eq_funci(r->ptr, doRCmul, *S->qss);
    CRarg.r = 0;
    CRarg.c = 0;
    QDP_reset_R(b->ptr);
    QDP_reset_C(a->ptr);

    return 1;
}

static int
q_R_mul_C(lua_State *L)
{
    mLatReal *a = qlua_checkLatReal(L, 1, NULL);
    mLattice *S = qlua_ObjLattice(L, 1);
    mLatComplex *b = qlua_checkLatComplex(L, 2, S);
    mLatComplex *r = qlua_newLatComplex(L, lua_gettop(L));

    CALL_QDP(L);
    CRarg.r = QDP_expose_R(a->ptr);
    CRarg.c = QDP_expose_C(b->ptr);
    QDP_C_eq_funci(r->ptr, doRCmul, *S->qss);
    CRarg.r = 0;
    CRarg.c = 0;
    QDP_reset_R(a->ptr);
    QDP_reset_C(b->ptr);
    
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
    mLatReal *b = qlua_checkLatReal(L, 2, NULL);
    mLattice *S = qlua_ObjLattice(L, 2);
    mLatComplex *r = qlua_newLatComplex(L, lua_gettop(L));

    CALL_QDP(L);
    cRarg.a = a;
    cRarg.b = QDP_expose_R(b->ptr);
    QDP_C_eq_funci(r->ptr, docRmul, *S->qss);
    cRarg.a = 0;
    cRarg.b = 0;
    QDP_reset_R(b->ptr);

    return 1;
}

static int
q_R_mul_c(lua_State *L)
{
    mLatReal *a = qlua_checkLatReal(L, 1, NULL);
    mLattice *S = qlua_ObjLattice(L, 1);
    QLA_Complex *b = qlua_checkComplex(L, 2);
    mLatComplex *r = qlua_newLatComplex(L, lua_gettop(L));

    CALL_QDP(L);
    cRarg.a = b;
    cRarg.b = QDP_expose_R(a->ptr);
    QDP_C_eq_funci(r->ptr, docRmul, *S->qss);
    cRarg.a = 0;
    cRarg.b = 0;
    QDP_reset_R(a->ptr);

    return 1;
}

static int
q_C_div_C(lua_State *L)
{
    mLatComplex *a = qlua_checkLatComplex(L, 1, NULL);
    mLattice *S = qlua_ObjLattice(L, 1);
    mLatComplex *b = qlua_checkLatComplex(L, 2, S);
    mLatComplex *c = qlua_newLatComplex(L, lua_gettop(L));

    CALL_QDP(L);
    QDP_C_eq_C_divide_C(c->ptr, a->ptr, b->ptr, *S->qss);

    return 1;
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
    mLatComplex *a = qlua_checkLatComplex(L, 1, NULL);
    mLattice *S = qlua_ObjLattice(L, 1);
    mLatReal *b = qlua_checkLatReal(L, 2, S);
    mLatComplex *r = qlua_newLatComplex(L, lua_gettop(L));
    
    CALL_QDP(L);
    CRarg.r = QDP_expose_R(b->ptr);
    CRarg.c = QDP_expose_C(a->ptr);
    QDP_C_eq_funci(r->ptr, doCRdiv, *S->qss);
    CRarg.r = 0;
    CRarg.c = 0;
    QDP_reset_R(b->ptr);
    QDP_reset_C(a->ptr);

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

static int
q_R_div_C(lua_State *L)
{
    mLatReal *a = qlua_checkLatReal(L, 1, NULL);
    mLattice *S = qlua_ObjLattice(L, 1);
    mLatComplex *b = qlua_checkLatComplex(L, 2, S);
    mLatComplex *r = qlua_newLatComplex(L, lua_gettop(L));

    CALL_QDP(L);
    CRarg.r = QDP_expose_R(a->ptr);
    CRarg.c = QDP_expose_C(b->ptr);
    QDP_C_eq_funci(r->ptr, doRCdiv, *S->qss);
    CRarg.r = 0;
    CRarg.c = 0;
    QDP_reset_R(a->ptr);
    QDP_reset_C(b->ptr);

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
    mLatComplex *b = qlua_checkLatComplex(L, 2, NULL);
    mLattice *S = qlua_ObjLattice(L, 2);
    mLatComplex *r = qlua_newLatComplex(L, lua_gettop(L));

    CALL_QDP(L);
    CRarg.r = &a;
    CRarg.c = QDP_expose_C(b->ptr);
    QDP_C_eq_funci(r->ptr, dorCdiv, *S->qss);
    CRarg.r = 0;
    CRarg.c = 0;
    QDP_reset_C(b->ptr);

    return 1;
}

static int
q_C_div_r(lua_State *L)
{
    mLatComplex *a = qlua_checkLatComplex(L, 1, NULL);
    mLattice *S = qlua_ObjLattice(L, 1);
    QLA_Real b = 1 / luaL_checknumber(L, 2);
    mLatComplex *r = qlua_newLatComplex(L, lua_gettop(L));

    CALL_QDP(L);
    QDP_C_eq_r_times_C(r->ptr, &b, a->ptr, *S->qss);

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
    mLatComplex *b = qlua_checkLatComplex(L, 2, NULL);
    mLattice *S = qlua_ObjLattice(L, 2);
    mLatComplex *r = qlua_newLatComplex(L, lua_gettop(L));

    CALL_QDP(L);
    cCarg.a = a;
    cCarg.b = QDP_expose_C(b->ptr);
    QDP_C_eq_funci(r->ptr, docCdiv, *S->qss);
    cCarg.a = 0;
    cCarg.b = 0;
    QDP_reset_C(b->ptr);

    return 1;
}

static int
q_C_div_c(lua_State *L)
{
    mLatComplex *a = qlua_checkLatComplex(L, 1, NULL);
    mLattice *S = qlua_ObjLattice(L, 1);
    QLA_Complex *b = qlua_checkComplex(L, 2);
    mLatComplex *r = qlua_newLatComplex(L, lua_gettop(L));
    double h = 1/hypot(QLA_real(*b), QLA_imag(*b));
    QLA_Complex v;

    QLA_real(v) = (QLA_real(*b) * h) * h;
    QLA_imag(v) = -(QLA_imag(*b) * h) * h;

    QDP_C_eq_c_times_C(r->ptr, &v, a->ptr, *S->qss);

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
    mLatReal *b = qlua_checkLatReal(L, 2, NULL);
    mLattice *S = qlua_ObjLattice(L, 2);
    mLatComplex *r = qlua_newLatComplex(L, lua_gettop(L));

    CALL_QDP(L);
    cRarg.a = a;
    cRarg.b = QDP_expose_R(b->ptr);
    QDP_C_eq_funci(r->ptr, docRdiv, *S->qss);
    cRarg.a = 0;
    cRarg.b = 0;
    QDP_reset_R(b->ptr);

    return 1;
}

static int
q_R_div_c(lua_State *L)
{
    mLatReal *a = qlua_checkLatReal(L, 1, NULL);
    mLattice *S = qlua_ObjLattice(L, 1);
    QLA_Complex *b = qlua_checkComplex(L, 2);
    mLatComplex *r = qlua_newLatComplex(L, lua_gettop(L));
    double h = 1/hypot(QLA_real(*b), QLA_imag(*b));
    QLA_Complex v;

    CALL_QDP(L);
    QLA_real(v) = (QLA_real(*b) * h) * h;
    QLA_imag(v) = -(QLA_imag(*b) * h) * h;
    cRarg.a = &v;
    cRarg.b = QDP_expose_R(a->ptr);
    QDP_C_eq_funci(r->ptr, docRmul, *S->qss);
    cRarg.a = 0;
    cRarg.b = 0;
    QDP_reset_R(a->ptr);

    return 1;
}

static int
q_C_project(lua_State *L)
{
    mLatComplex *a = qlua_checkLatComplex(L, 1, NULL);
    mLattice *S = qlua_ObjLattice(L, 1);
    mLatComplex *b = qlua_checkLatComplex(L, 2, S);
    mLatMulti *m = qlua_checkLatMulti(L, 3, S);
    int size = m->size;
    mVecComplex *r = qlua_newVecComplex(L, size);
    int sites = QDP_sites_on_node_L(S->lat);
    QLA_Real rr[2 * size];
    QLA_Int *ii = m->idx;
    QLA_Complex *aa;
    QLA_Complex *bb;
    int k;
    
    for (k = 0; k < 2 * size; k++)
        rr[k] = 0;
    CALL_QDP(L);
    aa = QDP_expose_C(a->ptr);
    bb = QDP_expose_C(b->ptr);
    for (k = 0; k < sites; k++, ii++, aa++, bb++) {
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
    mLatComplex *a = qlua_checkLatComplex(L, 1, NULL);
    int argc = lua_gettop(L);
    mLattice *S = qlua_ObjLattice(L, 1);
    int Sidx = lua_gettop(L);

    switch (argc) {
    case 1: {
        QLA_Complex *s = qlua_newComplex(L);

        CALL_QDP(L);
        if (S->lss.mask) {
            mLatComplex *b = qlua_newLatComplex(L, Sidx);
            QDP_C_eq_zero(b->ptr, *S->qss);
            QDP_C_eq_C_mask_I(b->ptr, a->ptr, S->lss.mask, *S->qss);
            QDP_c_eq_sum_C(s, b->ptr, *S->qss);
            lua_pop(L, 1);
        } else {
            QDP_c_eq_sum_C(s, a->ptr, *S->qss);
        }

        return 1;
    }
    case 2: {
        mLatMulti *m = qlua_checkLatMulti(L, 2, S);
        int size = m->size;
        QLA_Int *ii = m->idx;
        mVecComplex *r = qlua_newVecComplex(L, size);
        int sites = QDP_sites_on_node_L(S->lat);
        int k;
        QLA_Complex *xx;
        QLA_Real rr[2 * size];
        
        for (k = 0; k < 2 * size; k++)
            rr[k] = 0;

        CALL_QDP(L);
        xx = QDP_expose_C(a->ptr);
        for (k = 0; k < sites; k++, xx++, ii++) {
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

        return 1;
    }
    }
    return luaL_error(L, "bad arguments for Complex:sum()");
}

static int
q_C_norm2(lua_State *L)
{
    mLatComplex *a = qlua_checkLatComplex(L, 1, NULL);
    mLattice *S = qlua_ObjLattice(L, 1);
    QLA_Real n;

    CALL_QDP(L);
    if (S->lss.mask) {
        mLatComplex *b = qlua_newLatComplex(L, lua_gettop(L));
        QDP_C_eq_C_mask_I(b->ptr, a->ptr, S->lss.mask, *S->qss);
        QDP_r_eq_norm2_C(&n, b->ptr, *S->qss);
    } else {
        QDP_r_eq_norm2_C(&n, a->ptr, *S->qss);
    }
    lua_pushnumber(L, n);
    
    return 1;
}

static int
q_C_shift(lua_State *L)
{
    mLatComplex *a = qlua_checkLatComplex(L, 1, NULL);
    mLattice *S = qlua_ObjLattice(L, 1);
    QDP_Shift shift = qlua_checkShift(L, 2, S);
    QDP_ShiftDir dir = qlua_checkShiftDir(L, 3);
    mLatComplex *r = qlua_newLatComplex(L, lua_gettop(L));

    CALL_QDP(L);
    QDP_C_eq_sC(r->ptr, a->ptr, shift, dir, *S->qss);

    return 1;
}

static int
q_C_conj(lua_State *L)
{
    mLatComplex *a = qlua_checkLatComplex(L, 1, NULL);
    mLattice *S = qlua_ObjLattice(L, 1);
    mLatComplex *r = qlua_newLatComplex(L, lua_gettop(L));

    CALL_QDP(L);
    QDP_C_eq_conj_C(r->ptr, a->ptr, *S->qss);

    return 1;
}

static int
q_C_abs(lua_State *L)
{
    mLatComplex *a = qlua_checkLatComplex(L, 1, NULL);
    mLattice *S = qlua_ObjLattice(L, 1);
    mLatReal *r = qlua_newLatReal(L, lua_gettop(L));

    CALL_QDP(L);
    QDP_R_eq_norm_C(r->ptr, a->ptr, *S->qss);

    return 1;
}

static int
q_C_arg(lua_State *L)
{
    mLatComplex *a = qlua_checkLatComplex(L, 1, NULL);
    mLattice *S = qlua_ObjLattice(L, 1);
    mLatReal *r = qlua_newLatReal(L, lua_gettop(L));

    CALL_QDP(L);
    QDP_R_eq_arg_C(r->ptr, a->ptr, *S->qss);

    return 1;
}

static int
q_C_sqrt(lua_State *L)
{
    mLatComplex *a = qlua_checkLatComplex(L, 1, NULL);
    mLattice *S = qlua_ObjLattice(L, 1);
    mLatComplex *r = qlua_newLatComplex(L, lua_gettop(L));

    CALL_QDP(L);
    QDP_C_eq_csqrt_C(r->ptr, a->ptr, *S->qss);

    return 1;
}

static int
q_C_exp(lua_State *L)
{
    mLatComplex *a = qlua_checkLatComplex(L, 1, NULL);
    mLattice *S = qlua_ObjLattice(L, 1);
    mLatComplex *r = qlua_newLatComplex(L, lua_gettop(L));

    CALL_QDP(L);
    QDP_C_eq_cexp_C(r->ptr, a->ptr, *S->qss);

    return 1;
}

static int
q_C_log(lua_State *L)
{
    mLatComplex *a = qlua_checkLatComplex(L, 1, NULL);
    mLattice *S = qlua_ObjLattice(L, 1);
    mLatComplex *r = qlua_newLatComplex(L, lua_gettop(L));

    CALL_QDP(L);
    QDP_C_eq_clog_C(r->ptr, a->ptr, *S->qss);

    return 1;
}

static int
q_C_set(lua_State *L)
{
    mLatComplex *r = qlua_checkLatComplex(L, 1, NULL);
    mLattice *S = qlua_ObjLattice(L, 1);
    mLatComplex *a = qlua_checkLatComplex(L, 2, S);

    CALL_QDP(L);
    if (S->lss.mask)
        QDP_C_eq_C_mask_I(r->ptr, a->ptr, S->lss.mask, *S->qss);
    else
        QDP_C_eq_C(r->ptr, a->ptr, *S->qss);
    lua_pop(L, 1);

    return 1;
}

int
q_C_gaussian(lua_State *L)
{
    mLatRandom *a = qlua_checkLatRandom(L, 1, NULL);
    mLattice *S = qlua_ObjLattice(L, 1);
    mLatComplex *r = qlua_newLatComplex(L, lua_gettop(L));

    CALL_QDP(L);
    QDP_C_eq_gaussian_S(r->ptr, a->ptr, *S->qss);

    return 1;
}

static int
q_C_dot(lua_State *L)
{
    mLatComplex *a = qlua_checkLatComplex(L, 1, NULL);
    mLattice *S = qlua_ObjLattice(L, 1);
    mLatComplex *b = qlua_checkLatComplex(L, 2, S);
    mLatComplex *s = qlua_newLatComplex(L, lua_gettop(L));

    CALL_QDP(L);
    QDP_C_eq_C_dot_C(s->ptr, a->ptr, b->ptr, *S->qss);

    return 1;
}

static int
q_latcomplex(lua_State *L)
{
    mLattice *S = qlua_checkLattice(L, 1);

    switch (lua_gettop(L)) {
    case 1: {
        mLatComplex *v = qlua_newLatComplex(L, 1);
        
        CALL_QDP(L);
        QDP_C_eq_zero(v->ptr, *S->qss);

        return 1;
    }
    case 2:
        switch (qlua_qtype(L, 2)) {
        case qReal: {
            QLA_Real d = luaL_checknumber(L, 2);
            QLA_Complex z;
            mLatComplex *v = qlua_newLatComplex(L, 1);
            
            CALL_QDP(L);
            QLA_real(z) = d;
            QLA_imag(z) = 0.0;
            QDP_C_eq_c(v->ptr, &z, *S->qss);

            return 1;
        }
        case qComplex: {
            QLA_Complex *z = qlua_checkComplex(L, 2);
            mLatComplex *v = qlua_newLatComplex(L, 1);

            CALL_QDP(L);
            QDP_C_eq_c(v->ptr, z, *S->qss);
            
            return 1;
        }
        case qLatReal: {
            mLatReal *d = qlua_checkLatReal(L, 2, S);
            mLatComplex *v = qlua_newLatComplex(L, 1);
            
            CALL_QDP(L);
            QDP_C_eq_R(v->ptr, d->ptr, *S->qss);

            return 1;
        }
        case qLatComplex: {
            mLatComplex *d = qlua_checkLatComplex(L, 2, S);
            mLatComplex *v = qlua_newLatComplex(L, 1);
            
            CALL_QDP(L);
            QDP_C_eq_C(v->ptr, d->ptr, *S->qss);

            return 1;
        }
        default:
            break;
        }
    case 3: {
        mLatReal *a = qlua_checkLatReal(L, 2, S);
        mLatReal *b = qlua_checkLatReal(L, 3, S);
        mLatComplex *c = qlua_newLatComplex(L, 1);

        CALL_QDP(L);
        QDP_C_eq_R_plus_i_R(c->ptr, a->ptr, b->ptr, *S->qss);

        return 1;
    }
    }
    return qlua_badconstr(L, "Complex");
}

static struct luaL_Reg mtLatComplex[] = {
    { "__tostring",        q_C_fmt      },
    { "__gc",              q_C_gc       },
    { "__index",           q_C_get      },
    { "__newindex",        q_C_put      },
    { "__unm",             q_C_neg      },
    { "__add",             qlua_add     },
    { "__sub",             qlua_sub     },
    { "__mul",             qlua_mul     },
    { "__div",             qlua_div     },
    { "real",              q_C_real     },
    { "imag",              q_C_imag     },
    { "sum",               q_C_sum      },
    { "norm2",             q_C_norm2    },
    { "shift",             q_C_shift    },
    { "conj",              q_C_conj     },
    { "abs",               q_C_abs      },
    { "arg",               q_C_arg      },
    { "sqrt",              q_C_sqrt     },
    { "exp",               q_C_exp      },
    { "log",               q_C_log      },
    { "set",               q_C_set      },
    { "project",           q_C_project  },
    /* "lattice" */
    /* "a-type" */
    { NULL,                NULL         }
};

mLatComplex *
qlua_newLatComplex(lua_State *L, int Sidx)
{
    mLattice *S = qlua_checkLattice(L, Sidx);
    QDP_Complex *v = QDP_create_C_L(S->lat);
    mLatComplex *hdr;

    if (v == 0) {
        lua_gc(L, LUA_GCCOLLECT, 0);
        v = QDP_create_C_L(S->lat);
        if (v == 0)
            luaL_error(L, "not enough memory (QDP_Complex)");
    }
    hdr = lua_newuserdata(L, sizeof (mLatComplex));
    hdr->ptr = v;
    qlua_createLatticeTable(L, Sidx, mtLatComplex, qLatComplex, LatComplexName);
    lua_setmetatable(L, -2);

    return hdr;
}

mLatComplex *
qlua_checkLatComplex(lua_State *L, int idx, mLattice *S)
{
    void *v = qlua_checkLatticeType(L, idx, qLatComplex, LatComplexName);
    
    if (S) {
        mLattice *S1 = qlua_ObjLattice(L, idx);
        if (S1->id != S->id)
            luaL_error(L, "%s on a wrong lattice", LatComplexName);
        lua_pop(L, 1);
    }

    return (mLatComplex *)v;
}

static struct luaL_Reg fLatComplex[] = {
    { "Complex",     q_latcomplex  },
    { NULL,          NULL }
};

int
init_latcomplex(lua_State *L)
{
    static const QLUA_Op2 ops[] = {
        { qlua_add_table, qLatComplex,  qLatComplex,  q_C_add_C },
        { qlua_add_table, qLatComplex,  qLatReal,     q_C_add_R },
        { qlua_add_table, qLatReal,     qLatComplex,  q_R_add_C },
        { qlua_add_table, qLatComplex,  qReal,        q_C_add_r },
        { qlua_add_table, qReal,        qLatComplex,  q_r_add_C },
        { qlua_add_table, qLatComplex,  qComplex,     q_C_add_c },
        { qlua_add_table, qComplex,     qLatComplex,  q_c_add_C },
        { qlua_add_table, qComplex,     qLatReal,     q_c_add_R },
        { qlua_add_table, qLatReal,     qComplex,     q_R_add_c },
        { qlua_sub_table, qLatComplex,  qLatComplex,  q_C_sub_C },
        { qlua_sub_table, qLatComplex,  qLatReal,     q_C_sub_R },
        { qlua_sub_table, qLatReal,     qLatComplex,  q_R_sub_C },
        { qlua_sub_table, qLatComplex,  qReal,        q_C_sub_r },
        { qlua_sub_table, qReal,        qLatComplex,  q_r_sub_C },
        { qlua_sub_table, qLatComplex,  qComplex,     q_C_sub_c },
        { qlua_sub_table, qComplex,     qLatComplex,  q_c_sub_C },
        { qlua_sub_table, qComplex,     qLatReal,     q_c_sub_R },
        { qlua_sub_table, qLatReal,     qComplex,     q_R_sub_c },
        { qlua_mul_table, qLatComplex,  qLatComplex,  q_C_mul_C },
        { qlua_mul_table, qLatComplex,  qLatReal,     q_C_mul_R },
        { qlua_mul_table, qLatReal,     qLatComplex,  q_R_mul_C },
        { qlua_mul_table, qLatComplex,  qReal,        q_C_mul_r },
        { qlua_mul_table, qReal,        qLatComplex,  q_r_mul_C },
        { qlua_mul_table, qLatComplex,  qComplex,     q_C_mul_c },
        { qlua_mul_table, qComplex,     qLatComplex,  q_c_mul_C },
        { qlua_mul_table, qComplex,     qLatReal,     q_c_mul_R },
        { qlua_mul_table, qLatReal,     qComplex,     q_R_mul_c },
        { qlua_div_table, qLatComplex,  qLatComplex,  q_C_div_C },
        { qlua_div_table, qLatComplex,  qLatReal,     q_C_div_R },
        { qlua_div_table, qLatReal,     qLatComplex,  q_R_div_C },
        { qlua_div_table, qLatComplex,  qReal,        q_C_div_r },
        { qlua_div_table, qReal,        qLatComplex,  q_r_div_C },
        { qlua_div_table, qLatComplex,  qComplex,     q_C_div_c },
        { qlua_div_table, qComplex,     qLatComplex,  q_c_div_C },
        { qlua_div_table, qComplex,     qLatReal,     q_c_div_R },
        { qlua_div_table, qLatReal,     qComplex,     q_R_div_c },
        { NULL,           qNoType,      qNoType,      NULL }
    };

    luaL_getmetatable(L, opLattice);
    luaL_register(L, NULL, fLatComplex);
    lua_pop(L, 1);
    qlua_reg_op2(ops);

    qlua_reg_dot(qLatComplex, q_C_dot);

    return 0;
}

int
fini_latcomplex(lua_State *L)
{
    return 0;
}
