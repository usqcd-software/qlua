#include <qlua.h>
#include <latcolmat.h>
#include <qcomplex.h>
#include <latint.h>
#include <latrandom.h>
#include <latcomplex.h>
#include <latcolvec.h>
#include <qmp.h>

const char *mtnLatColMat = "qcd.lattice.ColorMatrix";
static char *opLatColMat = "qcd.lattice.ColorMatrix.op";

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
    
    luaL_argcheck(L, v != 0, idx, "qcd.ColorMatrix expected");
    
    return v;
}

int
qLatColMat_fmt(lua_State *L)
{
    char fmt[72];
    mLatColMat *b = qlua_checkLatColMat(L, 1);

    sprintf(fmt, "LatColMat(%p)", b->ptr);
    lua_pushstring(L, fmt);

    return 1;
}

int
qLatColMat_gc(lua_State *L)
{
    mLatColMat *b = qlua_checkLatColMat(L, 1);

    QDP_destroy_M(b->ptr);
    b->ptr = 0;

    return 0;
}

int
qLatColMat_get(lua_State *L)
{
    switch (qlua_gettype(L, 2)) {
    case qTable: {
        mLatColMat *V = qlua_checkLatColMat(L, 1);
        int b = qlua_checkrightindex(L, 2, QDP_Nc);
        int a = qlua_leftindex(L, 2, QDP_Nc);
        int *idx = qlua_latcoord(L, 2);
        if (idx == NULL) {
            if (a == -1) {
                mLatColVec *r = qlua_newLatColVec(L);
                
                QDP_V_eq_colorvec_M(r->ptr, V->ptr, b, QDP_all);
            } else {
                mLatComplex *r = qlua_newLatComplex(L);
                
                QDP_C_eq_elem_M(r->ptr, V->ptr, a, b, QDP_all);
            }
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
        qlua_free(L, idx);
        return 1;
    }
    case qString:
        return qlua_lookup(L, 2, opLatColMat);
    }
    return luaL_error(L, "bad index");
}

int
qLatColMat_put(lua_State *L)
{
    mLatColMat *V = qlua_checkLatColMat(L, 1);
    int b = qlua_checkrightindex(L, 2, QDP_Nc);
    int a = qlua_leftindex(L, 2, QDP_Nc);
    int *idx = qlua_latcoord(L, 2);

    if (idx == NULL) {
        if (a == -1) {
            mLatColVec *z = qlua_checkLatColVec(L, 3);

            QDP_M_eq_colorvec_V(V->ptr, z->ptr, b, QDP_all);
        } else {
            mLatComplex *z = qlua_checkLatComplex(L, 3);
            
            QDP_M_eq_elem_C(V->ptr, z->ptr, a, b, QDP_all);
        }
    } else {
        QLA_Complex *z = qlua_checkComplex(L, 3);
        if (QDP_node_number(idx) == QDP_this_node) {
            QLA_ColorMatrix *locked = QDP_expose_M(V->ptr);
            QLA_Complex *zz = &QLA_elem_M(locked[QDP_index(idx)], a, b);
            QLA_c_eq_c(*z, *zz);
            QDP_reset_M(V->ptr);
        }
    }
    qlua_free(L, idx);
    return 0;
}

/* index
 *  {n0, ..., nd, a=n, b=m} -- complex value at the index
 *  {a=n, b=m} -- Complex value at the index
 *  {b=n} -- ColorVector value at the index
 */

int
q_M_norm2(lua_State *L)
{
    mLatColMat *a = qlua_checkLatColMat(L, 1);
    QLA_Real n;

    QDP_r_eq_norm2_M(&n, a->ptr, QDP_all);
    lua_pushnumber(L, n);
    
    return 1;
}

int
q_M_shift(lua_State *L)
{
    mLatColMat *a = qlua_checkLatColMat(L, 1);
    QDP_Shift shift = qlua_checkShift(L, 2);
    QDP_ShiftDir dir = qlua_checkShiftDir(L, 3);
    mLatColMat *r = qlua_newLatColMat(L);

    QDP_M_eq_sM(r->ptr, a->ptr, shift, dir, QDP_all);

    return 1;
}

int
q_M_conj(lua_State *L)
{
    mLatColMat *a = qlua_checkLatColMat(L, 1);
    mLatColMat *r = qlua_newLatColMat(L);

    QDP_M_eq_conj_M(r->ptr, a->ptr, QDP_all);

    return 1;
}

int
q_M_trans(lua_State *L)
{
    mLatColMat *a = qlua_checkLatColMat(L, 1);
    mLatColMat *r = qlua_newLatColMat(L);

    QDP_M_eq_transpose_M(r->ptr, a->ptr, QDP_all);

    return 1;
}

int
q_M_adjoin(lua_State *L)
{
    mLatColMat *a = qlua_checkLatColMat(L, 1);
    mLatColMat *r = qlua_newLatColMat(L);

    QDP_M_eq_Ma(r->ptr, a->ptr, QDP_all);

    return 1;
}

int
q_M_trace(lua_State *L)
{
    mLatColMat *a = qlua_checkLatColMat(L, 1);
    mLatComplex *r = qlua_newLatComplex(L);

    QDP_C_eq_trace_M(r->ptr, a->ptr, QDP_all);

    return 1;
}

int
q_M_dot(lua_State *L)
{
    mLatColMat *a = qlua_checkLatColMat(L, 1);
    mLatColMat *b = qlua_checkLatColMat(L, 2);
    mLatComplex *s = qlua_newLatComplex(L);

    QDP_C_eq_M_dot_M(s->ptr, a->ptr, b->ptr, QDP_all);

    return 1;
}

int
q_M_gaussian(lua_State *L)
{
    mLatRandom *a = qlua_checkLatRandom(L, 1);
    mLatColMat *r = qlua_newLatColMat(L);

    QDP_M_eq_gaussian_S(r->ptr, a->ptr, QDP_all);

    return 1;
}


int
q_M_add_M(lua_State *L)
{
    mLatColMat *a = qlua_checkLatColMat(L, 1);
    mLatColMat *b = qlua_checkLatColMat(L, 2);
    mLatColMat *c = qlua_newLatColMat(L);

    QDP_M_eq_M_plus_M(c->ptr, a->ptr, b->ptr, QDP_all);

    return 1;
}

int
q_M_sub_M(lua_State *L)
{
    mLatColMat *a = qlua_checkLatColMat(L, 1);
    mLatColMat *b = qlua_checkLatColMat(L, 2);
    mLatColMat *c = qlua_newLatColMat(L);

    QDP_M_eq_M_minus_M(c->ptr, a->ptr, b->ptr, QDP_all);

    return 1;
}

int
q_M_mul_M(lua_State *L)
{
    mLatColMat *a = qlua_checkLatColMat(L, 1);
    mLatColMat *b = qlua_checkLatColMat(L, 2);
    mLatColMat *c = qlua_newLatColMat(L);

    QDP_M_eq_M_times_M(c->ptr, a->ptr, b->ptr, QDP_all);

    return 1;
}

int
q_M_mul_V(lua_State *L)
{
    mLatColMat *a = qlua_checkLatColMat(L, 1);
    mLatColVec *b = qlua_checkLatColVec(L, 2);
    mLatColVec *c = qlua_newLatColVec(L);

    QDP_V_eq_M_times_V(c->ptr, a->ptr, b->ptr, QDP_all);

    return 1;
}

int
q_M_mul_r(lua_State *L)
{
    QLA_Real a = luaL_checknumber(L, 1);
    mLatColMat *b = qlua_checkLatColMat(L, 2);
    mLatColMat *c = qlua_newLatColMat(L);

    QDP_M_eq_r_times_M(c->ptr, &a, b->ptr, QDP_all);

    return 1;
}

int
q_r_mul_M(lua_State *L)
{
    mLatColMat *a = qlua_checkLatColMat(L, 1);
    QLA_Real b = luaL_checknumber(L, 2);
    mLatColMat *c = qlua_newLatColMat(L);

    QDP_M_eq_r_times_M(c->ptr, &b, a->ptr, QDP_all);

    return 1;
}

int
q_M_mul_c(lua_State *L)
{
    QLA_Complex *a = qlua_checkComplex(L, 1);
    mLatColMat *b = qlua_checkLatColMat(L, 2);
    mLatColMat *c = qlua_newLatColMat(L);

    QDP_M_eq_c_times_M(c->ptr, a, b->ptr, QDP_all);

    return 1;
}

int
q_c_mul_M(lua_State *L)
{
    mLatColMat *a = qlua_checkLatColMat(L, 1);
    QLA_Complex *b = qlua_checkComplex(L, 2);
    mLatColMat *c = qlua_newLatColMat(L);

    QDP_M_eq_c_times_M(c->ptr, b, a->ptr, QDP_all);

    return 1;
}

int
q_neg_M(lua_State *L)
{
    mLatColMat *a = qlua_checkLatColMat(L, 1);
    mLatColMat *r = qlua_newLatColMat(L);
    QLA_Real m1 = -1;

    QDP_M_eq_r_times_M(r->ptr, &m1, a->ptr, QDP_all);

    return 1;
}

int
q_latcolmat(lua_State *L)
{
    switch (lua_gettop(L)) {
    case 0: {
        mLatColMat *v = qlua_newLatColMat(L);

        QDP_M_eq_zero(v->ptr, QDP_all);
        
        return 1;
    }
    case 1: {
        switch (qlua_gettype(L, 1)) {
        case qReal: {
            QLA_Real x = luaL_checknumber(L, 1);
            mLatColMat *v = qlua_newLatColMat(L);
            QLA_Complex z;

            QLA_real(z) = x;
            QLA_imag(z) = 0;
            QDP_M_eq_c(v->ptr, &z, QDP_all);

            return 1;
        }
        case qComplex: {
            QLA_Complex *z = qlua_checkComplex(L, 1);
            mLatColMat *v = qlua_newLatColMat(L);

            QDP_M_eq_c(v->ptr, z, QDP_all);

            return 1;
        }
        case qLatColMat: {
            mLatColMat *w = qlua_checkLatColMat(L, 1);
            mLatColMat *v = qlua_newLatColMat(L);

            QDP_M_eq_M(v->ptr, w->ptr, QDP_all);

            return 1;
        }
        }
        break;
    }
    case 2: {
        mLatColVec *x = qlua_checkLatColVec(L, 1);

        switch (qlua_gettype(L, 2)) {
        case qReal: {
            int b = luaL_checkint(L, 2);
            if ((b >= 0) && (b < QDP_Nc)) {
                mLatColMat *v = qlua_newLatColMat(L);
                
                QDP_M_eq_zero(v->ptr, QDP_all);
                QDP_M_eq_colorvec_V(v->ptr, x->ptr, b, QDP_all);
                
                return 1;
            }
        }
        case qLatColVec: {
            mLatColVec *y = qlua_checkLatColVec(L, 2);
            mLatColMat *v = qlua_newLatColMat(L);

            QDP_M_eq_V_times_Va(v->ptr, x->ptr, y->ptr, QDP_all);

            return 1;
        }
        }
        break;
    }
    case 3: {
        mLatComplex *z = qlua_checkLatComplex(L, 1);
        int a = luaL_checkint(L, 2);
        int b = luaL_checkint(L, 3);

        if ((a >= 0) && (a < QDP_Nc) && (b >= 0) && (b < QDP_Nc)) {
            mLatColMat *v = qlua_newLatColMat(L);

            QDP_M_eq_zero(v->ptr, QDP_all);
            QDP_M_eq_elem_C(v->ptr, z->ptr, a, b, QDP_all);

            return 1;
        }
        break;
    }
    }
    return luaL_error(L, "bad arguments");
}
/*
   0 - zero
   1 r -- M1 * r
     c -- M1 * c
     M -- copy
   2 V V -- V * Va
     V j -- elem_V
   3 C i j  elem_C 
 */

static struct luaL_Reg LatColMatMethods[] = {
    { "norm2",      q_M_norm2 },
    { "shift",      q_M_shift },
    { "conj",       q_M_conj },
    { "transpose",  q_M_trans },
    { "adjoin",     q_M_adjoin },
    { "trace",      q_M_trace },
    { NULL,         NULL }
};

static struct luaL_Reg mtLatColMat[] = {
    { "__tostring",        qLatColMat_fmt },
    { "__gc",              qLatColMat_gc },
    { "__index",           qLatColMat_get },
    { "__newindex",        qLatColMat_put },
    { "__umn",             q_neg_M },
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
    luaL_register(L, qcdlib, fLatColMat);
    qlua_metatable(L, mtnLatColMat, mtLatColMat);
    qlua_metatable(L, opLatColMat, LatColMatMethods);

    return 0;
}

int
fini_latcolmat(lua_State *L)
{
    return 0;
}
