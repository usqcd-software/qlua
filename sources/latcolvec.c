#include <qlua.h>                                                    /* DEPS */
#include <lattice.h>                                                 /* DEPS */
#include <latcolvec.h>                                               /* DEPS */
#include <qcomplex.h>                                                /* DEPS */
#include <latint.h>                                                  /* DEPS */
#include <latcomplex.h>                                              /* DEPS */
#include <latrandom.h>                                               /* DEPS */
#include <qmp.h>

const char mtnLatColVec[] = "qcd.lattice.colvec";
static const char opLatColVec[] = "qcd.lattice.colvec.op";

mLatColVec *
qlua_newLatColVec(lua_State *L)
{
    QDP_ColorVector *v = QDP_create_V();
    mLatColVec *hdr;

    if (v == 0) {
        lua_gc(L, LUA_GCCOLLECT, 0);
        v = QDP_create_V();
        if (v == 0)
            luaL_error(L, "not enough memory (QDP_ColorVector)");
    }
    hdr = lua_newuserdata(L, sizeof (mLatColVec));
    hdr->ptr = v;
    luaL_getmetatable(L, mtnLatColVec);
    lua_setmetatable(L, -2);

    return hdr;
}

mLatColVec *
qlua_checkLatColVec(lua_State *L, int idx)
{
    void *v = luaL_checkudata(L, idx, mtnLatColVec);
    
    luaL_argcheck(L, v != 0, idx, "qcd.ColorVector expected");
    
    return v;
}

static int
q_V_fmt(lua_State *L)
{
    char fmt[72];
    mLatColVec *b = qlua_checkLatColVec(L, 1);

    sprintf(fmt, "QDP:ColorVector(%p)", b->ptr);
    lua_pushstring(L, fmt);

    return 1;
}

static int
q_V_gc(lua_State *L)
{
    mLatColVec *b = qlua_checkLatColVec(L, 1);

    QDP_destroy_V(b->ptr);
    b->ptr = 0;

    return 0;
}

static int
q_V_get(lua_State *L)
{
    switch (qlua_gettype(L, 2)) {
    case qTable: {
        mLatColVec *V = qlua_checkLatColVec(L, 1);
        int c = qlua_checkcolorindex(L, 2);
        int *idx = qlua_latcoord(L, 2);

        if (idx == NULL) {
            mLatComplex *r = qlua_newLatComplex(L);
            
            QDP_C_eq_elem_V(r->ptr, V->ptr, c, qCurrent);
        } else {
            QLA_Complex *W = qlua_newComplex(L);
            QLA_ColorVector *locked;
            double z_re, z_im;

            locked = QDP_expose_V(V->ptr);
            if (QDP_node_number(idx) == QDP_this_node) {
                QLA_Complex *zz = &QLA_elem_V(locked[QDP_index(idx)], c);
                z_re = QLA_real(*zz);
                z_im = QLA_imag(*zz);
            } else {
                z_re = 0;
                z_im = 0;
            }
            QDP_reset_V(V->ptr);
            QMP_sum_double(&z_re);
            QMP_sum_double(&z_im);
            QLA_real(*W) = z_re;
            QLA_imag(*W) = z_im;
            if (QLA_imag(*W) == 0)
                lua_pushnumber(L, QLA_real(*W));
        }
        qlua_free(L, idx);
        return 1;
    }
    case qString:
        return qlua_lookup(L, 2, opLatColVec);
    }
    return qlua_badindex(L, "ColorVector");
}

static int
q_V_put(lua_State *L)
{
    mLatColVec *V = qlua_checkLatColVec(L, 1);
    int c = qlua_checkcolorindex(L, 2);
    int *idx = qlua_latcoord(L, 2);

    if (idx == NULL) {
        mLatComplex *z = qlua_checkLatComplex(L, 3);
        
        QDP_V_eq_elem_C(V->ptr, z->ptr, c, qCurrent);
    } else {
        QLA_Complex *z = qlua_checkComplex(L, 3);
        if (QDP_node_number(idx) == QDP_this_node) {
            QLA_ColorVector *locked = QDP_expose_V(V->ptr);
            QLA_Complex *zz = &QLA_elem_V(locked[QDP_index(idx)], c);
            QLA_c_eq_c(*z, *zz);
            QDP_reset_V(V->ptr);
        }
    }
    qlua_free(L, idx);
    return 0;
}

static int
q_V_dot(lua_State *L)
{
    mLatColVec *a = qlua_checkLatColVec(L, 1);
    mLatColVec *b = qlua_checkLatColVec(L, 2);
    mLatComplex *s = qlua_newLatComplex(L);

    QDP_C_eq_V_dot_V(s->ptr, a->ptr, b->ptr, qCurrent);

    return 1;
}

int
q_V_gaussian(lua_State *L)
{
    mLatRandom *a = qlua_checkLatRandom(L, 1);
    mLatColVec *r = qlua_newLatColVec(L);

    QDP_V_eq_gaussian_S(r->ptr, a->ptr, qCurrent);

    return 1;
}


static int
q_V_add_V(lua_State *L)
{
    mLatColVec *a = qlua_checkLatColVec(L, 1);
    mLatColVec *b = qlua_checkLatColVec(L, 2);
    mLatColVec *c = qlua_newLatColVec(L);

    QDP_V_eq_V_plus_V(c->ptr, a->ptr, b->ptr, qCurrent);

    return 1;
}

static int
q_V_sub_V(lua_State *L)
{
    mLatColVec *a = qlua_checkLatColVec(L, 1);
    mLatColVec *b = qlua_checkLatColVec(L, 2);
    mLatColVec *c = qlua_newLatColVec(L);

    QDP_V_eq_V_minus_V(c->ptr, a->ptr, b->ptr, qCurrent);

    return 1;
}

static int
q_V_mul_r(lua_State *L)
{
    mLatColVec *a = qlua_checkLatColVec(L, 1);
    QLA_Real b = luaL_checknumber(L, 2);
    mLatColVec *c = qlua_newLatColVec(L);

    QDP_V_eq_r_times_V(c->ptr, &b, a->ptr, qCurrent);

    return 1;
}

static int
q_r_mul_V(lua_State *L)
{
    QLA_Real a = luaL_checknumber(L, 1);
    mLatColVec *b = qlua_checkLatColVec(L, 2);
    mLatColVec *c = qlua_newLatColVec(L);

    QDP_V_eq_r_times_V(c->ptr, &a, b->ptr, qCurrent);

    return 1;
}

static int
q_V_mul_c(lua_State *L)
{
    mLatColVec *a = qlua_checkLatColVec(L, 1);
    QLA_Complex *b = qlua_checkComplex(L, 2);
    mLatColVec *c = qlua_newLatColVec(L);

    QDP_V_eq_c_times_V(c->ptr, b, a->ptr, qCurrent);

    return 1;
}

static int
q_c_mul_V(lua_State *L)
{
    QLA_Complex *a = qlua_checkComplex(L, 1);
    mLatColVec *b = qlua_checkLatColVec(L, 2);
    mLatColVec *c = qlua_newLatColVec(L);

    QDP_V_eq_c_times_V(c->ptr, a, b->ptr, qCurrent);

    return 1;
}

static int
q_V_norm2(lua_State *L)
{
    mLatColVec *a = qlua_checkLatColVec(L, 1);
    QLA_Real n;

    QDP_r_eq_norm2_V(&n, a->ptr, qCurrent);
    lua_pushnumber(L, n);
    
    return 1;
}

static int
q_V_shift(lua_State *L)
{
    mLatColVec *a = qlua_checkLatColVec(L, 1);
    QDP_Shift shift = qlua_checkShift(L, 2);
    QDP_ShiftDir dir = qlua_checkShiftDir(L, 3);
    mLatColVec *r = qlua_newLatColVec(L);

    QDP_V_eq_sV(r->ptr, a->ptr, shift, dir, qCurrent);

    return 1;
}

static int
q_V_conj(lua_State *L)
{
    mLatColVec *a = qlua_checkLatColVec(L, 1);
    mLatColVec *r = qlua_newLatColVec(L);

    QDP_V_eq_conj_V(r->ptr, a->ptr, qCurrent);

    return 1;
}

static int
q_V_set(lua_State *L)
{
    mLatColVec *r = qlua_checkLatColVec(L, 1);
    mLatColVec *a = qlua_checkLatColVec(L, 2);

    QDP_V_eq_V(r->ptr, a->ptr, qCurrent);
    lua_pop(L, 1);

    return 1;
}

static int
q_V_neg(lua_State *L)
{
    mLatColVec *a = qlua_checkLatColVec(L, 1);
    mLatColVec *r = qlua_newLatColVec(L);
    QLA_Real m1 = -1;

    QDP_V_eq_r_times_V(r->ptr, &m1, a->ptr, qCurrent);

    return 1;
}
static int
q_V_div_r(lua_State *L)
{
    mLatColVec *a = qlua_checkLatColVec(L, 1);
    QLA_Real b = 1 / luaL_checknumber(L, 2);
    mLatColVec *c = qlua_newLatColVec(L);

    QDP_V_eq_r_times_V(c->ptr, &b, a->ptr, qCurrent);

    return 1;
}

static int
q_V_div_c(lua_State *L)
{
    mLatColVec *a = qlua_checkLatColVec(L, 1);
    QLA_Complex *b = qlua_checkComplex(L, 2);
    mLatColVec *c = qlua_newLatColVec(L);
    double n = 1 / (QLA_real(*b) * QLA_real(*b) + QLA_imag(*b) * QLA_imag(*b));
    QLA_Complex s;

    QLA_real(s) = n * QLA_real(*b);
    QLA_imag(s) = -n * QLA_imag(*b);
    QDP_V_eq_c_times_V(c->ptr, &s, a->ptr, qCurrent);

    return 1;
}


static int
q_latcolvec(lua_State *L)
{
    switch (lua_gettop(L)) {
    case 1: {
        mLatColVec *v = qlua_newLatColVec(L);

        QDP_V_eq_zero(v->ptr, qCurrent);

        return 1;
    }
    case 2: {
        mLatColVec *a = qlua_checkLatColVec(L, 2);
        mLatColVec *r = qlua_newLatColVec(L);
        
        QDP_V_eq_V(r->ptr, a->ptr, qCurrent);
        
        return 1;
    }
    case 3: {
        mLatComplex *c = qlua_checkLatComplex(L, 2);
        int a = luaL_checkint(L, 3);
        mLatColVec *r = qlua_newLatColVec(L);

        QDP_V_eq_elem_C(r->ptr, c->ptr, a, qCurrent);

        return 1;
    }
    }
    return qlua_badconstr(L, "ColorVector");
}

static struct luaL_Reg LatColVecMethods[] = {
    { "norm2",      q_V_norm2 },
    { "shift",      q_V_shift },
    { "conj",       q_V_conj },
    { "set",        q_V_set },
    { NULL,         NULL }
};

static struct luaL_Reg mtLatColVec[] = {
    { "__tostring",        q_V_fmt },
    { "__gc",              q_V_gc },
    { "__index",           q_V_get },
    { "__newindex",        q_V_put },
    { "__unm",             q_V_neg },
    { "__add",             qlua_add },
    { "__sub",             qlua_sub },
    { "__mul",             qlua_mul },
    { "__div",             qlua_div },
    { NULL,                NULL }
};

static struct luaL_Reg fLatColVec[] = {
    { "ColorVector",       q_latcolvec },
    { NULL,                NULL }
};

int
init_latcolvec(lua_State *L)
{
    luaL_getmetatable(L, opLattice);
    luaL_register(L, NULL, fLatColVec);
    lua_pop(L, 1);
    qlua_metatable(L, mtnLatColVec, mtLatColVec);
    qlua_metatable(L, opLatColVec, LatColVecMethods);
    qlua_reg_add(qLatColVec, qLatColVec, q_V_add_V);
    qlua_reg_sub(qLatColVec, qLatColVec, q_V_sub_V);
    qlua_reg_mul(qReal,      qLatColVec, q_r_mul_V);
    qlua_reg_mul(qLatColVec, qReal,      q_V_mul_r);
    qlua_reg_mul(qComplex,   qLatColVec, q_c_mul_V);
    qlua_reg_mul(qLatColVec, qComplex,   q_V_mul_c);
    qlua_reg_div(qLatColVec, qReal,      q_V_div_r);
    qlua_reg_div(qLatColVec, qComplex,   q_V_div_c);
    qlua_reg_dot(qLatColVec, q_V_dot);

    return 0;
}

int
fini_latcolvec(lua_State *L)
{
    return 0;
}
