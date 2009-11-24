#include <qlua.h>                                                    /* DEPS */
#include <lattice.h>                                                 /* DEPS */
#include <latreal.h>                                                 /* DEPS */
#include <latcolvec.h>                                               /* DEPS */
#include <qcomplex.h>                                                /* DEPS */
#include <latint.h>                                                  /* DEPS */
#include <latcomplex.h>                                              /* DEPS */
#include <latrandom.h>                                               /* DEPS */
#include <qmp.h>

const char mtnLatColVec[] = "lattice.ColorVector";
static const char opLatColVec[] = "lattice.ColorVector.ops";

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
    
    luaL_argcheck(L, v != 0, idx, "lattice.ColorVector expected");
    
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

            CALL_QDP(L);
            QDP_C_eq_elem_V(r->ptr, V->ptr, c, *qCurrent);
        } else {
            QLA_Complex *W = qlua_newComplex(L);
            QLA_ColorVector *locked;
            double zri[2];

            qlua_verifylatcoord(L, idx);
            CALL_QDP(L);
            locked = QDP_expose_V(V->ptr);
            if (QDP_node_number(idx) == QDP_this_node) {
                QLA_Complex *zz = &QLA_elem_V(locked[QDP_index(idx)], c);
                zri[0] = QLA_real(*zz);
                zri[1] = QLA_imag(*zz);
            } else {
                zri[0] = 0;
                zri[1] = 0;
            }
            QDP_reset_V(V->ptr);
            QMP_sum_double_array(zri, 2);
            QLA_c_eq_r_plus_ir(*W, zri[0], zri[1]);
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
        
        CALL_QDP(L);
        QDP_V_eq_elem_C(V->ptr, z->ptr, c, *qCurrent);
    } else {
        QLA_Complex *z = qlua_checkComplex(L, 3);
        qlua_verifylatcoord(L, idx);
        if (QDP_node_number(idx) == QDP_this_node) {
            QLA_ColorVector *locked;
            QLA_Complex *zz;
            CALL_QDP(L);
            locked = QDP_expose_V(V->ptr);
            zz = &QLA_elem_V(locked[QDP_index(idx)], c);
            QLA_c_eq_c(*zz, *z);
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

    CALL_QDP(L);
    QDP_C_eq_V_dot_V(s->ptr, a->ptr, b->ptr, *qCurrent);

    return 1;
}

int
q_V_gaussian(lua_State *L)
{
    mLatRandom *a = qlua_checkLatRandom(L, 1);
    mLatColVec *r = qlua_newLatColVec(L);

    CALL_QDP(L);
    QDP_V_eq_gaussian_S(r->ptr, a->ptr, *qCurrent);

    return 1;
}


static int
q_V_add_V(lua_State *L)
{
    mLatColVec *a = qlua_checkLatColVec(L, 1);
    mLatColVec *b = qlua_checkLatColVec(L, 2);
    mLatColVec *c = qlua_newLatColVec(L);

    CALL_QDP(L);
    QDP_V_eq_V_plus_V(c->ptr, a->ptr, b->ptr, *qCurrent);

    return 1;
}

static int
q_V_sub_V(lua_State *L)
{
    mLatColVec *a = qlua_checkLatColVec(L, 1);
    mLatColVec *b = qlua_checkLatColVec(L, 2);
    mLatColVec *c = qlua_newLatColVec(L);

    CALL_QDP(L);
    QDP_V_eq_V_minus_V(c->ptr, a->ptr, b->ptr, *qCurrent);

    return 1;
}

static int
q_V_mul_r(lua_State *L)
{
    mLatColVec *a = qlua_checkLatColVec(L, 1);
    QLA_Real b = luaL_checknumber(L, 2);
    mLatColVec *c = qlua_newLatColVec(L);

    CALL_QDP(L);
    QDP_V_eq_r_times_V(c->ptr, &b, a->ptr, *qCurrent);

    return 1;
}

static int
q_r_mul_V(lua_State *L)
{
    QLA_Real a = luaL_checknumber(L, 1);
    mLatColVec *b = qlua_checkLatColVec(L, 2);
    mLatColVec *c = qlua_newLatColVec(L);

    CALL_QDP(L);
    QDP_V_eq_r_times_V(c->ptr, &a, b->ptr, *qCurrent);

    return 1;
}

static int
q_V_mul_c(lua_State *L)
{
    mLatColVec *a = qlua_checkLatColVec(L, 1);
    QLA_Complex *b = qlua_checkComplex(L, 2);
    mLatColVec *c = qlua_newLatColVec(L);

    CALL_QDP(L);
    QDP_V_eq_c_times_V(c->ptr, b, a->ptr, *qCurrent);

    return 1;
}

static int
q_c_mul_V(lua_State *L)
{
    QLA_Complex *a = qlua_checkComplex(L, 1);
    mLatColVec *b = qlua_checkLatColVec(L, 2);
    mLatColVec *c = qlua_newLatColVec(L);

    CALL_QDP(L);
    QDP_V_eq_c_times_V(c->ptr, a, b->ptr, *qCurrent);

    return 1;
}

static struct {
    QLA_Real *a;
    QLA_ColorVector *b;
} RVmul_args; /* YYY global state */

static void
do_RVmul(QLA_ColorVector *r, int idx)
{
    QLA_V_eq_r_times_V(r, &RVmul_args.a[idx], &RVmul_args.b[idx]);
}

static void
X_V_eq_R_times_V(QDP_ColorVector *r, QDP_Real *a, QDP_ColorVector *b,
                 QDP_Subset s)
{
    RVmul_args.a = QDP_expose_R(a);
    RVmul_args.b = QDP_expose_V(b);
    QDP_V_eq_funci(r, do_RVmul, s);
    QDP_reset_R(a);
    QDP_reset_V(b);
    RVmul_args.a = 0;
    RVmul_args.b = 0;
}

static int
q_R_mul_V(lua_State *L)
{
    mLatReal *a = qlua_checkLatReal(L, 1);
    mLatColVec *b = qlua_checkLatColVec(L, 2);
    mLatColVec *r = qlua_newLatColVec(L);

    CALL_QDP(L);
    X_V_eq_R_times_V(r->ptr, a->ptr, b->ptr, *qCurrent);

    return 1;
}

static int
q_V_mul_R(lua_State *L)
{
    mLatColVec *a = qlua_checkLatColVec(L, 1);
    mLatReal *b = qlua_checkLatReal(L, 2);
    mLatColVec *r = qlua_newLatColVec(L);

    CALL_QDP(L);
    X_V_eq_R_times_V(r->ptr, b->ptr, a->ptr, *qCurrent);

    return 1;
}

static struct {
    QLA_Complex *a;
    QLA_ColorVector *b;
} CVmul_args; /* YYY global state */

static void
do_CVmul(QLA_ColorVector *r, int idx)
{
    QLA_V_eq_c_times_V(r, &CVmul_args.a[idx], &CVmul_args.b[idx]);
}

static void
X_V_eq_C_times_V(QDP_ColorVector *r, QDP_Complex *a, QDP_ColorVector *b,
                 QDP_Subset s)
{
    CVmul_args.a = QDP_expose_C(a);
    CVmul_args.b = QDP_expose_V(b);
    QDP_V_eq_funci(r, do_CVmul, s);
    QDP_reset_C(a);
    QDP_reset_V(b);
    CVmul_args.a = 0;
    CVmul_args.b = 0;
}

static int
q_C_mul_V(lua_State *L)
{
    mLatComplex *a = qlua_checkLatComplex(L, 1);
    mLatColVec *b = qlua_checkLatColVec(L, 2);
    mLatColVec *r = qlua_newLatColVec(L);

    CALL_QDP(L);
    X_V_eq_C_times_V(r->ptr, a->ptr, b->ptr, *qCurrent);

    return 1;
}

static int
q_V_mul_C(lua_State *L)
{
    mLatColVec *a = qlua_checkLatColVec(L, 1);
    mLatComplex *b = qlua_checkLatComplex(L, 2);
    mLatColVec *r = qlua_newLatColVec(L);

    CALL_QDP(L);
    X_V_eq_C_times_V(r->ptr, b->ptr, a->ptr, *qCurrent);

    return 1;
}

static int
q_V_norm2(lua_State *L)
{
    mLatColVec *a = qlua_checkLatColVec(L, 1);
    QLA_Real n;

    CALL_QDP(L);
    QDP_r_eq_norm2_V(&n, a->ptr, *qCurrent);
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

    CALL_QDP(L);
    QDP_V_eq_sV(r->ptr, a->ptr, shift, dir, *qCurrent);

    return 1;
}

static int
q_V_conj(lua_State *L)
{
    mLatColVec *a = qlua_checkLatColVec(L, 1);
    mLatColVec *r = qlua_newLatColVec(L);

    CALL_QDP(L);
    QDP_V_eq_conj_V(r->ptr, a->ptr, *qCurrent);

    return 1;
}

static int
q_V_set(lua_State *L)
{
    mLatColVec *r = qlua_checkLatColVec(L, 1);
    mLatColVec *a = qlua_checkLatColVec(L, 2);

    CALL_QDP(L);
    QDP_V_eq_V(r->ptr, a->ptr, *qCurrent);
    lua_pop(L, 1);

    return 1;
}

static int
q_V_neg(lua_State *L)
{
    mLatColVec *a = qlua_checkLatColVec(L, 1);
    mLatColVec *r = qlua_newLatColVec(L);
    QLA_Real m1 = -1;

    CALL_QDP(L);
    QDP_V_eq_r_times_V(r->ptr, &m1, a->ptr, *qCurrent);

    return 1;
}
static int
q_V_div_r(lua_State *L)
{
    mLatColVec *a = qlua_checkLatColVec(L, 1);
    QLA_Real b = 1 / luaL_checknumber(L, 2);
    mLatColVec *c = qlua_newLatColVec(L);

    CALL_QDP(L);
    QDP_V_eq_r_times_V(c->ptr, &b, a->ptr, *qCurrent);

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

    CALL_QDP(L);
    QLA_real(s) = n * QLA_real(*b);
    QLA_imag(s) = -n * QLA_imag(*b);
    QDP_V_eq_c_times_V(c->ptr, &s, a->ptr, *qCurrent);

    return 1;
}


static int
q_latcolvec(lua_State *L)
{
    switch (lua_gettop(L)) {
    case 1: {
        mLatColVec *v = qlua_newLatColVec(L);

        CALL_QDP(L);
        QDP_V_eq_zero(v->ptr, *qCurrent);

        return 1;
    }
    case 2: {
        mLatColVec *a = qlua_checkLatColVec(L, 2);
        mLatColVec *r = qlua_newLatColVec(L);
        
        CALL_QDP(L);
        QDP_V_eq_V(r->ptr, a->ptr, *qCurrent);
        
        return 1;
    }
    case 3: {
        mLatComplex *c = qlua_checkLatComplex(L, 2);
        int a = qlua_checkcolorindex(L, 3);
        mLatColVec *r = qlua_newLatColVec(L);

        CALL_QDP(L);
        QDP_V_eq_elem_C(r->ptr, c->ptr, a, *qCurrent);

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
    qlua_reg_add(qLatColVec,  qLatColVec,  q_V_add_V);
    qlua_reg_sub(qLatColVec,  qLatColVec,  q_V_sub_V);
    qlua_reg_mul(qReal,       qLatColVec,  q_r_mul_V);
    qlua_reg_mul(qLatColVec,  qReal,       q_V_mul_r);
    qlua_reg_mul(qComplex,    qLatColVec,  q_c_mul_V);
    qlua_reg_mul(qLatColVec,  qComplex,    q_V_mul_c);
    qlua_reg_mul(qLatReal,    qLatColVec,  q_R_mul_V);
    qlua_reg_mul(qLatColVec,  qLatReal,    q_V_mul_R);
    qlua_reg_mul(qLatComplex, qLatColVec,  q_C_mul_V);
    qlua_reg_mul(qLatColVec,  qLatComplex, q_V_mul_C);
    qlua_reg_div(qLatColVec,  qReal,       q_V_div_r);
    qlua_reg_div(qLatColVec,  qComplex,    q_V_div_c);
    qlua_reg_dot(qLatColVec,  q_V_dot);

    return 0;
}

int
fini_latcolvec(lua_State *L)
{
    return 0;
}
