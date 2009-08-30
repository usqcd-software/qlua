#include <qlua.h>
#include <qcomplex.h>
#include <qgamma.h>
#include <qdp.h>

const char mtnGamma[] = "qcd.mtGamma";

enum {
    qG_z,  /* zero */
    qG_p,  /* +1 */
    qG_m,  /* -1 */
    qG_r,  /* real */
    qG_c,  /* complex */
    qG_t
};
#define Gi(a,b)  ((a)*qG_t+(b))

typedef struct mGamma_s {
    int t;
    QLA_Real r;
    QLA_Complex c;
} mGamma;

typedef struct mClifford_s {
    mGamma g[16];
} mClifford;

static const char *gn[] = {
    NULL, "g0",  "g1",  "g01",  "g2",  "g02", "g12",   "g012",
    "g3", "g03", "g13", "g013", "g23", "g023", "g123", "g5"};

/* generated from QDP gamma operations */
#define gm_k(i,j)  (15 & gm[(i)][(j)])
#define gm_s(i,j)  (gm[(i)][(j)] > 15)
static const char gm[16][16] = {
  { 0,  1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11, 12, 13, 14, 15},
  { 1,  0,  3,  2,  5,  4,  7,  6,  9,  8, 11, 10, 13, 12, 15, 14},
  { 2, 19,  0, 17,  6, 23,  4, 21, 10, 27,  8, 25, 14, 31, 12, 29},
  { 3, 18,  1, 16,  7, 22,  5, 20, 11, 26,  9, 24, 15, 30, 13, 28},
  { 4, 21, 22,  7,  0, 17, 18,  3, 12, 29, 30, 15,  8, 25, 26, 11},
  { 5, 20, 23,  6,  1, 16, 19,  2, 13, 28, 31, 14,  9, 24, 27, 10},
  { 6,  7, 20, 21,  2,  3, 16, 17, 14, 15, 28, 29, 10, 11, 24, 25},
  { 7,  6, 21, 20,  3,  2, 17, 16, 15, 14, 29, 28, 11, 10, 25, 24},
  { 8, 25, 26, 11, 28, 13, 14, 31,  0, 17, 18,  3, 20,  5,  6, 23},
  { 9, 24, 27, 10, 29, 12, 15, 30,  1, 16, 19,  2, 21,  4,  7, 22},
  {10, 11, 24, 25, 30, 31, 12, 13,  2,  3, 16, 17, 22, 23,  4,  5},
  {11, 10, 25, 24, 31, 30, 13, 12,  3,  2, 17, 16, 23, 22,  5,  4},
  {12, 13, 14, 15, 24, 25, 26, 27,  4,  5,  6,  7, 16, 17, 18, 19},
  {13, 12, 15, 14, 25, 24, 27, 26,  5,  4,  7,  6, 17, 16, 19, 18},
  {14, 31, 12, 29, 26, 11, 24,  9,  6, 23,  4, 21, 18,  3, 16,  1},
  {15, 30, 13, 28, 27, 10, 25,  8,  7, 22,  5, 20, 19,  2, 17,  0}};

#define g_OP(t,a,b) ((t)[(a)*qG_t + (b)])

#if 0
static void
gn_r(mGamma *r)
{
    if (r->r == 0)
        r->t = qG_z;
    else if (r->r == 1)
        r->t = qG_p;
    else if (r->r == -1)
        r->t = qG_m;
    else
        r->t = qG_r;
}

static void
gn_c(mGamma *r)
{
    if (QLA_imag(r->c) == 0) {
        r->r = QLA_real(r->c);
    } else {
        r->t = qG_c;
    }
}
#endif

/* additions */
typedef void (*gadd)(mGamma *r, const mGamma *a, const mGamma *b);

static void gadd_z_z(mGamma *r, const mGamma *a, const mGamma *b) {r->t = qG_z;}
static void gadd_z_p(mGamma *r, const mGamma *a, const mGamma *b) {r->t = qG_p;}
static void gadd_z_m(mGamma *r, const mGamma *a, const mGamma *b) {r->t = qG_m;}
static void gadd_z_r(mGamma *r, const mGamma *a, const mGamma *b) {*r = *b;}
static void gadd_z_c(mGamma *r, const mGamma *a, const mGamma *b) {*r = *b;}

static void gadd_p_z(mGamma *r, const mGamma *a, const mGamma *b) {r->t = qG_p;}
static void gadd_p_p(mGamma *r, const mGamma *a, const mGamma *b)
{
    r->r = 2;
    r->t = qG_r;
}
static void gadd_p_m(mGamma *r, const mGamma *a, const mGamma *b) {r->t = qG_z;}
static void gadd_p_r(mGamma *r, const mGamma *a, const mGamma *b)
{
    r->r = 1 + b->r;
    r->t = qG_r;
}
static void gadd_p_c(mGamma *r, const mGamma *a, const mGamma *b)
{
    QLA_real(r->c) = 1 + QLA_real(b->c);
    QLA_imag(r->c) = QLA_imag(b->c);
    r->t = qG_c;
}

static void gadd_m_z(mGamma *r, const mGamma *a, const mGamma *b) {r->t = qG_m;}
static void gadd_m_p(mGamma *r, const mGamma *a, const mGamma *b) {r->t = qG_z;}
static void gadd_m_m(mGamma *r, const mGamma *a, const mGamma *b)
{
    r->r = -2;
    r->t = qG_r;
}
static void gadd_m_r(mGamma *r, const mGamma *a, const mGamma *b)
{
    r->r = -1 + b->r;
    r->t = qG_r;
}
static void gadd_m_c(mGamma *r, const mGamma *a, const mGamma *b)
{
    QLA_real(r->c) = -1 + QLA_real(b->c);
    QLA_imag(r->c) = QLA_imag(b->c);
    r->t = qG_c;
}

static void gadd_r_z(mGamma *r, const mGamma *a, const mGamma *b) {*r = *a;}
static void gadd_r_p(mGamma *r, const mGamma *a, const mGamma *b)
{
    r->r = a->r + 1;
    r->t = qG_r;
}
static void gadd_r_m(mGamma *r, const mGamma *a, const mGamma *b)
{
    r->r = a->r - 1;
    r->t = qG_r;
}
static void gadd_r_r(mGamma *r, const mGamma *a, const mGamma *b)
{
    r->r = a->r + b->r;
    r->t = qG_r;
}
static void gadd_r_c(mGamma *r, const mGamma *a, const mGamma *b)
{
    QLA_real(r->c) = a->r + QLA_real(b->c);
    QLA_imag(r->c) = QLA_imag(b->c);
    r->t = qG_c;
}

static void gadd_c_z(mGamma *r, const mGamma *a, const mGamma *b) {*r = *a;}
static void gadd_c_p(mGamma *r, const mGamma *a, const mGamma *b)
{
    QLA_real(r->c) = QLA_real(a->c) + 1;
    QLA_imag(r->c) = QLA_imag(a->c);
    r->t = qG_c;
}
static void gadd_c_m(mGamma *r, const mGamma *a, const mGamma *b)
{
    QLA_real(r->c) = QLA_real(a->c) - 1;
    QLA_imag(r->c) = QLA_imag(a->c);
    r->t = qG_c;
}
static void gadd_c_r(mGamma *r, const mGamma *a, const mGamma *b)
{
    QLA_real(r->c) = QLA_real(a->c) + b->r;
    QLA_imag(r->c) = QLA_imag(a->c);
    r->t = qG_c;
}
static void gadd_c_c(mGamma *r, const mGamma *a, const mGamma *b)
{
    QLA_real(r->c) = QLA_real(a->c) + QLA_real(b->c);
    QLA_imag(r->c) = QLA_imag(a->c) + QLA_imag(b->c);
    r->t = qG_c;
}

static const gadd t_add[] = {
    gadd_z_z, gadd_z_p, gadd_z_m, gadd_z_r, gadd_z_c,
    gadd_p_z, gadd_p_p, gadd_p_m, gadd_p_r, gadd_p_c,
    gadd_m_z, gadd_m_p, gadd_m_m, gadd_m_r, gadd_m_c,
    gadd_r_z, gadd_r_p, gadd_r_m, gadd_r_r, gadd_r_c,
    gadd_c_z, gadd_c_p, gadd_c_m, gadd_c_r, gadd_c_c
};

/* multiplications */
typedef void (*gmul)(mGamma *r, const mGamma *a, const mGamma *b);
static void gmul_z_z(mGamma *r, const mGamma *a, const mGamma *b) {r->t = qG_z;}
static void gmul_z_p(mGamma *r, const mGamma *a, const mGamma *b) {r->t = qG_z;}
static void gmul_z_m(mGamma *r, const mGamma *a, const mGamma *b) {r->t = qG_z;}
static void gmul_z_r(mGamma *r, const mGamma *a, const mGamma *b) {r->t = qG_z;}
static void gmul_z_c(mGamma *r, const mGamma *a, const mGamma *b) {r->t = qG_z;}

static void gmul_p_z(mGamma *r, const mGamma *a, const mGamma *b) {r->t = qG_z;}
static void gmul_p_p(mGamma *r, const mGamma *a, const mGamma *b) {r->t = qG_p;}
static void gmul_p_m(mGamma *r, const mGamma *a, const mGamma *b) {r->t = qG_m;}
static void gmul_p_r(mGamma *r, const mGamma *a, const mGamma *b) {*r = *b;}
static void gmul_p_c(mGamma *r, const mGamma *a, const mGamma *b) {*r = *b;}

static void gmul_m_z(mGamma *r, const mGamma *a, const mGamma *b) {r->t = qG_z;}
static void gmul_m_p(mGamma *r, const mGamma *a, const mGamma *b) {r->t = qG_m;}
static void gmul_m_m(mGamma *r, const mGamma *a, const mGamma *b) {r->t = qG_p;}
static void gmul_m_r(mGamma *r, const mGamma *a, const mGamma *b)
{
    r->r = -b->r;
    r->t = qG_r;
}
static void gmul_m_c(mGamma *r, const mGamma *a, const mGamma *b)
{
    QLA_real(r->c) = -QLA_real(b->c);
    QLA_imag(r->c) = -QLA_imag(b->c);
    r->t = qG_c;
}

static void gmul_r_z(mGamma *r, const mGamma *a, const mGamma *b) {r->t = qG_z;}
static void gmul_r_p(mGamma *r, const mGamma *a, const mGamma *b) {*r = *a;}
static void gmul_r_m(mGamma *r, const mGamma *a, const mGamma *b)
{
    r->r = -a->r;
    r->t = qG_r;
}
static void gmul_r_r(mGamma *r, const mGamma *a, const mGamma *b)
{
    r->r = a->r * b->r;
    r->t = qG_r;
}
static void gmul_r_c(mGamma *r, const mGamma *a, const mGamma *b)
{
    QLA_c_eq_c_times_r(r->c, b->c, a->r);
    r->t = qG_c;
}

static void gmul_c_z(mGamma *r, const mGamma *a, const mGamma *b) {r->t = qG_z;}
static void gmul_c_p(mGamma *r, const mGamma *a, const mGamma *b) {*r = *a;}
static void gmul_c_m(mGamma *r, const mGamma *a, const mGamma *b)
{
    QLA_real(r->c) = -QLA_real(a->c);
    QLA_imag(r->c) = -QLA_imag(a->c);
    r->t = qG_c;
}
static void gmul_c_r(mGamma *r, const mGamma *a, const mGamma *b)
{
    QLA_c_eq_c_times_r(r->c, a->c, b->r);
    r->t = qG_c;
}
static void gmul_c_c(mGamma *r, const mGamma *a, const mGamma *b)
{
    QLA_c_eq_c_times_c(r->c, a->c, b->c);
    r->t = qG_c;
}

static const gmul t_mul[] = {
    gmul_z_z, gmul_z_p, gmul_z_m, gmul_z_r, gmul_z_c,
    gmul_p_z, gmul_p_p, gmul_p_m, gmul_p_r, gmul_p_c,
    gmul_m_z, gmul_m_p, gmul_m_m, gmul_m_r, gmul_m_c,
    gmul_r_z, gmul_r_p, gmul_r_m, gmul_r_r, gmul_r_c,
    gmul_c_z, gmul_c_p, gmul_c_m, gmul_c_r, gmul_c_c
};

mClifford *
qlua_newClifford(lua_State *L)
{
    mClifford *g;
    int i;

    g = lua_newuserdata(L, sizeof (mClifford));
    luaL_getmetatable(L, mtnGamma);
    lua_setmetatable(L, -2);

    for (i = 0; i < 16; i++)
        g->g[i].t = qG_z;

    return g;
}

mClifford *
qlua_checkClifford(lua_State *L, int idx)
{
    void *v = luaL_checkudata(L, idx, mtnGamma);
    
    luaL_argcheck(L, v != 0, idx, "qcd.gamma expected");
    
    return v;
}

static int
q_g_fmt(lua_State *L)
{
    luaL_Buffer b;
    mClifford *v = qlua_checkClifford(L, 1);
    char fmt[72];
    int fst, i;

    luaL_buffinit(L, &b);
    luaL_addstring(&b, "gamma[");
    switch (v->g[0].t) {
    case qG_z:
        fst = 1;
        break;
    case qG_p:
        fst = 0;
        luaL_addstring(&b, "1");
        break;
    case qG_m:
        fst = 0;
        luaL_addstring(&b, "-1");
        break;
    case qG_r:
        fst = 0;
        sprintf(fmt, "%g", v->g[0].r);
        luaL_addstring(&b, fmt);
        break;
    case qG_c:
        fst = 0;
        sprintf(fmt, "complex(%g,%g)",
                QLA_real(v->g[0].c), QLA_imag(v->g[0].c));
        luaL_addstring(&b, fmt);
        break;
    default:
        return luaL_error(L, "can't happen");
    }
    for (i = 1; i < 16; i++) {
        switch (v->g[i].t) {
        case qG_z:
            continue;
        case qG_p:
            if (fst == 0) {
                luaL_addstring(&b, "+");
            }
            fst = 0;
            break;
        case qG_m:
            luaL_addstring(&b, "-");
            fst = 0;
            break;
        case qG_r:
            if (fst == 0) {
                sprintf(fmt, "%+g*", v->g[i].r);
            } else {
                sprintf(fmt, "%g*", v->g[i].r);
            }
            luaL_addstring(&b, fmt);
            fst = 0;
            break;
        case qG_c:
            if (fst == 0) {
                sprintf(fmt, "+complex(%g,%g)*", 
                        QLA_real(v->g[i].c), QLA_imag(v->g[i].c));
            } else {
                sprintf(fmt, "complex(%g,%g)*", 
                        QLA_real(v->g[i].c), QLA_imag(v->g[i].c));
            }
            luaL_addstring(&b, fmt);
            fst = 0;
            break;
        }
        luaL_addstring(&b, gn[i]);
    }

    luaL_addstring(&b, "]");
    luaL_pushresult(&b);
    qlua_checkClifford(L, 1);

    return 1;
}

static void
g_norm(mGamma *r)
{
    switch (r->t) {
    case qG_z: case qG_p: case qG_m: break;
    case qG_c:
        if (QLA_imag(r->c) != 0) break;
        r->r = QLA_real(r->c);
        r->t = qG_r;
        /* through */
    case qG_r:
        if (r->r == 0)
            r->t = qG_z;
        else if (r->r == 1)
            r->t = qG_p;
        else if (r->r == -1)
            r->t = qG_m;
        break;
    }
}

static int
c_norm(lua_State *L, mClifford *x)
{
    int i, m;

    for (m = -1, i = 16; i--;) {
        g_norm(&x->g[i]);
        if ((m == -1) && x->g[i].t != qG_z)
            m = i;
    }
    if (m < 1)
        switch (x->g[0].t) {
        case qG_z:
            lua_pushnumber(L, 0);
            break;
        case qG_p:
            lua_pushnumber(L, 1);
            break;
        case qG_m:
            lua_pushnumber(L, -1);
            break;
        case qG_r:
            lua_pushnumber(L, x->g[0].r);
        case qG_c: {
            QLA_Complex *z = qlua_newComplex(L);
            QLA_c_eq_c(*z, x->g[0].c);
            break;
        }
        }
    return 1;
}

static int
q_g_add_g(lua_State *L)
{
    mClifford *x = qlua_checkClifford(L, 1);
    mClifford *y = qlua_checkClifford(L, 2);
    mClifford *r = qlua_newClifford(L);
    int i;

    for (i = 0; i < 16; i++)
        g_OP(t_add, x->g[i].t, y->g[i].t)(&r->g[i], &x->g[i], &y->g[i]);

    return c_norm(L, r);
}

static int
q_r_add_g(lua_State *L)
{
    QLA_Real b = luaL_checknumber(L, 1);
    mClifford *x = qlua_checkClifford(L, 2);
    mClifford *r = qlua_newClifford(L);
    int i;
    mGamma v;
    
    v.t = qG_r;
    v.r = b;
    g_OP(t_add, x->g[0].t, v.t)(&r->g[0], &x->g[0], &v);
    g_norm(&r->g[0]);
    for (i = 1; i < 16; i++)
        r->g[i] = x->g[i];

    return 1;
}

static int
q_g_add_r(lua_State *L)
{
    mClifford *x = qlua_checkClifford(L, 1);
    QLA_Real b = luaL_checknumber(L, 2);
    mClifford *r = qlua_newClifford(L);
    int i;
    mGamma v;
    
    v.t = qG_r;
    v.r = b;
    g_OP(t_add, x->g[0].t, v.t)(&r->g[0], &x->g[0], &v);
    g_norm(&r->g[0]);
    for (i = 1; i < 16; i++)
        r->g[i] = x->g[i];

    return 1;
}

static int
q_c_add_g(lua_State *L)
{
    QLA_Complex *b = qlua_checkComplex(L, 1);
    mClifford *x = qlua_checkClifford(L, 2);
    mClifford *r = qlua_newClifford(L);
    int i;
    mGamma v;
    
    v.t = qG_c;
    QLA_c_eq_c(v.c, *b);
    g_OP(t_add, x->g[0].t, v.t)(&r->g[0], &x->g[0], &v);
    g_norm(&r->g[0]);
    for (i = 1; i < 16; i++)
        r->g[i] = x->g[i];

    return 1;
}

static int
q_g_add_c(lua_State *L)
{
    mClifford *x = qlua_checkClifford(L, 1);
    QLA_Complex *b = qlua_checkComplex(L, 2);
    mClifford *r = qlua_newClifford(L);
    int i;
    mGamma v;
    
    v.t = qG_c;
    QLA_c_eq_c(v.c, *b);
    g_OP(t_add, x->g[0].t, v.t)(&r->g[0], &x->g[0], &v);
    g_norm(&r->g[0]);
    for (i = 1; i < 16; i++)
        r->g[i] = x->g[i];

    return 1;
}

static void
g_neg(mGamma *r, const mGamma *a)
{
    switch (a->t) {
    case qG_z: r->t = qG_z; break;
    case qG_p: r->t = qG_m; break;
    case qG_m: r->t = qG_p; break;
    case qG_r: r->t = qG_r; r->r = -a->r; break;
    case qG_c:
        r->t = qG_c;
        QLA_real(r->c) = -QLA_real(a->c);
        QLA_imag(r->c) = -QLA_imag(a->c);
        break;
    }
}

static int
q_g_sub_g(lua_State *L)
{
    mClifford *x = qlua_checkClifford(L, 1);
    mClifford *y = qlua_checkClifford(L, 2);
    mClifford *r = qlua_newClifford(L);
    mGamma t;
    int i;

    for (i = 0; i < 16; i++) {
        g_neg(&t, &y->g[i]);
        g_OP(t_add, x->g[i].t, t.t)(&r->g[i], &x->g[i], &t);
    }

    return c_norm(L, r);
}

static int
q_r_sub_g(lua_State *L)
{
    QLA_Real b = luaL_checknumber(L, 1);
    mClifford *x = qlua_checkClifford(L, 2);
    mClifford *r = qlua_newClifford(L);
    int i;
    mGamma v, t;
    
    v.t = qG_r;
    v.r = b;
    g_neg(&t, &x->g[0]);
    g_OP(t_add, t.t, v.t)(&r->g[0], &t, &v);
    g_norm(&r->g[0]);
    for (i = 1; i < 16; i++)
        g_neg(&r->g[i], &x->g[i]);

    return 1;
}

static int
q_g_sub_r(lua_State *L)
{
    mClifford *x = qlua_checkClifford(L, 1);
    QLA_Real b = luaL_checknumber(L, 2);
    mClifford *r = qlua_newClifford(L);
    int i;
    mGamma v;
    
    v.t = qG_r;
    v.r = -b;
    g_OP(t_add, x->g[0].t, v.t)(&r->g[0], &x->g[0], &v);
    g_norm(&r->g[0]);
    for (i = 1; i < 16; i++)
        r->g[i] = x->g[i];

    return 1;
}

static int
q_c_sub_g(lua_State *L)
{
    QLA_Complex *b = qlua_checkComplex(L, 1);
    mClifford *x = qlua_checkClifford(L, 2);
    mClifford *r = qlua_newClifford(L);
    int i;
    mGamma v, t;
    
    v.t = qG_c;
    QLA_c_eq_c(v.c, *b);
    g_neg(&t, &x->g[0]);
    g_OP(t_add, t.t, v.t)(&r->g[0], &t, &v);
    g_norm(&r->g[0]);
    for (i = 1; i < 16; i++)
        g_neg(&r->g[i], &x->g[i]);

    return 1;
}

static int q_g_sub_c(lua_State *L)
{
    mClifford *x = qlua_checkClifford(L, 1);
    QLA_Complex *b = qlua_checkComplex(L, 2);
    mClifford *r = qlua_newClifford(L);
    int i;
    mGamma v;
    
    v.t = qG_c;
    QLA_c_eqm_c(v.c, *b);
    g_OP(t_add, x->g[0].t, v.t)(&r->g[0], &x->g[0], &v);
    g_norm(&r->g[0]);
    for (i = 1; i < 16; i++)
        r->g[i] = x->g[i];

    return 1;
}

static int
q_g_mul_g(lua_State *L)
{
    mClifford *x = qlua_checkClifford(L, 1);
    mClifford *y = qlua_checkClifford(L, 2);
    mClifford *r = qlua_newClifford(L);
    int i, j;

    for (i = 0; i < 16; i++) {
        for (j = 0; j < 16; j++) {
            int s = gm_s(i, j);
            int k = gm_k(i, j);
            mGamma t;
            g_OP(t_mul, x->g[i].t, y->g[j].t)(&t, &x->g[i], &y->g[j]);
            if (s) {
                g_neg(&t, &t);
            }
            g_OP(t_add, r->g[k].t, t.t)(&r->g[k], &r->g[k], &t);
        }
    }
    
    return c_norm(L, r);
}

static int
q_r_mul_g(lua_State *L)
{
    QLA_Real a = luaL_checknumber(L, 1);
    mClifford *b = qlua_checkClifford(L, 2);
    mClifford *r = qlua_newClifford(L);
    int i;
    mGamma v;

    v.t = qG_r;
    v.r = a;
    for (i = 0; i < 16; i++)
        g_OP(t_mul, v.t, b->g[i].t)(&r->g[i], &v, &b->g[i]);
    
    return c_norm(L, r);
}

static int q_g_mul_r(lua_State *L)
{
    mClifford *b = qlua_checkClifford(L, 1);
    QLA_Real a = luaL_checknumber(L, 2);
    mClifford *r = qlua_newClifford(L);
    int i;
    mGamma v;

    v.t = qG_r;
    v.r = a;
    for (i = 0; i < 16; i++)
        g_OP(t_mul, b->g[i].t, v.t)(&r->g[i], &b->g[i], &v);
    
    return c_norm(L, r);
}

static int q_c_mul_g(lua_State *L)
{
    QLA_Complex *a = qlua_checkComplex(L, 1);
    mClifford *b = qlua_checkClifford(L, 2);
    mClifford *r = qlua_newClifford(L);
    int i;
    mGamma v;

    v.t = qG_c;
    QLA_c_eq_c(v.c, *a);
    for (i = 0; i < 16; i++)
        g_OP(t_mul, v.t, b->g[i].t)(&r->g[i], &v, &b->g[i]);
    
    return c_norm(L, r);
}

static int
q_g_mul_c(lua_State *L)
{
    mClifford *b = qlua_checkClifford(L, 1);
    QLA_Complex *a = qlua_checkComplex(L, 2);
    mClifford *r = qlua_newClifford(L);
    int i;
    mGamma v;

    v.t = qG_c;
    QLA_c_eq_c(v.c, *a);
    for (i = 0; i < 16; i++)
        g_OP(t_mul, b->g[i].t, v.t)(&r->g[i], &b->g[i], &v);
    
    return c_norm(L, r);
}

/* XXX */ static int q_g_div_r(lua_State *L) { return 0; }
/* XXX */ static int q_g_div_c(lua_State *L) { return 0; }

/* XXX */ static int q_g_mul_D(lua_State *L) { return 0; }
/* XXX */ static int q_g_mul_P(lua_State *L) { return 0; }
/* XXX */ static int q_P_mul_g(lua_State *L) { return 0; }

static int
q_g_neg(lua_State *L)
{
    mClifford *x = qlua_checkClifford(L, 1);
    mClifford *r = qlua_newClifford(L);
    int i;

    for (i = 0; i < 16; i++) {
        g_neg(&r->g[i], &x->g[i]);
    }
    
    return 1;
}

static int
q_gamma(lua_State *L)
{
    switch (qlua_gettype(L, 1)) {
    case qTable: {
        int mu = qlua_gammaindex(L, 1);
        int n = qlua_gammabinary(L, 1);

        if ((n == -1) && (mu != -1)) {
            mClifford *x = qlua_newClifford(L);
            
            n = mu == 5? 15: (1 << mu);
            x->g[n].t = qG_p;

            return 1;
        } else if ((n != -1) && (mu == -1)) {
            mClifford *x = qlua_newClifford(L);
            
            x->g[n].t = qG_p;

            return c_norm(L, x);
        }
        break;
    }
    }
    return luaL_error(L, "bad arguments");
}

static struct luaL_Reg mtGamma[] = {
    { "__tostring",        q_g_fmt },
    { "__unm",             q_g_neg },
    { "__add",             qlua_add },
    { "__sub",             qlua_sub },
    { "__mul",             qlua_mul },
    { NULL,                NULL }
};

static struct luaL_Reg fGamma[] = {
    { "gamma",             q_gamma },
    { NULL,                NULL }
};

int
init_gamma(lua_State *L)
{
    lua_getglobal(L, "_G");
    luaL_register(L, NULL, fGamma);
    qlua_metatable(L, mtnGamma, mtGamma);
    qlua_reg_add(qGamma,      qGamma,      q_g_add_g);
    qlua_reg_add(qReal,       qGamma,      q_r_add_g);
    qlua_reg_add(qGamma,      qReal,       q_g_add_r);
    qlua_reg_add(qComplex,    qGamma,      q_c_add_g);
    qlua_reg_add(qGamma,      qComplex,    q_g_add_c);
    qlua_reg_sub(qGamma,      qGamma,      q_g_sub_g);
    qlua_reg_sub(qReal,       qGamma,      q_r_sub_g);
    qlua_reg_sub(qGamma,      qReal,       q_g_sub_r);
    qlua_reg_sub(qComplex,    qGamma,      q_c_sub_g);
    qlua_reg_sub(qGamma,      qComplex,    q_g_sub_c);
    qlua_reg_mul(qGamma,      qGamma,      q_g_mul_g);
    qlua_reg_mul(qReal,       qGamma,      q_r_mul_g);
    qlua_reg_mul(qGamma,      qReal,       q_g_mul_r);
    qlua_reg_mul(qComplex,    qGamma,      q_c_mul_g);
    qlua_reg_mul(qGamma,      qComplex,    q_g_mul_c);
    qlua_reg_mul(qGamma,      qLatDirFerm, q_g_mul_D);
    qlua_reg_mul(qGamma,      qLatDirProp, q_g_mul_P);
    qlua_reg_mul(qLatDirProp, qGamma,      q_P_mul_g);
    qlua_reg_div(qGamma,      qReal,       q_g_div_r);
    qlua_reg_div(qGamma,      qComplex,    q_g_div_c);

    return 0;
}

int
fini_gamma(lua_State *L)
{
    return 0;
}

