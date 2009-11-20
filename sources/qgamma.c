#include <qlua.h>                                                    /* DEPS */
#include <qcomplex.h>                                                /* DEPS */
#include <qgamma.h>                                                  /* DEPS */
#include <lattice.h>                                                 /* DEPS */
#include <latdirferm.h>                                              /* DEPS */
#include <latdirprop.h>                                              /* DEPS */
#include <qdp.h>
#include <math.h>
#include <string.h>

const char mtnGamma[] = "qlua.mtGamma";

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

static char gconj[] = {0, 0, 0, 1, 0, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 0};

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

#define add_r(x,y,z) \
   QLA_real(x->c) = QLA_real(y->c) + z; QLA_imag(x->c) = QLA_imag(y->c)
#define add_c(x,y,z) QLA_c_eq_c_plus_c(x->c, y->c, z->c)

#define g_ID(a,b) ((a)*qG_t + (b))

static void
g_add(mGamma *r, const mGamma *a, const mGamma *b)
{
    switch (g_ID(a->t, b->t)) {
    case g_ID(qG_z, qG_z): r->t = qG_z; break;
    case g_ID(qG_z, qG_p): r->t = qG_p; break;
    case g_ID(qG_z, qG_m): r->t = qG_m; break;
    case g_ID(qG_z, qG_r): r->r = b->r; r->t = qG_r; break;
    case g_ID(qG_z, qG_c): QLA_c_eq_c(r->c, b->c); r->t = qG_c; break;
    case g_ID(qG_p, qG_z): r->t = qG_p; break;
    case g_ID(qG_p, qG_p): r->r = 2; r->t = qG_r; break;
    case g_ID(qG_p, qG_m): r->t = qG_z; break;
    case g_ID(qG_p, qG_r): r->r = 1 + b->r; r->t = qG_r; break;
    case g_ID(qG_p, qG_c): add_r(r, b, 1); r->t = qG_c; break;
    case g_ID(qG_m, qG_z): r->t = qG_m; break;
    case g_ID(qG_m, qG_p): r->t = qG_z; break;
    case g_ID(qG_m, qG_m): r->r = -2; r->t = qG_r; break;
    case g_ID(qG_m, qG_r): r->r = b->r - 1; r->t = qG_r; break;
    case g_ID(qG_m, qG_c): add_r(r, b, -1); r->t = qG_c; break;
    case g_ID(qG_r, qG_z): r->r = a->r; r->t = qG_r; break;
    case g_ID(qG_r, qG_p): r->r = a->r + 1; r->t = qG_r; break;
    case g_ID(qG_r, qG_m): r->r = a->r - 1; r->t = qG_r; break;
    case g_ID(qG_r, qG_r): r->r = a->r + b->r; r->t = qG_r; break;
    case g_ID(qG_r, qG_c): add_r(r, b, a->r); r->t=qG_c; break;
    case g_ID(qG_c, qG_z): QLA_c_eq_c(r->c, a->c); r->t = qG_c; break;
    case g_ID(qG_c, qG_p): add_r(r, b,  1); r->t = qG_c; break;
    case g_ID(qG_c, qG_m): add_r(r, b, -1); r->t = qG_c; break;
    case g_ID(qG_c, qG_r): add_r(r, a, b->r); r->t = qG_c; break;
    case g_ID(qG_c, qG_c): add_c(r, a, b); r->t = qG_c; break;
    }
}

#define mul_r(x,y,z) QLA_c_eq_c_times_r(x->c, y->c, z)
#define mul_c(x,y,z) QLA_c_eq_c_times_c(x->c, y->c, z->c)

static void
g_mul(mGamma *r, const mGamma *a, const mGamma *b)
{
    switch (g_ID(a->t, b->t)) {
    case g_ID(qG_z, qG_z): r->t = qG_z; break;
    case g_ID(qG_z, qG_p): r->t = qG_z; break;
    case g_ID(qG_z, qG_m): r->t = qG_z; break;
    case g_ID(qG_z, qG_r): r->t = qG_z; break;
    case g_ID(qG_z, qG_c): r->t = qG_z; break;
    case g_ID(qG_p, qG_z): r->t = qG_z; break;
    case g_ID(qG_p, qG_p): r->t = qG_p; break;
    case g_ID(qG_p, qG_m): r->t = qG_m; break;
    case g_ID(qG_p, qG_r): r->r = b->r; r->t = qG_r; break;
    case g_ID(qG_p, qG_c): QLA_c_eq_c(r->c, b->c); r->t = qG_c; break;
    case g_ID(qG_m, qG_z): r->t = qG_z; break;
    case g_ID(qG_m, qG_p): r->t = qG_m; break;
    case g_ID(qG_m, qG_m): r->t = qG_p; break;
    case g_ID(qG_m, qG_r): r->r = -b->r; r->t = qG_r; break;
    case g_ID(qG_m, qG_c): QLA_c_eqm_c(r->c, b->c); r->t = qG_c; break;
    case g_ID(qG_r, qG_z): r->t = qG_z; break;
    case g_ID(qG_r, qG_p): r->r = a->r; r->t = qG_r; break;
    case g_ID(qG_r, qG_m): r->r = -a->r; r->t = qG_r; break;
    case g_ID(qG_r, qG_r): r->r = a->r * b->r; r->t = qG_r; break;
    case g_ID(qG_r, qG_c): mul_r(r, b, a->r); r->t = qG_c; break;
    case g_ID(qG_c, qG_z): r->t = qG_z; break;
    case g_ID(qG_c, qG_p): QLA_c_eq_c(r->c, a->c); r->t = qG_c; break;
    case g_ID(qG_c, qG_m): QLA_c_eqm_c(r->c, a->c); r->t = qG_c; break;
    case g_ID(qG_c, qG_r): mul_r(r, a, b->r); r->t = qG_c; break;
    case g_ID(qG_c, qG_c): mul_c(r, a, b); r->t = qG_c; break;
    }
}

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
    
    luaL_argcheck(L, v != 0, idx, "gamma expected");
    
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

    if (fst)
        luaL_addstring(&b, "0]");
    else
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
    int i;

    for (i = 0; i < 16; i++) {
        g_norm(&x->g[i]);
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
        g_add(&r->g[i], &x->g[i], &y->g[i]);

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
    g_add(&r->g[0], &x->g[0], &v);
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
    g_add(&r->g[0], &x->g[0], &v);
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
    g_add(&r->g[0], &x->g[0], &v);
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
    g_add(&r->g[0], &x->g[0], &v);
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
        g_add(&r->g[i], &x->g[i], &t);
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
    g_add(&r->g[0], &t, &v);
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
    g_add(&r->g[0], &x->g[0], &v);
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
    g_add(&r->g[0], &t, &v);
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
    g_add(&r->g[0], &x->g[0], &v);
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
            g_mul(&t, &x->g[i], &y->g[j]);
            if (s) {
                g_neg(&t, &t);
            }
            g_add(&r->g[k], &r->g[k], &t);
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
    QLA_real(v.c) = 0;
    QLA_imag(v.c) = 0;
    for (i = 0; i < 16; i++)
        g_mul(&r->g[i], &v, &b->g[i]);
    
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
    QLA_real(v.c) = 0;
    QLA_imag(v.c) = 0;
    for (i = 0; i < 16; i++)
        g_mul(&r->g[i], &b->g[i], &v);
    
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
    v.r = 0;
    QLA_c_eq_c(v.c, *a);
    for (i = 0; i < 16; i++)
        g_mul(&r->g[i], &v, &b->g[i]);
    
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
    v.r = 0;
    QLA_c_eq_c(v.c, *a);
    for (i = 0; i < 16; i++)
        g_mul(&r->g[i], &b->g[i], &v);
    
    return c_norm(L, r);
}

static int
q_g_div_r(lua_State *L)
{
    mClifford *b = qlua_checkClifford(L, 1);
    QLA_Real a = luaL_checknumber(L, 2);
    mClifford *r = qlua_newClifford(L);
    int i;
    mGamma v;

    v.t = qG_r;
    v.r = 1/a;
    QLA_real(v.c) = 0;
    QLA_imag(v.c) = 0;
    for (i = 0; i < 16; i++)
        g_mul(&r->g[i], &b->g[i], &v);
    
    return c_norm(L, r);
}

static int
q_g_div_c(lua_State *L)
{
    mClifford *b = qlua_checkClifford(L, 1);
    QLA_Complex *a = qlua_checkComplex(L, 2);
    mClifford *r = qlua_newClifford(L);
    int i;
    mGamma v;
    double n = 1 / (QLA_real(*a) * QLA_real(*a) + QLA_imag(*a) * QLA_imag(*a));

    v.t = qG_c;
    v.r = 0;
    QLA_real(v.c) = n * QLA_real(*a);
    QLA_imag(v.c) = -n * QLA_imag(*a);
    for (i = 0; i < 16; i++)
        g_mul(&r->g[i], &b->g[i], &v);
    
    return c_norm(L, r);
}

static int
q_g_mul_D(lua_State *L)
{
    mClifford *m = qlua_checkClifford(L, 1);
    mLatDirFerm *f = qlua_checkLatDirFerm(L, 2);
    mLatDirFerm *mf = qlua_newLatDirFerm(L);
    mLatDirFerm *r = qlua_newLatDirFerm(L);
    int i;

    CALL_QDP(L);
    QDP_D_eq_zero(r->ptr, *qCurrent);
    for (i = 0; i < 16; i++) {
        switch (m->g[i].t) {
        case qG_z: continue;
        case qG_p:
            QDP_D_eq_gamma_times_D(mf->ptr, f->ptr, i, *qCurrent);
            QDP_D_peq_D(r->ptr, mf->ptr, *qCurrent);
            break;
        case qG_m:
            QDP_D_eq_gamma_times_D(mf->ptr, f->ptr, i, *qCurrent);
            QDP_D_meq_D(r->ptr, mf->ptr, *qCurrent);
            break;
        case qG_r:
            QDP_D_eq_gamma_times_D(mf->ptr, f->ptr, i, *qCurrent);
            QDP_D_peq_r_times_D(r->ptr, &m->g[i].r, mf->ptr, *qCurrent);
            break;
        case qG_c:
            QDP_D_eq_gamma_times_D(mf->ptr, f->ptr, i, *qCurrent);
            QDP_D_peq_c_times_D(r->ptr, &m->g[i].c, mf->ptr, *qCurrent);
            break;
        }
    }
    return 1;
}

static int
q_g_mul_P(lua_State *L)
{
    mClifford *m = qlua_checkClifford(L, 1);
    mLatDirProp *f = qlua_checkLatDirProp(L, 2);
    mLatDirProp *mf = qlua_newLatDirProp(L);
    mLatDirProp *r = qlua_newLatDirProp(L);
    int i;

    CALL_QDP(L);
    QDP_P_eq_zero(r->ptr, *qCurrent);
    for (i = 0; i < 16; i++) {
        switch (m->g[i].t) {
        case qG_z: continue;
        case qG_p:
            QDP_P_eq_gamma_times_P(mf->ptr, f->ptr, i, *qCurrent);
            QDP_P_peq_P(r->ptr, mf->ptr, *qCurrent);
            break;
        case qG_m:
            QDP_P_eq_gamma_times_P(mf->ptr, f->ptr, i, *qCurrent);
            QDP_P_meq_P(r->ptr, mf->ptr, *qCurrent);
            break;
        case qG_r:
            QDP_P_eq_gamma_times_P(mf->ptr, f->ptr, i, *qCurrent);
            QDP_P_peq_r_times_P(r->ptr, &m->g[i].r, mf->ptr, *qCurrent);
            break;
        case qG_c:
            QDP_P_eq_gamma_times_P(mf->ptr, f->ptr, i, *qCurrent);
            QDP_P_peq_c_times_P(r->ptr, &m->g[i].c, mf->ptr, *qCurrent);
            break;
        }
    }
    return 1;
}

static int
q_P_mul_g(lua_State *L)
{
    mLatDirProp *f = qlua_checkLatDirProp(L, 1);
    mClifford *m = qlua_checkClifford(L, 2);
    mLatDirProp *mf = qlua_newLatDirProp(L);
    mLatDirProp *r = qlua_newLatDirProp(L);
    int i;

    CALL_QDP(L);
    QDP_P_eq_zero(r->ptr, *qCurrent);
    for (i = 0; i < 16; i++) {
        switch (m->g[i].t) {
        case qG_z: continue;
        case qG_p:
            QDP_P_eq_P_times_gamma(mf->ptr, f->ptr, i, *qCurrent);
            QDP_P_peq_P(r->ptr, mf->ptr, *qCurrent);
            break;
        case qG_m:
            QDP_P_eq_P_times_gamma(mf->ptr, f->ptr, i, *qCurrent);
            QDP_P_meq_P(r->ptr, mf->ptr, *qCurrent);
            break;
        case qG_r:
            QDP_P_eq_P_times_gamma(mf->ptr, f->ptr, i, *qCurrent);
            QDP_P_peq_r_times_P(r->ptr, &m->g[i].r, mf->ptr, *qCurrent);
            break;
        case qG_c:
            QDP_P_eq_P_times_gamma(mf->ptr, f->ptr, i, *qCurrent);
            QDP_P_peq_c_times_P(r->ptr, &m->g[i].c, mf->ptr, *qCurrent);
            break;
        }
    }
    return 1;
}

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
q_g_conj(lua_State *L)
{
    mClifford *x = qlua_checkClifford(L, 1);
    mClifford *r = qlua_newClifford(L);
    int i;

    for (i = 0; i < 16; i++) {
        if (gconj[i])
            g_neg(&r->g[i], &x->g[i]);
        else
            r->g[i] = x->g[i];
        if (r->g[i].t == qG_c)
            QLA_c_eq_ca(r->g[i].c, r->g[i].c);
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

            return 1;
        }
        break;
    }
    }
    return qlua_badconstr(L, "gamma");
}

typedef void X_project(lua_State *L,
                       QLA_DiracFermion *r[QDP_Nc][QDP_Ns/2],
                       int mu, int sign,
                       QLA_DiracPropagator *f);

static void
q_P_left_proj(lua_State *L,
              QLA_DiracFermion *r[QDP_Nc][QDP_Ns/2],
              int mu,
              int sign,
              QLA_DiracPropagator *f)
{
    int count = QDP_sites_on_node;
    int k, ic, is;

    for (k = 0; k < count; k++) {
        for (ic = 0; ic < QDP_Nc; ic++) {
            for (is = 0; is < QDP_Ns; is++) {
                int jc, js;
                QLA_DiracFermion fk;
                QLA_HalfFermion hk;
                for (jc = 0; jc < QDP_Nc; jc++) {
                    for (js = 0; js < QDP_Ns; js++) {
                        QLA_c_eq_c(QLA_elem_D(fk, jc, js),
                                   QLA_elem_P(*f, jc, js, ic, is));
                    }
                }
                QLA_H_eq_spproj_D(&hk, &fk, mu, sign);
                for (jc = 0; jc < QDP_Nc; jc++) {
                    for (js = 0; js < QDP_Ns/2; js++) {
                        QLA_c_eq_c(QLA_elem_D(*r[jc][js], ic, is),
                                   QLA_elem_H(hk, jc, js));
                    }
                }
            }
        }
        f++;
        for (ic = 0; ic < QDP_Nc; ic++)
            for (is = 0; is < QDP_Ns/2; is++)
                r[ic][is]++;
    }
}

static void
q_P_right_proj(lua_State *L,
               QLA_DiracFermion *r[QDP_Nc][QDP_Ns/2],
               int mu,
               int sign,
               QLA_DiracPropagator *f)
{
    int count = QDP_sites_on_node;
    int k, ic, is;

    if (mu == 0 || mu == 2) sign = 1 - sign; /* AAA gamma basis dependent */

    for (k = 0; k < count; k++) {
        for (ic = 0; ic < QDP_Nc; ic++) {
            for (is = 0; is < QDP_Ns; is++) {
                int jc, js;
                QLA_DiracFermion fk;
                QLA_HalfFermion hk;
                for (jc = 0; jc < QDP_Nc; jc++) {
                    for (js = 0; js < QDP_Ns; js++) {
                        QLA_c_eq_c(QLA_elem_D(fk, jc, js),
                                   QLA_elem_P(*f, ic, is, jc, js));
                    }
                }
                QLA_H_eq_spproj_D(&hk, &fk, mu, sign);
                for (jc = 0; jc < QDP_Nc; jc++) {
                    for (js = 0; js < QDP_Ns/2; js++) {
                        QLA_c_eq_c(QLA_elem_D(*r[jc][js], ic, is),
                                   QLA_elem_H(hk, jc, js));
                    }
                }
            }
        }
        f++;
        for (ic = 0; ic < QDP_Nc; ic++)
            for (is = 0; is < QDP_Ns/2; is++)
                r[ic][is]++;
    }
}

static int
q_P_project(lua_State *L, X_project *op)
{
    const char *sign = luaL_checkstring(L, 1);
    int mu = qlua_checkgammaindex(L, 2);
    mLatDirProp *f = qlua_checkLatDirProp(L, 3);
    int isign = 0;
    int ic, is;
    QLA_DiracPropagator *ff;
    mLatDirFerm *r[QDP_Nc][QDP_Ns/2];
    QLA_DiracFermion *rr[QDP_Nc][QDP_Ns/2];

    if (strcmp(sign, "plus") == 0)
        isign = 1; /* XXX check it */
    else if (strcmp(sign, "minus") == 0)
        isign = 0; /* XXX check it */
    else
        return luaL_error(L, "bad sign parameter");
    
    lua_createtable(L, QDP_Nc, 0);
    for (ic = 0; ic < QDP_Nc; ic++) {
        lua_createtable(L, QDP_Ns/2, 0);
        for (is = 0; is < QDP_Ns/2; is++) {
            r[ic][is] = qlua_newLatDirFerm(L);
            lua_rawseti(L, -2, is + 1);
        }
        lua_rawseti(L, -2, ic + 1);
    }

    if (mu == 5) mu = 4; /* [sic] -- funny numbering in QLA's spproj/sprecon */

    CALL_QDP(L);
    for (ic = 0; ic < QDP_Nc; ic++)
        for (is = 0; is < QDP_Ns/2; is++) 
            rr[ic][is] = QDP_expose_D(r[ic][is]->ptr);
    ff = QDP_expose_P(f->ptr);

    op(L, rr, mu, isign, ff);

    QDP_reset_P(f->ptr);
    for (ic = 0; ic < QDP_Nc; ic++)
        for (is = 0; is < QDP_Ns/2; is++) 
            QDP_reset_D(r[ic][is]->ptr);
        
    return 1;
}

static int
q_left_project(lua_State *L)
{
    return q_P_project(L, q_P_left_proj);
}

static int
q_right_project(lua_State *L)
{
    return q_P_project(L, q_P_right_proj);
}

typedef void X_reconstruct(lua_State *L,
                           QLA_DiracPropagator *r,
                           int mu, int sign,
                           QLA_DiracFermion *a[QDP_Nc][QDP_Ns/2]);

static void
q_P_left_recon(lua_State *L,
               QLA_DiracPropagator *r,
               int mu, int sign,
               QLA_DiracFermion *a[QDP_Nc][QDP_Ns/2])
{
    int count = QDP_sites_on_node;
    int k, ic, is;

    for (k = 0; k < count; k++) {
        for (ic = 0; ic < QDP_Nc; ic++) {
            for (is = 0; is < QDP_Ns; is++) {
                int jc, js;
                QLA_DiracFermion fk;
                QLA_HalfFermion hk;
                for (jc = 0; jc < QDP_Nc; jc++) {
                    for (js = 0; js < QDP_Ns/2; js++) {
                        QLA_c_eq_c(QLA_elem_H(hk, jc, js),
                                   QLA_elem_D(*a[jc][js], ic, is));
                    }
                }
                QLA_D_eq_sprecon_H(&fk, &hk, mu, sign);
                for (jc = 0; jc < QDP_Nc; jc++) {
                    for (js = 0; js < QDP_Ns; js++) {
                        QLA_c_eq_c(QLA_elem_P(*r, jc, js, ic, is),
                                   QLA_elem_D(fk, jc, js));
                    }
                }
            }
        }
        r++;
        for (ic = 0; ic < QDP_Nc; ic++)
            for (is = 0; is < QDP_Ns/2; is++)
                a[ic][is]++;
    }
    
}

static void
q_P_right_recon(lua_State *L,
                QLA_DiracPropagator *r,
                int mu, int sign,
                QLA_DiracFermion *a[QDP_Nc][QDP_Ns/2])
{
    int count = QDP_sites_on_node;
    int k, ic, is;

    if (mu == 0 || mu == 2) sign = 1 - sign; /* AAA gamma basis dependent */

    for (k = 0; k < count; k++) {
        for (ic = 0; ic < QDP_Nc; ic++) {
            for (is = 0; is < QDP_Ns; is++) {
                int jc, js;
                QLA_DiracFermion fk;
                QLA_HalfFermion hk;
                for (jc = 0; jc < QDP_Nc; jc++) {
                    for (js = 0; js < QDP_Ns/2; js++) {
                        QLA_c_eq_c(QLA_elem_H(hk, jc, js),
                                   QLA_elem_D(*a[jc][js], ic, is));
                    }
                }
                QLA_D_eq_sprecon_H(&fk, &hk, mu, sign);
                for (jc = 0; jc < QDP_Nc; jc++) {
                    for (js = 0; js < QDP_Ns; js++) {
                        QLA_c_eq_c(QLA_elem_P(*r, ic, is, jc, js),
                                   QLA_elem_D(fk, jc, js));
                    }
                }
            }
        }
        r++;
        for (ic = 0; ic < QDP_Nc; ic++)
            for (is = 0; is < QDP_Ns/2; is++)
                a[ic][is]++;
    }
    
}

static int
q_P_reconstruct(lua_State *L, X_reconstruct *op)
{
    const char *sign = luaL_checkstring(L, 1);
    int mu = qlua_checkgammaindex(L, 2);
    mLatDirProp *r = qlua_newLatDirProp(L);
    int isign = 0;
    int ic, is;
    mLatDirFerm *a[QDP_Nc][QDP_Ns/2];
    QLA_DiracFermion *aa[QDP_Nc][QDP_Ns/2];
    QLA_DiracPropagator *rr;

    luaL_checktype(L, 3, LUA_TTABLE);
    for (ic = 0; ic < QDP_Nc; ic++) {
        lua_pushnumber(L, ic + 1);
        lua_gettable(L, 3);
        luaL_checktype(L, -1, LUA_TTABLE);
        for (is = 0; is < QDP_Ns/2; is++) {
            lua_pushnumber(L, is + 1);
            lua_gettable(L, -2);
            a[ic][is] = qlua_checkLatDirFerm(L, -1);
            lua_pop(L, 1);
        }
        lua_pop(L, 1);
    }

    if (strcmp(sign, "plus") == 0)
        isign = 1; /* XXX check it */
    else if (strcmp(sign, "minus") == 0)
        isign = 0; /* XXX check it */
    else
        return luaL_error(L, "bad sign parameter");

    if (mu == 5) mu = 4; /* [sic] -- funny numbering in QLA's spproj/sprecon */

    CALL_QDP(L);
    /* ZZZ will break if any of a[][] is aliased */
    for (ic = 0; ic < QDP_Nc; ic++)
        for (is = 0; is < QDP_Ns/2; is++)
            aa[ic][is] = QDP_expose_D(a[ic][is]->ptr);
    rr = QDP_expose_P(r->ptr);

    op(L, rr, mu, isign, aa);

    QDP_reset_P(r->ptr);
    for (ic = 0; ic < QDP_Nc; ic++)
        for (is = 0; is < QDP_Ns/2; is++)
            QDP_reset_D(a[ic][is]->ptr);

    return 1;
}

static int
q_left_reconstruct(lua_State *L)
{
    return q_P_reconstruct(L, q_P_left_recon);
}

static int
q_right_reconstruct(lua_State *L)
{
    return q_P_reconstruct(L, q_P_right_recon);
}

static struct luaL_Reg mtGamma[] = {
    { "__tostring",        q_g_fmt },
    { "__unm",             q_g_neg },
    { "__add",             qlua_add },
    { "__sub",             qlua_sub },
    { "__mul",             qlua_mul },
    { "__div",             qlua_div },
    { "conj",              q_g_conj },
    { NULL,                NULL }
};

static struct luaL_Reg fGamma[] = {
    { "gamma",              q_gamma },
    { NULL,                 NULL }
};

static struct luaL_Reg fProj[] = {
    { "left_project",       q_left_project },
    { "left_reconstruct",   q_left_reconstruct },
    { "right_project",      q_right_project },
    { "right_reconstruct",  q_right_reconstruct },
    { NULL,                 NULL }
};

int
init_gamma(lua_State *L)
{
    luaL_register(L, qcdlib, fProj);
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

