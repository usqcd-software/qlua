#include "modules.h"                                                 /* DEPS */
#include "qlua.h"                                                    /* DEPS */
#include "qcomplex.h"                                                /* DEPS */
#ifdef HAS_GSL
#include "qmatrix.h"                                                 /* DEPS */
#endif
#include "qgamma.h"                                                  /* DEPS */
#include "lattice.h"                                                 /* DEPS */
#include "latdirferm.h"                                              /* DEPS */
#include "latdirprop.h"                                              /* DEPS */
#include "seqdirferm.h"                                              /* DEPS */
#include "seqdirprop.h"                                              /* DEPS */
#include <string.h>

static const char mtnGamma[] = "qlua.mtGamma";

#define Gi(a,b)  ((a)*qG_t+(b))

/* generated from QDP gamma operations */
static char gconj[] =  {0, 1, 0, 1, 1, 0, 1, 0, 0, 1, 0, 1, 1, 0, 1, 0};
static char gtrans[] = {0, 1, 0, 0, 1, 1, 0, 1, 0, 0, 1, 0, 0, 1, 0, 0};
static char gadj[] =   {0, 0, 0, 1, 0, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 0};

static const char *gn[] = {
    NULL, "g0",  "g1",  "g01",  "g2",  "g02", "g12",   "g012",
    "g3", "g03", "g13", "g013", "g23", "g023", "g123", "g5"};

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

#ifdef HAS_GSL
/* 0,1: l, r
 * 2: im, im, re, re
 * 3: i, s
 */
static const signed char gmv[4][4][4][2] = {
	{{{ 3, -1}, {12,  1}, { 0,  1}, {15,  1}},
	 {{ 6, -1}, { 9,  1}, { 5,  1}, {10,  1}},
	 {{ 4, -1}, {11, -1}, { 7, -1}, { 8,  1}},
	 {{ 1, -1}, {14, -1}, { 2, -1}, {13,  1}}},
	{{{ 6, -1}, { 9,  1}, { 5, -1}, {10, -1}},
	 {{ 3,  1}, {12, -1}, { 0,  1}, {15,  1}},
	 {{ 1, -1}, {14, -1}, { 2,  1}, {13, -1}},
	 {{ 4,  1}, {11,  1}, { 7, -1}, { 8,  1}}},
	{{{ 4,  1}, {11, -1}, { 7,  1}, { 8,  1}},
	 {{ 1,  1}, {14, -1}, { 2,  1}, {13,  1}},
	 {{ 3, -1}, {12, -1}, { 0,  1}, {15, -1}},
	 {{ 6, -1}, { 9, -1}, { 5,  1}, {10, -1}}},
	{{{ 1,  1}, {14, -1}, { 2, -1}, {13, -1}},
	 {{ 4, -1}, {11,  1}, { 7,  1}, { 8,  1}},
	 {{ 6, -1}, { 9, -1}, { 5, -1}, {10,  1}},
	 {{ 3,  1}, {12,  1}, { 0,  1}, {15, -1}}}};
#endif

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
            sprintf(fmt, "%scomplex(%g,%g)*", fst ? "" : "+",
                    QLA_real(v->g[i].c), QLA_imag(v->g[i].c));
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
    QLA_D_Real b = luaL_checknumber(L, 1);
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
    QLA_D_Real b = luaL_checknumber(L, 2);
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
    QLA_D_Complex *b = qlua_checkComplex(L, 1);
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
    QLA_D_Complex *b = qlua_checkComplex(L, 2);
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
    QLA_D_Real b = luaL_checknumber(L, 1);
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
    QLA_D_Real b = luaL_checknumber(L, 2);
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
    QLA_D_Complex *b = qlua_checkComplex(L, 1);
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
    QLA_D_Complex *b = qlua_checkComplex(L, 2);
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
    QLA_D_Real a = luaL_checknumber(L, 1);
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
    QLA_D_Real a = luaL_checknumber(L, 2);
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
    QLA_D_Complex *a = qlua_checkComplex(L, 1);
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
    QLA_D_Complex *a = qlua_checkComplex(L, 2);
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
    QLA_D_Real a = luaL_checknumber(L, 2);
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
    QLA_D_Complex *a = qlua_checkComplex(L, 2);
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
q_g_adjoin(lua_State *L)
{
    mClifford *x = qlua_checkClifford(L, 1);
    mClifford *r = qlua_newClifford(L);
    int i;

    for (i = 0; i < 16; i++) {
        if (gadj[i])
            g_neg(&r->g[i], &x->g[i]);
        else
            r->g[i] = x->g[i];
        if (r->g[i].t == qG_c)
            QLA_c_eq_ca(r->g[i].c, r->g[i].c);
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
q_g_transpose(lua_State *L)
{
    mClifford *x = qlua_checkClifford(L, 1);
    mClifford *r = qlua_newClifford(L);
    int i;

    for (i = 0; i < 16; i++) {
        if (gtrans[i])
            g_neg(&r->g[i], &x->g[i]);
        else
            r->g[i] = x->g[i];
    }

    return 1;
}

#ifdef HAS_GSL
static gsl_complex
get_gmv(int i, int j, const mClifford *g, int k)
{
	int g_i = gmv[i][j][k][0];
	int g_s = gmv[i][j][k][1];
	gsl_complex r;
	GSL_REAL(r) = 0;
	GSL_IMAG(r) = 0;
	switch (g->g[g_i].t) {
	case qG_z:
		break;
	case qG_p:
		GSL_REAL(r) = g_s;
		GSL_IMAG(r) = 0;
		break;
	case qG_m:
		GSL_REAL(r) = -g_s;
		GSL_IMAG(r) = 0;
		break;
	case qG_r:
		GSL_REAL(r) = g_s * g->g[g_i].r;
		GSL_IMAG(r) = 0;
		break;
	case qG_c:
		GSL_REAL(r) = QLA_real(g->g[g_i].c) * g_s;
		GSL_IMAG(r) = QLA_imag(g->g[g_i].c) * g_s;
		break;
	}
	return r;
}

mMatComplex *
gamma2matrix(lua_State *L, int idx)
{
	mClifford *g = qlua_checkClifford(L, idx);
	mMatComplex *m = qlua_newMatComplex(L, 4, 4);
	int i, j;

	for (i = 0; i < 4; i++) {
		for (j = 0; j < 4; j++) {
			gsl_complex xi0 = get_gmv(i, j, g, 0);
			gsl_complex xi1 = get_gmv(i, j, g, 1);
			gsl_complex xr0 = get_gmv(i, j, g, 2);
			gsl_complex xr1 = get_gmv(i, j, g, 3);
			gsl_complex v;
			GSL_REAL(v) = GSL_REAL(xr0) + GSL_REAL(xr1) - GSL_IMAG(xi0) - GSL_IMAG(xi1);
			GSL_IMAG(v) = GSL_IMAG(xr0) + GSL_IMAG(xr1) + GSL_REAL(xi0) + GSL_REAL(xi1);
			gsl_matrix_complex_set(m->m, j, i, v);
		}
	}
	return m;
}

static int
q_g_matrix(lua_State *L)
{
	gamma2matrix(L, 1);
	return 1;
}
#endif

static int
q_gamma(lua_State *L)
{
    switch (qlua_qtype(L, 1)) {
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
    default:
        break;
    }
    return qlua_badconstr(L, "gamma");
}

#if USE_Nc2
#define QNc  '2'
#define Qcolors "2"
#define Qs(a)   a ## 2
#define Qx(a,b)  a ## 2 ## b
#define QC(x)    2
#define QNC(x)
#include "qgamma-x.c"                                               /* DEPS */
#endif

#if USE_Nc3
#define QNc  '3'
#define Qcolors "3"
#define Qs(a)   a ## 3
#define Qx(a,b)  a ## 3 ## b
#define QC(x)    3
#define QNC(x)
#include "qgamma-x.c"                                               /* DEPS */
#endif

#if USE_NcN
#define QNc  'N'
#define Qcolors "N"
#define Qs(a)   a ## N
#define Qx(a,b)  a ## N ## b
#define QC(x)    (x)->nc
#define QNC(x)   (x), 
#include "qgamma-x.c"                                               /* DEPS */
#endif

static int
q_left_project(lua_State *L)
{
    switch (qlua_qtype(L, 3)) {
#if USE_Nc2
    case qLatDirProp2: return q_P_project2(L, q_P_left_proj2);
#endif
#if USE_Nc3
    case qLatDirProp3: return q_P_project3(L, q_P_left_proj3);
#endif
#if USE_NcN
    case qLatDirPropN: return q_P_projectN(L, q_P_left_projN);
#endif
    default:
        return luaL_error(L, "unexpected type in g:left_project()");
    }
}

static int
q_right_project(lua_State *L)
{
    switch (qlua_qtype(L, 3)) {
#if USE_Nc2
    case qLatDirProp2: return q_P_project2(L, q_P_right_proj2);
#endif
#if USE_Nc3
    case qLatDirProp3: return q_P_project3(L, q_P_right_proj3);
#endif
#if USE_NcN
    case qLatDirPropN: return q_P_projectN(L, q_P_right_projN);
#endif
    default:
        return luaL_error(L, "unexpected type in g:right_project()");
    }
}

static int
q_left_reconstruct(lua_State *L)
{
    int n = lua_objlen(L, 3);
    if (n < 2)
        return luaL_error(L, "bad arg size for g:left_reconstruct");

    switch (n) {
#if USE_Nc2
    case 2: return q_P_reconstruct2(L, q_P_left_recon2, n);
#endif
#if USE_Nc3
    case 3: return q_P_reconstruct3(L, q_P_left_recon3, n);
#endif
#if USE_NcN
    default: return q_P_reconstructN(L, q_P_left_reconN, n);
#else
    default: return luaL_error(L, "bad number of colors");
#endif
    }
}

static int
q_right_reconstruct(lua_State *L)
{
    int n = lua_objlen(L, 3);
    if (n < 2)
        return luaL_error(L, "bad arg size for g:right_reconstruct");

    switch (n) {
#if USE_Nc2
    case 2: return q_P_reconstruct2(L, q_P_right_recon2, n);
#endif
#if USE_Nc3
    case 3: return q_P_reconstruct3(L, q_P_right_recon3, n);
#endif
#if USE_NcN
    default: return q_P_reconstructN(L, q_P_right_reconN, n);
#else
    default: return luaL_error(L, "bad number of colors");
#endif
    }
}

static struct luaL_Reg mtGamma[] = {
    { "__tostring",        q_g_fmt },
    { "__unm",             q_g_neg },
    { "__add",             qlua_add },
    { "__sub",             qlua_sub },
    { "__mul",             qlua_mul },
    { "__div",             qlua_div },
    { "adjoin",            q_g_adjoin },
    { "conj",              q_g_conj },
    { "transpose",         q_g_transpose },
#ifdef HAS_GSL
	{ "matrix",            q_g_matrix },
#endif
    { NULL,                NULL }
};

static struct luaL_Reg fGamma[] = {
    { "gamma",              q_gamma },
    { NULL,                 NULL }
};

static struct luaL_Reg fProj[] = {
    { "left_project",       q_left_project },
    { "right_project",      q_right_project },
    { "left_reconstruct",   q_left_reconstruct },
    { "right_reconstruct",  q_right_reconstruct },
    { NULL,                 NULL }
};

int
init_gamma(lua_State *L)
{
    static const QLUA_Op2 ops[] = {
        { qlua_add_table, qGamma,        qGamma,        q_g_add_g    },
        { qlua_add_table, qReal,         qGamma,        q_r_add_g    },
        { qlua_add_table, qGamma,        qReal,         q_g_add_r    },
        { qlua_add_table, qComplex,      qGamma,        q_c_add_g    },
        { qlua_add_table, qGamma,        qComplex,      q_g_add_c    },
        { qlua_sub_table, qGamma,        qGamma,        q_g_sub_g    },
        { qlua_sub_table, qReal,         qGamma,        q_r_sub_g    },
        { qlua_sub_table, qGamma,        qReal,         q_g_sub_r    },
        { qlua_sub_table, qComplex,      qGamma,        q_c_sub_g    },
        { qlua_sub_table, qGamma,        qComplex,      q_g_sub_c    },
        { qlua_mul_table, qGamma,        qGamma,        q_g_mul_g    },
        { qlua_mul_table, qReal,         qGamma,        q_r_mul_g    },
        { qlua_mul_table, qGamma,        qReal,         q_g_mul_r    },
        { qlua_mul_table, qComplex,      qGamma,        q_c_mul_g    },
        { qlua_mul_table, qGamma,        qComplex,      q_g_mul_c    },
        { qlua_div_table, qGamma,        qReal,         q_g_div_r    },
        { qlua_div_table, qGamma,        qComplex,      q_g_div_c    },
        { NULL,           qNoType,       qNoType,       NULL         }
    };
    luaL_register(L, qcdlib, fProj);
    lua_getglobal(L, "_G");
    luaL_register(L, NULL, fGamma);
    qlua_metatable(L, mtnGamma, mtGamma, qGamma);
    qlua_reg_op2(ops);
#if USE_Nc2
    qlua_reg_op2(ops2);
#endif
#if USE_Nc3
    qlua_reg_op2(ops3);
#endif
#if USE_NcN
    qlua_reg_op2(opsN);
#endif

    return 0;
}

int
fini_gamma(lua_State *L)
{
    return 0;
}
