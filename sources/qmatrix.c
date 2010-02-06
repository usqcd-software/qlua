/* Rings over matrices, real and complex, not mixed.
 * Modules are not implemented because there are too many of them
 */
#include <qlua.h>                                                    /* DEPS */
#include <qcomplex.h>                                                /* DEPS */
#include <qmatrix.h>                                                 /* DEPS */
#include <qvector.h>                                                 /* DEPS */
#include <matrix.h>                                                  /* DEPS */
#include <string.h>
#include <math.h>

const char mtnMatReal[]           = "matrix.mtReal";
const char mtnMatComplex[]        = "matrix.mtComplex";
static const char opMatReal[]     = "matrix.mtReal.ops";
static const char opMatComplex[]  = "matrix.mtComplex.ops";

static char matrix_ns[] = "matrix";

#define MSIZE(s,a,b,t) (sizeof (s) + ((a)*(b)-1)*sizeof (t))
#define RM(a,i,j)  ((a)->val + (i) + (j) * (a)->size_l)
#define CM(a,i,j)  ((a)->val + 2 * ((i) + (j) * (a)->size_l))

static void *
qlua_checkMatrix(lua_State *L, int idx, const char *mt, const char *name)
{
    void *v = luaL_checkudata(L, idx, mt);

    if (v == 0) {
        char fmt[72];
        sprintf(fmt, "matrix.%s expected", name);
        luaL_argcheck(L, v != NULL, idx, name);
    }

    return v;
}

static void *
qlua_newMatrix(lua_State *L, int s, const char *mt)
{
    void *v = lua_newuserdata(L, s);

    luaL_getmetatable(L, mt);
    lua_setmetatable(L, -2);
    
    return v;
}

/* real matrices */
mMatReal *
qlua_checkMatReal(lua_State *L, int idx)
{
    return (mMatReal *)qlua_checkMatrix(L, idx, mtnMatReal, "Real");
}

mMatReal *
qlua_newMatReal(lua_State *L, int sl, int sr)
{
    return qlua_newMatrix(L, MSIZE(mMatReal, sl, sr, double), mtnMatReal);
}

static int
md_fmt(lua_State *L)                                           /* (-1,+1,e) */
{
    char fmt[72];
    mMatReal *v = qlua_checkMatReal(L, 1);

    sprintf(fmt, "matrix.real[{%d,%d}]", v->size_l, v->size_r);
    lua_pushstring(L, fmt);

    return 1;
}

static int
md_inverse(lua_State *L)
{
    mMatReal *m = qlua_checkMatReal(L, 1);
    mMatReal *t;
    mMatReal *r;

    if (m->size_l != m->size_r)
        return luaL_error(L, "square matrix expected");
    
    t = qlua_newMatReal(L, m->size_l, m->size_l);
    r = qlua_newMatReal(L, m->size_l, m->size_l);
    memcpy(t->val, m->val, m->size_l * m->size_l * sizeof (double));
    r->size_l = m->size_l;
    r->size_r = m->size_l;
    if (matrix_rinverse(m->size_l, t->val, r->val))
        return luaL_error(L, "inverting singular matrix");

    return 1;
}

static int
md_eigen(lua_State *L)                                         /* (-1,+2,e) */
{
    mMatReal *m = qlua_checkMatReal(L, 1);
    mVecReal *lambda;
    mMatReal *trans;

    if (m->size_l != m->size_r)
        return luaL_error(L, "matrix:eigen() expects square matrix");
    
    lambda = qlua_newVecReal(L, m->size_l);
    lambda->size = m->size_l;
    trans = qlua_newMatReal(L, m->size_l, m->size_r);
    trans->size_l = m->size_l;
    trans->size_r = m->size_r;

    matrix_reigenvec(L, m->size_l, m->val, lambda->val, trans->val);

    return 2;
}

static int
md_transpose(lua_State *L)
{
    mMatReal *a = qlua_checkMatReal(L, 1);
    mMatReal *b = qlua_newMatReal(L, a->size_r, a->size_l);
    int l, r;

    b->size_l = a->size_r;
    b->size_r = a->size_l;
    for (l = 0; l < b->size_l; l++) {
        for (r = 0; r < b->size_r; r++) {
            RM(b,l,r)[0] = RM(a,r,l)[0];
        }
    }
    return 1;
}

static int
md_trace(lua_State *L)                                         /* (-1,+2,e) */
{
    mMatReal *m = qlua_checkMatReal(L, 1);
    int i;
    double tr;
    
    if (m->size_l != m->size_r)
        return luaL_error(L, "matrix:trace() expects square matrix");
    
    for (tr = 0, i = 0; i < m->size_l; i++)
        tr += RM(m,i,i)[0];
    lua_pushnumber(L, tr);

    return 1;
}

static int
md_dims(lua_State *L)                                          /* (-1,+2,e) */
{
    mMatReal *m = qlua_checkMatReal(L, 1);

    lua_pushnumber(L, m->size_l);
    lua_pushnumber(L, m->size_r);

    return 2;
}

static int
md_get(lua_State *L)                                           /* (-2,+1,e) */
{
    switch (qlua_gettype(L, 2)) {
    case qTable: {
        mMatReal *v = qlua_checkMatReal(L, 1);
        int sl, sr;

        qlua_checkindex2(L, 2, "matrix get", &sl, &sr);

        if ((sl >= 0) && (sl < v->size_l) &&
            (sr >= 0) && (sr < v->size_r)) {
            lua_pushnumber(L, RM(v, sl, sr)[0]);
            return 1;
        }
        break;
    }
    case qString:
        return qlua_lookup(L, 2, opMatReal);
    }
    return qlua_badindex(L, "matrix.real[]");
}


static int
md_put(lua_State *L)                                           /* (-3,+0,e) */
{
    mMatReal *v = qlua_checkMatReal(L, 1);
    int sl, sr;
    double x = luaL_checknumber(L, 3);

    qlua_checkindex2(L, 2, "matrix set", &sl, &sr);

    if ((sl >= 0) && (sl < v->size_l) &&
        (sr >= 0) && (sr < v->size_r)) {
        RM(v, sl, sr)[0] = x;
        return 0;
    }

    return qlua_badindex(L, "matrix.real[]");
}

static int
rm_add_rm(lua_State *L)
{
    mMatReal *a = qlua_checkMatReal(L, 1);
    mMatReal *b = qlua_checkMatReal(L, 2);
    int al = a->size_l;
    int ar = a->size_r;
    mMatReal *r = qlua_newMatReal(L, al, ar);
    int i;

    if ((al != b->size_l) || (ar != b->size_r))
        return luaL_error(L, "matrix sizes mismatch in m + m");

    r->size_l = al;
    r->size_r = ar;
    for (i = 0; i < al * ar; i++)
        r->val[i] = a->val[i] + b->val[i];

    return 1;
}

static int
rm_sub_rm(lua_State *L)
{
    mMatReal *a = qlua_checkMatReal(L, 1);
    mMatReal *b = qlua_checkMatReal(L, 2);
    int al = a->size_l;
    int ar = a->size_r;
    mMatReal *r = qlua_newMatReal(L, al, ar);
    int i;

    if ((al != b->size_l) || (ar != b->size_r))
        return luaL_error(L, "matrix sizes mismatch in m - m");

    r->size_l = al;
    r->size_r = ar;
    for (i = 0; i < al * ar; i++)
        r->val[i] = a->val[i] - b->val[i];

    return 1;
}

static int
rm_mul_rm(lua_State *L)
{
    mMatReal *a = qlua_checkMatReal(L, 1);
    mMatReal *b = qlua_checkMatReal(L, 2);
    int al = a->size_l;
    int ar = a->size_r;
    int bl = b->size_l;
    int br = b->size_r;
    mMatReal *r = qlua_newMatReal(L, al, br);
    int i, j, k;

    if (ar != bl)
        return luaL_error(L, "matrix sizes mismatch in m * m");

    r->size_l = al;
    r->size_r = br;
    for (i = 0; i < al; i++) {
        for (j = 0; j < br; j++) {
            double v = 0;
            for (k = 0; k < ar; k++)
                v += RM(a,i,k)[0] * RM(b, k, j)[0];
            RM(r, i, j)[0] = v;
        }
    }

    return 1;
}

static int
do_rrmul(lua_State *L, mMatReal *a, double b)
{
    mMatReal *r = qlua_newMatReal(L, a->size_l, a->size_r);
    int n = a->size_l * a->size_r;
    int i;

    r->size_l = a->size_l;
    r->size_r = a->size_r;

    for (i = 0; i < n; i++)
        r->val[i] = a->val[i] * b;

    return 1;
}

static int
rm_mul_r(lua_State *L)
{
    mMatReal *a = qlua_checkMatReal(L, 1);
    double b = luaL_checknumber(L, 2);

    return do_rrmul(L, a, b);
}

static int
r_mul_rm(lua_State *L)
{
    double b = luaL_checknumber(L, 1);
    mMatReal *a = qlua_checkMatReal(L, 2);

    return do_rrmul(L, a, b);
}

static int
rm_div_r(lua_State *L)
{
    mMatReal *a = qlua_checkMatReal(L, 1);
    double b = luaL_checknumber(L, 2);

    return do_rrmul(L, a, 1/b);
}

static int
m_real(lua_State *L)                                         /* (-1,+1,e) */
{
    int sl, sr;
    mMatReal *m;

    qlua_checkindex2(L, 1, "matrix size", &sl, &sr);

    m = qlua_newMatReal(L, sl, sr);
    m->size_l = sl;
    m->size_r = sr;
    memset(m->val, 0, sl * sr * sizeof (double));

    return 1;
}

static const luaL_Reg MatRealMethods[] = {
    { "dims",                 md_dims       },
    { "transpose",            md_transpose  },
    { "trace",                md_trace      },
    { "inverse",              md_inverse    },
    { "symmetric_eigen",      md_eigen      },
    { NULL, NULL}
};

static const luaL_Reg mtMatReal[] = {
    { "__tostring",     md_fmt    },
    { "__index",        md_get    },
    { "__newindex",     md_put    },
    { "__add",          qlua_add  },
    { "__sub",          qlua_sub  },
    { "__mul",          qlua_mul  },
    { "__div",          qlua_div  },
    { NULL,             NULL      }
};

/* complex matrices */
mMatComplex *
qlua_checkMatComplex(lua_State *L, int idx)
{
    return (mMatComplex *)qlua_checkMatrix(L, idx, mtnMatComplex, "Complex");
}

mMatComplex *
qlua_newMatComplex(lua_State *L, int sl, int sr)
{
    typedef double xc[2];

    return qlua_newMatrix(L, MSIZE(mMatComplex, sl, sr, xc), mtnMatComplex);
}

static int
mc_dims(lua_State *L)                                          /* (-1,+2,e) */
{
    mMatComplex *m = qlua_checkMatComplex(L, 1);

    lua_pushnumber(L, m->size_l);
    lua_pushnumber(L, m->size_r);

    return 2;
}

static int
mc_transpose(lua_State *L)
{
    mMatComplex *a = qlua_checkMatComplex(L, 1);
    mMatComplex *b = qlua_newMatComplex(L, a->size_r, a->size_l);
    int l, r;

    b->size_l = a->size_r;
    b->size_r = a->size_l;
    for (l = 0; l < b->size_l; l++) {
        for (r = 0; r < b->size_r; r++) {
            CM(b,l,r)[0] = CM(a,r,l)[0];
            CM(b,l,r)[1] = CM(a,r,l)[1];
        }
    }
    return 1;
}

static int
mc_adjoin(lua_State *L)
{
    mMatComplex *a = qlua_checkMatComplex(L, 1);
    mMatComplex *b = qlua_newMatComplex(L, a->size_r, a->size_l);
    int l, r;

    b->size_l = a->size_r;
    b->size_r = a->size_l;
    for (l = 0; l < b->size_l; l++) {
        for (r = 0; r < b->size_r; r++) {
            CM(b,l,r)[0] = CM(a,r,l)[0];
            CM(b,l,r)[1] = -CM(a,r,l)[1];
        }
    }
    return 1;
}

static int
mc_conj(lua_State *L)
{
    mMatComplex *a = qlua_checkMatComplex(L, 1);
    mMatComplex *b = qlua_newMatComplex(L, a->size_l, a->size_r);
    int l, r;

    b->size_l = a->size_r;
    b->size_r = a->size_l;
    for (l = 0; l < b->size_l; l++) {
        for (r = 0; r < b->size_r; r++) {
            CM(b,l,r)[0] = CM(a,l,r)[0];
            CM(b,l,r)[1] = -CM(a,l,r)[1];
        }
    }
    return 1;
}

static int
mc_trace(lua_State *L)                                         /* (-1,+2,e) */
{
    mMatComplex *m = qlua_checkMatComplex(L, 1);
    int i;
    double tr;
    double ti;
    QLA_Complex *t = qlua_newComplex(L);
    
    if (m->size_l != m->size_r)
        return luaL_error(L, "matrix:trace() expects square matrix");
    
    for (ti = tr = 0, i = 0; i < m->size_l; i++) {
        tr += RM(m,i,i)[0];
        ti += RM(m,i,i)[1];
    }
    QLA_real(*t) = tr;
    QLA_imag(*t) = ti;

    return 1;
}

static int
mc_get(lua_State *L)                                           /* (-2,+1,e) */
{
    switch (qlua_gettype(L, 2)) {
    case qTable: {
        mMatComplex *v = qlua_checkMatComplex(L, 1);
        int sl, sr;

        qlua_checkindex2(L, 2, "matrix get", &sl, &sr);

        if ((sl >= 0) && (sl < v->size_l) &&
            (sr >= 0) && (sr < v->size_r)) {
            QLA_Complex *z = qlua_newComplex(L);

            QLA_real(*z) = CM(v, sl, sr)[0];
            QLA_imag(*z) = CM(v, sl, sr)[1];
            return 1;
        }
        break;
    }
    case qString:
        return qlua_lookup(L, 2, opMatComplex);
    }
    return qlua_badindex(L, "matrix.complex[]");
}

static int
mc_put(lua_State *L)                                           /* (-3,+0,e) */
{
    mMatComplex *v = qlua_checkMatComplex(L, 1);
    int sl, sr;
    qlua_checkindex2(L, 2, "matrix set", &sl, &sr);

    if ((sl >= 0) && (sl < v->size_l) &&
        (sr >= 0) && (sr < v->size_r)) {
        switch (qlua_gettype(L, 3)) {
        case qReal: {
            double x = luaL_checknumber(L, 3);
            CM(v, sl, sr)[0] = x;
            CM(v, sl, sr)[1] = 0;
            return 0;
        }
        case qComplex: {
            QLA_Complex *z = qlua_checkComplex(L, 3);
            CM(v, sl, sr)[0] = QLA_real(*z);
            CM(v, sl, sr)[1] = QLA_imag(*z);
            return 0;
        }
        }
    }

    return qlua_badindex(L, "matrix.complex[]");
}

static int
mc_fmt(lua_State *L)                                           /* (-1,+1,e) */
{
    char fmt[72];
    mMatComplex *v = qlua_checkMatComplex(L, 1);

    sprintf(fmt, "matrix.complex[{%d,%d}]", v->size_l, v->size_r);
    lua_pushstring(L, fmt);

    return 1;
}

static int
mc_inverse(lua_State *L)
{
    mMatComplex *m = qlua_checkMatComplex(L, 1);
    mMatComplex *t;
    mMatComplex *r;

    if (m->size_l != m->size_r)
        return luaL_error(L, "square matrix expected");
    
    t = qlua_newMatComplex(L, m->size_l, m->size_l);
    r = qlua_newMatComplex(L, m->size_l, m->size_l);
    memcpy(t->val, m->val, 2 * m->size_l * m->size_l * sizeof (double));
    r->size_l = m->size_l;
    r->size_r = m->size_l;
    if (matrix_cinverse(m->size_l, t->val, r->val))
        return luaL_error(L, "inverting singular matrix");

    return 1;
}

static int
mc_eigen(lua_State *L)                                         /* (-1,+2,e) */
{
    mMatComplex *m = qlua_checkMatComplex(L, 1);
    mVecReal *lambda;
    mMatComplex *trans;

    if (m->size_l != m->size_r)
        return luaL_error(L, "matrix:eigen() expects square matrix");
    
    lambda = qlua_newVecReal(L, m->size_l);
    lambda->size = m->size_l;
    trans = qlua_newMatComplex(L, m->size_l, m->size_r);
    trans->size_l = m->size_l;
    trans->size_r = m->size_r;

    matrix_ceigenvec(L, m->size_l, m->val, lambda->val, trans->val);

    return 2;
}

static int
cm_add_cm(lua_State *L)
{
    mMatComplex *a = qlua_checkMatComplex(L, 1);
    mMatComplex *b = qlua_checkMatComplex(L, 2);
    int al = a->size_l;
    int ar = a->size_r;
    mMatComplex *r = qlua_newMatComplex(L, al, ar);
    int i;

    if ((al != b->size_l) || (ar != b->size_r))
        return luaL_error(L, "matrix sizes mismatch in m + m");

    r->size_l = al;
    r->size_r = ar;
    for (i = 0; i < 2 * al * ar; i++)
        r->val[i] = a->val[i] + b->val[i];

    return 1;
}

static int
cm_sub_cm(lua_State *L)
{
    mMatComplex *a = qlua_checkMatComplex(L, 1);
    mMatComplex *b = qlua_checkMatComplex(L, 2);
    int al = a->size_l;
    int ar = a->size_r;
    mMatComplex *r = qlua_newMatComplex(L, al, ar);
    int i;

    if ((al != b->size_l) || (ar != b->size_r))
        return luaL_error(L, "matrix sizes mismatch in m - m");

    r->size_l = al;
    r->size_r = ar;
    for (i = 0; i < 2 * al * ar; i++)
        r->val[i] = a->val[i] - b->val[i];

    return 1;
}

static int
do_crmul(lua_State *L, mMatComplex *a, double b)
{
    mMatComplex *r = qlua_newMatComplex(L, a->size_l, a->size_r);
    int n = 2 * a->size_l * a->size_r;
    int i;

    r->size_l = a->size_l;
    r->size_r = a->size_r;

    for (i = 0; i < n; i++)
        r->val[i] = a->val[i] * b;

    return 1;
}

static int
cm_mul_r(lua_State *L)
{
    mMatComplex *a = qlua_checkMatComplex(L, 1);
    double b = luaL_checknumber(L, 2);

    return do_crmul(L, a, b);
}

static int
r_mul_cm(lua_State *L)
{
    double b = luaL_checknumber(L, 1);
    mMatComplex *a = qlua_checkMatComplex(L, 2);

    return do_crmul(L, a, b);
}

static int
cm_div_r(lua_State *L)
{
    mMatComplex *a = qlua_checkMatComplex(L, 1);
    double b = luaL_checknumber(L, 2);

    return do_crmul(L, a, 1/b);
}

static int
do_ccmul(lua_State *L, mMatComplex *a, double br, double bi)
{
    mMatComplex *r = qlua_newMatComplex(L, a->size_l, a->size_r);
    int n = a->size_l * a->size_r;
    int i;

    r->size_l = a->size_l;
    r->size_r = a->size_r;

    for (i = 0; i < n; i++) {
        r->val[2*i] = a->val[2*i] * br - a->val[2*i+1] * bi;
        r->val[2*i+1] = a->val[2*i] * bi + a->val[2*i+1] * br;
    }

    return 1;
}

static int
cm_mul_c(lua_State *L)
{
    mMatComplex *a = qlua_checkMatComplex(L, 1);
    QLA_Complex *b = qlua_checkComplex(L, 2);

    return do_ccmul(L, a, QLA_real(*b), QLA_imag(*b));
}

static int
c_mul_cm(lua_State *L)
{
    QLA_Complex *b = qlua_checkComplex(L, 1);
    mMatComplex *a = qlua_checkMatComplex(L, 2);

    return do_ccmul(L, a, QLA_real(*b), QLA_imag(*b));
}

static int
cm_div_c(lua_State *L)
{
    mMatComplex *a = qlua_checkMatComplex(L, 1);
    QLA_Complex *b = qlua_checkComplex(L, 2);
    double br = QLA_real(*b);
    double bi = QLA_imag(*b);
    double h = hypot(br, bi);
    double n = 1/h;

    return do_ccmul(L, a, br * n * n, -bi * n * n);
}

static int
cm_mul_cm(lua_State *L)
{
    mMatComplex *a = qlua_checkMatComplex(L, 1);
    mMatComplex *b = qlua_checkMatComplex(L, 2);
    int al = a->size_l;
    int ar = a->size_r;
    int bl = b->size_l;
    int br = b->size_r;
    mMatComplex *r = qlua_newMatComplex(L, al, br);
    int i, j, k;

    if (ar != bl)
        return luaL_error(L, "matrix sizes mismatch in m * m");

    r->size_l = al;
    r->size_r = br;
    for (i = 0; i < al; i++) {
        for (j = 0; j < br; j++) {
            double vr = 0;
            double vi = 0;
            for (k = 0; k < ar; k++) {
                double ar = CM(a,i,k)[0];
                double ai = CM(a,i,k)[1];
                double br = CM(b,k,j)[0];
                double bi = CM(b,k,j)[1];
                vr += ar * br - ai * bi;
                vi += ar * bi + ai * br;
            }
            CM(r, i, j)[0] = vr;
            CM(r, i, j)[1] = vi;
        }
    }

    return 1;
}

static int
m_complex(lua_State *L)
{
    int sl, sr;
    mMatComplex *m;

    qlua_checkindex2(L, 1, "matrix size", &sl, &sr);

    m = qlua_newMatComplex(L, sl, sr);
    m->size_l = sl;
    m->size_r = sr;
    memset(m->val, 0, sl * sr * 2 * sizeof (double));

    return 1;
}

static const luaL_Reg MatComplexMethods[] = {
    { "dims",                 mc_dims       },
    { "transpose",            mc_transpose  },
    { "trace",                mc_trace      },
    { "conj",                 mc_conj       },
    { "adjoin",               mc_adjoin     },
    { "inverse",              mc_inverse    },
    { "hermitian_eigen",      mc_eigen      },
    { NULL, NULL}
};

static const luaL_Reg mtMatComplex[] = {
    { "__tostring",     mc_fmt    },
    { "__index",        mc_get    },
    { "__newindex",     mc_put    },
    { "__add",          qlua_add  },
    { "__sub",          qlua_sub  },
    { "__mul",          qlua_mul  },
    { "__div",          qlua_div  },
    { NULL,             NULL      }
};

/* matrix constructors */
static const luaL_Reg fMatrix[] = {
    { "real",          m_real    },
    { "complex",       m_complex },
    { NULL,            NULL      }
};

int
init_matrix(lua_State *L)
{
    luaL_register(L, matrix_ns,      fMatrix);
    qlua_metatable(L, mtnMatReal,    mtMatReal);
    qlua_metatable(L, opMatReal,     MatRealMethods);
    qlua_metatable(L, mtnMatComplex, mtMatComplex);
    qlua_metatable(L, opMatComplex,  MatComplexMethods);

    qlua_reg_add(qMatComplex, qMatComplex, cm_add_cm);
    qlua_reg_add(qMatReal,    qMatReal,    rm_add_rm);
    qlua_reg_div(qMatComplex, qComplex,    cm_div_c);
    qlua_reg_div(qMatComplex, qReal,       cm_div_r);
    qlua_reg_div(qMatReal,    qReal,       rm_div_r);
    qlua_reg_mul(qComplex,    qMatComplex, c_mul_cm);
    qlua_reg_mul(qMatComplex, qComplex,    cm_mul_c);
    qlua_reg_mul(qMatComplex, qMatComplex, cm_mul_cm);
    qlua_reg_mul(qMatComplex, qReal,       cm_mul_r);
    qlua_reg_mul(qMatReal,    qMatReal,    rm_mul_rm);
    qlua_reg_mul(qMatReal,    qReal,       rm_mul_r);
    qlua_reg_mul(qReal,       qMatComplex, r_mul_cm);
    qlua_reg_mul(qReal,       qMatReal,    r_mul_rm);
    qlua_reg_sub(qMatComplex, qMatComplex, cm_sub_cm);
    qlua_reg_sub(qMatReal,    qMatReal,    rm_sub_rm);

    return 0;
}

int
fini_matrix(lua_State *L)
{
    return 0;
}
