#include <qlua.h>                                                    /* DEPS */
#include <qmatrix.h>                                                 /* DEPS */
#include <qvector.h>                                                 /* DEPS */
#include <matrix.h>                                                  /* DEPS */
#include <string.h>

const char mtnMatReal[]    = "matrix.mtReal";
static const char opMatReal[] = "matrix.mtReal.ops";

static char matrix_ns[] = "matrix";

#define MSIZE(s,a,b,t) (sizeof (s) + ((a)*(b)-1)*sizeof (t))

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

mMatReal *
qlua_checkMatReal(lua_State *L, int idx)
{
    return (mMatReal *)qlua_checkMatrix(L, idx, mtnMatReal, "Real");
}

static void *
qlua_newMatrix(lua_State *L, int s, const char *mt)
{
    void *v = lua_newuserdata(L, s);

    luaL_getmetatable(L, mt);
    lua_setmetatable(L, -2);
    
    return v;
}

mMatReal *
qlua_newMatReal(lua_State *L, int sl, int sr)
{
    return qlua_newMatrix(L, MSIZE(mMatReal, sl, sr, double), mtnMatReal);
}

/* real matrices */
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
            b->val[l + r * b->size_l] = a->val[r + l * a->size_l];
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
        tr += m->val[i * (1 + m->size_l)];
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
            lua_pushnumber(L, v->val[sl + sr * v->size_l]);
            return 1;
        }
        break;
    }
    case qString:
        return qlua_lookup(L, 2, opMatReal);
    }
    return qlua_badindex(L, "vector.real[]");
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
        v->val[sl + sr * v->size_l] = x;
        return 0;
    }

    return qlua_badindex(L, "matrix.real[]");
}

static int
m_add_m(lua_State *L)
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
m_sub_m(lua_State *L)
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
m_mul_m(lua_State *L)
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
                v += a->val[i + k * al] * b->val[k + j * bl];
            r->val[i + j * al] = v;
        }
    }

    return 1;
}

static int
do_mul(lua_State *L, mMatReal *a, double b)
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
m_mul_r(lua_State *L)
{
    mMatReal *a = qlua_checkMatReal(L, 1);
    double b = luaL_checknumber(L, 2);

    return do_mul(L, a, b);
}

static int
r_mul_m(lua_State *L)
{
    double b = luaL_checknumber(L, 1);
    mMatReal *a = qlua_checkMatReal(L, 2);

    return do_mul(L, a, b);
}

static int
m_div_r(lua_State *L)
{
    mMatReal *a = qlua_checkMatReal(L, 2);
    double b = luaL_checknumber(L, 1);

    return do_mul(L, a, 1/b);
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
    { "dims",       md_dims       },
    { "eigen",      md_eigen      },
    { "transpose",  md_transpose  },
    { "trace",      md_trace      },
    { "inverse",    md_inverse    },
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

/* vector constructors */
static const luaL_Reg fMatrix[] = {
    { "real",      m_real    },
    { NULL,            NULL      }
};

int
init_matrix(lua_State *L)
{
    luaL_register(L, matrix_ns,      fMatrix);
    qlua_metatable(L, mtnMatReal,    mtMatReal);
    qlua_metatable(L, opMatReal,     MatRealMethods);
    qlua_reg_add(qMatReal, qMatReal, m_add_m);
    qlua_reg_sub(qMatReal, qMatReal, m_sub_m);
    qlua_reg_mul(qMatReal, qMatReal, m_mul_m);
    qlua_reg_mul(qReal,    qMatReal, r_mul_m);
    qlua_reg_mul(qMatReal, qReal,    m_mul_r);
    qlua_reg_div(qMatReal, qReal,    m_div_r);
    return 0;
}

int
fini_matrix(lua_State *L)
{
    return 0;
}

