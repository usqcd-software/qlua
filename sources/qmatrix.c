/* Rings over matrices, real and complex, not mixed.
 * Modules are not implemented because there are too many of them
 */
#include "qlua.h"                                                    /* DEPS */
#include "qcomplex.h"                                                /* DEPS */
#include "qmatrix.h"                                                 /* DEPS */
#include "qvector.h"                                                 /* DEPS */
#include <string.h>
#include <math.h>

const char mtnMatReal[]           = "matrix.mtReal";
const char mtnMatComplex[]        = "matrix.mtComplex";
static const char opMatReal[]     = "matrix.mtReal.ops";
static const char opMatComplex[]  = "matrix.mtComplex.ops";

static char matrix_ns[] = "matrix";

static gsl_permutation *
alloc_GSLPermutation(lua_State *L, int size)
{
    gsl_permutation *p = gsl_permutation_alloc(size);
    if (p == 0) {
        lua_gc(L, LUA_GCCOLLECT, 0);
        p = gsl_permutation_alloc(size);
        if (p == 0)
            luaL_error(L, "not enough memory");
    }
	return p;
}

static gsl_vector *
alloc_GSLRealVector(lua_State *L, int size)
{
	gsl_vector *w = gsl_vector_alloc(size);
	if (w == 0) {
		lua_gc(L, LUA_GCCOLLECT, 0);
		w = gsl_vector_alloc(size);
		if (w == 0)
			luaL_error(L, "not enough memory");
	}
	return w;
}

static gsl_vector *
qlua2gsl_real_vector(lua_State *L, int idx, int e_size)
{
	mVecReal *v = qlua_checkVecReal(L, idx);
	gsl_vector *w;
	int i;

	if ((e_size >= 0) && (v->size != e_size))
		luaL_error(L, "unexpected vector length");

	w = alloc_GSLRealVector(L, v->size);
	for (i = 0; i < v->size; i++) {
		gsl_vector_set(w, i, v->val[i]);
	}
	return w;
}

static mVecReal *
gsl2qlua_real_vector(lua_State *L, int size, gsl_vector *w)
{
	mVecReal *v = qlua_newVecReal(L, size);
	int i;

	for (i = 0; i < size; i++) {
		v->val[i] = gsl_vector_get(w, i);
	}
	gsl_vector_free(w);
	return v;
}

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
    mMatReal *m = (mMatReal *)qlua_newMatrix(L, sizeof (mMatReal), mtnMatReal);

    m->l_size = sl;
    m->r_size = sr;
    m->m = gsl_matrix_calloc(sl, sr);
    if (m->m == 0) {
        lua_gc(L, LUA_GCCOLLECT, 0);
        m->m = gsl_matrix_calloc(sl, sr);
        if (m->m == 0)
            luaL_error(L, "not enough memory");
    }
    return m;
}

static int
md_fmt(lua_State *L)                                           /* (-1,+1,e) */
{
    char fmt[72];
    mMatReal *v = qlua_checkMatReal(L, 1);

    sprintf(fmt, "matrix.real[{%d,%d}]", v->l_size, v->r_size);
    lua_pushstring(L, fmt);

    return 1;
}

static int
md_gc(lua_State *L)
{
    mMatReal *v = qlua_checkMatReal(L, 1);
    
    if (v->m)
        gsl_matrix_free(v->m);
    v->m = 0;

    return 0;
}

static int
md_inverse(lua_State *L)
{
    mMatReal *m = qlua_checkMatReal(L, 1);
    mMatReal *lu;
    mMatReal *r;
    gsl_permutation *p;
    int signum;

    if (m->l_size != m->r_size)
        return luaL_error(L, "square matrix expected");
    
    lu = qlua_newMatReal(L, m->l_size, m->l_size);
    r = qlua_newMatReal(L, m->l_size, m->l_size);
    gsl_matrix_memcpy(lu->m, m->m);
    p = alloc_GSLPermutation(L, m->l_size);
    gsl_linalg_LU_decomp(lu->m, p, &signum);
    if (gsl_linalg_LU_invert(lu->m, p, r->m))
        luaL_error(L, "matrix:inverse() failed");
    
    gsl_permutation_free(p);
    return 1;
}

static int
md_det(lua_State *L)                                            /* (-1,+1,e) */
{
    mMatReal *m = qlua_checkMatReal(L, 1);
    mMatReal *lu;
    gsl_permutation *p;
    int signum;
    double d;

    if (m->l_size != m->r_size)
        return luaL_error(L, "square matrix expected");
    
    lu = qlua_newMatReal(L, m->l_size, m->l_size);
    gsl_matrix_memcpy(lu->m, m->m);
    p = gsl_permutation_alloc(m->l_size);
    if (p == 0) {
        lua_gc(L, LUA_GCCOLLECT, 0);
        p = gsl_permutation_alloc(m->l_size);
        if (p == 0)
            luaL_error(L, "not enough memory");
    }
    gsl_linalg_LU_decomp(lu->m, p, &signum);
    d = gsl_linalg_LU_det(lu->m, signum);
    gsl_permutation_free(p);

    lua_pushnumber(L, d);
    return 1;
}

static int
md_qr(lua_State *L)                                            /* (-1,+2,e) */
{
    mMatReal *m = qlua_checkMatReal(L, 1);
    mMatReal *qr = qlua_newMatReal(L, m->l_size, m->r_size);
    mMatReal *q = qlua_newMatReal(L, m->l_size, m->l_size);
    mMatReal *r = qlua_newMatReal(L, m->l_size, m->r_size);
    int nm = m->l_size < m->r_size? m->l_size: m->r_size;
    gsl_vector *tau;

    gsl_matrix_memcpy(qr->m, m->m);
    tau = gsl_vector_alloc(nm);
    if (tau == 0) {
        lua_gc(L, LUA_GCCOLLECT, 0);
        tau = gsl_vector_alloc(nm);
        if (tau == 0)
            luaL_error(L, "not enough memory");
    }
    if (gsl_linalg_QR_decomp(qr->m, tau))
        luaL_error(L, "matrix:qr() failed");
    
    if (gsl_linalg_QR_unpack(qr->m, tau, q->m, r->m))
        luaL_error(L, "matrix:qr() failed");
    gsl_vector_free(tau);
    
    return 2;
}

static int
md_solve(lua_State *L)
{
	mMatReal *m = qlua_checkMatReal(L, 1);
	mMatReal *lu;
	gsl_permutation *p;
	gsl_vector *ss, *x;
	
	int signum;

	if (m->l_size != m->r_size)
		return luaL_error(L, "matrix:solve() expects a square matrix");
	ss = qlua2gsl_real_vector(L, 2, m->l_size);
	lu = qlua_newMatReal(L, m->l_size, m->r_size);
	gsl_matrix_memcpy(lu->m, m->m);
	p = alloc_GSLPermutation(L, m->l_size);
	x = alloc_GSLRealVector(L, m->l_size);
	gsl_linalg_LU_decomp(lu->m, p, &signum);
	if (gsl_linalg_LU_solve(lu->m, p, ss, x)) {
		luaL_error(L, "matrix:solve() failed");
	}
	gsl_permutation_free(p);
	gsl2qlua_real_vector(L, m->l_size, x);
	return 1;
}

static int
slice_out(lua_State *L)
{
    return luaL_error(L, "matrix:eigen() matrix slice out of bounds");
}

static int
md_eigen(lua_State *L)                                         /* (-1,+2,e) */
{
    mMatReal *m = qlua_checkMatReal(L, 1);
    gsl_matrix_view mx;
    gsl_eigen_symmv_workspace *w;
    gsl_vector *ev;
    mVecReal *lambda;
    mMatReal *trans;
    mMatReal *tmp;
    int n;
    int i;
    int lo, hi;

    switch (lua_gettop(L)) {
    case 1:
        if (m->l_size != m->r_size)
            return luaL_error(L, "matrix:eigen() expects square matrix");
        lo = 0;
        hi = m->l_size;
        break;
    case 2:
        lo = 0;
        hi = luaL_checkint(L, 2);
        if ((hi > m->l_size) || (hi > m->r_size))
            return slice_out(L);
        break;
    case 3:
        lo = luaL_checkint(L, 2);
        hi = luaL_checkint(L, 3);
        if ((lo >= hi) ||
            (lo > m->l_size) || (lo > m->r_size) ||
            (hi > m->l_size) || (hi > m->r_size))
            return slice_out(L);
        break;
    default:
        return luaL_error(L, "matrix:eigen(): illegal arguments");
    }

    n = hi - lo;
    mx = gsl_matrix_submatrix(m->m, lo, lo, n, n);
    tmp = qlua_newMatReal(L, n, n);
    gsl_matrix_memcpy(tmp->m, &mx.matrix);
    lambda = qlua_newVecReal(L, n);
    trans = qlua_newMatReal(L, n, n);

    ev = gsl_vector_alloc(n);
    if (ev == 0) {
        lua_gc(L, LUA_GCCOLLECT, 0);
        ev = gsl_vector_alloc(n);
        if (ev == 0)
            luaL_error(L, "not enough memory");
    }

    w = gsl_eigen_symmv_alloc(n);
    if (w == 0) {
        lua_gc(L, LUA_GCCOLLECT, 0);
        w = gsl_eigen_symmv_alloc(n);
        if (w == 0)
            luaL_error(L, "not enough memory");
    }
    
    if (gsl_eigen_symmv(tmp->m, ev, trans->m, w))
        luaL_error(L, "matrix:eigen() failed");

    if (gsl_eigen_symmv_sort(ev, trans->m, GSL_EIGEN_SORT_VAL_ASC))
        luaL_error(L, "matrix:eigen() eigenvalue ordering failed");

    for (i = 0; i < n; i++)
        lambda->val[i] = gsl_vector_get(ev, i);

    gsl_vector_free(ev);
    gsl_eigen_symmv_free(w);

    return 2;
}

static int
md_transpose(lua_State *L)
{
    mMatReal *a = qlua_checkMatReal(L, 1);
    mMatReal *b = qlua_newMatReal(L, a->r_size, a->r_size);
    int l, r;

    for (l = 0; l < b->l_size; l++) {
        for (r = 0; r < b->r_size; r++) {
            gsl_matrix_set(b->m, l, r, gsl_matrix_get(a->m, r, l));
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
    
    if (m->l_size != m->r_size)
        return luaL_error(L, "matrix:trace() expects square matrix");
    
    for (tr = 0, i = 0; i < m->l_size; i++)
        tr += gsl_matrix_get(m->m,i,i);
    lua_pushnumber(L, tr);

    return 1;
}

static int
md_dims(lua_State *L)                                          /* (-1,+2,e) */
{
    mMatReal *m = qlua_checkMatReal(L, 1);

    lua_pushnumber(L, m->l_size);
    lua_pushnumber(L, m->r_size);

    return 2;
}

static int
md_get(lua_State *L)                                           /* (-2,+1,e) */
{
    switch (qlua_qtype(L, 2)) {
    case qTable: {
        mMatReal *v = qlua_checkMatReal(L, 1);
        int sl, sr;

        qlua_checkindex2(L, 2, "matrix get", &sl, &sr);

        if ((sl >= 0) && (sl < v->l_size) &&
            (sr >= 0) && (sr < v->r_size)) {
            lua_pushnumber(L, gsl_matrix_get(v->m, sl, sr));
            return 1;
        }
        break;
    }
    case qString:
        return qlua_lookup(L, 2, opMatReal);
    default:
        break;
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

    if ((sl >= 0) && (sl < v->l_size) &&
        (sr >= 0) && (sr < v->r_size)) {
        gsl_matrix_set(v->m, sl, sr, x);
        return 0;
    }

    return qlua_badindex(L, "matrix.real[]");
}

static int
rm_add_rm(lua_State *L)
{
    mMatReal *a = qlua_checkMatReal(L, 1);
    mMatReal *b = qlua_checkMatReal(L, 2);
    int al = a->l_size;
    int ar = a->r_size;
    mMatReal *r = qlua_newMatReal(L, al, ar);

    if ((al != b->l_size) || (ar != b->r_size))
        return luaL_error(L, "matrix sizes mismatch in m + m");

    gsl_matrix_memcpy(r->m, a->m);
    gsl_matrix_add(r->m, b->m);

    return 1;
}

static int
rm_sub_rm(lua_State *L)
{
    mMatReal *a = qlua_checkMatReal(L, 1);
    mMatReal *b = qlua_checkMatReal(L, 2);
    int al = a->l_size;
    int ar = a->r_size;
    mMatReal *r = qlua_newMatReal(L, al, ar);

    if ((al != b->l_size) || (ar != b->r_size))
        return luaL_error(L, "matrix sizes mismatch in m - m");

    gsl_matrix_memcpy(r->m, a->m);
    gsl_matrix_sub(r->m, b->m);

    return 1;
}

static int
rm_mul_rm(lua_State *L)
{
    mMatReal *a = qlua_checkMatReal(L, 1);
    mMatReal *b = qlua_checkMatReal(L, 2);
    int al = a->l_size;
    int ar = a->r_size;
    int bl = b->l_size;
    int br = b->r_size;
    mMatReal *r = qlua_newMatReal(L, al, br);

    if (ar != bl)
        return luaL_error(L, "matrix sizes mismatch in m * m");

    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, a->m, b->m, 0.0, r->m);

    return 1;
}

static int
do_rrmul(lua_State *L, mMatReal *a, double b)
{
    mMatReal *r = qlua_newMatReal(L, a->l_size, a->r_size);

    gsl_matrix_memcpy(r->m, a->m);
    gsl_matrix_scale(r->m, b);

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

    qlua_checkindex2(L, 1, "matrix size", &sl, &sr);

    qlua_newMatReal(L, sl, sr);

    return 1;
}

static const luaL_Reg MatRealMethods[] = {
    { "dims",                 md_dims       },
    { "transpose",            md_transpose  },
    { "trace",                md_trace      },
    { "inverse",              md_inverse    },
    { "symmetric_eigen",      md_eigen      },
    { "qr",                   md_qr         },
    { "det",                  md_det        },
	{ "solve",                md_solve      },
    { NULL, NULL}
};

static const luaL_Reg mtMatReal[] = {
    { "__tostring",     md_fmt    },
    { "__gc",           md_gc     },
    { "__index",        md_get    },
    { "__newindex",     md_put    },
    { "__add",          qlua_add  },
    { "__sub",          qlua_sub  },
    { "__mul",          qlua_mul  },
    { "__div",          qlua_div  },
    { NULL,             NULL      }
};

static gsl_vector_complex *
alloc_GSLComplexVector(lua_State *L, int size)
{
	gsl_vector_complex *w = gsl_vector_complex_alloc(size);
	if (w == 0) {
		lua_gc(L, LUA_GCCOLLECT, 0);
		w = gsl_vector_complex_alloc(size);
		if (w == 0)
			luaL_error(L, "not enough memory");
	}
	return w;
}

static gsl_vector_complex *
qlua2gsl_complex_vector(lua_State *L, int idx, int e_size)
{
	mVecComplex *v = qlua_checkVecComplex(L, idx);
	gsl_vector_complex *w;
	int i;

	if ((e_size >= 0) && (v->size != e_size))
		luaL_error(L, "unexpected vector length");

	w = alloc_GSLComplexVector(L, v->size);
	for (i = 0; i < v->size; i++) {
		gsl_complex zz;
		QLA_D_Complex *vv = &v->val[i];

		GSL_SET_COMPLEX(&zz, QLA_real(*vv), QLA_imag(*vv));
		gsl_vector_complex_set(w, i, zz);
	}
	return w;
}

static mVecComplex *
gsl2qlua_complex_vector(lua_State *L, int size, gsl_vector_complex *w)
{
	mVecComplex *v = qlua_newVecComplex(L, size);
	int i;

	for (i = 0; i < size; i++) {
		gsl_complex zz = gsl_vector_complex_get(w, i);
		QLA_real(v->val[i]) = GSL_REAL(zz);
		QLA_imag(v->val[i]) = GSL_IMAG(zz);
	}
	gsl_vector_complex_free(w);
	return v;
}

/* complex matrices */
mMatComplex *
qlua_checkMatComplex(lua_State *L, int idx)
{
    return (mMatComplex *)qlua_checkMatrix(L, idx, mtnMatComplex, "Complex");
}

mMatComplex *
qlua_newMatComplex(lua_State *L, int sl, int sr)
{
    mMatComplex *m = qlua_newMatrix(L, sizeof (mMatComplex), mtnMatComplex);
    m->l_size = sl;
    m->r_size = sr;
    m->m = gsl_matrix_complex_calloc(sl, sr);
    if (m->m == 0) {
        lua_gc(L, LUA_GCCOLLECT, 0);
        m->m = gsl_matrix_complex_calloc(sl, sr);
        if (m->m == 0)
            luaL_error(L, "not enough memory");
    }
    return m;
}

static int
mc_dims(lua_State *L)                                          /* (-1,+2,e) */
{
    mMatComplex *m = qlua_checkMatComplex(L, 1);

    lua_pushnumber(L, m->l_size);
    lua_pushnumber(L, m->r_size);

    return 2;
}

static int
mc_transpose(lua_State *L)
{
    mMatComplex *a = qlua_checkMatComplex(L, 1);
    mMatComplex *b = qlua_newMatComplex(L, a->r_size, a->l_size);
    int l, r;

    for (l = 0; l < b->l_size; l++) {
        for (r = 0; r < b->r_size; r++) {
            gsl_matrix_complex_set(b->m, l, r,
                                   gsl_matrix_complex_get(a->m, r, l));
        }
    }
    return 1;
}

static int
mc_adjoin(lua_State *L)
{
    mMatComplex *a = qlua_checkMatComplex(L, 1);
    mMatComplex *b = qlua_newMatComplex(L, a->r_size, a->l_size);
    int l, r;

    for (l = 0; l < b->l_size; l++) {
        for (r = 0; r < b->r_size; r++) {
            gsl_complex q;
            gsl_complex z = gsl_matrix_complex_get(a->m, r, l);
            GSL_SET_COMPLEX(&q, GSL_REAL(z), -GSL_IMAG(z));
            gsl_matrix_complex_set(b->m, l, r, q);
        }
    }
    return 1;
}

static int
mc_conj(lua_State *L)
{
    mMatComplex *a = qlua_checkMatComplex(L, 1);
    mMatComplex *b = qlua_newMatComplex(L, a->l_size, a->r_size);
    int l, r;

    for (l = 0; l < b->l_size; l++) {
        for (r = 0; r < b->r_size; r++) {
            gsl_complex q;
            gsl_complex z = gsl_matrix_complex_get(a->m, l, r);
            GSL_SET_COMPLEX(&q, GSL_REAL(z), -GSL_IMAG(z));
            gsl_matrix_complex_set(b->m, l, r, q);
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
    QLA_D_Complex *t = qlua_newComplex(L);
    
    if (m->l_size != m->r_size)
        return luaL_error(L, "matrix:trace() expects square matrix");
    
    for (ti = tr = 0, i = 0; i < m->l_size; i++) {
        gsl_complex z = gsl_matrix_complex_get(m->m, i, i);
        tr += GSL_REAL(z);
        ti += GSL_IMAG(z);
    }
    QLA_real(*t) = tr;
    QLA_imag(*t) = ti;

    return 1;
}

static int
mc_get(lua_State *L)                                           /* (-2,+1,e) */
{
    switch (qlua_qtype(L, 2)) {
    case qTable: {
        mMatComplex *v = qlua_checkMatComplex(L, 1);
        int sl, sr;

        qlua_checkindex2(L, 2, "matrix get", &sl, &sr);

        if ((sl >= 0) && (sl < v->l_size) &&
            (sr >= 0) && (sr < v->r_size)) {
            QLA_D_Complex *z = qlua_newComplex(L);
            gsl_complex zz = gsl_matrix_complex_get(v->m, sl, sr);

            QLA_real(*z) = GSL_REAL(zz);
            QLA_imag(*z) = GSL_IMAG(zz);
            return 1;
        }
        break;
    }
    case qString:
        return qlua_lookup(L, 2, opMatComplex);
    default:
        break;
    }
    return qlua_badindex(L, "matrix.complex[]");
}

static int
mc_put(lua_State *L)                                           /* (-3,+0,e) */
{
    mMatComplex *v = qlua_checkMatComplex(L, 1);
    int sl, sr;
    

    qlua_checkindex2(L, 2, "matrix set", &sl, &sr);

    if ((sl >= 0) && (sl < v->l_size) &&
        (sr >= 0) && (sr < v->r_size)) {
        switch (qlua_qtype(L, 3)) {
        case qReal: {
            gsl_complex z;
            double x = luaL_checknumber(L, 3);
            GSL_SET_COMPLEX(&z, x, 0);
            gsl_matrix_complex_set(v->m, sl, sr, z);
            return 0;
        }
        case qComplex: {
            gsl_complex zz;
            QLA_D_Complex *z = qlua_checkComplex(L, 3);
            GSL_SET_COMPLEX(&zz, QLA_real(*z), QLA_imag(*z));
            gsl_matrix_complex_set(v->m, sl, sr, zz);
            return 0;
        }
        default:
            break;
        }
    }
    return qlua_badindex(L, "matrix.complex[]");
}

static int
mc_fmt(lua_State *L)                                           /* (-1,+1,e) */
{
    char fmt[72];
    mMatComplex *v = qlua_checkMatComplex(L, 1);

    sprintf(fmt, "matrix.complex[{%d,%d}]", v->l_size, v->r_size);
    lua_pushstring(L, fmt);

    return 1;
}

static int
mc_gc(lua_State *L)
{
    mMatComplex *v = qlua_checkMatComplex(L, 1);

    if (v->m)
        gsl_matrix_complex_free(v->m);
    v->m = 0;

    return 0;
}

static int
mc_inverse(lua_State *L)
{
    mMatComplex *m = qlua_checkMatComplex(L, 1);
    mMatComplex *lu;
    mMatComplex *r;
    gsl_permutation *p;
    int signum;

    if (m->l_size != m->r_size)
        return luaL_error(L, "square matrix expected");
    
    lu = qlua_newMatComplex(L, m->l_size, m->l_size);
    r = qlua_newMatComplex(L, m->l_size, m->l_size);
    gsl_matrix_complex_memcpy(lu->m, m->m);
    p = gsl_permutation_alloc(m->l_size);
    if (p == 0) {
        lua_gc(L, LUA_GCCOLLECT, 0);
        p = gsl_permutation_alloc(m->l_size);
        if (p == 0)
            luaL_error(L, "not enough memory");
    }
    gsl_linalg_complex_LU_decomp(lu->m, p, &signum);
    if (gsl_linalg_complex_LU_invert(lu->m, p, r->m))
        luaL_error(L, "matrix:inverse() failed");
    
    gsl_permutation_free(p);
    return 1;
}

static int
mc_det(lua_State *L)
{
    mMatComplex *m = qlua_checkMatComplex(L, 1);
    mMatComplex *lu;
    mMatComplex *r;
    gsl_permutation *p;
    int signum;
    gsl_complex d;
    QLA_D_Complex *z;

    if (m->l_size != m->r_size)
        return luaL_error(L, "square matrix expected");
    
    lu = qlua_newMatComplex(L, m->l_size, m->l_size);
    r = qlua_newMatComplex(L, m->l_size, m->l_size);
    gsl_matrix_complex_memcpy(lu->m, m->m);
    p = gsl_permutation_alloc(m->l_size);
    if (p == 0) {
        lua_gc(L, LUA_GCCOLLECT, 0);
        p = gsl_permutation_alloc(m->l_size);
        if (p == 0)
            luaL_error(L, "not enough memory");
    }
    gsl_linalg_complex_LU_decomp(lu->m, p, &signum);
    d = gsl_linalg_complex_LU_det(lu->m, signum);
    gsl_permutation_free(p);
    z = qlua_newComplex(L);
    QLA_real(*z) = GSL_REAL(d);
    QLA_imag(*z) = GSL_IMAG(d);

    return 1;
}

static int
mc_qr(lua_State *L)                                            /* (-1,+2,e) */
{
    mMatComplex *m = qlua_checkMatComplex(L, 1);
    mMatComplex *qr = qlua_newMatComplex(L, m->l_size, m->r_size);
    mMatComplex *q = qlua_newMatComplex(L, m->l_size, m->l_size);
    mMatComplex *r = qlua_newMatComplex(L, m->l_size, m->r_size);
    int nm = m->l_size < m->r_size? m->l_size: m->r_size;
    gsl_vector_complex *tau;

    gsl_matrix_complex_memcpy(qr->m, m->m);
    tau = gsl_vector_complex_alloc(nm);
    if (tau == 0) {
        lua_gc(L, LUA_GCCOLLECT, 0);
        tau = gsl_vector_complex_alloc(nm);
        if (tau == 0)
            luaL_error(L, "not enough memory");
    }
    if (gsl_linalg_complex_QR_decomp(qr->m, tau))
        luaL_error(L, "matrix:qr() failed");
    
    if (gsl_linalg_complex_QR_unpack(qr->m, tau, q->m, r->m))
        luaL_error(L, "matrix:qr() failed");
    gsl_vector_complex_free(tau);
    
    return 2;
}
static int
mc_solve(lua_State *L)
{
	mMatComplex *m = qlua_checkMatComplex(L, 1);
	mMatComplex *lu;
	gsl_permutation *p;
	gsl_vector_complex *ss, *x;
	
	int signum;

	if (m->l_size != m->r_size)
		return luaL_error(L, "matrix:solve() expects a square matrix");
	ss = qlua2gsl_complex_vector(L, 2, m->l_size);
	lu = qlua_newMatComplex(L, m->l_size, m->r_size);
	gsl_matrix_complex_memcpy(lu->m, m->m);
	p = alloc_GSLPermutation(L, m->l_size);
	x = alloc_GSLComplexVector(L, m->l_size);
	gsl_linalg_complex_LU_decomp(lu->m, p, &signum);
	if (gsl_linalg_complex_LU_solve(lu->m, p, ss, x)) {
		luaL_error(L, "matrix:solve() failed");
	}
	gsl_permutation_free(p);
	gsl2qlua_complex_vector(L, m->l_size, x);
	return 1;
}

static int
mc_eigen(lua_State *L)                                         /* (-1,+2,e) */
{
    mMatComplex *m = qlua_checkMatComplex(L, 1);
    gsl_matrix_complex_view mx;
    gsl_eigen_hermv_workspace *w;
    gsl_vector *ev;
    mVecReal *lambda;
    mMatComplex *trans;
    mMatComplex *tmp;
    int n;
    int i;
    int lo, hi;

    switch (lua_gettop(L)) {
    case 1:
        if (m->l_size != m->r_size)
            return luaL_error(L, "matrix:eigen() expects square matrix");
        lo = 0;
        hi = m->l_size;
        break;
    case 2:
        lo = 0;
        hi = luaL_checkint(L, 2);
        if ((hi > m->l_size) || (hi > m->r_size))
            return slice_out(L);
        break;
    case 3:
        lo = luaL_checkint(L, 2);
        hi = luaL_checkint(L, 3);
        if ((lo >= hi) ||
            (lo > m->l_size) || (lo > m->r_size) ||
            (hi > m->l_size) || (hi > m->r_size))
            return slice_out(L);
        break;
    default:
        return luaL_error(L, "matrix:eigen(): illegal arguments");
    }

    n = hi - lo;
    mx = gsl_matrix_complex_submatrix(m->m, lo, lo, n, n);
    tmp = qlua_newMatComplex(L, n, n);
    gsl_matrix_complex_memcpy(tmp->m, &mx.matrix);
    lambda = qlua_newVecReal(L, n);
    trans = qlua_newMatComplex(L, n, n);

    ev = gsl_vector_alloc(n);
    if (ev == 0) {
        lua_gc(L, LUA_GCCOLLECT, 0);
        ev = gsl_vector_alloc(n);
        if (ev == 0)
            luaL_error(L, "not enough memory");
    }

    w = gsl_eigen_hermv_alloc(n);
    if (w == 0) {
        lua_gc(L, LUA_GCCOLLECT, 0);
        w = gsl_eigen_hermv_alloc(n);
        if (w == 0)
            luaL_error(L, "not enough memory");
    }
    
    if (gsl_eigen_hermv(tmp->m, ev, trans->m, w))
        luaL_error(L, "matrix:eigen() failed");

    if (gsl_eigen_hermv_sort(ev, trans->m, GSL_EIGEN_SORT_VAL_ASC))
        luaL_error(L, "matrix:eigen() eigenvalue ordering failed");

    for (i = 0; i < n; i++)
        lambda->val[i] = gsl_vector_get(ev, i);

    gsl_vector_free(ev);
    gsl_eigen_hermv_free(w);

    return 2;
}

static int
cm_add_cm(lua_State *L)
{
    mMatComplex *a = qlua_checkMatComplex(L, 1);
    mMatComplex *b = qlua_checkMatComplex(L, 2);
    int al = a->l_size;
    int ar = a->r_size;
    mMatComplex *r = qlua_newMatComplex(L, al, ar);

    if ((al != b->l_size) || (ar != b->r_size))
        return luaL_error(L, "matrix sizes mismatch in m + m");

    gsl_matrix_complex_memcpy(r->m, a->m);
    gsl_matrix_complex_add(r->m, b->m);

    return 1;
}

static int
cm_sub_cm(lua_State *L)
{
    mMatComplex *a = qlua_checkMatComplex(L, 1);
    mMatComplex *b = qlua_checkMatComplex(L, 2);
    int al = a->l_size;
    int ar = a->r_size;
    mMatComplex *r = qlua_newMatComplex(L, al, ar);

    if ((al != b->l_size) || (ar != b->r_size))
        return luaL_error(L, "matrix sizes mismatch in m + m");

    gsl_matrix_complex_memcpy(r->m, a->m);
    gsl_matrix_complex_sub(r->m, b->m);

    return 1;
}

static int
do_ccmul(lua_State *L, mMatComplex *a, double b_re, double b_im)
{
    mMatComplex *r = qlua_newMatComplex(L, a->l_size, a->r_size);
    gsl_complex z;

    gsl_matrix_complex_memcpy(r->m, a->m);
    GSL_SET_COMPLEX(&z, b_re, b_im);
    gsl_matrix_complex_scale(r->m, z);

    return 1;
}

static int
cm_mul_r(lua_State *L)
{
    mMatComplex *a = qlua_checkMatComplex(L, 1);
    double b = luaL_checknumber(L, 2);

    return do_ccmul(L, a, b, 0);
}

static int
r_mul_cm(lua_State *L)
{
    double b = luaL_checknumber(L, 1);
    mMatComplex *a = qlua_checkMatComplex(L, 2);

    return do_ccmul(L, a, b, 0);
}

static int
cm_div_r(lua_State *L)
{
    mMatComplex *a = qlua_checkMatComplex(L, 1);
    double b = luaL_checknumber(L, 2);

    return do_ccmul(L, a, 1/b, 0);
}

static int
cm_mul_c(lua_State *L)
{
    mMatComplex *a = qlua_checkMatComplex(L, 1);
    QLA_D_Complex *b = qlua_checkComplex(L, 2);

    return do_ccmul(L, a, QLA_real(*b), QLA_imag(*b));
}

static int
c_mul_cm(lua_State *L)
{
    QLA_D_Complex *b = qlua_checkComplex(L, 1);
    mMatComplex *a = qlua_checkMatComplex(L, 2);

    return do_ccmul(L, a, QLA_real(*b), QLA_imag(*b));
}

static int
cm_div_c(lua_State *L)
{
    mMatComplex *a = qlua_checkMatComplex(L, 1);
    QLA_D_Complex *b = qlua_checkComplex(L, 2);
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
    int al = a->l_size;
    int ar = a->r_size;
    int bl = b->l_size;
    int br = b->r_size;
    mMatComplex *r = qlua_newMatComplex(L, al, br);
    gsl_complex z1, z0;

    if (ar != bl)
        return luaL_error(L, "matrix sizes mismatch in m * m");

    GSL_SET_COMPLEX(&z1, 1.0, 0.0);
    GSL_SET_COMPLEX(&z0, 0.0, 0.0);
    gsl_blas_zgemm(CblasNoTrans, CblasNoTrans, z1, a->m, b->m, z0, r->m);

    return 1;
}

static int
m_complex(lua_State *L)
{
    int sl, sr;
    
    qlua_checkindex2(L, 1, "matrix size", &sl, &sr);
    qlua_newMatComplex(L, sl, sr);

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
    { "qr",                   mc_qr         },
    { "det",                  mc_det        },
	{ "solve",                mc_solve      },
    { NULL, NULL}
};

static const luaL_Reg mtMatComplex[] = {
    { "__tostring",     mc_fmt    },
    { "__gc",           mc_gc     },
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
    static const QLUA_Op2 ops[] = {
        {qlua_add_table, qMatReal,    qMatReal,    rm_add_rm },
        {qlua_add_table, qMatComplex, qMatComplex, cm_add_cm },
        {qlua_sub_table, qMatReal,    qMatReal,    rm_sub_rm },
        {qlua_sub_table, qMatComplex, qMatComplex, cm_sub_cm },
        {qlua_mul_table, qMatReal,    qMatReal,    rm_mul_rm },
        {qlua_mul_table, qReal,       qMatReal,    r_mul_rm  },
        {qlua_mul_table, qMatReal,    qReal,       rm_mul_r  },
        {qlua_mul_table, qMatComplex, qMatComplex, cm_mul_cm },
        {qlua_mul_table, qReal,       qMatComplex, r_mul_cm  },
        {qlua_mul_table, qComplex,    qMatComplex, c_mul_cm  },
        {qlua_mul_table, qMatComplex, qComplex,    cm_mul_c  },
        {qlua_mul_table, qMatComplex, qReal,       cm_mul_r  },
        {qlua_div_table, qMatReal,    qReal,       rm_div_r  },
        {qlua_div_table, qMatComplex, qComplex,    cm_div_c  },
        {qlua_div_table, qMatComplex, qReal,       cm_div_r  },
        {NULL,           qOther,      qOther,      NULL      }
    };

    gsl_set_error_handler_off();
    luaL_register(L, matrix_ns,      fMatrix);
    qlua_metatable(L, mtnMatReal,    mtMatReal,          qMatReal);
    qlua_metatable(L, opMatReal,     MatRealMethods,     qNoType);
    qlua_metatable(L, mtnMatComplex, mtMatComplex,       qMatComplex);
    qlua_metatable(L, opMatComplex,  MatComplexMethods,  qNoType);
    qlua_reg_op2(ops);

    return 0;
}

int
fini_matrix(lua_State *L)
{
    return 0;
}
