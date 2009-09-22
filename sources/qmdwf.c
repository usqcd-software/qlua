#include <qlua.h>                                                    /* DEPS */
#include <qmdwf.h>                                                   /* DEPS */
#include <qcomplex.h>                                                /* DEPS */
#include <lattice.h>                                                 /* DEPS */
#include <latcolmat.h>                                               /* DEPS */
#include <latdirferm.h>                                              /* DEPS */
#include <latdirprop.h>                                              /* DEPS */
#define QOP_MDWF_DEFAULT_PRECISION QDP_Precision
#include <qop-mdwf3.h>
#include <qmp.h>

/* NB: MDWF operator definition follows the docs, must agree with Chroma too */

static char mtnMDWF[] = "qcd.mdwf";

typedef struct {
    struct QOP_MDWF_State       *state;
    struct QOP_MDWF_Parameters  *params;
    struct QOP_MDWF_Gauge       *gauge;
    int                          Ls;
} mMDWF;

static mMDWF *
qlua_newMDWF(lua_State *L)
{
    mMDWF *c = lua_newuserdata(L, sizeof (mMDWF));
    
    c->state = 0;
    c->params = 0;
    c->gauge = 0;
    c->Ls = 0;
    luaL_getmetatable(L, mtnMDWF);
    lua_setmetatable(L, -2);

    return c;
}

static mMDWF *
qlua_checkMDWF(lua_State *L, int idx)
{
    void *v = luaL_checkudata(L, idx, mtnMDWF);

    luaL_argcheck(L, v != 0, idx, "qcd.MDWF expected");
    
    return v;
}

static int
q_MDWF_fmt(lua_State *L)
{
    char fmt[72];
    mMDWF *c = qlua_checkMDWF(L, 1);
    
    sprintf(fmt, "MDWF[%d]", c->Ls);
    lua_pushstring(L, fmt);

    return 1;
}

static int
q_MDWF_gc(lua_State *L)
{
    mMDWF *c = qlua_checkMDWF(L, 1);

    if (c->gauge) {
        QOP_MDWF_free_gauge(&c->gauge);
        c->gauge = 0;
    }
    if (c->params) {
        QOP_MDWF_free_parameters(&c->params);
        c->params = 0;
    }
    if (c->state) {
        QOP_MDWF_fini(&c->state);
        c->state = 0;
    }
    c->Ls = 0;

    return 0;
}

static int
q_MDWF_close(lua_State *L)
{
    mMDWF *c = qlua_checkMDWF(L, 1);

    if (c->gauge) {
        QOP_MDWF_free_gauge(&c->gauge);
        c->gauge = 0;
    }
    if (c->params) {
        QOP_MDWF_free_parameters(&c->params);
        c->params = 0;
    }
    if (c->state) {
        QOP_MDWF_fini(&c->state);
        c->state = 0;
    }
    c->Ls = 0;

    return 0;
}

static double
f5_reader(const int pos[QOP_MDWF_DIM + 1], int c, int d, int re_im, void *env)
{
    QLA_DiracFermion **in = env;
    int idx = QDP_index(pos);
    QLA_Complex *z = &QLA_elem_D(in[pos[QOP_MDWF_DIM]][idx], c, d);

    if (re_im == 0) {
        return QLA_real(*z);
    } else {
        return QLA_imag(*z);
    }
}

static void
f5_writer(const int pos[QOP_MDWF_DIM + 1], int c, int d, int re_im,
          double v, void *env)
{
    QLA_DiracFermion **out = env;
    int idx = QDP_index(pos);
    
    if (re_im == 0) {
        QLA_real(QLA_elem_D(out[pos[QOP_MDWF_DIM]][idx], c, d)) = v;
    } else {
        QLA_imag(QLA_elem_D(out[pos[QOP_MDWF_DIM]][idx], c, d)) = v;
    }
}

static int
MDWF_op(lua_State *L,
        const char *name,
        int (*op)(struct QOP_D3_MDWF_Fermion *result,
                  const struct QOP_MDWF_Parameters *params,
                  const struct QOP_D3_MDWF_Gauge *gauge,
                  const struct QOP_D3_MDWF_Fermion *fermion))
{
    mMDWF *c = qlua_checkMDWF(L, 1);
    int i;
    int len;
    QDP_DiracFermion **in;
    QLA_DiracFermion **x_in;
    struct QOP_MDWF_Fermion *m_in;
    QDP_DiracFermion **out;
    QLA_DiracFermion **x_out;
    struct QOP_MDWF_Fermion *m_out;

    luaL_checktype(L, 2, LUA_TTABLE);
    len = lua_objlen(L, 1);
    if (len != c->Ls) 
        return luaL_error(L, "%s(): wrong size of s-dimension", name);
    in = qlua_malloc(L, c->Ls * sizeof (QDP_DiracFermion *));
    x_in = qlua_malloc(L, c->Ls * sizeof (QLA_DiracFermion *));
    out = qlua_malloc(L, c->Ls * sizeof (QDP_DiracFermion *));
    x_out = qlua_malloc(L, c->Ls * sizeof (QLA_DiracFermion *));
    for (i = 0; i < c->Ls; i++) {
        in[i] = QDP_create_D();
        lua_pushnumber(L, i + 1); /* [sic] lua indexing */
        lua_gettable(L, 1);
        QDP_D_eq_D(in[i], qlua_checkLatDirFerm(L, -1)->ptr, QDP_all);
        x_in[i] = QDP_expose_D(in[i]);
    }
    lua_pop(L, 1);
    lua_createtable(L, c->Ls, 0);
    for (i = 0; i < c->Ls; i++) {
        mLatDirFerm *x = qlua_newLatDirFerm(L);
        lua_rawseti(L, -1, i + 1); /* [sic] lua indexing */
        out[i] = x->ptr;
        x_out[i] = QDP_expose_D(out[i]);
    }

    lua_gc(L, LUA_GCCOLLECT, 0);
    if (QOP_MDWF_import_fermion(&m_in, c->state, f5_reader, x_in))
        return luaL_error(L, "%s(): input conversion failed", name);
    if (QOP_MDWF_allocate_fermion(&m_out, c->state))
        return luaL_error(L, "%s(): output allocation failed", name);
    if (op(m_out, c->params, c->gauge, m_in))
        return luaL_error(L, "%s(): operator failed");
    if (QOP_MDWF_export_fermion(f5_writer, x_out, m_out))
        return luaL_error(L, "%s(): output conversion failed", name);

    for (i = 0; i < c->Ls; i++) {
        QDP_reset_D(in[i]);
        QDP_destroy_D(in[i]);
        QDP_reset_D(out[i]);
    }
    qlua_free(L, x_out);
    qlua_free(L, out);
    qlua_free(L, x_in);
    qlua_free(L, in);
    
    return 1;
}

static int
q_MDWF_D(lua_State *L)
{
    return MDWF_op(L, "MDWF:D", QOP_MDWF_DDW_operator);
}

static int
q_MDWF_Dx(lua_State *L)
{
    return MDWF_op(L, "MDWF:Dx", QOP_MDWF_DDW_operator_conjugated);
}

static int
q_MDWF_make_solver(lua_State *L)
{
    return luaL_error(L, "XXX MDWF:solver not implemented");
}

static int
q_MDWF_make_mixed_solver(lua_State *L)
{
    return luaL_error(L, "XXX MDWF:mixed_solver not implemented");
}

typedef struct {
    int lattice[QOP_MDWF_DIM + 1];
    int network[QOP_MDWF_DIM];
    QLA_Complex bf[QOP_MDWF_DIM];
    QLA_ColorMatrix *uf[QOP_MDWF_DIM];
} QMArgs;

static void
q_mdwf_sublattice(int lo[], int hi[], const int node[], void *env)
{
    QMArgs *args = env;
    int i;

    for (i = 0; i < QOP_MDWF_DIM; i++) {
        lo[i] = (args->lattice[i] * node[i]) / args->network[i];
        hi[i] = (args->lattice[i] * (node[i] + 1)) / args->network[i];
    }
}

static void
get_vector(int v[], int def, int dim, const int d[])
{
    int i;

    for (i = 0; i < dim && i < QOP_MDWF_DIM; i++)
        v[i] = d[i];
    for (; i < QOP_MDWF_DIM; i++)
        v[i] = def;
}

static double
q_MDWF_u_reader(int d, const int p[], int a, int b, int re_im, void *env)
{
    QLA_Complex z;
    QMArgs *args = env;
    int i = QDP_index(p);

    if (p[d] == (args->lattice[d] - 1)) {
        QLA_c_eq_c_times_c(z, args->bf[d], QLA_elem_M(args->uf[d][i], a, b));
    } else {
        QLA_c_eq_c(z, QLA_elem_M(args->uf[d][i], a, b));
    }

    if (re_im == 0)
        return QLA_real(z);
    else
        return QLA_imag(z);
}

/* XXX handle all kinds of domain wall */
static int
q_MDWF(lua_State *L)
{
    mMDWF *c = qlua_newMDWF(L);
    int Ls = luaL_checkint(L, 2);
    QDP_ColorMatrix *UF[QOP_MDWF_DIM];
    int i;
    QMArgs args;
    int node[QOP_MDWF_DIM];

    /* collect BC factors */
    if (lua_gettop(L) == 4) { /* (U, Ls, bc), c */
        for (i = 0; i < QOP_MDWF_DIM; i++) {
            lua_pushnumber(L, i + 1);
            lua_gettable(L, 3);
            switch (qlua_gettype(L, -1)) {
            case qReal:
                QLA_c_eq_r_plus_ir(args.bf[i], luaL_checknumber(L, -1), 0);
                break;
            case qComplex:
                QLA_c_eq_c(args.bf[i], *qlua_checkComplex(L, -1));
                break;
            default:
                luaL_error(L, "bad MDWF boundary condition type");
            }
            lua_pop(L, 1);
        }
    } else {
        for (i = 0; i < QOP_MDWF_DIM; i++)
            QLA_c_eq_r_plus_ir(args.bf[i], 1.0, 0.0);
    }
    /* check argument types */
    luaL_checktype(L, 1, LUA_TTABLE);
    lua_gc(L, LUA_GCCOLLECT, 0);
    /* check restrictions */
    if (qRank != QOP_MDWF_DIM)
        return luaL_error(L, "MDWF not implemented for #L=%d", qRank);
    if (QDP_Nc != QOP_MDWF_COLORS)
        return luaL_error(L, "MDWF not implemented for Nc=%d", QDP_Nc);
    if (QDP_Ns != QOP_MDWF_FERMION_DIM)
        return luaL_error(L, "MDWF not implemented for Ns=%d", QDP_Ns);
    /* get the gauge field from arguments */
    for (i = 0; i < QOP_MDWF_DIM; i++) {
        UF[i] = QDP_create_M();
        lua_pushnumber(L, i + 1); /* [sic] lua idexing */
        lua_gettable(L, 1);
        QDP_M_eq_M(UF[i], qlua_checkLatColMat(L, -1)->ptr, QDP_all);
        args.uf[i] = QDP_expose_M(UF[i]);
        lua_pop(L, 1);
    }
    /* init the state */
    get_vector(args.network, 1, QMP_get_logical_number_of_dimensions(),
               QMP_get_logical_dimensions());
    get_vector(node, 0, QMP_get_logical_number_of_dimensions(),
               QMP_get_logical_coordinates());
    QDP_latsize(args.lattice);
    args.lattice[QOP_MDWF_DIM] = Ls;
    c->Ls = Ls;
    if (QOP_MDWF_init(&c->state, args.lattice, args.network, node,
                      QMP_is_primary_node(), q_mdwf_sublattice, &args))
        return luaL_error(L, "MDWF_init() failed");

    /* import gauge field */
    if (QOP_MDWF_import_gauge(&c->gauge, c->state, q_MDWF_u_reader, &args))
        return luaL_error(L, "MDWF_import_gauge() failed");

    /* free gauge field */
    for (i = 0; i < QOP_MDWF_DIM; i++) {
        QDP_reset_M(UF[i]);
        QDP_destroy_M(UF[i]);
    }

    return 1;
}

static struct luaL_Reg mtMDWF[] = {
    { "__tostring",    q_MDWF_fmt },
    { "__gc",          q_MDWF_gc },
    { "close",         q_MDWF_close },
    { "D",             q_MDWF_D },
    { "Dx",            q_MDWF_Dx },
    { "solver",        q_MDWF_make_solver },
    { "mixed_solver",  q_MDWF_make_mixed_solver },
    { NULL,            NULL }
};

static struct luaL_Reg fMDWF[] = {
    /* XXX allocators for all kinds of DWF */
    { "MDWF",      q_MDWF }, /* XXX this becomes an inner allocator */
    { NULL,        NULL }
};

int
init_mdwf(lua_State *L)
{
    luaL_register(L, qcdlib, fMDWF);
    qlua_metatable(L, mtnMDWF, mtMDWF);
    return 0;
}

int
fini_mdwf(lua_State *L)
{
    return 0;
}
