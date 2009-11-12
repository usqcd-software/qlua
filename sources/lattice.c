#include <qlua.h>                                                    /* DEPS */
#include <qcomplex.h>                                                /* DEPS */
#include <qvector.h>                                                 /* DEPS */
#include <lattice.h>                                                 /* DEPS */
#include <latint.h>                                                  /* DEPS */
#include <latcomplex.h>                                              /* DEPS */
#include <string.h>
#include <qmp.h>
#include <math.h>

/* NB: This code works only for a single lattice */

QDP_Subset *qCurrent;
const char opLattice[] = "lattice.ops";

typedef struct {
    int rank;
    int dim[1];
} mLattice;

static const char mtnLattice[] = "qcd.lattice";
int qRank = 0;
int *qDim = NULL;

static int
q_L_fmt(lua_State *L)
{
    luaL_Buffer b;
    char fmt[72];
    char *sep;
    int i;

    luaL_buffinit(L, &b);
    luaL_addstring(&b, "lattice(");
    for (i = 0, sep = ""; i < qRank; i++, sep = ", ") {
        sprintf(fmt, "%s%d", sep, qDim[i]);
        luaL_addstring(&b, fmt);
    }
    luaL_addstring(&b, ")");
    luaL_pushresult(&b);
    
    return 1;
}

static int
q_L_get(lua_State *L)
{
    switch (qlua_gettype(L, 2)) {
    case qReal: {
        int d = luaL_checkint(L, 2);

        if ((d >= 0) && (d < qRank)) {
            lua_pushnumber(L, qDim[d]);
            return 1;
        }
        break;
    }
    case qString:
        return qlua_lookup(L, 2, opLattice);
    }

    return qlua_badindex(L, "Lattice");
}

static int
q_L_dim(lua_State *L)
{
    lua_pushnumber(L, qRank);

    return 1;
}

static int pcoord_d = -1; /* YYY global state */
static void
pcoord_set(QLA_Int *dst, int coords[])
{
    *dst = coords[pcoord_d];
}

static int
q_pcoord(lua_State *L)
{
    int d = luaL_checkint(L, 2);
    mLatInt *v = qlua_newLatInt(L);
    
    if ((d < 0) || (d >= qRank))
        return luaL_error(L, "coordinate out of range");
    
    /* YYY global state */
    CALL_QDP(L);
    pcoord_d = d;
    QDP_I_eq_func(v->ptr, pcoord_set, *qCurrent);

    return 1;
}

static struct {
    int *s;
    int *p;
} PW_arg;

static void
pw_simple(QLA_Complex *dst, int coord[])
{
    int i;
    double ph;

    for (ph = 0, i = 0; i < qRank; i++) {
        double d = (coord[i]-PW_arg.s[i]) * ((double)PW_arg.p[i]);
        ph += 2 * M_PI * d / qDim[i];
    }
    *dst = QLA_cexpi(ph);
}

static int
q_planewave(lua_State *L)
{
    int *s = qlua_checklatcoord(L, 2);
    int *p = qlua_checkintarray(L, 3, qRank, NULL);
    mLatComplex *w = qlua_newLatComplex(L);

    /* YYY global state */
    PW_arg.s = s;
    PW_arg.p = p;
    CALL_QDP(L);
    QDP_C_eq_func(w->ptr, pw_simple, *qCurrent);

    PW_arg.s = 0;
    PW_arg.p = 0;
    qlua_free(L, s);
    qlua_free(L, p);
    return 1;
}

static void
mp_full(lua_State *L, mLatComplex *v, int *src, int *p)
{
    int *coord = qlua_malloc(L, qRank * sizeof (int));
    QLA_Complex *r = qlua_newComplex(L);
    int idx;
    double c2pt_mp[2]; /* re, im */
    QLA_Complex *z;

    CALL_QDP(L);
    c2pt_mp[0] = 0;
    c2pt_mp[1] = 0;
    z = QDP_expose_C(v->ptr);
    for (idx = 0; idx <  QDP_sites_on_node ; idx++) {
        double ph, re, im, c, s;
        int i;
        QDP_get_coords(coord, QDP_this_node, idx);
        for (i = 0, ph = 0.; i < qRank; i++) {
            double d = p[i] * (coord[i] - src[i]);
            ph += 2 * M_PI * d / qDim[i];
        }
        re   = QLA_real(z[idx]);
        im   = QLA_imag(z[idx]);
        c    = cos(ph);
        s    = sin(ph);
        c2pt_mp[0]    += re * c - im * s;
        c2pt_mp[1]    += re * s + im * c;
    }
    QDP_reset_C(v->ptr);
    QMP_sum_double_array(c2pt_mp, 2);
    QLA_c_eq_r_plus_ir(*r, c2pt_mp[0], c2pt_mp[1]);

    qlua_free(L, coord);
}

static void
mp_axis(lua_State *L, mLatComplex *v, int *src, int *p, int axis)
{
    int *coord = qlua_malloc(L, qRank * sizeof (int));
    double *c2pt = qlua_malloc(L, 2 * qDim[axis] * sizeof (double));
    mVecComplex *r = qlua_newVecComplex(L, qDim[axis]);
    int i, idx;
    QLA_Complex *z;

    CALL_QDP(L);
    p[axis] = 0; /* avoid phases along the axis */
    r->size = qDim[axis];
    for (i = 0; i < 2 * qDim[axis]; i++)
        c2pt[i] = 0;
    z = QDP_expose_C(v->ptr);
    for (idx = 0; idx <  QDP_sites_on_node ; idx++) {
        double ph, re, im, c, s;
        QDP_get_coords(coord, QDP_this_node, idx);
        for (i = 0, ph = 0.; i < qRank; i++) {
            double d = p[i] * (coord[i] - src[i]);
            ph += 2 * M_PI * d / qDim[i];
        }
        re   = QLA_real(z[idx]);
        im   = QLA_imag(z[idx]);
        c    = cos(ph);
        s    = sin(ph);
        c2pt[2 * coord[axis] + 0]    += re * c - im * s;
        c2pt[2 * coord[axis] + 1]    += re * s + im * c;
    }
    QDP_reset_C(v->ptr);
    QMP_sum_double_array(c2pt, 2 * qDim[axis]);
    for (i = 0; i < qDim[axis]; i++)
        QLA_c_eq_r_plus_ir(r->val[i], c2pt[2 * i + 0], c2pt[2 * i + 1]);

    qlua_free(L, coord);
    qlua_free(L, c2pt);
}

static int
q_momentum_project(lua_State *L)
{
    mLatComplex *v = qlua_checkLatComplex(L, 2);
    int *s = qlua_checklatcoord(L, 3);
    int *p = qlua_checkintarray(L, 4, qRank, NULL);

    switch (lua_type(L, 5)) {
    case LUA_TNUMBER: {
        int axis = luaL_checkint(L, 5);
        if (axis < 0 || axis >= qRank)
            return luaL_error(L, "axis out of range");
        mp_axis(L, v, s, p, axis);
        break;
    }
    case LUA_TNONE:
        mp_full(L, v, s, p);
    break;
    default:
        return luaL_error(L, "axis must be a number");
    }
    qlua_free(L, s);
    qlua_free(L, p);
    return 1;
}

static void
set_Nd(lua_State *L)
{
    lua_getglobal(L, qcdlib);
    lua_pushnumber(L, qRank);
    lua_setfield(L, -2, "Nd");
    lua_pop(L, 1);
}

static int
q_lattice(lua_State *L)
{
    int r, i;
    mLattice *res;

    if (qRank != 0)
        return luaL_error(L, "redefining the lattice");

    luaL_checktype(L, 1, LUA_TTABLE);
    r = lua_objlen(L, 1);
    if (r <= 0)
        return luaL_error(L, "Bad lattice rank");
    qRank = r;
    qDim = qlua_malloc(L, r * sizeof (int));
    for (i = 0; i < r; i++) {
        lua_pushnumber(L, i + 1);
        lua_gettable(L, 1);
        qDim[i] = luaL_checkint(L, -1);
    }
    CALL_QDP(L);
    QDP_set_latsize(qRank, qDim);
    if (QDP_create_layout()) {
        return luaL_error(L, "can not create lattice");
    }
    qCurrent = &QDP_all;

    res = lua_newuserdata(L, sizeof (mLattice) + (qRank - 1) * sizeof (int));
    res->rank = qRank;
    for (i = 0; i < qRank; i++) {
        res->dim[i] = qDim[i];
    }
    luaL_getmetatable(L, mtnLattice);
    lua_setmetatable(L, -2);

    set_Nd(L);
    
    return 1;
}

int *
qlua_intarray(lua_State *L, int n, int *dim)
{
    int d, i;
    int *idx;

    if (lua_type(L, n) != LUA_TTABLE)
        return NULL;

    d = lua_objlen(L, n);
    idx = qlua_malloc(L, d * sizeof (int));
    for (i = 0; i < d; i++) {
        lua_pushnumber(L, i + 1);
        lua_gettable(L, n);
        if (lua_type(L, -1) != LUA_TNUMBER) {
            qlua_free(L, idx);
            return NULL;
        }
        idx[i] = luaL_checkint(L, -1);
        lua_pop(L, 1);
    }
    *dim = d;
    return idx;
}

int *
qlua_checkintarray(lua_State *L, int n, int dim, int *out_dim)
{
    int d_dim;
    int *idx = qlua_intarray(L, n, &d_dim);

    if (idx == 0)
        luaL_error(L, "table of integers expected");

    if (out_dim)
        *out_dim = d_dim;
    else if (d_dim != dim)
        luaL_error(L, "table of integer has wrong size");

    return idx;
}

int *
qlua_latcoord(lua_State *L, int n)
{
    int d, i;
    int *idx;

    luaL_checktype(L, n, LUA_TTABLE);
    d = lua_objlen(L, n);
    if (d != qRank) {
        return NULL;
    }
    idx = qlua_malloc(L, d * sizeof (int));
    for (i = 0; i < d; i++) {
        lua_pushnumber(L, i + 1);
        lua_gettable(L, n);
        idx[i] = luaL_checkint(L, -1);
        lua_pop(L, 1);
        if ((idx[i] < 0) || (idx[i] >= qDim[i])) {
            qlua_free(L, idx);        
            return NULL;
        }
    }
    
    return idx;
}

static const char *bad_lat_coords = "bad lattice coordinates";

void
qlua_verifylatcoord(lua_State *L, int *coord)
{
    int i;

    for (i = 0; i < qRank; i++) {
        if ((coord[i] < 0) || (coord[i] >= qDim[i]))
            luaL_error(L, bad_lat_coords);
    }
}

int *
qlua_checklatcoord(lua_State *L, int n)
{
    int *idx = qlua_latcoord(L, n);

    if (idx == 0)
        luaL_error(L, bad_lat_coords);

    return idx;
}

QDP_Shift
qlua_checkShift(lua_State *L, int idx)
{
    int d = luaL_checkint(L, idx);

    if ((d < 0) || (d >= qRank))
        luaL_error(L, "bad shift dimension");

    return QDP_neighbor[d];
}

QDP_ShiftDir
qlua_checkShiftDir(lua_State *L, int idx)
{
    static const struct {
        char *name;
        QDP_ShiftDir dir;
    } t[] = {
        { "from_forward",  QDP_forward },
        { "from_backward", QDP_backward },
        { "to_forward",    QDP_backward },
        { "to_backward",   QDP_forward },
        { NULL,            QDP_forward }
    };
    int i;
    const char *d = luaL_checkstring(L, idx);

    for (i = 0; t[i].name; i++) {
        if (strcmp(d, t[i].name) == 0)
            return t[i].dir;
    }
    luaL_error(L, "bad shift direction");
    /* NEVER HAPPENS */
    return QDP_forward;
}

static int
q_network(lua_State *L)
{
    int n = QMP_get_logical_number_of_dimensions();
    const int *d = QMP_get_logical_dimensions();
    int i;

    lua_pushnumber(L, QMP_get_number_of_nodes());
    lua_createtable(L, n, 0);
    for (i = 0; i < n; i++) {
        lua_pushnumber(L, d[i]);
        lua_rawseti(L, -2, i + 1);
    }

    return 2;
}

static struct luaL_Reg LatticeMethods[] = {
    { "pcoord",           q_pcoord },
    { "planewave",        q_planewave },
    { "momentum_project", q_momentum_project },
    { NULL,           NULL },
};

static struct luaL_Reg mtLattice[] = {
    { "__tostring",   q_L_fmt   },
    { "__index",      q_L_get   },
    { "__len",        q_L_dim   },
    { NULL,           NULL      }
};

static struct luaL_Reg fLattice[] = {
    { "lattice",            q_lattice },
    { "network",            q_network },
    { NULL, NULL}
};

int
init_lattice(lua_State *L)
{
    luaL_register(L, qcdlib, fLattice);
    qlua_metatable(L, mtnLattice, mtLattice);
    qlua_metatable(L, opLattice, LatticeMethods);

    set_Nd(L);

    return 0;
}

int
fini_lattice(lua_State *L)
{
    qlua_free(L, qDim);
    qDim = 0;
    qRank = 0;
    return 0;
}
