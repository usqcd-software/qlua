#include "qlua.h"                                                    /* DEPS */
#include "qcomplex.h"                                                /* DEPS */
#include "qvector.h"                                                 /* DEPS */
#include "lattice.h"                                                 /* DEPS */
#include "qlayout.h"                                                 /* DEPS */
#include "latsubset.h"                                               /* DEPS */
#include "latint.h"                                                  /* DEPS */
#include "latcomplex.h"                                              /* DEPS */
#include "qmp.h"
#include "qla_types.h"
#include <string.h>

/* NB: This code is not tested for multiple lattices */
const char opLattice[] = "lattice.ops";

mLattice *
qlua_checkLattice(lua_State *L, int idx)
{
    void *v = lua_touserdata(L, idx);

    if (qlua_qtype(L, idx) != qLattice)
        luaL_error(L, "lattice expected");

    return v;
}

mLattice *
qlua_ObjLattice(lua_State *L, int idx)
{
    qlua_getlattice(L, idx);
    return qlua_checkLattice(L, -1);
}

static int
q_L_fmt(lua_State *L)
{
    mLattice *S = qlua_checkLattice(L, 1);
    luaL_Buffer b;
    char fmt[72];
    char *sep;
    int i;

    luaL_buffinit(L, &b);
    luaL_addstring(&b, "lattice(");
    for (i = 0, sep = ""; i < S->rank; i++, sep = ", ") {
        sprintf(fmt, "%s%d", sep, S->dim[i]);
        luaL_addstring(&b, fmt);
    }
    luaL_addstring(&b, ")");
    luaL_pushresult(&b);
    
    return 1;
}

static int
q_L_gc(lua_State *L)
{
    mLattice *S = qlua_checkLattice(L, 1);

    if (S->neighbor_up)
        qlua_free(L, S->neighbor_up);
    S->neighbor_up = 0;

    if (S->neighbor_down)
        qlua_free(L, S->neighbor_down);
    S->neighbor_down = 0;

    if (S->net)
        qlua_free(L, S->net);
    S->net = 0;

    if (S->dim)
        qlua_free(L, S->dim);
    S->dim = 0;

    if (S->lss.cl > qss_last_static)
        QDP_destroy_subset(S->qss);
    S->lss.cl = qss_none;

    if (S->lss.mask)
        QDP_destroy_I(S->lss.mask);
    S->lss.mask = NULL;

    if (S->none)
        QDP_destroy_subset(S->none);
    S->none = NULL;
    /* standard subsets are freed by the lattice */

    if (S->lat)
        QDP_destroy_lattice(S->lat);
    S->lat = NULL;

    return 0;
}

static int
q_L_get(lua_State *L)
{
    mLattice *S = qlua_checkLattice(L, 1);

    switch (qlua_qtype(L, 2)) {
    case qReal: {
        int d = luaL_checkint(L, 2);

        if ((d >= 0) && (d < S->rank)) {
            lua_pushnumber(L, S->dim[d]);
            return 1;
        }
        break;
    }
    case qString:
        return qlua_lookup(L, 2, opLattice);
    default:
        break;
    }

    return qlua_badindex(L, "Lattice");
}

static int
q_L_dim(lua_State *L)
{
    mLattice *S = qlua_checkLattice(L, 1);

    lua_pushnumber(L, S->rank);

    return 1;
}

static void
pcoord_set(QLA_Int *dst, int coords[], void *env)
{
    int *pcoord_d = env;
    *dst = coords[*pcoord_d];
}

static int
q_pcoord(lua_State *L)
{
    mLattice *S = qlua_checkLattice(L, 1);
    QDP_Subset *qCurrent = S->qss;
    int d = luaL_checkint(L, 2);
    mLatInt *v = qlua_newLatInt(L, 1);
    
    if ((d < 0) || (d >= S->rank))
        return luaL_error(L, "coordinate out of range");
    
    CALL_QDP(L);
    QDP_I_eq_funca(v->ptr, pcoord_set, &d, *qCurrent);

    return 1;
}

static void
add_default(lua_State *L, const char *key, int value)
{
    lua_pushnumber(L, value);
    lua_setfield(L, -2, key);
}

static int
set_default(lua_State *L, const char *key, int def)
{
    int v = def;

    lua_getfield(L, 2, key);
    if (qlua_qtype(L, -1) == qReal) {
        v = luaL_checkint(L, -1);
    }
    lua_pop(L, 1);
    return v;
}

static int
q_defaults(lua_State *L)
{
    static const char *colors_key = "colors";

    mLattice *S = qlua_checkLattice(L, 1);

    switch (lua_gettop(L)) {
    case 1: {
        lua_createtable(L, 0, 1);
        add_default(L, colors_key, S->nc);
        return 1;
    }
    case 2: {
        if (qlua_qtype(L, 2) == qTable) {
            int nc = set_default(L, colors_key, S->nc);
            if (nc > 0)
                S->nc = nc;
            else
                return luaL_error(L, "bad default number of colors %d", nc);
            return 1;
        }
        break;
    }
    default:
        break;
    }
    return luaL_error(L, "bad parameters for L:defaults()");
}

static int
q_latnet(lua_State *L)
{
    mLattice *S = qlua_checkLattice(L, 1);
    
    lua_createtable(L, S->rank, 0);
    for (int i = 0; i < S->rank; i++) {
        lua_pushnumber(L, S->net[i]);
        lua_rawseti(L, -2, i + 1);
    }

    return 1;
}

typedef struct {
    int *s;
    int *p;
    int rank;
    int *dim;
} PW_arg;

static void
pw_simple(QLA_Complex *dst, int coord[], void *env)
{
    int i;
    double ph;
    PW_arg *arg = env;

    for (ph = 0, i = 0; i < arg->rank; i++) {
        double d = (coord[i] - arg->s[i]) * ((double)arg->p[i]);
        ph += 2 * M_PI * d / arg->dim[i];
    }
    *dst = QLA_cexpi(ph);
}

static int
q_planewave(lua_State *L)
{
    mLattice *S = qlua_checkLattice(L, 1);
    QDP_Subset *qCurrent = S->qss;
    int *s = qlua_checklatcoord(L, 2, S);
    int *p = qlua_checkintarray(L, 3, S->rank, NULL);
    mLatComplex *w = qlua_newLatComplex(L, 1);
    PW_arg arg;

    arg.s = s;
    arg.p = p;
    arg.rank = S->rank;
    arg.dim = S->dim;
    CALL_QDP(L);
    QDP_D_C_eq_funca(w->ptr, pw_simple, &arg, *qCurrent);

    qlua_free(L, s);
    qlua_free(L, p);
    return 1;
}

static int
subset_none(QDP_Lattice *S, int *coord, void *arg)
{
    return 2;
}

static int
q_lattice(lua_State *L)
{
    int r, i;
    static int lat_id = 0;
    mLattice *S = lua_newuserdata(L, sizeof (mLattice));
    static const struct luaL_Reg mtLattice[] = {
        { "__tostring",   q_L_fmt        },
        { "__gc",         q_L_gc         },
        { "__index",      q_L_get        },
        { "__newindex",   qlua_nowrite   },
        { "__len",        q_L_dim        },
        { NULL,           NULL           }
    };
    
    if (lat_id) /* XXX need more work for multiple lattices */
        return luaL_error(L, "multiple lattices not supported yet");

    luaL_checktype(L, 1, LUA_TTABLE);
    r = lua_objlen(L, 1);
    if (r <= 0)
        return luaL_error(L, "Bad lattice rank");
	if (r > QLUA_MAX_LATTICE_RANK)
		return luaL_error(L, "latice rank is too large");
    S->rank = r;
    S->node = 0;
    S->neighbor_up = qlua_malloc(L, r * sizeof (int));
    S->neighbor_down = qlua_malloc(L, r * sizeof (int));
    S->net = qlua_malloc(L, r * sizeof (int));
    S->dim = qlua_malloc(L, r * sizeof (int));
    for (i = 0; i < r; i++) {
        lua_pushnumber(L, i + 1);
        lua_gettable(L, 1);
        S->dim[i] = qlua_checkint(L, -1, "lattice dim %d", i);
        lua_pop(L, 1);
    }
    CALL_QDP(L);
    S->lat = QDP_create_lattice(&qlua_layout, S, S->rank, S->dim);
    if (S->lat == 0)
        return luaL_error(L, "can not create lattice");
    S->all = QDP_all_L(S->lat);
    S->even = QDP_even_L(S->lat);
    S->odd = QDP_odd_L(S->lat);
    S->none = QDP_create_subset_L(S->lat, subset_none, 0, 0, 1);
    S->qss = &S->all;
    S->lss.cl = qss_all;
    S->lss.mask = NULL;
    S->id = lat_id;
    S->nc = 3;
    lat_id++;

    qlua_selftable(L, mtLattice, qLattice);
    lua_setmetatable(L, -2);
    
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
        idx[i] = qlua_checkint(L, -1, "array element %d", i + 1);
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
qlua_latcoord(lua_State *L, int n, mLattice *S)
{
    int d, i;
    int *idx;

    luaL_checktype(L, n, LUA_TTABLE);
    d = lua_objlen(L, n);
    if (d != S->rank) {
        return NULL;
    }
    idx = qlua_malloc(L, d * sizeof (int));
    for (i = 0; i < d; i++) {
        lua_pushnumber(L, i + 1);
        lua_gettable(L, n);
        idx[i] = qlua_checkint(L, -1, "lattice coord %d", i);
        lua_pop(L, 1);
        if ((idx[i] < 0) || (idx[i] >= S->dim[i])) {
            qlua_free(L, idx);        
            return NULL;
        }
    }
    
    return idx;
}

static const char *bad_lat_coords = "bad lattice coordinates";

void
qlua_verifylatcoord(lua_State *L, int *coord, mLattice *S)
{
    int i;

    for (i = 0; i < S->rank; i++) {
        if ((coord[i] < 0) || (coord[i] >= S->dim[i]))
            luaL_error(L, bad_lat_coords);
    }
}

int *
qlua_checklatcoord(lua_State *L, int n, mLattice *S)
{
    int *idx = qlua_latcoord(L, n, S);

    if (idx == 0)
        luaL_error(L, bad_lat_coords);

    return idx;
}

QDP_Shift
qlua_checkShift(lua_State *L, int idx, mLattice *S)
{
    int d = luaL_checkint(L, idx);

    if ((d < 0) || (d >= S->rank))
        luaL_error(L, "bad shift dimension");

    return QDP_neighbor_L(S->lat)[d];
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
        { NULL,            QDP_forward } /* never used */
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

static int
q_volume(lua_State *L)
{
  mLattice *S = qlua_checkLattice(L, 1);
  int vol = QDP_volume_L(S->lat);
  lua_pushnumber(L, vol);
  return 1;
}

static struct luaL_Reg LatticeMethods[] = {
    { "pcoord",           q_pcoord        },
    { "planewave",        q_planewave     },
    { "defaults",         q_defaults      },
    { "network",          q_latnet        },
    { "everywhere",       qlua_everywhere },
    { "volume",           q_volume        },
    { NULL,               NULL            }
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
    qlua_metatable(L, opLattice, LatticeMethods, qNoType);

    return 0;
}

int
fini_lattice(lua_State *L)
{
    return 0;
}
