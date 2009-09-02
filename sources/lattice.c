#include <qlua.h>                                                    /* DEPS */
#include <lattice.h>                                                 /* DEPS */
#include <latint.h>                                                  /* DEPS */
#include <string.h>

/* NB: This code works only for a single lattice */

QDP_Subset qCurrent;
const char opLattice[] = "lattice.ops";

typedef struct {
    int rank;
    int dim[1];
} mLattice;

static const char mtnLattice[] = "qcd.lattice";
static int qRank = 0;
static int *qDim = NULL;

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

    return luaL_error(L, "bad index");
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
    pcoord_d = d;
    QDP_I_eq_func(v->ptr, pcoord_set, qCurrent);

    return 1;
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
    QDP_set_latsize(qRank, qDim);
    if (QDP_create_layout()) {
        return luaL_error(L, "can not create lattice");
    }
    qCurrent = QDP_all;

    res = lua_newuserdata(L, sizeof (mLattice) + (qRank - 1) * sizeof (int));
    res->rank = qRank;
    for (i = 0; i < qRank; i++) {
        res->dim[i] = qDim[i];
    }
    luaL_getmetatable(L, mtnLattice);
    lua_setmetatable(L, -2);
    
    return 1;
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
        if ((idx[i] < 0) || (idx[i] >= qDim[i])) {
            qlua_free(L, idx);        
            return NULL;
        }
    }
    
    return idx;
}

int *
qlua_checklatcoord(lua_State *L, int n)
{
    int *idx = qlua_latcoord(L, n);

    if (idx == 0)
        luaL_error(L, "bad lattice coordinates");

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

#if 0 /* XXX */

static int
q_dims(lua_State *L)
{
    int i;

    lua_createtable(L, qRank, 0);
    for (i = 0; i < qRank; i++) {
        lua_pushnumber(L, qDim[i]);
        lua_rawseti(L, -2, i + 1);
    }
    return 1;
}
#endif

static struct luaL_Reg LatticeMethods[] = {
    { "pcoord",       q_pcoord },
    { NULL,           NULL },
};

static struct luaL_Reg mtLattice[] = {
    { "__tostring",   q_L_fmt   },
    { "__index",      q_L_get   },
    { "__len",        q_L_dim   },
    { NULL,           NULL      }
};

static struct luaL_Reg fLattice[] = {
    { "lattice", q_lattice },
    { NULL, NULL}
};

int
init_lattice(lua_State *L)
{
    luaL_register(L, qcdlib, fLattice);
    qlua_metatable(L, mtnLattice, mtLattice);
    qlua_metatable(L, opLattice, LatticeMethods);

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
