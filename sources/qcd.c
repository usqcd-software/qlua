#include <qlua.h>
#include <stdlib.h>
#include <qcd.h>

static int qRank = 0;
static int *qDim = NULL;

/* lattice integers */
const char *mtnLatInt = "qcd.lattice.int";

mLatInt *
qlua_newLatInt(lua_State *L)
{
    QDP_Int *v = QDP_create_I();
    mLatInt *hdr;
    
    if (v == 0) {
        lua_gc(L, LUA_GCCOLLECT, 0);
        v = QDP_create_I();
        if (v == 0)
            luaL_error(L, "not enough memory (QDP_Int)");
    }
    hdr = lua_newuserdata(L, sizeof (mLatInt));
    hdr->ptr = v;
    luaL_getmetatable(L, mtnLatInt);
    lua_setmetatable(L, -2);

    return hdr;
}

mLatInt *
qlua_checkLatInt(lua_State *L, int idx)
{
    void *v = luaL_checkudata(L, idx, mtnLatInt);

    luaL_argcheck(L, v != 0, idx, "qcd.LatInt expected");
    
    return v;
}

static int
qLatInt_fmt(lua_State *L)
{
    char fmt[72];
    mLatInt *b = qlua_checkLatInt(L, 1);

    sprintf(fmt, "LatInt(%p)", b->ptr);
    lua_pushstring(L, fmt);

    return 1;
}

static int
qLatInt_gc(lua_State *L)
{
    mLatInt *b = qlua_checkLatInt(L, 1);

    QDP_destroy_I(b->ptr);
    b->ptr = 0;

    return 0;
}

static struct luaL_reg mtLatInt[] = {
    { "__tostring", qLatInt_fmt },
    { "__gc",       qLatInt_gc },
    { NULL,         NULL}
};

/* lattice definition */
static int
qcd_lattice(lua_State *L)
{
    int r, i;

    if (qRank != 0)
        return luaL_error(L, "redefining the lattice");
    r = lua_gettop(L);
    if (r <= 0)
        return luaL_error(L, "Bad lattice rank");
    qRank = r;
    qDim = qlua_malloc(L, r * sizeof (int));
    for (i = 0; i < r; i++) {
        qDim[i] = luaL_checkint(L, i + 1);
    }
    QDP_set_latsize(qRank, qDim);
    if (QDP_create_layout()) {
        return luaL_error(L, "can not create lattice");
    }
    
    return 0;
}

static int
qcd_rank(lua_State *L)
{
    lua_pushnumber(L, qRank);
    return 1;
}

static int
qcd_dim(lua_State *L)
{
    int d = luaL_checkint(L, 1);

    if ((d < 0) || (d >= qRank))
        return luaL_error(L, "bad dimension number");
    lua_pushnumber(L, qDim[d]);

    return 1;
}

/* pcoords */
static int pcoord_d = -1;
static void
pcoord_set(QLA_Int *dst, int coords[])
{
    *dst = coords[pcoord_d];
}

static int
qcd_pcoord(lua_State *L)
{
    int d = luaL_checkint(L, 1);
    mLatInt *v = qlua_newLatInt(L);
    
    if ((d < 0) || (d >= qRank))
        return luaL_error(L, "coordinate out of range");
    
    /* YYY global state */
    pcoord_d = d;
    QDP_I_eq_func(v->ptr, pcoord_set, QDP_all);

    return 1;
}

static struct luaL_Reg fQcd[] = {
    { "lattice", qcd_lattice },
    { "rank",    qcd_rank },
    { "dim",     qcd_dim },
    { "pcoord",  qcd_pcoord },
    { NULL, NULL}
};

int
init_qcd(lua_State *L)
{
    luaL_register(L, qcdlib, fQcd);
    qlua_metatable(L, mtnLatInt, mtLatInt);

    return 0;
}
