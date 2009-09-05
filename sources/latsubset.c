#include <qlua.h>                                                    /* DEPS */
#include <lattice.h>                                                 /* DEPS */
#include <latsubset.h>                                               /* DEPS */
#include <string.h>

const char mtnLatSubset[] = "lattice.subset";
static const char LatSubsetStack[] = "lattice.subsetstack";

typedef struct {
    int heap;
    char *name;
    QDP_Subset *s;
} mLatSubset;

static mLatSubset *
qlua_newLatSubset(lua_State *L, const char *name)
{
    mLatSubset *v = lua_newuserdata(L, sizeof (mLatSubset));
    
    v->heap = 0;
    v->s = NULL;
    v->name = qlua_malloc(L, strlen(name) + 1);
    strcpy(v->name, name);
    luaL_getmetatable(L, mtnLatSubset);
    lua_setmetatable(L, -2);

    return v;
}

static mLatSubset *
qlua_checkLatSubset(lua_State *L, int idx)
{
    void *v = luaL_checkudata(L, idx, mtnLatSubset);

    luaL_argcheck(L, v != 0, idx, "lattice.Subset expected");

    return v;
}

static int
q_U_fmt(lua_State *L)
{
    mLatSubset *m = qlua_checkLatSubset(L, 1);
    char fmt[72];

    sprintf(fmt, "Subset[%s]", m->name);

    lua_pushstring(L, fmt);

    return 1;
}

static int
q_U_gc(lua_State *L)
{
    mLatSubset *v = qlua_checkLatSubset(L, 1);

    if (v->heap) {
        QDP_destroy_subset(v->s);
    }
    v->s = 0;
    v->heap = 0;
    if (v->name) {
        qlua_free(L, v->name);
    }
    v->name = 0;
    
    return 0;
}

static int
q_U_where(lua_State *L)
{
    mLatSubset *v = qlua_checkLatSubset(L, 1);
    int argc = lua_gettop(L) - 2;
    int resc;
    QDP_Subset old_set = qCurrent;

    qCurrent = *v->s;

    /* save old context from GC on the subset stack */
    luaL_getmetatable(L, LatSubsetStack);
    lua_createtable(L, 2, 0);
    lua_rawgeti(L, -2, 1);
    lua_rawseti(L, -2, 1);
    lua_rawgeti(L, -2, 2);
    lua_rawseti(L, -2, 2);
    lua_rawseti(L, -2, 2);
    lua_pushvalue(L, 1);
    lua_rawseti(L, -2, 1);
    lua_pop(L, 1);

    if (lua_pcall(L, argc, LUA_MULTRET, 0))
        return luaL_error(L, luaL_checkstring(L, -1));
    resc = lua_gettop(L) - 1;

    /* pop old subset from the subset stack */
    luaL_getmetatable(L, LatSubsetStack);
    lua_rawgeti(L, -1, 2);
    lua_rawgeti(L, -1, 1);
    lua_rawseti(L, -3, 1);
    lua_rawgeti(L, -1, 2);
    lua_rawseti(L, -3, 2);
    lua_pop(L, 2);

    qCurrent = old_set;

    return resc;
}

typedef struct {
    int axis;
    int pos;
} Subset_Arg;

static int
subset_func(int *coord, void *arg)
{
    Subset_Arg *x = arg;

    return coord[x->axis] != x->pos;
}

static int
q_subset(lua_State *L)
{
    int d = qlua_checkindex(L, 2, "axis", qRank);
    int p = qlua_checkindex(L, 2, "position", qDim[d]);
    char fmt[72];
    Subset_Arg arg;
    mLatSubset *m;

    sprintf(fmt, "{axis=%d,position=%d}", d, p);
    m = qlua_newLatSubset(L, fmt);
    m->heap = 1;
    arg.axis = d;
    arg.pos = p;
    lua_gc(L, LUA_GCCOLLECT, 0);
    m->s = QDP_create_subset(subset_func, &arg, sizeof (arg), 1);

    if (m->s == 0)
        luaL_error(L, "subset creation failure");

    luaL_getmetatable(L, mtnLatSubset);
    lua_setmetatable(L, -2);

    return 1;
}

static void
add_static(lua_State *L, QDP_Subset *s, const char *name)
{
    mLatSubset *m = qlua_newLatSubset(L, name);

    m->s = s;
    lua_setfield(L, -2, name);
}

static struct luaL_Reg mtLatSubset[] = {
    { "__tostring",     q_U_fmt   },
    { "__gc",           q_U_gc    },
    { "where",          q_U_where },
    { NULL,             NULL      }
};

static struct luaL_Reg fLatSubset[] = {
    { "Subset",    q_subset },
    { NULL,        NULL     }
};

int
init_latsubset(lua_State *L)
{    
    luaL_newmetatable(L, LatSubsetStack);
    qlua_metatable(L, mtnLatSubset, mtLatSubset);
    luaL_getmetatable(L, opLattice);
    luaL_register(L, NULL, fLatSubset);
    add_static(L, &QDP_all, "all");
    add_static(L, &QDP_even, "even");
    add_static(L, &QDP_odd, "odd");
    lua_pop(L, 1);
    
    return 0;
}

int
fini_latsubset(lua_State *L)
{
    return 0;
}
