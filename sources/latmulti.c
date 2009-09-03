#include <qlua.h>                                                    /* DEPS */
#include <lattice.h>                                                 /* DEPS */
#include <latmulti.h>                                                /* DEPS */

/* only axis subsets are implemented */

const char mtnLatMulti[] = "lattice.multi";

mLatMulti *
qlua_checkLatMulti(lua_State *L, int idx)
{
    void *v = luaL_checkudata(L, idx, mtnLatMulti);
    
    luaL_argcheck(L, v != 0, idx, "lattice.MultiSet expected");
    
    return v;
}

static int
q_multi_fmt(lua_State *L)
{
    mLatMulti *b = qlua_checkLatMulti(L, 1);
    char fmt[72];

    sprintf(fmt, "QDP:MultiSet(d=%d,n=%d,%p)", b->axis, b->count, b->subset);
    lua_pushstring(L, fmt);
    
    return 1;
}

static int
q_multi_gc(lua_State *L)
{
    mLatMulti *b = qlua_checkLatMulti(L, 1);

    if (b->arg) {
        qlua_free(L, b->arg);
        QDP_destroy_subset(b->subset);
    }
    b->arg = 0;
    b->subset = NULL;

    return 0;
}

static int
multi_func(int *coord, void *arg)
{
    int *axis = arg;

    return coord[*axis];
}

static int
q_latmulti(lua_State *L)
{
    int axis = qlua_checkindex(L, 2, "d", qRank);
    int *arg = qlua_malloc(L, sizeof (int)); /* just in case LUA moves udata */
    int count = qDim[axis];
    mLatMulti *v = lua_newuserdata(L, sizeof (mLatMulti));

    *arg = axis;
    v->axis = axis;
    v->arg = arg;
    v->count = count;
    lua_gc(L, LUA_GCCOLLECT, 0);
    v->subset = QDP_create_subset(multi_func, arg, sizeof (int), count);

    if (v->subset == 0)
        return luaL_error(L, "multiset creation failure");

    luaL_getmetatable(L, mtnLatMulti);
    lua_setmetatable(L, -2);

    return 1;
}

static struct luaL_Reg mtLatMulti[] = {
    { "__tostring",     q_multi_fmt },
    { "__gc",           q_multi_gc  },
    { NULL,             NULL        }
};

static struct luaL_Reg fLatMulti[] = {
    { "MultiSet",    q_latmulti },
    { NULL,          NULL       }
};

int
init_latmulti(lua_State *L)
{
    luaL_getmetatable(L, opLattice);
    luaL_register(L, NULL, fLatMulti);
    lua_pop(L, 1);
    qlua_metatable(L, mtnLatMulti, mtLatMulti);

    return 0;
}

int
fini_latmulti(lua_State *L)
{
    return 0;
}
