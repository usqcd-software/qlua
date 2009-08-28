#include <qlua.h>
#include <latint.h>
#include <latrandom.h>
#include <latreal.h>

const char *mtnLatRandom = "qcd.lattice.random";

mLatRandom *
qlua_newLatRandom(lua_State *L)
{
    QDP_RandomState *v = QDP_create_S();
    mLatRandom *hdr;

    if (v == 0) {
        lua_gc(L, LUA_GCCOLLECT, 0);
        v = QDP_create_S();
        if (v == 0)
            luaL_error(L, "not enough memory (QDP_RandomState)");
    }
    hdr = lua_newuserdata(L, sizeof (mLatRandom));
    hdr->ptr = v;
    luaL_getmetatable(L, mtnLatRandom);
    lua_setmetatable(L, -2);

    return hdr;
}

mLatRandom *
qlua_checkLatRandom(lua_State *L, int idx)
{
    void *v = luaL_checkudata(L, idx, mtnLatRandom);

    luaL_argcheck(L, v != 0, idx, "qcd.LatRandom expected");
    
    return v;
}

static int
qLatRandom_fmt(lua_State *L)
{
    char fmt[72];
    mLatRandom *b = qlua_checkLatRandom(L, 1);

    sprintf(fmt, "LatRandom(%p)", b->ptr);
    lua_pushstring(L, fmt);
   
    return 1;
}

static int
qLatRandom_gc(lua_State *L)
{
    mLatRandom *b = qlua_checkLatRandom(L, 1);

    QDP_destroy_S(b->ptr);
    b->ptr = 0;
    
    return 0;
}

static int
q_S_eq_S(lua_State *L)
{
    mLatRandom *a = qlua_checkLatRandom(L, 1);
    mLatRandom *b = qlua_newLatRandom(L);

    QDP_S_eq_S(b->ptr, a->ptr, QDP_all);
    return 1;
}

static struct luaL_Reg mtLatRandom[] = {
    { "__tostring",       qLatRandom_fmt },
    { "__gc",             qLatRandom_gc },
    { "random_real",      q_R_random },
    { "gaussian_real",    q_R_gaussian },
    /* supported random generators for QDP types */
    { NULL,          NULL}
};

static int
q_latrandom(lua_State *L)
{
    int b = lua_gettop(L);

    switch (b) {
    case 1:
        return q_S_eq_S(L);
    case 2: {
        QLA_Int seed_i = luaL_checkinteger(L, 1);
        mLatInt *seed_I = qlua_checkLatInt(L, 2);
        mLatRandom *state = qlua_newLatRandom(L);
        QDP_S_eq_seed_i_I(state->ptr, seed_i, seed_I->ptr, QDP_all);
        break;
    }
    default:
        return luaL_error(L, "bad arguments");
    }
    return 1;
}

static struct luaL_Reg fLatRandom[] = {
    { "lat_random", q_latrandom},
    { NULL, NULL }
};

int
init_latrandom(lua_State *L)
{
    luaL_register(L, qcdlib, fLatRandom);
    qlua_metatable(L, mtnLatRandom, mtLatRandom);

    return 0;
}

int
fini_latrandom(lua_State *L)
{
    return 0;
}
