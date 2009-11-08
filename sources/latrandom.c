#include <qlua.h>                                                    /* DEPS */
#include <lattice.h>                                                 /* DEPS */
#include <latint.h>                                                  /* DEPS */
#include <latrandom.h>                                               /* DEPS */
#include <latreal.h>                                                 /* DEPS */
#include <latcomplex.h>                                              /* DEPS */
#include <latcolvec.h>                                               /* DEPS */
#include <latcolmat.h>                                               /* DEPS */
#include <latdirferm.h>                                              /* DEPS */
#include <latdirprop.h>                                              /* DEPS */
/* ZZZ other packages */

const char mtnLatRandom[] = "lattice.RandomState";

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

    luaL_argcheck(L, v != 0, idx, "lattice.RandomState expected");
    
    return v;
}

static int
q_S_fmt(lua_State *L)
{
    char fmt[72];
    mLatRandom *b = qlua_checkLatRandom(L, 1);

    sprintf(fmt, "QDP:RandomState(%p)", b->ptr);
    lua_pushstring(L, fmt);
   
    return 1;
}

static int
q_S_gc(lua_State *L)
{
    mLatRandom *b = qlua_checkLatRandom(L, 1);

    QDP_destroy_S(b->ptr);
    b->ptr = 0;
    
    return 0;
}

static int
q_S_set(lua_State *L)
{
    mLatRandom *r = qlua_checkLatRandom(L, 1);
    mLatRandom *a = qlua_checkLatRandom(L, 2);

    CALL_QDP(L);
    QDP_S_eq_S(r->ptr, a->ptr, *qCurrent);
    lua_pop(L, 1);

    return 1;
}

static int
q_latrandom(lua_State *L)
{
    int b = lua_gettop(L);

    switch (b) {
    case 2: {
        mLatRandom *a = qlua_checkLatRandom(L, 2);
        mLatRandom *b = qlua_newLatRandom(L);

        CALL_QDP(L);
        QDP_S_eq_S(b->ptr, a->ptr, *qCurrent);

        return 1;
    }
    case 3: {
        QLA_Int seed_i = luaL_checkint(L, 2);
        mLatInt *seed_I = qlua_checkLatInt(L, 3);
        mLatRandom *state = qlua_newLatRandom(L);

        CALL_QDP(L);
        QDP_S_eq_seed_i_I(state->ptr, seed_i, seed_I->ptr, *qCurrent);

        return 1;
    }
    }
    return qlua_badconstr(L, "RandomState");
}

static struct luaL_Reg mtLatRandom[] = {
    { "__tostring",               q_S_fmt },
    { "__gc",                     q_S_gc },
    { "random_Real",              q_R_random },
    { "gaussian_Real",            q_R_gaussian },
    { "gaussian_Complex",         q_C_gaussian },
    { "gaussian_ColorVector",     q_V_gaussian },
    { "gaussian_ColorMatrix",     q_M_gaussian },
    { "gaussian_DiracFermion",    q_D_gaussian },
    { "gaussian_DiracPropagator", q_P_gaussian },
    { "set",                      q_S_set },
    /* ZZZ other gaussian randoms */
    { NULL,          NULL}
};

static struct luaL_Reg fLatRandom[] = {
    { "RandomState", q_latrandom},
    { NULL,          NULL }
};

int
init_latrandom(lua_State *L)
{
    luaL_getmetatable(L, opLattice);
    luaL_register(L, NULL, fLatRandom);
    lua_pop(L, 1);
    qlua_metatable(L, mtnLatRandom, mtLatRandom);

    return 0;
}

int
fini_latrandom(lua_State *L)
{
    return 0;
}
