#include "modules.h"                                                 /* DEPS */
#include "qlua.h"                                                    /* DEPS */
#include "qcomplex.h"                                                /* DEPS */
#include "seqrandom.h"                                               /* DEPS */
#include "seqcolvec.h"                                               /* DEPS */
#include "seqcolmat.h"                                               /* DEPS */
#include "seqdirferm.h"                                              /* DEPS */
#include "seqdirprop.h"                                              /* DEPS */

#include <string.h>
/* ZZZ other packages */
static const char mtnSeqRandom[] = "qlua.mtSeqRandomState";
static int
q_s_fmt(lua_State *L)
{
    char fmt[72];
    mSeqRandom *b = qlua_checkSeqRandom(L, 1);

    sprintf(fmt, "qcd.RandomState(%p)", b);
    lua_pushstring(L, fmt);
   
    return 1;
}

static int
q_s_set(lua_State *L)
{
    mSeqRandom *r = qlua_checkSeqRandom(L, 1);
    mSeqRandom *a = qlua_checkSeqRandom(L, 2);

    CALL_QDP(L);
    QLA_S_eq_S(&r->state, &a->state);
    lua_pop(L, 1);

    return 1;
}

static int
q_seqrandom(lua_State *L)
{
    int b = lua_gettop(L);

    switch (b) {
    case 1: {
        mSeqRandom *a = qlua_checkSeqRandom(L, 1);
        mSeqRandom *b = qlua_newSeqRandom(L);

        CALL_QDP(L);
        QLA_S_eq_S(&b->state, &a->state);

        return 1;
    }
    case 2: {
        QLA_Int seed_i = luaL_checkint(L, 1);
        QLA_Int seed_I = luaL_checkint(L, 2);
        mSeqRandom *state = qlua_newSeqRandom(L);

        QLA_S_eq_seed_i_I(&state->state, seed_i, &seed_I);

        return 1;
    }
    }
    return qlua_badconstr(L, "RandomState");
}

static int
q_r_random(lua_State *L)
{
    mSeqRandom *a = qlua_checkSeqRandom(L, 1);
    QLA_D_Real value;

    QLA_D_R_eq_random_S(&value, &a->state);
    lua_pushnumber(L, value);
    return 1;
}

static int
q_r_gaussian(lua_State *L)
{
    mSeqRandom *a = qlua_checkSeqRandom(L, 1);
    QLA_D_Real value;

    QLA_D_R_eq_gaussian_S(&value, &a->state);
    lua_pushnumber(L, value);
    return 1;
}

static struct luaL_Reg mtSeqRandom[] = {
    { "__tostring",               q_s_fmt        },
    { "__newindex",               qlua_nowrite   },
    { "random_Real",              q_r_random     },
    { "gaussian_Real",            q_r_gaussian   },
    { "gaussian_Complex",         q_c_gaussian   },
    { "gaussian_ColorVectorN",    q_v_gaussian_N },
    { "gaussian_ColorMatrixN",    q_m_gaussian_N },
    { "gaussian_DiracFermionN",   q_d_gaussian_N },
    { "gaussian_DiracPropagatorN",q_p_gaussian_N },
    { "set",                      q_s_set        },
    /* ZZZ other gaussian randoms */
    { NULL,                       NULL           }
};

mSeqRandom *
qlua_newSeqRandom(lua_State *L)
{
    mSeqRandom *m;

    m = lua_newuserdata(L, sizeof (mSeqRandom));
    memset(m, 0, sizeof (mSeqRandom));
    luaL_getmetatable(L, mtnSeqRandom);
    lua_setmetatable(L, -2);

    return m;
}

mSeqRandom *
qlua_checkSeqRandom(lua_State *L, int idx)
{
    void *ud = luaL_checkudata(L, idx, mtnSeqRandom);

    luaL_argcheck(L, ud != NULL, idx, "RandomState expected");

    return (mSeqRandom *)ud;
}

static struct luaL_Reg fSeqRandom[] = {
    { "RandomState", q_seqrandom},
    { NULL,          NULL }
};

int
init_seqrandom(lua_State *L)
{
    luaL_register(L, qcdlib, fSeqRandom);
    qlua_metatable(L, mtnSeqRandom, mtSeqRandom, qSeqRandom);

    return 0;
}

int
fini_seqrandom(lua_State *L)
{
    return 0;
}
