#include "qlua.h"                                                    /* DEPS */
#include "lattice.h"                                                 /* DEPS */
#include "latsubset.h"                                               /* DEPS */
#include "latint.h"                                                  /* DEPS */
#include "latrandom.h"                                               /* DEPS */
#include "latreal.h"                                                 /* DEPS */
#include "latcomplex.h"                                              /* DEPS */
#include "latcolvec.h"                                               /* DEPS */
#include "latcolmat.h"                                               /* DEPS */
#if 0 /* XXX includes */
#include "latdirferm.h"                                              /* DEPS */
#include "latdirprop.h"                                              /* DEPS */
#endif /* XXX includes */
/* ZZZ other packages */

static const char LatRandomName[] = "lattice.RandomState";

static int
q_S_fmt(lua_State *L)
{
    char fmt[72];
    mLatRandom *b = qlua_checkLatRandom(L, 1, NULL);

    sprintf(fmt, "QDP:RandomState(%p)", b->ptr);
    lua_pushstring(L, fmt);
   
    return 1;
}

static int
q_S_gc(lua_State *L)
{
    mLatRandom *b = qlua_checkLatRandom(L, 1, NULL);

    QDP_destroy_S(b->ptr);
    b->ptr = 0;
    
    return 0;
}

static int
q_S_set(lua_State *L)
{
    mLatRandom *r = qlua_checkLatRandom(L, 1, NULL);
    mLattice *S = qlua_ObjLattice(L, 1);
    mLatRandom *a = qlua_checkLatRandom(L, 2, S);

    CALL_QDP(L);
    QDP_S_eq_S(r->ptr, a->ptr, *S->qss);
    lua_pop(L, 1);

    return 1;
}

static int
q_latrandom(lua_State *L)
{
    mLattice *S = qlua_checkLattice(L, 1);
    int b = lua_gettop(L);

    switch (b) {
    case 2: {
        mLatRandom *a = qlua_checkLatRandom(L, 2, S);
        mLatRandom *b = qlua_newLatRandom(L, 1);

        CALL_QDP(L);
        QDP_S_eq_S(b->ptr, a->ptr, *S->qss);

        return 1;
    }
    case 3: {
        QLA_Int seed_i = luaL_checkint(L, 2);
        mLatInt *seed_I = qlua_checkLatInt(L, 3, S);
        mLatRandom *state = qlua_newLatRandom(L, 1);

        CALL_QDP(L);
        QDP_S_eq_seed_i_I(state->ptr, seed_i, seed_I->ptr, *S->qss);

        return 1;
    }
    }
    return qlua_badconstr(L, "RandomState");
}

static struct luaL_Reg mtLatRandom[] = {
    { "__tostring",               q_S_fmt        },
    { "__gc",                     q_S_gc         },
    { "__newindex",               qlua_nowrite   },
    { "random_Real",              q_R_random     },
    { "gaussian_Real",            q_R_gaussian   },
    { "gaussian_Complex",         q_C_gaussian   },
    { "gaussian_ColorVector",     q_V_gaussian   },
    { "gaussian_ColorVectorN",    q_V_gaussian_N },
    { "gaussian_ColorMatrix",     q_M_gaussian   },
    { "gaussian_ColorMatrixN",    q_M_gaussian_N },
#if 0 /* XXX random objects */
    { "gaussian_DiracFermion",    q_D_gaussian   },
    { "gaussian_DiracPropagator", q_P_gaussian   },
#endif /* XXX random objects */
    { "set",                      q_S_set        },
    /* ZZZ other gaussian randoms */
    { NULL,                       NULL           }
};

mLatRandom *
qlua_newLatRandom(lua_State *L, int Sidx)
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
    qlua_createLatticeTable(L, Sidx, mtLatRandom, qLatRandom, LatRandomName);
    lua_setmetatable(L, -2);

    return hdr;
}

mLatRandom *
qlua_checkLatRandom(lua_State *L, int idx, mLattice *S)
{
    void *v = qlua_checkLatticeType(L, idx, qLatRandom, LatRandomName);
    
    if (S) {
        mLattice *S1 = qlua_ObjLattice(L, idx);
        if (S1->id != S->id)
            luaL_error(L, "%s on a wrong lattice", LatRandomName);
        lua_pop(L, 1);
    }

    return (mLatRandom *)v;
}

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

    return 0;
}

int
fini_latrandom(lua_State *L)
{
    return 0;
}
