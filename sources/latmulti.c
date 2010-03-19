#include "qlua.h"                                                    /* DEPS */
#include "latsubset.h"                                               /* DEPS */
#include "lattice.h"                                                 /* DEPS */
#include "latint.h"                                                  /* DEPS */
#include "latmulti.h"                                                /* DEPS */

/* QLUA multisets are not QDP multisets! */

static const char LatMultiName[] = "lattice.multi";

static int
q_multi_fmt(lua_State *L)
{
    mLatMulti *v = qlua_checkLatMulti(L, 1, NULL);
    char fmt[72];

    sprintf(fmt, "MultiSet(%d,...)", v->size);
    lua_pushstring(L, fmt);
    
    return 1;
}

static int
q_multi_gc(lua_State *L)
{
    mLatMulti *v = qlua_checkLatMulti(L, 1, NULL);

    if (v->idx)
        qlua_free(L, v->idx);
    v->idx = 0;

    return 0;
}

static struct luaL_Reg mtLatMulti[] = {
    { "__tostring",     q_multi_fmt  },
    { "__gc",           q_multi_gc   },
    { "__newindex",     qlua_nowrite },
    /* "lattice" is inserted when the table is creates */
    /* "a-type"  is inserted as well */
    { NULL,             NULL        }
};

mLatMulti *
qlua_newLatMulti(lua_State *L, int Sidx)
{
    mLatMulti *v = lua_newuserdata(L, sizeof (mLatMulti));
    
    v->size = 0;
    v->idx = qlua_malloc(L, QDP_sites_on_node * sizeof (int));
    qlua_createLatticeTable(L, Sidx, mtLatMulti, qLatMulti, LatMultiName);
    lua_setmetatable(L, -2);

    return v;
}

mLatMulti *
qlua_checkLatMulti(lua_State *L, int idx, mLattice *S)
{
    void *v = qlua_checkLatticeType(L, idx, qLatMulti, LatMultiName);
    
    if (S) {
        mLattice *S1 = qlua_ObjLattice(L, idx);
        if (S1->id != S->id)
            luaL_error(L, "%s on a wrong lattice", LatMultiName);
    }

    return (mLatMulti *)v;
}

static int
q_latmulti(lua_State *L)
{
    int size = luaL_checkint(L, 2);
    mLatInt *m = qlua_checkLatInt(L, 3, NULL);
    mLatMulti *v = qlua_newLatMulti(L, lua_gettop(L));
    QLA_Int *mm;
    int k;

    v->size = size;
    CALL_QDP(L);
    mm = QDP_expose_I(m->ptr);
    for (k = 0; k < QDP_sites_on_node; k++)
        v->idx[k] = mm[k];
    QDP_reset_I(m->ptr);

    return 1;
}

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

    return 0;
}

int
fini_latmulti(lua_State *L)
{
    return 0;
}
