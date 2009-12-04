#include <qlua.h>                                                   /* DEPS */
#include <lattice.h>                                                /* DEPS */
#include <aff_io.h>                                                 /* DEPS */
#include <latdirferm.h>                                             /* DEPS */
#include <latdirprop.h>                                             /* DEPS */
#include <latcolmat.h>                                              /* DEPS */
#include <latcolvec.h>                                              /* DEPS */
#include <extras.h>                                                 /* DEPS */

static int
q_save_bb(lua_State *L)
{
    mAffWriter *aff_w = qlua_checkAffWriter(L, 1);
    const char *key_path = luaL_checkstring(L, 2);
    mLatDirProp *F = qlua_checkLatDirProp(L, 3);
    mLatDirProp *B = qlua_checkLatDirProp(L, 4);
    int *csrc = qlua_checkintarray(L, 5, qRank, NULL);
    int tsnk = luaL_checkint(L, 6);
    int n_qext; /* extracted from #7 */
    int *qext; /* #7 */
    int time_rev = luaL_checkint(L, 8);
    int t_axis = luaL_checkint(L, 9);
    double bc_baryon = luaL_checknumber(L, 10);
    int i, j, k;
    const char *status = NULL;

    if (csrc == NULL)
        return luaL_error(L, "bad value for coord_src");

    luaL_checktype(L, 7, LUA_TTABLE);
    n_qext = lua_objlen(L, 7);
    qext = qlua_malloc(L, n_qext * qRank * sizeof (int));
    for (k = i = 0; i < n_qext; i++) {
        lua_pushnumber(L, i + 1);
        lua_gettable(L, 7);
        qlua_checktable(L, -1, "momentum at #7[%d]", i + 1);
        for (j = 0; j < qRank; j++, k++) {
            lua_pushnumber(L, j + 1);
            lua_gettable(L, -2);
            qext[k] = qlua_checkint(L, -1, "momentum component at #7[%d][%d]",
                                    i + 1, j + 1);
            lua_pop(L, 1);
        }
        lua_pop(L, 1);
    }

    qlua_Aff_enter(L);
    CALL_QDP(L);
    
    status = save_bb(L, aff_w, key_path, F->ptr, B->ptr, csrc, tsnk, n_qext, qext,
                     time_rev, t_axis, bc_baryon);
    qlua_Aff_leave();
    
    qlua_free(L, csrc);
    qlua_free(L, qext);

    if (status)
        luaL_error(L, status);

    return 0;
}

static int
q_laplacian(lua_State *L)
{
    double a = luaL_checknumber(L, 1); /* [1] : a */
    double b = luaL_checknumber(L, 2); /* [2] : b */
    QDP_ColorMatrix **g; /* [3]: { U[0], ... } */
    /* [4]: input */
    int skip; /* [5]: skip axis (if present) */
    int i;
    const char *status = NULL;

    switch (lua_type(L, 5)) {
    case LUA_TNUMBER:
        skip = luaL_checkint(L, 5);
        break;
    case LUA_TNONE:
        skip = -1;
        break;
    default:
        return luaL_error(L, "axis must be a number");
    }

    g = qlua_malloc(L, qRank * sizeof (QDP_ColorMatrix *));
    qlua_checktable(L, 3, "gauge field expected");
    for (i = 0; i < qRank; i++) {
        lua_pushnumber(L, i + 1);
        lua_gettable(L, 3);
        g[i] = qlua_checkLatColMat(L, -1)->ptr;
        lua_pop(L, 1);
    }

    switch (qlua_gettype(L, 4)) {
    case qLatColVec: {
        mLatColVec *x = qlua_checkLatColVec(L, 4);
        mLatColVec *r = qlua_newLatColVec(L);
        CALL_QDP(L);
        status = gen_laplace_V(L, r->ptr, a, b, g, x->ptr, skip);
        break;
    }
    case qLatColMat: {
        mLatColMat *x = qlua_checkLatColMax(L, 4);
        mLatColMat *r = qlua_newLatColMax(L);
        CALL_QDP(L);
        status = gen_laplace_M(L, r->ptr, a, b, g, x->ptr, skip);
        break;
    }
    case qLatDirFerm: {
        mLatDirFerm *x = qlua_checkLatDirFerm(L, 4);
        mLatDirFerm *r = qlua_newLatDirFerm(L);
        CALL_QDP(L);
        status = gen_laplace_D(L, r->ptr, a, b, g, x->ptr, skip);
        break;
    }
    case qLatDirProp: {
        mLatDirProp *x = qlua_checkLatDirProp(L, 4);
        mLatDirProp *r = qlua_newLatDirProp(L);
        CALL_QDP(L);
        status = gen_laplace_P(L, r->ptr, a, b, g, x->ptr, skip);
        break;
    }
    default:
        return luaL_error(L, "arg #4 must be colored");
    }

    qlua_free(L, g);
    if (status)
        return luaL_error(L, status);
    return 1;
}

static struct luaL_Reg fExtra[] = {
    { "save_bb",    q_save_bb },
    { "laplacian",  q_laplacian },
    { NULL,         NULL }
};

int
init_extras(lua_State *L)
{
    luaL_register(L, qcdlib, fExtra);
    return 0;
}

int
fini_extras(lua_State *L)
{
    return 0;
}
