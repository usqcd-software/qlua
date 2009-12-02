#include <qlua.h>                                                   /* DEPS */
#include <lattice.h>                                                /* DEPS */
#include <aff_io.h>                                                 /* DEPS */
#include <latdirprop.h>                                             /* DEPS */
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
        luaL_checktype(L, -1, LUA_TTABLE);
        for (j = 0; j < qRank; j++, k++) {
            lua_pushnumber(L, j + 1);
            lua_gettable(L, -2);
            qext[k] = luaL_checkint(L, -1);
            lua_pop(L, 1);
        }
        lua_pop(L, 1);
    }

    qlua_Aff_enter(L);
    CALL_QDP(L);
    
    status = save_bb(aff_w, key_path, F->ptr, B->ptr, csrc, tsnk, n_qext, qext,
                     time_rev, t_axis, bc_baryon);
    qlua_Aff_leave();
    
    qlua_free(L, csrc);
    qlua_free(L, qext);

    if (status)
        luaL_error(L, status);

    return 0;
}

static struct luaL_Reg fExtra[] = {
    { "save_bb",    q_save_bb },
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
