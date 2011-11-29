
#include "laph_common.h"

#include "modules.h"                                                /* DEPS */
#include "latdirferm.h"                                             /* DEPS */
#include <assert.h>

/**/
#define DEBUG_LEVEL 1
#include "debug.h"

/* calculate subgrid dimensions and "lowest corner" coordinate for `lat'
   TODO check that the subgrid is rectilinear 
 */
int 
calc_subgrid(QDP_Lattice *lat, int dims[], int cmin[])
{
    assert(LDIM == QDP_ndim_L(lat));
    int n_site = QDP_sites_on_node_L(lat);
    assert(0 < n_site);

    QDP_get_coords_L(lat, cmin, QDP_this_node, 0);
    int cmax[LDIM];
    for (int k = 0 ; k < LDIM ; k++)
        cmax[k] = cmin[k];
    
    int coord[LDIM];
    for (int i_site = 1; i_site < n_site ; i_site++) {
        QDP_get_coords_L(lat, coord, QDP_this_node, i_site);
        for (int k = 0 ; k < LDIM ; k++) {
            if (coord[k] < cmin[k])
                cmin[k] = coord[k];
            if (cmax[k] < coord[k])
                cmax[k] = coord[k];
        }
    }
    
    int vol = 1;
    for (int k = 0; k < LDIM; k++) {
        dims[k] = 1 + cmax[k] - cmin[k];
        assert(0 < dims[k]);
        vol *= dims[k];
    }
    assert(vol == n_site);

    return 0;
}

/* parse a list of {{tsrc, jvec, jspin, sol}} 
    tab_idx     table index in the stack
    n_sol, tsrc, jvec, jspin, sol
                location to put the size and pointers to arrays
    result:
        0   OK  ; tsrc, jvec, jspin, sol arrays are allocated [n_sol]
        1   Fail; no memory allocated on exit
 */
int 
qlua_check_laph_sol_list(lua_State *L, int tab_idx, mLattice *S,
               int *n_sol, int **tsrc, int **jvec, int **jspin,
               QDP_D3_DiracFermion ***sol, mLattice **have_S)
{
    assert(NULL != n_sol && NULL != tsrc && NULL != jvec 
            && NULL != jspin && NULL != sol);
    if (LUA_TTABLE != lua_type(L, tab_idx)) {
        luaL_error(L, "list of {tsrc, jvec, jspin, sol} objects expected");
        return 1;
    }

    *n_sol = lua_objlen(L, tab_idx);
    *tsrc  = qlua_malloc(L, sizeof((*tsrc)[0]) * (*n_sol));
    *jvec  = qlua_malloc(L, sizeof((*jvec)[0]) * (*n_sol));
    *jspin = qlua_malloc(L, sizeof((*jspin)[0]) * (*n_sol));
    *sol   = qlua_malloc(L, sizeof((*sol)[0]) * (*n_sol));
    assert(NULL != *tsrc && NULL != *jvec && NULL != *jspin && NULL != *sol);

    for (int isol = 0 ; isol < (*n_sol) ; isol++) {
        lua_pushnumber(L, isol + 1);
        lua_gettable(L, tab_idx);
        int sol_idx = lua_gettop(L);

        if (LUA_TTABLE != lua_type(L, sol_idx) || 4 != lua_objlen(L, sol_idx)){
            lua_pop(L, 1);
            luaL_error(L, "list of {tsrc, jvec, jspin, sol} objects expected");
            goto clearerr_1;
        }
        
        lua_pushnumber(L, 1);
        lua_gettable(L, sol_idx);
        (*tsrc)[isol]  = luaL_checkint(L, -1);
        lua_pop(L, 1);

        lua_pushnumber(L, 2); 
        lua_gettable(L, sol_idx);
        (*jvec)[isol]  = luaL_checkint(L, -1);
        lua_pop(L, 1);

        lua_pushnumber(L, 3); 
        lua_gettable(L, sol_idx);
        (*jspin)[isol]  = luaL_checkint(L, -1);
        lua_pop(L, 1);

        lua_pushnumber(L, 4); 
        lua_gettable(L, sol_idx);
        (*sol)[isol]  = qlua_checkLatDirFerm3(L, -1, S, NCOLOR)->ptr;
        
        mLattice *S1 = qlua_ObjLattice(L, -1);
        lua_pop(L, 1);
        assert(NULL != S1);
        if (NULL == S && 0 == isol) 
            S = S1;
        else {
            if (S1->id != S->id) {
                luaL_error(L, "expect array of same lattice fields");
                goto clearerr_1;
            }
        }
        lua_pop(L, 1);

        lua_pop(L, 1);
    }

    if (NULL != have_S)
        *have_S = S;

    return 0;

clearerr_1:
    qlua_free(L, *tsrc);
    qlua_free(L, *jvec);
    qlua_free(L, *jspin);
    qlua_free(L, *sol);

    return 1;
}

int
is_masternode()
{
    return (0 == QDP_this_node ||       /* masternode */
            QDP_this_node < 0);         /* qdp not initialized */
}
