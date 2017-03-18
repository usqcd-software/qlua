#include <string.h>
#include <assert.h>

#include "modules.h"                                                /* DEPS */
#include "qparam.h"
#include "latcolvec.h"                                              /* DEPS */

/* parse an 1d list of integers 
    result:
        !=NULL  OK  ; return allocated array [dim1]
        ==NULL  Fail; return NULL; no memory allocated
 */
int *
qlua_check_array1d_int(lua_State *L, int idx, 
        int need_dim1, 
        int *have_dim1)
{
    if (LUA_TTABLE != lua_type(L, idx))
        luaL_error(L, "array1d(int) expected");

    int d = lua_objlen(L, idx);
    if (0 <= need_dim1 && d != need_dim1) {
        luaL_error(L, "array1d(int)[%d] expected", need_dim1);
        return NULL;
    }
    if (NULL != have_dim1)
        *have_dim1 = d;

    int *res = qlua_malloc(L, d * sizeof(int));
    for (int i = 0; i < d; i++) {
        lua_pushnumber(L, i + 1);
        lua_gettable(L, idx);
        if (lua_type(L, -1) != LUA_TNUMBER) {
            qlua_free(L, res);
            luaL_error(L, "non-int encountered in array1d(int)");
            return NULL;
        }
        res[i] = qlua_checkint(L, -1, "array element %d", i + 1);
        lua_pop(L, 1);
    }
    return res;
}


/* parse a 2d list of integers 
    result:
        !=NULL  OK  ; return allocated array [dim1][dim2]
        ==NULL  Fail; return NULL; no memory allocated
 */
int *
qlua_check_array2d_int(lua_State *L, int idx, 
        int need_dim1, int need_dim2,
        int *have_dim1, int *have_dim2)
{
    if (LUA_TTABLE != lua_type(L, idx)) {
        luaL_error(L, "array2d(int) expected");
        return NULL;
    }

    int d1 = lua_objlen(L, idx);
    if (0 <= need_dim1 && d1 != need_dim1) {
        luaL_error(L, "array2d(int)[%d,..] expected", need_dim1);
        return NULL;
    }
    if (NULL != have_dim1)
        *have_dim1 = d1;
    
    lua_pushnumber(L, 1);
    lua_gettable(L, idx);
    int d2;
    int *res0;
    if (NULL == (res0 = qlua_check_array1d_int(L, lua_gettop(L), need_dim2, &d2))) {
        luaL_error(L, "array2d(int) expected");
        return NULL;
    }

    lua_pop(L, 1);
    int *res    = qlua_malloc(L, d1 * d2 * sizeof(int));
    assert(NULL != res);
    memcpy(res, res0, d2 * sizeof(int));
    qlua_free(L, res0);

    if (NULL != have_dim2)
        *have_dim2 = d2;

    for (int i = 1; i < d1; i++) {
        lua_pushnumber(L, i + 1);
        lua_gettable(L, idx);
        if (NULL == (res0 = qlua_check_array1d_int(L, lua_gettop(L), d2, NULL))) {
            qlua_free(L, res);
            luaL_error(L, "array2d(int) expected");
            return NULL;
        }
        lua_pop(L, 1);
        memcpy(res + i * d2, res0, d2 * sizeof(int));
        qlua_free(L, res0);
    }
//    lua_pop(L, 1); // <---- why did I put it here???

    return res;
}

/* parse an 1d list of LatColVec 
    result:
        !=NULL  OK  ; return allocated array [dim1]
        ==NULL  Fail; no memory allocated
 */
QDP_D3_ColorVector **
qlua_check_latcolvec_table(lua_State *L, int idx, 
        mLattice *S, int need_dim1, 
        mLattice **have_S, int *have_dim1)
{
    if (LUA_TTABLE != lua_type(L, idx)) {
        luaL_error(L, "array1d(LatColVec) expected");
        return NULL;
    }
    int d = lua_objlen(L, idx);
    if (0 <= need_dim1 && d != need_dim1) {
        luaL_error(L, "array1d(LatColVec)[%d] expected", need_dim1);
        return NULL;
    }
    if (NULL != have_dim1)
        *have_dim1 = d;

    QDP_D3_ColorVector **res = qlua_malloc(L, d * sizeof(QDP_D3_ColorVector *));
    for (int i = 0; i < d; i++) {
        lua_pushnumber(L, i + 1);
        lua_gettable(L, idx);
        mLatColVec3 *ch = qlua_checkLatColVec3(L, lua_gettop(L), S, 3);
        assert(NULL != ch);
        res[i] = ch->ptr;

        mLattice *S1 = qlua_ObjLattice(L, lua_gettop(L));
        lua_pop(L, 1);
        if (NULL == S && 0 == i)
            S = S1;
        else {
            if (S1->id != S->id) {
                qlua_free(L, res);
                luaL_error(L, "expect array of same lattice fields");
                return NULL;
            }
        }

        lua_pop(L, 1);
    }
    if (NULL != have_S)
        *have_S = S;

    return res;
}

int 
qlua_table_objcount(lua_State *L, int idx) 
{
    if (LUA_TTABLE != lua_type(L, idx)) 
        return 0;

    int cnt = 0;
    lua_pushnil(L);
    while (0 != lua_next(L, idx)) {
        cnt ++;
        lua_pop(L, 1);
    }
    return cnt;
}


int
qlua_check_latcomplex_strmap(lua_State *L, int idx, mLattice **S, 
        int *len, const char ***strkey, QDP_D_Complex ***latcomplex)
/* 
    *S      check that all lattice objects have this lattice type; 
            if ==NULL, initialize from input
    *len    check that the number of elements is <= *len
            if ==-1, initialize from data; also, *strkey and *latcomplex must 
            be NULL to avoid memleaks
    *strkey array for strkeys; malloc if ==NULL or *len==-1
    *latcomplex array for fields; malloc if ==NULL or *len==-1
*/            
{
    assert(NULL != len);
    assert(NULL != strkey);
    assert(NULL != latcomplex);
    assert(NULL != S);
    if (LUA_TTABLE != lua_type(L, idx)) {
        luaL_error(L, "strmap(LatComplex) expected");
        return 1;
    }

    int d = qlua_table_objcount(L, idx);
    if (*len < 0) {
        assert(NULL == *strkey);
        *strkey = qlua_malloc(L, sizeof(char *) * d);
        assert(NULL == *latcomplex);
        *latcomplex = qlua_malloc(L, sizeof(QDP_D3_ColorVector *) * d);
    } else {
        assert(d <= *len);
    }
    *len = d;
    
    int cnt = 0;
    lua_pushnil(L);
    while (0 != lua_next(L, idx)) {
        assert(cnt < d);
        (*strkey)[cnt] = luaL_checkstring(L, -2);
        (*latcomplex)[cnt]  = qlua_checkLatComplex(L, -1, *S)->ptr;
//        /**/printf("[%s(%s)] : %s\n", (*strkey)[cnt], 
//                        lua_typename(L, lua_type(L,-2)), 
//                        lua_typename(L, lua_type(L,-1)));
        if (NULL == *S) {
            mLattice *S1 = qlua_ObjLattice(L, lua_gettop(L));
            lua_pop(L, 1);
            *S = S1;
        }
        lua_pop(L, 1);
        cnt++;
    }

    return 0;   /* OK */
}

