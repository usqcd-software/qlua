#ifndef QPARAM_H_ipqo8ybrpPCnj8kyzdob
#define QPARAM_H_ipqo8ybrpPCnj8kyzdob

#include "qlua.h"                                                   /* DEPS */
#include "lattice.h"                                                /* DEPS */
#include "latcomplex.h"
#include "latdirprop.h"

int *qlua_check_array1d_int(
        lua_State *L, int idx,
        int need_dim1,
        int *have_dim1);

int *qlua_check_array2d_int(
        lua_State *L, int idx,
        int need_dim1, int need_dim2,
        int *have_dim1, int *have_dim2);

QDP_D3_ColorVector **
qlua_check_latcolvec_table(lua_State *L, int idx,
        mLattice *S, int need_dim1,
        mLattice **have_S, int *have_dim1);

int qlua_table_objcount(lua_State *L, int idx);

int qlua_check_latcomplex_strmap(
        lua_State *L, int idx, mLattice **S, 
        int *len, const char ***strkey, QDP_D_Complex ***latcomplex);

#endif/*QPARAM_H_ipqo8ybrpPCnj8kyzdob*/
