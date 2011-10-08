#include "modules.h"                                                /* DEPS */
#include "qlua.h"                                                   /* DEPS */
#ifdef HAS_GSL
#include "qmatrix.h"                                                /* DEPS */
#endif /* defined(HAS_GSL) */
#include "qgamma.h"                                                 /* DEPS */
#include "lattice.h"                                                /* DEPS */
#include "aff_io.h"                                                 /* DEPS */
#include "latdirferm.h"                                             /* DEPS */
#include "latdirprop.h"                                             /* DEPS */
#include "latcolmat.h"                                              /* DEPS */
#include "latcolvec.h"                                              /* DEPS */
#include "latcomplex.h"                                             /* DEPS */
#include "extras.h"                                                 /* DEPS */

#ifdef HAS_GSL
static int
q_baryon_duu(lua_State *L)
{
    mLatDirProp3 *Pd = qlua_checkLatDirProp3(L, 1, NULL, 3);
    mLattice *S = qlua_ObjLattice(L, 1);
	int Sidx = lua_gettop(L);
    mLatDirProp3 *Pu11 = qlua_checkLatDirProp3(L, 2, S, 3);
    mLatDirProp3 *Pu12 = qlua_checkLatDirProp3(L, 3, S, 3);
    mLatDirProp3 *Pu21 = qlua_checkLatDirProp3(L, 4, S, 3);
    mLatDirProp3 *Pu22 = qlua_checkLatDirProp3(L, 5, S, 3);
	mMatComplex *Sf = gamma2matrix(L, 6);
	mMatComplex *Si_bar = gamma2matrix(L, 7);
	mMatComplex *RTR = gamma2matrix(L, 8);
	mLatComplex *B = qlua_newLatComplex(L, Sidx);

	if ((Sf->l_size != QDP_Ns) || (Sf->r_size != QDP_Ns))
		return luaL_error(L, "bad size of Sf (%d, %d)", Sf->l_size, Sf->r_size);

	if ((Si_bar->l_size != QDP_Ns) || (Si_bar->r_size != QDP_Ns))
		return luaL_error(L, "bad size of Si_bar (%d, %d)", Si_bar->l_size, Si_bar->r_size);

	if ((RTR->l_size != QDP_Ns) || (RTR->r_size != QDP_Ns))
		return luaL_error(L, "bad size of RTR (%d, %d)", RTR->l_size, RTR->r_size);

	CALL_QDP(L);
	baryon_duu(S, B->ptr,
			   Pd->ptr, Pu11->ptr, Pu12->ptr, Pu21->ptr, Pu22->ptr,
			   Sf->m, Si_bar->m, RTR->m);

	return 1;
}
#endif /* defined(HAS_GSL) */

static int
q_save_bb(lua_State *L)
{
    mAffWriter *aff_w = qlua_checkAffWriter(L, 1);
    const char *key_path = luaL_checkstring(L, 2);
    mLatDirProp3 *F = qlua_checkLatDirProp3(L, 3, NULL, 3);
    mLattice *S = qlua_ObjLattice(L, 3);
    int Sidx = lua_gettop(L);
    mLatDirProp3 *B = qlua_checkLatDirProp3(L, 4, S, 3);
    int *csrc = qlua_checkintarray(L, 5, S->rank, NULL);
    int tsnk = luaL_checkint(L, 6);
    int time_rev = luaL_checkint(L, 8);
    int t_axis = luaL_checkint(L, 9);
    double bc_baryon = luaL_checknumber(L, 10);
    int i, j, k;
    const char *status = NULL;

    if (Sidx < 11)
        return luaL_error(L, "bad arguments");

    if (csrc == NULL)
        return luaL_error(L, "bad value for coord_src");

    luaL_checktype(L, 7, LUA_TTABLE);
    int n_qext = lua_objlen(L, 7);
    int qext[n_qext * S->rank];
    for (k = i = 0; i < n_qext; i++) {
        lua_pushnumber(L, i + 1);
        lua_gettable(L, 7);
        qlua_checktable(L, -1, "momentum at #7[%d]", i + 1);
        for (j = 0; j < S->rank; j++, k++) {
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
    
    status = save_bb(L, S, aff_w, key_path, F->ptr, B->ptr, csrc, tsnk, n_qext, qext,
                     time_rev, t_axis, bc_baryon);
    qlua_Aff_leave();
    
    qlua_free(L, csrc);

    if (status)
        luaL_error(L, status);

    return 0;
}

static int
q_laplacian(lua_State *L)
{
    double a = luaL_checknumber(L, 1); /* [1] : a */
    double b = luaL_checknumber(L, 2); /* [2] : b */
    /* [3]: { U[0], ... } */
    /* [4]: input */
    int skip; /* [5]: skip axis (if present) */
    int i;
    const char *status = NULL;

    switch (lua_type(L, 5)) {
    case LUA_TNUMBER:
        skip = luaL_checkint(L, 5);
        break;
    case LUA_TNONE:
    case LUA_TNIL:
        skip = -1;
        break;
    default:
        return luaL_error(L, "axis must be a number");
    }

    qlua_checktable(L, 3, "gauge field expected");
    lua_pushnumber(L, 1);
    lua_gettable(L, 3);
    mLattice *S = qlua_ObjLattice(L, -1);
    int Sidx = lua_gettop(L);
    if (Sidx < 6)
        return luaL_error(L, "bad arguments");
    
    QDP_D3_ColorMatrix *g[S->rank];
    for (i = 0; i < S->rank; i++) {
        lua_pushnumber(L, i + 1);
        lua_gettable(L, 3);
        g[i] = qlua_checkLatColMat3(L, -1, S, 3)->ptr;
        lua_pop(L, 1);
    }

    switch (qlua_qtype(L, 4)) {
    case qLatColVec3: {
        mLatColVec3 *x = qlua_checkLatColVec3(L, 4, S, 3);
        mLatColVec3 *r = qlua_newLatColVec3(L, Sidx, 3);
        CALL_QDP(L);
        status = gen_laplace_V(L, S, r->ptr, a, b, g, x->ptr, skip);
        break;
    }
    case qLatColMat3: {
        mLatColMat3 *x = qlua_checkLatColMat3(L, 4, S, 3);
        mLatColMat3 *r = qlua_newLatColMat3(L, Sidx, 3);
        CALL_QDP(L);
        status = gen_laplace_M(L, S, r->ptr, a, b, g, x->ptr, skip);
        break;
    }
    case qLatDirFerm3: {
        mLatDirFerm3 *x = qlua_checkLatDirFerm3(L, 4, S, 3);
        mLatDirFerm3 *r = qlua_newLatDirFerm3(L, Sidx, 3);
        CALL_QDP(L);
        status = gen_laplace_D(L, S, r->ptr, a, b, g, x->ptr, skip);
        break;
    }
    case qLatDirProp3: {
        mLatDirProp3 *x = qlua_checkLatDirProp3(L, 4, S, 3);
        mLatDirProp3 *r = qlua_newLatDirProp3(L, Sidx, 3);
        CALL_QDP(L);
        status = gen_laplace_P(L, S, r->ptr, a, b, g, x->ptr, skip);
        break;
    }
    default:
        return luaL_error(L, "arg #4 must be colored");
    }

    if (status)
        return luaL_error(L, status);
    return 1;
}

extern int q_save_squark_wf(lua_State *L);
extern int q_save_q2pt(lua_State *L);
extern int q_save_q2pt_list(lua_State *L);

static struct luaL_Reg fExtra[] = {
    { "save_bb",    q_save_bb },
    { "laplacian",  q_laplacian },
#ifdef HAS_GSL
	{ "baryon_duu", q_baryon_duu },
#endif /* defined(HAS_GSL) */
    { "save_squark_wf",     q_save_squark_wf },
    { "save_q2pt",          q_save_q2pt },
    { "save_q2pt_list",     q_save_q2pt_list },
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
