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


static struct luaL_Reg fExtra[] = {
    { "save_bb",                    q_save_bb },
    { "laplacian",                  q_laplacian },
/* (latcomplex x, sign, dir)
    dir     (= 0 .. (NDIM-1)) : FT direction
            nil : do all directions
    sign    +1: y_k = \sum_x exp(+ikx) x_k
            -1: y_k = \sum_x exp(-ikx) x_k

    RETURN latcomplex y on the same lattice
 */
static int
q_fourier_transf(lua_State *L)
{
    const char *status = NULL;
    mLattice *S = qlua_ObjLattice(L, 1);
    int ndim = QDP_ndim_L(S->lat);
    int Sidx = lua_gettop(L);
    int ft_sign = luaL_checkint(L, 2);
    int ft_dir = -1;    /* all directions by default */
    if (4 <= Sidx)
        ft_dir  = luaL_checkint(L, 3);
    if (ndim <= ft_dir) {
        luaL_error(L, "bad direction in arg #3");
    }
    
    CALL_QDP(L);
    /* XXX the code relies on particular organization of QLA structures: 
       for each site, an array of double complex numbers, no padding */
    switch (qlua_qtype(L, 1)) {
    case qLatComplex: {
        mLatComplex *x = qlua_checkLatComplex(L, 1, S);
        mLatComplex *r = qlua_newLatComplex(L, Sidx);
        QLA_D_Complex *qla_x = QDP_D_expose_C(x->ptr),
                      *qla_r = QDP_D_expose_C(r->ptr);
        CALL_QDP(L);
        status = extra_lat_fourier_qla(L, S, qla_r, qla_x, 1, ft_sign, ft_dir);
        QDP_D_reset_C(x->ptr);
        QDP_D_reset_C(r->ptr);
        break;
    }
    case qLatColVec3: {
        mLatColVec3 *x = qlua_checkLatColVec3(L, 1, S, 3);
        mLatColVec3 *r = qlua_newLatColVec3(L, Sidx, 3);
        QLA_D_Complex *qla_x = QDP_D_expose_V(x->ptr),
                      *qla_r = QDP_D_expose_V(r->ptr);
        CALL_QDP(L);
        status = extra_lat_fourier_qla(L, S, (QLA_D_Complex *)qla_r, qla_x, 3, ft_sign, ft_dir);
        QDP_D_reset_V(x->ptr);
        QDP_D_reset_V(r->ptr);
        break;
    }
    case qLatColMat3: {
        mLatColMat3 *x = qlua_checkLatColMat3(L, 1, S, 3);
        mLatColMat3 *r = qlua_newLatColMat3(L, Sidx, 3);
        QLA_D_Complex *qla_x = QDP_D_expose_M(x->ptr),
                      *qla_r = QDP_D_expose_M(r->ptr);
        CALL_QDP(L);
        status = extra_lat_fourier_qla(L, S, qla_r, qla_x, 9, ft_sign, ft_dir);
        QDP_D_reset_M(x->ptr);
        QDP_D_reset_M(r->ptr);
        break;
    }
    case qLatDirFerm3: {
        mLatDirFerm3 *x = qlua_checkLatDirFerm3(L, 1, S, 3);
        mLatDirFerm3 *r = qlua_newLatDirFerm3(L, Sidx, 3);
        QLA_D_Complex *qla_x = QDP_D_expose_D(x->ptr),
                      *qla_r = QDP_D_expose_D(r->ptr);
        CALL_QDP(L);
        status = extra_lat_fourier_qla(L, S, qla_r, qla_x, 12, ft_sign, ft_dir);
        QDP_D_reset_D(x->ptr);
        QDP_D_reset_D(r->ptr);
        break;
    }
    case qLatDirProp3: {
        mLatDirProp3 *x = qlua_checkLatDirProp3(L, 1, S, 3);
        mLatDirProp3 *r = qlua_newLatDirProp3(L, Sidx, 3);
        QLA_D_Complex *qla_x = QDP_D_expose_P(x->ptr),
                      *qla_r = QDP_D_expose_P(r->ptr);
        CALL_QDP(L);
        status = extra_lat_fourier_qla(L, S, qla_r, qla_x, 144, ft_sign, ft_dir);
        QDP_D_reset_P(x->ptr);
        QDP_D_reset_P(r->ptr);
        break;
    }
    default: 
        return luaL_error(L, "arg #1 must be complex");
    }

    if (NULL != status)
        return luaL_error(L, status);
    return 1;

}



static struct luaL_Reg fExtra[] = {
    { "save_bb",    q_save_bb },
    { "laplacian",  q_laplacian },
    { "fourier_transf", q_fourier_transf },
#ifdef HAS_GSL
    { "baryon_duu",                 q_baryon_duu },
#endif /* defined(HAS_GSL) */
    { "save_laph_wf_baryon_pwave",  q_save_laph_wf_baryon_pwave },
    { "save_q2pt",                  q_save_q2pt },
    { "save_q2pt_list",             q_save_q2pt_list },
    { "save_q3pt_selectspin",       q_save_q3pt_0deriv_selectspin },
    { "save_npr_prop",              q_save_npr_prop },
    { "save_npr_2qvertex",          q_save_npr_2qvertex },
    { NULL,         NULL }
};

int
init_extras(lua_State *L)
{
    luaL_register(L, qcdlib, fExtra);
    return 0;
}

void
fini_extras(void)
{
}
