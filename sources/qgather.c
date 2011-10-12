#include "modules.h"                                                 /* DEPS */
#include "qlua.h"                                                    /* DEPS */
#include "lattice.h"                                                 /* DEPS */
#include "latsubset.h"                                               /* DEPS */
#include "latcolmat.h"                                               /* DEPS */
#include "latcolvec.h"                                               /* DEPS */
#include "latcomplex.h"                                              /* DEPS */
#include "latdirferm.h"                                              /* DEPS */
#include "latdirprop.h"                                              /* DEPS */
#include "latint.h"                                                  /* DEPS */
#include "latreal.h"                                                 /* DEPS */
#include "qgather.h"                                                 /* DEPS */
#include "qdp.h"

static const char GatherName[] = "qlua.gather";
static const char src_lattice_key[]  = "gather.src.lattice";
static const char dst_lattice_key[]  = "gather.dst.lattice";

static int
q_gather_fmt(lua_State *L)
{
    char fmt[72];
    mGather *b = qlua_checkGather(L, 1);

    sprintf(fmt, "QDP:Gather(%p)", b->ptr);
    lua_pushstring(L, fmt);

    return 1;
}

static int
q_gather_gc(lua_State *L)
{
    mGather *b = qlua_checkGather(L, 1);

    QDP_destroy_collect(b->ptr);
    b->ptr = 0;

    return 0;
}

mGather *
qlua_newGather(lua_State *L)
{
    mGather *g;

    g = lua_newuserdata(L, sizeof (mGather));
    g->ptr = 0;

    return g;
}

mGather *
qlua_checkGather(lua_State *L, int idx)
{
    void *v = qlua_checkLatticeType(L, idx, qGather, GatherName);
    return v;
}


static void
gather_and_I(QLA_Int *res, QLA_Int *val, void *env)
{
    QLA_Int r;
    QLA_I_eq_I_and_I(&r, res, val);
    QLA_I_eq_I(res, &r);
}

static int
q_gather_and(lua_State *L)
{
    int top = lua_gettop(L);
    mGather *g = qlua_checkGather(L, 1);
    mLattice *dstLattice;
    mLattice *srcLattice;
    
    luaL_getmetafield(L, 1, dst_lattice_key);
    dstLattice = qlua_checkLattice(L, top + 1);
    luaL_getmetafield(L, 1, src_lattice_key);
    srcLattice = qlua_checkLattice(L, top + 2);

    switch (qlua_qtype(L, 2)) {
    case qLatInt: {
        QLA_Int zval = -1; /* assume two-complement encoding if QLA_Int */
        mLatInt *src = qlua_checkLatInt(L, 2, srcLattice);
        mLatInt *dst = qlua_newZeroLatInt(L, top + 1);
        CALL_QDP(L);
        QDP_I_eq_collect_I_funca(dst->ptr, g->ptr, src->ptr,
                                 &zval, gather_and_I, NULL,
                                 *srcLattice->qss);
        break;
    }
    default:
        luaL_argcheck(L, 0, 2, "Unexpected gather argument");
        break;
    }
    return 1;
}

static void
gather_or_I(QLA_Int *res, QLA_Int *val, void *env)
{
    QLA_Int r;
    QLA_I_eq_I_or_I(&r, res, val);
    QLA_I_eq_I(res, &r);
}

static int
q_gather_or(lua_State *L)
{
    int top = lua_gettop(L);
    mGather *g = qlua_checkGather(L, 1);
    mLattice *dstLattice;
    mLattice *srcLattice;
    
    luaL_getmetafield(L, 1, dst_lattice_key);
    dstLattice = qlua_checkLattice(L, top + 1);
    luaL_getmetafield(L, 1, src_lattice_key);
    srcLattice = qlua_checkLattice(L, top + 2);

    switch (qlua_qtype(L, 2)) {
    case qLatInt: {
        QLA_Int zval = 0;
        mLatInt *src = qlua_checkLatInt(L, 2, srcLattice);
        mLatInt *dst = qlua_newZeroLatInt(L, top + 1);
        CALL_QDP(L);
        QDP_I_eq_collect_I_funca(dst->ptr, g->ptr, src->ptr,
                                 &zval, gather_or_I, NULL,
                                 *srcLattice->qss);
        break;
    }
    default:
        luaL_argcheck(L, 0, 2, "Unexpected gather argument");
        break;
    }
    return 1;
}

static void
gather_xor_I(QLA_Int *res, QLA_Int *val, void *env)
{
    QLA_Int r;
    QLA_I_eq_I_xor_I(&r, res, val);
    QLA_I_eq_I(res, &r);
}

static int
q_gather_xor(lua_State *L)
{
    int top = lua_gettop(L);
    mGather *g = qlua_checkGather(L, 1);
    mLattice *dstLattice;
    mLattice *srcLattice;
    
    luaL_getmetafield(L, 1, dst_lattice_key);
    dstLattice = qlua_checkLattice(L, top + 1);
    luaL_getmetafield(L, 1, src_lattice_key);
    srcLattice = qlua_checkLattice(L, top + 2);

    switch (qlua_qtype(L, 2)) {
    case qLatInt: {
        QLA_Int zval = 0;
        mLatInt *src = qlua_checkLatInt(L, 2, srcLattice);
        mLatInt *dst = qlua_newZeroLatInt(L, top + 1);
        CALL_QDP(L);
        QDP_I_eq_collect_I_funca(dst->ptr, g->ptr, src->ptr,
                                 &zval, gather_xor_I, NULL,
                                 *srcLattice->qss);
        break;
    }
    default:
        luaL_argcheck(L, 0, 2, "Unexpected gather argument");
        break;
    }
    return 1;
}

static void
gather_min_I(QLA_Int *res, QLA_Int *val, void *env)
{
    QLA_Int r;
    QLA_I_eq_I_min_I(&r, res, val);
    QLA_I_eq_I(res, &r);
}

static void
gather_min_R(QLA_Real *res, QLA_Real *val, void *env)
{
    QLA_Real r;
    QLA_R_eq_R_min_R(&r, res, val);
    QLA_R_eq_R(res, &r);
}

static int
q_gather_min(lua_State *L)
{
    int top = lua_gettop(L);
    mGather *g = qlua_checkGather(L, 1);
    mLattice *dstLattice;
    mLattice *srcLattice;
    
    luaL_getmetafield(L, 1, dst_lattice_key);
    dstLattice = qlua_checkLattice(L, top + 1);
    luaL_getmetafield(L, 1, src_lattice_key);
    srcLattice = qlua_checkLattice(L, top + 2);

    switch (qlua_qtype(L, 2)) {
    case qLatInt: {
        QLA_Int zval = QLA_INT_MAX;
        mLatInt *src = qlua_checkLatInt(L, 2, srcLattice);
        mLatInt *dst = qlua_newZeroLatInt(L, top + 1);
        CALL_QDP(L);
        QDP_I_eq_collect_I_funca(dst->ptr, g->ptr, src->ptr,
                                 &zval, gather_min_I, NULL,
                                 *srcLattice->qss);
        break;
    }
    case qLatReal: {
        QLA_D_Real zval = QLA_D_REAL_MAX;
        mLatReal *src = qlua_checkLatReal(L, 2, srcLattice);
        mLatReal *dst = qlua_newZeroLatReal(L, top + 1);
        CALL_QDP(L);
        QDP_D_R_eq_collect_R_funca(dst->ptr, g->ptr, src->ptr,
                                   &zval, gather_min_R, NULL,
                                   *srcLattice->qss);
        break;
    }
    default:
        luaL_argcheck(L, 0, 2, "Unexpected gather argument");
        break;
    }
    return 1;
}

static void
gather_max_I(QLA_Int *res, QLA_Int *val, void *env)
{
    QLA_Int r;
    QLA_I_eq_I_max_I(&r, res, val);
    QLA_I_eq_I(res, &r);
}

static void
gather_max_R(QLA_Real *res, QLA_Real *val, void *env)
{
    QLA_Real r;
    QLA_R_eq_R_max_R(&r, res, val);
    QLA_R_eq_R(res, &r);
}

static int
q_gather_max(lua_State *L)
{
    int top = lua_gettop(L);
    mGather *g = qlua_checkGather(L, 1);
    mLattice *dstLattice;
    mLattice *srcLattice;
    
    luaL_getmetafield(L, 1, dst_lattice_key);
    dstLattice = qlua_checkLattice(L, top + 1);
    luaL_getmetafield(L, 1, src_lattice_key);
    srcLattice = qlua_checkLattice(L, top + 2);

    switch (qlua_qtype(L, 2)) {
    case qLatInt: {
        QLA_Int zval = QLA_INT_MIN;
        mLatInt *src = qlua_checkLatInt(L, 2, srcLattice);
        mLatInt *dst = qlua_newZeroLatInt(L, top + 1);
        CALL_QDP(L);
        QDP_I_eq_collect_I_funca(dst->ptr, g->ptr, src->ptr,
                                 &zval, gather_max_I, NULL,
                                 *srcLattice->qss);
        break;
    }
    case qLatReal: {
        QLA_D_Real zval = QLA_D_REAL_MIN;
        mLatReal *src = qlua_checkLatReal(L, 2, srcLattice);
        mLatReal *dst = qlua_newZeroLatReal(L, top + 1);
        CALL_QDP(L);
        QDP_D_R_eq_collect_R_funca(dst->ptr, g->ptr, src->ptr,
                                   &zval, gather_max_R, NULL,
                                   *srcLattice->qss);
        break;
    }
    default:
        luaL_argcheck(L, 0, 2, "Unexpected gather argument");
        break;
    }
    return 1;
}

static void
gather_mul_I(QLA_Int *res, QLA_Int *val, void *env)
{
    QLA_Int r;
    QLA_I_eq_I_times_I(&r, res, val);
    QLA_I_eq_I(res, &r);
}

static void
gather_mul_R(QLA_Real *res, QLA_Real *val, void *env)
{
    QLA_Real r;
    QLA_R_eq_R_times_R(&r, res, val);
    QLA_R_eq_R(res, &r);
}

static void
gather_mul_C(QLA_Complex *res, QLA_Complex *val, void *env)
{
    QLA_Complex r;
    QLA_C_eq_C_times_C(&r, res, val);
    QLA_C_eq_C(res, &r);
}

static int
q_gather_mul(lua_State *L)
{
    int top = lua_gettop(L);
    mGather *g = qlua_checkGather(L, 1);
    mLattice *dstLattice;
    mLattice *srcLattice;
    
    luaL_getmetafield(L, 1, dst_lattice_key);
    dstLattice = qlua_checkLattice(L, top + 1);
    luaL_getmetafield(L, 1, src_lattice_key);
    srcLattice = qlua_checkLattice(L, top + 2);

    switch (qlua_qtype(L, 2)) {
    case qLatInt: {
        QLA_Int zval = 1;
        mLatInt *src = qlua_checkLatInt(L, 2, srcLattice);
        mLatInt *dst = qlua_newZeroLatInt(L, top + 1);
        CALL_QDP(L);
        QDP_I_eq_collect_I_funca(dst->ptr, g->ptr, src->ptr,
                                 &zval, gather_mul_I, NULL,
                                 *srcLattice->qss);
        break;
    }
    case qLatReal: {
        QLA_Real zval = 1.0;
        mLatReal *src = qlua_checkLatReal(L, 2, srcLattice);
        mLatReal *dst = qlua_newZeroLatReal(L, top + 1);
        CALL_QDP(L);
        QDP_D_R_eq_collect_R_funca(dst->ptr, g->ptr, src->ptr,
                                   &zval, gather_mul_R, NULL,
                                   *srcLattice->qss);
        break;
    }
    case qLatComplex: {
        QLA_Real rone = 1.0;
        QLA_Complex zval;
        mLatComplex *src = qlua_checkLatComplex(L, 2, srcLattice);
        mLatComplex *dst = qlua_newZeroLatComplex(L, top + 1);
        QLA_C_eq_R(&zval, &rone);
        CALL_QDP(L);
        QDP_D_C_eq_collect_C_funca(dst->ptr, g->ptr, src->ptr,
                                   &zval, gather_mul_C, NULL,
                                   *srcLattice->qss);
        break;
    }
    default:
        luaL_argcheck(L, 0, 2, "Unexpected gather argument");
        break;
    }
    return 1;
}

static void
gather_add_I(QLA_Int *res, QLA_Int *val, void *env)
{
    QLA_Int r;
    QLA_I_eq_I_plus_I(&r, res, val);
    QLA_I_eq_I(res, &r);
}

static void
gather_add_R(QLA_Real *res, QLA_Real *val, void *env)
{
    QLA_Real r;
    QLA_R_eq_R_plus_R(&r, res, val);
    QLA_R_eq_R(res, &r);
}

static void
gather_add_C(QLA_Complex *res, QLA_Complex *val, void *env)
{
    QLA_Complex r;
    QLA_C_eq_C_plus_C(&r, res, val);
    QLA_C_eq_C(res, &r);
}

#if USE_Nc2
static void
gather_add_V2(QLA_D2_ColorVector *res, QLA_D2_ColorVector *val, void *env)
{
    QLA_D2_ColorVector r;
    QLA_D2_V_eq_V_plus_V(&r, res, val);
    QLA_D2_V_eq_V(res, &r);
}

static void
gather_add_M2(QLA_D2_ColorMatrix *res, QLA_D2_ColorMatrix *val, void *env)
{
    QLA_D2_ColorMatrix r;
    QLA_D2_M_eq_M_plus_M(&r, res, val);
    QLA_D2_M_eq_M(res, &r);
}

static void
gather_add_D2(QLA_D2_DiracFermion *res, QLA_D2_DiracFermion *val, void *env)
{
    QLA_D2_DiracFermion r;
    QLA_D2_D_eq_D_plus_D(&r, res, val);
    QLA_D2_D_eq_D(res, &r);
}

static void
gather_add_P2(QLA_D2_DiracPropagator *res, QLA_D2_DiracPropagator *val, void *env)
{
    QLA_D2_DiracPropagator r;
    QLA_D2_P_eq_P_plus_P(&r, res, val);
    QLA_D2_P_eq_P(res, &r);
}
#endif /* USE_Nc2 */

#if USE_Nc3
static void
gather_add_V3(QLA_D3_ColorVector *res, QLA_D3_ColorVector *val, void *env)
{
    QLA_D3_ColorVector r;
    QLA_D3_V_eq_V_plus_V(&r, res, val);
    QLA_D3_V_eq_V(res, &r);
}

static void
gather_add_M3(QLA_D3_ColorMatrix *res, QLA_D3_ColorMatrix *val, void *env)
{
    QLA_D3_ColorMatrix r;
    QLA_D3_M_eq_M_plus_M(&r, res, val);
    QLA_D3_M_eq_M(res, &r);
}

static void
gather_add_D3(QLA_D3_DiracFermion *res, QLA_D3_DiracFermion *val, void *env)
{
    QLA_D3_DiracFermion r;
    QLA_D3_D_eq_D_plus_D(&r, res, val);
    QLA_D3_D_eq_D(res, &r);
}

static void
gather_add_P3(QLA_D3_DiracPropagator *res, QLA_D3_DiracPropagator *val, void *env)
{
    QLA_D3_DiracPropagator r;
    QLA_D3_P_eq_P_plus_P(&r, res, val);
    QLA_D3_P_eq_P(res, &r);
}
#endif /* USE_Nc3 */

#if USE_NcN
static void
gather_add_VN(int nc, QLA_DN_ColorVector(nc,(*res)), QLA_DN_ColorVector(nc,(*val)), void *env)
{
    QLA_DN_ColorVector(nc,(r));
    QLA_DN_V_eq_V_plus_V(nc, &r, res, val);
    QLA_DN_V_eq_V(nc, res, &r);
}

static void
gather_add_MN(int nc, QLA_DN_ColorMatrix(nc,(*res)), QLA_DN_ColorMatrix(nc,(*val)), void *env)
{
    QLA_DN_ColorMatrix(nc,(r));
    QLA_DN_M_eq_M_plus_M(nc, &r, res, val);
    QLA_DN_M_eq_M(nc, res, &r);
}

static void
gather_add_DN(int nc, QLA_DN_DiracFermion(nc,(*res)), QLA_DN_DiracFermion(nc,(*val)), void *env)
{
    QLA_DN_DiracFermion(nc,(r));
    QLA_DN_D_eq_D_plus_D(nc, &r, res, val);
    QLA_DN_D_eq_D(nc, res, &r);
}

static void
gather_add_PN(int nc, QLA_DN_DiracPropagator(nc,(*res)), QLA_DN_DiracPropagator(nc,(*val)), void *env)
{
    QLA_DN_DiracPropagator(nc,(r));
    QLA_DN_P_eq_P_plus_P(nc, &r, res, val);
    QLA_DN_P_eq_P(nc, res, &r);
}
#endif /* USE_NcN */

static int
q_gather_add(lua_State *L)
{
    int top = lua_gettop(L);
    mGather *g = qlua_checkGather(L, 1);
    mLattice *dstLattice;
    mLattice *srcLattice;
    
    luaL_getmetafield(L, 1, dst_lattice_key);
    dstLattice = qlua_checkLattice(L, top + 1);
    luaL_getmetafield(L, 1, src_lattice_key);
    srcLattice = qlua_checkLattice(L, top + 2);

    switch (qlua_qtype(L, 2)) {
    case qLatInt: {
        mLatInt *src = qlua_checkLatInt(L, 2, srcLattice);
        mLatInt *dst = qlua_newZeroLatInt(L, top + 1);
        QLA_Int zval = 0;
        CALL_QDP(L);
        QDP_I_eq_collect_I_funca(dst->ptr, g->ptr, src->ptr,
                                 &zval, gather_add_I, NULL,
                                 *srcLattice->qss);
        break;
    }
    case qLatReal: {
        mLatReal *src = qlua_checkLatReal(L, 2, srcLattice);
        mLatReal *dst = qlua_newZeroLatReal(L, top + 1);
        QLA_Real zval = 0;
        CALL_QDP(L);
        QDP_D_R_eq_collect_R_funca(dst->ptr, g->ptr, src->ptr,
                                   &zval, gather_add_R, NULL,
                                   *srcLattice->qss);
        break;
    }
    case qLatComplex: {
        mLatComplex *src = qlua_checkLatComplex(L, 2, srcLattice);
        mLatComplex *dst = qlua_newZeroLatComplex(L, top + 1);
        QLA_Complex zval;
        QLA_C_eq_zero(&zval);
        CALL_QDP(L);
        QDP_D_C_eq_collect_C_funca(dst->ptr, g->ptr, src->ptr,
                                   &zval, gather_add_C, NULL,
                                   *srcLattice->qss);
        break;
    }
#if USE_Nc2
    case qLatColVec2: {
        mLatColVec2 *src = qlua_checkLatColVec2(L, 2, srcLattice, 2);
        mLatColVec2 *dst = qlua_newZeroLatColVec2(L, top + 1, 2);
        QLA_D2_ColorVector zval;
        QLA_D2_V_eq_zero(&zval);
        CALL_QDP(L);
        QDP_D2_V_eq_collect_V_funca(dst->ptr, g->ptr, src->ptr,
                                    &zval, gather_add_V2, NULL,
                                    *srcLattice->qss);
        break;
    }
    case qLatColMat2: {
        mLatColMat2 *src = qlua_checkLatColMat2(L, 2, srcLattice, 2);
        mLatColMat2 *dst = qlua_newZeroLatColMat2(L, top + 1, 2);
        QLA_D2_ColorMatrix zval;
        QLA_D2_M_eq_zero(&zval);
        CALL_QDP(L);
        QDP_D2_M_eq_collect_M_funca(dst->ptr, g->ptr, src->ptr,
                                    &zval, gather_add_M2, NULL,
                                    *srcLattice->qss);
        break;
    }
    case qLatDirFerm2: {
        mLatDirFerm2 *src = qlua_checkLatDirFerm2(L, 2, srcLattice, 2);
        mLatDirFerm2 *dst = qlua_newZeroLatDirFerm2(L, top + 1, 2);
        QLA_D2_DiracFermion zval;
        QLA_D2_D_eq_zero(&zval);
        CALL_QDP(L);
        QDP_D2_D_eq_collect_D_funca(dst->ptr, g->ptr, src->ptr,
                                    &zval, gather_add_D2, NULL,
                                    *srcLattice->qss);
        break;
    }
    case qLatDirProp2: {
        mLatDirProp2 *src = qlua_checkLatDirProp2(L, 2, srcLattice, 2);
        mLatDirProp2 *dst = qlua_newZeroLatDirProp2(L, top + 1, 2);
        QLA_D2_DiracPropagator zval;
        QLA_D2_P_eq_zero(&zval);
        CALL_QDP(L);
        QDP_D2_P_eq_collect_P_funca(dst->ptr, g->ptr, src->ptr,
                                    &zval, gather_add_P2, NULL,
                                    *srcLattice->qss);
        break;
    }
#endif /* USE_Nc2 */
#if USE_Nc3
    case qLatColVec3: {
        mLatColVec3 *src = qlua_checkLatColVec3(L, 2, srcLattice, 3);
        mLatColVec3 *dst = qlua_newZeroLatColVec3(L, top + 1, 3);
        QLA_D3_ColorVector zval;
        QLA_D3_V_eq_zero(&zval);
        CALL_QDP(L);
        QDP_D3_V_eq_collect_V_funca(dst->ptr, g->ptr, src->ptr,
                                    &zval, gather_add_V3, NULL,
                                    *srcLattice->qss);
        break;
    }
    case qLatColMat3: {
        mLatColMat3 *src = qlua_checkLatColMat3(L, 2, srcLattice, 3);
        mLatColMat3 *dst = qlua_newZeroLatColMat3(L, top + 1, 3);
        QLA_D3_ColorMatrix zval;
        QLA_D3_M_eq_zero(&zval);
        CALL_QDP(L);
        QDP_D3_M_eq_collect_M_funca(dst->ptr, g->ptr, src->ptr,
                                    &zval, gather_add_M3, NULL,
                                    *srcLattice->qss);
        break;
    }
    case qLatDirFerm3: {
        mLatDirFerm3 *src = qlua_checkLatDirFerm3(L, 2, srcLattice, 3);
        mLatDirFerm3 *dst = qlua_newZeroLatDirFerm3(L, top + 1, 3);
        QLA_D3_DiracFermion zval;
        QLA_D3_D_eq_zero(&zval);
        CALL_QDP(L);
        QDP_D3_D_eq_collect_D_funca(dst->ptr, g->ptr, src->ptr,
                                    &zval, gather_add_D3, NULL,
                                    *srcLattice->qss);
        break;
    }
    case qLatDirProp3: {
        mLatDirProp3 *src = qlua_checkLatDirProp3(L, 2, srcLattice, 3);
        mLatDirProp3 *dst = qlua_newZeroLatDirProp3(L, top + 1, 3);
        QLA_D3_DiracPropagator zval;
        QLA_D3_P_eq_zero(&zval);
        CALL_QDP(L);
        QDP_D3_P_eq_collect_P_funca(dst->ptr, g->ptr, src->ptr,
                                    &zval, gather_add_P3, NULL,
                                    *srcLattice->qss);
        break;
    }
#endif /* USE_Nc3 */
#if USE_NcN
    case qLatColVecN: {
        mLatColVecN *src = qlua_checkLatColVecN(L, 2, srcLattice, -1);
        int nc = src->nc;
        mLatColVecN *dst = qlua_newZeroLatColVecN(L, top + 1, nc);
        QLA_DN_ColorVector(nc, (zval));
        QLA_DN_V_eq_zero(nc, &zval);
        CALL_QDP(L);
        QDP_DN_V_eq_collect_V_funca(dst->ptr, g->ptr, src->ptr,
                                    &zval, gather_add_VN, NULL,
                                    *srcLattice->qss);
        break;
    }
    case qLatColMatN: {
        mLatColMatN *src = qlua_checkLatColMatN(L, 2, srcLattice, -1);
        int nc = src->nc;
        mLatColMatN *dst = qlua_newZeroLatColMatN(L, top + 1, nc);
        QLA_DN_ColorMatrix(nc,(zval));
        QLA_DN_M_eq_zero(nc, &zval);
        CALL_QDP(L);
        QDP_DN_M_eq_collect_M_funca(dst->ptr, g->ptr, src->ptr,
                                    &zval, gather_add_MN, NULL,
                                    *srcLattice->qss);
        break;
    }
    case qLatDirFermN: {
        mLatDirFermN *src = qlua_checkLatDirFermN(L, 2, srcLattice, -1);
        int nc = src->nc;
        mLatDirFermN *dst = qlua_newZeroLatDirFermN(L, top + 1, nc);
        QLA_DN_DiracFermion(nc,(zval));
        QLA_DN_D_eq_zero(nc, &zval);
        CALL_QDP(L);
        QDP_DN_D_eq_collect_D_funca(dst->ptr, g->ptr, src->ptr,
                                    &zval, gather_add_DN, NULL,
                                    *srcLattice->qss);
        break;
    }
    case qLatDirPropN: {
        mLatDirPropN *src = qlua_checkLatDirPropN(L, 2, srcLattice, -1);
        int nc = src->nc;
        mLatDirPropN *dst = qlua_newZeroLatDirPropN(L, top + 1, nc);
        QLA_DN_DiracPropagator(nc,(zval));
        QLA_DN_P_eq_zero(nc, &zval);
        CALL_QDP(L);
        QDP_DN_P_eq_collect_P_funca(dst->ptr, g->ptr, src->ptr,
                                    &zval, gather_add_PN, NULL,
                                    *srcLattice->qss);
        break;
    }
#endif /* USE_NcN */
    default:
        luaL_argcheck(L, 0, 2, "Unexpected gather argument");
        break;
    }
    return 1;
}

static struct luaL_Reg GatherMethods[] = {
    { "__tostring", q_gather_fmt  },
    { "__gc",       q_gather_gc   },
    { "or",         q_gather_or },
    { "xor",        q_gather_xor },
    { "and",        q_gather_and },
    { "min",        q_gather_min },
    { "max",        q_gather_max },
    { "mul",        q_gather_mul },
    { "add",        q_gather_add },
    { NULL, NULL }
};

static int
q_gather(lua_State *L)
{
    mLattice *dstLattice = qlua_checkLattice(L, 1);
    mLattice *srcLattice = qlua_checkLattice(L, 2);
    int dstRank = dstLattice->rank;
    QDP_Int **addr = qlua_malloc(L, dstRank * sizeof (QDP_Int *));
    int i;
    mGather *g;

    qlua_checktable(L, 3, "gather address");
    for (i = 0; i < dstRank; i++) {
        mLatInt *ai;
        lua_pushnumber(L, i + 1);
        lua_gettable(L, 3);
        ai = qlua_checkLatInt(L, -1, srcLattice);
        addr[i] = ai->ptr;
        lua_pop(L, 1);
    }
    g = qlua_newGather(L);
    CALL_QDP(L);
    g->ptr = QDP_create_collect(dstLattice->lat, srcLattice->lat, addr);
    qlua_free(L, addr);
    
    lua_createtable(L, 0, 12);
    qlua_fillmeta(L, GatherMethods, qGather);
    /* Must keep references to both lattices to avoid Lua gc while gather is alive */
    lua_pushvalue(L, 2);
    lua_setfield(L, -2, src_lattice_key);
    lua_pushvalue(L, 1);
    lua_setfield(L, -2, dst_lattice_key);
    lua_setmetatable(L, -2);

    return 1;
}

static struct luaL_reg fGather[] = {
    { "gather",    q_gather },
    { NULL, NULL}
};

int
init_gather(lua_State *L)
{
    luaL_register(L, qcdlib, fGather);
    return 0;
}

int
fini_gather(lua_State *L)
{
    return 0;
}
