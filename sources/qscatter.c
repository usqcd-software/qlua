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
#include "qscatter.h"                                                /* DEPS */
#include "qdp.h"

static const char ScatterName[] = "qlua.scatter";
static const char src_lattice_key[]  = "scatter.src.lattice";
static const char dst_lattice_key[]  = "scatter.dst.lattice";

static int
q_scatter_fmt(lua_State *L)
{
    char fmt[72];
    mScatter *b = qlua_checkScatter(L, 1);

    sprintf(fmt, "QDP:Scatter(%p)", b->ptr);
    lua_pushstring(L, fmt);

    return 1;
}

static int
q_scatter_gc(lua_State *L)
{
    mScatter *b = qlua_checkScatter(L, 1);

    QDP_destroy_scatter(b->ptr);
    b->ptr = 0;

    return 0;
}

mScatter *
qlua_newScatter(lua_State *L)
{
    mScatter *s;

    s = lua_newuserdata(L, sizeof (mScatter));
    s->ptr = 0;

    return s;
}

mScatter *
qlua_checkScatter(lua_State *L, int idx)
{
    void *v = qlua_checkLatticeType(L, idx, qScatter, ScatterName);
    return v;
}

static int
q_scatter_call(lua_State *L)
{
    int top = lua_gettop(L);
    mScatter *s = qlua_checkScatter(L, 1);
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
        CALL_QDP(L);
        QDP_I_eq_scatter_I(dst->ptr, s->ptr, src->ptr, *dstLattice->qss);
        break;
    }
    case qLatReal: {
        mLatReal *src = qlua_checkLatReal(L, 2, srcLattice);
        mLatReal *dst = qlua_newZeroLatReal(L, top + 1);
        CALL_QDP(L);
        QDP_R_eq_scatter_R(dst->ptr, s->ptr, src->ptr, *dstLattice->qss);
        break;
    }
    case qLatComplex: {
        mLatComplex *src = qlua_checkLatComplex(L, 2, srcLattice);
        mLatComplex *dst = qlua_newZeroLatComplex(L, top + 1);
        CALL_QDP(L);
        QDP_C_eq_scatter_C(dst->ptr, s->ptr, src->ptr, *dstLattice->qss);
        break;
    }
#if USE_Nc2
    case qLatColVec2: {
        mLatColVec2 *src = qlua_checkLatColVec2(L, 2, srcLattice, 2);
        mLatColVec2 *dst = qlua_newZeroLatColVec2(L, top + 1, 2);
        CALL_QDP(L);
        QDP_D2_V_eq_scatter_V(dst->ptr, s->ptr, src->ptr, *dstLattice->qss);
        break;
    }
    case qLatColMat2: {
        mLatColMat2 *src = qlua_checkLatColMat2(L, 2, srcLattice, 2);
        mLatColMat2 *dst = qlua_newZeroLatColMat2(L, top + 1, 2);
        CALL_QDP(L);
        QDP_D2_M_eq_scatter_M(dst->ptr, s->ptr, src->ptr, *dstLattice->qss);
        break;
    }
    case qLatDirFerm2: {
        mLatDirFerm2 *src = qlua_checkLatDirFerm2(L, 2, srcLattice, 2);
        mLatDirFerm2 *dst = qlua_newZeroLatDirFerm2(L, top + 1, 2);
        CALL_QDP(L);
        QDP_D2_D_eq_scatter_D(dst->ptr, s->ptr, src->ptr, *dstLattice->qss);
        break;
    }
    case qLatDirProp2: {
        mLatDirProp2 *src = qlua_checkLatDirProp2(L, 2, srcLattice, 2);
        mLatDirProp2 *dst = qlua_newZeroLatDirProp2(L, top + 1, 2);
        CALL_QDP(L);
        QDP_D2_P_eq_scatter_P(dst->ptr, s->ptr, src->ptr, *dstLattice->qss);
        break;
    }
#endif
#if USE_Nc3
    case qLatColVec3: {
        mLatColVec3 *src = qlua_checkLatColVec3(L, 2, srcLattice, 3);
        mLatColVec3 *dst = qlua_newZeroLatColVec3(L, top + 1, 3);
        CALL_QDP(L);
        QDP_D3_V_eq_scatter_V(dst->ptr, s->ptr, src->ptr, *dstLattice->qss);
        break;
    }
    case qLatColMat3: {
        mLatColMat3 *src = qlua_checkLatColMat3(L, 2, srcLattice, 3);
        mLatColMat3 *dst = qlua_newZeroLatColMat3(L, top + 1, 3);
        CALL_QDP(L);
        QDP_D3_M_eq_scatter_M(dst->ptr, s->ptr, src->ptr, *dstLattice->qss);
        break;
    }
    case qLatDirFerm3: {
        mLatDirFerm3 *src = qlua_checkLatDirFerm3(L, 2, srcLattice, 3);
        mLatDirFerm3 *dst = qlua_newZeroLatDirFerm3(L, top + 1, 3);
        CALL_QDP(L);
        QDP_D3_D_eq_scatter_D(dst->ptr, s->ptr, src->ptr, *dstLattice->qss);
        break;
    }
    case qLatDirProp3: {
        mLatDirProp3 *src = qlua_checkLatDirProp3(L, 2, srcLattice, 3);
        mLatDirProp3 *dst = qlua_newZeroLatDirProp3(L, top + 1, 3);
        CALL_QDP(L);
        QDP_D3_P_eq_scatter_P(dst->ptr, s->ptr, src->ptr, *dstLattice->qss);
        break;
    }
#endif
#if USE_NcN
    case qLatColVecN: {
        mLatColVecN *src = qlua_checkLatColVecN(L, 2, srcLattice, -1);
        int nc = src->nc;
        mLatColVecN *dst = qlua_newZeroLatColVecN(L, top + 1, nc);
        CALL_QDP(L);
        QDP_DN_V_eq_scatter_V(dst->ptr, s->ptr, src->ptr, *dstLattice->qss);
        break;
    }
    case qLatColMatN: {
        mLatColMatN *src = qlua_checkLatColMatN(L, 2, srcLattice, -1);
        int nc = src->nc;
        mLatColMatN *dst = qlua_newZeroLatColMatN(L, top + 1, nc);
        CALL_QDP(L);
        QDP_DN_M_eq_scatter_M(dst->ptr, s->ptr, src->ptr, *dstLattice->qss);
        break;
    }
    case qLatDirFermN: {
        mLatDirFermN *src = qlua_checkLatDirFermN(L, 2, srcLattice, -1);
        int nc = src->nc;
        mLatDirFermN *dst = qlua_newZeroLatDirFermN(L, top + 1, nc);
        CALL_QDP(L);
        QDP_DN_D_eq_scatter_D(dst->ptr, s->ptr, src->ptr, *dstLattice->qss);
        break;
    }
    case qLatDirPropN: {
        mLatDirPropN *src = qlua_checkLatDirPropN(L, 2, srcLattice, -1);
        int nc = src->nc;
        mLatDirPropN *dst = qlua_newZeroLatDirPropN(L, top + 1, nc);
        CALL_QDP(L);
        QDP_DN_P_eq_scatter_P(dst->ptr, s->ptr, src->ptr, *dstLattice->qss);
        break;
    }
#endif
    default:
        luaL_argcheck(L, 0, 2, "Unexpected scatter argument");
        break;
    }
    return 1;
}

static struct luaL_Reg ScatterMethods[] = {
    { "__tostring", q_scatter_fmt  },
    { "__gc",       q_scatter_gc   },
    { "__call",     q_scatter_call },
    { NULL, NULL }
};


static int
q_scatter(lua_State *L)
{
    mLattice *dstLattice = qlua_checkLattice(L, 1);
    mLattice *srcLattice = qlua_checkLattice(L, 2);
    int srcRank = srcLattice->rank;
    QDP_Int **addr = qlua_malloc(L, srcRank * sizeof (QDP_Int *));
    int i;
    mScatter *s;

    qlua_checktable(L, 3, "scatter address");
    for (i = 0; i < srcRank; i++) {
        mLatInt *ai;
        lua_pushnumber(L, i + 1);
        lua_gettable(L, 3);
        ai = qlua_checkLatInt(L, -1, dstLattice);
        addr[i] = ai->ptr;
        lua_pop(L, 1);
    }
    s = qlua_newScatter(L);
    CALL_QDP(L);
    s->ptr = QDP_create_scatter(dstLattice->lat, srcLattice->lat, addr);
    qlua_free(L, addr);
    
    lua_createtable(L, 0, 6);
    qlua_fillmeta(L, ScatterMethods, qScatter);
    /* Must keep references to both lattices to avoid Lua gc while scatter is alive */
    lua_pushvalue(L, 2);
    lua_setfield(L, -2, src_lattice_key);
    lua_pushvalue(L, 1);
    lua_setfield(L, -2, dst_lattice_key);
    lua_setmetatable(L, -2);

    return 1;
}

static struct luaL_Reg fScatter[] = {
    { "scatter",            q_scatter },
    { NULL, NULL}
};


int
init_scatter(lua_State *L)
{
    luaL_register(L, qcdlib, fScatter);
    return 0;
}

void
fini_scatter(void)
{
}
