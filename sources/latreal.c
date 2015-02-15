#include "qlua.h"                                                    /* DEPS */
#include "qvector.h"                                                 /* DEPS */
#include "lattice.h"                                                 /* DEPS */
#include "latsubset.h"                                               /* DEPS */
#include "latmulti.h"                                                /* DEPS */
#include "latreal.h"                                                 /* DEPS */
#include "latint.h"                                                  /* DEPS */
#include "latrandom.h"                                               /* DEPS */
#include "latcomplex.h"                                              /* DEPS */
#include "qmp.h"

const char LatRealName[] = "lattice.Real";

static int
q_R_fmt(lua_State *L)
{
    char fmt[72];
    mLatReal *b = qlua_checkLatReal(L, 1, NULL);

    sprintf(fmt, "LatReal(%p)", b->ptr);
    lua_pushstring(L, fmt);

    return 1;
}

static int
q_R_gc(lua_State *L)
{
    mLatReal *b = qlua_checkLatReal(L, 1, NULL);

    QDP_destroy_R(b->ptr);
    b->ptr = 0;
    qlua_qdp_memuse(L, "Real", -1);

    return 0;
}

static int
q_R_neg(lua_State *L)
{
    mLatReal *a = qlua_checkLatReal(L, 1, NULL);
    mLattice *S = qlua_ObjLattice(L, 1);
    mLatReal *res = qlua_newLatReal(L, lua_gettop(L));

    CALL_QDP(L);
    QDP_R_eqm_R(res->ptr, a->ptr, *S->qss);
    return 1;
}

static int
q_R_sum(lua_State *L)
{
    mLatReal *a = qlua_checkLatReal(L, 1, NULL);
    int argc = lua_gettop(L);
    mLattice *S = qlua_ObjLattice(L, 1);
    int Sidx = lua_gettop(L);
    
    switch (argc) {
    case 1: {
        QLA_Real sum;

        CALL_QDP(L);
        if (S->lss.mask) {
            mLatReal *b = qlua_newZeroLatReal(L, Sidx);
            QDP_R_eq_R_mask_I(b->ptr, a->ptr, S->lss.mask, *S->qss);
            QDP_r_eq_sum_R(&sum, b->ptr, *S->qss);
            lua_pop(L, 1);
        } else {
            QDP_r_eq_sum_R(&sum, a->ptr, *S->qss);
        }
        lua_pushnumber(L, sum);

        return 1;
    }
    case 2: {
        mLatMulti *m = qlua_checkLatMulti(L, 2, S);
        int size = m->size;
        QLA_Int *ii = m->idx;
        mVecReal *r = qlua_newVecReal(L, size);
        int sites = QDP_sites_on_node_L(S->lat);
        int k;
        QLA_Real *xx;
        
        for (k = 0; k < size; k++)
            r->val[k] = 0;

        CALL_QDP(L);
        xx = QDP_expose_R(a->ptr);
        for (k = 0; k < sites; k++, xx++, ii++) {
            int t = *ii;
            if ((t < 0) || (t >= size))
                continue;
            r->val[t] += *xx;
        }
        QDP_reset_R(a->ptr);
        QMP_sum_double_array(r->val, size);

        return 1;
    }
    }
    return luaL_error(L, "bad arguments for Real:sum()");
}

static int
q_R_norm2(lua_State *L)
{
    mLatReal *a = qlua_checkLatReal(L, 1, NULL);
    mLattice *S = qlua_ObjLattice(L, 1);
    QLA_Real n;

    CALL_QDP(L);
    if (S->lss.mask) {
        mLatReal *b = qlua_newLatReal(L, lua_gettop(L));
        QDP_R_eq_R_mask_I(b->ptr, a->ptr, S->lss.mask, *S->qss);
        QDP_r_eq_norm2_R(&n, b->ptr, *S->qss);
    } else {
        QDP_r_eq_norm2_R(&n, a->ptr, *S->qss);
    }
    lua_pushnumber(L, n);

    return 1;
}

static int
q_R_shift(lua_State *L)
{
    mLatReal *a = qlua_checkLatReal(L, 1, NULL);
    mLattice *S = qlua_ObjLattice(L, 1);
    int Sidx = lua_gettop(L);
    QDP_Shift shift = qlua_checkShift(L, 2, S);
    QDP_ShiftDir dir = qlua_checkShiftDir(L, 3);
    mLatReal *r = qlua_newLatReal(L, Sidx);

    CALL_QDP(L);
    QDP_R_eq_sR(r->ptr, a->ptr, shift, dir, *S->qss);

    return 1;
}

static int
q_R_sin(lua_State *L)
{
    mLatReal *a = qlua_checkLatReal(L, 1, NULL);
    mLattice *S = qlua_ObjLattice(L, 1);
    mLatReal *r = qlua_newLatReal(L, lua_gettop(L));

    CALL_QDP(L);
    QDP_R_eq_sin_R(r->ptr, a->ptr, *S->qss);

    return 1;
}

static int
q_R_cos(lua_State *L)
{
    mLatReal *a = qlua_checkLatReal(L, 1, NULL);
    mLattice *S = qlua_ObjLattice(L, 1);
    mLatReal *r = qlua_newLatReal(L, lua_gettop(L));

    CALL_QDP(L);
    QDP_R_eq_cos_R(r->ptr, a->ptr, *S->qss);

    return 1;
}

static int
q_R_tan(lua_State *L)
{
    mLatReal *a = qlua_checkLatReal(L, 1, NULL);
    mLattice *S = qlua_ObjLattice(L, 1);
    mLatReal *r = qlua_newLatReal(L, lua_gettop(L));

    CALL_QDP(L);
    QDP_R_eq_tan_R(r->ptr, a->ptr, *S->qss);

    return 1;
}

static int
q_R_asin(lua_State *L)
{
    mLatReal *a = qlua_checkLatReal(L, 1, NULL);
    mLattice *S = qlua_ObjLattice(L, 1);
    mLatReal *r = qlua_newLatReal(L, lua_gettop(L));

    CALL_QDP(L);
    QDP_R_eq_asin_R(r->ptr, a->ptr, *S->qss);

    return 1;
}

static int
q_R_acos(lua_State *L)
{
    mLatReal *a = qlua_checkLatReal(L, 1, NULL);
    mLattice *S = qlua_ObjLattice(L, 1);
    mLatReal *r = qlua_newLatReal(L, lua_gettop(L));

    CALL_QDP(L);
    QDP_R_eq_acos_R(r->ptr, a->ptr, *S->qss);

    return 1;
}

static int
q_R_atan(lua_State *L)
{
    mLatReal *a = qlua_checkLatReal(L, 1, NULL);
    mLattice *S = qlua_ObjLattice(L, 1);
    mLatReal *r = qlua_newLatReal(L, lua_gettop(L));

    CALL_QDP(L);
    QDP_R_eq_atan_R(r->ptr, a->ptr, *S->qss);

    return 1;
}

static int
q_R_sqrt(lua_State *L)
{
    mLatReal *a = qlua_checkLatReal(L, 1, NULL);
    mLattice *S = qlua_ObjLattice(L, 1);
    mLatReal *r = qlua_newLatReal(L, lua_gettop(L));

    CALL_QDP(L);
    QDP_R_eq_sqrt_R(r->ptr, a->ptr, *S->qss);

    return 1;
}

static int
q_R_abs(lua_State *L)
{
    mLatReal *a = qlua_checkLatReal(L, 1, NULL);
    mLattice *S = qlua_ObjLattice(L, 1);
    mLatReal *r = qlua_newLatReal(L, lua_gettop(L));

    CALL_QDP(L);
    QDP_R_eq_fabs_R(r->ptr, a->ptr, *S->qss);

    return 1;
}

static int
q_R_exp(lua_State *L)
{
    mLatReal *a = qlua_checkLatReal(L, 1, NULL);
    mLattice *S = qlua_ObjLattice(L, 1);
    mLatReal *r = qlua_newLatReal(L, lua_gettop(L));

    CALL_QDP(L);
    QDP_R_eq_exp_R(r->ptr, a->ptr, *S->qss);

    return 1;
}

static int
q_R_expi(lua_State *L)
{
    mLatReal *a = qlua_checkLatReal(L, 1, NULL);
    mLattice *S = qlua_ObjLattice(L, 1);
    mLatComplex *r = qlua_newLatComplex(L, lua_gettop(L));

    CALL_QDP(L);
    QDP_C_eq_cexpi_R(r->ptr, a->ptr, *S->qss);

    return 1;
}

static int
q_R_log(lua_State *L)
{
    mLatReal *a = qlua_checkLatReal(L, 1, NULL);
    mLattice *S = qlua_ObjLattice(L, 1);
    mLatReal *r = qlua_newLatReal(L, lua_gettop(L));

    CALL_QDP(L);
    QDP_R_eq_log_R(r->ptr, a->ptr, *S->qss);

    return 1;

}

static int
q_R_sign(lua_State *L)
{
    mLatReal *a = qlua_checkLatReal(L, 1, NULL);
    mLattice *S = qlua_ObjLattice(L, 1);
    mLatReal *r = qlua_newLatReal(L, lua_gettop(L));

    CALL_QDP(L);
    QDP_R_eq_sign_R(r->ptr, a->ptr, *S->qss);

    return 1;
}

static int
q_R_floor(lua_State *L)
{
    mLatReal *a = qlua_checkLatReal(L, 1, NULL);
    mLattice *S = qlua_ObjLattice(L, 1);
    mLatReal *r = qlua_newLatReal(L, lua_gettop(L));

    CALL_QDP(L);
    QDP_R_eq_floor_R(r->ptr, a->ptr, *S->qss);

    return 1;
}

static int
q_R_ceil(lua_State *L)
{
    mLatReal *a = qlua_checkLatReal(L, 1, NULL);
    mLattice *S = qlua_ObjLattice(L, 1);
    mLatReal *r = qlua_newLatReal(L, lua_gettop(L));

    CALL_QDP(L);
    QDP_R_eq_ceil_R(r->ptr, a->ptr, *S->qss);

    return 1;
}

static int
q_R_sinh(lua_State *L)
{
    mLatReal *a = qlua_checkLatReal(L, 1, NULL);
    mLattice *S = qlua_ObjLattice(L, 1);
    mLatReal *r = qlua_newLatReal(L, lua_gettop(L));

    CALL_QDP(L);
    QDP_R_eq_sinh_R(r->ptr, a->ptr, *S->qss);

    return 1;
}

static int
q_R_cosh(lua_State *L)
{
    mLatReal *a = qlua_checkLatReal(L, 1, NULL);
    mLattice *S = qlua_ObjLattice(L, 1);
    mLatReal *r = qlua_newLatReal(L, lua_gettop(L));

    CALL_QDP(L);
    QDP_R_eq_cosh_R(r->ptr, a->ptr, *S->qss);

    return 1;
}

static int
q_R_tanh(lua_State *L)
{
    mLatReal *a = qlua_checkLatReal(L, 1, NULL);
    mLattice *S = qlua_ObjLattice(L, 1);
    mLatReal *r = qlua_newLatReal(L, lua_gettop(L));

    CALL_QDP(L);
    QDP_R_eq_tanh_R(r->ptr, a->ptr, *S->qss);

    return 1;

}

static int
q_R_log10(lua_State *L)
{
    mLatReal *a = qlua_checkLatReal(L, 1, NULL);
    mLattice *S = qlua_ObjLattice(L, 1);
    mLatReal *r = qlua_newLatReal(L, lua_gettop(L));

    CALL_QDP(L);
    QDP_R_eq_log10_R(r->ptr, a->ptr, *S->qss);

    return 1;
}

static int
q_R_trunc(lua_State *L)
{
    mLatReal *a = qlua_checkLatReal(L, 1, NULL);
    mLattice *S = qlua_ObjLattice(L, 1);
    mLatInt *r = qlua_newLatInt(L, lua_gettop(L));

    CALL_QDP(L);
    QDP_I_eq_trunc_R(r->ptr, a->ptr, *S->qss);

    return 1;
}

static int
q_R_round(lua_State *L)
{
    mLatReal *a = qlua_checkLatReal(L, 1, NULL);
    mLattice *S = qlua_ObjLattice(L, 1);
    mLatInt *r = qlua_newLatInt(L, lua_gettop(L));

    CALL_QDP(L);
    QDP_I_eq_round_R(r->ptr, a->ptr, *S->qss);

    return 1;
}

static int
q_R_min(lua_State *L)
{
    mLatReal *a = qlua_checkLatReal(L, 1, NULL);
    mLattice *S = qlua_ObjLattice(L, 1);
    mLatReal *b = qlua_checkLatReal(L, 2, S);
    mLatReal *r = qlua_newLatReal(L, lua_gettop(L));

    CALL_QDP(L);
    QDP_R_eq_R_min_R(r->ptr, a->ptr, b->ptr, *S->qss);
    return 1;
}

static int
q_R_max(lua_State *L)
{
    mLatReal *a = qlua_checkLatReal(L, 1, NULL);
    mLattice *S = qlua_ObjLattice(L, 1);
    mLatReal *b = qlua_checkLatReal(L, 2, S);
    mLatReal *r = qlua_newLatReal(L, lua_gettop(L));

    CALL_QDP(L);
    QDP_R_eq_R_max_R(r->ptr, a->ptr, b->ptr, *S->qss);
    return 1;
}

static int
q_R_set(lua_State *L)
{
    mLatReal *r = qlua_checkLatReal(L, 1, NULL);
    mLattice *S = qlua_ObjLattice(L, 1);
    mLatReal *a = qlua_checkLatReal(L, 2, S);

    CALL_QDP(L);
    if (S->lss.mask) {
        QDP_R_eq_R_mask_I(r->ptr, a->ptr, S->lss.mask, *S->qss);
    } else {
        QDP_R_eq_R(r->ptr, a->ptr, *S->qss);
    }
    lua_pop(L, 1);

    return 1;
}

static int
q_R_get(lua_State *L)
{
    switch (qlua_qtype(L, 2)) {
    case qTable: {
        mLatReal *V = qlua_checkLatReal(L, 1, NULL);
        mLattice *S = qlua_ObjLattice(L, 1);
        QLA_Real *locked;
        int *idx = 0;
        double z;

        idx = qlua_checklatcoord(L, 2, S);
        CALL_QDP(L);
        locked = QDP_expose_R(V->ptr);
        int site_node = QDP_node_number_L(S->lat, idx);
        if (site_node == QDP_this_node) {
            z = QLA_elem_R(locked[QDP_index_L(S->lat, idx)]);
        }
        QDP_reset_R(V->ptr);
        qlua_free(L, idx);
        XMP_dist_double_array(site_node, 1, &z);
        lua_pushnumber(L, z);

        return 1;
    }
    case qString:
        return qlua_selflookup(L, 1, luaL_checkstring(L, 2));
    default:
        break;
    }
    return qlua_badindex(L, "Real");
}


static int
q_R_put(lua_State *L)
{
    mLatReal *V = qlua_checkLatReal(L, 1, NULL);
    mLattice *S = qlua_ObjLattice(L, 1);
    QLA_Real *locked;
    int *idx = 0;
    double z = luaL_checknumber(L, 3);

    idx = qlua_checklatcoord(L, 2, S);
    CALL_QDP(L);
    locked = QDP_expose_R(V->ptr);
    if (QDP_node_number_L(S->lat, idx) == QDP_this_node) {
        QLA_elem_R(locked[QDP_index_L(S->lat, idx)]) = z;
    }
    QDP_reset_R(V->ptr);
    qlua_free(L, idx);

    return 0;
}

static int
q_R_add_R(lua_State *L)
{
    mLatReal *a = qlua_checkLatReal(L, 1, NULL);
    mLattice *S = qlua_ObjLattice(L, 1);
    mLatReal *b = qlua_checkLatReal(L, 2, S);
    mLatReal *c = qlua_newLatReal(L, lua_gettop(L));

    CALL_QDP(L);
    QDP_R_eq_R_plus_R(c->ptr, a->ptr, b->ptr, *S->qss);

    return 1;
}

static int
q_r_add_R(lua_State *L)
{
    QLA_D_Real a = luaL_checknumber(L, 1);
    mLatReal *b = qlua_checkLatReal(L, 2, NULL);
    mLattice *S = qlua_ObjLattice(L, 2);
    mLatReal *c = qlua_newLatReal(L, lua_gettop(L));

    CALL_QDP(L);
    QDP_R_eq_r(c->ptr, &a, *S->qss);
    QDP_R_peq_R(c->ptr, b->ptr, *S->qss);

    return 1;
}

static int
q_R_add_r(lua_State *L)
{
    mLatReal *a = qlua_checkLatReal(L, 1, NULL);
    mLattice *S = qlua_ObjLattice(L, 1);
    QLA_D_Real b = luaL_checknumber(L, 2);
    mLatReal *c = qlua_newLatReal(L, lua_gettop(L));

    CALL_QDP(L);
    QDP_R_eq_r(c->ptr, &b, *S->qss);
    QDP_R_peq_R(c->ptr, a->ptr, *S->qss);

    return 1;
}

static int
q_R_sub_R(lua_State *L)
{
    mLatReal *a = qlua_checkLatReal(L, 1, NULL);
    mLattice *S = qlua_ObjLattice(L, 1);
    mLatReal *b = qlua_checkLatReal(L, 2, S);
    mLatReal *c = qlua_newLatReal(L, lua_gettop(L));

    CALL_QDP(L);
    QDP_R_eq_R_minus_R(c->ptr, a->ptr, b->ptr, *S->qss);

    return 1;
}

static int
q_r_sub_R(lua_State *L)
{
    QLA_D_Real a = luaL_checknumber(L, 1);
    mLatReal *b = qlua_checkLatReal(L, 2, NULL);
    mLattice *S = qlua_ObjLattice(L, 2);
    mLatReal *c = qlua_newLatReal(L, lua_gettop(L));

    CALL_QDP(L);
    QDP_R_eq_r(c->ptr, &a, *S->qss);
    QDP_R_meq_R(c->ptr, b->ptr, *S->qss);

    return 1;
}

static int
q_R_sub_r(lua_State *L)
{
    mLatReal *a = qlua_checkLatReal(L, 1, NULL);
    mLattice *S = qlua_ObjLattice(L, 1);
    QLA_D_Real b = -luaL_checknumber(L, 2);
    mLatReal *c = qlua_newLatReal(L, lua_gettop(L));

    CALL_QDP(L);
    QDP_R_eq_r(c->ptr, &b, *S->qss);
    QDP_R_peq_R(c->ptr, a->ptr, *S->qss);

    return 1;

}

static int
q_R_mul_R(lua_State *L)
{
    mLatReal *a = qlua_checkLatReal(L, 1, NULL);
    mLattice *S = qlua_ObjLattice(L, 1);
    mLatReal *b = qlua_checkLatReal(L, 2, S);
    mLatReal *c = qlua_newLatReal(L, lua_gettop(L));

    CALL_QDP(L);
    QDP_R_eq_R_times_R(c->ptr, a->ptr, b->ptr, *S->qss);

    return 1;
}

static int
q_r_mul_R(lua_State *L)
{
    QLA_Real b = luaL_checknumber(L, 1);
    mLatReal *a = qlua_checkLatReal(L, 2, NULL);
    mLattice *S = qlua_ObjLattice(L, 2);
    mLatReal *c = qlua_newLatReal(L, lua_gettop(L));

    CALL_QDP(L);
    QDP_R_eq_r_times_R(c->ptr, &b, a->ptr, *S->qss);

    return 1;
}

static int
q_R_mul_r(lua_State *L)
{
    mLatReal *a = qlua_checkLatReal(L, 1, NULL);
    mLattice *S = qlua_ObjLattice(L, 1);
    QLA_Real b = luaL_checknumber(L, 2);
    mLatReal *c = qlua_newLatReal(L, lua_gettop(L));

    CALL_QDP(L);
    QDP_R_eq_r_times_R(c->ptr, &b, a->ptr, *S->qss);

    return 1;
}

static int
q_R_div_R(lua_State *L)
{
    mLatReal *a = qlua_checkLatReal(L, 1, NULL);
    mLattice *S = qlua_ObjLattice(L, 1);
    mLatReal *b = qlua_checkLatReal(L, 2, S);
    mLatReal *c = qlua_newLatReal(L, lua_gettop(L));

    CALL_QDP(L);
    QDP_R_eq_R_divide_R(c->ptr, a->ptr, b->ptr, *S->qss);

    return 1;
}

static int
q_R_div_r(lua_State *L)
{
    mLatReal *a = qlua_checkLatReal(L, 1, NULL);
    mLattice *S = qlua_ObjLattice(L, 1);
    QLA_Real b = 1/luaL_checknumber(L, 2);
    mLatReal *c = qlua_newLatReal(L, lua_gettop(L));

    CALL_QDP(L);
    QDP_R_eq_r_times_R(c->ptr, &b, a->ptr, *S->qss);

    return 1;
}

static int
q_r_div_R(lua_State *L)
{
    QLA_D_Real a = luaL_checknumber(L, 1);
    mLatReal *b = qlua_checkLatReal(L, 2, NULL);
    mLattice *S = qlua_ObjLattice(L, 2);
    int Sidx = lua_gettop(L);
    mLatReal *z = qlua_newLatReal(L, Sidx);
    mLatReal *c = qlua_newLatReal(L, Sidx);

    CALL_QDP(L);
    QDP_R_eq_r(z->ptr, &a, *S->qss);
    QDP_R_eq_R_divide_R(c->ptr, z->ptr, b->ptr, *S->qss);

    return 1;
}

static int
q_r_min_R(lua_State *L)
{
    QLA_D_Real a = luaL_checknumber(L, 1);
    mLatReal *b = qlua_checkLatReal(L, 2, NULL);
    mLattice *S = qlua_ObjLattice(L, 2);
    int Sidx = lua_gettop(L);
    mLatReal *c = qlua_newLatReal(L, Sidx);
    mLatReal *r = qlua_newLatReal(L, Sidx);
    
    CALL_QDP(L);
    QDP_R_eq_r(c->ptr, &a, *S->qss);
    QDP_R_eq_R_min_R(r->ptr, c->ptr, b->ptr, *S->qss);

    return 1;
}

static int
q_R_min_r(lua_State *L)
{
    mLatReal *a = qlua_checkLatReal(L, 1, NULL);
    mLattice *S = qlua_ObjLattice(L, 1);
    int Sidx = lua_gettop(L);
    QLA_D_Real b = luaL_checknumber(L, 2);
    mLatReal *c = qlua_newLatReal(L, Sidx);
    mLatReal *r = qlua_newLatReal(L, Sidx);
    
    CALL_QDP(L);
    QDP_R_eq_r(c->ptr, &b, *S->qss);
    QDP_R_eq_R_min_R(r->ptr, a->ptr, c->ptr, *S->qss);

    return 1;
}

static int
q_R_min_R(lua_State *L)
{
    mLatReal *a = qlua_checkLatReal(L, 1, NULL);
    mLattice *S = qlua_ObjLattice(L, 1);
    mLatReal *b = qlua_checkLatReal(L, 2, S);
    mLatReal *r = qlua_newLatReal(L, lua_gettop(L));
    
    CALL_QDP(L);
    QDP_R_eq_R_min_R(r->ptr, a->ptr, b->ptr, *S->qss);

    return 1;
}

static int
q_r_max_R(lua_State *L)
{
    QLA_D_Real a = luaL_checknumber(L, 1);
    mLatReal *b = qlua_checkLatReal(L, 2, NULL);
    mLattice *S = qlua_ObjLattice(L, 2);
    int Sidx = lua_gettop(L);
    mLatReal *c = qlua_newLatReal(L, Sidx);
    mLatReal *r = qlua_newLatReal(L, Sidx);
    
    CALL_QDP(L);
    QDP_R_eq_r(c->ptr, &a, *S->qss);
    QDP_R_eq_R_max_R(r->ptr, c->ptr, b->ptr, *S->qss);

    return 1;
}

static int
q_R_max_r(lua_State *L)
{
    mLatReal *a = qlua_checkLatReal(L, 1, NULL);
    mLattice *S = qlua_ObjLattice(L, 1);
    int Sidx = lua_gettop(L);
    QLA_D_Real b = luaL_checknumber(L, 2);
    mLatReal *c = qlua_newLatReal(L, Sidx);
    mLatReal *r = qlua_newLatReal(L, Sidx);
    
    CALL_QDP(L);
    QDP_R_eq_r(c->ptr, &b, *S->qss);
    QDP_R_eq_R_max_R(r->ptr, a->ptr, c->ptr, *S->qss);

    return 1;
}

static int
q_R_max_R(lua_State *L)
{
    mLatReal *a = qlua_checkLatReal(L, 1, NULL);
    mLattice *S = qlua_ObjLattice(L, 1);
    mLatReal *b = qlua_checkLatReal(L, 2, S);
    mLatReal *r = qlua_newLatReal(L, lua_gettop(L));
    
    CALL_QDP(L);
    QDP_R_eq_R_max_R(r->ptr, a->ptr, b->ptr, *S->qss);

    return 1;
}

static int
q_r_eq_R(lua_State *L)
{
    QLA_D_Real a = luaL_checknumber(L, 1);
    mLatReal *b = qlua_checkLatReal(L, 2, NULL);
    mLattice *S = qlua_ObjLattice(L, 2);
    int Sidx = lua_gettop(L);
    mLatReal *c = qlua_newLatReal(L, Sidx);
    mLatInt *r = qlua_newLatInt(L, Sidx);
    
    CALL_QDP(L);
    QDP_R_eq_r(c->ptr, &a, *S->qss);
    QDP_I_eq_R_eq_R(r->ptr, c->ptr, b->ptr, *S->qss);

    return 1;
}

static int
q_R_eq_r(lua_State *L)
{
    mLatReal *a = qlua_checkLatReal(L, 1, NULL);
    mLattice *S = qlua_ObjLattice(L, 1);
    int Sidx = lua_gettop(L);
    QLA_D_Real b = luaL_checknumber(L, 2);
    mLatReal *c = qlua_newLatReal(L, Sidx);
    mLatInt *r = qlua_newLatInt(L, Sidx);
    
    CALL_QDP(L);
    QDP_R_eq_r(c->ptr, &b, *S->qss);
    QDP_I_eq_R_eq_R(r->ptr, a->ptr, c->ptr, *S->qss);

    return 1;
}

static int
q_R_eq_R(lua_State *L)
{
    mLatReal *a = qlua_checkLatReal(L, 1, NULL);
    mLattice *S = qlua_ObjLattice(L, 1);
    mLatReal *b = qlua_checkLatReal(L, 2, S);
    mLatInt *r = qlua_newLatInt(L, lua_gettop(L));
    
    CALL_QDP(L);
    QDP_I_eq_R_eq_R(r->ptr, a->ptr, b->ptr, *S->qss);

    return 1;
}

static int
q_r_ne_R(lua_State *L)
{
    QLA_D_Real a = luaL_checknumber(L, 1);
    mLatReal *b = qlua_checkLatReal(L, 2, NULL);
    mLattice *S = qlua_ObjLattice(L, 2);
    int Sidx = lua_gettop(L);
    mLatReal *c = qlua_newLatReal(L, Sidx);
    mLatInt *r = qlua_newLatInt(L, Sidx);
    
    CALL_QDP(L);
    QDP_R_eq_r(c->ptr, &a, *S->qss);
    QDP_I_eq_R_ne_R(r->ptr, c->ptr, b->ptr, *S->qss);

    return 1;
}

static int
q_R_ne_r(lua_State *L)
{
    mLatReal *a = qlua_checkLatReal(L, 1, NULL);
    mLattice *S = qlua_ObjLattice(L, 1);
    int Sidx = lua_gettop(L);
    QLA_D_Real b = luaL_checknumber(L, 2);
    mLatReal *c = qlua_newLatReal(L, Sidx);
    mLatInt *r = qlua_newLatInt(L, Sidx);
    
    CALL_QDP(L);
    QDP_R_eq_r(c->ptr, &b, *S->qss);
    QDP_I_eq_R_ne_R(r->ptr, a->ptr, c->ptr, *S->qss);

    return 1;
}

static int
q_R_ne_R(lua_State *L)
{
    mLatReal *a = qlua_checkLatReal(L, 1, NULL);
    mLattice *S = qlua_ObjLattice(L, 1);
    mLatReal *b = qlua_checkLatReal(L, 2, S);
    mLatInt *r = qlua_newLatInt(L, lua_gettop(L));
    
    CALL_QDP(L);
    QDP_I_eq_R_ne_R(r->ptr, a->ptr, b->ptr, *S->qss);

    return 1;
}

static int
q_r_lt_R(lua_State *L)
{
    QLA_D_Real a = luaL_checknumber(L, 1);
    mLatReal *b = qlua_checkLatReal(L, 2, NULL);
    mLattice *S = qlua_ObjLattice(L, 2);
    int Sidx = lua_gettop(L);
    mLatReal *c = qlua_newLatReal(L, Sidx);
    mLatInt *r = qlua_newLatInt(L, Sidx);
    
    CALL_QDP(L);
    QDP_R_eq_r(c->ptr, &a, *S->qss);
    QDP_I_eq_R_lt_R(r->ptr, c->ptr, b->ptr, *S->qss);

    return 1;
}

static int
q_R_lt_r(lua_State *L)
{
    mLatReal *a = qlua_checkLatReal(L, 1, NULL);
    mLattice *S = qlua_ObjLattice(L, 1);
    int Sidx = lua_gettop(L);
    QLA_D_Real b = luaL_checknumber(L, 2);
    mLatReal *c = qlua_newLatReal(L, Sidx);
    mLatInt *r = qlua_newLatInt(L, Sidx);
    
    CALL_QDP(L);
    QDP_R_eq_r(c->ptr, &b, *S->qss);
    QDP_I_eq_R_lt_R(r->ptr, a->ptr, c->ptr, *S->qss);

    return 1;
}

static int
q_R_lt_R(lua_State *L)
{
    mLatReal *a = qlua_checkLatReal(L, 1, NULL);
    mLattice *S = qlua_ObjLattice(L, 1);
    mLatReal *b = qlua_checkLatReal(L, 2, S);
    mLatInt *r = qlua_newLatInt(L, lua_gettop(L));
    
    CALL_QDP(L);
    QDP_I_eq_R_lt_R(r->ptr, a->ptr, b->ptr, *S->qss);

    return 1;
}

static int
q_r_le_R(lua_State *L)
{
    QLA_D_Real a = luaL_checknumber(L, 1);
    mLatReal *b = qlua_checkLatReal(L, 2, NULL);
    mLattice *S = qlua_ObjLattice(L, 2);
    int Sidx = lua_gettop(L);
    mLatReal *c = qlua_newLatReal(L, Sidx);
    mLatInt *r = qlua_newLatInt(L, Sidx);
    
    CALL_QDP(L);
    QDP_R_eq_r(c->ptr, &a, *S->qss);
    QDP_I_eq_R_le_R(r->ptr, c->ptr, b->ptr, *S->qss);

    return 1;
}

static int
q_R_le_r(lua_State *L)
{
    mLatReal *a = qlua_checkLatReal(L, 1, NULL);
    mLattice *S = qlua_ObjLattice(L, 1);
    int Sidx = lua_gettop(L);
    QLA_D_Real b = luaL_checknumber(L, 2);
    mLatReal *c = qlua_newLatReal(L, Sidx);
    mLatInt *r = qlua_newLatInt(L, Sidx);
    
    CALL_QDP(L);
    QDP_R_eq_r(c->ptr, &b, *S->qss);
    QDP_I_eq_R_le_R(r->ptr, a->ptr, c->ptr, *S->qss);

    return 1;
}

static int
q_R_le_R(lua_State *L)
{
    mLatReal *a = qlua_checkLatReal(L, 1, NULL);
    mLattice *S = qlua_ObjLattice(L, 1);
    mLatReal *b = qlua_checkLatReal(L, 2, S);
    mLatInt *r = qlua_newLatInt(L, lua_gettop(L));
    
    CALL_QDP(L);
    QDP_I_eq_R_le_R(r->ptr, a->ptr, b->ptr, *S->qss);

    return 1;
}

static int
q_r_gt_R(lua_State *L)
{
    QLA_D_Real a = luaL_checknumber(L, 1);
    mLatReal *b = qlua_checkLatReal(L, 2, NULL);
    mLattice *S = qlua_ObjLattice(L, 2);
    int Sidx = lua_gettop(L);
    mLatReal *c = qlua_newLatReal(L, Sidx);
    mLatInt *r = qlua_newLatInt(L, Sidx);
    
    CALL_QDP(L);
    QDP_R_eq_r(c->ptr, &a, *S->qss);
    QDP_I_eq_R_gt_R(r->ptr, c->ptr, b->ptr, *S->qss);

    return 1;
}

static int
q_R_gt_r(lua_State *L)
{
    mLatReal *a = qlua_checkLatReal(L, 1, NULL);
    mLattice *S = qlua_ObjLattice(L, 1);
    int Sidx = lua_gettop(L);
    QLA_D_Real b = luaL_checknumber(L, 2);
    mLatReal *c = qlua_newLatReal(L, Sidx);
    mLatInt *r = qlua_newLatInt(L, Sidx);
    
    CALL_QDP(L);
    QDP_R_eq_r(c->ptr, &b, *S->qss);
    QDP_I_eq_R_gt_R(r->ptr, a->ptr, c->ptr, *S->qss);

    return 1;
}

static int
q_R_gt_R(lua_State *L)
{
    mLatReal *a = qlua_checkLatReal(L, 1, NULL);
    mLattice *S = qlua_ObjLattice(L, 1);
    mLatReal *b = qlua_checkLatReal(L, 2, S);
    mLatInt *r = qlua_newLatInt(L, lua_gettop(L));
    
    CALL_QDP(L);
    QDP_I_eq_R_gt_R(r->ptr, a->ptr, b->ptr, *S->qss);

    return 1;
}

static int
q_r_ge_R(lua_State *L)
{
    QLA_D_Real a = luaL_checknumber(L, 1);
    mLatReal *b = qlua_checkLatReal(L, 2, NULL);
    mLattice *S = qlua_ObjLattice(L, 2);
    int Sidx = lua_gettop(L);
    mLatReal *c = qlua_newLatReal(L, Sidx);
    mLatInt *r = qlua_newLatInt(L, Sidx);
    
    CALL_QDP(L);
    QDP_R_eq_r(c->ptr, &a, *S->qss);
    QDP_I_eq_R_ge_R(r->ptr, c->ptr, b->ptr, *S->qss);

    return 1;
}

static int
q_R_ge_r(lua_State *L)
{
    mLatReal *a = qlua_checkLatReal(L, 1, NULL);
    mLattice *S = qlua_ObjLattice(L, 1);
    int Sidx = lua_gettop(L);
    QLA_D_Real b = luaL_checknumber(L, 2);
    mLatReal *c = qlua_newLatReal(L, Sidx);
    mLatInt *r = qlua_newLatInt(L, Sidx);
    
    CALL_QDP(L);
    QDP_R_eq_r(c->ptr, &b, *S->qss);
    QDP_I_eq_R_ge_R(r->ptr, a->ptr, c->ptr, *S->qss);

    return 1;
}

static int
q_R_ge_R(lua_State *L)
{
    mLatReal *a = qlua_checkLatReal(L, 1, NULL);
    mLattice *S = qlua_ObjLattice(L, 1);
    mLatReal *b = qlua_checkLatReal(L, 2, S);
    mLatInt *r = qlua_newLatInt(L, lua_gettop(L));
    
    CALL_QDP(L);
    QDP_I_eq_R_ge_R(r->ptr, a->ptr, b->ptr, *S->qss);

    return 1;
}



static int
q_latreal(lua_State *L)
{
    mLattice *S = qlua_checkLattice(L, 1);

    switch (lua_gettop(L)) {
    case 1: {
        qlua_newZeroLatReal(L, 1);
        return 1;
    }
    case 2:
        switch (qlua_qtype(L, 2)) {
        case qReal: {
            QLA_Real d = luaL_checknumber(L, 2);
            mLatReal *v = qlua_newLatReal(L, 1);
            
            CALL_QDP(L);
            QDP_R_eq_r(v->ptr, &d, *S->qss);
            
            return 1;
        }
        case qLatInt: {
            mLatInt *d = qlua_checkLatInt(L, 2, S);
            mLatReal *v = qlua_newLatReal(L, 1);
            
            CALL_QDP(L);
            QDP_R_eq_I(v->ptr, d->ptr, *S->qss);
            
            return 1;
        }
        case qLatReal: { /* no mixing constructor */
            mLatReal *d = qlua_checkLatReal(L, 2, S);
            mLatReal *v = qlua_newLatReal(L, 1);
            
            CALL_QDP(L);
            QDP_R_eq_R(v->ptr, d->ptr, *S->qss);
            
            return 1;
        }
        default:
            break;
        }
    }
    return qlua_badconstr(L, "Real");
}

static struct luaL_Reg mtLatReal[] = {
    { "__tostring",   q_R_fmt      },
    { "__gc",         q_R_gc       },
    { "__index",      q_R_get      },
    { "__newindex",   q_R_put      },
    { "__unm",        q_R_neg      },
    { "__add",        qlua_add     },
    { "__sub",        qlua_sub     },
    { "__mul",        qlua_mul     },
    { "__div",        qlua_div     },
    { "sum",          q_R_sum      },
    { "norm2",        q_R_norm2    },
    { "shift",        q_R_shift    },
    { "sin",          q_R_sin      },
    { "cos",          q_R_cos      },
    { "tan",          q_R_tan      },
    { "asin",         q_R_asin     },
    { "acos",         q_R_acos     },
    { "atan",         q_R_atan     },
    { "sqrt",         q_R_sqrt     },
    { "abs",          q_R_abs      },
    { "exp",          q_R_exp      },
    { "log",          q_R_log      },
    { "sign",         q_R_sign     },
    { "ceil",         q_R_ceil     },
    { "floor",        q_R_floor    },
    { "sinh",         q_R_sinh     },
    { "cosh",         q_R_cosh     },
    { "tanh",         q_R_tanh     },
    { "log10",        q_R_log10    },
    { "expi",         q_R_expi     },
    { "trunc",        q_R_trunc    },
    { "round",        q_R_round    },
    { "set",          q_R_set      },
    { "min",          q_R_min      },
    { "max",          q_R_max      },
    /* "lattice" */
    /* "a-type"  */
    { NULL,           NULL         }
};

mLatReal *
qlua_newLatReal(lua_State *L, int Sidx)
{
    mLattice *S = qlua_checkLattice(L, Sidx);
    QDP_Real *v = QDP_create_R_L(S->lat);
    mLatReal *hdr;

    if (v == 0) {
        lua_gc(L, LUA_GCCOLLECT, 0);
        v = QDP_D_create_R_L(S->lat);
        if (v == 0)
            luaL_error(L, "not enough memory (QDP_Real)");
    }
    hdr = lua_newuserdata(L, sizeof (mLatReal));
    hdr->ptr = v;
    qlua_createLatticeTable(L, Sidx, mtLatReal, qLatReal, LatRealName);
    lua_setmetatable(L, -2);
    qlua_qdp_memuse(L, "Real", 1);

    return hdr;
}

mLatReal *
qlua_newZeroLatReal(lua_State *L, int Sidx)
{
        mLatReal *v = qlua_newLatReal(L, Sidx);
        mLattice *S = qlua_checkLattice(L, Sidx);
        QDP_R_eq_zero(v->ptr, S->all);
        return v;
}

mLatReal *
qlua_checkLatReal(lua_State *L, int idx, mLattice *S)
{
    void *v = qlua_checkLatticeType(L, idx, qLatReal, LatRealName);
    
    if (S) {
        mLattice *S1 = qlua_ObjLattice(L, idx);
        if (S1->id != S->id)
            luaL_error(L, "%s on a wrong lattice", LatRealName);
        lua_pop(L, 1);
    }

    return (mLatReal *)v;
}

int
q_R_random(lua_State *L)
{
    mLatRandom *a = qlua_checkLatRandom(L, 1, NULL);
    mLattice *S = qlua_ObjLattice(L, 1);
    mLatReal *r = qlua_newLatReal(L, lua_gettop(L));

    CALL_QDP(L);
    QDP_R_eq_random_S(r->ptr, a->ptr, *S->qss);

    return 1;
}

int
q_R_gaussian(lua_State *L)
{
    mLatRandom *a = qlua_checkLatRandom(L, 1, NULL);
    mLattice *S = qlua_ObjLattice(L, 1);
    mLatReal *r = qlua_newLatReal(L, lua_gettop(L));

    CALL_QDP(L);
    QDP_R_eq_gaussian_S(r->ptr, a->ptr, *S->qss);

    return 1;
}

static struct luaL_Reg fLatReal[] = {
    { "Real",       q_latreal },
    { NULL,         NULL }
};

int
init_latreal(lua_State *L)
{
    static const QLUA_Op2 ops[] = {
        {qlua_add_table, qLatReal,  qLatReal, q_R_add_R  },
        {qlua_add_table, qReal,     qLatReal, q_r_add_R  },
        {qlua_add_table, qLatReal,  qReal,    q_R_add_r  },
        {qlua_sub_table, qLatReal,  qLatReal, q_R_sub_R  },
        {qlua_sub_table, qReal,     qLatReal, q_r_sub_R  },
        {qlua_sub_table, qLatReal,  qReal,    q_R_sub_r  },
        {qlua_mul_table, qLatReal,  qLatReal, q_R_mul_R  },
        {qlua_mul_table, qReal,     qLatReal, q_r_mul_R  },
        {qlua_mul_table, qLatReal,  qReal,    q_R_mul_r  },
        {qlua_div_table, qLatReal,  qLatReal, q_R_div_R  },
        {qlua_div_table, qReal,     qLatReal, q_r_div_R  },
        {qlua_div_table, qLatReal,  qReal,    q_R_div_r  },
        {NULL,           qNoType,   qNoType,  NULL        }
    };
    static const QLUA_ZOp2 zops[] = {
        {qlua_min_table,  zReal,    zLatReal, q_r_min_R },
        {qlua_min_table,  zLatReal, zReal,    q_R_min_r },
        {qlua_min_table,  zLatReal, zLatReal, q_R_min_R },
        {qlua_max_table,  zReal,    zLatReal, q_r_max_R },
        {qlua_max_table,  zLatReal, zReal,    q_R_max_r },
        {qlua_max_table,  zLatReal, zLatReal, q_R_max_R },
        {qlua_eq_table,   zReal,    zLatReal, q_r_eq_R  },
        {qlua_eq_table,   zLatReal, zReal,    q_R_eq_r  },
        {qlua_eq_table,   zLatReal, zLatReal, q_R_eq_R  },
        {qlua_ne_table,   zReal,    zLatReal, q_r_ne_R  },
        {qlua_ne_table,   zLatReal, zReal,    q_R_ne_r  },
        {qlua_ne_table,   zLatReal, zLatReal, q_R_ne_R  },
        {qlua_lt_table,   zReal,    zLatReal, q_r_lt_R  },
        {qlua_lt_table,   zLatReal, zReal,    q_R_lt_r  },
        {qlua_lt_table,   zLatReal, zLatReal, q_R_lt_R  },
        {qlua_le_table,   zReal,    zLatReal, q_r_le_R  },
        {qlua_le_table,   zLatReal, zReal,    q_R_le_r  },
        {qlua_le_table,   zLatReal, zLatReal, q_R_le_R  },
        {qlua_gt_table,   zReal,    zLatReal, q_r_gt_R  },
        {qlua_gt_table,   zLatReal, zReal,    q_R_gt_r  },
        {qlua_gt_table,   zLatReal, zLatReal, q_R_gt_R  },
        {qlua_ge_table,   zReal,    zLatReal, q_r_ge_R  },
        {qlua_ge_table,   zLatReal, zReal,    q_R_ge_r  },
        {qlua_ge_table,   zLatReal, zLatReal, q_R_ge_R  },
        {NULL,            zNoType, zNoType, NULL      }
    };

    luaL_getmetatable(L, opLattice);
    luaL_register(L, NULL, fLatReal);
    lua_pop(L, 1);
    qlua_reg_op2(ops);
    qlua_reg_zop2(zops);
    qlua_reg_dot(qLatReal, q_R_mul_R);

    return 0;
}

void
fini_latreal(void)
{
}
