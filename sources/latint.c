#include "qlua.h"                                                    /* DEPS */
#include "qvector.h"                                                 /* DEPS */
#include "lattice.h"                                                 /* DEPS */
#include "latsubset.h"                                               /* DEPS */
#include "latint.h"                                                  /* DEPS */
#include "latmulti.h"                                                /* DEPS */
#include <stdlib.h>  
#include "qmp.h"
#include <string.h>

/* lattice integers */
static const char LatIntName[] = "lattice.Int";

static int
q_I_fmt(lua_State *L)
{
    char fmt[72];
    mLatInt *b = qlua_checkLatInt(L, 1, NULL);

    sprintf(fmt, "QDP:Integer(%p)", b->ptr);
    lua_pushstring(L, fmt);

    return 1;
}

static int
q_I_gc(lua_State *L)
{
    mLatInt *b = qlua_checkLatInt(L, 1, NULL);

    QDP_destroy_I(b->ptr);
    b->ptr = 0;
    qlua_qdp_memuse(L, "Integer", -1);

    return 0;
}

static int
q_I_sum(lua_State *L)
{
    mLatInt *a = qlua_checkLatInt(L, 1, NULL);
    int argc = lua_gettop(L);
    mLattice *S = qlua_ObjLattice(L, 1);
    int Sidx = lua_gettop(L);

    switch (argc) {
    case 1: {
        QLA_Real sum;

        CALL_QDP(L);
        if (S->lss.mask) {
            mLatInt *b = qlua_newZeroLatInt(L, Sidx);
            QDP_I_eq_I_mask_I(b->ptr, a->ptr, S->lss.mask, *S->qss);
            QDP_r_eq_sum_I(&sum, b->ptr, *S->qss);
            lua_pop(L, 1);
        } else {
            QDP_r_eq_sum_I(&sum, a->ptr, *S->qss);
        }
        lua_pushnumber(L, sum);
        
        return 1;
    }
    case 2: { /* NB: does not use the subset mask */
        mLatMulti *m = qlua_checkLatMulti(L, 2, S);
        int size = m->size;
        QLA_Int *ii = m->idx;
        mVecReal *r = qlua_newVecReal(L, size);
        int sites = QDP_sites_on_node_L(S->lat);
        int k;
        QLA_Int *xx;
        
        for (k = 0; k < size; k++)
            r->val[k] = 0;

        CALL_QDP(L);
        xx = QDP_expose_I(a->ptr);
        for (k = 0; k < sites; k++, xx++, ii++) {
            int t = *ii;
            if ((t < 0) || (t >= size))
                continue;
            r->val[t] += *xx;
        }
        QDP_reset_I(a->ptr);
        QMP_sum_double_array(r->val, size);

        return 1;
    }
    }
    return luaL_error(L, "bad arguments for Int:sum()");
}

static int
q_I_set(lua_State *L)
{
    mLatInt *r = qlua_checkLatInt(L, 1, NULL);
    mLattice *S = qlua_ObjLattice(L, 1);
    mLatInt *a = qlua_checkLatInt(L, 2, S);

    CALL_QDP(L);
    if (S->lss.mask) {
        QDP_I_eq_I_mask_I(r->ptr, a->ptr, S->lss.mask, *S->qss);
    } else {
        QDP_I_eq_I(r->ptr, a->ptr, *S->qss);
    }
    lua_pop(L, 1);

    return 1;
}

static int
q_I_norm2(lua_State *L)
{
    mLatInt *a = qlua_checkLatInt(L, 1, NULL);
    mLattice *S = qlua_ObjLattice(L, 1);
    QLA_Real n;

    CALL_QDP(L);
    if (S->lss.mask) {
        mLatInt *b = qlua_newLatInt(L, lua_gettop(L));
        QDP_I_eq_I_mask_I(b->ptr, a->ptr, S->lss.mask, *S->qss);
        QDP_r_eq_norm2_I(&n, b->ptr, *S->qss);
    } else {
        QDP_r_eq_norm2_I(&n, a->ptr, *S->qss);
    }
    lua_pushnumber(L, n);

    return 1;
}

static int
q_I_xor(lua_State *L)
{
    mLatInt *a = qlua_checkLatInt(L, 1, NULL);
    mLattice *S = qlua_ObjLattice(L, 1);
    mLatInt *b = qlua_checkLatInt(L, 2, S);
    mLatInt *r = qlua_newLatInt(L, lua_gettop(L));

    CALL_QDP(L);
    QDP_I_eq_I_xor_I(r->ptr, a->ptr, b->ptr, *S->qss);
    return 1;
}

static int
q_I_or(lua_State *L)
{
    mLatInt *a = qlua_checkLatInt(L, 1, NULL);
    mLattice *S = qlua_ObjLattice(L, 1);
    mLatInt *b = qlua_checkLatInt(L, 2, S);
    mLatInt *r = qlua_newLatInt(L, lua_gettop(L));

    CALL_QDP(L);
    QDP_I_eq_I_or_I(r->ptr, a->ptr, b->ptr, *S->qss);
    return 1;
}

static int
q_I_and(lua_State *L)
{
    mLatInt *a = qlua_checkLatInt(L, 1, NULL);
    mLattice *S = qlua_ObjLattice(L, 1);
    mLatInt *b = qlua_checkLatInt(L, 2, S);
    mLatInt *r = qlua_newLatInt(L, lua_gettop(L));

    CALL_QDP(L);
    QDP_I_eq_I_and_I(r->ptr, a->ptr, b->ptr, *S->qss);
    return 1;
}

static int
q_I_min(lua_State *L)
{
    mLatInt *a = qlua_checkLatInt(L, 1, NULL);
    mLattice *S = qlua_ObjLattice(L, 1);
    mLatInt *b = qlua_checkLatInt(L, 2, S);
    mLatInt *r = qlua_newLatInt(L, lua_gettop(L));

    CALL_QDP(L);
    QDP_I_eq_I_min_I(r->ptr, a->ptr, b->ptr, *S->qss);
    return 1;
}

static int
q_I_max(lua_State *L)
{
    mLatInt *a = qlua_checkLatInt(L, 1, NULL);
    mLattice *S = qlua_ObjLattice(L, 1);
    mLatInt *b = qlua_checkLatInt(L, 2, S);
    mLatInt *r = qlua_newLatInt(L, lua_gettop(L));

    CALL_QDP(L);
    QDP_I_eq_I_max_I(r->ptr, a->ptr, b->ptr, *S->qss);
    return 1;
}

static int
q_I_shift(lua_State *L)
{
    mLatInt *a = qlua_checkLatInt(L, 1, NULL);
    mLattice *S = qlua_ObjLattice(L, 1);
    QDP_Shift shift = qlua_checkShift(L, 2, S);
    QDP_ShiftDir dir = qlua_checkShiftDir(L, 3);
    mLatInt *b = qlua_newLatInt(L, lua_gettop(L));

    CALL_QDP(L);
    QDP_I_eq_sI(b->ptr, a->ptr, shift, dir, *S->qss);

    return 1;
}

static int
q_I_get(lua_State *L)
{
    switch (qlua_qtype(L, 2)) {
    case qTable: {
        mLatInt *V = qlua_checkLatInt(L, 1, NULL);
        mLattice *S = qlua_ObjLattice(L, 1);
        QLA_Int *locked;
        int *idx = 0;
        int z;

        idx = qlua_checklatcoord(L, 2, S);
        CALL_QDP(L);
        locked = QDP_expose_I(V->ptr);
        int site_node = QDP_node_number_L(S->lat, idx);
        if (site_node == QDP_this_node) {
            z = QLA_elem_I(locked[QDP_index_L(S->lat, idx)]);
        }
        QDP_reset_I(V->ptr);
        qlua_free(L, idx);
        XMP_dist_int_array(site_node, 1, &z);
        lua_pushnumber(L, z);

        return 1;
    }
    case qString:
        return qlua_selflookup(L, 1, luaL_checkstring(L, 2));
    default:
        break;
    }

    return qlua_badindex(L, "Int");
}

static int
q_I_put(lua_State *L)
{
    mLatInt *V = qlua_checkLatInt(L, 1, NULL);
    mLattice *S = qlua_ObjLattice(L, 1);
    QLA_Int *locked;
    int *idx = 0;
    int z = luaL_checkint(L, 3);

    idx = qlua_checklatcoord(L, 2, S);
    CALL_QDP(L);
    locked = QDP_expose_I(V->ptr);
    if (QDP_node_number_L(S->lat, idx) == QDP_this_node) {
        QLA_elem_I(locked[QDP_index_L(S->lat, idx)]) = z;
    }
    QDP_reset_I(V->ptr);
    qlua_free(L, idx);

    return 0;
}

static int
q_I_add_I(lua_State *L)
{
    mLatInt *a = qlua_checkLatInt(L, 1, NULL);
    mLattice *S = qlua_ObjLattice(L, 1);
    mLatInt *b = qlua_checkLatInt(L, 2, S);
    mLatInt *r = qlua_newLatInt(L, lua_gettop(L));

    CALL_QDP(L);
    QDP_I_eq_I_plus_I(r->ptr, a->ptr, b->ptr, *S->qss);

    return 1;
}

static int
q_i_add_I(lua_State *L)
{
    QLA_Int x = luaL_checkint(L, 1);
    mLatInt *a = qlua_checkLatInt(L, 2, NULL);
    mLattice *S = qlua_ObjLattice(L, 2);
    mLatInt *r = qlua_newLatInt(L, lua_gettop(L));

    CALL_QDP(L);
    QDP_I_eq_i(r->ptr, &x, *S->qss);
    QDP_I_peq_I(r->ptr, a->ptr, *S->qss);

    return 1;
}

static int
q_I_add_i(lua_State *L)
{
    mLatInt *a = qlua_checkLatInt(L, 1, NULL);
    mLattice *S = qlua_ObjLattice(L, 1);
    QLA_Int x = luaL_checkint(L, 2);
    mLatInt *r = qlua_newLatInt(L, lua_gettop(L));

    CALL_QDP(L);
    QDP_I_eq_i(r->ptr, &x, *S->qss);
    QDP_I_peq_I(r->ptr, a->ptr, *S->qss);

    return 1;
}

static int
q_I_sub_I(lua_State *L)
{
    mLatInt *a = qlua_checkLatInt(L, 1, NULL);
    mLattice *S = qlua_ObjLattice(L, 1);
    mLatInt *b = qlua_checkLatInt(L, 2, S);
    mLatInt *r = qlua_newLatInt(L, lua_gettop(L));

    CALL_QDP(L);
    QDP_I_eq_I_minus_I(r->ptr, a->ptr, b->ptr, *S->qss);

    return 1;
}

static int
q_i_sub_I(lua_State *L)
{
    QLA_Int x = luaL_checkint(L, 1);
    mLatInt *a = qlua_checkLatInt(L, 2, NULL);
    mLattice *S = qlua_ObjLattice(L, 2);
    mLatInt *r = qlua_newLatInt(L, lua_gettop(L));

    CALL_QDP(L);
    QDP_I_eq_i(r->ptr, &x, *S->qss);
    QDP_I_meq_I(r->ptr, a->ptr, *S->qss);

    return 1;
}

static int
q_I_sub_i(lua_State *L)
{
    mLatInt *a = qlua_checkLatInt(L, 1, NULL);
    mLattice *S = qlua_ObjLattice(L, 1);
    QLA_Int x = -luaL_checkint(L, 2);
    mLatInt *r = qlua_newLatInt(L, lua_gettop(L));

    CALL_QDP(L);
    QDP_I_eq_i(r->ptr, &x, *S->qss);
    QDP_I_peq_I(r->ptr, a->ptr, *S->qss);

    return 1;
}

static int
q_i_mul_I(lua_State *L)
{
    QLA_Int a = luaL_checkint(L, 1);
    mLatInt *b = qlua_checkLatInt(L, 2, NULL);
    mLattice *S = qlua_ObjLattice(L, 2);
    mLatInt *r = qlua_newLatInt(L, lua_gettop(L));

    CALL_QDP(L);
    QDP_I_eq_i_times_I(r->ptr, &a, b->ptr, *S->qss);

    return 1;
}

static int
q_I_mul_i(lua_State *L)
{
    mLatInt *b = qlua_checkLatInt(L, 1, NULL);
    mLattice *S = qlua_ObjLattice(L, 1);
    QLA_Int a = luaL_checkint(L, 2);
    mLatInt *r = qlua_newLatInt(L, lua_gettop(L));

    CALL_QDP(L);
    QDP_I_eq_i_times_I(r->ptr, &a, b->ptr, *S->qss);

    return 1;
}

static int
q_I_mul_I(lua_State *L)
{
    mLatInt *a = qlua_checkLatInt(L, 1, NULL);
    mLattice *S = qlua_ObjLattice(L, 1);
    mLatInt *b = qlua_checkLatInt(L, 2, S);
    mLatInt *r = qlua_newLatInt(L, lua_gettop(L));

    CALL_QDP(L);
    QDP_I_eq_I_times_I(r->ptr, a->ptr, b->ptr, *S->qss);

    return 1;
}

static int
q_I_div_I(lua_State *L)
{
    mLatInt *a = qlua_checkLatInt(L, 1, NULL);
    mLattice *S = qlua_ObjLattice(L, 1);
    mLatInt *b = qlua_checkLatInt(L, 2, S);
    mLatInt *r = qlua_newLatInt(L, lua_gettop(L));

    CALL_QDP(L);
    QDP_I_eq_I_divide_I(r->ptr, a->ptr, b->ptr, *S->qss);

    return 1;
}

static int
q_i_div_I(lua_State *L)
{
    QLA_Int x = luaL_checkint(L, 1);
    mLatInt *a = qlua_checkLatInt(L, 2, NULL);
    mLattice *S = qlua_ObjLattice(L, 2);
    int Sidx = lua_gettop(L);
    mLatInt *z = qlua_newLatInt(L, Sidx);
    mLatInt *r = qlua_newLatInt(L, Sidx);

    CALL_QDP(L);
    QDP_I_eq_i(z->ptr, &x, *S->qss);
    QDP_I_eq_I_divide_I(r->ptr, z->ptr, a->ptr, *S->qss);

    return 1;
}

static int
q_I_div_i(lua_State *L)
{
    mLatInt *a = qlua_checkLatInt(L, 1, NULL);
    mLattice *S = qlua_ObjLattice(L, 1);
    int Sidx = lua_gettop(L);
    QLA_Int x = luaL_checkint(L, 2);
    mLatInt *z = qlua_newLatInt(L, Sidx);
    mLatInt *r = qlua_newLatInt(L, Sidx);

    CALL_QDP(L);
    QDP_I_eq_i(z->ptr, &x, *S->qss);
    QDP_I_eq_I_divide_I(r->ptr, a->ptr, z->ptr, *S->qss);

    return 1;
}

static int
q_I_mod_I(lua_State *L)
{
    mLatInt *a = qlua_checkLatInt(L, 1, NULL);
    mLattice *S = qlua_ObjLattice(L, 1);
    mLatInt *b = qlua_checkLatInt(L, 2, S);
    mLatInt *r = qlua_newLatInt(L, lua_gettop(L));

    CALL_QDP(L);
    QDP_I_eq_I_mod_I(r->ptr, a->ptr, b->ptr, *S->qss);

    return 1;
}

static int
q_i_mod_I(lua_State *L)
{
    QLA_Int x = luaL_checkint(L, 1);
    mLatInt *a = qlua_checkLatInt(L, 2, NULL);
    mLattice *S = qlua_ObjLattice(L, 2);
    int Sidx = lua_gettop(L);
    mLatInt *z = qlua_newLatInt(L, Sidx);
    mLatInt *r = qlua_newLatInt(L, Sidx);

    CALL_QDP(L);
    QDP_I_eq_i(z->ptr, &x, *S->qss);
    QDP_I_eq_I_mod_I(r->ptr, z->ptr, a->ptr, *S->qss);

    return 1;
}

static int
q_I_mod_i(lua_State *L)
{
    mLatInt *a = qlua_checkLatInt(L, 1, NULL);
    mLattice *S = qlua_ObjLattice(L, 1);
    int Sidx = lua_gettop(L);
    int x = luaL_checkint(L, 2);
    mLatInt *z = qlua_newLatInt(L, Sidx);
    mLatInt *r = qlua_newLatInt(L, Sidx);

    CALL_QDP(L);
    QDP_I_eq_i(z->ptr, &x, *S->qss);
    QDP_I_eq_I_mod_I(r->ptr, a->ptr, z->ptr, *S->qss);

    return 1;
}

static int
q_I_neg(lua_State *L)
{
    mLatInt *a = qlua_checkLatInt(L, 1, NULL);
    mLattice *S = qlua_ObjLattice(L, 1);
    mLatInt *r = qlua_newLatInt(L, lua_gettop(L));
    QLA_Int m1 = -1;

    CALL_QDP(L);
    QDP_I_eq_i_times_I(r->ptr, &m1, a->ptr, *S->qss);

    return 1;
}

static int
q_i_min_I(lua_State *L)
{
    QLA_Int a = luaL_checkint(L, 1);
    mLatInt *b = qlua_checkLatInt(L, 2, NULL);
    mLattice *S = qlua_ObjLattice(L, 2);
    int Sidx = lua_gettop(L);
    mLatInt *c = qlua_newLatInt(L, Sidx);
    mLatInt *r = qlua_newLatInt(L, Sidx);
    
    CALL_QDP(L);
    QDP_I_eq_i(c->ptr, &a, *S->qss);
    QDP_I_eq_I_min_I(r->ptr, c->ptr, b->ptr, *S->qss);

    return 1;
}

static int
q_I_min_i(lua_State *L)
{
    mLatInt *a = qlua_checkLatInt(L, 1, NULL);
    mLattice *S = qlua_ObjLattice(L, 1);
    int Sidx = lua_gettop(L);
    QLA_Int b = luaL_checkint(L, 2);
    mLatInt *c = qlua_newLatInt(L, Sidx);
    mLatInt *r = qlua_newLatInt(L, Sidx);
    
    CALL_QDP(L);
    QDP_I_eq_i(c->ptr, &b, *S->qss);
    QDP_I_eq_I_min_I(r->ptr, a->ptr, c->ptr, *S->qss);

    return 1;
}

static int
q_I_min_I(lua_State *L)
{
    mLatInt *a = qlua_checkLatInt(L, 1, NULL);
    mLattice *S = qlua_ObjLattice(L, 1);
    mLatInt *b = qlua_checkLatInt(L, 2, S);
    mLatInt *r = qlua_newLatInt(L, lua_gettop(L));
    
    CALL_QDP(L);
    QDP_I_eq_I_min_I(r->ptr, a->ptr, b->ptr, *S->qss);

    return 1;
}

static int
q_i_max_I(lua_State *L)
{
    QLA_Int a = luaL_checkint(L, 1);
    mLatInt *b = qlua_checkLatInt(L, 2, NULL);
    mLattice *S = qlua_ObjLattice(L, 2);
    int Sidx = lua_gettop(L);
    mLatInt *c = qlua_newLatInt(L, Sidx);
    mLatInt *r = qlua_newLatInt(L, Sidx);
    
    CALL_QDP(L);
    QDP_I_eq_i(c->ptr, &a, *S->qss);
    QDP_I_eq_I_max_I(r->ptr, c->ptr, b->ptr, *S->qss);

    return 1;
}

static int
q_I_max_i(lua_State *L)
{
    mLatInt *a = qlua_checkLatInt(L, 1, NULL);
    mLattice *S = qlua_ObjLattice(L, 1);
    int Sidx = lua_gettop(L);
    QLA_Int b = luaL_checkint(L, 2);
    mLatInt *c = qlua_newLatInt(L, Sidx);
    mLatInt *r = qlua_newLatInt(L, Sidx);
    
    CALL_QDP(L);
    QDP_I_eq_i(c->ptr, &b, *S->qss);
    QDP_I_eq_I_max_I(r->ptr, a->ptr, c->ptr, *S->qss);

    return 1;
}

static int
q_I_max_I(lua_State *L)
{
    mLatInt *a = qlua_checkLatInt(L, 1, NULL);
    mLattice *S = qlua_ObjLattice(L, 1);
    mLatInt *b = qlua_checkLatInt(L, 2, S);
    mLatInt *r = qlua_newLatInt(L, lua_gettop(L));
    
    CALL_QDP(L);
    QDP_I_eq_I_max_I(r->ptr, a->ptr, b->ptr, *S->qss);

    return 1;
}

static int
q_i_eq_I(lua_State *L)
{
    QLA_Int a = luaL_checkint(L, 1);
    mLatInt *b = qlua_checkLatInt(L, 2, NULL);
    mLattice *S = qlua_ObjLattice(L, 2);
    int Sidx = lua_gettop(L);
    mLatInt *c = qlua_newLatInt(L, Sidx);
    mLatInt *r = qlua_newLatInt(L, Sidx);
    
    CALL_QDP(L);
    QDP_I_eq_i(c->ptr, &a, *S->qss);
    QDP_I_eq_I_eq_I(r->ptr, c->ptr, b->ptr, *S->qss);

    return 1;
}

static int
q_I_eq_i(lua_State *L)
{
    mLatInt *a = qlua_checkLatInt(L, 1, NULL);
    mLattice *S = qlua_ObjLattice(L, 1);
    int Sidx = lua_gettop(L);
    QLA_Int b = luaL_checkint(L, 2);
    mLatInt *c = qlua_newLatInt(L, Sidx);
    mLatInt *r = qlua_newLatInt(L, Sidx);
    
    CALL_QDP(L);
    QDP_I_eq_i(c->ptr, &b, *S->qss);
    QDP_I_eq_I_eq_I(r->ptr, a->ptr, c->ptr, *S->qss);

    return 1;
}

static int
q_I_eq_I(lua_State *L)
{
    mLatInt *a = qlua_checkLatInt(L, 1, NULL);
    mLattice *S = qlua_ObjLattice(L, 1);
    mLatInt *b = qlua_checkLatInt(L, 2, S);
    mLatInt *r = qlua_newLatInt(L, lua_gettop(L));
    
    CALL_QDP(L);
    QDP_I_eq_I_eq_I(r->ptr, a->ptr, b->ptr, *S->qss);

    return 1;
}

static int
q_i_ne_I(lua_State *L)
{
    QLA_Int a = luaL_checkint(L, 1);
    mLatInt *b = qlua_checkLatInt(L, 2, NULL);
    mLattice *S = qlua_ObjLattice(L, 2);
    int Sidx = lua_gettop(L);
    mLatInt *c = qlua_newLatInt(L, Sidx);
    mLatInt *r = qlua_newLatInt(L, Sidx);
    
    CALL_QDP(L);
    QDP_I_eq_i(c->ptr, &a, *S->qss);
    QDP_I_eq_I_ne_I(r->ptr, c->ptr, b->ptr, *S->qss);

    return 1;
}

static int
q_I_ne_i(lua_State *L)
{
    mLatInt *a = qlua_checkLatInt(L, 1, NULL);
    mLattice *S = qlua_ObjLattice(L, 1);
    int Sidx = lua_gettop(L);
    QLA_Int b = luaL_checkint(L, 2);
    mLatInt *c = qlua_newLatInt(L, Sidx);
    mLatInt *r = qlua_newLatInt(L, Sidx);
    
    CALL_QDP(L);
    QDP_I_eq_i(c->ptr, &b, *S->qss);
    QDP_I_eq_I_ne_I(r->ptr, a->ptr, c->ptr, *S->qss);

    return 1;
}

static int
q_I_ne_I(lua_State *L)
{
    mLatInt *a = qlua_checkLatInt(L, 1, NULL);
    mLattice *S = qlua_ObjLattice(L, 1);
    mLatInt *b = qlua_checkLatInt(L, 2, S);
    mLatInt *r = qlua_newLatInt(L, lua_gettop(L));
    
    CALL_QDP(L);
    QDP_I_eq_I_ne_I(r->ptr, a->ptr, b->ptr, *S->qss);

    return 1;
}

static int
q_i_lt_I(lua_State *L)
{
    QLA_Int a = luaL_checkint(L, 1);
    mLatInt *b = qlua_checkLatInt(L, 2, NULL);
    mLattice *S = qlua_ObjLattice(L, 2);
    int Sidx = lua_gettop(L);
    mLatInt *c = qlua_newLatInt(L, Sidx);
    mLatInt *r = qlua_newLatInt(L, Sidx);
    
    CALL_QDP(L);
    QDP_I_eq_i(c->ptr, &a, *S->qss);
    QDP_I_eq_I_lt_I(r->ptr, c->ptr, b->ptr, *S->qss);

    return 1;
}

static int
q_I_lt_i(lua_State *L)
{
    mLatInt *a = qlua_checkLatInt(L, 1, NULL);
    mLattice *S = qlua_ObjLattice(L, 1);
    int Sidx = lua_gettop(L);
    QLA_Int b = luaL_checkint(L, 2);
    mLatInt *c = qlua_newLatInt(L, Sidx);
    mLatInt *r = qlua_newLatInt(L, Sidx);
    
    CALL_QDP(L);
    QDP_I_eq_i(c->ptr, &b, *S->qss);
    QDP_I_eq_I_lt_I(r->ptr, a->ptr, c->ptr, *S->qss);

    return 1;
}

static int
q_I_lt_I(lua_State *L)
{
    mLatInt *a = qlua_checkLatInt(L, 1, NULL);
    mLattice *S = qlua_ObjLattice(L, 1);
    mLatInt *b = qlua_checkLatInt(L, 2, S);
    mLatInt *r = qlua_newLatInt(L, lua_gettop(L));
    
    CALL_QDP(L);
    QDP_I_eq_I_lt_I(r->ptr, a->ptr, b->ptr, *S->qss);

    return 1;
}

static int
q_i_le_I(lua_State *L)
{
    QLA_Int a = luaL_checkint(L, 1);
    mLatInt *b = qlua_checkLatInt(L, 2, NULL);
    mLattice *S = qlua_ObjLattice(L, 2);
    int Sidx = lua_gettop(L);
    mLatInt *c = qlua_newLatInt(L, Sidx);
    mLatInt *r = qlua_newLatInt(L, Sidx);
    
    CALL_QDP(L);
    QDP_I_eq_i(c->ptr, &a, *S->qss);
    QDP_I_eq_I_le_I(r->ptr, c->ptr, b->ptr, *S->qss);

    return 1;
}

static int
q_I_le_i(lua_State *L)
{
    mLatInt *a = qlua_checkLatInt(L, 1, NULL);
    mLattice *S = qlua_ObjLattice(L, 1);
    int Sidx = lua_gettop(L);
    QLA_Int b = luaL_checkint(L, 2);
    mLatInt *c = qlua_newLatInt(L, Sidx);
    mLatInt *r = qlua_newLatInt(L, Sidx);
    
    CALL_QDP(L);
    QDP_I_eq_i(c->ptr, &b, *S->qss);
    QDP_I_eq_I_le_I(r->ptr, a->ptr, c->ptr, *S->qss);

    return 1;
}

static int
q_I_le_I(lua_State *L)
{
    mLatInt *a = qlua_checkLatInt(L, 1, NULL);
    mLattice *S = qlua_ObjLattice(L, 1);
    mLatInt *b = qlua_checkLatInt(L, 2, S);
    mLatInt *r = qlua_newLatInt(L, lua_gettop(L));
    
    CALL_QDP(L);
    QDP_I_eq_I_le_I(r->ptr, a->ptr, b->ptr, *S->qss);

    return 1;
}

static int
q_i_gt_I(lua_State *L)
{
    QLA_Int a = luaL_checkint(L, 1);
    mLatInt *b = qlua_checkLatInt(L, 2, NULL);
    mLattice *S = qlua_ObjLattice(L, 2);
    int Sidx = lua_gettop(L);
    mLatInt *c = qlua_newLatInt(L, Sidx);
    mLatInt *r = qlua_newLatInt(L, Sidx);
    
    CALL_QDP(L);
    QDP_I_eq_i(c->ptr, &a, *S->qss);
    QDP_I_eq_I_gt_I(r->ptr, c->ptr, b->ptr, *S->qss);

    return 1;
}

static int
q_I_gt_i(lua_State *L)
{
    mLatInt *a = qlua_checkLatInt(L, 1, NULL);
    mLattice *S = qlua_ObjLattice(L, 1);
    int Sidx = lua_gettop(L);
    QLA_Int b = luaL_checkint(L, 2);
    mLatInt *c = qlua_newLatInt(L, Sidx);
    mLatInt *r = qlua_newLatInt(L, Sidx);
    
    CALL_QDP(L);
    QDP_I_eq_i(c->ptr, &b, *S->qss);
    QDP_I_eq_I_gt_I(r->ptr, a->ptr, c->ptr, *S->qss);

    return 1;
}

static int
q_I_gt_I(lua_State *L)
{
    mLatInt *a = qlua_checkLatInt(L, 1, NULL);
    mLattice *S = qlua_ObjLattice(L, 1);
    mLatInt *b = qlua_checkLatInt(L, 2, S);
    mLatInt *r = qlua_newLatInt(L, lua_gettop(L));
    
    CALL_QDP(L);
    QDP_I_eq_I_gt_I(r->ptr, a->ptr, b->ptr, *S->qss);

    return 1;
}

static int
q_i_ge_I(lua_State *L)
{
    QLA_Int a = luaL_checkint(L, 1);
    mLatInt *b = qlua_checkLatInt(L, 2, NULL);
    mLattice *S = qlua_ObjLattice(L, 2);
    int Sidx = lua_gettop(L);
    mLatInt *c = qlua_newLatInt(L, Sidx);
    mLatInt *r = qlua_newLatInt(L, Sidx);
    
    CALL_QDP(L);
    QDP_I_eq_i(c->ptr, &a, *S->qss);
    QDP_I_eq_I_ge_I(r->ptr, c->ptr, b->ptr, *S->qss);

    return 1;
}

static int
q_I_ge_i(lua_State *L)
{
    mLatInt *a = qlua_checkLatInt(L, 1, NULL);
    mLattice *S = qlua_ObjLattice(L, 1);
    int Sidx = lua_gettop(L);
    QLA_Int b = luaL_checkint(L, 2);
    mLatInt *c = qlua_newLatInt(L, Sidx);
    mLatInt *r = qlua_newLatInt(L, Sidx);
    
    CALL_QDP(L);
    QDP_I_eq_i(c->ptr, &b, *S->qss);
    QDP_I_eq_I_ge_I(r->ptr, a->ptr, c->ptr, *S->qss);

    return 1;
}

static int
q_I_ge_I(lua_State *L)
{
    mLatInt *a = qlua_checkLatInt(L, 1, NULL);
    mLattice *S = qlua_ObjLattice(L, 1);
    mLatInt *b = qlua_checkLatInt(L, 2, S);
    mLatInt *r = qlua_newLatInt(L, lua_gettop(L));
    
    CALL_QDP(L);
    QDP_I_eq_I_ge_I(r->ptr, a->ptr, b->ptr, *S->qss);

    return 1;
}

static int
q_latint(lua_State *L)
{
    mLattice *S = qlua_checkLattice(L, 1);

    switch (lua_gettop(L)) {
    case 1:{
        qlua_newZeroLatInt(L, 1);
        return 1;
    }
    case 2:
        switch (qlua_qtype(L, 2)) {
        case qReal: {
            QLA_Int d = luaL_checkint(L, 2);
            mLatInt *v = qlua_newLatInt(L, 1);
            
            CALL_QDP(L);
            QDP_I_eq_i(v->ptr, &d, *S->qss);
            
            return 1;
        }
        case qLatInt: {
            mLatInt *a = qlua_checkLatInt(L, 2, S);
            mLatInt *res = qlua_newLatInt(L, 1);
            
            CALL_QDP(L);
            QDP_I_eq_I(res->ptr, a->ptr, *S->qss);
            
            return 1;
        }
        default:
            break;
        }
    }

    return qlua_badconstr(L, "Int");
}

static struct luaL_Reg mtLatInt[] = {
    { "__tostring",   q_I_fmt },
    { "__gc",         q_I_gc },
    { "__index",      q_I_get },
    { "__newindex",   q_I_put },
    { "__unm",        q_I_neg },
    { "__add",        qlua_add },
    { "__sub",        qlua_sub },
    { "__mul",        qlua_mul },
    { "__div",        qlua_div },
    { "__mod",        qlua_mod },
    { "norm2",        q_I_norm2 },
    { "shift",        q_I_shift },
    { "sum",          q_I_sum },
    { "set",          q_I_set },
    { "xor",          q_I_xor },
    { "or",           q_I_or },
    { "and",          q_I_and },
    { "max",          q_I_max },
    { "min",          q_I_min },
    /* "lattice" is inserted when the table is creates */
    /* "a-type"  is inserted as well */
    { NULL,           NULL}
};

mLatInt *
qlua_newLatInt(lua_State *L, int Sidx)
{
    mLattice *S = qlua_checkLattice(L, Sidx);
    QDP_Int *v = QDP_create_I_L(S->lat);
    mLatInt *hdr;
    
    if (v == 0) {
        lua_gc(L, LUA_GCCOLLECT, 0);
        v = QDP_create_I_L(S->lat);
        if (v == 0)
            luaL_error(L, "not enough memory (QDP_Int)");
    }
    hdr = lua_newuserdata(L, sizeof (mLatInt));
    hdr->ptr = v;
    qlua_createLatticeTable(L, Sidx, mtLatInt, qLatInt, LatIntName);
    lua_setmetatable(L, -2);
    qlua_qdp_memuse(L, "Integer", 1);

    return hdr;
}

mLatInt *
qlua_newZeroLatInt(lua_State *L, int Sidx)
{
        mLatInt *v = qlua_newLatInt(L, Sidx);
        mLattice *S = qlua_checkLattice(L, Sidx);
        
        QDP_I_eq_zero(v->ptr, S->all);
        return v;
}

mLatInt *
qlua_checkLatInt(lua_State *L, int idx, mLattice *S)
{
    void *v = qlua_checkLatticeType(L, idx, qLatInt, LatIntName);
    
    if (S) {
        mLattice *S1 = qlua_ObjLattice(L, idx);
        if (S1->id != S->id)
            luaL_error(L, "%s on a wrong lattice", LatIntName);
        lua_pop(L, 1);
    }

    return (mLatInt *)v;
}

mLatInt *
qlua_tabkey_LatInt(lua_State *L, int idx, const char *key, mLattice *S)
{
  mLatInt *v;

  if (!qlua_tabpushopt_key(L, idx, key))
    luaL_error(L, "expecting LatInt in { %s = ...}", key);
  v = qlua_checkLatInt(L, -1, S);
  lua_pop(L, 1); /* expect the user not to drop the object */

  return v;
}

mLatInt *
qlua_tabidx_LatInt(lua_State *L, int idx, int subidx, mLattice *S)
{
  mLatInt *v;

  if (!qlua_tabpushopt_idx(L, idx, subidx))
    luaL_error(L, "expecting LatInt in { [%d] = ...}", subidx);
  v = qlua_checkLatInt(L, -1, S);
  lua_pop(L, 1); /* expect the user not to drop the object */

  return v;
}

static struct luaL_Reg fLatInt[] = {
    { "Int",     q_latint },
    { NULL, NULL}
};

int
init_latint(lua_State *L)
{
    static const QLUA_Op2 ops[] = {
        {qlua_add_table, qLatInt, qLatInt, q_I_add_I},
        {qlua_add_table, qReal,   qLatInt, q_i_add_I},
        {qlua_add_table, qLatInt, qReal,   q_I_add_i},
        {qlua_sub_table, qLatInt, qLatInt, q_I_sub_I},
        {qlua_sub_table, qReal,   qLatInt, q_i_sub_I},
        {qlua_sub_table, qLatInt, qReal,   q_I_sub_i},
        {qlua_mul_table, qLatInt, qLatInt, q_I_mul_I},
        {qlua_mul_table, qReal,   qLatInt, q_i_mul_I},
        {qlua_mul_table, qLatInt, qReal,   q_I_mul_i},
        {qlua_div_table, qLatInt, qLatInt, q_I_div_I},
        {qlua_div_table, qReal,   qLatInt, q_i_div_I},
        {qlua_div_table, qLatInt, qReal,   q_I_div_i},
        {qlua_mod_table, qLatInt, qLatInt, q_I_mod_I},
        {qlua_mod_table, qReal,   qLatInt, q_i_mod_I},
        {qlua_mod_table, qLatInt, qReal,   q_I_mod_i},
        {NULL,           qOther,  qOther,  NULL     }
    };
    static const QLUA_ZOp2 zops[] = {
        {qlua_min_table,  zReal,   zLatInt, q_i_min_I },
        {qlua_min_table,  zLatInt, zReal,   q_I_min_i },
        {qlua_min_table,  zLatInt, zLatInt, q_I_min_I },
        {qlua_max_table,  zReal,   zLatInt, q_i_max_I },
        {qlua_max_table,  zLatInt, zReal,   q_I_max_i },
        {qlua_max_table,  zLatInt, zLatInt, q_I_max_I },
        {qlua_eq_table,   zReal,   zLatInt, q_i_eq_I  },
        {qlua_eq_table,   zLatInt, zReal,   q_I_eq_i  },
        {qlua_eq_table,   zLatInt, zLatInt, q_I_eq_I  },
        {qlua_ne_table,   zReal,   zLatInt, q_i_ne_I  },
        {qlua_ne_table,   zLatInt, zReal,   q_I_ne_i  },
        {qlua_ne_table,   zLatInt, zLatInt, q_I_ne_I  },
        {qlua_lt_table,   zReal,   zLatInt, q_i_lt_I  },
        {qlua_lt_table,   zLatInt, zReal,   q_I_lt_i  },
        {qlua_lt_table,   zLatInt, zLatInt, q_I_lt_I  },
        {qlua_le_table,   zReal,   zLatInt, q_i_le_I  },
        {qlua_le_table,   zLatInt, zReal,   q_I_le_i  },
        {qlua_le_table,   zLatInt, zLatInt, q_I_le_I  },
        {qlua_gt_table,   zReal,   zLatInt, q_i_gt_I  },
        {qlua_gt_table,   zLatInt, zReal,   q_I_gt_i  },
        {qlua_gt_table,   zLatInt, zLatInt, q_I_gt_I  },
        {qlua_ge_table,   zReal,   zLatInt, q_i_ge_I  },
        {qlua_ge_table,   zLatInt, zReal,   q_I_ge_i  },
        {qlua_ge_table,   zLatInt, zLatInt, q_I_ge_I  },
        {NULL,            zNoType, zNoType, NULL      }
    };

    luaL_getmetatable(L, opLattice);
    luaL_register(L, NULL, fLatInt);
    lua_pop(L, 1);
    qlua_reg_op2(ops);
    qlua_reg_zop2(zops);

    return 0;
}

int
fini_latint(lua_State *L)
{
    return 0;
}
