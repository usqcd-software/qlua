#include <qlua.h>                                                    /* DEPS */
#include <qvector.h>                                                 /* DEPS */
#include <lattice.h>                                                 /* DEPS */
#include <latint.h>                                                  /* DEPS */
#include <latmulti.h>                                                /* DEPS */
#include <stdlib.h>  
#include <qmp.h>
#include <string.h>

/* lattice integers */
const char mtnLatInt[] = "lattice.Int";
static const char opLatInt[] = "lattice.Int.ops";

mLatInt *
qlua_newLatInt(lua_State *L)
{
    QDP_Int *v = QDP_create_I();
    mLatInt *hdr;
    
    if (v == 0) {
        lua_gc(L, LUA_GCCOLLECT, 0);
        v = QDP_create_I();
        if (v == 0)
            luaL_error(L, "not enough memory (QDP_Int)");
    }
    hdr = lua_newuserdata(L, sizeof (mLatInt));
    hdr->ptr = v;
    luaL_getmetatable(L, mtnLatInt);
    lua_setmetatable(L, -2);

    return hdr;
}

mLatInt *
qlua_checkLatInt(lua_State *L, int idx)
{
    void *v = luaL_checkudata(L, idx, mtnLatInt);

    luaL_argcheck(L, v != 0, idx, "lattice.Int expected");
    
    return v;
}

static int
q_I_fmt(lua_State *L)
{
    char fmt[72];
    mLatInt *b = qlua_checkLatInt(L, 1);

    sprintf(fmt, "QDP:Integer(%p)", b->ptr);
    lua_pushstring(L, fmt);

    return 1;
}

static int
q_I_gc(lua_State *L)
{
    mLatInt *b = qlua_checkLatInt(L, 1);

    QDP_destroy_I(b->ptr);
    b->ptr = 0;

    return 0;
}

static int
q_I_sum(lua_State *L)
{
    mLatInt *a = qlua_checkLatInt(L, 1);

    switch (lua_gettop(L)) {
    case 1: {
        QLA_Real sum;

        CALL_QDP(L);
        QDP_r_eq_sum_I(&sum, a->ptr, *qCurrent);
        lua_pushnumber(L, sum);
        
        return 1;
    }
    case 2: {
        mLatMulti *m = qlua_checkLatMulti(L, 2);
        int size = m->size;
        QLA_Int *ii = m->idx;
        mVecReal *r = qlua_newVecReal(L, size);
        int k;
        QLA_Int *xx;
        
        r->size = size;
        for (k = 0; k < size; k++)
            r->val[k] = 0;

        CALL_QDP(L);
        xx = QDP_expose_I(a->ptr);
        for (k = 0; k < QDP_sites_on_node; k++, xx++, ii++) {
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
    mLatInt *r = qlua_checkLatInt(L, 1);
    mLatInt *a = qlua_checkLatInt(L, 2);

    CALL_QDP(L);
    QDP_I_eq_I(r->ptr, a->ptr, *qCurrent);
    lua_pop(L, 1);

    return 1;
}

static int
q_I_norm2(lua_State *L)
{
    mLatInt *a = qlua_checkLatInt(L, 1);
    QLA_Real sum;

    CALL_QDP(L);
    QDP_r_eq_norm2_I(&sum, a->ptr, *qCurrent);
    lua_pushnumber(L, sum);

    return 1;
}

static int
q_I_shift(lua_State *L)
{
    mLatInt *a = qlua_checkLatInt(L, 1);
    QDP_Shift shift = qlua_checkShift(L, 2);
    QDP_ShiftDir dir = qlua_checkShiftDir(L, 3);
    mLatInt *b = qlua_newLatInt(L);

    CALL_QDP(L);
    QDP_I_eq_sI(b->ptr, a->ptr, shift, dir, *qCurrent);

    return 1;
}

static int
q_I_get(lua_State *L)
{
    switch (qlua_gettype(L, 2)) {
    case qTable: {
        mLatInt *V = qlua_checkLatInt(L, 1);
        QLA_Int *locked;
        int *idx = 0;
        int z;

        idx = qlua_checklatcoord(L, 2);
        CALL_QDP(L);
        locked = QDP_expose_I(V->ptr);
        if (QDP_node_number(idx) == QDP_this_node) {
            z = QLA_elem_I(locked[QDP_index(idx)]);
        } else {
            z = 0;
        }
        QDP_reset_I(V->ptr);
        qlua_free(L, idx);
        QMP_sum_int(&z);
        lua_pushnumber(L, z);

        return 1;
    }
    case qString:
        return qlua_lookup(L, 2, opLatInt);
    }

    return qlua_badindex(L, "Int");
}

static int
q_I_put(lua_State *L)
{
    mLatInt *V = qlua_checkLatInt(L, 1);
    QLA_Int *locked;
    int *idx = 0;
    int z = luaL_checkint(L, 3);

    idx = qlua_checklatcoord(L, 2);
    CALL_QDP(L);
    locked = QDP_expose_I(V->ptr);
    if (QDP_node_number(idx) == QDP_this_node) {
        QLA_elem_I(locked[QDP_index(idx)]) = z;
    }
    QDP_reset_I(V->ptr);
    qlua_free(L, idx);

    return 0;
}

static struct {
    QLA_Int *a;
    QLA_Int *b;
    int x;
} Iargs; /* YYY global state */

static int
q_I_add_I(lua_State *L)
{
    mLatInt *a = qlua_checkLatInt(L, 1);
    mLatInt *b = qlua_checkLatInt(L, 2);
    mLatInt *r = qlua_newLatInt(L);

    CALL_QDP(L);
    QDP_I_eq_I_plus_I(r->ptr, a->ptr, b->ptr, *qCurrent);

    return 1;
}

static void
do_add_iI(QLA_Int *r, int idx)
{
    *r = Iargs.x + Iargs.a[idx];
}

static int
q_i_add_I(lua_State *L)
{
    int x = luaL_checkint(L, 1);
    mLatInt *a = qlua_checkLatInt(L, 2);
    mLatInt *r = qlua_newLatInt(L);

    CALL_QDP(L);
    Iargs.x = x;
    Iargs.a = QDP_expose_I(a->ptr);
    QDP_I_eq_funci(r->ptr, do_add_iI, *qCurrent);
    Iargs.a = 0;
    QDP_reset_I(a->ptr);

    return 1;
}

static int
q_I_add_i(lua_State *L)
{
    mLatInt *a = qlua_checkLatInt(L, 1);
    int x = luaL_checkint(L, 2);
    mLatInt *r = qlua_newLatInt(L);

    CALL_QDP(L);
    Iargs.x = x;
    Iargs.a = QDP_expose_I(a->ptr);
    QDP_I_eq_funci(r->ptr, do_add_iI, *qCurrent);
    Iargs.a = 0;
    QDP_reset_I(a->ptr);

    return 1;
}

static int
q_I_sub_I(lua_State *L)
{
    mLatInt *a = qlua_checkLatInt(L, 1);
    mLatInt *b = qlua_checkLatInt(L, 2);
    mLatInt *r = qlua_newLatInt(L);

    CALL_QDP(L);
    QDP_I_eq_I_minus_I(r->ptr, a->ptr, b->ptr, *qCurrent);

    return 1;
}

static void
do_sub_iI(QLA_Int *r, int idx)
{
    *r = Iargs.x - Iargs.a[idx];
}

static void
do_sub_Ii(QLA_Int *r, int idx)
{
    *r = Iargs.a[idx] - Iargs.x;
}

static int
q_i_sub_I(lua_State *L)
{
    int x = luaL_checkint(L, 1);
    mLatInt *a = qlua_checkLatInt(L, 2);
    mLatInt *r = qlua_newLatInt(L);

    CALL_QDP(L);
    Iargs.x = x;
    Iargs.a = QDP_expose_I(a->ptr);
    QDP_I_eq_funci(r->ptr, do_sub_iI, *qCurrent);
    Iargs.a = 0;
    QDP_reset_I(a->ptr);

    return 1;
}

static int
q_I_sub_i(lua_State *L)
{
    mLatInt *a = qlua_checkLatInt(L, 1);
    int x = luaL_checkint(L, 2);
    mLatInt *r = qlua_newLatInt(L);

    CALL_QDP(L);
    Iargs.x = x;
    Iargs.a = QDP_expose_I(a->ptr);
    QDP_I_eq_funci(r->ptr, do_sub_Ii, *qCurrent);
    Iargs.a = 0;
    QDP_reset_I(a->ptr);

    return 1;
}

static int
q_i_mul_I(lua_State *L)
{
    QLA_Int a = luaL_checkint(L, 1);
    mLatInt *b = qlua_checkLatInt(L, 2);
    mLatInt *r = qlua_newLatInt(L);

    CALL_QDP(L);
    QDP_I_eq_i_times_I(r->ptr, &a, b->ptr, *qCurrent);

    return 1;
}

static int
q_I_mul_i(lua_State *L)
{
    mLatInt *b = qlua_checkLatInt(L, 1);
    QLA_Int a = luaL_checkint(L, 2);
    mLatInt *r = qlua_newLatInt(L);

    CALL_QDP(L);
    QDP_I_eq_i_times_I(r->ptr, &a, b->ptr, *qCurrent);

    return 1;
}

static int
q_I_mul_I(lua_State *L)
{
    mLatInt *a = qlua_checkLatInt(L, 1);
    mLatInt *b = qlua_checkLatInt(L, 2);
    mLatInt *r = qlua_newLatInt(L);

    CALL_QDP(L);
    QDP_I_eq_I_times_I(r->ptr, a->ptr, b->ptr, *qCurrent);

    return 1;
}

static int
q_I_div_I(lua_State *L)
{
    mLatInt *a = qlua_checkLatInt(L, 1);
    mLatInt *b = qlua_checkLatInt(L, 2);
    mLatInt *r = qlua_newLatInt(L);

    CALL_QDP(L);
    QDP_I_eq_I_divide_I(r->ptr, a->ptr, b->ptr, *qCurrent);

    return 1;
}

static void
do_div_iI(QLA_Int *r, int idx)
{
    *r = Iargs.x / Iargs.a[idx];
}

static void
do_div_Ii(QLA_Int *r, int idx)
{
    *r = Iargs.a[idx] / Iargs.x;
}

static int
q_i_div_I(lua_State *L)
{
    int x = luaL_checkint(L, 1);
    mLatInt *a = qlua_checkLatInt(L, 2);
    mLatInt *r = qlua_newLatInt(L);

    CALL_QDP(L);
    Iargs.x = x;
    Iargs.a = QDP_expose_I(a->ptr);
    QDP_I_eq_funci(r->ptr, do_div_iI, *qCurrent);
    Iargs.a = 0;
    QDP_reset_I(a->ptr);

    return 1;
}

static int
q_I_div_i(lua_State *L)
{
    mLatInt *a = qlua_checkLatInt(L, 1);
    int x = luaL_checkint(L, 2);
    mLatInt *r = qlua_newLatInt(L);

    CALL_QDP(L);
    Iargs.x = x;
    Iargs.a = QDP_expose_I(a->ptr);
    QDP_I_eq_funci(r->ptr, do_div_Ii, *qCurrent);
    Iargs.a = 0;
    QDP_reset_I(a->ptr);

    return 1;
}

static void
do_mod_II(QLA_Int *r, int idx)
{
    *r = Iargs.a[idx] % Iargs.b[idx];
}

static void
do_mod_iI(QLA_Int *r, int idx)
{
    *r = Iargs.x % Iargs.a[idx];
}

static void
do_mod_Ii(QLA_Int *r, int idx)
{
    *r = Iargs.a[idx] % Iargs.x;
}

static int
q_I_mod_I(lua_State *L)
{
    mLatInt *a = qlua_checkLatInt(L, 1);
    mLatInt *b = qlua_checkLatInt(L, 2);
    mLatInt *r = qlua_newLatInt(L);

    CALL_QDP(L);
    Iargs.a = QDP_expose_I(a->ptr);
    Iargs.b = QDP_expose_I(b->ptr);
    QDP_I_eq_funci(r->ptr, do_mod_II, *qCurrent);
    Iargs.b = 0;
    Iargs.a = 0;
    QDP_reset_I(a->ptr);
    QDP_reset_I(b->ptr);

    return 1;
}

static int
q_i_mod_I(lua_State *L)
{
    int x = luaL_checkint(L, 1);
    mLatInt *a = qlua_checkLatInt(L, 2);
    mLatInt *r = qlua_newLatInt(L);

    CALL_QDP(L);
    Iargs.x = x;
    Iargs.a = QDP_expose_I(a->ptr);
    QDP_I_eq_funci(r->ptr, do_mod_iI, *qCurrent);
    Iargs.a = 0;
    QDP_reset_I(a->ptr);

    return 1;
}

static int
q_I_mod_i(lua_State *L)
{
    mLatInt *a = qlua_checkLatInt(L, 1);
    int x = luaL_checkint(L, 2);
    mLatInt *r = qlua_newLatInt(L);

    CALL_QDP(L);
    Iargs.x = x;
    Iargs.a = QDP_expose_I(a->ptr);
    QDP_I_eq_funci(r->ptr, do_mod_Ii, *qCurrent);
    Iargs.a = 0;
    QDP_reset_I(a->ptr);

    return 1;
}

static int
q_I_neg(lua_State *L)
{
    mLatInt *a = qlua_checkLatInt(L, 1);
    mLatInt *r = qlua_newLatInt(L);
    QLA_Int m1 = -1;

    CALL_QDP(L);
    QDP_I_eq_i_times_I(r->ptr, &m1, a->ptr, *qCurrent);

    return 1;
}

static int
q_latint(lua_State *L)
{
    switch (lua_gettop(L)) {
    case 1:{
        mLatInt *v = qlua_newLatInt(L);

        CALL_QDP(L);
        QDP_I_eq_zero(v->ptr, *qCurrent);

        return 1;
    }
    case 2:
        switch (qlua_gettype(L, 2)) {
        case qReal: {
            QLA_Int d = luaL_checkint(L, 2);
            mLatInt *v = qlua_newLatInt(L);
            
            CALL_QDP(L);
            QDP_I_eq_i(v->ptr, &d, *qCurrent);
            
            return 1;
        }
        case qLatInt: {
            mLatInt *res = qlua_newLatInt(L);
            mLatInt *a = qlua_checkLatInt(L, 2);
            
            CALL_QDP(L);
            QDP_I_eq_I(res->ptr, a->ptr, *qCurrent);
            
            return 1;
        }
        }
    }

    return qlua_badconstr(L, "Int");
}

static struct luaL_Reg LatIntMethods[] = {
    { "norm2",  q_I_norm2 },
    { "shift",  q_I_shift },
    { "sum",    q_I_sum },
    { "set",    q_I_set },
    { NULL,     NULL}
};

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
    { NULL,           NULL}
};

static struct luaL_Reg fLatInt[] = {
    { "Int",     q_latint },
    { NULL, NULL}
};

int
init_latint(lua_State *L)
{
    luaL_getmetatable(L, opLattice);
    luaL_register(L, NULL, fLatInt);
    lua_pop(L, 1);
    qlua_metatable(L, mtnLatInt, mtLatInt);
    qlua_metatable(L, opLatInt, LatIntMethods);
    qlua_reg_add(qLatInt, qLatInt, q_I_add_I);
    qlua_reg_add(qReal,   qLatInt, q_i_add_I);
    qlua_reg_add(qLatInt, qReal,   q_I_add_i);
    qlua_reg_sub(qLatInt, qLatInt, q_I_sub_I);
    qlua_reg_sub(qReal,   qLatInt, q_i_sub_I);
    qlua_reg_sub(qLatInt, qReal,   q_I_sub_i);
    qlua_reg_mul(qLatInt, qLatInt, q_I_mul_I);
    qlua_reg_mul(qReal,   qLatInt, q_i_mul_I);
    qlua_reg_mul(qLatInt, qReal,   q_I_mul_i);
    qlua_reg_div(qLatInt, qLatInt, q_I_div_I);
    qlua_reg_div(qReal,   qLatInt, q_i_div_I);
    qlua_reg_div(qLatInt, qReal,   q_I_div_i);
    qlua_reg_mod(qLatInt, qLatInt, q_I_mod_I);
    qlua_reg_mod(qReal,   qLatInt, q_i_mod_I);
    qlua_reg_mod(qLatInt, qReal,   q_I_mod_i);

    return 0;
}

int
fini_latint(lua_State *L)
{
    return 0;
}
