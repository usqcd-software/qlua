#include <qlua.h>                                                    /* DEPS */
#include <qvector.h>                                                 /* DEPS */
#include <lattice.h>                                                 /* DEPS */
#include <latint.h>                                                  /* DEPS */
#include <latmulti.h>                                                /* DEPS */
#include <stdlib.h>  
#include <qmp.h>
#include <string.h>

/* lattice integers */
const char mtnLatInt[] = "qcd.lattice.int";
static const char opLatInt[] = "qcd.lattice.int.op";

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

        QDP_r_eq_sum_I(&sum, a->ptr, qCurrent);
        lua_pushnumber(L, sum);
        
        return 1;
    }
    case 2: {
        mLatMulti *m = qlua_checkLatMulti(L, 2);
        mVecReal *r = qlua_newVecReal(L, m->count);

        r->size = m->count;
        QDP_r_eq_sum_I_multi(r->val, a->ptr, m->subset, m->count);

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

    QDP_I_eq_I(r->ptr, a->ptr, qCurrent);
    lua_pop(L, 1);

    return 1;
}

static int
q_I_norm2(lua_State *L)
{
    mLatInt *a = qlua_checkLatInt(L, 1);
    QLA_Real sum;

    QDP_r_eq_norm2_I(&sum, a->ptr, qCurrent);
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

    QDP_I_eq_sI(b->ptr, a->ptr, shift, dir, qCurrent);

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
    locked = QDP_expose_I(V->ptr);
    if (QDP_node_number(idx) == QDP_this_node) {
        QLA_elem_I(locked[QDP_index(idx)]) = z;
    }
    QDP_reset_I(V->ptr);
    qlua_free(L, idx);

    return 0;
}

static int
q_I_add_I(lua_State *L)
{
    mLatInt *res = qlua_newLatInt(L);
    mLatInt *a = qlua_checkLatInt(L, 1);
    mLatInt *b = qlua_checkLatInt(L, 2);

    QDP_I_eq_I_plus_I(res->ptr, a->ptr, b->ptr, qCurrent);

    return 1;
}

static int
q_I_sub_I(lua_State *L)
{
    mLatInt *res = qlua_newLatInt(L);
    mLatInt *a = qlua_checkLatInt(L, 1);
    mLatInt *b = qlua_checkLatInt(L, 2);

    QDP_I_eq_I_minus_I(res->ptr, a->ptr, b->ptr, qCurrent);

    return 1;
}

static int
q_i_mul_I(lua_State *L)
{
    mLatInt *res = qlua_newLatInt(L);
    QLA_Int a = luaL_checkint(L, 1);
    mLatInt *b = qlua_checkLatInt(L, 2);

    QDP_I_eq_i_times_I(res->ptr, &a, b->ptr, qCurrent);

    return 1;
}

static int
q_I_mul_i(lua_State *L)
{
    mLatInt *res = qlua_newLatInt(L);
    mLatInt *b = qlua_checkLatInt(L, 1);
    QLA_Int a = luaL_checkint(L, 2);

    QDP_I_eq_i_times_I(res->ptr, &a, b->ptr, qCurrent);

    return 1;
}

static int
q_I_mul_I(lua_State *L)
{
    mLatInt *res = qlua_newLatInt(L);
    mLatInt *a = qlua_checkLatInt(L, 1);
    mLatInt *b = qlua_checkLatInt(L, 2);

    QDP_I_eq_I_times_I(res->ptr, a->ptr, b->ptr, qCurrent);

    return 1;
}

static int
q_I_div_I(lua_State *L)
{
    mLatInt *res = qlua_newLatInt(L);
    mLatInt *a = qlua_checkLatInt(L, 1);
    mLatInt *b = qlua_checkLatInt(L, 2);

    QDP_I_eq_I_divide_I(res->ptr, a->ptr, b->ptr, qCurrent);

    return 1;
}

static int
q_I_neg(lua_State *L)
{
    mLatInt *a = qlua_checkLatInt(L, 1);
    mLatInt *res = qlua_newLatInt(L);
    QLA_Int m1 = -1;

    QDP_I_eq_i_times_I(res->ptr, &m1, a->ptr, qCurrent);

    return 1;
}

static int
q_I_dot(lua_State *L)
{
    mLatInt *a = qlua_checkLatInt(L, 1);
    mLatInt *b = qlua_checkLatInt(L, 2);
    QLA_Real s;

    QDP_r_eq_I_dot_I(&s, a->ptr, b->ptr, qCurrent);
    lua_pushnumber(L, s);

    return 1;
}

static int
q_latint(lua_State *L)
{
    switch (qlua_gettype(L, 2)) {
    case qReal: {
        QLA_Int d = luaL_checkint(L, 2);
        mLatInt *v = qlua_newLatInt(L);

        QDP_I_eq_i(v->ptr, &d, qCurrent);

        return 1;
    }
    case qLatInt: {
        mLatInt *res = qlua_newLatInt(L);
        mLatInt *a = qlua_checkLatInt(L, 2);

        QDP_I_eq_I(res->ptr, a->ptr, qCurrent);

        return 1;
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
    qlua_reg_sub(qLatInt, qLatInt, q_I_sub_I);
    qlua_reg_mul(qLatInt, qLatInt, q_I_mul_I);
    qlua_reg_mul(qReal,   qLatInt, q_i_mul_I);
    qlua_reg_mul(qLatInt, qReal,   q_I_mul_i);
    qlua_reg_div(qLatInt, qLatInt, q_I_div_I);
    qlua_reg_dot(qLatInt, q_I_dot);

    return 0;
}

int
fini_latint(lua_State *L)
{
    return 0;
}
