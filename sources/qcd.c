#include <qlua.h>
#include <stdlib.h>
#include <qcd.h>

static int qRank = 0;
static int *qDim = NULL;

/* lattice integers */
const char *mtnLatInt = "qcd.lattice.int";

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

    luaL_argcheck(L, v != 0, idx, "qcd.LatInt expected");
    
    return v;
}

static int
qLatInt_fmt(lua_State *L)
{
    char fmt[72];
    mLatInt *b = qlua_checkLatInt(L, 1);

    sprintf(fmt, "LatInt(%p)", b->ptr);
    lua_pushstring(L, fmt);

    return 1;
}

static int
qLatInt_gc(lua_State *L)
{
    mLatInt *b = qlua_checkLatInt(L, 1);

    QDP_destroy_I(b->ptr);
    b->ptr = 0;

    return 0;
}

static int
qLatInt_get(lua_State *L)
{
    switch (qlua_gettype(L, 2)) {
    case qTable: {
        mLatInt *V = qlua_checkLatInt(L, 1);
        QLA_Int *locked;
        int *idx = 0;
        int i, d;
        int z;
        
        d = lua_objlen(L, 2);
        if (d != qRank)
            return luaL_error(L, "parallel index rank mismatch");
        idx = qlua_malloc(L, d * sizeof (int));
        for (i = 0; i < d; i++) {
            lua_pushnumber(L, i + 1);
            lua_gettable(L, 2);
            idx[i] = luaL_checkinteger(L, -1);
            if ((idx[i] < 0) || (idx[i] >= qDim[i]))
                return luaL_error(L, "parallel index out of range");
        }
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
        lua_pushstring(L, "copy");
        if (lua_equal(L, 2, -1)) {
            lua_pushcfunction(L, q_I_eq_I);
            return 1;
        }
        /* THROUGH */
    default:
        return luaL_error(L, "bad index");
    }

}

static int
qLatInt_put(lua_State *L)
{
    mLatInt *V = qlua_checkLatInt(L, 1);
    QLA_Int *locked;
    int *idx = 0;
    int i, d;
    int z = luaL_checkinteger(L, 3);

    luaL_checktype(L, 2, LUA_TTABLE);
    d = lua_objlen(L, 2);
    if (d != qRank)
        return luaL_error(L, "parallel index rank mismatch");
    idx = qlua_malloc(L, d * sizeof (int));
    for (i = 0; i < d; i++) {
        lua_pushnumber(L, i + 1);
        lua_gettable(L, 2);
        idx[i] = luaL_checkinteger(L, -1);
        if ((idx[i] < 0) || (idx[i] >= qDim[i]))
            return luaL_error(L, "parallel index out of range");
    }
    locked = QDP_expose_I(V->ptr);
    if (QDP_node_number(idx) == QDP_this_node) {
        QLA_elem_I(locked[QDP_index(idx)]) = z;
    }
    QDP_reset_I(V->ptr);
    qlua_free(L, idx);

    return 0;
}

static int
qmk_latint(lua_State *L, QLA_Int value)
{
    mLatInt *v = qlua_newLatInt(L);

    QDP_I_eq_i(v->ptr, &value, QDP_all);

    return 1;
}

static int
q_copy_latint(lua_State *L, int idx)
{
    mLatInt *dst = qlua_newLatInt(L);
    mLatInt *src = qlua_checkLatInt(L, idx);

    QDP_I_eq_I(dst->ptr, src->ptr, QDP_all);

    return 1;
}

static int
qcd_latint(lua_State *L)
{
    switch (qlua_gettype(L, 1)) {
    case qReal:
        return qmk_latint(L, luaL_checkinteger(L, 1));
    case qLatInt:
        return q_copy_latint(L, 1);
    default:
        return luaL_error(L, "bad argument");
    }
}

int
q_I_eq_I(lua_State *L)
{
    mLatInt *res = qlua_newLatInt(L);
    mLatInt *a = qlua_checkLatInt(L, 1);

    QDP_I_eq_I(res->ptr, a->ptr, QDP_all);

    return 1;
}

int
q_I_add_I(lua_State *L)
{
    mLatInt *res = qlua_newLatInt(L);
    mLatInt *a = qlua_checkLatInt(L, 1);
    mLatInt *b = qlua_checkLatInt(L, 2);

    QDP_I_eq_I_plus_I(res->ptr, a->ptr, b->ptr, QDP_all);

    return 1;
}

int
q_I_sub_I(lua_State *L)
{
    mLatInt *res = qlua_newLatInt(L);
    mLatInt *a = qlua_checkLatInt(L, 1);
    mLatInt *b = qlua_checkLatInt(L, 2);

    QDP_I_eq_I_minus_I(res->ptr, a->ptr, b->ptr, QDP_all);

    return 1;
}

int
q_i_mul_I(lua_State *L)
{
    mLatInt *res = qlua_newLatInt(L);
    QLA_Int a = luaL_checkinteger(L, 1);
    mLatInt *b = qlua_checkLatInt(L, 2);

    QDP_I_eq_i_times_I(res->ptr, &a, b->ptr, QDP_all);

    return 1;
}

int
q_I_mul_i(lua_State *L)
{
    mLatInt *res = qlua_newLatInt(L);
    mLatInt *b = qlua_checkLatInt(L, 1);
    QLA_Int a = luaL_checkinteger(L, 2);

    QDP_I_eq_i_times_I(res->ptr, &a, b->ptr, QDP_all);

    return 1;
}

int
q_I_mul_I(lua_State *L)
{
    mLatInt *res = qlua_newLatInt(L);
    mLatInt *a = qlua_checkLatInt(L, 1);
    mLatInt *b = qlua_checkLatInt(L, 2);

    QDP_I_eq_I_times_I(res->ptr, a->ptr, b->ptr, QDP_all);

    return 1;
}

int
q_I_div_I(lua_State *L)
{
    mLatInt *res = qlua_newLatInt(L);
    mLatInt *a = qlua_checkLatInt(L, 1);
    mLatInt *b = qlua_checkLatInt(L, 2);

    QDP_I_eq_I_divide_I(res->ptr, a->ptr, b->ptr, QDP_all);

    return 1;
}

int
q_neg_I(lua_State *L)
{
    mLatInt *res = qlua_newLatInt(L);
    mLatInt *a = qlua_checkLatInt(L, 1);
    QLA_Int m1 = -1;

    QDP_I_eq_i_times_I(res->ptr, &m1, a->ptr, QDP_all);

    return 1;
}

static struct luaL_reg mtLatInt[] = {
    { "__tostring",   qLatInt_fmt },
    { "__gc",         qLatInt_gc },
    { "__index",      qLatInt_get },
    { "__newindex",   qLatInt_put },
    { "__umn",        q_neg_I },
    { "__add",        qlua_add },
    { "__sub",        qlua_sub },
    { "__mul",        qlua_mul },
    { "__div",        qlua_div },
    /* handled by __index { "copy",         q_I_eq_I }, */
    { NULL,           NULL}
};

/* lattice definition */
static int
qcd_lattice(lua_State *L)
{
    int r, i;

    if (qRank != 0)
        return luaL_error(L, "redefining the lattice");

    luaL_checktype(L, 1, LUA_TTABLE);
    r = lua_objlen(L, 1);
    if (r <= 0)
        return luaL_error(L, "Bad lattice rank");
    qRank = r;
    qDim = qlua_malloc(L, r * sizeof (int));
    for (i = 0; i < r; i++) {
        lua_pushnumber(L, i + 1);
        lua_gettable(L, 1);
        qDim[i] = luaL_checkinteger(L, -1);
    }
    QDP_set_latsize(qRank, qDim);
    if (QDP_create_layout()) {
        return luaL_error(L, "can not create lattice");
    }
    
    return 0;
}

static int
qcd_dims(lua_State *L)
{
    int i;

    lua_createtable(L, qRank, 0);
    for (i = 0; i < qRank; i++) {
        lua_pushnumber(L, qDim[i]);
        lua_rawseti(L, -2, i + 1);
    }
    return 1;
}

/* pcoords */
static int pcoord_d = -1;
static void
pcoord_set(QLA_Int *dst, int coords[])
{
    *dst = coords[pcoord_d];
}

static int
qcd_pcoord(lua_State *L)
{
    int d = luaL_checkint(L, 1);
    mLatInt *v = qlua_newLatInt(L);
    
    if ((d < 0) || (d >= qRank))
        return luaL_error(L, "coordinate out of range");
    
    /* YYY global state */
    pcoord_d = d;
    QDP_I_eq_func(v->ptr, pcoord_set, QDP_all);

    return 1;
}

static struct luaL_Reg fQcd[] = {
    { "lattice", qcd_lattice },
    { "dims",    qcd_dims },
    { "pcoord",  qcd_pcoord },
    { "lat_int", qcd_latint },
    { NULL, NULL}
};

int
init_qcd(lua_State *L)
{
    luaL_register(L, qcdlib, fQcd);
    qlua_metatable(L, mtnLatInt, mtLatInt);

    return 0;
}

int
fini_qcd(lua_State *L)
{
    qlua_free(L, qDim);
    qDim = 0;
    qRank = 0;
}
