#include <qlua.h>
#include <latint.h>
#include <stdlib.h>
#include <qmp.h>
#include <string.h>

static int qRank = 0;
static int *qDim = NULL;

/* lattice integers */
const char *mtnLatInt = "qcd.lattice.int";
static const char *opLatInt = "qcd.lattice.int.op";

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
q_I_sum(lua_State *L)
{
    mLatInt *a = qlua_checkLatInt(L, 1);
    QLA_D_Real sum;


    QDP_D_r_eq_sum_I(&sum, a->ptr, QDP_all);
    lua_pushnumber(L, sum);

    return 1;
}

static int
q_I_norm2(lua_State *L)
{
    mLatInt *a = qlua_checkLatInt(L, 1);
    QLA_D_Real sum;

    QDP_D_r_eq_norm2_I(&sum, a->ptr, QDP_all);
    lua_pushnumber(L, sum);

    return 1;
}


QDP_Shift
qlua_checkShift(lua_State *L, int idx)
{
    int d = luaL_checkinteger(L, idx);

    if ((d < 0) || (d >= qRank))
        luaL_error(L, "bad shift dimension");

    return QDP_neighbor[d];
}

QDP_ShiftDir
qlua_checkShiftDir(lua_State *L, int idx)
{
    static const struct {
        char *name;
        QDP_ShiftDir dir;
    } t[] = {
        { "from_forward",  QDP_forward },
        { "from_backward", QDP_backward },
        { "to_forward",    QDP_backward },
        { "to_backward",   QDP_forward },
        { NULL,            QDP_forward }
    };
    int i;
    const char *d = luaL_checkstring(L, idx);

    for (i = 0; t[i].name; i++) {
        if (strcmp(d, t[i].name) == 0)
            return t[i].dir;
    }
    luaL_error(L, "bad shift direction");
    /* NEVER HAPPENS */
    return QDP_forward;
}

static int
q_I_shift(lua_State *L)
{
    mLatInt *a = qlua_checkLatInt(L, 1);
    QDP_Shift shift = qlua_checkShift(L, 2);
    QDP_ShiftDir dir = qlua_checkShiftDir(L, 3);
    mLatInt *b = qlua_newLatInt(L);

    QDP_I_eq_sI(b->ptr, a->ptr, shift, dir, QDP_all);

    return 1;
}

int *
qlua_lattice_coord(lua_State *L, int n)
{
    int d, i;
    int *idx;

    luaL_checktype(L, n, LUA_TTABLE);
    d = lua_objlen(L, n);
    if (d != qRank) {
        luaL_error(L, "parallel index rank mismatch");
        return NULL;
    }
    idx = qlua_malloc(L, d * sizeof (int));
    for (i = 0; i < d; i++) {
        lua_pushnumber(L, i + 1);
        lua_gettable(L, n);
        idx[i] = luaL_checkinteger(L, -1);
        if ((idx[i] < 0) || (idx[i] >= qDim[i])) {
            qlua_free(L, idx);
            luaL_error(L, "parallel index out of range");
            return NULL;
        }
    }
    
    return idx;
}

static int
qLatInt_get(lua_State *L)
{
    switch (qlua_gettype(L, 2)) {
    case qTable: {
        mLatInt *V = qlua_checkLatInt(L, 1);
        QLA_Int *locked;
        int *idx = 0;
        int z;

        idx = qlua_lattice_coord(L, 2);
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

    return luaL_error(L, "bad index");
}

static struct luaL_Reg LatIntMethods[] = {
    { "norm2",  q_I_norm2 },
    { "shift",  q_I_shift },
    { "sum",    q_I_sum },
    { NULL,     NULL}
};

static int
qLatInt_put(lua_State *L)
{
    mLatInt *V = qlua_checkLatInt(L, 1);
    QLA_Int *locked;
    int *idx = 0;
    int z = luaL_checkinteger(L, 3);

    idx = qlua_lattice_coord(L, 2);
    locked = QDP_expose_I(V->ptr);
    if (QDP_node_number(idx) == QDP_this_node) {
        QLA_elem_I(locked[QDP_index(idx)]) = z;
    }
    QDP_reset_I(V->ptr);
    qlua_free(L, idx);

    return 0;
}

static int
q_latint(lua_State *L)
{
    switch (qlua_gettype(L, 1)) {
    case qReal: {
        QLA_Int d = luaL_checkinteger(L, 1);
        mLatInt *v = qlua_newLatInt(L);

        QDP_I_eq_i(v->ptr, &d, QDP_all);
        break;
    }
    case qLatInt: {
        mLatInt *res = qlua_newLatInt(L);
        mLatInt *a = qlua_checkLatInt(L, 1);

        QDP_I_eq_I(res->ptr, a->ptr, QDP_all);
        break;
    }
    default:
        return luaL_error(L, "bad argument");
    }
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
    mLatInt *a = qlua_checkLatInt(L, 1);
    mLatInt *res = qlua_newLatInt(L);
    QLA_Int m1 = -1;

    QDP_I_eq_i_times_I(res->ptr, &m1, a->ptr, QDP_all);

    return 1;
}

int
q_I_dot(lua_State *L)
{
    mLatInt *a = qlua_checkLatInt(L, 1);
    mLatInt *b = qlua_checkLatInt(L, 2);
    QLA_D_Real s;

    QDP_D_r_eq_I_dot_I(&s, a->ptr, b->ptr, QDP_all);
    lua_pushnumber(L, s);

    return 1;
}

static struct luaL_Reg mtLatInt[] = {
    { "__tostring",   qLatInt_fmt },
    { "__gc",         qLatInt_gc },
    { "__index",      qLatInt_get },
    { "__newindex",   qLatInt_put },
    { "__umn",        q_neg_I },
    { "__add",        qlua_add },
    { "__sub",        qlua_sub },
    { "__mul",        qlua_mul },
    { "__div",        qlua_div },
    { NULL,           NULL}
};

/* lattice definition */
static int
q_lattice(lua_State *L)
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
q_dims(lua_State *L)
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
q_pcoord(lua_State *L)
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

static struct luaL_Reg fLatInt[] = {
    { "lattice", q_lattice },
    { "dims",    q_dims },
    { "pcoord",  q_pcoord },
    { "lat_int", q_latint },
    { NULL, NULL}
};

int
init_latint(lua_State *L)
{
    luaL_register(L, qcdlib, fLatInt);
    qlua_metatable(L, mtnLatInt, mtLatInt);
    qlua_metatable(L, opLatInt, LatIntMethods);

    return 0;
}

int
fini_latint(lua_State *L)
{
    qlua_free(L, qDim);
    qDim = 0;
    qRank = 0;

    return 0;
}
