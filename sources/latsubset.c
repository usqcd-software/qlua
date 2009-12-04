#include <qlua.h>                                                    /* DEPS */
#include <lattice.h>                                                 /* DEPS */
#include <latsubset.h>                                               /* DEPS */
#include <string.h>

const char mtnLatSubset[] = "lattice.subset";
static const char LatSubsetStack[] = "lattice.subsetstack";

static mLatSubset current_set;

static mLatSubset *
qlua_newLatSubset(lua_State *L)
{
    mLatSubset *v = lua_newuserdata(L, sizeof (mLatSubset));
    
    v->cl = qss_all;
    luaL_getmetatable(L, mtnLatSubset);
    lua_setmetatable(L, -2);

    return v;
}

static mLatSubset *
qlua_checkLatSubset(lua_State *L, int idx)
{
    void *v = luaL_checkudata(L, idx, mtnLatSubset);

    luaL_argcheck(L, v != 0, idx, "lattice.Subset expected");

    return v;
}

static int
q_U_fmt(lua_State *L)
{
    mLatSubset *m = qlua_checkLatSubset(L, 1);
    char fmt[72];

    switch (m->cl) {
    case qss_all:
        sprintf(fmt, "Subset[all]");
        break;
    case qss_even:
        sprintf(fmt, "Subset[even]"); 
        break;
    case qss_odd: 
        sprintf(fmt, "Subset[odd]");
        break;
    case qss_slice:
        sprintf(fmt, "Subset[axis=%d,position=%d]", m->axis, m->position);
        break;
    case qss_upper:
        sprintf(fmt, "Subset[axis=%d,position=%d,semispace=upper]",
                m->axis, m->position);
        break;
    case qss_lower:
        sprintf(fmt, "Subset[axis=%d,position=%d,semispace=lower]",
                m->axis, m->position);
        break;
    default:
        return luaL_error(L, "unknown subset class %d", m->cl);
    }
    lua_pushstring(L, fmt);

    return 1;
}

/* subset functions return 0 on positions belonging to the subset */
static int
subset_slice(int *coord, void *arg)
{
    mLatSubset *x = arg;

    return coord[x->axis] != x->position;
}

static int
subset_upper(int *coord, void *arg)
{
    mLatSubset *x = arg;

    return coord[x->axis] < x->position;
}

static int
subset_lower(int *coord, void *arg)
{
    mLatSubset *x = arg;

    return coord[x->axis] >= x->position;
}

static void
switch_subset(lua_State *L, mLatSubset *s)
{
    if (current_set.cl > qss_dynamic)
        QDP_destroy_subset(qCurrent);

    switch (s->cl) {
    case qss_all:  qCurrent = &QDP_all;  break;
    case qss_even: qCurrent = &QDP_even; break;
    case qss_odd:  qCurrent = &QDP_odd;  break;
    case qss_slice:
        qCurrent = QDP_create_subset(subset_slice, s, sizeof (*s), 1);
        if (qCurrent == 0) {
            lua_gc(L, LUA_GCCOLLECT, 0);
            qCurrent = QDP_create_subset(subset_slice, s, sizeof (*s), 1);
            if (qCurrent == 0)
                luaL_error(L, "QDP_create_subset() failed");
        }
        break;
    case qss_upper:
        qCurrent = QDP_create_subset(subset_upper, s, sizeof (*s), 1);
        if (qCurrent == 0) {
            lua_gc(L, LUA_GCCOLLECT, 0);
            qCurrent = QDP_create_subset(subset_upper, s, sizeof (*s), 1);
            if (qCurrent == 0)
                luaL_error(L, "QDP_create_subset() failed");
        }
        break;
    case qss_lower:
        qCurrent = QDP_create_subset(subset_lower, s, sizeof (*s), 1);
        if (qCurrent == 0) {
            lua_gc(L, LUA_GCCOLLECT, 0);
            qCurrent = QDP_create_subset(subset_lower, s, sizeof (*s), 1);
            if (qCurrent == 0)
                luaL_error(L, "QDP_create_subset() failed");
        }
        break;
    default:
        luaL_error(L, "unknown subset class %d", s->cl);
    }
    current_set = *s;
}

static int
q_U_where(lua_State *L)
{
    mLatSubset *v = qlua_checkLatSubset(L, 1);
    int argc = lua_gettop(L) - 2;
    int resc;
    mLatSubset old_subset = current_set;

    switch_subset(L, v);

    if (lua_pcall(L, argc, LUA_MULTRET, 0))
        return luaL_error(L, lua_tostring(L, -1));
    resc = lua_gettop(L) - 1;

    switch_subset(L, &old_subset);

    return resc;
}

static int
q_subset(lua_State *L)
{
    int d = qlua_checkindex(L, 2, "axis", qRank);
    int p = qlua_checkindex(L, 2, "position", qDim[d]);
    mLatSubset *m = qlua_newLatSubset(L);
    const char *sub = 0;

    lua_getfield(L, 2, "semispace");
    if (lua_type(L, -1) == LUA_TSTRING)
        sub = lua_tostring(L, -1);
    lua_pop(L, 1);

    if (sub == 0)
        m->cl = qss_slice;
    else if (strcmp(sub, "upper") == 0)
        m->cl = qss_upper;
    else if (strcmp(sub, "lower") == 0)
        m->cl = qss_lower;
    else
        luaL_error(L, "bad semispace specifier");

    m->axis = d;
    m->position = p;
    luaL_getmetatable(L, mtnLatSubset);
    lua_setmetatable(L, -2);

    return 1;
}

static void
add_static(lua_State *L, qSubsetClass cl, const char *name)
{
    mLatSubset *m = qlua_newLatSubset(L);

    m->cl = cl;
    luaL_getmetatable(L, mtnLatSubset);
    lua_setmetatable(L, -2);

    lua_setfield(L, -2, name);
}

static struct luaL_Reg mtLatSubset[] = {
    { "__tostring",     q_U_fmt   },
    { "where",          q_U_where },
    { NULL,             NULL      }
};

static struct luaL_Reg fLatSubset[] = {
    { "Subset",    q_subset },
    { NULL,        NULL     }
};

int
init_latsubset(lua_State *L)
{    
    luaL_newmetatable(L, LatSubsetStack);
    qlua_metatable(L, mtnLatSubset, mtLatSubset);
    luaL_getmetatable(L, opLattice);
    luaL_register(L, NULL, fLatSubset);
    add_static(L, qss_all, "all");
    add_static(L, qss_even, "even");
    add_static(L, qss_odd, "odd");
    lua_pop(L, 1);

    current_set.cl = qss_all;
    
    
    return 0;
}

int
fini_latsubset(lua_State *L)
{
    return 0;
}
