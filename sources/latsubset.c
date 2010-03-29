#include "qlua.h"                                                    /* DEPS */
#include "lattice.h"                                                 /* DEPS */
#include "latsubset.h"                                               /* DEPS */
#include <string.h>

static const char LatSubsetName[] = "lattice.Subset";

static int
q_U_fmt(lua_State *L)
{
    mLatSubset *m = qlua_checkLatSubset(L, 1, NULL);
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
switch_subset(lua_State *L, mLattice *S, mLatSubset *v)
{
    if (S->lss.cl > qss_dynamic)
        QDP_destroy_subset(S->qss);

    S->lss = *v;
    switch (S->lss.cl) {
    case qss_all:  S->qss = &QDP_all;  break;
    case qss_even: S->qss = &QDP_even; break;
    case qss_odd:  S->qss = &QDP_odd;  break;
    case qss_slice:
        S->qss = QDP_create_subset(subset_slice, v, sizeof (*v), 1);
        if (S->qss == 0) {
            lua_gc(L, LUA_GCCOLLECT, 0);
            S->qss = QDP_create_subset(subset_slice, v, sizeof (*v), 1);
            if (S->qss == 0)
                luaL_error(L, "QDP_create_subset() failed");
        }
        break;
    case qss_upper:
        S->qss = QDP_create_subset(subset_upper, v, sizeof (*v), 1);
        if (S->qss == 0) {
            lua_gc(L, LUA_GCCOLLECT, 0);
            S->qss = QDP_create_subset(subset_upper, v, sizeof (*v), 1);
            if (S->qss == 0)
                luaL_error(L, "QDP_create_subset() failed");
        }
        break;
    case qss_lower:
        S->qss = QDP_create_subset(subset_lower, v, sizeof (*v), 1);
        if (S->qss == 0) {
            lua_gc(L, LUA_GCCOLLECT, 0);
            S->qss = QDP_create_subset(subset_lower, v, sizeof (*v), 1);
            if (S->qss == 0)
                luaL_error(L, "QDP_create_subset() failed");
        }
        break;
    default:
        luaL_error(L, "unknown subset class %d", v->cl);
    }
}

static int
q_U_where(lua_State *L)
{
    mLatSubset *v = qlua_checkLatSubset(L, 1, NULL);
    int argc = lua_gettop(L) - 2;
    int resc;
    mLattice *S = qlua_ObjLattice(L, 1);
    mLatSubset old_subset = S->lss;
    
    switch_subset(L, S, v);
    lua_pop(L, 1);
    if (lua_pcall(L, argc, LUA_MULTRET, 0))
        return luaL_error(L, lua_tostring(L, -1));
    resc = lua_gettop(L) - 1;
    switch_subset(L, S, &old_subset);

    return resc;
}

int
qlua_everywhere(lua_State *L)
{
    mLattice *S = qlua_checkLattice(L, 1);
    mLatSubset old_subset = S->lss;
    int argc = lua_gettop(L) - 1;
    int resc;
    mLatSubset all;

    all.cl = qss_all;
    switch_subset(L, S, &all);
    if (lua_pcall(L, argc, LUA_MULTRET, 0))
        return luaL_error(L, lua_tostring(L, -1));
    resc = lua_gettop(L) - 1;
    switch_subset(L, S, &old_subset);

    return resc;
}

static int
q_subset(lua_State *L)
{
    mLattice *S = qlua_checkLattice(L, 1);
    mLatSubset *m = qlua_newLatSubset(L, 1);

    switch (qlua_qtype(L, 2)) {
    case qTable: {
        int d = qlua_checkindex(L, 2, "axis", S->rank);
        int p = qlua_checkindex(L, 2, "position", S->dim[d]);
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
        return 1;
    }
    case qString: {
        const char *sset = lua_tostring(L, 2);
        if (strcmp(sset, "all") == 0)
            m->cl = qss_all;
        else if (strcmp(sset, "even") == 0)
            m->cl = qss_even;
        else if (strcmp(sset, "odd") == 0)
            m->cl = qss_odd;
        else
            break;
        return 1;
    }
    default:
        break;
    }
    return qlua_badconstr(L, "Subset");
}

static struct luaL_Reg mtLatSubset[] = {
    { "__tostring",     q_U_fmt       },
    { "__newindex",     qlua_nowrite  },
    { "where",          q_U_where     },
    /* "lattice" */
    /* *a-type" */
    { NULL,             NULL          }
};

mLatSubset *
qlua_newLatSubset(lua_State *L, int Sidx)
{
    mLatSubset *v = lua_newuserdata(L, sizeof (mLatSubset));
    
    v->cl = qss_all;
    qlua_createLatticeTable(L, Sidx, mtLatSubset, qLatSubset, LatSubsetName);
    lua_setmetatable(L, -2);

    return v;
}

mLatSubset *
qlua_checkLatSubset(lua_State *L, int idx, mLattice *S)
{
    void *v = qlua_checkLatticeType(L, idx, qLatSubset, LatSubsetName);
    
    if (S) {
        mLattice *S1 = qlua_ObjLattice(L, idx);
        if (S1->id != S->id)
            luaL_error(L, "%s on a wrong lattice", LatSubsetName);
        lua_pop(L, 1);
    }

    return (mLatSubset *)v;
}

static struct luaL_Reg fLatSubset[] = {
    { "Subset",    q_subset },
    { NULL,        NULL     }
};

int
init_latsubset(lua_State *L)
{    
    luaL_getmetatable(L, opLattice);
    luaL_register(L, NULL, fLatSubset);
    lua_pop(L, 1);

    return 0;
}

int
fini_latsubset(lua_State *L)
{
    return 0;
}
