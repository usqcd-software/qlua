#include "qlua.h"                                                    /* DEPS */
#include "lattice.h"                                                 /* DEPS */
#include "latsubset.h"                                               /* DEPS */
#include "latint.h"                                                  /* DEPS */
#include <string.h>

static const char LatSubsetName[] = "lattice.Subset";

static int
q_U_gc(lua_State *L)
{
    mLatSubset *m = qlua_checkLatSubset(L, 1, NULL);

    if (m->mask)
        QDP_destroy_I(m->mask);
    m->mask = 0;

    return 0;
}

static int
q_U_fmt(lua_State *L)
{
    mLatSubset *m = qlua_checkLatSubset(L, 1, NULL);
    char fmt[72];

    switch (m->cl) {
    case qss_all:
        sprintf(fmt, "Subset[all%s]", m->mask? ",mask": "");
        break;
    case qss_even:
        sprintf(fmt, "Subset[even%s]", m->mask? ",mask": ""); 
        break;
    case qss_odd: 
        sprintf(fmt, "Subset[odd%s]", m->mask? ",mask": "");
        break;
    case qss_none:
        sprintf(fmt, "Subset[none]");
        break;
    case qss_slice:
        sprintf(fmt, "Subset[axis=%d,position=%d%s]",
                m->axis, m->position, m->mask? ",mask": "");
        break;
    case qss_upper:
        sprintf(fmt, "Subset[axis=%d,position=%d,semispace=upper%s]",
                m->axis, m->position, m->mask? ",mask": "");
        break;
    case qss_lower:
        sprintf(fmt, "Subset[axis=%d,position=%d,semispace=lower%s]",
                m->axis, m->position, m->mask? ",mask": "");
        break;
    }
    lua_pushstring(L, fmt);

    return 1;
}

/* local subset function returns 0 on positions belonging to the subset */
static int
subset_local(QDP_Lattice *lat, int *coord, void *arg)
{
    mLatSubset *x = arg;

    switch (x->cl) {
    case qss_slice:
        return coord[x->axis] != x->position;
    case qss_upper:
        return coord[x->axis] < x->position;
    case qss_lower:
        return coord[x->axis] >= x->position;
    default:
        return -1;
    }
}

static QDP_Int *
subset_mask(lua_State *L, mLattice *S, QDP_Int *w)
{
    QDP_Int *v = QDP_create_I_L(S->lat);
    if (v == 0) {
        CALL_QDP(L);
        v = QDP_create_I_L(S->lat);
        if (v == 0)
            luaL_error(L, "not enough memory (subset mask)");
    }
    QDP_I_eq_zero(v, S->all);
    if (w)
        QDP_I_eq_I(v, w, *S->qss);

    return v;
}

static void
subset_copy(lua_State *L, mLattice *S, mLatSubset *dst, const mLatSubset *src)
{
    *dst = *src;
    dst->mask = 0;
    if (src->mask) {
        dst->mask = subset_mask(L, S, src->mask);
    }
}

typedef struct {
    QLA_Int *a_mask;
    QLA_Int *b_mask;
    mLatSubset b;
    QDP_Lattice *lat;
    int rank;
} sjArg;

static void
build_mask(QLA_Int *r, int coord[], void *env)
{
    sjArg *arg = env;
    int idx = QDP_index_L(arg->lat, coord);
    int i;
    int p;
    int v = 0;

    for (i = 0, p = 0; i < arg->rank; i++)
        p += coord[i];

    switch (arg->b.cl) {
    case qss_slice:
        v = (coord[arg->b.axis] == arg->b.position);
        break;
    case qss_upper:
        v = (coord[arg->b.axis] >= arg->b.position);
        break;
    case qss_lower:
        v = (coord[arg->b.axis] < arg->b.position);
        break;
    case qss_even:
        v = ((p & 1) == 0);
        break;
    case qss_odd:
        v = ((p & 1) == 1);
        break;
    case qss_all:
        v = 1;
        break;
    case qss_none:
        v = 0;
        break;
    }
    if (arg->a_mask)
        v = v & arg->a_mask[idx];
    if (arg->b_mask)
        v = v & arg->b_mask[idx];

    *r = v;
}

static void
subset_join(lua_State *L,
            mLatSubset *dst,
            mLattice *S,  /* use S->lss as src1 */
            const mLatSubset *b) /* src2 */
{
    int refine = ((b->mask != NULL) || (S->lss.mask != NULL))? 1: 0;
    QLA_D_Real count = 0;
    sjArg arg;

    if (b->cl == qss_none) {
        subset_copy(L, S, dst, b);
        return;
    }

    subset_copy(L, S, dst, &S->lss);
    switch (S->lss.cl) {
    case qss_none:
        return;
    case qss_all:
        if (S->lss.mask == 0) {
            subset_copy(L, S, dst, b);
            return;
        }
        break;
    case qss_even:
        if (b->cl == qss_odd) {
            if (dst->mask)
                QDP_destroy_I(dst->mask);
            dst->mask = NULL;
            dst->cl = qss_none;
            return;
        }
        break;
    case qss_odd:
        if (b->cl == qss_even) {
            if (dst->mask)
                QDP_destroy_I(dst->mask);
            dst->mask = NULL;
            dst->cl = qss_none;
            return;
        }
        break;
    case qss_slice:
    case qss_lower:
    case qss_upper:
        refine = 1;
        break;
    }
    if (!refine)
        return;
    arg.lat = S->lat;
    arg.rank = S->rank;
    arg.a_mask = S->lss.mask ? QDP_expose_I(S->lss.mask): NULL;
    arg.b_mask = b->mask ? QDP_expose_I(b->mask): NULL;
    arg.b = *b;
    if (dst->mask == 0)
        dst->mask = subset_mask(L, S, NULL);
    QDP_I_eq_funca(dst->mask, build_mask, &arg, *S->qss);
    if (S->lss.mask) QDP_reset_I(S->lss.mask);
    if (b->mask) QDP_reset_I(b->mask);
    
    QDP_r_eq_sum_I(&count, dst->mask, *S->qss);
    if (count == 0) {
        QDP_destroy_I(dst->mask);
        dst->cl = qss_none;
    }
}

static void
switch_subset(lua_State *L, mLattice *S, const mLatSubset *v)
{
    if (S->lss.cl > qss_last_static)
        QDP_destroy_subset(S->qss);

    S->lss = *v;
    switch (S->lss.cl) {
    case qss_none: S->qss = S->none; break;
    case qss_all:  S->qss = &S->all;  break;
    case qss_even: S->qss = &S->even; break;
    case qss_odd:  S->qss = &S->odd;  break;
    case qss_lower:
    case qss_upper:
    case qss_slice:
        S->qss = QDP_create_subset_L(S->lat, subset_local,
                                     &S->lss, sizeof (S->lss), 1);
        if (S->qss == 0) {
            lua_gc(L, LUA_GCCOLLECT, 0);
            S->qss = QDP_create_subset_L(S->lat, subset_local,
                                         &S->lss, sizeof (S->lss), 1);
            if (S->qss == 0)
                luaL_error(L, "QDP_create_subset() failed");
        }
        break;
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
    mLatSubset new_subset;

    subset_join(L, &new_subset, S, v);
    switch_subset(L, S, &new_subset);
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
    static const mLatSubset all = { qss_all, 0, 0, NULL };

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
    case qLatInt: {
        mLatInt *mask = qlua_checkLatInt(L, 2, S);
        mLatSubset *x = qlua_newLatSubset(L, 1);
        m->mask = subset_mask(L, S, mask->ptr);
        subset_join(L, x, S, m);
        return 1;
    }
    default:
        break;
    }
    return qlua_badconstr(L, "Subset");
}

static struct luaL_Reg mtLatSubset[] = {
    { "__gc",           q_U_gc        },
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
    v->mask = NULL;
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
