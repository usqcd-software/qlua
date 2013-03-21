#include "modules.h"                        /* DEPS */
#include "qlua.h"                           /* DEPS */
#include "qhql.h"                           /* DEPS */
#include "lattice.h"                        /* DEPS */

static const char *HQLGridName = "hql";

static int
q_g_fmt(lua_State *L)
{
  qlua_checkHQLGrid(L, 1);
  lua_pushstring(L, "HQLGrid{}");
  return 1;
}

static int
q_g_gc(lua_State *L)
{
  qlua_checkHQLGrid(L, 1);
  return 0;
}

static int
q_g_get(lua_State *L)
{
  return qlua_selflookup(L, 1, luaL_checkstring(L, 2));
}

static struct luaL_Reg mtHQLGrid[] = {
    { "__tostring",        q_g_fmt      },
    { "__gc",              q_g_gc       },
    { "__index",           q_g_get      },
    { "__newindex",        qlua_nowrite },
    { "matrix",            qhql_matrix  },
    /* "lattice" */
    /* "a-type" */
    { NULL,                NULL         }
};

mHQLGrid *
qlua_newHQLGrid(lua_State *L, int Sidx, int flavor_dim, int spin_dim)
{
    mHQLGrid *grid = lua_newuserdata(L, sizeof (mHQLGrid));

    grid->flavor_dim = flavor_dim;
    grid->spin_dim = spin_dim;

    qlua_checkLattice(L, Sidx);
    qlua_createLatticeTable(L, Sidx, mtHQLGrid, qHQLGrid, HQLGridName);
    lua_setmetatable(L, -2);

    return grid;
}

mHQLGrid *
qlua_checkHQLGrid(lua_State *L, int idx)
{
  void *v = qlua_checkLatticeType(L, idx, qHQLGrid, HQLGridName);
    
  return (mHQLGrid *)v;
}

static int
q_get_intopt(lua_State *L, int tidx, const char *name, int def)
{
  lua_pushstring(L, name);
  lua_gettable(L, tidx);
  return luaL_optint(L, -1, def);
}

static int
q_hql(lua_State *L)
{
  int Sidx;
  int flavor_dim;
  int spin_dim;

  if (lua_istable(L, 1) == 0)
    luaL_error(L, "qcd.hql expects table");
  lua_pushstring(L, "Lattice");
  lua_gettable(L, 1);
  Sidx = lua_gettop(L);
  if (qlua_qtype(L, Sidx) != qLattice)
    luaL_error(L, "qcd.hql expects .Lattice");
  flavor_dim = q_get_intopt(L, 1, "Flavor", 1);
  spin_dim = q_get_intopt(L, 1, "Spin", 4);
  if (flavor_dim < 1)
    luaL_error(L, "qcd.hql expects .Flavor >= 1");
  if ((spin_dim != 1) && (spin_dim != 4))
    luaL_error(L, "qcd.hql expects .Spin to be 1 or 4");

  qlua_newHQLGrid(L, Sidx, flavor_dim, spin_dim);
  return 1;
}

static const struct luaL_Reg fHQL[] = {
  { "hql", q_hql },
  { NULL, NULL }
};

int
init_hql(lua_State *L)
{
  luaL_register(L, qcdlib, fHQL);
  /* XXX install hql matrix ops */
  /* XXX install hql vector ops */
  return 0;
}

int
fini_hql(lua_State *L)
{
  return 0;
}
