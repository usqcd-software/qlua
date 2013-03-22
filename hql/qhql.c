#include "modules.h"                        /* DEPS */
#include "qlua.h"                           /* DEPS */
#include "lattice.h"                        /* DEPS */
#include "qhql.h"                           /* DEPS */

static const char *HQLGridName = "hql";
static const char *grid_key = "grid";

static int
q_g_fmt(lua_State *L)
{
  char fmt[72];
  mHQLGrid *grid = qlua_checkHQLGrid(L, 1);

  sprintf(fmt, "HQLGrid(%s)", grid->grid? "" : "template");
  lua_pushstring(L, fmt);
  return 1;
}

static int
q_g_gc(lua_State *L)
{
  mHQLGrid *grid = qlua_checkHQLGrid(L, 1);

  if (grid->grid)
    HQL_GridFree(&grid->grid);
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
qlua_newHQLGrid(lua_State *L, int Sidx, int flavor_dim, int spin_dim, int colors)
{
  mHQLGrid *grid = lua_newuserdata(L, sizeof (mHQLGrid));
  
  grid->flavor_dim = flavor_dim;
  grid->spin_dim = spin_dim;
  grid->colors = colors;
  grid->grid = NULL;
  
  qlua_checkLattice(L, Sidx);
  qlua_createLatticeTable(L, Sidx, mtHQLGrid, qHQLGrid, HQLGridName);
  lua_setmetatable(L, -2);
  
  return grid;
}

mHQLGrid *
qlua_copyHQLGrid(lua_State *L, int idx)
{
  mHQLGrid *templ = qlua_checkHQLGridTemplate(L, idx);
  mHQLGrid *grid = lua_newuserdata(L, sizeof (mHQLGrid));

  grid->flavor_dim = templ->flavor_dim;
  grid->spin_dim = templ->spin_dim;
  grid->colors = templ->colors;
  grid->grid = NULL;

  lua_getmetatable(L, idx);
  lua_setmetatable(L, -2);

  return grid;
}

mHQLGrid *
qlua_checkHQLGridTemplate(lua_State *L, int idx)
{
  mHQLGrid *grid = qlua_checkLatticeType(L, idx, qHQLGrid, HQLGridName);
  
  if (grid->grid)
    luaL_error(L, "assembled grid instead of template");

  return grid;
}

mHQLGrid *
qlua_checkHQLGrid(lua_State *L, int idx)
{
  void *v = qlua_checkLatticeType(L, idx, qHQLGrid, HQLGridName);
    
  return (mHQLGrid *)v;
}

void
qlua_createHQLTable(lua_State *L,
                         int Sidx,
                         int Gidx,
                         const struct luaL_Reg *ft,
                         QLUA_Type t_id,
                         const char *name)
{
  if (luaL_getmetafield(L, Gidx, name) != 0)
    goto end;
  lua_getmetatable(L, Gidx);
  qlua_latticetable(L, ft, t_id, Sidx);
  lua_pushstring(L, grid_key);
  lua_pushvalue(L, Gidx);
  lua_settable(L, -3);
  lua_setfield(L, -2, name);
  lua_pop(L, 1);
  luaL_getmetafield(L, Gidx, name);
 end:
  return;
}

mHQLGrid *
qlua_ObjGrid(lua_State *L, int idx)
{
  if (luaL_getmetafield(L, idx, grid_key) == 0)
    luaL_error(L, "hql object expected");
  return qlua_checkHQLGrid(L, -1);
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
  int colors;

  if (lua_istable(L, 1) == 0)
    luaL_error(L, "qcd.hql expects table");
  lua_pushstring(L, "Lattice");
  lua_gettable(L, 1);
  Sidx = lua_gettop(L);
  if (qlua_qtype(L, Sidx) != qLattice)
    luaL_error(L, "qcd.hql expects .Lattice");
  lua_pushstring(L, "Colors");
  lua_gettable(L, 1);
  colors = luaL_checkint(L, -1);
  lua_pop(L, 1);
  if (colors < 1)
    luaL_error(L, "qcd.hql expects .Colors >= 1");
  flavor_dim = q_get_intopt(L, 1, "Flavor", 1);
  spin_dim = q_get_intopt(L, 1, "Spin", 4);
  if (flavor_dim < 1)
    luaL_error(L, "qcd.hql expects .Flavor >= 1");
  if ((spin_dim != 1) && (spin_dim != 4))
    luaL_error(L, "qcd.hql expects .Spin to be 1 or 4");

  qlua_newHQLGrid(L, Sidx, flavor_dim, spin_dim, colors);
  return 1;
}

static const struct luaL_Reg fHQL[] = {
  { "hql", q_hql },
  { NULL, NULL }
};

int
init_hql(lua_State *L)
{
    static const QLUA_Op2 ops[] = {
      { qlua_mul_table, qHQLMatrix,  qHQLMatrix,  qhql_m_mul_m },
      { qlua_mul_table, qHQLMatrix,  qHQLVector,  qhql_m_mul_v },
      { qlua_mul_table, qHQLMatrix,  qComplex,    qhql_m_mul_c },
      { qlua_mul_table, qHQLMatrix,  qReal,       qhql_m_mul_r },
      { qlua_mul_table, qComplex,    qHQLMatrix,  qhql_c_mul_m },
      { qlua_mul_table, qReal,       qHQLMatrix,  qhql_r_mul_m },
      { qlua_div_table, qHQLMatrix,  qComplex,    qhql_m_div_c },
      { qlua_div_table, qHQLMatrix,  qReal,       qhql_m_div_r },
      { qlua_add_table, qHQLMatrix,  qHQLMatrix,  qhql_m_add_m },
      { qlua_sub_table, qHQLMatrix,  qHQLMatrix,  qhql_m_sub_m },
      { qlua_mul_table, qHQLVector,  qComplex,    qhql_v_mul_c },
      { qlua_mul_table, qHQLVector,  qReal,       qhql_v_mul_r },
      { qlua_mul_table, qComplex,    qHQLVector,  qhql_c_mul_v },
      { qlua_mul_table, qReal,       qHQLVector,  qhql_r_mul_v },
      { qlua_div_table, qHQLVector,  qComplex,    qhql_v_div_c },
      { qlua_div_table, qHQLVector,  qReal,       qhql_v_div_r },
      { qlua_add_table, qHQLVector,  qHQLVector,  qhql_v_add_v },
      { qlua_sub_table, qHQLVector,  qHQLVector,  qhql_v_sub_v },
      { NULL,           qNoType,     qNoType,     NULL         }
    };

  luaL_register(L, qcdlib, fHQL);
  qlua_reg_op2(ops);
  qlua_reg_dot(qHQLVector, qhql_v_dot_v);

  return 0;
}

int
fini_hql(lua_State *L)
{
  return 0;
}
