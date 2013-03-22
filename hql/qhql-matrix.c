#include "modules.h"                        /* DEPS */
#include "qlua.h"                           /* DEPS */
#include "lattice.h"                        /* DEPS */
#include "qhql.h"                           /* DEPS */
#include "hql.h"
#include <string.h>

static const char *HQLMatrixName = "hql.matrix";

static int
q_m_fmt(lua_State *L)
{
  char fmt[72];
  
  qlua_checkHQLMatrix(L, 1, NULL, NULL);
  sprintf(fmt, "HQLMatrix()");
  lua_pushstring(L, fmt);
  return 1;
}

static int
q_m_gc(lua_State *L)
{
  mHQLMatrix *matrix = qlua_checkHQLMatrix(L, 1, NULL, NULL);

  HQL_OperatorFree(&matrix->op);
  return 0;
}

static int
q_m_get(lua_State *L)
{
  return qlua_selflookup(L, 1, luaL_checkstring(L, 2));
}

static struct luaL_Reg mtHQLMatrix[] = {
  { "__tostring",   q_m_fmt       },
  { "__gc",         q_m_gc        },
  { "__index",      q_m_get       },
  { "__newindex",   qlua_nowrite  },
  { "__add",        qlua_add      },
  { "__sub",        qlua_sub      },
  { "__mul",        qlua_mul      },
  { "__div",        qlua_div      },
  { "__unm",        qhql_m_neg    },
  { "vector",       qhql_vector   },
  { "apply",        qlua_mul      },
  /* "lattice" */
  /* "grid" */
  /* "a-type" */
  { NULL,            NULL}
};

int qhql_m_mul_m(lua_State *L) { /* XXX */ return 0; }
int qhql_m_mul_v(lua_State *L) { /* XXX */ return 0; }
int qhql_m_mul_c(lua_State *L) { /* XXX */ return 0; }
int qhql_m_mul_r(lua_State *L) { /* XXX */ return 0; }
int qhql_r_mul_m(lua_State *L) { /* XXX */ return 0; }
int qhql_c_mul_m(lua_State *L) { /* XXX */ return 0; }
int qhql_m_add_m(lua_State *L) { /* XXX */ return 0; }
int qhql_m_sub_m(lua_State *L) { /* XXX */ return 0; }
int qhql_m_div_c(lua_State *L) { /* XXX */ return 0; }
int qhql_m_div_r(lua_State *L) { /* XXX */ return 0; }
int qhql_m_neg(lua_State *L)   { /* XXX */ return 0; }

typedef struct {
  HQL_Stencil_t *stencil;
#if USE_Nc2
  QDP_D2_ColorMatrix *U2;
#endif
#if USE_Nc3
  QDP_D3_ColorMatrix *U3;
#endif
#if USE_NcN
  QDP_DN_ColorMatrix *UN;
  int nc;
#endif
} qStencil;

mHQLMatrix *
qlua_newHQLMatrix(lua_State *L, int Sidx, int Gidx, HQL_Operator_t *op, mHQLGrid *grid)
{
  mHQLMatrix *v = lua_newuserdata(L, sizeof (mHQLMatrix));

  v->grid = grid;
  v->op = op;

  qlua_createHQLTable(L, Sidx, Gidx, mtHQLMatrix, qHQLMatrix, HQLMatrixName);
  lua_setmetatable(L, -2);

  return v;
}

mHQLMatrix *
qlua_checkHQLMatrix(lua_State *L, int idx, mLattice *S, mHQLGrid *G)
{
  void *v = qlua_checkLatticeType(L, idx, qHQLMatrix, HQLMatrixName);
  
  if (S) {
    mLattice *S1 = qlua_ObjLattice(L, idx);
    if (S1->id != S->id)
      luaL_error(L, "%s on a wrong lattice", HQLMatrixName);
    lua_pop(L, 1);
  }
  if (G) {
    mHQLGrid *G1 = qlua_ObjGrid(L, idx);
    if (G != G1)
      luaL_error(L, "%s on a wrong grid", HQLMatrixName);
    lua_pop(L, 1);
  }

  return (mHQLMatrix *)v;
}

static void
check_status(HQL_err_t status, lua_State *L, const char *msg)
{
  if (status != 0)
    luaL_error(L, "hql.matrix failed to %s: %s", msg, HQL_error_msg(status));
}

static void
collect_stencil(mHQLGrid *grid,
                qStencil *stencil,
                int idx,
                lua_State *L,
                int *offset,
                mLattice *S)
{
  memset(stencil, 0, sizeof (qStencil));
  /* XXX collect_stencil */
}

static void
setup_matrix(HQL_Operator_t *ho,
             qStencil *stencil,
             lua_State *L,
             int *offset,
             mLattice *S,
             mHQLGrid *grid)
{
  /* XXX setup_matrix */
}

int
qhql_matrix(lua_State *L)
{
  int stencil_size = lua_objlen(L, 2);
  mLattice *S = qlua_ObjLattice(L, 1);
  int Sidx = lua_gettop(L);
  mHQLGrid *grid = qlua_copyHQLGrid(L, 1);
  int Gidx = lua_gettop(L);
  MPI_Comm comm = qlua_latticeMPI(S);
  HQL_Operator_t *ho = NULL;
  int *offset = qlua_malloc(L, S->rank * sizeof (int));
  qStencil *stencil = qlua_malloc(L, stencil_size * sizeof (qStencil));
  int i;

  CALL_QDP(L);
  check_status(HQL_GridCreate(&grid->grid,
                              comm,
                              S->rank,
                              grid->flavor_dim,
                              grid->spin_dim,
                              grid->colors,
                              S->dim,
                              S->box_low,
                              S->box_high),
               L, "create grid");
  for (i = 0; i < stencil_size; i++)
    collect_stencil(grid, &stencil[i], i, L, offset, S);
  check_status(HQL_GridAssemble(grid->grid), L, "assemble grid");
  check_status(HQL_OperatorCreate(&ho, grid->grid), L, "create matrix");
  for (i = 0; i < stencil_size; i++)
    setup_matrix(ho, &stencil[i], L, offset, S, grid);
  check_status(HQL_OperatorAssemble(ho), L, "assemble matrix");
  qlua_free(L, offset);
  qlua_free(L, stencil);
  qlua_newHQLMatrix(L, Sidx, Gidx, ho, grid);

  return 1;
}
