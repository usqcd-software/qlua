#include "modules.h"                        /* DEPS */
#include "qlua.h"                           /* DEPS */
#include "lattice.h"                        /* DEPS */
#include "qmatrix.h"                        /* DEPS */
#include "qgamma.h"                         /* DEPS */
#include "latcolmat.h"                      /* DEPS */
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
  { "apply",        qhql_m_mul_v  },
  /* "lattice" */
  /* "grid" */
  /* "a-type" */
  { NULL,            NULL}
};

int
qhql_m_mul_v(lua_State *L)
{
  mHQLMatrix *m = qlua_checkHQLMatrix(L, 1, NULL, NULL);
  mLattice *S = qlua_ObjLattice(L, 1);
  int Sidx = lua_gettop(L);
  mHQLGrid *G = qlua_ObjGrid(L, 1);
  int Gidx = lua_gettop(L);
  mHQLVector *v = qlua_checkHQLVector(L, 2, S, G);
  mHQLVector *r = qlua_newHQLVector(L, Sidx, Gidx, NULL, G);

  qhql_check_status(HQL_OperatorApply(r->vec, m->op, v->vec),
                    L, "apply operator");

  return 1;
}

/* other operations of HQL matrices to be implemented later */
int qhql_m_mul_m(lua_State *L) { luaL_error(L, "/* YYY qhql_m_mul_m() */\n"); return 0; }
int qhql_m_mul_c(lua_State *L) { luaL_error(L, "/* YYY qhql_m_mul_c() */\n"); return 0; }
int qhql_m_mul_r(lua_State *L) { luaL_error(L, "/* YYY qhql_m_mul_r() */\n"); return 0; }
int qhql_r_mul_m(lua_State *L) { luaL_error(L, "/* YYY qhql_r_mul_m() */\n"); return 0; }
int qhql_c_mul_m(lua_State *L) { luaL_error(L, "/* YYY qhql_c_mul_m() */\n"); return 0; }
int qhql_m_add_m(lua_State *L) { luaL_error(L, "/* YYY qhql_m_add_m() */\n"); return 0; }
int qhql_m_sub_m(lua_State *L) { luaL_error(L, "/* YYY qhql_m_sub_m() */\n"); return 0; }
int qhql_m_div_c(lua_State *L) { luaL_error(L, "/* YYY qhql_m_div_c() */\n"); return 0; }
int qhql_m_div_r(lua_State *L) { luaL_error(L, "/* YYY qhql_m_div_r() */\n"); return 0; }
int qhql_m_neg(lua_State *L)   { luaL_error(L, "/* YYY qhql_m_neg() */\n"); return 0; }

typedef struct {
  HQL_Stencil_t *stencil;
#if USE_Nc2
  mLatColMat2 *U2;
#endif
#if USE_Nc3
  mLatColMat3 *U3;
#endif
#if USE_NcN
  mLatColMatN *UN;
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

void
qhql_check_status(HQL_err_t status, lua_State *L, const char *msg)
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
                QLA_D_Complex *flavor,
                mLattice *S)
{
  int i;

  memset(stencil, 0, sizeof (qStencil));
  /* default offset is [0,...] */
  memset(offset, 0, S->rank * sizeof (int));
  lua_pushnumber(L, idx + 1);
  lua_gettable(L, 2);
  /* if .offset is defined, load it */
  lua_getfield(L, -1, "offset");
  if (!lua_isnil(L, -1)) {
    for (i = 0; i < S->rank; i++) {
      lua_pushnumber(L, i + 1);
      lua_gettable(L, -2);
      offset[i] = luaL_checkint(L, -1);
      lua_pop(L, 1);
    }
  }
  lua_pop(L, 1);
  qhql_check_status(HQL_StencilCreate(&stencil->stencil, grid->grid, offset),
                    L, "create stencil");
  /* if .gamma is defined, load it to the stencil */
  lua_getfield(L, -1, "gamma");
  if (!lua_isnil(L, -1)) {
    int a, b;
    mMatComplex *m = gamma2matrix(L, -1);
    QLA_D_Complex g[QDP_Ns * QDP_Ns];
    if (grid->spin_dim != QDP_Ns)
      luaL_error(L, "no gamma structure defined in the grid");
    for (a = 0; a < QDP_Ns; a++) {
      for (b = 0; b < QDP_Ns; b++) {
        gsl_complex zz = gsl_matrix_complex_get(m->m, a, b);
        QLA_real(g[a * QDP_Ns + b]) = GSL_REAL(zz);
        QLA_imag(g[a * QDP_Ns + b]) = GSL_IMAG(zz);
      }
    }
    qhql_check_status(HQL_StencilSetSpin(stencil->stencil, g),
                      L, "set gamma");
    lua_pop(L, 1);
  }
  lua_pop(L, 1);
  /* if .flavor is defined, load it */
  lua_getfield(L, -1, "flavor");
  if (!lua_isnil(L, -1)) {
    int a, b;
    mMatComplex *m = qlua_checkMatComplex(L, -1);
    if ((m->l_size != grid->flavor_dim) || (m->r_size != grid->flavor_dim))
      luaL_error(L, "flavor matrix shape mismatch");
    for (a = 0; a < grid->flavor_dim; a++) {
      for (b = 0; b < grid->flavor_dim; b++) {
        gsl_complex zz = gsl_matrix_complex_get(m->m, a, b);
        QLA_real(flavor[a * grid->flavor_dim + b]) = GSL_REAL(zz);
        QLA_imag(flavor[a * grid->flavor_dim + b]) = GSL_IMAG(zz);
      }
    }
    qhql_check_status(HQL_StencilSetFlavor(stencil->stencil, flavor),
                      L, "set flavor");
  }
  lua_pop(L, 1);
  /* if .U is defined, store it into appropriate part of stencil for now */
  lua_getfield(L, -1, "U");
  if (!lua_isnil(L, -1)) {
    switch (qlua_qtype(L, -1)) {
#if USE_Nc2
    case qLatColMat2:
      stencil->U2 = qlua_checkLatColMat2(L, -1, S, 2);
      if (grid->colors != 2)
        luaL_error(L, "Nc mismatch (%d vs %d)", 2, grid->colors);
      break;
#endif
#if USE_Nc3
    case qLatColMat3:
      stencil->U3 = qlua_checkLatColMat3(L, -1, S, 3);
      if (grid->colors != 3)
        luaL_error(L, "Nc mismatch (%d vs %d)", 3, grid->colors);
      break;
#endif
#if USE_NcN
    case qLatColMatN:
      stencil->UN = qlua_checkLatColMatN(L, -1, S, grid->colors);
      if (grid->colors != stencil->UN->nc)
        luaL_error(L, "Nc mismatch (%d vs %d)", stencil->UN->nc, grid->colors);
      break;
#endif
    default:
      luaL_error(L, "stencil[%d].U of unexpected type", idx);
      break;
    }
  }
  lua_pop(L, 2);
}

static void
setup_matrix(HQL_Operator_t *ho,
             qStencil *stencil,
             lua_State *L,
             mLattice *S,
             mHQLGrid *grid)
{
#if USE_Nc2
  if (stencil->U2) {
    qhql_check_status(HQL_OperatorSetColorMatrix2(ho,
                                                  stencil->stencil,
                                                  stencil->U2->ptr),
                 L, "set color Nc=2");
    return;
  }
#endif
#if USE_Nc3
  if (stencil->U3) {
    qhql_check_status(HQL_OperatorSetColorMatrix3(ho,
                                                  stencil->stencil,
                                                  stencil->U3->ptr),
                 L, "set color Nc=3");
    return;
  }
#endif
#if USE_NcN
  if (stencil->UN) {
    qhql_check_status(HQL_OperatorSetColorMatrixN(ho,
                                                  stencil->stencil,
                                                  stencil->UN->ptr),
                 L, "set color Nc=N");
    return;
  }
#endif
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
  QLA_D_Complex *flavor = qlua_malloc(L,
                                      grid->flavor_dim
                                      * grid->flavor_dim
                                      * sizeof (QLA_D_Complex));
  int i;

  CALL_QDP(L);
  qhql_check_status(HQL_GridCreate(&grid->grid,
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
    collect_stencil(grid, &stencil[i], i, L, offset, flavor, S);
  qhql_check_status(HQL_GridAssemble(grid->grid), L, "assemble grid");
  qhql_check_status(HQL_OperatorCreate(&ho, grid->grid), L, "create matrix");
  for (i = 0; i < stencil_size; i++)
    setup_matrix(ho, &stencil[i], L, S, grid);
  qhql_check_status(HQL_OperatorAssemble(ho), L, "assemble matrix");
  qlua_free(L, flavor);
  qlua_free(L, offset);
  qlua_free(L, stencil);
  qlua_newHQLMatrix(L, Sidx, Gidx, ho, grid);

  return 1;
}
