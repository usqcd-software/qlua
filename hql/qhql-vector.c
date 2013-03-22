#include "modules.h"                        /* DEPS */
#include "qlua.h"                           /* DEPS */
#include "lattice.h"                        /* DEPS */
#include "latcolvec.h"                      /* DEPS */
#include "latdirferm.h"                     /* DEPS */
#include "qcomplex.h"                       /* DEPS */
#include "qhql.h"                           /* DEPS */

static const char *HQLVectorName = "hql.vector";

static int
q_v_fmt(lua_State *L)
{
  char fmt[72];
  
  qlua_checkHQLVector(L, 1, NULL, NULL);
  sprintf(fmt, "HQLVector()");
  lua_pushstring(L, fmt);
  return 1;
}

static int
q_v_gc(lua_State *L)
{
  mHQLVector *matrix = qlua_checkHQLVector(L, 1, NULL, NULL);

  HQL_VectorFree(&matrix->vec);
  return 0;
}

static int
q_v_get(lua_State *L)
{
  return qlua_selflookup(L, 1, luaL_checkstring(L, 2));
}

static struct luaL_Reg mtHQLVector[] = {
  { "__tostring",   q_v_fmt       },
  { "__gc",         q_v_gc        },
  { "__index",      q_v_get       },
  { "__newindex",   qlua_nowrite  },
  { "__add",        qlua_add      },
  { "__sub",        qlua_sub      },
  { "__mul",        qlua_mul      },
  { "__div",        qlua_div      },
  { "__unm",        qhql_m_neg    },
  { "export",       qhql_v_export },
  /* "lattice" */
  /* "grid" */
  /* "a-type" */
  { NULL,            NULL}
};

static int
v_scale(lua_State *L, int Vidx, QLA_D_Complex *c)
{
  mHQLVector *v = qlua_checkHQLVector(L, Vidx, NULL, NULL);
  qlua_ObjLattice(L, Vidx);
  {
    int Sidx = lua_gettop(L);
    mHQLGrid *G = qlua_ObjGrid(L, Vidx);
    int Gidx = lua_gettop(L);
    mHQLVector *r = qlua_newHQLVector(L, Sidx, Gidx, NULL, G);
    
    CALL_QDP(L);
    qhql_check_status(HQL_VectorScale(r->vec, c, v->vec),
                      L, "scale vector");
  }
  return 1;
}

int
qhql_v_mul_c(lua_State *L)
{
  QLA_D_Complex *c = qlua_checkComplex(L, 2);
  return v_scale(L, 1, c);
}

int
qhql_v_mul_r(lua_State *L)
{
  double r = luaL_checknumber(L, 2);
  QLA_D_Complex c;
  QLA_real(c) = r;
  QLA_imag(c) = 0.0;
  return v_scale(L, 1, &c);
   
}

int
qhql_c_mul_v(lua_State *L)
{
  QLA_D_Complex *c = qlua_checkComplex(L, 1);
  return v_scale(L, 2, c);
}

int
qhql_r_mul_v(lua_State *L)
{
  double r = luaL_checknumber(L, 1);
  QLA_D_Complex c;
  QLA_real(c) = r;
  QLA_imag(c) = 0.0;
  return v_scale(L, 2, &c);
}


int
qhql_v_div_c(lua_State *L)
{
  QLA_D_Complex *c = qlua_checkComplex(L, 2);
  double h = 1/hypot(QLA_real(*c), QLA_imag(*c));
  QLA_D_Complex v;
  QLA_real(v) = (QLA_real(*c) * h) * h;
  QLA_imag(v) = -(QLA_imag(*c) * h) * h;
  return v_scale(L, 1, &v);
}

int
qhql_v_div_r(lua_State *L)
{
  double r = luaL_checknumber(L, 2);
  QLA_D_Complex c;
  QLA_real(c) = 1/r;
  QLA_imag(c) = 0.0;
  return v_scale(L, 1, &c);
}

int
qhql_v_neg(lua_State *L)
{
  QLA_D_Complex c;
  QLA_real(c) = -1;
  QLA_imag(c) = 0.0;
  return v_scale(L, 1, &c);
}

int
qhql_v_add_v(lua_State *L)
{
  mHQLVector *a = qlua_checkHQLVector(L, 1, NULL, NULL);
  mLattice *S = qlua_ObjLattice(L, 1);
  int Sidx = lua_gettop(L);
  mHQLGrid *G = qlua_ObjGrid(L, 1);
  int Gidx = lua_gettop(L);
  mHQLVector *b = qlua_checkHQLVector(L, 2, S, G);
  mHQLVector *r = qlua_newHQLVector(L, Sidx, Gidx, NULL, G);
    
  CALL_QDP(L);
  qhql_check_status(HQL_VectorAdd(r->vec, a->vec, b->vec),
                    L, "add vectors");
  return 1;
}

int
qhql_v_sub_v(lua_State *L)
{
  mHQLVector *a = qlua_checkHQLVector(L, 1, NULL, NULL);
  mLattice *S = qlua_ObjLattice(L, 1);
  int Sidx = lua_gettop(L);
  mHQLGrid *G = qlua_ObjGrid(L, 1);
  int Gidx = lua_gettop(L);
  mHQLVector *b = qlua_checkHQLVector(L, 2, S, G);
  mHQLVector *r = qlua_newHQLVector(L, Sidx, Gidx, NULL, G);
    
  CALL_QDP(L);
  qhql_check_status(HQL_VectorSub(r->vec, a->vec, b->vec),
                    L, "sub vectors");
  return 1;
}

int
qhql_v_dot_v(lua_State *L)
{
  mHQLVector *a = qlua_checkHQLVector(L, 1, NULL, NULL);
  mLattice *S = qlua_ObjLattice(L, 1);
  mHQLGrid *G = qlua_ObjGrid(L, 1);
  mHQLVector *b = qlua_checkHQLVector(L, 2, S, G);
  QLA_D_Complex *c = qlua_newComplex(L);
    
  CALL_QDP(L);
  qhql_check_status(HQL_VectorDot(c, a->vec, b->vec),
                    L, "dot vectors");
  return 1;
}

int
qhql_v_norm2(lua_State *L)
{
  mHQLVector *a = qlua_checkHQLVector(L, 1, NULL, NULL);
  double n = 0.0;
    
  CALL_QDP(L);
  qhql_check_status(HQL_VectorNorm2(&n, a->vec),
                    L, "norm2 vector");
  lua_pushnumber(L, n);
  return 1;
}

static void
vec_export(lua_State *L,
           int fidx,
           mHQLVector *vec,
           mLattice *S,
           int Sidx,
           mHQLGrid *G)
{
  switch (G->spin_dim) {
  case 1: {
#if USE_Nc2
    if (G->colors == 2) {
      mLatColVec2 *r = qlua_newZeroLatColVec2(L, Sidx, 2);
      CALL_QDP(L);
      qhql_check_status(HQL_VectorReadColorVector2(r->ptr, vec->vec, fidx),
                        L, "export Nc=2");
      return;
    }
#endif
#if USE_Nc3
    if (G->colors == 3) {
      mLatColVec3 *r = qlua_newZeroLatColVec3(L, Sidx, 3);
      CALL_QDP(L);
      qhql_check_status(HQL_VectorReadColorVector3(r->ptr, vec->vec, fidx),
                        L, "export Nc=3");
      return;
    }
#endif
#if USE_NcN
    {
      mLatColVecN *r = qlua_newZeroLatColVecN(L, Sidx, G->colors);
      CALL_QDP(L);
      qhql_check_status(HQL_VectorReadColorVectorN(r->ptr, vec->vec, fidx),
                        L, "export Nc=N");
      return;
    }
#endif
    luaL_error(L, "Unexpected Nc = %d", G->colors);
  } break;
  case QDP_Ns: {
#if USE_Nc2
    if (G->colors == 2) {
      mLatDirFerm2 *r = qlua_newZeroLatDirFerm2(L, Sidx, 2);
      CALL_QDP(L);
      qhql_check_status(HQL_VectorReadDiracFermion2(r->ptr, vec->vec, fidx),
                        L, "export Nc=2");
      return;
    }
#endif
#if USE_Nc3
    if (G->colors == 3) {
      mLatDirFerm3 *r = qlua_newZeroLatDirFerm3(L, Sidx, 3);
      CALL_QDP(L);
      qhql_check_status(HQL_VectorReadDiracFermion3(r->ptr, vec->vec, fidx),
                        L, "export Nc=3");
      return;
    }
#endif
#if USE_NcN
    {
      mLatDirFermN *r = qlua_newZeroLatDirFermN(L, Sidx, G->colors);
      CALL_QDP(L);
      qhql_check_status(HQL_VectorReadDiracFermionN(r->ptr, vec->vec, fidx),
                        L, "export Nc=N");
      return;
    }
#endif
    luaL_error(L, "Unexpected Nc = %d", G->colors);
  } break;
  default:
    luaL_error(L, "Unexpected spin dimension %d", G->spin_dim);
    break;
  }
}

int
qhql_v_export(lua_State *L)
{
  mHQLVector *a = qlua_checkHQLVector(L, 1, NULL, NULL);
  mLattice *S = qlua_ObjLattice(L, 1);
  int Sidx = lua_gettop(L);
  mHQLGrid *G = qlua_ObjGrid(L, 1);

  if (G->flavor_dim == 1) {
    vec_export(L, 0, a, S, Sidx, G);
  } else {
    int i;
    lua_createtable(L, G->flavor_dim, 0);
    for (i = 0; i < G->flavor_dim; i++) {
      vec_export(L, i, a, S, Sidx, G);
      lua_pushnumber(L, i+1);
      lua_settable(L, -3);
    }
  }
  return 1;
}

mHQLVector *
qlua_newHQLVector(lua_State *L, int Sidx, int Gidx, HQL_Vector_t *vec, mHQLGrid *grid)
{
  mHQLVector *v = lua_newuserdata(L, sizeof (mHQLVector));

  v->grid = grid;
  if (vec == NULL) {
    CALL_QDP(L);
    qhql_check_status(HQL_VectorCreate(&vec, grid->grid),
                      L, "create vector");
    qhql_check_status(HQL_VectorAssemble(vec),
                      L, "assemble vector");
  }
  v->vec = vec;

  qlua_createHQLTable(L, Sidx, Gidx, mtHQLVector, qHQLVector, HQLVectorName);
  lua_setmetatable(L, -2);

  return v;
}

mHQLVector *
qlua_checkHQLVector(lua_State *L, int idx, mLattice *S, mHQLGrid *G)
{
  void *v = qlua_checkLatticeType(L, idx, qHQLVector, HQLVectorName);
  
  if (S) {
    mLattice *S1 = qlua_ObjLattice(L, idx);
    if (S1->id != S->id)
      luaL_error(L, "%s on a wrong lattice", HQLVectorName);
    lua_pop(L, 1);
  }
  if (G) {
    mHQLGrid *G1 = qlua_ObjGrid(L, idx);
    if (G != G1)
      luaL_error(L, "%s on a wrong grid", HQLVectorName);
    lua_pop(L, 1);
  }

  return (mHQLVector *)v;
}

static void
setup_vector(HQL_Vector_t *vx,
             int fidx,
             int vidx,
             lua_State *L,
             mLattice *S,
             mHQLGrid *G)
{
  switch (qlua_qtype(L, vidx)) {
#if USE_Nc2
  case qLatColVec2: {
    mLatColVec2 *v = qlua_checkLatColVec2(L, vidx, S, 2);
    qhql_check_status(HQL_VectorSetColorVector2(vx, fidx, v->ptr),
                      L, "set vector (color vector) Nc=2");
  } break;
  case qLatDirFerm2: {
    mLatDirFerm2 *v = qlua_checkLatDirFerm2(L, vidx, S, 2);
    qhql_check_status(HQL_VectorSetDiracFermion2(vx, fidx, v->ptr),
                      L, "set vector (dirac fermion) Nc=2");
  } break;
#endif
#if USE_Nc3
  case qLatColVec3: {
    mLatColVec3 *v = qlua_checkLatColVec3(L, vidx, S, 3);
    qhql_check_status(HQL_VectorSetColorVector3(vx, fidx, v->ptr),
                      L, "set vector (color vector) Nc=3");
  } break;
  case qLatDirFerm3: {
    mLatDirFerm3 *v = qlua_checkLatDirFerm3(L, vidx, S, 3);
    qhql_check_status(HQL_VectorSetDiracFermion3(vx, fidx, v->ptr),
                      L, "set vector (dirac fermion) Nc=3");
  } break;
#endif
#if USE_NcN
  case qLatColVecN: {
    mLatColVecN *v = qlua_checkLatColVecN(L, vidx, S, G->colors);
    qhql_check_status(HQL_VectorSetColorVectorN(vx, fidx, v->ptr),
                      L, "set vector (color vector) Nc=N");
  } break;
  case qLatDirFermN: {
    mLatDirFermN *v = qlua_checkLatDirFermN(L, vidx, S, G->colors);
    qhql_check_status(HQL_VectorSetDiracFermionN(vx, fidx, v->ptr),
                      L, "set vector (dirac fermion) Nc=N");
  } break;
#endif
  default:
    luaL_error(L, "qcd vector[%d] of unexpected type", fidx);
    break;
  }
}

int
qhql_vector(lua_State *L)
{
  qlua_checkHQLMatrix(L, 1, NULL, NULL);
  {
    mLattice *S = qlua_ObjLattice(L, 1);
    int Sidx = lua_gettop(L);
    mHQLGrid *G = qlua_ObjGrid(L, 1);
    int Gidx = lua_gettop(L);
    HQL_Vector_t *vx = NULL;
    
    CALL_QDP(L);
    qhql_check_status(HQL_VectorCreate(&vx, G->grid),
                      L, "create vector");
    if (qlua_qtype(L, 2) == qTable) {
      int i;
      for (i = 0; i < G->flavor_dim; i++) {
        lua_pushnumber(L, i + 1);
        lua_gettable(L, 2);
        setup_vector(vx, i, -1, L, S, G);
        lua_pop(L, 1);
      }
    } else {
      setup_vector(vx, 0, 2, L, S, G);
    }
    qhql_check_status(HQL_VectorAssemble(vx),
                      L, "assemble vector");
    qlua_newHQLVector(L, Sidx, Gidx, vx, G);
  }
  return 1;
}
