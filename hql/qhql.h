#ifndef MARK_22D544B0_05A3_46A9_910F_CA4C1194932B
#define MARK_22D544B0_05A3_46A9_910F_CA4C1194932B
#include "hql.h"

typedef struct {
  int         flavor_dim;
  int         spin_dim;
  int         colors;
  HQL_Grid_t *grid;
  /* Lattice is in the metatable */
} mHQLGrid;

typedef struct {
  mHQLGrid        *grid; /* duplicate of the .Grid from the metatable */
  HQL_Operator_t  *op;
  /* Lattice and Grid are in the metatable */
} mHQLMatrix;

typedef struct {
  mHQLGrid        *grid; /* duplicate of the .Grid from the metatable */
  HQL_Operator_t  *op;
  /* Lattice and Grid are in the metatable */
} mHQLVector;

mHQLGrid *qlua_ObjGrid(lua_State *L, int idx);

void qlua_createHQLTable(lua_State *L,
                         int Sidx,
                         int Gidx,
                         const struct luaL_Reg *ft,
                         QLUA_Type t_id,
                         const char *name);
mHQLGrid *qlua_newHQLGrid(lua_State *L,
                          int S_idx,
                          int flavor_dim,
                          int spin_dim,
                          int colors);
mHQLGrid *qlua_copyHQLGrid(lua_State *L,int idx);
mHQLGrid *qlua_checkHQLGridTemplate(lua_State *L, int idx);
mHQLGrid *qlua_checkHQLGrid(lua_State *L, int idx);

mHQLMatrix *qlua_newHQLMatrix(lua_State *L,
                              int Sidx,
                              int Gidx,
                              HQL_Operator_t *op,
                              mHQLGrid *grid);
mHQLMatrix *qlua_checkHQLMatrix(lua_State *L,
                                int idx,
                                mLattice *S,
                                mHQLGrid *G);

mHQLVector *qlua_newHQLVector(lua_State *L,
                              int Sidx,
                              int Gidx,
                              HQL_Vector_t *vec,
                              mHQLGrid *grid);
mHQLVector *qlua_checkHQLVector(lua_State *L,
                                int idx,
                                mLattice *S,
                                mHQLGrid *G);

int qhql_matrix(lua_State *L);
int qhql_vector(lua_State *L);

int qhql_m_mul_m(lua_State *L);
int qhql_m_mul_v(lua_State *L);
int qhql_m_mul_c(lua_State *L);
int qhql_m_mul_r(lua_State *L);
int qhql_r_mul_m(lua_State *L);
int qhql_c_mul_m(lua_State *L);
int qhql_m_add_m(lua_State *L);
int qhql_m_sub_m(lua_State *L);
int qhql_m_div_c(lua_State *L);
int qhql_m_div_r(lua_State *L);
int qhql_m_neg(lua_State *L);

int qhql_v_mul_c(lua_State *L);
int qhql_v_mul_r(lua_State *L);
int qhql_c_mul_v(lua_State *L);
int qhql_r_mul_v(lua_State *L);
int qhql_v_div_c(lua_State *L);
int qhql_v_div_r(lua_State *L);
int qhql_v_add_v(lua_State *L);
int qhql_v_sub_v(lua_State *L);
int qhql_v_dot_v(lua_State *L);
int qhql_v_norm2(lua_State *L);
int qhql_v_neg(lua_State *L);

int init_hql(lua_State *L);
int fini_hql(lua_State *L);

#endif /* !defined(MARK_22D544B0_05A3_46A9_910F_CA4C1194932B) */
