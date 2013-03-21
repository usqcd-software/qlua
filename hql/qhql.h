#ifndef MARK_22D544B0_05A3_46A9_910F_CA4C1194932B
#define MARK_22D544B0_05A3_46A9_910F_CA4C1194932B

typedef struct {
  int flavor_dim;
  int spin_dim;
  /* Lattice is in the metatable */
} mHQLGrid;

typedef struct {
  /* XXX */
} mHQLMatrix;

typedef struct {
  /* XXX */
} mHQLVector;

mHQLGrid *qlua_newHQLGrid(lua_State *L, int S_idx, int flavor_dim, int spin_dim);
mHQLGrid *qlua_checkHQLGrid(lua_State *L, int idx);
mHQLMatrix *qlua_newHQLMatrix(lua_State *L); /* XXX */
mHQLMatrix *qlua_checkHQLMatrix(lua_State *L); /* XXX */
mHQLVector *qlua_newHQLVector(lua_State *L); /* XXX */
mHQLVector *qlua_checkHQLVector(lua_State *L); /* XXX */

int init_hql(lua_State *L);
int fini_hql(lua_State *L);

int qhql_matrix(lua_State *L);
int qhql_vector(lua_State *L);

#endif /* !defined(MARK_22D544B0_05A3_46A9_910F_CA4C1194932B) */
