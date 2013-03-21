#include "modules.h"                        /* DEPS */
#include "qlua.h"                           /* DEPS */
#include "qhql.h"                           /* DEPS */
#include "lattice.h"                        /* DEPS */

int
qhql_matrix(lua_State *L)
{
  mHQLGrid *grid = qlua_checkHQLGrid(L, 1);
  int stencil_size = lua_objlen(L, 2);
  mLattice *S = qlua_ObjLattice(L, 1);
  /* XXX */
  printf("XXX qhql_matrix: WRITE ME, |s| = %d, S=%p\n", stencil_size, S);
  lua_pushnil(L);
  return 1;
}
