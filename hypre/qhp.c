#include "qlua.h"                           /* DEPS */
#include "qhp.h"                            /* DEPS */
#include "qhp-i.h"                          /* DEPS */

int
init_qhp(lua_State *L)
{
  /* XXX */
  qhp_init_mv(L);
  qhp_init_solvers(L);
  return 0;
}

int
fini_qhp(lua_State *L)
{
  qhp_fini_solvers(L);
  qhp_fini_mv(L);
  /* XXX */
  return 0;
}
