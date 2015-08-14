#include "modules.h"                                                 /* DEPS */
#include "qlua.h"                                                    /* DEPS */
#include "qqopqdp.h"                                                 /* DEPS */
#include <qop_qdp.h>

static const char qopqdp[] = "qop";

int qqq_inited = 0;

/* names and routines for qcd.qdpc table */
static const struct luaL_Reg fQOP[] = {
    { "WilsonMG",   qqq_wmg},
    { NULL,       NULL }
};

int
init_qopqdp(lua_State *L)
{
  lua_getglobal(L, qcdlib);
  lua_newtable(L);
  luaL_register(L, NULL, fQOP);
  lua_setfield(L, -2, qopqdp);
  lua_pop(L, 1);
    
  return 0;
}

void
fini_qopqdp(void)
{
  if (qqq_inited)
      QOP_finalize();
  qqq_inited = 0;
}
