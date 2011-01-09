#define QOP_Precision 2
#include "modules.h"                                                 /* DEPS */
#include "qlua.h"                                                    /* DEPS */
#include "qclover.h"                                                 /* DEPS */
#include "qcomplex.h"                                                /* DEPS */
#include "qvector.h"                                                 /* DEPS */
#include "lattice.h"                                                 /* DEPS */
#include "qlayout.h"                                                 /* DEPS */
#include "latcolvec.h"                                               /* DEPS */
#include "latcolmat.h"                                               /* DEPS */
#include "latdirferm.h"                                              /* DEPS */
#include "latdirprop.h"                                              /* DEPS */
#include "qop_qdp.h"
#include "qdp_f3.h"
#include "qdp_df3.h"

#define printf0 if(QDP_this_node==0) printf

#if USE_Nc3

static void
Init(QDP_Lattice *lat)
{
  static int inited=0;

  if(!QOP_is_initialized()) {
    QOP_layout_t qoplayout;
    qoplayout.latdim = QDP_ndim_L(lat);
    qoplayout.latsize = (int *) malloc(QDP_ndim_L(lat)*sizeof(int));
    qoplayout.machdim = -1;
    QDP_latsize_L(lat, qoplayout.latsize);
    QDP_set_default_lattice(lat);
    printf0("begin QOP init... ");
    QOP_init(&qoplayout);
    printf0("done\n");
  }
  inited = 1;
}

static void
add_default(lua_State *L, const char *key, double value)
{
  lua_pushnumber(L, value);
  lua_setfield(L, -2, key);
}

static double
set_default(lua_State *L, const char *key, double def)
{
  double v = def;

  lua_getfield(L, -1, key);
  if (qlua_qtype(L, -1) == qReal) {
    v = luaL_checknumber(L, -1);
  }
  lua_pop(L, 1);
  return v;
}

#include "qqopqdp-asqtad.c"
#include "qqopqdp-clover.c"
#include "qqopqdp-hisq.c"

static struct luaL_Reg fQOPQDP[] = {
    { "Asqtad",       q_asqtad },
    { "Clover",       q_clover },
    { "Hisq",         q_hisq },
    { NULL,           NULL }
};

int
init_qopqdp(lua_State *L)
{
  lua_getglobal(L, qcdlib);
  lua_newtable(L);
  luaL_register(L, NULL, fQOPQDP);
  lua_setfield(L, -2, "qopqdp");
  lua_pop(L, 1);
  return 0;
}
#else /* USE_Nc3 */
int
init_qopqdp(lua_State *L)
{
  return 0;
}
#endif /* USE_Nc3 */

int
fini_qopqdp(lua_State *L)
{
  return 0;
}
