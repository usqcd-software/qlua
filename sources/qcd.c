#include "modules.h"                                                 /* DEPS */
#include "qlua.h"                                                    /* DEPS */
#include "qio.h"
#include "qdp.h"
#ifdef HAS_QOPQDP
#include "qop.h"
#endif

static void
add_default(lua_State *L, const char *key, int value)
{
  lua_pushnumber(L, value);
  lua_setfield(L, -2, key);
}

static int
set_default(lua_State *L, const char *key, int def)
{
  int v = def;

  lua_getfield(L, -1, key);
  if (qlua_qtype(L, -1) == qReal) {
    v = luaL_checkint(L, -1);
  }
  lua_pop(L, 1);
  return v;
}

static int
q_defaults(lua_State *L)
{
  struct keystruct {
    const char *key;
    int (*func)(int);
  };
  static struct keystruct keys[] = {
    { "qioVerbose", QIO_verbose },
    { "qdpProfcontrol", QDP_profcontrol },
#ifdef HAS_QOPQDP
    { "qopqdpProfcontrol", QOP_profcontrol },
    { "qopqdpVerbose", QOP_verbose },
#endif
    { "", NULL }
  };
  int nkeys = sizeof(keys)/sizeof(struct keystruct) - 1;
  //printf("nkeys = %i\n", nkeys);

  switch (lua_gettop(L)) {
  case 0: {
    lua_createtable(L, 0, nkeys);
    for(int i=0; i<nkeys; i++) {
      int val = keys[i].func(0);
      keys[i].func(val);
      add_default(L, keys[i].key, val);
    }
    return 1;
  }
  case 1: {
    if (qlua_qtype(L, 1) == qTable) {
      for(int i=0; i<nkeys; i++) {
	int val = keys[i].func(0);
	val = set_default(L, keys[i].key, val);
	keys[i].func(val);
      }
      return 0;
    }
    break;
  }
  default:
    break;
  }
  return luaL_error(L, "bad parameters for qcd.defaults()");
}

static struct luaL_Reg fQcd[] = {
  { "defaults",           q_defaults },
  { NULL, NULL}
};

int
init_qcd(lua_State *L)
{
  luaL_register(L, qcdlib, fQcd);
  lua_pop(L, 1);
  return 0;
}

int
fini_qcd(lua_State *L)
{
  return 0;
}
