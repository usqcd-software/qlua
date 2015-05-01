#ifdef _OPENMP
#include "modules.h"                                           /* DEPS */
#include "qlua.h"                                              /* DEPS */
#include "qomp.h"                                              /* DEPS */
#include <string.h>
#include <omp.h>

static const char qomp_name[] = "openmp";

/* OpenMP controls */

/*    ICV                     ENV                         C access
 *
 *    dynamic                 OMP_DYNAMIC                 bool         r/w *_dynamic()
 *    nested   	              OMP_NESTED                  bool         r/w *_nested()
 *    num_threads	      OMP_NUM_THREADS             int          r/w *_{max,num}_threads()
 *    schedule		      OMP_SCHEDULE                string,int   r/w *_schedule()
 *    bind		      OMP_PROC_BIND               string       r/  *_proc_bind()
 *    thread+limit	      OMP_THREAD_LIMIT            int          r/  *_thread_limit()
 *    max_active_levels	      OMP_MAX_ACTIVE_LEVELS       int          r/w *_max_active_levels()
 *    default-device	      OMP_DEFAULT_DEVICE          int          r/w *_default_device()
 *    num_devices             -                           int          r/  *_num_devices()
 *    version                 -                           string       r/  _OPENMP
 *    
 */

static int
qomp_dynamic(lua_State *L)
{
  switch (lua_gettop(L)) {
  case 0: {
    int n = omp_get_dynamic();
    lua_pushboolean(L, n);
    return 1;
  };
  case 1: {
    int v = lua_toboolean(L, 1);
    omp_set_dynamic(v);
  } break;
  default:
    luaL_error(L, "illegal arguments");
    break;
  }
  return 0;
}

static int
qomp_nested(lua_State *L)
{
  switch (lua_gettop(L)) {
  case 0: {
    int n = omp_get_nested();
    lua_pushboolean(L, n);
    return 1;
  };
  case 1: {
    int v = lua_toboolean(L, 1);
    omp_set_nested(v);
  } break;
  default:
    luaL_error(L, "illegal arguments");
    break;
  }
  return 0;
}

static int
qomp_num_threads(lua_State *L)
{
  switch (lua_gettop(L)) {
  case 0: {
    int n = omp_get_max_threads();
    lua_pushinteger(L, n);
    return 1;
  };
  case 1: {
    int v = lua_tointeger(L, 1);
    omp_set_num_threads(v);
  } break;
  default:
    luaL_error(L, "illegal arguments");
    break;
  }
  return 0;
}

static int
qomp_levels(lua_State *L)
{
  switch (lua_gettop(L)) {
  case 0: {
    int n = omp_get_max_active_levels();
    lua_pushinteger(L, n);
    return 1;
  };
  case 1: {
    int v = lua_tointeger(L, 1);
    omp_set_max_active_levels(v);
  } break;
  default:
    luaL_error(L, "illegal arguments");
    break;
  }
  return 0;
}

static int
qomp_thread_limit(lua_State *L)
{
  switch (lua_gettop(L)) {
  case 0: {
    int n = omp_get_thread_limit();
    lua_pushinteger(L, n);
    return 1;
  };
  default:
    luaL_error(L, "illegal arguments");
    break;
  }
  return 0;
}

static int
qomp_default_dev(lua_State *L)
{
  switch (lua_gettop(L)) {
  case 0: {
    int n = omp_get_default_device();
    lua_pushinteger(L, n);
    return 1;
  };
  case 1: {
    int v = lua_tointeger(L, 1);
    omp_set_default_device(v);
  } break;
  default:
    luaL_error(L, "illegal arguments");
    break;
  }
  return 0;
}

static int
qomp_num_devices(lua_State *L)
{
  switch (lua_gettop(L)) {
  case 0: {
    int n = omp_get_num_devices();
    lua_pushinteger(L, n);
    return 1;
  };
  default:
    luaL_error(L, "illegal arguments");
    break;
  }
  return 0;
}

static struct {
  omp_sched_t   sched;
  const char *  name;
} sched[] = {
  { omp_sched_static,    "static" },
  { omp_sched_dynamic,   "dynamic" },
  { omp_sched_guided,    "guided" },
  { omp_sched_auto,      "auto" },
  { -1,                  NULL }
};

static int
qomp_schedule(lua_State *L)
{
  switch (lua_gettop(L)) {
  case 0: {
    int mod = 0;
    omp_sched_t kind = omp_sched_auto;
    int i;
    const char *name = NULL;

    omp_get_schedule(&kind, &mod);
    for (i = 0; sched[i].name; i++) {
      if (sched[i].sched == kind) {
	name = sched[i].name;
	break;
      }
    }
    qlua_assert(name != NULL, "Unknown schedule kind");
    lua_pushstring(L, name);
    lua_pushinteger(L, mod);
    return 2;
  };
  case 2: {
    const char *name = lua_tostring(L, 1);
    omp_sched_t kind = omp_sched_auto;
    int mod = lua_tointeger(L, 2);
    int i;

    for (i = 0; sched[i].name; i++) {
      if (strcmp(sched[i].name, name) == 0) {
	kind = sched[i].sched;
	break;
      }
    }
    qlua_assert(sched[i].name != NULL, "Unknown schedule kind");
    omp_set_schedule(kind, mod);
  } break;
  default:
    luaL_error(L, "illegal arguments");
    break;
  }
  return 0;
}

static struct {
  omp_proc_bind_t   bind;
  const char       *name;
} binds[] = {
  /* { omp_proc_bind_false, "false" } */
  /* { omp_proc_bind_true,   "true" } */
  { omp_proc_bind_master,   "master" },
  { omp_proc_bind_close,    "close" },
  { omp_proc_bind_spread,   "spread" },
  { omp_proc_bind_false,     NULL }
};

static int
qomp_bind(lua_State *L)
{
  switch (lua_gettop(L)) {
  case 0: {
    omp_proc_bind_t kind = omp_get_proc_bind();

    switch (kind) {
    case omp_proc_bind_false:
      lua_pushboolean(L, 0);
      break;
    case omp_proc_bind_true:
      lua_pushboolean(L, 1);
      break;
    default: {
      int i;
      const char *name = NULL;
      for (i = 0; binds[i].name; i++) {
	if (binds[i].bind == kind) {
	  name = binds[i].name;
	  break;
	}
      }
      qlua_assert(name != NULL, "Unknown bind value");
      lua_pushstring(L, name);
      break;
    }
    }
    return 1;
  }
  default:
    luaL_error(L, "illegal arguments");
    break;
  }
  return 0;
}

static const luaL_Reg mtOpenMP[] = {
  { "dynamic",              qomp_dynamic},       /*  bool         r/w */
  { "nested",               qomp_nested},        /*  bool         r/w */
  { "num_threads",          qomp_num_threads},   /*  int          r/w */
  { "thread_limit",         qomp_thread_limit},  /*  int          r/- */
  { "levels",               qomp_levels},        /*  int          r/w */
  { "default_dev",          qomp_default_dev},   /*  int          r/w */
  { "num_dev",              qomp_num_devices},   /*  int          r/- */
  { "schedule",             qomp_schedule},      /*  string,int   r/w */
  { "binding",              qomp_bind},          /*  string       r/- */
  { NULL,                   NULL}
};

int
init_qomp(lua_State *L)
{
  char version[64];

  sprintf(version, "%d", _OPENMP);
  lua_getglobal(L, "_G"); 
  lua_newtable(L);
  lua_pushstring(L, version);
  lua_setfield(L, -2, "version");
  luaL_register(L, NULL, mtOpenMP);
  lua_setfield(L, -2, qomp_name);
  lua_pop(L, 1);
  return 0;
}

void fini_qomp(void)
{
}
#endif /* defined(_OPENMP) */
