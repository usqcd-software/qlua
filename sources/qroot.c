/* GSL root finders */
#include "qlua.h"                                                           /* DEPS */
#include "qroot.h"                                                          /* DEPS */
#include "qmatrix.h"                                                        /* DEPS */

enum {
  QROOT_dnewton,
  QROOT_broyden,
  QROOT_hybrid,
  QROOT_hybrids,
};

/* root finder parameters names */
static const char *root_params[] = {
  "MaxIterations",  /* maximal number of iterations to use */
  "Logging",        /* if true, record minimization progress */
  "RelativeError",  /* stop if the relative error is below this number */
  "AbsoluteError",  /* stop if the absolute error is below this number */
  NULL
};

/* root finder paramters indices in the closure
 * These must agree with root_params order 
 */
enum {
  QS_kind = 1,
  QS_max_iter,
  QS_logging,
  QS_rel_err,
  QS_abs_err,
  QS_param_count
};

struct frootN_params {
  int ndim;
  lua_State *L;
  void *func;
};

static int
frootN_func(const gsl_vector *x, void *p, gsl_vector *f)
{
  struct frootN_params *params = p;
  lua_State *L = params->L;
  int ndim = params->ndim;
  int i, j;

  lua_pushlightuserdata(L, params);
  lua_gettable(L, LUA_REGISTRYINDEX);
  for (i = 0; i < ndim; i++) {
    lua_pushinteger(L, i+1);
    lua_gettable(L, -2);
    for (j = 0; j < ndim; j++) {
      double xj = gsl_vector_get(x, j);
      lua_pushnumber(L, xj);
    }
    lua_call(L, ndim, 1);
    gsl_vector_set(f, i, luaL_checknumber(L, -1));
    lua_pop(L, 1);
  }
  lua_pop(L, 1);

  return GSL_SUCCESS;
}

static int
frootN(lua_State *L, int idx_x)
{
  const gsl_multiroot_fsolver_type *T = NULL;
  gsl_multiroot_fsolver *s = NULL;
  char *name = NULL;
  struct frootN_params params;
  gsl_multiroot_function func;
  int ndim = lua_objlen(L, 1);
  int max_iter;
  int logging;
  double rel_err;
  double abs_err;
  gsl_vector *x = NULL;
  int iter;
  int i;
  int status = 1;

  if (ndim < 1)
    luaL_error(L, "Dimension too small in solver");

  if (lua_objlen(L, idx_x) != ndim)
    luaL_error(L, "wrong number of variables in solver");

  switch (luaL_checkint(L, lua_upvalueindex(QS_kind))) {
  case QROOT_dnewton:
    T = gsl_multiroot_fsolver_dnewton;
    name = "dnewton";
    break;
  case QROOT_broyden:
    T = gsl_multiroot_fsolver_broyden;
    name = "broyden";
    break;
  case QROOT_hybrid:
    T = gsl_multiroot_fsolver_hybrid;
    name = "hybrid";
    break;
  case QROOT_hybrids:
    T = gsl_multiroot_fsolver_hybrids;
    name = "hybrids";
    break;
  default:
    luaL_error(L, "internal error: unexpected solver");
    return 0;
  }
  x = new_gsl_vector(L, ndim);
  for (i = 0; i < ndim; i++) {
    double xi;
    lua_pushnumber(L, i + 1);
    lua_gettable(L, idx_x);
    xi = luaL_checknumber(L, -1);
    lua_pop(L, 1);
    gsl_vector_set(x, i, xi);
  }

  max_iter = luaL_optint(L, lua_upvalueindex(QS_max_iter), 100);
  rel_err = luaL_optnumber(L, lua_upvalueindex(QS_rel_err), 0.0);
  abs_err = luaL_optnumber(L, lua_upvalueindex(QS_abs_err), 0.0);
  logging = lua_toboolean(L, lua_upvalueindex(QS_logging));

  params.func = &params;
  params.ndim = ndim;
  params.L = L;
  lua_pushlightuserdata(L, &params);
  lua_pushvalue(L, 1);
  lua_settable(L, LUA_REGISTRYINDEX);
  func.params = &params;
  func.n = ndim;
  func.f = frootN_func;
  
  s = gsl_multiroot_fsolver_alloc (T, ndim);
  if (s == NULL) {
    lua_gc(L, LUA_GCCOLLECT, 0);
    s = gsl_multiroot_fsolver_alloc (T, ndim);
    if (s == NULL)
      luaL_error(L, "not enough memory");
  }
  gsl_multiroot_fsolver_set(s, &func, x);

  lua_pushnil(L);
  lua_createtable(L, 0, logging?4:3);
  lua_pushstring(L, name);
  lua_setfield(L, -2, "Name");
  if (logging) {
    lua_createtable(L, 0, 2); /* Logs */
    lua_createtable(L, max_iter, 0); /* X */
    lua_createtable(L, max_iter, 0); /* f */
  }

  for (iter = 1; iter < max_iter; iter++) {
    if (gsl_multiroot_fsolver_iterate(s) != 0) {
      iter --;
      status = 2;
      break;
    }
    if (logging) {
      lua_createtable(L, ndim, 0); /* x */
      lua_createtable(L, ndim, 0); /* f */
      for (i = 0; i < ndim; i++) {
        lua_pushnumber(L, gsl_vector_get(s->x, i));
        lua_rawseti(L, -3, i+1);
        lua_pushnumber(L, gsl_vector_get(s->f, i));
        lua_rawseti(L, -2, i+1);
      }
      lua_rawseti(L, -3, iter);
      lua_rawseti(L, -3, iter);
    }
    if (gsl_multiroot_test_delta(s->dx, s->x, abs_err, rel_err) == GSL_SUCCESS) {
      status = 0;
      break;
    }
  }

  if (logging) {
    lua_setfield(L, -3, "f");
    lua_setfield(L, -2, "x");
    lua_setfield(L, -2, "Logs");
  }

  lua_pushstring(L, status == 0? "OK": "FAILED");
  lua_setfield(L, -2, "Status");
  lua_pushinteger(L, iter);
  lua_setfield(L, -2, "Iterations");

  lua_createtable(L, ndim, 0);
  for (i = 0; i < ndim; i++) {
    lua_pushnumber(L, gsl_vector_get(s->x, i));
    lua_rawseti(L, -2, i+1);
  }
  lua_replace(L, -3);

  lua_pushlightuserdata(L, &params);
  lua_pushnil(L);
  lua_settable(L, LUA_REGISTRYINDEX);
  gsl_multiroot_fsolver_free(s);
  gsl_vector_free(x);

  return 2;
}


static int
frooter(lua_State *L)
{
  int narg = lua_gettop(L);
  
  switch (narg) {
  case 2:
    return frootN(L, 2);
  default:
    luaL_error(L, "bad number of arguments");
    break;
  }
  return 0;
}

static int
build_rooter(lua_State *L, int id)
{
  int i;
  lua_pushinteger(L, id);
  for (i = 0; root_params[i]; i++)
    lua_getfield(L, 1, root_params[i]);
  lua_pushcclosure(L, frooter, i + 1);
  return 1;
}

#define ROOT_IF(r) static int froot_##r(lua_State *L) \
  { return build_rooter(L, QROOT_##r); }

ROOT_IF(dnewton)
ROOT_IF(broyden)
ROOT_IF(hybrid)
ROOT_IF(hybrids)

static const luaL_Reg fRoot[] = {
  { "dnewton",  froot_dnewton  },
  { "broyden",  froot_broyden  },
  { "hybrid",   froot_hybrid   },
  { "hybrids",  froot_hybrids  },
  { NULL,       NULL           }
};

int
init_root(lua_State *L)
{
  gsl_set_error_handler_off();
  lua_getglobal(L, "gsl");
  lua_newtable(L);
  luaL_register(L, NULL, fRoot);
  lua_setfield(L, -2, "root");
  lua_pop(L, 1);
  return 0;
}

int
fini_root(lua_State *L)
{
  return 0;
}
