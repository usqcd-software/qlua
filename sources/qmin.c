/* GSL minimizers */
#include "qlua.h"                                                           /* DEPS */
#include "qmin.h"                                                           /* DEPS */
#include "qmatrix.h"                                                        /* DEPS */

enum {
#if 0
  QMIN_steepest_descent,   /* Nd,df */
  QMIN_conjugate_fr,       /* Nd,df */
  QMIN_conjugate_pr,       /* Nd,df */
  QMIN_vector_bfgs,        /* Nd,df */
  QMIN_vector_bfgs2,       /* Nd,df */
#endif
  QMIN_nmsimplex,          /* Nd,f */
  QMIN_nmsimplex2,         /* Nd,f */
  QMIN_nmsimplex2rand      /* Nd,f */
};

/* minimizer paramters names */
static const char *min_params[] = {
  "MaxIterations",  /* maximal number of iterations to use */
  "Logging",        /* if true, record minimization progress */
  "RelativeError",  /* stop if the relative error is below this number */
  "AbsoluteError",  /* stop if the absolute error is below this number */
  "StepSize",       /* step size (n-d for df minimizers, 1-d for f-minimizers */
  NULL
};

/* minimizer paramters indices in the closure
 * These must agree with min_params order 
 */
enum {
  QM_kind = 1,
  QM_max_iter,
  QM_logging,
  QM_rel_err,
  QM_abs_err,
  QM_step_size,
  QM_param_count
};

double
vnorm2_gsl(gsl_vector *v, int dim)
{
  double n = 0;
  int i;
  for (i = 0; i < dim; i++) {
    double x = gsl_vector_get(v, i);
    n = n + x * x;
  }
  return sqrt(n);
}

struct fminN_params {
  int ndim;
  lua_State *L;
  void *func;
};

static double
fminN_func(const gsl_vector *x, void *params_v)
{
  struct fminN_params *params = params_v;
  lua_State *L = params->L;
  int ndim = params->ndim;
  double r;
  int i;

  lua_pushlightuserdata(L, params);
  lua_gettable(L, LUA_REGISTRYINDEX);
  lua_createtable(L, ndim, 0);
  for (i = 0; i < ndim; i++) {
    double xi = gsl_vector_get(x, i);
    lua_pushnumber(L, xi);
    lua_rawseti(L, -2, i+1);
  }
  lua_call(L, 1, 1);
  r = luaL_checknumber(L, -1);
  lua_pop(L, 1);

  return r;
}

static int
fminN(lua_State *L, int idx_x)
{
  const gsl_multimin_fminimizer_type *T = NULL;
  gsl_multimin_fminimizer *s = NULL;
  gsl_vector *x = NULL;
  gsl_vector *v = NULL;
  char *name = "";
  struct fminN_params params;
  gsl_multimin_function func;
  int ndim = lua_objlen(L, idx_x);
  int max_iter;
  int logging;
  double rel_err;
  double abs_err;
  double step_size;
  int iter;
  int status = 1;
  double size = 0.0;
  int i;

  if (ndim < 1)
    luaL_error(L, "Dimension too small in minimizer");

  switch (luaL_checkint(L, lua_upvalueindex(QM_kind))) {
  case QMIN_nmsimplex:
    T = gsl_multimin_fminimizer_nmsimplex;
    name = "nmsimplex";
    break;
  case QMIN_nmsimplex2:
    T = gsl_multimin_fminimizer_nmsimplex2;
    name = "nmsimplex2";
    break;
  case QMIN_nmsimplex2rand:
    T = gsl_multimin_fminimizer_nmsimplex2rand;
    name = "nmsimplex2rand";
    break;
  default:
    luaL_error(L, "internal error: unexpected minimizer");
    return 0;
  }

  max_iter = luaL_optint(L, lua_upvalueindex(QM_max_iter), 100);
  rel_err = luaL_optnumber(L, lua_upvalueindex(QM_rel_err), 0.0);
  abs_err = luaL_optnumber(L, lua_upvalueindex(QM_abs_err), 0.0);
  step_size = luaL_optnumber(L, lua_upvalueindex(QM_step_size), 0.01);
  logging = lua_toboolean(L, lua_upvalueindex(QM_logging));
  x = new_gsl_vector(L, ndim);
  v = new_gsl_vector(L, ndim);
  for (i = 0; i < ndim; i++) {
    double xi;
    lua_pushnumber(L, i + 1);
    lua_gettable(L, idx_x);
    xi = luaL_checknumber(L, -1);
    lua_pop(L, 1);
    gsl_vector_set(x, i, xi);
  }
  gsl_vector_set_all(v, step_size);

  params.func = &params;
  params.ndim = ndim;
  params.L = L;
  lua_pushlightuserdata(L, &params);
  lua_pushvalue(L, 1);
  lua_settable(L, LUA_REGISTRYINDEX);
  func.params = &params;
  func.n = ndim;
  func.f = fminN_func;

  s = gsl_multimin_fminimizer_alloc (T, ndim);
  if (s == NULL) {
    lua_gc(L, LUA_GCCOLLECT, 0);
    s = gsl_multimin_fminimizer_alloc (T, ndim);
    if (s == NULL)
      luaL_error(L, "not enough memory");
  }
  gsl_multimin_fminimizer_set(s, &func, x, v);

  lua_pushnil(L);
  lua_createtable(L, 0, logging?5:4);
  lua_pushstring(L, name);
  lua_setfield(L, -2, "Name");
  if (logging) {
    lua_createtable(L, 0, 3); /* Logs */
    lua_createtable(L, max_iter, 0); /* X */
    lua_createtable(L, max_iter, 0); /* f */
    lua_createtable(L, max_iter, 0); /* size */
  }
  
  for (iter = 1; iter < max_iter; iter++) {
    if (gsl_multimin_fminimizer_iterate(s) != 0) {
      status = 2;
      break;
    }
    size = gsl_multimin_fminimizer_size(s);
    if (logging) {
      lua_pushnumber(L, fabs(size));
      lua_rawseti(L, -2, iter);
      lua_pushnumber(L, s->fval);
      lua_rawseti(L, -3, iter);
      lua_createtable(L, ndim, 0);
      for (i = 0; i < ndim; i++) {
        lua_pushnumber(L, gsl_vector_get(s->x, i));
        lua_rawseti(L, -2, i+1);
      }
      lua_rawseti(L, -4, iter);
    }
    if (fabs(size) < abs_err + rel_err * vnorm2_gsl(s->x, ndim)) {
      status = 0;
      break;
    }
  }
  if (logging) {
    lua_setfield(L, -4, "size");
    lua_setfield(L, -3, "f");
    lua_setfield(L, -2, "x");
    lua_setfield(L, -2, "Logs");
  }
  lua_pushstring(L, status == 0? "OK": "FAILED");
  lua_setfield(L, -2, "Status");
  lua_pushinteger(L, iter);
  lua_setfield(L, -2, "Iterations");
  lua_pushnumber(L, fabs(size));
  lua_setfield(L, -2, "Size");

  lua_createtable(L, ndim, 0);
  for (i = 0; i < ndim; i++) {
    lua_pushnumber(L, gsl_vector_get(s->x, i));
    lua_rawseti(L, -2, i+1);
  }
  lua_replace(L, -3);

  lua_pushlightuserdata(L, &params);
  lua_pushnil(L);
  lua_settable(L, LUA_REGISTRYINDEX);
  gsl_multimin_fminimizer_free(s);
  gsl_vector_free(x);
  gsl_vector_free(v);

  return 2;
}

static int
fdfminN(lua_State *L)
{
  switch (luaL_checkint(L, lua_upvalueindex(QM_kind))) {
#if 0
  case QMIN_steepest_descent:
  case QMIN_conjugate_fr:
  case QMIN_conjugate_pr:
  case QMIN_vector_bfgs:
  case QMIN_vector_bfgs2:
    /* ... */
    break;
#endif
  case QMIN_nmsimplex:
  case QMIN_nmsimplex2:
  case QMIN_nmsimplex2rand:
    return fminN(L, 3);
  default:
    luaL_error(L, "internal error: unexpected minimizer\n");
    return 0;
  }
#if 0
  luaL_error(L, "f/df min N not implemented XXX");
  /* ... */
  return 0;
#endif
}

static int
fminimizer(lua_State *L)
{
  int narg = lua_gettop(L);

  switch (narg) {
  case 2:
    return fminN(L, 2);
  case 3:
    return fdfminN(L);
  default:
    luaL_error(L, "bad number of arguments");
    break;
  }
  return 0;
}

static int
build_minimizer(lua_State *L, int id)
{
  int i;
  lua_pushinteger(L, id);
  for (i = 0; min_params[i]; i++)
    lua_getfield(L, 1, min_params[i]);
  lua_pushcclosure(L, fminimizer, i + 1);
  return 1;
}

#define MIN_IF(r) static int fmin_##r(lua_State *L) \
   { return build_minimizer(L, QMIN_##r); }

#if 0
MIN_IF(steepest_descent)
MIN_IF(conjugate_pr)
MIN_IF(conjugate_fr)
MIN_IF(vector_bfgs)
MIN_IF(vector_bfgs2)
#endif
MIN_IF(nmsimplex)
MIN_IF(nmsimplex2)
MIN_IF(nmsimplex2rand)

static const luaL_Reg fMin[] = {
#if 0
  {"steepest_descent",   fmin_steepest_descent },
  {"conjugate_pr",       fmin_conjugate_pr     },
  {"conjugate_fr",       fmin_conjugate_fr     },
  {"vector_bfgs",        fmin_vector_bfgs      },
  {"vector_bfgs2",       fmin_vector_bfgs2     },
#endif
  {"nmsimplex",          fmin_nmsimplex        },
  {"nmsimplex2",         fmin_nmsimplex2       },
  {"nmsimplex2rand",     fmin_nmsimplex2rand   },
  { NULL,                NULL                  }
};

int
init_min(lua_State *L)
{
  gsl_set_error_handler_off();
  lua_getglobal(L, "gsl");
  lua_newtable(L);
  luaL_register(L, NULL, fMin);
  lua_setfield(L, -2, "min");
  lua_pop(L, 1);
  return 0;
}

int
fini_min(lua_State *L)
{
  return 0;
}
