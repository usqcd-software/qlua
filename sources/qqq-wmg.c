#include "modules.h"                                                 /* DEPS */
#include "qlua.h"                                                    /* DEPS */
#include "lattice.h"                                                 /* DEPS */
#include "latcolmat.h"                                               /* DEPS */
#include "qqopqdp.h"                                                 /* DEPS */
#include <string.h>
#ifdef USE_Nc3
#define QOP_PRECISION 'D'
#define QOP_Colors 3
#include <qdp.h>
#include <qdp.h>
#include <qop_qdp.h>

static const char QOPwmgStateName[]         = "qop.WilsonMG";

#define Nc 3       /* Only Nc = 3, see FIXED Nc comments */
#define MG_DIM 4   /* Only dim = 4, see FIXED dim comments */

/* default values -- borrowed from James's example */
#define QQQ_RSQMIN_DEFAULT        1e-8
#define QQQ_MAX_ITER_DEFAULT      600
#define QQQ_RESTART_DEFAULT       200
#define QQQ_MAX_RESTARTS_DEFAULT  5
#define QQQ_EVENODD_DEFAULT       QOP_EVENODD

typedef struct {
  char                       *name;
  int                         nc;
  int                         state;
  QOP_layout_t                layout;
  QOP_info_t                  info;
  QOP_invert_arg_t            inv_arg;
  QOP_resid_arg_t             res_arg;
  QOP_FermionLinksWilson     *flw;
  QOP_WilsonMg               *wilmg;
} mQOPwmgState;

static mQOPwmgState *
qlua_checkQOPwmgState(lua_State *L, int idx, mLattice *S, int nc, int live)
{
  mQOPwmgState *v = qlua_checkLatticeType(L, idx, qQOPwmgState, QOPwmgStateName);
  
  if (S) {
    mLattice *S1 = qlua_ObjLattice(L, idx);
    if (S1->id != S->id)
      luaL_error(L, "%s on a wrong lattice", QOPwmgStateName);
    lua_pop(L, 1);
  }
  if (nc != -1) {
    if (v->nc != nc)
      luaL_error(L, "Wrong number of colors");
  }
  if (live && (v->state == 0)) {
    luaL_error(L, "Expecting live %s", QOPwmgStateName);
  }
  return v;
}

static int
wmg_fmt(lua_State *L)
{
  mQOPwmgState *v = qlua_checkQOPwmgState(L, 1, NULL, -1, 0);

  if (v->state) {
    if (v->name) {
      lua_pushfstring(L, "%s(%s)", QOPwmgStateName, v->name);
    } else {
      lua_pushfstring(L, "%s(%p)", QOPwmgStateName, v);
    }
  } else {
    lua_pushfstring(L, "%s(closed)", QOPwmgStateName);
  }
  return 1;
}

static void
wmg_do_close(lua_State *L, mQOPwmgState *v)
{
  QLUA_ASSERT(v->state != 0);
  if (v->wilmg)
    QOP_wilsonMgFree(v->wilmg);
  if (v->flw)
    QOP_wilson_destroy_L(v->flw);
  if (v->layout.latsize)
    qlua_free(L, v->layout.latsize);

  v->layout.latsize = NULL;
  v->flw = NULL;
  v->wilmg = NULL;
  v->state = 0;
}

static int
wmg_gc(lua_State *L)
{
  mQOPwmgState *v = qlua_checkQOPwmgState(L, 1, NULL, -1, 0);
  if (v->state)
    wmg_do_close(L, v);
  return 0;
}

static int
wmg_info(lua_State *L)
{
  mQOPwmgState *wmg = qlua_checkQOPwmgState(L, 1, NULL, -1, 1);
  lua_newtable(L);
  qlua_push_key_number(L, -1, "final_sec", wmg->info.final_sec);
  qlua_push_key_number(L, -1, "final_flops", wmg->info.final_flop);
  qlua_push_key_number(L, -1, "count1", wmg->info.count1);
  qlua_push_key_number(L, -1, "count2", wmg->info.count2);
  return 1;
}

static int
wmg_inverter(lua_State *L)
{
  mQOPwmgState *wmg = qlua_checkQOPwmgState(L, 1, NULL, -1, 1);
  lua_newtable(L);
  qlua_push_key_number(L, -1, "mixed_rsq", wmg->inv_arg.mixed_rsq);
  qlua_push_key_number(L, -1, "max_iter", wmg->inv_arg.max_iter);
  qlua_push_key_number(L, -1, "max_restarts", wmg->inv_arg.max_restarts);
  qlua_push_key_number(L, -1, "restart", wmg->inv_arg.restart);
  const char *eo = NULL;
  switch (wmg->inv_arg.evenodd) {
  case QOP_EVENODD: eo = "evenodd"; break;
  case QOP_EVEN: eo = "even"; break;
  case QOP_ODD: eo = "odd"; break;
  }
  if (eo)
    qlua_push_key_string(L, -1, "evenodd", eo);
  return 1;
}

static int
wmg_result(lua_State *L)
{
  lua_pushstring(L, "XXX wmg method called result");
  return 1;
}

static int
wmg_status(lua_State *L)
{
  lua_pushstring(L, "XXX wmg method called status");
  return 1;
}

static int
wmg_solve(lua_State *L)
{
  lua_pushstring(L, "XXX wmg method called solve");
  return 1;
}

static int
wmg_close(lua_State *L)
{
  mQOPwmgState *v = qlua_checkQOPwmgState(L, 1, NULL, -1, 1);

  wmg_do_close(L, v);
  return 0;
}

static int
wmg_colors(lua_State *L)
{
  mQOPwmgState *v = qlua_checkQOPwmgState(L, 1, NULL, -1, 0);
  lua_pushnumber(L, v->nc);
  return 1;
}

static const struct luaL_Reg mtQOPwmgState[] = {
  { "__tostring",       wmg_fmt       },
  { "__gc",             wmg_gc        },
  { "info",             wmg_info      },
  { "result",           wmg_result    },
  { "inverter",         wmg_inverter  },
  { "status",           wmg_status    },
  { "solve",            wmg_solve     },
  { "close",            wmg_close     },
  { "colors",           wmg_colors    },
  /* "lattice" */
  /* "a-type" */
  { NULL,               NULL          }
};

static mQOPwmgState *
qlua_newQOPwmgState(lua_State *L, int Sidx)
{
    mQOPwmgState *c = lua_newuserdata(L, sizeof (mQOPwmgState));
    QOP_layout_t layz = QOP_LAYOUT_ZERO;
    QOP_info_t infoz = QOP_INFO_ZERO;
    QOP_invert_arg_t invz = QOP_INVERT_ARG_DEFAULT;
    QOP_resid_arg_t resz = QOP_RESID_ARG_DEFAULT;

    c->name = NULL;
    c->nc = 0;
    c->state = 0;
    c->layout = layz;
    c->info = infoz;
    c->inv_arg = invz;
    c->inv_arg.max_iter = QQQ_MAX_ITER_DEFAULT;
    c->inv_arg.restart = QQQ_RESTART_DEFAULT;
    c->inv_arg.max_restarts = QQQ_MAX_RESTARTS_DEFAULT;
    c->inv_arg.evenodd = QQQ_EVENODD_DEFAULT;
    c->res_arg = resz;
    c->flw = NULL;
    c->wilmg = NULL;

    qlua_createLatticeTable(L, Sidx, mtQOPwmgState, qQOPwmgState, QOPwmgStateName);
    lua_setmetatable(L, -2);

    return c;
}

static int
load_inverter_args(lua_State *L, int tidx, QOP_invert_arg_t *inv_arg, QOP_resid_arg_t *res_arg)
{
  QOP_invert_arg_t invz = QOP_INVERT_ARG_DEFAULT;
  QOP_resid_arg_t resz = QOP_RESID_ARG_DEFAULT;
  
  *inv_arg = invz;
  *res_arg = resz;

  res_arg->rsqmin = qlua_tabkey_doubleopt(L, tidx, "rsqmin", QQQ_RSQMIN_DEFAULT);

  if (qlua_tabkey_tableopt(L, tidx, "inverter")) {
    inv_arg->max_iter = qlua_tabkey_intopt(L, tidx, "max_iter", QQQ_MAX_ITER_DEFAULT);
    inv_arg->restart = qlua_tabkey_intopt(L, tidx, "restart", QQQ_RESTART_DEFAULT);
    inv_arg->max_restarts = qlua_tabkey_intopt(L, tidx, "max_restarts", QQQ_MAX_RESTARTS_DEFAULT);
    const char *eo = qlua_tabkey_stringopt(L, tidx, "evenodd", NULL);
    if (eo) {
      if (!strcmp(eo, "evenodd")) inv_arg->evenodd = QOP_EVENODD;
      if (!strcmp(eo, "even")) inv_arg->evenodd = QOP_EVEN;
      if (!strcmp(eo, "odd")) inv_arg->evenodd = QOP_ODD;
    } else {
      inv_arg->evenodd = QOP_EVENODD;
    }
    lua_pop(L, 1);
  }
  return 0;
}

static int
load_globals(lua_State *L, int tidx, mQOPwmgState *wmg)
{
  int verbose = qlua_tabkey_intopt(L, tidx, "verbose", 0);
  int nlevels = 0;
  if (qlua_tabkey_tableopt(L, tidx, "multigrid")) {
    nlevels = lua_objlen(L, -1);
    lua_pop(L, 1);
  }
  QOP_wilsonMgSet(wmg->wilmg, -1, "nlevels", nlevels);
  QOP_wilsonMgSet(wmg->wilmg, -2, "verbose", verbose);
  QOP_wilsonMgSet(wmg->wilmg, -1, "nc", wmg->nc);

  /* insert all elements of arg[2].global */
  if (qlua_tabkey_tableopt(L, tidx, "global")) {
    int gidx = lua_gettop(L);
    lua_pushnil(L);
    while (lua_next(L, gidx)) {
      const char *key = lua_tostring(L, -2);
      double val = lua_tonumber(L, -1);
      QOP_wilsonMgSet(wmg->wilmg, -1, (char *)key, val); /* drop const qualifier in the key */
      lua_pop(L, 1);
    }
    lua_pop(L, 1);
  }
  return 0;
}

static int
load_mglevel(lua_State *L, int tidx, int level,  mQOPwmgState *wmg)
{
  if (!qlua_tabpushopt_idx(L, tidx, level + 1))
    luaL_error(L, "missing multigrid level %d", level);
  int lidx = lua_gettop(L);
  lua_pushnil(L);
  while (lua_next(L, lidx)) {
    const char *key = lua_tostring(L, -2);
    if (strcmp(key, "lattice") == 0) {
      int dim = 0;
      double *sizes = qlua_numberarray(L, lua_gettop(L), &dim);
      if (dim != MG_DIM)
	luaL_error(L, "bad multigrid lattice rank %d at level %d", dim, level);
      QOP_wilsonMgSetArray(wmg->wilmg, level, "lattice", sizes, dim);
      qlua_free(L, sizes);
    } else {
      double val = lua_tonumber(L, lua_gettop(L));
      QOP_wilsonMgSet(wmg->wilmg, level, (char *)key, val); /* drop const qualifier in the key */
    }
    lua_pop(L, 1);
  }
  lua_pop(L, 1);
  return 0;
}

int
qqq_wmg(lua_State *L)
{
  int i;
  double t0, t1;

  t0 = QDP_time();
  /* check arguments */
  luaL_checktype(L, 1, LUA_TTABLE);
  luaL_checktype(L, 2, LUA_TTABLE);
  /* extract gauge fields and the lattice */
  int dim = lua_objlen(L, 1);
  if (dim != MG_DIM) /* FIXED dim */
    luaL_error(L, "Gauge table must have %d elements", MG_DIM);
  lua_pushnumber(L, 1);
  lua_gettable(L, 1);
  int U0idx = lua_gettop(L);
  qlua_checkLatColMat3(L, U0idx, NULL, Nc); /* FIXED Nc */
  mLattice *S = qlua_ObjLattice(L, U0idx);
  int Sidx = lua_gettop(L);
  if (dim != S->rank)
    luaL_error(L, "Gauge table size=%d, lattice rank=%d", dim, S->rank);
  QDP_F3_ColorMatrix *Ugauge[MG_DIM]; /* FIXED Nc, FIXED dim */
  for (i = 0; i < MG_DIM; i++) {
    lua_pushnumber(L, i + 1);
    lua_gettable(L, 1);
    QDP_D3_ColorMatrix *ui = qlua_checkLatColMat3(L, -1, S, Nc)->ptr;
    Ugauge[i] = QDP_F3_create_M_L(S->lat);
    QDP_FD3_M_eq_M(Ugauge[i], ui, S->all);
    lua_pop(L, 1);
  }
  CALL_QDP(L);
  /* construct the QOPwmgStat object */
  mQOPwmgState *wmg = qlua_newQOPwmgState(L, Sidx);
  wmg->nc = Nc;
  wmg->layout.latdim = dim;
  wmg->layout.latsize = NULL;
  wmg->layout.machdim = -1;
  if (QOP_init(&wmg->layout) != QOP_SUCCESS)
    luaL_error(L, "QOP init failed");
  QOP_verbose(0);

  QOP_wilson_coeffs_t coeffs = QOP_WILSON_COEFFS_ZERO;
  if (qlua_tabkey_tableopt(L, 2, "clover")) {
    coeffs.clov_s = qlua_tabkey_doubleopt(L, -1, "clov_s", 0.0);
    coeffs.clov_t = qlua_tabkey_doubleopt(L, -1, "clov_t", coeffs.clov_s);
    coeffs.aniso = 1;
    lua_pop(L,1);
  }

  QOP_GaugeField *gf = QOP_create_G_from_qdp(Ugauge);
  wmg->flw = QOP_wilson_create_L_from_G(&wmg->info, &coeffs, gf);
  if (wmg->flw == NULL)
    luaL_error(L, "Wilson MG link create failed");
  QOP_destroy_G(gf);

  wmg->wilmg = QOP_3_wilsonMgNew(S->lat);
  if (wmg->wilmg == NULL)
    luaL_error(L, "Wilson MG new failed");
  load_inverter_args(L, 2, &wmg->inv_arg, &wmg->res_arg);
  load_globals(L, 2, wmg);
  if (qlua_tabkey_tableopt(L, 2, "multigrid")) {
    int nlevels = lua_objlen(L, -1);
    int i;
    for (i = 0; i < nlevels; i++)
      load_mglevel(L, -1, i, wmg);
    lua_pop(L, 1);
  }
  QOP_wilsonMgSetLinks(wmg->wilmg, wmg->flw);
  QOP_wilsonMgSetup(wmg->wilmg);
  const char *name = qlua_tabkey_stringopt(L, 2, "name", NULL);
  if (name)
    wmg->name = qlua_strdup(L, name);

  t1 = QDP_time();
  wmg->info.final_sec = t1 - t0;

  wmg->state = 1;
  qqq_inited = 1;
  return 1;
}
#endif /* USE_Nc3 */
