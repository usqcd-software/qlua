#include "modules.h"                                                 /* DEPS */
#include "qlua.h"                                                    /* DEPS */
#include "lattice.h"                                                 /* DEPS */
#include "qqopqdp.h"                                                 /* DEPS */
#include "latcolmat.h"                                               /* DEPS */
#include "latdirferm.h"                                              /* DEPS */

#include <string.h>
#if USE_Nc3
#define QOP_PRECISION 'D'
#define QOP_Colors 3
#include <qdp.h>
#include <qdp.h>
#include <qop_qdp.h>

static const char QOPwmgStateName[]         = "qop.WilsonMG";

#define Nc 3       /* Only Nc = 3 is supported */
#define MG_DIM 4   /* Only dim = 4 is supported */
#define Qs(a)   a ## 3

/* default values -- borrowed from James's example */
#define QQQ_RSQMIN_DEFAULT        1e-8 /* 0 to ignore */
#define QQQ_RELMIN_DEFAULT        0 /* 0 to ignore */
#define QQQ_MAX_ITER_DEFAULT      600
#define QQQ_RESTART_DEFAULT       200
#define QQQ_MAX_RESTARTS_DEFAULT  5
#define QQQ_EVENODD_DEFAULT       QOP_EVENODD

typedef struct {
  char                       *name;
  int                         nc;
  int                         state;
  double                      kappa;
  double                      c_sw;
  int                         has_bc;
  QOP_Complex                 phase[MG_DIM];
  QOP_wilson_coeffs_t         coeffs;
  QOP_layout_t                layout;
  QOP_info_t                  info;
  QOP_invert_arg_t            inv_arg;
  QOP_resid_arg_t             res_in;
  QOP_resid_arg_t             res_out;
  QOP_D3_FermionLinksWilson  *flw;
  QOP_3_WilsonMg             *wilmg;
} mQOPwmgState;

static int
push_info(lua_State *L, QOP_info_t *info)
{
  lua_newtable(L);
  qlua_push_key_number(L, -1, "final_sec", info->final_sec);
  qlua_push_key_number(L, -1, "final_flops", info->final_flop);
  qlua_push_key_number(L, -1, "count1", info->count1);
  qlua_push_key_number(L, -1, "count2", info->count2);
  return 1;
}

static int
push_residual(lua_State *L, QOP_resid_arg_t *res)
{
  lua_newtable(L);
  qlua_push_key_number(L, -1, "rsqmin", res->rsqmin);
  qlua_push_key_number(L, -1, "final_rsq", res->final_rsq);
  qlua_push_key_number(L, -1, "relmin", res->relmin);
  qlua_push_key_number(L, -1, "final_rel", res->final_rel);
  qlua_push_key_number(L, -1, "final_iter", res->final_iter);
  qlua_push_key_number(L, -1, "final_restart", res->final_restart);
  return 1;
}

static int
load_residual_args(lua_State *L, int tidx, QOP_resid_arg_t *res)
{
  QOP_resid_arg_t resz = QOP_RESID_ARG_DEFAULT;
  *res = resz;

  if (qlua_tabkey_tableopt(L, tidx, "residual")) {
    int ridx = lua_gettop(L);
    res->rsqmin = qlua_tabkey_doubleopt(L, ridx, "rsqmin", QQQ_RSQMIN_DEFAULT);
    res->relmin = qlua_tabkey_doubleopt(L, ridx, "relmin", QQQ_RELMIN_DEFAULT);
    lua_pop(L, 1);
    return 1;
  }
  return 0;
}

static int
load_inverter_args(lua_State *L, int tidx, QOP_invert_arg_t *inv_arg)
{
  QOP_invert_arg_t invz = QOP_INVERT_ARG_DEFAULT;
  *inv_arg = invz;

  if (qlua_tabkey_tableopt(L, tidx, "inverter")) {
    int iidx = lua_gettop(L);
    inv_arg->max_iter = qlua_tabkey_intopt(L, iidx, "max_iter", QQQ_MAX_ITER_DEFAULT);
    inv_arg->restart = qlua_tabkey_intopt(L, iidx, "restart", QQQ_RESTART_DEFAULT);
    inv_arg->max_restarts = qlua_tabkey_intopt(L, iidx, "max_restarts", QQQ_MAX_RESTARTS_DEFAULT);
    const char *eo = qlua_tabkey_stringopt(L, iidx, "evenodd", NULL);
    if (eo) {
      if (!strcmp(eo, "evenodd")) inv_arg->evenodd = QOP_EVENODD;
      if (!strcmp(eo, "even")) inv_arg->evenodd = QOP_EVEN;
      if (!strcmp(eo, "odd")) inv_arg->evenodd = QOP_ODD;
    } else {
      inv_arg->evenodd = QOP_EVENODD;
    }
    lua_pop(L, 1);
    return 1;
  }
  return 0;
}

static int
load_action(lua_State *L, int tidx, mQOPwmgState *wmg)
{
  if (!qlua_tabkey_tableopt(L, tidx, "action"))
    luaL_error(L, "no action record found");

  wmg->kappa = qlua_tabkey_double(L, -1, "kappa");
  wmg->c_sw = qlua_tabkey_doubleopt(L, -1, "c_sw", 0.0);
  wmg->coeffs.clov_s = wmg->c_sw;
  wmg->coeffs.clov_t = wmg->c_sw;
  wmg->coeffs.aniso = 1;
  if (qlua_tabkey_tableopt(L, -1, "boundary")) {
    int i;
    QLA_D_Complex qla_phase[MG_DIM];
    qlua_checkcomplexarray(L, -1, MG_DIM, qla_phase);
    lua_pop(L, 1);
    wmg->has_bc = 1;
    for (i = 0; i < MG_DIM; i++) {
      wmg->phase[i].re = QLA_real(qla_phase[i]);
      wmg->phase[i].im = QLA_imag(qla_phase[i]);
    }
  }
  lua_pop(L, 1);
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
  QOP_3_wilsonMgSet(wmg->wilmg, -1, "nlevels", nlevels);
  QOP_3_wilsonMgSet(wmg->wilmg, -2, "verbose", verbose);
  QOP_3_wilsonMgSet(wmg->wilmg, -1, "nc", wmg->nc);

  /* insert all elements of arg[2].global */
  if (qlua_tabkey_tableopt(L, tidx, "global")) {
    int gidx = lua_gettop(L);
    lua_pushnil(L);
    while (lua_next(L, gidx)) {
      const char *key = lua_tostring(L, -2);
      double val = lua_tonumber(L, -1);
      QOP_3_wilsonMgSet(wmg->wilmg, -1, (char *)key, val); /* drop const qualifier in the key */
      lua_pop(L, 1);
    }
    lua_pop(L, 1);
  }
  QOP_3_wilsonMgSet(wmg->wilmg, -1, "kappa", wmg->kappa);
  QOP_3_wilsonMgSet(wmg->wilmg, -1, "kappanv", wmg->kappa);
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
      QOP_3_wilsonMgSetArray(wmg->wilmg, level, "lattice", sizes, dim);
      qlua_free(L, sizes);
    } else {
      double val = lua_tonumber(L, lua_gettop(L));
      QOP_3_wilsonMgSet(wmg->wilmg, level, (char *)key, val); /* drop const qualifier in the key */
    }
    lua_pop(L, 1);
  }
  lua_pop(L, 1);
  return 0;
}

static void
wmg_do_close(lua_State *L, mQOPwmgState *wmg)
{
  QLUA_ASSERT(wmg->state != 0);
  if (wmg->wilmg)
    QOP_3_wilsonMgFree(wmg->wilmg);
  if (wmg->flw)
    QOP_D3_wilson_destroy_L(wmg->flw);
  if (wmg->layout.latsize)
    qlua_free(L, wmg->layout.latsize);

  wmg->layout.latsize = NULL;
  wmg->flw = NULL;
  wmg->wilmg = NULL;
  wmg->state = 0;
}

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
  mQOPwmgState *wmg = qlua_checkQOPwmgState(L, 1, NULL, -1, 0);

  if (wmg->state) {
    if (wmg->name) {
      lua_pushfstring(L, "%s(%s)", QOPwmgStateName, wmg->name);
    } else {
      lua_pushfstring(L, "%s(%p)", QOPwmgStateName, wmg);
    }
  } else {
    lua_pushfstring(L, "%s(closed)", QOPwmgStateName);
  }
  return 1;
}

static int
wmg_gc(lua_State *L)
{
  mQOPwmgState *wmg = qlua_checkQOPwmgState(L, 1, NULL, -1, 0);
  if (wmg->state)
    wmg_do_close(L, wmg);
  return 0;
}

static int
wmg_info(lua_State *L)
{
  mQOPwmgState *wmg = qlua_checkQOPwmgState(L, 1, NULL, -1, 1);
  push_info(L, &wmg->info);
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
wmg_status(lua_State *L)
{
  mQOPwmgState *wmg = qlua_checkQOPwmgState(L, 1, NULL, -1, 0);

  lua_newtable(L);
  push_info(L, &wmg->info);
  qlua_push_key_object(L, -2, "info");
  push_residual(L, &wmg->res_out);
  qlua_push_key_object(L, -2, "residual");
  return 1;
}

static int
wmg_operator(lua_State *L)
{
  mQOPwmgState *wmg = qlua_checkQOPwmgState(L, 1, NULL, -1, 1);
  mLattice *S = qlua_ObjLattice(L, 1);
  int Sidx = lua_gettop(L);
  QLA_D_Real scale;

  Qs(mLatDirFerm) *rhs = Qs(qlua_checkLatDirFerm)(L, 2, S, wmg->nc);
  Qs(mLatDirFerm) *res = Qs(qlua_newZeroLatDirFerm)(L, Sidx, wmg->nc);

  wmg->info.final_flop = 0;
  wmg->info.final_sec = 0;
  scale = 2 * wmg->kappa;
  QMP_barrier();
  QOP_D3_wilson_dslash_qdp(&wmg->info, wmg->flw, wmg->kappa, +1,
			   res->ptr, rhs->ptr,
			   QOP_EVENODD, QOP_EVENODD);
  QMP_barrier();
  QDP_D3_D_eq_r_times_D(res->ptr, &scale, res->ptr, S->all);

  return 1;
}

static int
wmg_solve(lua_State *L)
{
  int narg = lua_gettop(L);
  mQOPwmgState *wmg = qlua_checkQOPwmgState(L, 1, NULL, -1, 1);
  mLattice *S = qlua_ObjLattice(L, 1);
  int Sidx = lua_gettop(L);
  Qs(mLatDirFerm) *rhs = Qs(qlua_checkLatDirFerm)(L, 2, S, wmg->nc);
  Qs(mLatDirFerm) *tmp = Qs(qlua_newLatDirFerm)(L, Sidx, wmg->nc);
  Qs(mLatDirFerm) *lhs = Qs(qlua_newZeroLatDirFerm)(L, Sidx, wmg->nc);
  QOP_invert_arg_t inv_arg = wmg->inv_arg;
  QLA_Real scale = 1/(2 * wmg->kappa);

  wmg->res_out = wmg->res_in;
  if (narg > 2) {
    qlua_checktable(L, 3, "wilsonMG:solve()");
    QOP_resid_arg_t res_c;
    QOP_invert_arg_t inv_c;
    if (load_residual_args(L, 3, &res_c)) wmg->res_out= res_c;
    if (load_inverter_args(L, 3, &inv_c)) inv_arg = inv_c;
  }
  wmg->info.final_flop = 0;
  wmg->info.final_sec = 0;
  QDP_D3_D_eq_r_times_D(tmp->ptr, &scale, rhs->ptr, S->all);
  QMP_barrier();
  QOP_D3_wilsonMgSolve(&wmg->info, wmg->wilmg, wmg->flw, &inv_arg, &wmg->res_out, wmg->kappa,
		       lhs->ptr, tmp->ptr);
  QMP_barrier();
  lua_newtable(L);
  push_info(L, &wmg->info);
  qlua_push_key_object(L, -2, "info");
  push_residual(L, &wmg->res_out);
  qlua_push_key_object(L, -2, "residual");
  return 2;
}

static int
wmg_close(lua_State *L)
{
  mQOPwmgState *wmg = qlua_checkQOPwmgState(L, 1, NULL, -1, 1);

  wmg_do_close(L, wmg);
  return 0;
}

static int
wmg_colors(lua_State *L)
{
  mQOPwmgState *wmg = qlua_checkQOPwmgState(L, 1, NULL, -1, 0);
  lua_pushnumber(L, wmg->nc);
  return 1;
}

static const struct luaL_Reg mtQOPwmgState[] = {
  { "__tostring",       wmg_fmt       },
  { "__gc",             wmg_gc        },
  { "info",             wmg_info      },
  { "inverter",         wmg_inverter  },
  { "status",           wmg_status    },
  { "operator",         wmg_operator  },
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
    mQOPwmgState *wmg = lua_newuserdata(L, sizeof (mQOPwmgState));
    QOP_layout_t layz = QOP_LAYOUT_ZERO;
    QOP_info_t infoz = QOP_INFO_ZERO;
    QOP_invert_arg_t invz = QOP_INVERT_ARG_DEFAULT;
    QOP_resid_arg_t resz = QOP_RESID_ARG_DEFAULT;
    QOP_wilson_coeffs_t coeffz = QOP_WILSON_COEFFS_ZERO;

    wmg->name = NULL;
    wmg->nc = 0;
    wmg->state = 0;
    wmg->kappa = 0.0;
    wmg->c_sw = 0.0;
    wmg->has_bc = 0;
    wmg->layout = layz;
    wmg->coeffs = coeffz;
    wmg->info = infoz;
    wmg->inv_arg = invz;
    wmg->inv_arg.max_iter = QQQ_MAX_ITER_DEFAULT;
    wmg->inv_arg.restart = QQQ_RESTART_DEFAULT;
    wmg->inv_arg.max_restarts = QQQ_MAX_RESTARTS_DEFAULT;
    wmg->inv_arg.evenodd = QQQ_EVENODD_DEFAULT;
    wmg->res_in = resz;
    wmg->res_out = resz;
    wmg->flw = NULL;
    wmg->wilmg = NULL;

    qlua_createLatticeTable(L, Sidx, mtQOPwmgState, qQOPwmgState, QOPwmgStateName);
    lua_setmetatable(L, -2);

    return wmg;
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
  QDP_D3_ColorMatrix *Ugauge[MG_DIM]; /* FIXED Nc, FIXED dim */
  for (i = 0; i < MG_DIM; i++) {
    lua_pushnumber(L, i + 1);
    lua_gettable(L, 1);
    Ugauge[i] = qlua_checkLatColMat3(L, -1, S, Nc)->ptr;
    lua_pop(L, 1);
  }
  CALL_QDP(L);
  /* construct the QOPwmgStat object */
  mQOPwmgState *wmg = qlua_newQOPwmgState(L, Sidx);
  wmg->nc = Nc;
  wmg->layout.latdim = dim;
  wmg->layout.latsize = NULL;
  wmg->layout.machdim = -1;
  load_action(L, 2, wmg);
  if (QOP_init(&wmg->layout) != QOP_SUCCESS)
    luaL_error(L, "QOP init failed");
  QOP_verbose(0);

  QOP_D3_GaugeField *gf = QOP_D3_create_G_from_qdp(Ugauge);
  if (wmg->has_bc) {
    int i;
    int r0[MG_DIM];
    QOP_bc_t bc;
    QOP_staggered_sign_t sign;
    for (i = 0; i < MG_DIM; i++)
      r0[i] = 0;
    bc.phase = wmg->phase;
    sign.signmask = NULL;
    QOP_D3_rephase_G(gf, r0, &bc, &sign);
  }
  wmg->flw = QOP_D3_wilson_create_L_from_G(&wmg->info, &wmg->coeffs, gf);
  if (wmg->flw == NULL)
    luaL_error(L, "Wilson MG link create failed");
  QOP_D3_destroy_G(gf);

  wmg->wilmg = QOP_3_wilsonMgNew(S->lat);
  if (wmg->wilmg == NULL)
    luaL_error(L, "Wilson MG new failed");
  load_residual_args(L, 2, &wmg->res_in);
  load_inverter_args(L, 2, &wmg->inv_arg);
  load_globals(L, 2, wmg);
  if (qlua_tabkey_tableopt(L, 2, "multigrid")) {
    int nlevels = lua_objlen(L, -1);
    int i;
    for (i = 0; i < nlevels; i++)
      load_mglevel(L, -1, i, wmg);
    lua_pop(L, 1);
  }
  QOP_D3_wilsonMgSetLinks(wmg->wilmg, wmg->flw);
  QOP_3_wilsonMgSetup(wmg->wilmg);
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
