#include "modules.h"                                                 /* DEPS */
#include "qlua.h"                                                    /* DEPS */
#include "qtwisted.h"                                                /* DEPS */
#include "qcomplex.h"                                                /* DEPS */
#include "qvector.h"                                                 /* DEPS */
#include "lattice.h"                                                 /* DEPS */
#include "qlayout.h"                                                 /* DEPS */
#include "latcolmat.h"                                               /* DEPS */
#include "latdirferm.h"                                              /* DEPS */
#define QOP_TWISTED_DEFAULT_PRECISION QDP_Precision
#include QLUA_TWISTED_HDR
#include "qmp.h"
#include "qlanczos.h"

#include <math.h>
#include <string.h>

#define QNc(a,b) QNc_1(a,QLUA_TWISTED_NC,b)
#define QNc_1(a,b,c) QNc_2(a,b,c)
#define QNc_2(a,b,c) a ## b ## c
#define QNz(a) QNz_1(a,QLUA_TWISTED_NC)
#define QNz_1(a,b) QNz_2(a,b)
#define QNz_2(a,b) a ## b

#if QLUA_TWISTED_NC != QNc(QOP_, _TWISTED_COLORS)
#error "Twisted Nc mismatch"
#endif /*  QLUA_TWISTED_NC != QNc(QOP_, _TWISTED_COLORS) */

static const char TwistedName[]               = "lattice.Twisted";
static const char TwistedDeflatorName[]       = "lattice.Twisted.Deflator";
static const char TwistedDeflatorStateName[]  = "lattice.Twisted.DeflatorState";

typedef int TwistedInverter(lua_State *L,
                           struct QNc(QOP_,_TWISTED_Fermion) *solution,
                           int *out_iters,
                           double *out_epsilon,
                           const struct QNc(QOP_,_TWISTED_Fermion) *rhs,
                           int log_level);

typedef struct {
    TwistedInverter  *proc;
    const char      *name;
} TwistedSolver;

typedef int TwistedMxMInverter(lua_State *L,
                           struct QNc(QOP_,_TWISTED_HalfFermion) *solution,
                           int *out_iters,
                           double *out_epsilon,
                           const struct QNc(QOP_,_TWISTED_HalfFermion) *rhs,
                           int log_level);

typedef struct {
    TwistedMxMInverter  *proc;
    const char      *name;
} TwistedMxMSolver;

typedef struct {
  struct QNc(QOP_, _TWISTED_State) *state;
  struct QNc(QOP_, _TWISTED_Gauge) *gauge;
} mTwisted;

typedef struct {
  int nev;
  int umax;
  int vmax;
  struct QNc(QOP_F, _TWISTED_Deflator) *deflator;
} mDeflatorState;

static mTwisted *qlua_newTwisted(lua_State *L, int Sidx);
static mTwisted *qlua_checkTwisted(lua_State *L,
                                 int idx,
                                 mLattice *S,
                                 int live);

typedef struct {
  QDP_Lattice *lat;
  QNc(QLA_D, _DiracFermion) *f;
  double s;
} TW_D_env;

static double
q_TW_D_reader_scaled(const int p[], int c, int d, int re_im, void *e)
{
  TW_D_env *env = e;
  int i = QDP_index_L(env->lat, p);
  QNc(QLA_D, _DiracFermion) *f = env->f;
  double s = env->s;
  QLA_D_Real xx;
  
  if (re_im == 0) {
    QLA_r_eq_Re_c(xx, QLA_elem_D(f[i], c, d));
  } else {
    QLA_r_eq_Im_c(xx, QLA_elem_D(f[i], c, d));
  }
  
  return xx * s;
}

static void
q_TW_D_writer_scaled(const int p[], int c, int d, int re_im, double v, void *e)
{
  TW_D_env *env = e;
  int i = QDP_index_L(env->lat, p);
  QNc(QLA_D, _DiracFermion) *f = env->f;
  double s = env->s;
  
  v = v * s;
  if (re_im == 0) {
    QLA_real(QLA_elem_D(f[i], c, d)) = v;
  } else {
    QLA_imag(QLA_elem_D(f[i], c, d)) = v;
  }
}

static int
q_dirac_solver(lua_State *L)
{
  int relaxed_p = 0, log_p = 0;
  if (qlua_checkopt_table(L,2)) {
    relaxed_p = qlua_tabkey_boolopt(L, 2, "relaxed", 0);
    log_p = qlua_tabkey_boolopt(L, 2, "log", 0);
  }
  int log_level = 0;
  if (log_p)
    log_level = (QOP_TWISTED_LOG_CG_RESIDUAL |
		 QOP_TWISTED_LOG_EIG_POSTAMBLE |
		 QOP_TWISTED_LOG_EIG_UPDATE1);

  TwistedSolver *solver = lua_touserdata(L, lua_upvalueindex(1));
  mTwisted *c = qlua_checkTwisted(L, lua_upvalueindex(2), NULL, 1);
  mLattice *S = qlua_ObjLattice(L, lua_upvalueindex(2));
  int Sidx = lua_gettop(L);
  long long fl1;
  double t1;
  double out_eps;
  int out_iters;
  

  QNz(mLatDirFerm) *psi = QNz(qlua_checkLatDirFerm)(L, 1, S, QLUA_TWISTED_NC);
  QNz(mLatDirFerm) *eta = QNz(qlua_newZeroLatDirFerm)(L, Sidx, QLUA_TWISTED_NC);
  struct QNc(QOP_,_TWISTED_Fermion) *c_psi;
  struct QNc(QOP_,_TWISTED_Fermion) *c_eta;
  TW_D_env env;
  QLA_D_Real rhs_norm2 = 0;
  double rhs_n;
  int status;
  const char *err_str;
  
  CALL_QDP(L);
  QNc(QDP_D,_r_eq_norm2_D)(&rhs_norm2, psi->ptr, S->all);
  if (rhs_norm2 == 0) {
    lua_pushnumber(L, 0.0);
    lua_pushnumber(L, 0);
    lua_pushnumber(L, 0);
    lua_pushnumber(L, 0);
    return 5;
  }
  rhs_n = sqrt(rhs_norm2);
  env.lat = S->lat;
  env.f = QNc(QDP_D, _expose_D)(psi->ptr);
  env.s = 1 / rhs_n;
  if (QNc(QOP_,_TWISTED_import_fermion)(&c_psi, c->state, q_TW_D_reader_scaled,
					&env))
    return luaL_error(L, "TWISTED_import_fermion() failed");
  QNc(QDP_D, _reset_D)(psi->ptr);
  
  if (QNc(QOP_, _TWISTED_allocate_fermion)(&c_eta, c->state))
    return luaL_error(L, "TWISTED_allocate_fermion() failed");
  
  status = solver->proc(L, c_eta, &out_iters, &out_eps, c_psi, log_level);
  
  if (status)
    err_str = QNc(QOP_, _TWISTED_error)(c->state);

  QNc(QOP_, _TWISTED_performance)(&t1, &fl1, NULL, NULL, c->state);
  if (t1 == 0)
    t1 = -1;
  if (QDP_this_node == qlua_master_node)
    printf("TWISTED %s solver: status = %d,"
	   " eps = %.4e, iters = %d, time = %.3f sec,"
	   " perf = %.2f MFlops/sec\n",
	   solver->name, status,
	   out_eps, out_iters, t1, fl1 * 1e-6 / t1);
  
  env.lat = S->lat;
  env.f = QNc(QDP_D, _expose_D)(eta->ptr);
  env.s = rhs_n;
  QNc(QOP_, _TWISTED_export_fermion)(q_TW_D_writer_scaled, &env, c_eta);
  QNc(QDP_D, _reset_D)(eta->ptr);
  
  QNc(QOP_, _TWISTED_free_fermion)(&c_eta);
  QNc(QOP_, _TWISTED_free_fermion)(&c_psi);
  
  if (status && !relaxed_p)
    return luaL_error(L, QNc(QOP_, _TWISTED_error)(c->state));
  
  /* eta is on the stack already */
  lua_pushnumber(L, out_eps);
  lua_pushnumber(L, out_iters);
  lua_pushnumber(L, t1);
  lua_pushnumber(L, (double)fl1);

  return 5;
}

static int
q_mxm_solver(lua_State *L)
{
  int relaxed_p = 0, log_p = 0;
  if (qlua_checkopt_table(L, 2)) {
    relaxed_p = qlua_tabkey_boolopt(L, 2, "relaxed", 0);
    log_p = qlua_tabkey_boolopt(L, 2, "log", 0);
  }
  int log_level = 0;
  if (log_p)
    log_level = (QOP_TWISTED_LOG_CG_RESIDUAL |
		 QOP_TWISTED_LOG_EIG_POSTAMBLE |
		 QOP_TWISTED_LOG_EIG_UPDATE1);

  TwistedMxMSolver *solver = lua_touserdata(L, lua_upvalueindex(1));
  mTwisted *c = qlua_checkTwisted(L, lua_upvalueindex(2), NULL, 1);
  mLattice *S = qlua_ObjLattice(L, lua_upvalueindex(2));
  int Sidx = lua_gettop(L);
  long long fl1;
  double t1;
  double out_eps;
  int out_iters;
  

  QNz(mLatDirFerm) *psi = QNz(qlua_checkLatDirFerm)(L, 1, S, QLUA_TWISTED_NC);
  QNz(mLatDirFerm) *eta = QNz(qlua_newZeroLatDirFerm)(L, Sidx, QLUA_TWISTED_NC);
  struct QNc(QOP_,_TWISTED_HalfFermion) *c_psi;
  struct QNc(QOP_,_TWISTED_HalfFermion) *c_eta;
  TW_D_env env;
  QLA_D_Real rhs_norm2 = 0;
  double rhs_n;
  int status;
  const char *err_str;
  
  CALL_QDP(L);
  QNc(QDP_D,_r_eq_norm2_D)(&rhs_norm2, psi->ptr, S->all);
  if (rhs_norm2 == 0) {
    lua_pushnumber(L, 0.0);
    lua_pushnumber(L, 0);
    lua_pushnumber(L, 0);
    lua_pushnumber(L, 0);
    return 5;
  }
  rhs_n = sqrt(rhs_norm2);
  env.lat = S->lat;
  env.f = QNc(QDP_D, _expose_D)(psi->ptr);
  env.s = 1 / rhs_n;
  if (QNc(QOP_,_TWISTED_import_half_fermion)(&c_psi, c->state, q_TW_D_reader_scaled,
                                             &env))
    return luaL_error(L, "TWISTED_import_fermion() failed");
  QNc(QDP_D, _reset_D)(psi->ptr);
  
  if (QNc(QOP_, _TWISTED_allocate_half_fermion)(&c_eta, c->state))
    return luaL_error(L, "TWISTED_allocate_fermion() failed");
  
  status = solver->proc(L, c_eta, &out_iters, &out_eps, c_psi, log_level);

  if (status)
    err_str = QNc(QOP_, _TWISTED_error)(c->state);
  
  QNc(QOP_, _TWISTED_performance)(&t1, &fl1, NULL, NULL, c->state);
  if (t1 == 0)
    t1 = -1;
  if (QDP_this_node == qlua_master_node)
    printf("TWISTED %s solver: status = %d,"
	   " eps = %.4e, iters = %d, time = %.3f sec,"
	   " perf = %.2f MFlops/sec\n",
	   solver->name, status,
	   out_eps, out_iters, t1, fl1 * 1e-6 / t1);
  
  env.lat = S->lat;
  env.f = QNc(QDP_D, _expose_D)(eta->ptr);
  env.s = rhs_n;
  QNc(QOP_, _TWISTED_export_half_fermion)(q_TW_D_writer_scaled, &env, c_eta);
  QNc(QDP_D, _reset_D)(eta->ptr);
  
  QNc(QOP_, _TWISTED_free_half_fermion)(&c_eta);
  QNc(QOP_, _TWISTED_free_half_fermion)(&c_psi);
  
  if (status && !relaxed_p)
    return luaL_error(L, QNc(QOP_, _TWISTED_error)(c->state));
  
  /* eta is on the stack already */
  lua_pushnumber(L, out_eps);
  lua_pushnumber(L, out_iters);
  lua_pushnumber(L, t1);
  lua_pushnumber(L, (double)fl1);
  
  return 5;
}

/***** deflator state interface */
static mDeflatorState *q_checkDeflatorState(lua_State *L,
                                            int idx,
                                            mLattice *S,
                                            int live);

static int
q_TWS_gc(lua_State *L)
{
    mDeflatorState *d = q_checkDeflatorState(L, 1, NULL, 0);

    if (d->deflator)
      QNc(QOP_F, _TWISTED_free_deflator)(&d->deflator);
    d->deflator = 0;

    return 0;
}

static struct luaL_Reg mtTwistedDeflatorState[] = {
    { "__gc",         q_TWS_gc },
    { NULL, NULL},
};

static mDeflatorState *
q_newDeflatorState(lua_State *L, int Sidx)
{
  mDeflatorState *d= lua_newuserdata(L, sizeof (mDeflatorState));
  d->nev = 0;
  d->vmax = 0;
  d->umax = 0;
  d->deflator = 0;
  qlua_createLatticeTable(L, Sidx, mtTwistedDeflatorState, qTwistedDeflatorState,
			  TwistedDeflatorStateName);
  lua_setmetatable(L, -2);
  
  return d;
}

static mDeflatorState*
q_checkDeflatorState(lua_State *L, int idx, mLattice *S, int live)
{
  mDeflatorState *d = qlua_checkLatticeType(L, idx, qTwistedDeflatorState,
					    TwistedDeflatorStateName);
  
  if (S) {
    mLattice *S1 = qlua_ObjLattice(L, idx);
    if (S1->id != S->id)
      luaL_error(L, "%s on a wrong lattice", TwistedDeflatorStateName);
    lua_pop(L, 1);
  }
  
  if (live && (d->deflator == 0))
    luaL_error(L, "Using closed Twisted.DeflatorState");
  
  return d;
}

/***** delfator interface */
static mTwisted *
q_Deflator_get_Twisted(lua_State *L, int idx, mLattice *S, int live)
{
  mTwisted *c;
  qlua_checktable(L, idx, "");
  
  lua_rawgeti(L, idx, 1);
  c = qlua_checkTwisted(L, -1, S, live);
  
  return c;
}

static mDeflatorState *
q_Deflator_get_State(lua_State *L, int idx, mLattice *S, int live)
{
  mDeflatorState *d;
  qlua_checktable(L, idx, "");

  lua_rawgeti(L, idx, 2);
  d = q_checkDeflatorState(L, -1, S, live);

  return d;
}

static int
q_TWDF_fmt(lua_State *L)
{
    mDeflatorState *d = q_Deflator_get_State(L, 1, NULL, 0);
    char fmt[72];

    if (d->deflator)
        sprintf(fmt, "Twisted.Deflator(0x%p)", d->deflator);
    else
        sprintf(fmt, "Twisted.Deflator(closed)");

    lua_pushstring(L, fmt);

    return 1;
}

static int
q_TWDF_close(lua_State *L)
{
  mDeflatorState *d = q_Deflator_get_State(L, 1, NULL, 1);

  QNc(QOP_F, _TWISTED_free_deflator)(&d->deflator);

  return 0;
}

static int
q_TWDF_reset(lua_State *L)
{
  mDeflatorState *d = q_Deflator_get_State(L, 1, NULL, 1);

  QNc(QOP_F, _TWISTED_deflator_eigcg_reset)(d->deflator);

  return 0;
}

static int
q_TWDF_stop(lua_State *L)
{
  mDeflatorState *d = q_Deflator_get_State(L, 1, NULL, 1);
  
  QNc(QOP_F, _TWISTED_deflator_eigcg_stop)(d->deflator);
  
  return 0;
}

static int
q_TWDF_resume(lua_State *L)
{
  mDeflatorState *d = q_Deflator_get_State(L, 1, NULL, 1);

  QNc(QOP_F, _TWISTED_deflator_eigcg_resume)(d->deflator);

  return 0;
}

static int
q_TWDF_eigenvalues(lua_State *L)
{
  mDeflatorState *d = q_Deflator_get_State(L, 1, NULL, 1);
  int df_cur_dim = QNc(QOP_F, _TWISTED_deflator_current_dim)(d->deflator);
  mVecReal *v = qlua_newVecReal(L, df_cur_dim);
  double *t = qlua_malloc(L, df_cur_dim * sizeof (double));

  CALL_QDP(L);
  int status = QNc(QOP_F, _TWISTED_deflator_eigen)(df_cur_dim, t, d->deflator);

  if (status == 0) {
    int i;
    for (i = 0 ; i < df_cur_dim ; i++)
      v->val[i] = t[i];
  }

  qlua_free(L, t);

  if (status == 0)
    return 1;
  else
    return 0;
}

static int
q_TWDF_current_dim(lua_State *L)
{
  mDeflatorState *d = q_Deflator_get_State(L, 1, NULL, 1);
  lua_pushnumber(L, QNc(QOP_F, _TWISTED_deflator_current_dim)(d->deflator));
  return 1;
}

static int
q_TWDF_start_load(lua_State *L)
{
  mDeflatorState *d = q_Deflator_get_State(L, 1, NULL, 1);
  if (QNc(QOP_F, _TWISTED_deflator_start_load)(d->deflator))
    return luaL_error(L, "TWISTED_deflator_start_load() failed");
  lua_pushnumber(L, QNc(QOP_F, _TWISTED_deflator_current_dim)(d->deflator));
  return 1;
}

static int
q_TWDF_stop_load(lua_State *L)
{
  mDeflatorState *d = q_Deflator_get_State(L, 1, NULL, 1);
  if (QNc(QOP_F, _TWISTED_deflator_stop_load)(d->deflator))
    return luaL_error(L, "TWISTED_deflator_start_load() failed");
  lua_pushnumber(L, QNc(QOP_F, _TWISTED_deflator_current_dim)(d->deflator));
  return 1;
}

static int
q_TWDF_add_vector(lua_State *L)
{
  mDeflatorState *d = q_Deflator_get_State(L, 1, NULL, 1); /* push1 */
  mTwisted *c = q_Deflator_get_Twisted(L, 1, NULL, 1); /* push1 */
  mLattice *S = qlua_ObjLattice(L, -1);
  QNz(mLatDirFerm) *psi = QNz(qlua_checkLatDirFerm)(L, 2, S, QLUA_TWISTED_NC);
  int n_vec = QNc(QOP_F, _TWISTED_deflator_current_dim)(d->deflator);
  double norm, norm2;
  QNc(QLA_D, _DiracFermion) *e_psi;
  struct QNc(QOP_F,_TWISTED_HalfFermion) *c_psi;
  struct QNc(QOP_F,_TWISTED_Gauge) *gaugeF;

  if (d->umax <= n_vec)
    return luaL_error(L, "vector space is full");

  CALL_QDP(L);
  QNc(QDP_D, _r_eq_norm2_D)(&norm2, psi->ptr, S->all);
  norm = sqrt(norm2);

  e_psi = QNc(QDP_D, _expose_D)(psi->ptr);

  TW_D_env    r_env;
  r_env.lat   = S->lat;
  r_env.f     = e_psi;
  r_env.s     = 1. / norm;
  if (QNc(QOP_F, _TWISTED_import_half_fermion)(&c_psi, c->state, q_TW_D_reader_scaled, &r_env))
    return luaL_error(L, "TWISTED_import_fermion() failed");
  if (QNc(QOP_, _TWISTED_gauge_float_from_double)(&gaugeF, c->gauge))
    return luaL_error(L, "TWISTED_gauge_float_from_double() failed");
  if (QNc(QOP_F, _TWISTED_deflator_add_vector)(gaugeF, d->deflator, c_psi))
    return luaL_error(L, QNc(QOP_, _TWISTED_error)(c->state));
  
  /* cleanup */
  QNc(QOP_F, _TWISTED_free_gauge)(&gaugeF);
  QNc(QOP_F, _TWISTED_free_half_fermion)(&c_psi);

  QNc(QDP_D, _reset_D)(psi->ptr);

  /* normal return */
  lua_pushnumber(L, QNc(QOP_F, _TWISTED_deflator_current_dim)(d->deflator));
  return 1;
}

static int
q_TWDF_get_vector(lua_State *L)
{
  mDeflatorState *d = q_Deflator_get_State(L, 1, NULL, 1);
  mTwisted *c = q_Deflator_get_Twisted(L, 1, NULL, 1);
  mLattice *S = qlua_ObjLattice(L, -1);
  int Sidx = lua_gettop(L);
  int num_vec = QNc(QOP_F, _TWISTED_deflator_current_dim)(d->deflator);
  int idx_vec = qlua_checkint(L, 2, "expect vector index");

  if (idx_vec < 0 || num_vec <= idx_vec)
    return luaL_error(L, "expect vector index 0 <= i < dim");
  
  QNz(mLatDirFerm) *psi = QNz(qlua_newZeroLatDirFerm)(L, Sidx, QLUA_TWISTED_NC);
  
  CALL_QDP(L);
  QNc(QLA_D, _DiracFermion) *e_psi = QNc(QDP_D,_expose_D)(psi->ptr);

  TW_D_env w_env;
  w_env.lat = S->lat;
  w_env.f   = e_psi;
  w_env.s   = 1.;

  struct QNc(QOP_F, _TWISTED_HalfFermion) *c_psi;
  if (QNc(QOP_F,_TWISTED_allocate_half_fermion)(&c_psi, c->state))
    return luaL_error(L, "TWISTED_allocate_half_fermion() failed");
  if (QNc(QOP_F, _TWISTED_deflator_extract_vector)(c_psi, d->deflator, idx_vec))
    return luaL_error(L, QNc(QOP_, _TWISTED_error)(c->state));
  if (QNc(QOP_F, _TWISTED_export_half_fermion)(q_TW_D_writer_scaled, &w_env, c_psi))
    return luaL_error(L, "TWISTED_export_half_fermion() failed");
  
  QNc(QDP_D, _reset_D)(psi->ptr);
  
  /* clean up */
  QNc(QOP_F, _TWISTED_free_half_fermion)(&c_psi);
  return 1;
}

static int
q_TWDF_deflated_mixed_solver(lua_State *L,
			     struct QNc(QOP_, _TWISTED_Fermion) *solution,
			     int *out_iters,
			     double *out_epsilon,
			     const struct QNc(QOP_, _TWISTED_Fermion) *rhs,
			     int log_level)
{
  mTwisted        *c = qlua_checkTwisted(L, lua_upvalueindex(2), NULL, 1);
  mLattice       *S = qlua_ObjLattice(L, lua_upvalueindex(2));
  mDeflatorState *d = q_checkDeflatorState(L, lua_upvalueindex(3), S, 1);
  double      f_eps = luaL_checknumber(L, lua_upvalueindex(4));
  int   inner_iters = luaL_checkint(L, lua_upvalueindex(5));
  double        eps = luaL_checknumber(L, lua_upvalueindex(6));
  int     max_iters = luaL_checkint(L, lua_upvalueindex(7));
  
  lua_pop(L, 1);
  return QNc(QOP_, _TWISTED_deflated_mixed_D_CG)(solution, out_iters, out_epsilon,
						 rhs, c->gauge, rhs, d->deflator,
						 inner_iters, f_eps,
						 max_iters, eps,
						 log_level);
}

static TwistedSolver deflated_mixed_solver = {
    q_TWDF_deflated_mixed_solver, "eigCG"
};

static int
q_TWDF_make_mixed_solver(lua_State *L)
{
  double inner_eps = qlua_tabkey_double(L, 2, "inner_eps");
  int inner_iter = qlua_tabkey_int(L, 2, "inner_iter");
  double eps = qlua_tabkey_double(L, 2, "eps");
  int max_iter = qlua_tabkey_int(L, 2, "max_iter");
  
  lua_pushlightuserdata(L, &deflated_mixed_solver); /* clo[1]: solver */
  q_Deflator_get_Twisted(L, 1, NULL, 1);            /* clo[2]: Twisted */
  q_Deflator_get_State(L, 1, NULL, 1);              /* clo[3]: Deflator */
  lua_pushnumber(L, inner_eps);                     /* clo[4]: inner_eps */
  lua_pushnumber(L, inner_iter);                    /* clo[5]: inner_iter */
  lua_pushnumber(L, eps);                           /* clo[6]: epsilon */
  lua_pushnumber(L, max_iter);                      /* clo[7]: max_iter */
  lua_pushcclosure(L, q_dirac_solver, 7);

  return 1;
}

static struct luaL_Reg mtDeflator[] = {
    { "__tostring",         q_TWDF_fmt               },
    { "__newindex",         qlua_nowrite             },
    { "mixed_solver",       q_TWDF_make_mixed_solver },
    { "close",              q_TWDF_close             },
    { "reset",              q_TWDF_reset             },
    { "stop",               q_TWDF_stop              },
    { "resume",             q_TWDF_resume            },
    { "eigenvalues",        q_TWDF_eigenvalues       },
    { "current_dim",        q_TWDF_current_dim       },
    { "start_load",         q_TWDF_start_load        },
    { "stop_load",          q_TWDF_stop_load         },
    { "add_vector",         q_TWDF_add_vector        },
    { "get_vector",         q_TWDF_get_vector        },
    { NULL,                 NULL                     }
};

static int
q_TW_make_deflator(lua_State *L)
{
  int vmax = qlua_tabkey_int(L, 2, "Vmax");
  int nev = qlua_tabkey_int(L, 2, "Nev");
  double eps = qlua_tabkey_double(L, 2, "eps");
  int umax = qlua_tabkey_int(L, 2, "Umax");

  mTwisted *c = qlua_checkTwisted(L, 1, NULL, 1);
  qlua_ObjLattice(L, 1);
  int Sidx = lua_gettop(L);

  if (nev <= 0 || vmax <= 0 || umax <= 0)
    return luaL_error(L, "bad eigenspace size");
  if (2 * nev >= vmax)
    return luaL_error(L, "eigcg VMAX: must satisfy VMAX > 2*NEV");
  
  if ((c->state == 0) || (c->gauge == 0))
    return luaL_error(L, "closed Twisted used");
  
  lua_createtable(L, 2, 0);
  lua_pushvalue(L, 1);
  lua_rawseti(L, -2, 1);
  mDeflatorState *d = q_newDeflatorState(L, Sidx);
  d->nev = nev;
  d->vmax = vmax;
  d->umax = umax;
  
  CALL_QDP(L);
  if (QNc(QOP_F, _TWISTED_create_deflator)(&d->deflator, c->state,
					   vmax, nev, eps, umax))
    return luaL_error(L, "TWISTED_create_deflator() failed");
  
  lua_rawseti(L, -2, 2);
  qlua_createLatticeTable(L, Sidx, mtDeflator, qTwistedDeflator,
			  TwistedDeflatorName);
  lua_setmetatable(L, -2);
  
  return 1;
}

#ifdef HAS_ARPACK
/* operator for lanczos : function and oblique arg */
typedef struct {
  struct QNc(QOP_, _TWISTED_State)      *twisted_state;
  struct QNc(QOP_F, _TWISTED_Gauge)     *twisted_gauge;

  /* workspace: must be allocated before calling twisted_eoprec_op */
  struct QNc(QOP_F,_TWISTED_HalfFermion) *x, *y;

  /* polynomial acc parameters :
     compute orthogonal polynomial of 'poly_n'-degree of the Op
     using 3-term recursion
     p_{i+1}(Op).x := [(poly_a[i] + poly_b[i]*Op) . p_{i}(Op)
     + poly_c[i]*p_{i-1}(Op)] . x
     where n = 0 .. {poly_n-1}, p_0(Op).x := x, p_{-1}.x := 0
  */
  int        poly_n; /* degree */
  double    *poly_a,
            *poly_b,
            *poly_c;
} QNc(op_TWISTED_F,_eoprec_MdagM_arg_s);

void
QNc(op_TWISTED_F, _eoprec_MdagM_op)(int loc_dim,
                                   float complex *x,
                                   float complex *y,
                                   void *op_arg) /* x<-op(y) */
{
  long long fl1, fl2;
  double t1, t2;

  QNc(op_TWISTED_F, _eoprec_MdagM_arg_s) *a = (QNc(op_TWISTED_F, _eoprec_MdagM_arg_s) *)op_arg;
  QLUA_ASSERT(2 * loc_dim == QNc(QOP_, _TWISTED_half_fermion_size)(a->twisted_state));

  QNc(QOP_F, _TWISTED_half_fermion_from_blas)(a->y, (float *)y, 2 * loc_dim);
  if (0 < a->poly_n) {
    QNc(QOP_F, _TWISTED_MxM_poly)(a->x, NULL, a->twisted_gauge, a->y,
                                 a->poly_n, a->poly_a, a->poly_b, a->poly_c);
    QNc(QOP_, _TWISTED_performance)(&t1, &fl1, NULL, NULL, a->twisted_state);

    QNc(QOP_F, _TWISTED_blas_from_half_fermion)((float *)x, 2 * loc_dim, a->x);
  } else {
    QNc(QOP_F, _TWISTED_M_operator)(a->x, a->twisted_gauge, a->y);
    QNc(QOP_, _TWISTED_performance)(&t1, &fl1, NULL, NULL, a->twisted_state);

    QNc(QOP_F, _TWISTED_M_operator_conjugated)(a->y, a->twisted_gauge, a->x);
    QNc(QOP_, _TWISTED_performance)(&t2, &fl2, NULL, NULL, a->twisted_state);
    t1  += t2;
    fl1 += fl2;

    QNc(QOP_F, _TWISTED_blas_from_half_fermion)((float *)x, 2 * loc_dim, a->y);
  }

  if (t1 == 0)
    t1 = -1;
  if (QDP_this_node == qlua_master_node)
    printf("TWISTED MdagM(poly_n=%d): time = %.3f sec,"
           " perf = %.2f MFlops/sec\n",
           (0 < a->poly_n ? a->poly_n : 1),
           t1, fl1 * 1e-6 / t1);
}

/* double precision operator */
typedef struct {
  struct QNc(QOP_, _TWISTED_State)           *twisted_state;
  struct QNc(QOP_D, _TWISTED_Gauge)        *twisted_gauge;

  /* workspace: must be allocated before calling twisted_eoprec_op */
  struct QNc(QOP_D, _TWISTED_HalfFermion)  *x, *y;

  /* polynomial acc parameters :
     compute orthogonal polynomial of 'poly_n'-degree of the Op
       using 3-term recursion
       p_{i+1}(Op).x := [(poly_a[i] + poly_b[i]*Op) . p_{i}(Op)
                      + poly_c[i]*p_{i-1}(Op)] . x
       where n = 0 .. {poly_n-1}, p_0(Op).x := x, p_{-1}.x := 0
     */
  int      poly_n; /* degree */
  double  *poly_a,
          *poly_b,
          *poly_c;
} QNc(op_TWISTED_D, _eoprec_MdagM_arg_s);


void
QNc(op_TWISTED_F, _eoprec_MdagM_double_op)(int loc_dim,
                                          float complex *x,
                                          float complex *y,
                                          void *op_arg) /* x<-op(y) */
{
  long long fl1, fl2;
  double t1, t2;
  int i;
  double complex *d_buf = NULL;

  QNc(op_TWISTED_D, _eoprec_MdagM_arg_s) *a = (QNc(op_TWISTED_D, _eoprec_MdagM_arg_s) *)op_arg;
  QLUA_ASSERT(2 * loc_dim == QNc(QOP_, _TWISTED_half_fermion_size)(a->twisted_state));

  d_buf = malloc(sizeof(d_buf[0]) * loc_dim);
  QLUA_ASSERT(NULL != d_buf);

  for (i = 0 ; i < loc_dim ; i++)
    d_buf[i] = y[i];
  QNc(QOP_D, _TWISTED_half_fermion_from_blas)(a->y, (double *)d_buf, 2 * loc_dim);

  if (0 < a->poly_n) {
    QNc(QOP_D, _TWISTED_MxM_poly)(a->x, NULL, a->twisted_gauge, a->y,
                                 a->poly_n, a->poly_a, a->poly_b, a->poly_c);
    QNc(QOP_, _TWISTED_performance)(&t1, &fl1, NULL, NULL, a->twisted_state);

    QNc(QOP_D, _TWISTED_blas_from_half_fermion)((double *)d_buf, 2 * loc_dim, a->x);
    for (i = 0 ; i < loc_dim ; i++)
      x[i] = d_buf[i];
  } else {
    QNc(QOP_D, _TWISTED_M_operator)(a->x, a->twisted_gauge, a->y);
    QNc(QOP_, _TWISTED_performance)(&t1, &fl1, NULL, NULL, a->twisted_state);

    QNc(QOP_D, _TWISTED_M_operator_conjugated)(a->y, a->twisted_gauge, a->x);
    QNc(QOP_, _TWISTED_performance)(&t2, &fl2, NULL, NULL, a->twisted_state);
    t1  += t2;
    fl1 += fl2;

    QNc(QOP_D, _TWISTED_blas_from_half_fermion)((double *)d_buf, 2 * loc_dim, a->y);
    for (i = 0 ; i < loc_dim ; i++)
      x[i] = d_buf[i];
  }

  free(d_buf);

  if (t1 == 0)
    t1 = -1;
  if (QDP_this_node == qlua_master_node)
    printf("TWISTED MdagM(poly_n=%d): time = %.3f sec,"
           " perf = %.2f MFlops/sec\n",
           (0 < a->poly_n ? a->poly_n : 1),
           t1, fl1 * 1e-6 / t1);
}

/* TWISTED:eig_deflator_lanczos(nev, ncv, max_iter, tol, [param])
   param = {
   [cheb_accel]={n, a, b},
   [eigcg] = {vmax, nev, eps, umax}
   ... }
   return: deflator object, n_converged, n_iter

   notes:
   * first, Lanczos/Arnoldi is used to compute eigenvectors;
   * converged vectors are saved into a newly created EigCG deflator,
   which is set to frozen state;
   * only min(eigcg_umax, #converged) vectors will be saved;
   * the user may specify eigcg_umax, eigcg_vmax, eigcg_nev values
   to continue EigCG after Lanczos (which is perhaps pointless though);
   * if any of (eigcg_umax, eigcg_vmax, eigcg_nev) are set to invalid
   values, the code will figure out the best values to continue operation
*/

static int
q_TW_make_deflator_lanczos(lua_State *L)
{

#define LANCZOS_MXM_DOUBLE  1
  /* by default, search for ev with smallest real part */
  const char *lanczos_which= "SR";
  const char *arpack_logfile = NULL;
  struct QNc(QOP_F, _TWISTED_HalfFermionMat) *hfm = NULL;
  /* operator parameters, init to empty */
  struct QNc(QOP_F, _TWISTED_Gauge) *gaugeF = NULL;
#ifdef LANCZOS_MXM_DOUBLE
  QNc(op_TWISTED_D, _eoprec_MdagM_arg_s) op_arg;
#else
  QNc(op_TWISTED_F, _eoprec_MdagM_arg_s) op_arg;
#endif/*LANCZOS_MXM_DOUBLE*/
  op_arg.poly_n = -1;
  op_arg.poly_a = op_arg.poly_b = op_arg.poly_c = NULL;
  op_arg.x = op_arg.y = NULL;

  int status;
  int do_lanczos_inplace = 0;
  int loc_dim;
  int n_iters, nconv;
  int n_evecs;
  int i;
  int Sidx;
  float complex *evec = NULL,
    *eval = NULL;

  mTwisted *c = qlua_checkTwisted(L, 1, NULL, 1);
  if (NULL == c->state || NULL == c->gauge)
    return luaL_error(L, "closed TWISTED used");

  mDeflatorState *d = NULL;

  /* parse parameters */
  /* XXX convert parameters to a single table */
  int nev      = qlua_checkint(L, 2, "expect NEV at #2");
  int ncv      = qlua_checkint(L, 3, "expect NCV at #3");
  int max_iter = qlua_checkint(L, 4, "expect MAX_ITER at #4");
  double tol   = luaL_checknumber(L, 5);

  /* parse optional parameters */
  int eigcg_vmax  = 0;
  int eigcg_umax  = 0;
  int eigcg_nev   = 0;
  double eigcg_eps= 0.;

  if (qlua_checkopt_table(L, 6)) {
    if (qlua_tabpushopt_key(L, 6, "cheb_accel")) {
      /* Chebyshev acceleration parameters */
#ifdef LANCZOS_MXM_DOUBLE
      QNc(op_TWISTED_D, _eoprec_MdagM_arg_s) *a = &op_arg;
#else
      QNc(op_TWISTED_F, _eoprec_MdagM_arg_s) *a = &op_arg;
#endif
      int cheb_n = -1;
      double cheb_a = 0.,
        cheb_b = 0.,
        cheb_x0= 1.;
      int do_norm = 0;

      if (0 <= op_arg.poly_n)
        return luaL_error(L, "more than one poly.accel. parameter");

      cheb_n  = qlua_tabidx_int(L, -1, 1);
      if (cheb_n < 0)
        return luaL_error(L, "poly.degree must be positive");

      cheb_a  = qlua_tabidx_double(L, -1, 2);
      cheb_b  = qlua_tabidx_double(L, -1, 3);
      if (cheb_a == cheb_b)
        return luaL_error(L, "invalid segment [a;b]");

      if (qlua_tabpushopt_idx(L, -1, 4)) {
        do_norm = 1;
        cheb_x0 = luaL_checknumber(L, -1);
        lua_pop(L, 1);
      }

      /* set up poly parameters for P_n(x) = T_n(x'), with
         x:[a;b] -> x':[-1;+1] */
      a->poly_a = a->poly_b = a->poly_c = NULL;
      a->poly_a = qlua_malloc(L, cheb_n * sizeof(a->poly_a[0]));
      a->poly_b = qlua_malloc(L, cheb_n * sizeof(a->poly_b[0]));
      a->poly_c = qlua_malloc(L, cheb_n * sizeof(a->poly_c[0]));
      if (NULL == a->poly_a
          || NULL == a->poly_b
          || NULL == a->poly_c)
        return luaL_error(L, "not enough memory");

      a->poly_a[0] = (-cheb_b - cheb_a) / (cheb_b - cheb_a);
      a->poly_b[0] = 2. / (cheb_b - cheb_a);
      a->poly_c[0] = 1;
      for (i = 1; i < cheb_n ; i++) {
        a->poly_a[i] = 2 * a->poly_a[0];
        a->poly_b[i] = 2 * a->poly_b[0];
        a->poly_c[i] = -1;
      }

      if (do_norm)
        QNc(QOP_, _TWISTED_poly_normalize)(cheb_n,
                                          a->poly_a, a->poly_b, a->poly_c,
                                          cheb_x0, 1e-8);

      a->poly_n = cheb_n;

      lua_pop(L, 1);
    }

    if (qlua_tabpushopt_key(L, 6, "which")) {
      lanczos_which = luaL_checkstring(L, -1);
      if (NULL == lanczos_which ||
          (  strcmp("SR", lanczos_which)
             && strcmp("LR", lanczos_which)
             && strcmp("SI", lanczos_which)
             && strcmp("LI", lanczos_which)
             && strcmp("SM", lanczos_which)
             && strcmp("LM", lanczos_which)))
        return luaL_error(L, "invalid value for which='%s'",
                          NULL == lanczos_which ? "null" : lanczos_which);
      lua_pop(L, 1);
    }

    if (qlua_tabpushopt_key(L, 6, "eigcg")) {
      eigcg_vmax  = qlua_tabidx_int(L, -1, 1);
      eigcg_nev   = qlua_tabidx_int(L, -1, 2);
      eigcg_eps   = qlua_tabidx_double(L, -1, 3);
      eigcg_umax  = qlua_tabidx_int(L, -1, 4);
      lua_pop(L, 1);
    }

    if (qlua_tabpushopt_key(L, 6, "arpack_logfile")) {
      arpack_logfile = luaL_checkstring(L, -1);
      lua_pop(L, 1);
    }

    if (qlua_tabpushopt_key(L, 6, "inplace")) {
      if (lua_isboolean(L, -1)) {
        do_lanczos_inplace = lua_toboolean(L, -1);
      } else
        return luaL_error(L, "'inplace' : expect boolean");
      lua_pop(L, 1);
    }
  }

#if 0    // no auto-setting of eigcg
  if (eigcg_vmax <= 0)
    eigcg_vmax = ncv;
  if (eigcg_nev <= 0) {
    eigcg_nev = 2; /* quite arbitrary, honestly */
    if (eigcg_vmax < 2 * eigcg_nev)
      eigcg_nev = eigcg_vmax / 2;
  }
  if (eigcg_vmax < 2 * eigcg_nev)
    return luaL_error(L, "eigcg VMAX: must satisfy VMAX > 2*NEV");
#endif


  /* create Qlua deflator object {TWISTED, DeflatorState} and set its META */
  /* FIXME should be refactored into a separate function?;
     this code duplicates parts of q_TW_make_deflator(qtwisted.c) */
  qlua_ObjLattice(L, 1);
  Sidx = lua_gettop(L);
  lua_createtable(L, 2, 0);
  lua_pushvalue(L, 1);
  lua_rawseti(L, -2, 1);
  if (NULL == (d = q_newDeflatorState(L, Sidx)))
    return luaL_error(L, "cannot create deflator state");
  lua_rawseti(L, -2, 2);
  qlua_createLatticeTable(L, Sidx, mtDeflator, qTwistedDeflator,
                          TwistedDeflatorName);
  lua_setmetatable(L, -2);



  CALL_QDP(L);

  gaugeF = NULL;
  if (QNc(QOP_, _TWISTED_gauge_float_from_double)(&gaugeF, c->gauge))
    return luaL_error(L, "QOP_TWISTED_gauge_float_from_double() failed");

#ifdef LANCZOS_MXM_DOUBLE
  op_arg.twisted_state = c->state;
  op_arg.twisted_gauge = c->gauge;

  op_arg.x = op_arg.y = NULL;
  if (QNc(QOP_D, _TWISTED_allocate_half_fermion)(&op_arg.x, c->state)
      || QNc(QOP_D, _TWISTED_allocate_half_fermion)(&op_arg.y, c->state))
    return luaL_error(L, "cannot allocate HalfFermion");
#else
  op_arg.twisted_state   = c->state;
  if (QNc(QOP_, _TWISTED_gauge_float_from_double)(&(op_arg.twisted_gauge), c->gauge))
    return luaL_error(L, "TWISTED_gauge_float_from_double() failed");

  op_arg.x = op_arg.y = NULL;
  if (QNc(QOP_F, _TWISTED_allocate_half_fermion)(&op_arg.x, c->state)
      || QNc(QOP_F, _TWISTED_allocate_half_fermion)(&op_arg.y, c->state))
    return luaL_error(L, "cannot allocate HalfFermion");
#endif/*LANCZOS_MXM_DOUBLE*/

  MPI_Comm mpi_comm = MPI_COMM_WORLD; /* FIXME any better choice? */

  loc_dim = QNc(QOP_, _TWISTED_half_fermion_size)(c->state) / 2;
  QLUA_ASSERT(0 == QNc(QOP_, _TWISTED_half_fermion_size)(c->state) % 2);

  /* run Arnoldi/Lanczos iterations */
  n_iters = nconv = 0;
  evec = eval = NULL;
  if (do_lanczos_inplace) {
    int inplace_umax = ncv < eigcg_umax ? eigcg_umax : ncv;
    float *hfm_blas_ptr = NULL;
    int hfm_nrow_loc = 0,
      hfm_ncol     = 0,
      hfm_ld       = 0;

    if (NULL == (eval = qlua_malloc(L, sizeof(eval[0]) * nev)))
      return luaL_error(L, "not enough memory");
    if (QNc(QOP_F, _TWISTED_alloc_half_fermion_matrix)(&hfm, c->state, inplace_umax))
      return luaL_error(L, "TWISTED_alloc_half_fermion_matrix failed");
    if (QNc(QOP_F, _TWISTED_blas_view_half_fermion_matrix)(hfm, &hfm_nrow_loc, &hfm_ncol,
                                                          &hfm_blas_ptr, &hfm_ld))
      return luaL_error(L, QNc(QOP_, _TWISTED_error)(c->state));

    if (0 != (status = lanczos_inplace_float(
                                             L, mpi_comm,
#ifdef LANCZOS_MXM_DOUBLE
                                             QNc(op_TWISTED_F, _eoprec_MdagM_double_op),
#else
                                             QNc(op_TWISTED_F, _eoprec_MdagM_op),
#endif/*LANCZOS_MXM_DOUBLE*/
                                             &op_arg,
                                             lanczos_which, loc_dim, nev, ncv, max_iter, tol,
                                             NULL, /* initial vector - pending implementation */
                                             eval, (float complex *)hfm_blas_ptr, hfm_ld, hfm_ncol,
                                             &n_iters, &nconv, arpack_logfile))) {
      return luaL_error(L, "lanczos_float_inplace returned %d", status);
    }
    if (QNc(QOP_F, _TWISTED_create_deflator_inplace)(&(d->deflator), gaugeF, &hfm,
                                                    nconv, eigcg_vmax, eigcg_nev, eigcg_eps,
                                                    eigcg_umax))
      return luaL_error(L, QNc(QOP_, _TWISTED_error)(c->state));

  } else {
    if (0 != (status = lanczos_float(
                                     L, mpi_comm,
#ifdef LANCZOS_MXM_DOUBLE
                                     QNc(op_TWISTED_F, _eoprec_MdagM_double_op),
#else
                                     QNc(op_TWISTED_F, _eoprec_MdagM_op),
#endif/*LANCZOS_MXM_DOUBLE*/
                                     &op_arg,
                                     lanczos_which, loc_dim, nev, ncv, max_iter, tol,
                                     NULL, /* initial vector - pending implementation */
                                     &eval, &evec, &n_iters, &nconv, arpack_logfile)))
      return luaL_error(L, "lanczos_float returned %d", status);
    /* FIXME rewrite with clear explanation for the choice of
       default values */
    if (eigcg_umax <= 0 || eigcg_umax <= nconv)
      eigcg_umax = nconv;


    if (QNc(QOP_F, _TWISTED_create_deflator)(&d->deflator, c->state,
                                            eigcg_vmax, eigcg_nev, eigcg_eps, eigcg_umax))
      return luaL_error(L, "TWISTED_create_deflator() failed");

    /* fill deflator with e.vecs */
    if (QNc(QOP_F, _TWISTED_deflator_start_load)(d->deflator))
      return luaL_error(L, "TWISTED_deflator_start_load() failed");
    /* (have space for eigcg_umax) */
    n_evecs = (nconv <= eigcg_umax ? nconv : eigcg_umax);

    struct QNc(QOP_F, _TWISTED_HalfFermion) *hf_buf = NULL;
    if (QNc(QOP_F, _TWISTED_allocate_half_fermion)(&hf_buf, c->state))
      return luaL_error(L, "cannot allocate HalfFermion");

    for (i = 0 ; i < n_evecs ; i++) {
      QNc(QOP_F, _TWISTED_half_fermion_from_blas)(hf_buf,
                                                 (float *)(evec + i * loc_dim), 2 * loc_dim);
      if (QNc(QOP_F, _TWISTED_deflator_add_vector)(gaugeF, d->deflator, hf_buf)) {
        QNc(QOP_F, _TWISTED_free_half_fermion)(&hf_buf);
        return luaL_error(L, "TWISTED_deflator_add_vector() failed");
      }
    }

    if (QNc(QOP_F, _TWISTED_deflator_stop_load)(d->deflator))
      return luaL_error(L, "TWISTED_deflator_end_load() failed");
  }

  /* initialize TWISTED side of deflator */
  d->nev  = eigcg_nev;
  d->vmax = eigcg_vmax;
  d->umax = eigcg_umax;

  CALL_QDP(L);

  QNc(QOP_F,_TWISTED_deflator_eigcg_stop)(d->deflator);

  /* cleanup */
  if (NULL != evec) free(evec);
  if (NULL != eval) free(eval);
  if (NULL != gaugeF) QNc(QOP_F, _TWISTED_free_gauge)(&(gaugeF));
#ifdef LANCZOS_MXM_DOUBLE
  if (NULL != op_arg.y) QNc(QOP_D, _TWISTED_free_half_fermion)(&op_arg.y);
  if (NULL != op_arg.x) QNc(QOP_D, _TWISTED_free_half_fermion)(&op_arg.x);
#else
  if (NULL != op_arg.twisted_gauge) QNc(QOP_F, _TWISTED_free_gauge)(&(op_arg.twisted_gauge));
  if (NULL != op_arg.y) QNc(QOP_F, _TWISTED_free_half_fermion)(&op_arg.y);
  if (NULL != op_arg.x) QNc(QOP_F, _TWISTED_free_half_fermion)(&op_arg.x);
#endif/*LANCZOS_MXM_DOUBLE*/

  if (NULL != op_arg.poly_a) qlua_free(L, op_arg.poly_a);
  if (NULL != op_arg.poly_b) qlua_free(L, op_arg.poly_b);
  if (NULL != op_arg.poly_c) qlua_free(L, op_arg.poly_c);

  if (NULL != hfm)  QNc(QOP_F,_TWISTED_free_half_fermion_matrix)(&hfm);

  lua_pushnumber(L, nconv);
  lua_pushnumber(L, n_iters);
  return 3; /* deflator object, n_converged, n_iter */
}

#endif /* HAS_ARPACK */

/***** twisted interface */
static int
q_TW_fmt(lua_State *L)
{
    char fmt[72];
    mTwisted *c = qlua_checkTwisted(L, 1, NULL, 0);

    if (c->state)
        sprintf(fmt, "Twisted[%p]", c);
    else
        sprintf(fmt, "Twisted(closed)");

    lua_pushstring(L, fmt);

    return 1;
}

static int
q_TW_gc(lua_State *L)
{
    mTwisted *c = qlua_checkTwisted(L, 1, NULL, 0);

    if (c->gauge) {
      QNc(QOP_, _TWISTED_free_gauge)(&c->gauge);
        c->gauge = 0;
    }
    if (c->state) {
      QNc(QOP_, _TWISTED_fini)(&c->state);
        c->state = 0;
    }

    return 0;
}

static int
q_TW_close(lua_State *L)
{
    mTwisted *c = qlua_checkTwisted(L, 1, NULL, 1);

    if (c->gauge)
      QNc(QOP_, _TWISTED_free_gauge)(&c->gauge);
    if (c->state)
      QNc(QOP_, _TWISTED_fini)(&c->state);
    c->gauge = 0;
    c->state = 0;

    return 0;
}

static double
q_TW_D_reader(const int p[], int c, int d, int re_im, void *e)
{
    TW_D_env *env = e;
    int i = QDP_index_L(env->lat, p);
    QNc(QLA_D, _DiracFermion) *f = env->f;
    QLA_D_Real xx;

    if (re_im == 0) {
        QLA_r_eq_Re_c(xx, QLA_elem_D(f[i], c, d));
    } else {
        QLA_r_eq_Im_c(xx, QLA_elem_D(f[i], c, d));
    }

    return xx;
}

static void
q_TW_D_writer(const int p[], int c, int d, int re_im, double v, void *e)
{
    TW_D_env *env = e;
    int i = QDP_index_L(env->lat, p);
    QNc(QLA_D, _DiracFermion) *f = env->f;

    if (re_im == 0) {
        QLA_real(QLA_elem_D(f[i], c, d)) = v;
    } else {
        QLA_imag(QLA_elem_D(f[i], c, d)) = v;
    }
}

static int
q_TW_operator(lua_State *L,
              const char *name,
              int (*op)(struct QNc(QOP_D, _TWISTED_Fermion) *result,
                        const struct QNc(QOP_D, _TWISTED_Gauge) *gauge,
                        const struct QNc(QOP_D, _TWISTED_Fermion) *source))
{
    mTwisted *c = qlua_checkTwisted(L, 1, NULL, 1);
    mLattice *S = qlua_ObjLattice(L, 1);
    int Sidx = lua_gettop(L);
    QNz(mLatDirFerm) *psi = QNz(qlua_checkLatDirFerm)(L, 2, S, QLUA_TWISTED_NC);
    QNz(mLatDirFerm) *eta = QNz(qlua_newLatDirFerm)(L, Sidx, QLUA_TWISTED_NC);
    struct QNc(QOP_, _TWISTED_Fermion) *c_psi;
    struct QNc(QOP_, _TWISTED_Fermion) *c_eta;
    QNc(QLA_D, _DiracFermion) *e_psi;
    QNc(QLA_D, _DiracFermion) *e_eta;
    TW_D_env env;
    
    CALL_QDP(L);
    e_psi = QNc(QDP_D, _expose_D)(psi->ptr);
    env.lat = S->lat;
    env.f = e_psi;
    if (QNc(QOP_, _TWISTED_import_fermion)(&c_psi, c->state, q_TW_D_reader, &env))
      return luaL_error(L, "TWISTED_import_fermion() failed");
    QNc(QDP_D, _reset_D)(psi->ptr);
    
    if (QNc(QOP_, _TWISTED_allocate_fermion)(&c_eta, c->state))
      return luaL_error(L, "TWISTED_allocate_fermion() failed");
    
    (*op)(c_eta, c->gauge, c_psi);
    
    e_eta = QNc(QDP_D, _expose_D)(eta->ptr);
    env.f = e_eta;
    QNc(QOP_, _TWISTED_export_fermion)(q_TW_D_writer, &env, c_eta);
    QNc(QDP_D, _reset_D)(eta->ptr);
    
    QNc(QOP_, _TWISTED_free_fermion)(&c_eta);
    QNc(QOP_, _TWISTED_free_fermion)(&c_psi);
    
    return 1;
}

static int
q_TW_D(lua_State *L)
{
  return q_TW_operator(L, "D", QNc(QOP_D, _TWISTED_D_operator));
}

static int
q_TW_Dx(lua_State *L)
{
  return q_TW_operator(L, "Dx", QNc(QOP_D, _TWISTED_D_operator_conjugated));
}

/* the standard twisted solver */
static int
q_TW_std_solver(lua_State *L,
                struct QNc(QOP_, _TWISTED_Fermion) *solution,
                int *out_iters,
                double *out_epsilon,
                const struct QNc(QOP_, _TWISTED_Fermion) *rhs,
                int log_level)
{
    mTwisted *c = qlua_checkTwisted(L, lua_upvalueindex(2), NULL, 1);
    double eps = luaL_checknumber(L, lua_upvalueindex(3));
    int max_iters = luaL_checkint(L, lua_upvalueindex(4));

    return QNc(QOP_, _TWISTED_D_CG)(solution, out_iters, out_epsilon,
                           rhs, c->gauge, rhs, max_iters, eps,
                           log_level);
}

static TwistedSolver std_solver = { q_TW_std_solver, "CG" };

static int
q_TW_make_solver(lua_State *L)
{
  qlua_checkTwisted(L, 1, NULL, 1);   /* mTwisted */
  double eps = qlua_tabkey_double(L, 2, "eps");
  int max_iter = qlua_tabkey_int(L, 2, "max_iter");
  
  lua_pushlightuserdata(L, &std_solver);  /* cl[1]: solver */
  lua_pushvalue(L, 1);                    /* cl[2]: mTwisted */
  lua_pushnumber(L, eps);                 /* cl[3]: epsilon */
  lua_pushnumber(L, max_iter);            /* cl[4]: max_iters */
  lua_pushcclosure(L, q_dirac_solver, 4);
  
  return 1;
}

/* the standard twisted MxM solver */
static int
q_TW_mxm_std_solver(lua_State *L,
                struct QNc(QOP_, _TWISTED_HalfFermion) *solution,
                int *out_iters,
                double *out_epsilon,
                const struct QNc(QOP_, _TWISTED_HalfFermion) *rhs,
                int log_level)
{
    mTwisted *c = qlua_checkTwisted(L, lua_upvalueindex(2), NULL, 1);
    double eps = luaL_checknumber(L, lua_upvalueindex(3));
    int max_iters = luaL_checkint(L, lua_upvalueindex(4));

    return QNc(QOP_, _TWISTED_MxM_CG)(solution, out_iters, out_epsilon,
                           rhs, c->gauge, rhs, max_iters, eps,
                           log_level);
}

static TwistedMxMSolver mxm_std_solver = { q_TW_mxm_std_solver, "MxM_CG" };

static int
q_TW_make_mxm_solver(lua_State *L)
{
  qlua_checkTwisted(L, 1, NULL, 1);   /* mTwisted */
  double eps = qlua_tabkey_double(L, 2, "eps");
  int max_iter = qlua_tabkey_int(L, 2, "max_iter");
  
  lua_pushlightuserdata(L, &mxm_std_solver);  /* cl[1]: solver */
  lua_pushvalue(L, 1);                        /* cl[2]: mTwisted */
  lua_pushnumber(L, eps);                     /* cl[3]: epsilon */
  lua_pushnumber(L, max_iter);                /* cl[4]: max_iters */
  lua_pushcclosure(L, q_mxm_solver, 4);
  
  return 1;
}

/* the mixed twisted solver */
static int
q_TW_mixed_solver(lua_State *L,
                  struct QNc(QOP_, _TWISTED_Fermion) *solution,
                  int *out_iters,
                  double *out_epsilon,
                  const struct QNc(QOP_, _TWISTED_Fermion) *rhs,
                  int log_level)
{
    mTwisted *c = qlua_checkTwisted(L, lua_upvalueindex(2), NULL, 1);
    double f_eps    = luaL_checknumber(L, lua_upvalueindex(3));
    int inner_iters = luaL_checkint(L, lua_upvalueindex(4));
    double eps      = luaL_checknumber(L, lua_upvalueindex(5));
    int max_iters   = luaL_checkint(L, lua_upvalueindex(6));

    return QNc(QOP_, _TWISTED_mixed_D_CG)(solution, out_iters, out_epsilon,
                                 rhs, c->gauge, rhs,
                                 inner_iters, f_eps,
                                 max_iters, eps,
                                 log_level);
}

static TwistedSolver mixed_solver = { q_TW_mixed_solver, "mixedCG" };

static int
q_TW_make_mixed_solver(lua_State *L)
{
  qlua_checkTwisted(L, 1, NULL, 1);   /* mTwisted */
  double inner_eps = qlua_tabkey_double(L, 2, "inner_eps");
  int inner_iter = qlua_tabkey_int(L, 2, "inner_iter");
  double eps = qlua_tabkey_double(L, 2, "eps");
  int max_iter = qlua_tabkey_int(L, 2, "max_iter");

  lua_pushlightuserdata(L, &mixed_solver);  /* cl[1]: solver */
  lua_pushvalue(L, 1);                      /* cl[2]: mTwisted */
  lua_pushnumber(L, inner_eps);             /* cl[3]: inner_epsilon */
  lua_pushnumber(L, inner_iter);            /* cl[4]: inner_iters */
  lua_pushnumber(L, eps);                   /* cl[5]: epsilon */
  lua_pushnumber(L, max_iter);              /* cl[6]: max_iters */
  lua_pushcclosure(L, q_dirac_solver, 6);

  return 1;
}

#define Nu  QNc(QOP_, _TWISTED_DIM)

typedef struct {
  QDP_Lattice *lat;
  int lattice[QNc(QOP_, _TWISTED_DIM)];
  QLA_D_Complex mq;
  QLA_D_Complex mu;
  QLA_D_Complex bf[QNc(QOP_, _TWISTED_DIM)];
  QNc(QLA_D, _ColorMatrix) *uf[Nu];
} QTWArgs;

static double
q_TW_u_reader(int d, const int p[], int a, int b, int re_im, void *env)
{
    QLA_D_Complex z;
    QTWArgs *args = env;
    int i = QDP_index_L(args->lat, p);

    if (p[d] == (args->lattice[d] - 1)) {
        QLA_c_eq_c_times_c(z, args->bf[d], QLA_elem_M(args->uf[d][i], a, b));
    } else {
        QLA_c_eq_c(z, QLA_elem_M(args->uf[d][i], a, b));
    }

    if (re_im == 0)
        return QLA_real(z);
    else
        return QLA_imag(z);
}

static int
build_twisted(lua_State *L,
	      mLattice *S,
	      double mq_re,
	      double mq_im,
	      double mu_re,
	      double mu_im,
	      mTwisted *twisted,
	      QNc(QLA_D, _ColorMatrix) *uf[],
	      double (*u_reader)(int, const int [], int, int, int, void *),
	      void *args)
{
  QNc(QDP_D, _ColorMatrix) *UF[Nu];
  
  luaL_checktype(L, 1, LUA_TTABLE);
  CALL_QDP(L);
  
  /* create a temporary F */
  int i;
  for (i = QNc(QOP_, _TWISTED_DIM); i < Nu; i++)
    UF[i] = QNc(QDP_D, _create_M_L)(S->lat);
  
  /* extract U from the arguments */
  for (i = 0; i < QNc(QOP_, _TWISTED_DIM); i++) {
    lua_pushnumber(L, i + 1); /* [sic] lua indexing */
    lua_gettable(L, 1);
    UF[i] = QNz(qlua_checkLatColMat)(L, -1, S, QLUA_TWISTED_NC)->ptr;
    lua_pop(L, 1);
  }
  
  struct QNc(QOP_, _TWISTED_Config) cc;
  cc.self = S->node;
  cc.master_p = QMP_is_primary_node();
  cc.rank = S->rank;
  cc.lat = S->dim;
  cc.net = S->net;
  cc.neighbor_up = S->neighbor_up;
  cc.neighbor_down = S->neighbor_down;
  cc.sublattice = qlua_sublattice;
  cc.env = S;
  if (QNc(QOP_, _TWISTED_init)(&twisted->state, &cc))
    return luaL_error(L, "TWISTED_init() failed");
  
  /* import the gauge field */
  for (i = 0; i < Nu; i++) {
    uf[i] = QNc(QDP_D, _expose_M)(UF[i]);
  }
  
  if (QNc(QOP_, _TWISTED_import_gauge)(&twisted->gauge, twisted->state,
				       mq_re, mq_im, mu_re, mu_im, u_reader, args)) {
    return luaL_error(L, "TWISTED_import_gauge() failed");
  }

  for (i = 0; i < Nu; i++)
    QNc(QDP_D, _reset_M)(UF[i]);

    return 1;
}

static void
start_twisted(lua_State *L, mLattice **ptr_S, mTwisted **ptr_twisted)
{
  luaL_checktype(L, 1, LUA_TTABLE);
  lua_pushnumber(L, 1);
  lua_gettable(L, 1);
  QNz(qlua_checkLatColMat)(L, -1, NULL, QLUA_TWISTED_NC);
  *ptr_S = qlua_ObjLattice(L, -1);
  int Sidx = lua_gettop(L);
  *ptr_twisted = qlua_newTwisted(L, Sidx);
  if ((*ptr_S)->rank != QNc(QOP_, _TWISTED_DIM))
    luaL_error(L, "twisted is not implemented for #L=%d", (*ptr_S)->rank);
  if (QDP_Ns != QNc(QOP_,  _TWISTED_FERMION_DIM))
    luaL_error(L, "twisted does not support Ns=%d", QDP_Ns);
}

/*
 * qcd.Twisted(U,                 -- gauge field, { Ux,Uy,Uz,Ut }
 *             {mq = number,      -- complex quark mass
 *              mu = number,      -- complex twist
 *              boundary = { bx, by, bz, bt } })  -- complex boundary conditions
 */

static int
q_twisted(lua_State *L)
{
  mLattice     *S = NULL;
  mTwisted      *twisted = NULL;
  QTWArgs       args;
  double        mq_re, mq_im, mu_re, mu_im;

  qlua_tabkey_complex(L, 2, "mq", &mq_re, &mq_im, "mq");
  qlua_tabkey_complex(L, 2, "mu", &mu_re, &mu_im, "mu");
  if (qlua_tabkey_tableopt(L, 2, "boundary") == 0)
    luaL_error(L, "missing boundary values");
  qlua_get_complex_vector(L, lua_gettop(L),  QNc(QOP_, _TWISTED_DIM), args.bf, "boundary");
  lua_pop(L, 1);

  start_twisted(L, &S, &twisted);
  args.lat = S->lat;
  QDP_latsize_L(S->lat, args.lattice);

  return build_twisted(L, S, mq_re, mq_im, mu_re, mu_im, twisted, args.uf, q_TW_u_reader, &args);
}


static struct luaL_Reg mtTwisted[] = {
    { "__tostring",                 q_TW_fmt },
    { "__gc",                       q_TW_gc },
    { "close",                      q_TW_close },
    { "D",                          q_TW_D },
    { "Dx",                         q_TW_Dx },
    { "solver",                     q_TW_make_solver },
    { "mxm_solver",                 q_TW_make_mxm_solver },
    { "mixed_solver",               q_TW_make_mixed_solver },
    { "eig_deflator",               q_TW_make_deflator },
#ifdef HAS_ARPACK
    { "eig_deflator_lanczos",       q_TW_make_deflator_lanczos      },
#endif /* HAS_ARPACK */
    { NULL, NULL }
};

static mTwisted *
qlua_newTwisted(lua_State *L, int Sidx)
{
    mTwisted *c = lua_newuserdata(L, sizeof (mTwisted));

    c->state = 0;
    c->gauge = 0;
    qlua_createLatticeTable(L, Sidx, mtTwisted, qTwisted, TwistedName);
    lua_setmetatable(L, -2);

    return c;
}

static mTwisted *
qlua_checkTwisted(lua_State *L, int idx, mLattice *S, int live)
{
    mTwisted *c = qlua_checkLatticeType(L, idx, qTwisted, TwistedName);

    if (S) {
        mLattice *S1 = qlua_ObjLattice(L, idx);
        if (S1->id != S->id)
            luaL_error(L, "%s on a wrong lattice", TwistedName);
        lua_pop(L, 1);
    }

    if (live && (c->state == 0 || c->gauge == 0))
        luaL_error(L, "using closed qcd.Twisted");

    return c;
}

static struct luaL_Reg fTwisted[] = {
    { "Twisted",                 q_twisted },
    { NULL,                      NULL }
};

int
init_twisted(lua_State *L)
{
    luaL_register(L, qcdlib, fTwisted);
    return 0;
}

void
fini_twisted(void)
{
}
