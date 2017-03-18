#include "modules.h"                                                 /* DEPS */
#include "qlua.h"                                                    /* DEPS */
#include "qclover.h"                                                 /* DEPS */
#include "qcomplex.h"                                                /* DEPS */
#include "qvector.h"                                                 /* DEPS */
#include "lattice.h"                                                 /* DEPS */
#include "qlayout.h"                                                 /* DEPS */
#include "latcolmat.h"                                               /* DEPS */
#include "latdirferm.h"                                              /* DEPS */
#include "latdirprop.h"                                              /* DEPS */
#define QOP_CLOVER_DEFAULT_PRECISION QDP_Precision
#include QLUA_CLOVER_HDR
#include "qmp.h"
#include "qlanczos.h"

#include <math.h>
#include <string.h>

#define QNc(a,b) QNc_1(a,QLUA_CLOVER_NC,b)
#define QNc_1(a,b,c) QNc_2(a,b,c)
#define QNc_2(a,b,c) a ## b ## c
#define QNz(a) QNz_1(a,QLUA_CLOVER_NC)
#define QNz_1(a,b) QNz_2(a,b)
#define QNz_2(a,b) a ## b

#if QLUA_CLOVER_NC != QNc(QOP_, _CLOVER_COLORS)
#error "Clover Nc mismatch"
#endif /*  QLUA_CLOVER_NC != QNc(QOP_, _CLOVER_COLORS) */


static const char CloverName[]               = "lattice.Clover";
static const char CloverDeflatorName[]       = "lattice.Clover.Deflator";
static const char CloverDeflatorStateName[]  = "lattice.Clover.DeflatorState";

typedef int CloverInverter(lua_State *L,
                           struct QNc(QOP_,_CLOVER_Fermion) *solution,
                           int *out_iters,
                           double *out_epsilon,
                           const struct QNc(QOP_,_CLOVER_Fermion) *rhs,
                           int log_level);

typedef struct {
    CloverInverter  *proc;
    const char      *name;
} CloverSolver;

typedef int CloverMxMInverter(lua_State *L,
                           struct QNc(QOP_,_CLOVER_HalfFermion) *solution,
                           int *out_iters,
                           double *out_epsilon,
                           const struct QNc(QOP_,_CLOVER_HalfFermion) *rhs,
                           int log_level);

typedef struct {
    CloverMxMInverter  *proc;
    const char      *name;
} CloverMxMSolver;

typedef struct {
  struct QNc(QOP_, _CLOVER_State) *state;
  struct QNc(QOP_, _CLOVER_Gauge) *gauge;
} mClover;

typedef struct {
  int nev;
  int umax;
  int vmax;
  struct QNc(QOP_F, _CLOVER_Deflator) *deflator;
} mDeflatorState;

static mClover *qlua_newClover(lua_State *L, int Sidx);
static mClover *qlua_checkClover(lua_State *L,
                                 int idx,
                                 mLattice *S,
                                 int live);


/* The generic solver */
typedef struct {
  QDP_Lattice *lat;
  int c, d;
  QNc(QLA_D, _DiracPropagator) *in;
  QNc(QLA_D, _DiracPropagator) *out;
  double s;
} qCL_P_env;

typedef struct {
  QDP_Lattice *lat;
  QNc(QLA_D, _DiracFermion) *f;
  double s;
} CL_D_env;

static double
q_CL_D_reader_scaled(const int p[], int c, int d, int re_im, void *e)
{
    CL_D_env *env = e;
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
q_CL_D_writer_scaled(const int p[], int c, int d, int re_im, double v, void *e)
{
    CL_D_env *env = e;
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

static double
q_CL_P_reader_scaled(const int p[], int c, int d, int re_im, void *e)
{
    qCL_P_env *env = e;
    int i = QDP_index_L(env->lat, p);
    QLA_D_Real xx;

    if (re_im == 0) {
        QLA_r_eq_Re_c(xx, QLA_elem_P(env->in[i], c, d, env->c, env->d));
    } else {
        QLA_r_eq_Im_c(xx, QLA_elem_P(env->in[i], c, d, env->c, env->d));
    }

    return xx * env->s;
}

static void
q_CL_P_writer_scaled(const int p[], int c, int d, int re_im, double v, void *e)
{
    qCL_P_env *env = e;
    int i = QDP_index_L(env->lat, p);

    v = v * env->s;
    if (re_im == 0) {
        QLA_real(QLA_elem_P(env->out[i], c, d, env->c, env->d)) = v;
    } else {
        QLA_imag(QLA_elem_P(env->out[i], c, d, env->c, env->d)) = v;
    }
}

static int
q_dirac_solver(lua_State *L)
{
    CloverSolver *solver = lua_touserdata(L, lua_upvalueindex(1));
    mClover *c = qlua_checkClover(L, lua_upvalueindex(2), NULL, 1);
    mLattice *S = qlua_ObjLattice(L, lua_upvalueindex(2));
    int Sidx = lua_gettop(L);
    int relaxed_p;
    long long fl1;
    double t1;
    double out_eps;
    int out_iters;
    int log_level;

    switch (lua_type(L, 2)) {
    case LUA_TNONE:
    case LUA_TNIL:
        relaxed_p = 0;
        break;
    case LUA_TBOOLEAN:
        relaxed_p = lua_toboolean(L, 2);
        break;
    default:
        relaxed_p = 1;
        break;
    }

    if ((lua_type(L, 3) == LUA_TBOOLEAN) && (lua_toboolean(L, 3) != 0))
      log_level = (QOP_CLOVER_LOG_CG_RESIDUAL |
                   QOP_CLOVER_LOG_EIG_POSTAMBLE |
                   QOP_CLOVER_LOG_EIG_UPDATE1);
    else
        log_level = 0;

    switch (qlua_qtype(L, 1)) {
    case QNz(qLatDirFerm): {
        QNz(mLatDirFerm) *psi = QNz(qlua_checkLatDirFerm)(L, 1, S, QLUA_CLOVER_NC);
        QNz(mLatDirFerm) *eta = QNz(qlua_newZeroLatDirFerm)(L, Sidx, QLUA_CLOVER_NC);
        struct QNc(QOP_,_CLOVER_Fermion) *c_psi;
        struct QNc(QOP_,_CLOVER_Fermion) *c_eta;
        CL_D_env env;
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
        if (QNc(QOP_,_CLOVER_import_fermion)(&c_psi, c->state, q_CL_D_reader_scaled,
                                             &env))
            return luaL_error(L, "CLOVER_import_fermion() failed");
        QNc(QDP_D, _reset_D)(psi->ptr);

        if (QNc(QOP_, _CLOVER_allocate_fermion)(&c_eta, c->state))
            return luaL_error(L, "CLOVER_allocate_fermion() failed");

        status = solver->proc(L, c_eta, &out_iters, &out_eps, c_psi, log_level);

        if (status)
            err_str =  QNc(QOP_,_CLOVER_error)(c->state);

        QNc(QOP_, _CLOVER_performance)(&t1, &fl1, NULL, NULL, c->state);
        if (t1 == 0)
            t1 = -1;
        if (QDP_this_node == qlua_master_node)
            printf("CLOVER %s solver: status = %d,"
                   " eps = %.4e, iters = %d, time = %.3f sec,"
                   " perf = %.2f MFlops/sec\n",
                   solver->name, status,
                   out_eps, out_iters, t1, fl1 * 1e-6 / t1);

        env.lat = S->lat;
        env.f = QNc(QDP_D, _expose_D)(eta->ptr);
        env.s = rhs_n;
        QNc(QOP_, _CLOVER_export_fermion)(q_CL_D_writer_scaled, &env, c_eta);
        QNc(QDP_D, _reset_D)(eta->ptr);

        QNc(QOP_, _CLOVER_free_fermion)(&c_eta);
        QNc(QOP_, _CLOVER_free_fermion)(&c_psi);

        if (status && !relaxed_p)
            return luaL_error(L, err_str);

        /* eta is on the stack already */
        lua_pushnumber(L, out_eps);
        lua_pushnumber(L, out_iters);
        lua_pushnumber(L, t1);
        lua_pushnumber(L, (double)fl1);

        return 5;
    }

    case QNz(qLatDirProp): {
        QNz(mLatDirProp) *psi = QNz(qlua_checkLatDirProp)(L, 1, S, QLUA_CLOVER_NC);
        QNz(mLatDirProp) *eta = QNz(qlua_newZeroLatDirProp)(L, Sidx, QLUA_CLOVER_NC);
        struct QNc(QOP_, _CLOVER_Fermion) *c_psi;
        struct QNc(QOP_, _CLOVER_Fermion) *c_eta;
        int status;
        qCL_P_env env;
        QLA_D_Real rhs_norm2 = 0;
        double rhs_n;
        const char *err_str;

        lua_createtable(L, QNc(QOP_, _CLOVER_COLORS), 0);  /* eps */
        lua_createtable(L, QNc(QOP_, _CLOVER_COLORS), 0);  /* iters */
        CALL_QDP(L);
        if (QNc(QOP_, _CLOVER_allocate_fermion)(&c_eta, c->state))
            return luaL_error(L, "CLOVER_allocate_fermion() failed");

        QNc(QDP_D, _r_eq_norm2_P)(&rhs_norm2, psi->ptr, S->all);
        if (rhs_norm2 == 0) {
            return 3;
        }
        rhs_n = sqrt(rhs_norm2);
        env.lat = S->lat;
        env.in = QNc(QDP_D, _expose_P)(psi->ptr);
        env.out = QNc(QDP_D, _expose_P)(eta->ptr);
        for (env.c = 0; env.c < QNc(QOP_, _CLOVER_COLORS); env.c++) {
            lua_createtable(L, QNc(QOP_,  _CLOVER_FERMION_DIM), 0); /* eps.c */
            lua_createtable(L, QNc(QOP_,  _CLOVER_FERMION_DIM), 0); /* iters.c */

            for (env.d = 0; env.d < QNc(QOP_,  _CLOVER_FERMION_DIM); env.d++) {
                env.s = 1 / rhs_n;
                if (QNc(QOP_, _CLOVER_import_fermion)(&c_psi, c->state,
                                              q_CL_P_reader_scaled, &env))
                    return luaL_error(L, "CLOVER_import_fermion() failed");
                status = solver->proc(L, c_eta, &out_iters, &out_eps, c_psi,
                                      log_level);

                if (status)
                    err_str =  QNc(QOP_,_CLOVER_error)(c->state);

                QNc(QOP_, _CLOVER_performance)(&t1, &fl1, NULL, NULL, c->state);
                if (t1 == 0)
                    t1 = -1;
                if (QDP_this_node == qlua_master_node)
                    printf("CLOVER %s solver: status = %d, c = %d, d = %d,"
                           " eps = %.4e, iters = %d, time = %.3f sec,"
                           " perf = %.2f MFlops/sec\n",
                           solver->name, status,
                           env.c, env.d, out_eps, out_iters, t1,
                           fl1 * 1e-6 / t1);
                QNc(QOP_, _CLOVER_free_fermion)(&c_psi);
                if (status && !relaxed_p)
                    return luaL_error(L, err_str);

                env.s = rhs_n;
                QNc(QOP_, _CLOVER_export_fermion)(q_CL_P_writer_scaled, &env, c_eta);
                lua_pushnumber(L, out_eps);
                lua_rawseti(L, -3, env.d + 1);
                lua_pushnumber(L, out_iters);
                lua_rawseti(L, -2, env.d + 1);
            }
            lua_rawseti(L, -3, env.c + 1);
            lua_rawseti(L, -3, env.c + 1);
        }
        QNc(QDP_D, _reset_P)(psi->ptr);
        QNc(QDP_D, _reset_P)(eta->ptr);
        QNc(QOP_, _CLOVER_free_fermion)(&c_eta);

        return 3;
    }
    default:
        break;
    }
    return luaL_error(L, "bad argument to CLOVER solver");
}

static int
q_mxm_solver(lua_State *L)
{
    CloverMxMSolver *solver = lua_touserdata(L, lua_upvalueindex(1));
    mClover *c = qlua_checkClover(L, lua_upvalueindex(2), NULL, 1);
    mLattice *S = qlua_ObjLattice(L, lua_upvalueindex(2));
    int Sidx = lua_gettop(L);
    int relaxed_p;
    long long fl1;
    double t1;
    double out_eps;
    int out_iters;
    int log_level;

    switch (lua_type(L, 2)) {
    case LUA_TNONE:
    case LUA_TNIL:
        relaxed_p = 0;
        break;
    case LUA_TBOOLEAN:
        relaxed_p = lua_toboolean(L, 2);
        break;
    default:
        relaxed_p = 1;
        break;
    }

    if ((lua_type(L, 3) == LUA_TBOOLEAN) && (lua_toboolean(L, 3) != 0))
      log_level = (QOP_CLOVER_LOG_CG_RESIDUAL |
                   QOP_CLOVER_LOG_EIG_POSTAMBLE |
                   QOP_CLOVER_LOG_EIG_UPDATE1);
    else
        log_level = 0;

    switch (qlua_qtype(L, 1)) {
    case QNz(qLatDirFerm): {
        QNz(mLatDirFerm) *psi = QNz(qlua_checkLatDirFerm)(L, 1, S, QLUA_CLOVER_NC);
        QNz(mLatDirFerm) *eta = QNz(qlua_newZeroLatDirFerm)(L, Sidx, QLUA_CLOVER_NC);
        struct QNc(QOP_,_CLOVER_HalfFermion) *c_psi;
        struct QNc(QOP_,_CLOVER_HalfFermion) *c_eta;
        CL_D_env env;
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
        if (QNc(QOP_,_CLOVER_import_half_fermion)(&c_psi, c->state, q_CL_D_reader_scaled,
                                             &env))
            return luaL_error(L, "CLOVER_import_fermion() failed");
        QNc(QDP_D, _reset_D)(psi->ptr);

        if (QNc(QOP_, _CLOVER_allocate_half_fermion)(&c_eta, c->state))
            return luaL_error(L, "CLOVER_allocate_fermion() failed");

        status = solver->proc(L, c_eta, &out_iters, &out_eps, c_psi, log_level);

        if (status)
            err_str =  QNc(QOP_,_CLOVER_error)(c->state);

        QNc(QOP_, _CLOVER_performance)(&t1, &fl1, NULL, NULL, c->state);
        if (t1 == 0)
            t1 = -1;
        if (QDP_this_node == qlua_master_node)
            printf("CLOVER %s solver: status = %d,"
                   " eps = %.4e, iters = %d, time = %.3f sec,"
                   " perf = %.2f MFlops/sec\n",
                   solver->name, status,
                   out_eps, out_iters, t1, fl1 * 1e-6 / t1);

        env.lat = S->lat;
        env.f = QNc(QDP_D, _expose_D)(eta->ptr);
        env.s = rhs_n;
        QNc(QOP_, _CLOVER_export_half_fermion)(q_CL_D_writer_scaled, &env, c_eta);
        QNc(QDP_D, _reset_D)(eta->ptr);

        QNc(QOP_, _CLOVER_free_half_fermion)(&c_eta);
        QNc(QOP_, _CLOVER_free_half_fermion)(&c_psi);

        if (status && !relaxed_p)
            return luaL_error(L, err_str);

        /* eta is on the stack already */
        lua_pushnumber(L, out_eps);
        lua_pushnumber(L, out_iters);
        lua_pushnumber(L, t1);
        lua_pushnumber(L, (double)fl1);

        return 5;
    }

    case QNz(qLatDirProp): {
        QNz(mLatDirProp) *psi = QNz(qlua_checkLatDirProp)(L, 1, S, QLUA_CLOVER_NC);
        QNz(mLatDirProp) *eta = QNz(qlua_newZeroLatDirProp)(L, Sidx, QLUA_CLOVER_NC);
        struct QNc(QOP_, _CLOVER_HalfFermion) *c_psi;
        struct QNc(QOP_, _CLOVER_HalfFermion) *c_eta;
        int status;
        qCL_P_env env;
        QLA_D_Real rhs_norm2 = 0;
        double rhs_n;
        const char *err_str;

        lua_createtable(L, QNc(QOP_, _CLOVER_COLORS), 0);  /* eps */
        lua_createtable(L, QNc(QOP_, _CLOVER_COLORS), 0);  /* iters */
        CALL_QDP(L);
        if (QNc(QOP_, _CLOVER_allocate_half_fermion)(&c_eta, c->state))
            return luaL_error(L, "CLOVER_allocate_fermion() failed");

        QNc(QDP_D, _r_eq_norm2_P)(&rhs_norm2, psi->ptr, S->all);
        if (rhs_norm2 == 0) {
            return 3;
        }
        rhs_n = sqrt(rhs_norm2);
        env.lat = S->lat;
        env.in = QNc(QDP_D, _expose_P)(psi->ptr);
        env.out = QNc(QDP_D, _expose_P)(eta->ptr);
        for (env.c = 0; env.c < QNc(QOP_, _CLOVER_COLORS); env.c++) {
            lua_createtable(L, QNc(QOP_,  _CLOVER_FERMION_DIM), 0); /* eps.c */
            lua_createtable(L, QNc(QOP_,  _CLOVER_FERMION_DIM), 0); /* iters.c */

            for (env.d = 0; env.d < QNc(QOP_,  _CLOVER_FERMION_DIM); env.d++) {
                env.s = 1 / rhs_n;
                if (QNc(QOP_, _CLOVER_import_half_fermion)(&c_psi, c->state,
                                              q_CL_P_reader_scaled, &env))
                    return luaL_error(L, "CLOVER_import_fermion() failed");
                status = solver->proc(L, c_eta, &out_iters, &out_eps, c_psi,
                                      log_level);

                if (status)
                    err_str =  QNc(QOP_,_CLOVER_error)(c->state);

                QNc(QOP_, _CLOVER_performance)(&t1, &fl1, NULL, NULL, c->state);
                if (t1 == 0)
                    t1 = -1;
                if (QDP_this_node == qlua_master_node)
                    printf("CLOVER %s solver: status = %d, c = %d, d = %d,"
                           " eps = %.4e, iters = %d, time = %.3f sec,"
                           " perf = %.2f MFlops/sec\n",
                           solver->name, status,
                           env.c, env.d, out_eps, out_iters, t1,
                           fl1 * 1e-6 / t1);
                QNc(QOP_, _CLOVER_free_half_fermion)(&c_psi);
                if (status && !relaxed_p)
                    return luaL_error(L, err_str);

                env.s = rhs_n;
                QNc(QOP_, _CLOVER_export_half_fermion)(q_CL_P_writer_scaled, &env, c_eta);
                lua_pushnumber(L, out_eps);
                lua_rawseti(L, -3, env.d + 1);
                lua_pushnumber(L, out_iters);
                lua_rawseti(L, -2, env.d + 1);
            }
            lua_rawseti(L, -3, env.c + 1);
            lua_rawseti(L, -3, env.c + 1);
        }
        QNc(QDP_D, _reset_P)(psi->ptr);
        QNc(QDP_D, _reset_P)(eta->ptr);
        QNc(QOP_, _CLOVER_free_half_fermion)(&c_eta);

        return 3;
    }
    default:
        break;
    }
    return luaL_error(L, "bad argument to CLOVER solver");
}

/***** deflator state interface */
static mDeflatorState *q_checkDeflatorState(lua_State *L,
                                            int idx,
                                            mLattice *S,
                                            int live);

static int
q_DFS_gc(lua_State *L)
{
    mDeflatorState *d = q_checkDeflatorState(L, 1, NULL, 0);

    if (d->deflator)
      QNc(QOP_F, _CLOVER_free_deflator)(&d->deflator);
    d->deflator = 0;

    return 0;
}


static struct luaL_Reg mtCloverDeflatorState[] = {
    { "__gc",         q_DFS_gc },
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
    qlua_createLatticeTable(L, Sidx, mtCloverDeflatorState, qCloverDeflatorState,
                            CloverDeflatorStateName);
    lua_setmetatable(L, -2);

    return d;
}



static mDeflatorState*
q_checkDeflatorState(lua_State *L, int idx, mLattice *S, int live)
{
    mDeflatorState *d = qlua_checkLatticeType(L, idx, qCloverDeflatorState,
                                              CloverDeflatorStateName);

    if (S) {
        mLattice *S1 = qlua_ObjLattice(L, idx);
        if (S1->id != S->id)
            luaL_error(L, "%s on a wrong lattice", CloverDeflatorStateName);
        lua_pop(L, 1);
    }

    if (live && (d->deflator == 0))
        luaL_error(L, "Using closed Clover.DeflatorState");

    return d;
}



/***** delfator interface */
static mClover *
q_Deflator_get_Clover(lua_State *L, int idx, mLattice *S, int live)
{
    mClover *c;
    qlua_checktable(L, idx, "");

    lua_rawgeti(L, idx, 1);
    c = qlua_checkClover(L, -1, S, live);

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
q_DF_fmt(lua_State *L)
{
    mDeflatorState *d = q_Deflator_get_State(L, 1, NULL, 0);
    char fmt[72];

    if (d->deflator)
        sprintf(fmt, "Clover.Deflator(0x%p)", d->deflator);
    else
        sprintf(fmt, "Clover.Deflator(closed)");

    lua_pushstring(L, fmt);

    return 1;
}

static int
q_DF_close(lua_State *L)
{
    mDeflatorState *d = q_Deflator_get_State(L, 1, NULL, 1);

    QNc(QOP_F, _CLOVER_free_deflator)(&d->deflator);

    return 0;
}

static int
q_DF_reset(lua_State *L)
{
    mDeflatorState *d = q_Deflator_get_State(L, 1, NULL, 1);

    QNc(QOP_F, _CLOVER_deflator_eigcg_reset)(d->deflator);

    return 0;
}

static int
q_DF_stop(lua_State *L)
{
    mDeflatorState *d = q_Deflator_get_State(L, 1, NULL, 1);

    QNc(QOP_F, _CLOVER_deflator_eigcg_stop)(d->deflator);

    return 0;
}

static int
q_DF_resume(lua_State *L)
{
    mDeflatorState *d = q_Deflator_get_State(L, 1, NULL, 1);

    QNc(QOP_F, _CLOVER_deflator_eigcg_resume)(d->deflator);

    return 0;
}

static int
q_DF_eigenvalues(lua_State *L)
{
    mDeflatorState *d = q_Deflator_get_State(L, 1, NULL, 1);
    int df_cur_dim = QNc(QOP_F, _CLOVER_deflator_current_dim)(d->deflator);
    mVecReal *v = qlua_newVecReal(L, df_cur_dim);
    double *t = qlua_malloc(L, df_cur_dim * sizeof (double));

    CALL_QDP(L);
    int status = QNc(QOP_F, _CLOVER_deflator_eigen)(df_cur_dim, t, d->deflator);

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
q_DF_current_dim(lua_State *L)
{
    mDeflatorState *d = q_Deflator_get_State(L, 1, NULL, 1);
    lua_pushnumber(L, QNc(QOP_F, _CLOVER_deflator_current_dim)(d->deflator));
    return 1;
}

static int
q_DF_start_load(lua_State *L)
{
    mDeflatorState *d = q_Deflator_get_State(L, 1, NULL, 1);
    if (QNc(QOP_F, _CLOVER_deflator_start_load)(d->deflator))
        return luaL_error(L, "CLOVER_deflator_start_load() failed");
    lua_pushnumber(L, QNc(QOP_F, _CLOVER_deflator_current_dim)(d->deflator));
    return 1;
}

static int
q_DF_stop_load(lua_State *L)
{
    mDeflatorState *d = q_Deflator_get_State(L, 1, NULL, 1);
    if (QNc(QOP_F, _CLOVER_deflator_stop_load)(d->deflator))
        return luaL_error(L, "CLOVER_deflator_start_load() failed");
    lua_pushnumber(L, QNc(QOP_F, _CLOVER_deflator_current_dim)(d->deflator));
    return 1;
}

static int
q_DF_add_vector(lua_State *L)
{
    mDeflatorState *d = q_Deflator_get_State(L, 1, NULL, 1); /* push1 */
    mClover *c = q_Deflator_get_Clover(L, 1, NULL, 1); /* push1 */
    mLattice *S = qlua_ObjLattice(L, -1);
    QNz(mLatDirFerm) *psi = QNz(qlua_checkLatDirFerm)(L, 2, S, QLUA_CLOVER_NC);
    int n_vec = QNc(QOP_F, _CLOVER_deflator_current_dim)(d->deflator);
    double norm, norm2;
    QNc(QLA_D, _DiracFermion) *e_psi;
    struct QNc(QOP_F,_CLOVER_HalfFermion) *c_psi;
    struct QNc(QOP_F,_CLOVER_Gauge) *gaugeF;

    if (d->umax <= n_vec)
      return luaL_error(L, "vector space is full");

    CALL_QDP(L);
    QNc(QDP_D, _r_eq_norm2_D)(&norm2, psi->ptr, S->all);
    norm = sqrt(norm2);

    e_psi = QNc(QDP_D, _expose_D)(psi->ptr);

    CL_D_env    r_env;
    r_env.lat   = S->lat;
    r_env.f     = e_psi;
    r_env.s     = 1. / norm;
    if (QNc(QOP_F, _CLOVER_import_half_fermion)(&c_psi, c->state, q_CL_D_reader_scaled, &r_env))
      return luaL_error(L, "CLOVER_import_fermion() failed");
    if (QNc(QOP_, _CLOVER_gauge_float_from_double)(&gaugeF, c->gauge))
      return luaL_error(L, "CLOVER_gauge_float_from_double() failed");
    if (QNc(QOP_F, _CLOVER_deflator_add_vector)(gaugeF, d->deflator, c_psi))
      return luaL_error(L, QNc(QOP_, _CLOVER_error)(c->state));

    /* cleanup */
    QNc(QOP_F, _CLOVER_free_gauge)(&gaugeF);
    QNc(QOP_F, _CLOVER_free_half_fermion)(&c_psi);

    QNc(QDP_D, _reset_D)(psi->ptr);

    /* normal return */
    lua_pushnumber(L, QNc(QOP_F, _CLOVER_deflator_current_dim)(d->deflator));
    return 1;
}

static int
q_DF_get_vector(lua_State *L)
{
    mDeflatorState *d = q_Deflator_get_State(L, 1, NULL, 1);
    mClover *c = q_Deflator_get_Clover(L, 1, NULL, 1);
    mLattice *S = qlua_ObjLattice(L, -1);
    int Sidx = lua_gettop(L);
    int num_vec = QNc(QOP_F, _CLOVER_deflator_current_dim)(d->deflator);
    int idx_vec = qlua_checkint(L, 2, "expect vector index");

    if (idx_vec < 0 || num_vec <= idx_vec)
        return luaL_error(L, "expect vector index 0 <= i < dim");

    QNz(mLatDirFerm) *psi = QNz(qlua_newZeroLatDirFerm)(L, Sidx, QLUA_CLOVER_NC);

    CALL_QDP(L);
    QNc(QLA_D, _DiracFermion) *e_psi = QNc(QDP_D,_expose_D)(psi->ptr);

    CL_D_env w_env;
    w_env.lat = S->lat;
    w_env.f   = e_psi;
    w_env.s   = 1.;

    struct QNc(QOP_F, _CLOVER_HalfFermion) *c_psi;
    if (QNc(QOP_F,_CLOVER_allocate_half_fermion)(&c_psi, c->state))
      return luaL_error(L, "CLOVER_allocate_half_fermion() failed");
    if (QNc(QOP_F, _CLOVER_deflator_extract_vector)(c_psi, d->deflator, idx_vec))
      return luaL_error(L, QNc(QOP_, _CLOVER_error)(c->state));
    if (QNc(QOP_F, _CLOVER_export_half_fermion)(q_CL_D_writer_scaled, &w_env, c_psi))
      return luaL_error(L, "CLOVER_export_half_fermion() failed");

    QNc(QDP_D, _reset_D)(psi->ptr);

    /* clean up */
    QNc(QOP_F, _CLOVER_free_half_fermion)(&c_psi);
    return 1;
}

static int
q_DF_deflated_mixed_solver(lua_State *L,
                           struct QNc(QOP_, _CLOVER_Fermion) *solution,
                           int *out_iters,
                           double *out_epsilon,
                           const struct QNc(QOP_, _CLOVER_Fermion) *rhs,
                           int log_level)
{
    mClover        *c = qlua_checkClover(L, lua_upvalueindex(2), NULL, 1);
    mLattice       *S = qlua_ObjLattice(L, lua_upvalueindex(2));
    mDeflatorState *d = q_checkDeflatorState(L, lua_upvalueindex(3), S, 1);
    double      f_eps = luaL_checknumber(L, lua_upvalueindex(4));
    int   inner_iters = luaL_checkint(L, lua_upvalueindex(5));
    double        eps = luaL_checknumber(L, lua_upvalueindex(6));
    int     max_iters = luaL_checkint(L, lua_upvalueindex(7));

    lua_pop(L, 1);
    return QNc(QOP_, _CLOVER_deflated_mixed_D_CG)(solution, out_iters, out_epsilon,
                                          rhs, c->gauge, rhs, d->deflator,
                                          inner_iters, f_eps,
                                          max_iters, eps,
                                          log_level);
}

static CloverSolver deflated_mixed_solver = {
    q_DF_deflated_mixed_solver, "eigCG"
};

static int
q_DF_make_mixed_solver(lua_State *L)
{
    double inner_eps = luaL_checknumber(L, 2);
    int inner_iter = luaL_checkint(L, 3);
    double eps = luaL_checknumber(L, 4);
    int max_iter = luaL_checkint(L, 5);

    lua_pushlightuserdata(L, &deflated_mixed_solver); /* clo[1]: solver */
    q_Deflator_get_Clover(L, 1, NULL, 1);             /* clo[2]: Clover */
    q_Deflator_get_State(L, 1, NULL, 1);              /* clo[3]: Deflator */
    lua_pushnumber(L, inner_eps);                     /* clo[4]: inner_eps */
    lua_pushnumber(L, inner_iter);                    /* clo[5]: inner_iter */
    lua_pushnumber(L, eps);                           /* clo[6]: epsilon */
    lua_pushnumber(L, max_iter);                      /* clo[7]: max_iter */
    lua_pushcclosure(L, q_dirac_solver, 7);

    return 1;
}

static struct luaL_Reg mtDeflator[] = {
    { "__tostring",         q_DF_fmt               },
    { "__newindex",         qlua_nowrite           },
    { "mixed_solver",       q_DF_make_mixed_solver },
    { "close",              q_DF_close             },
    { "reset",              q_DF_reset             },
    { "stop",               q_DF_stop              },
    { "resume",             q_DF_resume            },
    { "eigenvalues",        q_DF_eigenvalues       },
    { "current_dim",        q_DF_current_dim       },
    { "start_load",         q_DF_start_load        },
    { "stop_load",          q_DF_stop_load         },
    { "add_vector",         q_DF_add_vector        },
    { "get_vector",         q_DF_get_vector        },
#if 0  /* XXX deflator extra methods */
    { "truncate",           q_DF_truncate          },
    { "get_counter",        q_DF_get_counter       },
    { "put_counter",        q_DF_put_counter       },
    { "write_eigspace",     q_DF_write             },
    { "read_eigspace",      q_DF_read              },
#endif  /* XXX deflator extra methods */
    { NULL,                 NULL                   }
};

static int
q_CL_make_deflator(lua_State *L)
{
    mClover *c = qlua_checkClover(L, 1, NULL, 1);
    qlua_ObjLattice(L, 1);
    int Sidx = lua_gettop(L);
    int vmax = luaL_checkint(L, 2);
    int nev = luaL_checkint(L, 3);
    double eps = luaL_checknumber(L, 4);
    int umax = luaL_checkint(L, 5);

    if (nev <= 0 || vmax <= 0 || umax <= 0)
        return luaL_error(L, "bad eigenspace size");
    if (2 * nev >= vmax)
        return luaL_error(L, "eigcg VMAX: must satisfy VMAX > 2*NEV");

    if ((c->state == 0) || (c->gauge == 0))
        return luaL_error(L, "closed Clover used");

    lua_createtable(L, 2, 0);
    lua_pushvalue(L, 1);
    lua_rawseti(L, -2, 1);
    mDeflatorState *d = q_newDeflatorState(L, Sidx);
    d->nev = nev;
    d->vmax = vmax;
    d->umax = umax;

    CALL_QDP(L);
    if (QNc(QOP_F, _CLOVER_create_deflator)(&d->deflator, c->state,
                                   vmax, nev, eps, umax))
        return luaL_error(L, "CLOVER_create_deflator() failed");

    lua_rawseti(L, -2, 2);
    qlua_createLatticeTable(L, Sidx, mtDeflator, qCloverDeflator,
                            CloverDeflatorName);
    lua_setmetatable(L, -2);

    return 1;
}

#ifdef HAS_ARPACK
/* operator for lanczos : function and oblique arg */
typedef struct {
  struct QNc(QOP_, _CLOVER_State)      *clover_state;
  struct QNc(QOP_F, _CLOVER_Gauge)     *clover_gauge;

  /* workspace: must be allocated before calling clover_eoprec_op */
  struct QNc(QOP_F,_CLOVER_HalfFermion) *x, *y;

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
} QNc(op_CLOVER_F,_eoprec_MdagM_arg_s);

void
QNc(op_CLOVER_F, _eoprec_MdagM_op)(int loc_dim,
                                   float complex *x,
                                   float complex *y,
                                   void *op_arg) /* x<-op(y) */
{
  long long fl1, fl2;
  double t1, t2;

  QNc(op_CLOVER_F, _eoprec_MdagM_arg_s) *a = (QNc(op_CLOVER_F, _eoprec_MdagM_arg_s) *)op_arg;
  QLUA_ASSERT(2 * loc_dim == QNc(QOP_, _CLOVER_half_fermion_size)(a->clover_state));

  QNc(QOP_F, _CLOVER_half_fermion_from_blas)(a->y, (float *)y, 2 * loc_dim);
  if (0 < a->poly_n) {
    QNc(QOP_F, _CLOVER_MxM_poly)(a->x, NULL, a->clover_gauge, a->y,
                                 a->poly_n, a->poly_a, a->poly_b, a->poly_c);
    QNc(QOP_, _CLOVER_performance)(&t1, &fl1, NULL, NULL, a->clover_state);

    QNc(QOP_F, _CLOVER_blas_from_half_fermion)((float *)x, 2 * loc_dim, a->x);
  } else {
    QNc(QOP_F, _CLOVER_M_operator)(a->x, a->clover_gauge, a->y);
    QNc(QOP_, _CLOVER_performance)(&t1, &fl1, NULL, NULL, a->clover_state);

    QNc(QOP_F, _CLOVER_M_operator_conjugated)(a->y, a->clover_gauge, a->x);
    QNc(QOP_, _CLOVER_performance)(&t2, &fl2, NULL, NULL, a->clover_state);
    t1  += t2;
    fl1 += fl2;

    QNc(QOP_F, _CLOVER_blas_from_half_fermion)((float *)x, 2 * loc_dim, a->y);
  }

  if (t1 == 0)
    t1 = -1;
  if (QDP_this_node == qlua_master_node)
    printf("CLOVER MdagM(poly_n=%d): time = %.3f sec,"
           " perf = %.2f MFlops/sec\n",
           (0 < a->poly_n ? a->poly_n : 1),
           t1, fl1 * 1e-6 / t1);
}

/* double precision operator */
typedef struct {
  struct QNc(QOP_, _CLOVER_State)           *clover_state;
  struct QNc(QOP_D, _CLOVER_Gauge)        *clover_gauge;

  /* workspace: must be allocated before calling clover_eoprec_op */
  struct QNc(QOP_D, _CLOVER_HalfFermion)  *x, *y;

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
} QNc(op_CLOVER_D, _eoprec_MdagM_arg_s);



void
QNc(op_CLOVER_F, _eoprec_MdagM_double_op)(int loc_dim,
                                          float complex *x,
                                          float complex *y,
                                          void *op_arg) /* x<-op(y) */
{
  long long fl1, fl2;
  double t1, t2;
  int i;
  double complex *d_buf = NULL;

  QNc(op_CLOVER_D, _eoprec_MdagM_arg_s) *a = (QNc(op_CLOVER_D, _eoprec_MdagM_arg_s) *)op_arg;
  QLUA_ASSERT(2 * loc_dim == QNc(QOP_, _CLOVER_half_fermion_size)(a->clover_state));

  d_buf = malloc(sizeof(d_buf[0]) * loc_dim);
  QLUA_ASSERT(NULL != d_buf);

  for (i = 0 ; i < loc_dim ; i++)
    d_buf[i] = y[i];
  QNc(QOP_D, _CLOVER_half_fermion_from_blas)(a->y, (double *)d_buf, 2 * loc_dim);

  if (0 < a->poly_n) {
    QNc(QOP_D, _CLOVER_MxM_poly)(a->x, NULL, a->clover_gauge, a->y,
                                 a->poly_n, a->poly_a, a->poly_b, a->poly_c);
    QNc(QOP_, _CLOVER_performance)(&t1, &fl1, NULL, NULL, a->clover_state);

    QNc(QOP_D, _CLOVER_blas_from_half_fermion)((double *)d_buf, 2 * loc_dim, a->x);
    for (i = 0 ; i < loc_dim ; i++)
      x[i] = d_buf[i];
  } else {
    QNc(QOP_D, _CLOVER_M_operator)(a->x, a->clover_gauge, a->y);
    QNc(QOP_, _CLOVER_performance)(&t1, &fl1, NULL, NULL, a->clover_state);

    QNc(QOP_D, _CLOVER_M_operator_conjugated)(a->y, a->clover_gauge, a->x);
    QNc(QOP_, _CLOVER_performance)(&t2, &fl2, NULL, NULL, a->clover_state);
    t1  += t2;
    fl1 += fl2;

    QNc(QOP_D, _CLOVER_blas_from_half_fermion)((double *)d_buf, 2 * loc_dim, a->y);
    for (i = 0 ; i < loc_dim ; i++)
      x[i] = d_buf[i];
  }

  free(d_buf);

  if (t1 == 0)
    t1 = -1;
  if (QDP_this_node == qlua_master_node)
    printf("CLOVER MdagM(poly_n=%d): time = %.3f sec,"
           " perf = %.2f MFlops/sec\n",
           (0 < a->poly_n ? a->poly_n : 1),
           t1, fl1 * 1e-6 / t1);
}

/* CLOVER:eig_deflator_lanczos(nev, ncv, max_iter, tol, [param])
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
q_CL_make_deflator_lanczos(lua_State *L)
{
#define LANCZOS_MXM_DOUBLE  1
  /* by default, search for ev with smallest real part */
  const char *lanczos_which= "SR";
  const char *arpack_logfile = NULL;
  struct QNc(QOP_F, _CLOVER_HalfFermionMat) *hfm = NULL;
  /* operator parameters, init to empty */
  struct QNc(QOP_F, _CLOVER_Gauge) *gaugeF = NULL;
#ifdef LANCZOS_MXM_DOUBLE
  QNc(op_CLOVER_D, _eoprec_MdagM_arg_s) op_arg;
#else
  QNc(op_CLOVER_F, _eoprec_MdagM_arg_s) op_arg;
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

  mClover *c = qlua_checkClover(L, 1, NULL, 1);
  if (NULL == c->state || NULL == c->gauge)
    return luaL_error(L, "closed CLOVER used");

  mDeflatorState *d = NULL;

  /* parse parameters */
  int nev, ncv, max_iter;
  double tol;
  nev     = qlua_checkint(L, 2, "expect NEV at #2");
  ncv     = qlua_checkint(L, 3, "expect NCV at #3");
  max_iter= qlua_checkint(L, 4, "expect MAX_ITER at #4");
  tol     = luaL_checknumber(L, 5);

  /* parse optional parameters */
  int eigcg_vmax  = 0,
    eigcg_umax  = 0,
    eigcg_nev   = 0;
  double eigcg_eps= 0.;

  if (qlua_checkopt_paramtable(L, 6)) {
    if (qlua_tabpushopt_key(L, 6, "cheb_accel")) {
      /* Chebyshev acceleration parameters */
#ifdef LANCZOS_MXM_DOUBLE
      QNc(op_CLOVER_D, _eoprec_MdagM_arg_s) *a = &op_arg;
#else
      QNc(op_CLOVER_F, _eoprec_MdagM_arg_s) *a = &op_arg;
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
        QNc(QOP_, _CLOVER_poly_normalize)(cheb_n,
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


  /* create Qlua deflator object {CLOVER, DeflatorState} and set its META */
  /* FIXME should be refactored into a separate function?;
     this code duplicates parts of q_CL_make_deflator(qclover.c) */
  qlua_ObjLattice(L, 1);
  Sidx = lua_gettop(L);
  lua_createtable(L, 2, 0);
  lua_pushvalue(L, 1);
  lua_rawseti(L, -2, 1);
  if (NULL == (d = q_newDeflatorState(L, Sidx)))
    return luaL_error(L, "cannot create deflator state");
  lua_rawseti(L, -2, 2);
  qlua_createLatticeTable(L, Sidx, mtDeflator, qCloverDeflator,
                          CloverDeflatorName);
  lua_setmetatable(L, -2);



  CALL_QDP(L);

  gaugeF = NULL;
  if (QNc(QOP_, _CLOVER_gauge_float_from_double)(&gaugeF, c->gauge))
    return luaL_error(L, "QOP_CLOVER_gauge_float_from_double() failed");

#ifdef LANCZOS_MXM_DOUBLE
  op_arg.clover_state = c->state;
  op_arg.clover_gauge = c->gauge;

  op_arg.x = op_arg.y = NULL;
  if (QNc(QOP_D, _CLOVER_allocate_half_fermion)(&op_arg.x, c->state)
      || QNc(QOP_D, _CLOVER_allocate_half_fermion)(&op_arg.y, c->state))
    return luaL_error(L, "cannot allocate HalfFermion");
#else
  op_arg.clover_state   = c->state;
  if (QNc(QOP_, _CLOVER_gauge_float_from_double)(&(op_arg.clover_gauge), c->gauge))
    return luaL_error(L, "CLOVER_gauge_float_from_double() failed");

  op_arg.x = op_arg.y = NULL;
  if (QNc(QOP_F, _CLOVER_allocate_half_fermion)(&op_arg.x, c->state)
      || QNc(QOP_F, _CLOVER_allocate_half_fermion)(&op_arg.y, c->state))
    return luaL_error(L, "cannot allocate HalfFermion");
#endif/*LANCZOS_MXM_DOUBLE*/

  MPI_Comm mpi_comm = MPI_COMM_WORLD; /* FIXME any better choice? */

  loc_dim = QNc(QOP_, _CLOVER_half_fermion_size)(c->state) / 2;
  QLUA_ASSERT(0 == QNc(QOP_, _CLOVER_half_fermion_size)(c->state) % 2);

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
    if (QNc(QOP_F, _CLOVER_alloc_half_fermion_matrix)(&hfm, c->state, inplace_umax))
      return luaL_error(L, "CLOVER_alloc_half_fermion_matrix failed");
    if (QNc(QOP_F, _CLOVER_blas_view_half_fermion_matrix)(hfm, &hfm_nrow_loc, &hfm_ncol,
                                                          &hfm_blas_ptr, &hfm_ld))
      return luaL_error(L, QNc(QOP_, _CLOVER_error)(c->state));

    if (0 != (status = lanczos_inplace_float(
                                             L, mpi_comm,
#ifdef LANCZOS_MXM_DOUBLE
                                             QNc(op_CLOVER_F, _eoprec_MdagM_double_op),
#else
                                             QNc(op_CLOVER_F, _eoprec_MdagM_op),
#endif/*LANCZOS_MXM_DOUBLE*/
                                             &op_arg,
                                             lanczos_which, loc_dim, nev, ncv, max_iter, tol, NULL/*TODO implement v0*/,
                                             eval, (float complex *)hfm_blas_ptr, hfm_ld, hfm_ncol,
                                             &n_iters, &nconv, arpack_logfile))) {
      return luaL_error(L, "lanczos_float_inplace returned %d", status);
    }
    if (QNc(QOP_F, _CLOVER_create_deflator_inplace)(&(d->deflator), gaugeF, &hfm,
                                                    nconv, eigcg_vmax, eigcg_nev, eigcg_eps,
                                                    eigcg_umax))
      return luaL_error(L, QNc(QOP_, _CLOVER_error)(c->state));

  } else {
    if (0 != (status = lanczos_float(
                                     L, mpi_comm,
#ifdef LANCZOS_MXM_DOUBLE
                                     QNc(op_CLOVER_F, _eoprec_MdagM_double_op),
#else
                                     QNc(op_CLOVER_F, _eoprec_MdagM_op),
#endif/*LANCZOS_MXM_DOUBLE*/
                                     &op_arg,
                                     lanczos_which, loc_dim, nev, ncv, max_iter, tol, NULL/*TODO implement v0*/,
                                     &eval, &evec, &n_iters, &nconv, arpack_logfile)))
      return luaL_error(L, "lanczos_float returned %d", status);
    /* FIXME rewrite with clear explanation for the choice of
       default values */
    if (eigcg_umax <= 0 || eigcg_umax <= nconv)
      eigcg_umax = nconv;


    if (QNc(QOP_F, _CLOVER_create_deflator)(&d->deflator, c->state,
                                            eigcg_vmax, eigcg_nev, eigcg_eps, eigcg_umax))
      return luaL_error(L, "CLOVER_create_deflator() failed");

    /* fill deflator with e.vecs */
    if (QNc(QOP_F, _CLOVER_deflator_start_load)(d->deflator))
      return luaL_error(L, "CLOVER_deflator_start_load() failed");
    /* (have space for eigcg_umax) */
    n_evecs = (nconv <= eigcg_umax ? nconv : eigcg_umax);

    struct QNc(QOP_F, _CLOVER_HalfFermion) *hf_buf = NULL;
    if (QNc(QOP_F, _CLOVER_allocate_half_fermion)(&hf_buf, c->state))
      return luaL_error(L, "cannot allocate HalfFermion");

    for (i = 0 ; i < n_evecs ; i++) {
      QNc(QOP_F, _CLOVER_half_fermion_from_blas)(hf_buf,
                                                 (float *)(evec + i * loc_dim), 2 * loc_dim);
      if (QNc(QOP_F, _CLOVER_deflator_add_vector)(gaugeF, d->deflator, hf_buf)) {
        QNc(QOP_F, _CLOVER_free_half_fermion)(&hf_buf);
        return luaL_error(L, "CLOVER_deflator_add_vector() failed");
      }
    }

    if (QNc(QOP_F, _CLOVER_deflator_stop_load)(d->deflator))
      return luaL_error(L, "CLOVER_deflator_end_load() failed");
  }

  /* initialize CLOVER side of deflator */
  d->nev  = eigcg_nev;
  d->vmax = eigcg_vmax;
  d->umax = eigcg_umax;

  CALL_QDP(L);

  QNc(QOP_F,_CLOVER_deflator_eigcg_stop)(d->deflator);

  /* cleanup */
  if (NULL != evec) free(evec);
  if (NULL != eval) free(eval);
  if (NULL != gaugeF) QNc(QOP_F, _CLOVER_free_gauge)(&(gaugeF));
#ifdef LANCZOS_MXM_DOUBLE
  if (NULL != op_arg.y) QNc(QOP_D, _CLOVER_free_half_fermion)(&op_arg.y);
  if (NULL != op_arg.x) QNc(QOP_D, _CLOVER_free_half_fermion)(&op_arg.x);
#else
  if (NULL != op_arg.clover_gauge) QNc(QOP_F, _CLOVER_free_gauge)(&(op_arg.clover_gauge));
  if (NULL != op_arg.y) QNc(QOP_F, _CLOVER_free_half_fermion)(&op_arg.y);
  if (NULL != op_arg.x) QNc(QOP_F, _CLOVER_free_half_fermion)(&op_arg.x);
#endif/*LANCZOS_MXM_DOUBLE*/

  if (NULL != op_arg.poly_a) qlua_free(L, op_arg.poly_a);
  if (NULL != op_arg.poly_b) qlua_free(L, op_arg.poly_b);
  if (NULL != op_arg.poly_c) qlua_free(L, op_arg.poly_c);

  if (NULL != hfm)  QNc(QOP_F,_CLOVER_free_half_fermion_matrix)(&hfm);

  lua_pushnumber(L, nconv);
  lua_pushnumber(L, n_iters);
  return 3; /* deflator object, n_converged, n_iter */
}

#endif /* HAS_ARPACK */
/* XXXXXXX */

/***** clover interface */
static int
q_CL_fmt(lua_State *L)
{
    char fmt[72];
    mClover *c = qlua_checkClover(L, 1, NULL, 0);

    if (c->state)
        sprintf(fmt, "Clover[%p]", c);
    else
        sprintf(fmt, "Clover(closed)");

    lua_pushstring(L, fmt);

    return 1;
}

static int
q_CL_gc(lua_State *L)
{
    mClover *c = qlua_checkClover(L, 1, NULL, 0);

    if (c->gauge) {
      QNc(QOP_, _CLOVER_free_gauge)(&c->gauge);
        c->gauge = 0;
    }
    if (c->state) {
      QNc(QOP_, _CLOVER_fini)(&c->state);
        c->state = 0;
    }

    return 0;
}

static int
q_CL_close(lua_State *L)
{
    mClover *c = qlua_checkClover(L, 1, NULL, 1);

    if (c->gauge)
      QNc(QOP_, _CLOVER_free_gauge)(&c->gauge);
    if (c->state)
      QNc(QOP_, _CLOVER_fini)(&c->state);
    c->gauge = 0;
    c->state = 0;

    return 0;
}

static double
q_CL_D_reader(const int p[], int c, int d, int re_im, void *e)
{
    CL_D_env *env = e;
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
q_CL_D_writer(const int p[], int c, int d, int re_im, double v, void *e)
{
    CL_D_env *env = e;
    int i = QDP_index_L(env->lat, p);
    QNc(QLA_D, _DiracFermion) *f = env->f;

    if (re_im == 0) {
        QLA_real(QLA_elem_D(f[i], c, d)) = v;
    } else {
        QLA_imag(QLA_elem_D(f[i], c, d)) = v;
    }
}

static double
q_CL_P_reader(const int p[], int c, int d, int re_im, void *e)
{
    qCL_P_env *env = e;
    int i = QDP_index_L(env->lat, p);
    QLA_D_Real xx;

    if (re_im == 0) {
        QLA_r_eq_Re_c(xx, QLA_elem_P(env->in[i], c, d, env->c, env->d));
    } else {
        QLA_r_eq_Im_c(xx, QLA_elem_P(env->in[i], c, d, env->c, env->d));
    }

    return xx;
}

static void
q_CL_P_writer(const int p[], int c, int d, int re_im, double v, void *e)
{
    qCL_P_env *env = e;
    int i = QDP_index_L(env->lat, p);

    if (re_im == 0) {
        QLA_real(QLA_elem_P(env->out[i], c, d, env->c, env->d)) = v;
    } else {
        QLA_imag(QLA_elem_P(env->out[i], c, d, env->c, env->d)) = v;
    }
}

static int
q_CL_operator(lua_State *L,
              const char *name,
              int (*op)(struct QNc(QOP_D, _CLOVER_Fermion) *result,
                        const struct QNc(QOP_D, _CLOVER_Gauge) *gauge,
                        const struct QNc(QOP_D, _CLOVER_Fermion) *source))
{
    mClover *c = qlua_checkClover(L, 1, NULL, 1);
    mLattice *S = qlua_ObjLattice(L, 1);
    int Sidx = lua_gettop(L);

    switch (qlua_qtype(L, 2)) {
    case QNz(qLatDirFerm): {
        QNz(mLatDirFerm) *psi = QNz(qlua_checkLatDirFerm)(L, 2, S, QLUA_CLOVER_NC);
        QNz(mLatDirFerm) *eta = QNz(qlua_newLatDirFerm)(L, Sidx, QLUA_CLOVER_NC);
        struct QNc(QOP_, _CLOVER_Fermion) *c_psi;
        struct QNc(QOP_, _CLOVER_Fermion) *c_eta;
        QNc(QLA_D, _DiracFermion) *e_psi;
        QNc(QLA_D, _DiracFermion) *e_eta;
        CL_D_env env;

        CALL_QDP(L);
        e_psi = QNc(QDP_D, _expose_D)(psi->ptr);
        env.lat = S->lat;
        env.f = e_psi;
        if (QNc(QOP_, _CLOVER_import_fermion)(&c_psi, c->state, q_CL_D_reader, &env))
            return luaL_error(L, "CLOVER_import_fermion() failed");
        QNc(QDP_D, _reset_D)(psi->ptr);

        if (QNc(QOP_, _CLOVER_allocate_fermion)(&c_eta, c->state))
            return luaL_error(L, "CLOVER_allocate_fermion() failed");

        (*op)(c_eta, c->gauge, c_psi);

        e_eta = QNc(QDP_D, _expose_D)(eta->ptr);
        env.f = e_eta;
        QNc(QOP_, _CLOVER_export_fermion)(q_CL_D_writer, &env, c_eta);
        QNc(QDP_D, _reset_D)(eta->ptr);

        QNc(QOP_, _CLOVER_free_fermion)(&c_eta);
        QNc(QOP_, _CLOVER_free_fermion)(&c_psi);

        return 1;
    }
    case QNz(qLatDirProp): {
        QNz(mLatDirProp) *psi = QNz(qlua_checkLatDirProp)(L, 2, S, QLUA_CLOVER_NC);
        QNz(mLatDirProp) *eta = QNz(qlua_newLatDirProp)(L, Sidx, QLUA_CLOVER_NC);
        struct QNc(QOP_, _CLOVER_Fermion) *c_psi;
        struct QNc(QOP_, _CLOVER_Fermion) *c_eta;
        qCL_P_env env;

        CALL_QDP(L);
        if (QNc(QOP_, _CLOVER_allocate_fermion)(&c_eta, c->state))
            return luaL_error(L, "CLOVER_allocate_fermion() failed");

        env.lat = S->lat;
        env.in = QNc(QDP_D, _expose_P)(psi->ptr);
        env.out = QNc(QDP_D, _expose_P)(eta->ptr);
        for (env.c = 0; env.c < QNc(QOP_, _CLOVER_COLORS); env.c++) {
            for (env.d = 0; env.d < QNc(QOP_,  _CLOVER_FERMION_DIM); env.d++) {
              if (QNc(QOP_, _CLOVER_import_fermion)(&c_psi, c->state,
                                              q_CL_P_reader, &env))
                    return luaL_error(L, "CLOVER_import_fermion failed");

                (*op)(c_eta, c->gauge, c_psi);

                QNc(QOP_, _CLOVER_free_fermion)(&c_psi);
                QNc(QOP_, _CLOVER_export_fermion)(q_CL_P_writer, &env, c_eta);
            }
        }
        QNc(QOP_, _CLOVER_free_fermion)(&c_eta);
        QNc(QDP_D, _reset_P)(psi->ptr);
        QNc(QDP_D, _reset_P)(eta->ptr);

        return 1;
    }
    default:
        break;
    }
    return luaL_error(L, "bad arguments in Clover:%s", name);
}

static int
q_CL_D(lua_State *L)
{
  return q_CL_operator(L, "D", QNc(QOP_D, _CLOVER_D_operator));
}

static int
q_CL_Dx(lua_State *L)
{
  return q_CL_operator(L, "Dx", QNc(QOP_D, _CLOVER_D_operator_conjugated));
}

static void
q_CL_P_writer_full(const int p[], int a, int i, int b, int j,
		   int re_im, double v, void *e)
{
    qCL_P_env *env = e;
    int x = QDP_index_L(env->lat, p);

    if (re_im == 0) {
        QLA_real(QLA_elem_P(env->out[x], a, i, b, j)) = v;
    } else {
        QLA_imag(QLA_elem_P(env->out[x], a, i, b, j)) = v;
    }
}

static int
q_CL_inv_clovterm(lua_State *L)
{
    mClover *c = qlua_checkClover(L, 1, NULL, 1);
    mLattice *S = qlua_ObjLattice(L, 1);
    int Sidx = lua_gettop(L);

    QNz(mLatDirProp) *eta = QNz(qlua_newLatDirProp)(L, Sidx, QLUA_CLOVER_NC);
    qCL_P_env env;

    CALL_QDP(L);

    env.lat = S->lat;
    env.out = QNc(QDP_D, _expose_P)(eta->ptr);

    QNc(QOP_D, _CLOVER_export_inv_clover)(q_CL_P_writer_full, &env, c->gauge);

    QNc(QDP_D, _reset_P)(eta->ptr);

    return 1;
}

/* the standard clover solver */
static int
q_CL_std_solver(lua_State *L,
                struct QNc(QOP_, _CLOVER_Fermion) *solution,
                int *out_iters,
                double *out_epsilon,
                const struct QNc(QOP_, _CLOVER_Fermion) *rhs,
                int log_level)
{
    mClover *c = qlua_checkClover(L, lua_upvalueindex(2), NULL, 1);
    double eps = luaL_checknumber(L, lua_upvalueindex(3));
    int max_iters = luaL_checkint(L, lua_upvalueindex(4));

    return QNc(QOP_, _CLOVER_D_CG)(solution, out_iters, out_epsilon,
                           rhs, c->gauge, rhs, max_iters, eps,
                           log_level);
}

static CloverSolver std_solver = { q_CL_std_solver, "CG" };

static int
q_CL_make_solver(lua_State *L)
{
    qlua_checkClover(L, 1, NULL, 1);   /* mClover */
    luaL_checknumber(L, 2);            /* double epsilon */
    (void)luaL_checkint(L, 3);         /* int max_iter */

    lua_pushlightuserdata(L, &std_solver);  /* cl[1]: solver */
    lua_pushvalue(L, 1);                    /* cl[2]: mClover */
    lua_pushvalue(L, 2);                    /* cl[3]: epsilon */
    lua_pushvalue(L, 3);                    /* cl[4]: max_iters */
    lua_pushcclosure(L, q_dirac_solver, 4);

    return 1;
}

/* the standard clover MxM solver */
static int
q_CL_mxm_std_solver(lua_State *L,
                struct QNc(QOP_, _CLOVER_HalfFermion) *solution,
                int *out_iters,
                double *out_epsilon,
                const struct QNc(QOP_, _CLOVER_HalfFermion) *rhs,
                int log_level)
{
    mClover *c = qlua_checkClover(L, lua_upvalueindex(2), NULL, 1);
    double eps = luaL_checknumber(L, lua_upvalueindex(3));
    int max_iters = luaL_checkint(L, lua_upvalueindex(4));

    return QNc(QOP_, _CLOVER_MxM_CG)(solution, out_iters, out_epsilon,
                           rhs, c->gauge, rhs, max_iters, eps,
                           log_level);
}

static CloverMxMSolver mxm_std_solver = { q_CL_mxm_std_solver, "MxM_CG" };

static int
q_CL_make_mxm_solver(lua_State *L)
{
    qlua_checkClover(L, 1, NULL, 1);   /* mClover */
    luaL_checknumber(L, 2);            /* double epsilon */
    (void)luaL_checkint(L, 3);         /* int max_iter */

    lua_pushlightuserdata(L, &mxm_std_solver);  /* cl[1]: solver */
    lua_pushvalue(L, 1);                    /* cl[2]: mClover */
    lua_pushvalue(L, 2);                    /* cl[3]: epsilon */
    lua_pushvalue(L, 3);                    /* cl[4]: max_iters */
    lua_pushcclosure(L, q_mxm_solver, 4);

    return 1;
}

/* the mixed clover solver */
static int
q_CL_mixed_solver(lua_State *L,
                  struct QNc(QOP_, _CLOVER_Fermion) *solution,
                  int *out_iters,
                  double *out_epsilon,
                  const struct QNc(QOP_, _CLOVER_Fermion) *rhs,
                  int log_level)
{
    mClover *c = qlua_checkClover(L, lua_upvalueindex(2), NULL, 1);
    double f_eps    = luaL_checknumber(L, lua_upvalueindex(3));
    int inner_iters = luaL_checkint(L, lua_upvalueindex(4));
    double eps      = luaL_checknumber(L, lua_upvalueindex(5));
    int max_iters   = luaL_checkint(L, lua_upvalueindex(6));

    return QNc(QOP_, _CLOVER_mixed_D_CG)(solution, out_iters, out_epsilon,
                                 rhs, c->gauge, rhs,
                                 inner_iters, f_eps,
                                 max_iters, eps,
                                 log_level);
}

static CloverSolver mixed_solver = { q_CL_mixed_solver, "mixedCG" };

static int
q_CL_make_mixed_solver(lua_State *L)
{
    qlua_checkClover(L, 1, NULL, 1);   /* mClover */
    luaL_checknumber(L, 2);            /* double inner_epsilon */
    (void)luaL_checkint(L, 3);         /* int inner_iter */
    luaL_checknumber(L, 4);            /* double epsilon */
    (void)luaL_checkint(L, 5);         /* int max_iter */

    lua_pushlightuserdata(L, &mixed_solver);  /* cl[1]: solver */
    lua_pushvalue(L, 1);                      /* cl[2]: mClover */
    lua_pushvalue(L, 2);                      /* cl[3]: inner_epsilon */
    lua_pushvalue(L, 3);                      /* cl[4]: inner_iters */
    lua_pushvalue(L, 4);                      /* cl[5]: epsilon */
    lua_pushvalue(L, 5);                      /* cl[6]: max_iters */
    lua_pushcclosure(L, q_dirac_solver, 6);

    return 1;
}

#define Nu  QNc(QOP_, _CLOVER_DIM)
#define Nf  ((QNc(QOP_, _CLOVER_DIM) * (QNc(QOP_, _CLOVER_DIM) - 1)) / 2)
#define Nt  (Nu + Nf)
#define Nz  (Nu + Nf + 6)

typedef struct {
  QDP_Lattice *lat;
  int lattice[QNc(QOP_, _CLOVER_DIM)];
  QLA_D_Complex bf[QNc(QOP_, _CLOVER_DIM)];
  QNc(QLA_D, _ColorMatrix) *uf[Nu + Nf];
} QCArgs;

static double
q_CL_u_reader(int d, const int p[], int a, int b, int re_im, void *env)
{
    QLA_D_Complex z;
    QCArgs *args = env;
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

static double
q_CL_f_reader(int mu, int nu, const int p[], int a, int b, int re_im, void *env)
{
    QLA_D_Real xx;
    QCArgs *args = env;
    int i = QDP_index_L(args->lat, p);
    int d, xm, xn;

    for (d = 0, xm = 0; xm < QNc(QOP_, _CLOVER_DIM); xm++) {
        for (xn = xm + 1; xn < QNc(QOP_, _CLOVER_DIM); xn++, d++) {
            if ((xn == nu) && (xm == mu))
                goto found;
        }
    }
    return 0.0; /* should never happen */

found:
    /* NB:m stores F * 8i */
    if (re_im == 0) {
        QLA_r_eq_Im_c(xx, QLA_elem_M(args->uf[Nu + d][i], a, b));
        xx = xx / 8;
    } else {
        QLA_r_eq_Re_c(xx, QLA_elem_M(args->uf[Nu + d][i], a, b));
        xx = -xx / 8;
    }

    return xx;
}

typedef struct {
  QDP_Lattice *lat;
  int lattice[QNc(QOP_, _CLOVER_DIM)];
  QLA_D_Complex kappa[QNc(QOP_, _CLOVER_DIM)];
  QLA_D_Complex c_sw[QNc(QOP_, _CLOVER_DIM)][QNc(QOP_, _CLOVER_DIM)];
  QLA_D_Complex boundary[QNc(QOP_, _CLOVER_DIM)];
  QNc(QLA_D, _ColorMatrix) *uf[Nu + Nf];
} QAniCloArgs;

static double
q_AniC_u_reader(int d, const int p[], int a, int b, int re_im, void *env)
{
  QLA_D_Complex z, w;
  QAniCloArgs *args = env;
  int i = QDP_index_L(args->lat, p);
  
  if (p[d] == (args->lattice[d] - 1)) {
    QLA_c_eq_c_times_c(z, args->boundary[d], QLA_elem_M(args->uf[d][i], a, b));
  } else {
    QLA_c_eq_c(z, QLA_elem_M(args->uf[d][i], a, b));
  }
  QLA_c_eq_c_times_c(w, z, args->kappa[d]);
  if (re_im == 0)
    return QLA_real(w);
  else
    return QLA_imag(w);
}

static double
q_AniC_f_reader(int mu, int nu, const int p[], int a, int b, int re_im, void *env)
{
  QLA_D_Real x;
  QLA_D_Complex z, w;
  QAniCloArgs *args = env;
  int i = QDP_index_L(args->lat, p);
  int d, xm, xn;

  for (d = 0, xm = 0; xm < QNc(QOP_, _CLOVER_DIM); xm++) {
    for (xn = xm + 1; xn < QNc(QOP_, _CLOVER_DIM); xn++, d++) {
      if ((xn == nu) && (xm == mu))
	goto found;
    }
  }
  return 0.0; /* should never happen */
  
 found:
  QLA_c_eq_c(z, QLA_elem_M(args->uf[Nu + d][i], a, b));
  QLA_c_eq_c_times_c(w, z, args->c_sw[mu][nu]);
  if (re_im == 0) {
    QLA_r_eq_Im_c(x, w);
    x = x / 8;
  } else {
    QLA_r_eq_Re_c(x, w);
    x = - x / 8;
  }
  return x;
}


static int
build_clover(lua_State *L,
	     mLattice *S,
	     double kappa,
	     double c_sw,
	     mClover *clover,
	     QNc(QLA_D, _ColorMatrix) *uf[],
	     double (*u_reader)(int, const int [], int, int, int, void *),
	     double (*f_reader)(int, int, const int [], int, int, int, void *),
	     void *args)
{
  QNc(QDP_D, _ColorMatrix) *UF[Nz];

  luaL_checktype(L, 1, LUA_TTABLE);
  CALL_QDP(L);
  
  /* create a temporary F, and temp M */
  int i;
  for (i = QNc(QOP_, _CLOVER_DIM); i < Nz; i++)
    UF[i] = QNc(QDP_D, _create_M_L)(S->lat);
  
  /* extract U from the arguments */
  for (i = 0; i < QNc(QOP_, _CLOVER_DIM); i++) {
    lua_pushnumber(L, i + 1); /* [sic] lua indexing */
    lua_gettable(L, 1);
    UF[i] = QNz(qlua_checkLatColMat)(L, -1, S, QLUA_CLOVER_NC)->ptr;
    lua_pop(L, 1);
  }

  int mu, nu;
  QDP_Shift *neighbor = QDP_neighbor_L(S->lat);
  CALL_QDP(L); /* just in case, because we touched LUA state above */
  /* compute 8i*F[mu,nu] in UF[Nf...] */
  for (i = 0, mu = 0; mu < QNc(QOP_, _CLOVER_DIM); mu++) {
    for (nu = mu + 1; nu < QNc(QOP_, _CLOVER_DIM); nu++, i++) {
      /* clover in [mu, nu] --> UF[Nu + i] */
      QNc(QDP_D, _M_eq_sM)(UF[Nt], UF[nu], neighbor[mu], QDP_forward, S->all);
      QNc(QDP_D, _M_eq_Ma_times_M)(UF[Nt+1], UF[nu], UF[mu], S->all);
      QNc(QDP_D, _M_eq_M_times_M)(UF[Nt+2], UF[Nt+1], UF[Nt], S->all);
      QNc(QDP_D, _M_eq_sM)(UF[Nt+3], UF[Nt+2], neighbor[nu], QDP_backward, S->all);
      QNc(QDP_D, _M_eq_M_times_Ma)(UF[Nt+4], UF[Nt+3], UF[mu], S->all);
      QNc(QDP_D, _M_eq_sM)(UF[Nt+1], UF[mu], neighbor[nu], QDP_forward, S->all);
      QNc(QDP_D, _M_eq_Ma_times_M)(UF[Nt+5], UF[mu], UF[Nt+3], S->all);
      QNc(QDP_D, _M_eq_M_times_Ma)(UF[Nt+2], UF[Nt], UF[Nt+1], S->all);
      QNc(QDP_D, _M_eq_M_times_Ma)(UF[Nt+3], UF[Nt+2], UF[nu], S->all);
      QNc(QDP_D, _M_peq_M_times_M)(UF[Nt+4], UF[mu], UF[Nt+3], S->all);
      QNc(QDP_D, _M_peq_M_times_M)(UF[Nt+5], UF[Nt+3], UF[mu], S->all);
      QNc(QDP_D, _M_eq_sM)(UF[Nt+2], UF[Nt+5], neighbor[mu], QDP_backward, S->all);
      QNc(QDP_D, _M_peq_M)(UF[Nt+4], UF[Nt+2], S->all);
      QNc(QDP_D, _M_eq_M)(UF[Nu+i], UF[Nt+4], S->all);
      QNc(QDP_D, _M_meq_Ma)(UF[Nu+i], UF[Nt+4], S->all);
    }
  }

  struct QNc(QOP_, _CLOVER_Config) cc;
  cc.self = S->node;
  cc.master_p = QMP_is_primary_node();
  cc.rank = S->rank;
  cc.lat = S->dim;
  cc.net = S->net;
  cc.neighbor_up = S->neighbor_up;
  cc.neighbor_down = S->neighbor_down;
  cc.sublattice = qlua_sublattice;
  cc.env = S;
  if (QNc(QOP_, _CLOVER_init)(&clover->state, &cc))
    return luaL_error(L, "CLOVER_init() failed");

  /* import the gauge field */
  for (i = 0; i < Nt; i++) {
    uf[i] = QNc(QDP_D, _expose_M)(UF[i]);
  }

  if (QNc(QOP_, _CLOVER_import_gauge)(&clover->gauge, clover->state, kappa, c_sw,
                              u_reader, f_reader, args)) {
      return luaL_error(L, "CLOVER_import_gauge() failed");
  }

  for (i = 0; i < Nt; i++)
    QNc(QDP_D, _reset_M)(UF[i]);

  /* clean up temporaries */
  for (i = QNc(QOP_, _CLOVER_DIM); i < Nz; i++)
    QNc(QDP_D, _destroy_M)(UF[i]);

  return 1;
}

static void
start_clover(lua_State *L, mLattice **ptr_S, mClover **ptr_clover)
{
  luaL_checktype(L, 1, LUA_TTABLE);
  lua_pushnumber(L, 1);
  lua_gettable(L, 1);
  QNz(qlua_checkLatColMat)(L, -1, NULL, QLUA_CLOVER_NC);
  *ptr_S = qlua_ObjLattice(L, -1);
  int Sidx = lua_gettop(L);
  *ptr_clover = qlua_newClover(L, Sidx);
  if ((*ptr_S)->rank != QNc(QOP_, _CLOVER_DIM))
    luaL_error(L, "clover is not implemented for #L=%d", (*ptr_S)->rank);
  if (QDP_Ns != QNc(QOP_,  _CLOVER_FERMION_DIM))
    luaL_error(L, "clover does not support Ns=%d", QDP_Ns);
}

/*
 *  qcd.Clover(U,         -- 1, {U0,U1,U2,U3}, a table of color matrices
 *             kappa,     -- 2, double, the hopping parameter
 *             c_sw,      -- 3, double, the clover term
 *             boundary)  -- 4, {r/c, ...}, a table of boundary phases
 */
static int
q_clover(lua_State *L)
{
    mLattice *S = NULL;
    mClover *clover = NULL;
    QCArgs args;
    double kappa = luaL_checknumber(L, 2);
    double c_sw = luaL_checknumber(L, 3);

    qlua_get_complex_vector(L, 4, QNc(QOP_, _CLOVER_DIM), args.bf, "boundary");
    start_clover(L, &S, &clover);
    args.lat = S->lat;
    QDP_latsize_L(S->lat, args.lattice);

    return build_clover(L, S, kappa, c_sw, clover, args.uf, q_CL_u_reader, q_CL_f_reader, &args);
}

/*
 * qcd.AnisotropicClover(U,                                -- gauge field, { Ux,Uy,Uz,Ut }
 *                       {kappa = { kx, ky, kz, kt },      -- complex kappa's in each direction
 *                        c_sw  = { {cxx,cxy,cxy,cxt },    -- complex clover factors
 *                                  {cyx,cyy,cyz,cyt },
 *                                  {czx,czy,czz,czt },
 *                                  {ctx,cty,ctz,ctt } },
 *                        boundary = { bx, by, bz, bt } })  -- complex boundary conditions
 */
static int
q_anisotropic_clover(lua_State *L)
{
  mLattice     *S = NULL;
  mClover      *clover = NULL;
  QAniCloArgs   args;
  int i, j;

  if (qlua_tabkey_tableopt(L, 2, "kappa") == 0)
    luaL_error(L, "missing kappa values");
  qlua_get_complex_vector(L, lua_gettop(L),  QNc(QOP_, _CLOVER_DIM), args.kappa, "kappa");
  lua_pop(L, 1);
  if (qlua_tabkey_tableopt(L, 2, "boundary") == 0)
    luaL_error(L, "missing boundary values");
  qlua_get_complex_vector(L, lua_gettop(L),  QNc(QOP_, _CLOVER_DIM), args.boundary, "boundary");
  lua_pop(L, 1);
  if (qlua_tabkey_tableopt(L, 2, "c_sw") == 0)
    luaL_error(L, "missing c_sw values");
  for (i = 0; i < QNc(QOP_, _CLOVER_DIM); i++) {
    if (qlua_tabidx_tableopt(L, lua_gettop(L), i + 1) == 0)
      luaL_error(L, "missing c_sw values for i = %d", i);
    int idx_i = lua_gettop(L);
    for (j = 0; j <  QNc(QOP_, _CLOVER_DIM); j++) {
      double rv, iv;
      qlua_tabidx_complex(L, idx_i, j + 1, &rv, &iv, "c_sw element");
      QLA_c_eq_r_plus_ir(args.c_sw[i][j], rv, iv);
    }
    lua_pop(L, 1);
  }
  lua_pop(L, 1);
  start_clover(L, &S, &clover);
  args.lat = S->lat;
  QDP_latsize_L(S->lat, args.lattice);

  return build_clover(L, S, 1.0, 1.0, clover, args.uf, q_AniC_u_reader, q_AniC_f_reader, &args);
}

static struct luaL_Reg mtClover[] = {
    { "__tostring",                 q_CL_fmt },
    { "__gc",                       q_CL_gc },
    { "close",                      q_CL_close },
    { "D",                          q_CL_D },
    { "Dx",                         q_CL_Dx },
    { "solver",                     q_CL_make_solver },
    { "mxm_solver",                 q_CL_make_mxm_solver },
    { "mixed_solver",               q_CL_make_mixed_solver },
    { "eig_deflator",               q_CL_make_deflator },
    { "inv_clovterm",               q_CL_inv_clovterm },
#ifdef HAS_ARPACK
    { "eig_deflator_lanczos",       q_CL_make_deflator_lanczos      },
#endif /* HAS_ARPACK */
    { NULL, NULL }
};

static mClover *
qlua_newClover(lua_State *L, int Sidx)
{
    mClover *c = lua_newuserdata(L, sizeof (mClover));

    c->state = 0;
    c->gauge = 0;
    qlua_createLatticeTable(L, Sidx, mtClover, qClover, CloverName);
    lua_setmetatable(L, -2);

    return c;
}

static mClover *
qlua_checkClover(lua_State *L, int idx, mLattice *S, int live)
{
    mClover *c = qlua_checkLatticeType(L, idx, qClover, CloverName);

    if (S) {
        mLattice *S1 = qlua_ObjLattice(L, idx);
        if (S1->id != S->id)
            luaL_error(L, "%s on a wrong lattice", CloverName);
        lua_pop(L, 1);
    }

    if (live && (c->state == 0 || c->gauge == 0))
        luaL_error(L, "using closed qcd.Clover");

    return c;
}

static struct luaL_Reg fClover[] = {
    { "Clover",                  q_clover },
    { "AnisotropicClover",       q_anisotropic_clover },
    { NULL,                      NULL }
};

int
init_clover(lua_State *L)
{
    luaL_register(L, qcdlib, fClover);
    return 0;
}

void
fini_clover(void)
{
}
