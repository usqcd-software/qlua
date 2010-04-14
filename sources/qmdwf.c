#include "qlua.h"                                                    /* DEPS */
#include "qmdwf.h"                                                   /* DEPS */
#include "qcomplex.h"                                                /* DEPS */
#include "qvector.h"                                                 /* DEPS */
#include "lattice.h"                                                 /* DEPS */
#include "qlayout.h"                                                 /* DEPS */
#include "latcolmat.h"                                               /* DEPS */
#include "latdirferm.h"                                              /* DEPS */
#include "latdirprop.h"                                              /* DEPS */
#define QOP_MDWF_DEFAULT_PRECISION QDP_Precision
#include "qop-mdwf3.h"
#include "qmp.h"

#include <math.h>

static const char mdwf_name[] = "MDWF";

static const char MDWFName[]         = "lattice.MDWF";
static const char MDWFDeflatorName[]       = "lattice.MDWF.Deflator";
static const char MDWFDeflatorStateName[]  = "lattice.MDWF.DeflatorState";

typedef int MDWFInverter(lua_State                      *L,
                         struct QOP_MDWF_Fermion        *solution,
                         int                            *out_iters,
                         double                         *out_epsilon,
                         const struct QOP_MDWF_Fermion  *rhs,
                         int                             log_level);

typedef struct {
    MDWFInverter    *proc;
    const char      *name;
} MDWFSolver;

typedef struct {
    struct QOP_MDWF_State        *state;
    struct QOP_MDWF_Parameters   *params;
    struct QOP_MDWF_Gauge        *gauge;
    int                           Ls;
    const char                   *name;
} mMDWF;

typedef struct {
    int nev;
    int umax;
    int vmax;
    struct QOP_MDWF_Deflator *deflator;
} mDeflatorState;

static mMDWF *qlua_newMDWF(lua_State *L, int Sidx);
static mMDWF *qlua_checkMDWF(lua_State *L,
                             int idx,
                             mLattice *S,
                             int live);

/* The generic solver */
typedef struct {
    QDP_Lattice *lat;
    int c, d;
    QLA_D3_DiracPropagator *in;
    QLA_D3_DiracPropagator *out;
    double s;
} DW_P_env;

typedef struct {
    QDP_Lattice *lat;
    QLA_D3_DiracFermion *f;
    double s;
} DW_D_env;

typedef struct {
    QDP_Lattice *lat;
    QLA_D3_DiracFermion **f;
    double s;
} DW_5_env;

static double
q_DW_5_reader_scaled(const int p[], int c, int d, int re_im, void *e)
{
    DW_5_env *env = e;
    int i = QDP_index_L(env->lat, p);
    QLA_D3_DiracFermion *f = env->f[p[QOP_MDWF_DIM]];
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
q_DW_5_writer_scaled(const int p[], int c, int d, int re_im, double v, void *e)
{
    DW_5_env *env = e;
    int i = QDP_index_L(env->lat, p);
    QLA_D3_DiracFermion *f = env->f[p[QOP_MDWF_DIM]];
    double s = env->s;
 
    v = v * s;
    if (re_im == 0) {
        QLA_real(QLA_elem_D(f[i], c, d)) = v;
    } else {
        QLA_imag(QLA_elem_D(f[i], c, d)) = v;
    }
}

#if 0 /* XXX scaled reader */
static double
q_CL_D_reader_scaled(const int p[], int c, int d, int re_im, void *e)
{
    CL_D_env *env = e;
    int i = QDP_index_L(env->lat, p);
    QLA_D3_DiracFermion *f = env->f;
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
    QLA_D3_DiracFermion *f = env->f;
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
#endif /* XXX scaled reader */

static int
q_dirac_solver(lua_State *L)
{
    MDWFSolver *solver = lua_touserdata(L, lua_upvalueindex(1));
    mMDWF *c = qlua_checkMDWF(L, lua_upvalueindex(2), NULL, 1);
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
        log_level = (QOP_MDWF_LOG_CG_RESIDUAL |
                     QOP_MDWF_LOG_EIG_POSTAMBLE |
                     QOP_MDWF_LOG_EIG_UPDATE1);
    else
        log_level = 0;

    switch (qlua_qtype(L, 1)) {
    case qTable: {
        mLatDirFerm3 *psi[c->Ls];
        mLatDirFerm3 *eta[c->Ls];
        struct QOP_MDWF_Fermion *c_psi;
        struct QOP_MDWF_Fermion *c_eta;
        QLA_D3_DiracFermion *e_psi[c->Ls];
        QLA_D3_DiracFermion *e_eta[c->Ls];
        DW_5_env env;
        double norm5;
        int status;
        int i;

        CALL_QDP(L);
        norm5 = 0;
        for (i = 0; i < c->Ls; i++) {
            QLA_D_Real normi;

            lua_pushnumber(L, i + 1); /* [sic] lua indexing */
            lua_gettable(L, 1);
            psi[i] = qlua_checkLatDirFerm3(L, -1, S, 3);
            QDP_D3_r_eq_norm2_D(&normi, psi[i]->ptr, S->all);
            norm5 += normi;
            lua_pop(L, 1);
        }
        if (norm5 == 0) {
            lua_createtable(L, c->Ls, 0);
            for (i = 0; i < c->Ls; i++) {
                qlua_newLatDirFerm3(L, Sidx, 3);
                lua_rawseti(L, -2, i + 1); /* [sic] Lua indexing */
            }
            lua_pushnumber(L, 0.0);
            lua_pushnumber(L, 0);
            return 3;
        }

        for (i = 0; i < c->Ls; i++) {
            e_psi[i] = QDP_D3_expose_D(psi[i]->ptr);
        }
        norm5 = sqrt(norm5);
        env.lat = S->lat;
        env.f = e_psi;
        env.s = 1 / norm5;
        if (QOP_MDWF_import_fermion(&c_psi, c->state, q_DW_5_reader_scaled,
                                    &env))
            return luaL_error(L, "MDWF_import_fermion() failed");
        for (i = 0; i < c->Ls; i++) {
            QDP_D3_reset_D(psi[i]->ptr);
        }
        if (QOP_MDWF_allocate_fermion(&c_eta, c->state))
            return luaL_error(L, "MDWF_allocate_fermion() failed");

        status = solver->proc(L, c_eta, &out_iters, &out_eps, c_psi, log_level);

        QOP_MDWF_performance(&t1, &fl1, NULL, NULL, c->state);
        if (t1 == 0)
            t1 = -1;
        if (qlua_primary_node)
            printf("CLOVER %s solver: status = %d,"
                   " eps = %.4e, iters = %d, time = %.3f sec,"
                   " perf = %.2f MFlops/sec\n",
                   solver->name, status,
                   out_eps, out_iters, t1, fl1 * 1e-6 / t1);
        QOP_MDWF_free_fermion(&c_psi);

        lua_createtable(L, c->Ls, 0);
        for (i = 0; i < c->Ls; i++) {
            eta[i] = qlua_newLatDirFerm3(L, Sidx, 3);
            e_eta[i] = QDP_D3_expose_D(eta[i]->ptr);
            lua_rawseti(L, -2, i + 1); /* [sic] Lua indexing */
        }

        env.lat = S->lat;
        env.f = e_eta;
        env.s = norm5;
        QOP_MDWF_export_fermion(q_DW_5_writer_scaled, &env, c_eta);
        for (i = 0; i < c->Ls; i++) {
            QDP_D3_reset_D(eta[i]->ptr);
        }
        
        lua_pushnumber(L, out_eps * norm5);
        lua_pushnumber(L, out_iters);
        return 3;
    }
#if 0 /* XXX solver */
    case qLatDirFerm3: {
        mLatDirFerm3 *psi = qlua_checkLatDirFerm3(L, 1, S, 3);
        mLatDirFerm3 *eta = qlua_newLatDirFerm3(L, Sidx, 3);
        struct QOP_CLOVER_Fermion *c_psi;
        struct QOP_CLOVER_Fermion *c_eta;
        CL_D_env env;
        QLA_D_Real rhs_norm2 = 0;
        double rhs_n;
        int status;

        CALL_QDP(L);
        QDP_D3_r_eq_norm2_D(&rhs_norm2, psi->ptr, S->all);
        if (rhs_norm2 == 0) {
            QDP_D3_D_eq_zero(eta->ptr, S->all);
            lua_pushnumber(L, 0.0);
            lua_pushnumber(L, 0);
            lua_pushnumber(L, 0);
            lua_pushnumber(L, 0);
            return 5;
        }
        rhs_n = sqrt(rhs_norm2);
        env.lat = S->lat;
        env.f = QDP_D3_expose_D(psi->ptr);
        env.s = 1 / rhs_n;
        if (QOP_CLOVER_import_fermion(&c_psi, c->state, q_CL_D_reader_scaled,
                                      &env))
            return luaL_error(L, "CLOVER_import_fermion() failed");
        QDP_D3_reset_D(psi->ptr);

        if (QOP_CLOVER_allocate_fermion(&c_eta, c->state))
            return luaL_error(L, "CLOVER_allocate_fermion() failed");

        status = solver->proc(L, c_eta, &out_iters, &out_eps, c_psi, log_level);

        QOP_CLOVER_performance(&t1, &fl1, NULL, NULL, c->state);
        if (t1 == 0)
            t1 = -1;
        if (qlua_primary_node)
            printf("CLOVER %s solver: status = %d,"
                   " eps = %.4e, iters = %d, time = %.3f sec,"
                   " perf = %.2f MFlops/sec\n",
                   solver->name, status,
                   out_eps, out_iters, t1, fl1 * 1e-6 / t1);

        env.lat = S->lat;
        env.f = QDP_D3_expose_D(eta->ptr);
        env.s = rhs_n;
        QOP_CLOVER_export_fermion(q_CL_D_writer_scaled, &env, c_eta);
        QDP_D3_reset_D(eta->ptr);
        
        QOP_CLOVER_free_fermion(&c_eta);
        QOP_CLOVER_free_fermion(&c_psi);

        if (status) {
            if (relaxed_p)
                return 0;
            else
                return luaL_error(L, QOP_CLOVER_error(c->state));
        }

        /* eta is on the stack already */
        lua_pushnumber(L, out_eps);
        lua_pushnumber(L, out_iters);
        lua_pushnumber(L, t1);
        lua_pushnumber(L, (double)fl1);

        return 5;
    }

    case qLatDirProp3: {
        mLatDirProp3 *psi = qlua_checkLatDirProp3(L, 1, S, 3);
        mLatDirProp3 *eta = qlua_newLatDirProp3(L, Sidx, 3);
        struct QOP_CLOVER_Fermion *c_psi;
        struct QOP_CLOVER_Fermion *c_eta;
        int gstatus = 0;
        int status;
        qCL_P_env env;
        QLA_D_Real rhs_norm2 = 0;
        double rhs_n;

        lua_createtable(L, QOP_CLOVER_COLORS, 0);  /* eps */
        lua_createtable(L, QOP_CLOVER_COLORS, 0);  /* iters */
        CALL_QDP(L);
        if (QOP_CLOVER_allocate_fermion(&c_eta, c->state))
            return luaL_error(L, "CLOVER_allocate_fermion() failed");

        QDP_D3_r_eq_norm2_P(&rhs_norm2, psi->ptr, S->all);
        if (rhs_norm2 == 0) {
            QDP_D3_P_eq_zero(eta->ptr, S->all);
            return 3;
        }
        rhs_n = sqrt(rhs_norm2);
        env.lat = S->lat;
        env.in = QDP_D3_expose_P(psi->ptr);
        env.out = QDP_D3_expose_P(eta->ptr);
        for (env.c = 0; env.c < QOP_CLOVER_COLORS; env.c++) {
            lua_createtable(L, QOP_CLOVER_FERMION_DIM, 0); /* eps.c */
            lua_createtable(L, QOP_CLOVER_FERMION_DIM, 0); /* iters.c */

            for (env.d = 0; env.d < QOP_CLOVER_FERMION_DIM; env.d++) {
                env.s = 1 / rhs_n;
                if (QOP_CLOVER_import_fermion(&c_psi, c->state,
                                              q_CL_P_reader_scaled, &env))
                    return luaL_error(L, "CLOVER_import_fermion() failed");
                status = solver->proc(L, c_eta, &out_iters, &out_eps, c_psi,
                                      log_level);

                QOP_CLOVER_performance(&t1, &fl1, NULL, NULL, c->state);
                if (t1 == 0)
                    t1 = -1;
                if (qlua_primary_node)
                    printf("CLOVER %s solver: status = %d, c = %d, d = %d,"
                           " eps = %.4e, iters = %d, time = %.3f sec,"
                           " perf = %.2f MFlops/sec\n",
                           solver->name, status,
                           env.c, env.d, out_eps, out_iters, t1,
                           fl1 * 1e-6 / t1);
                QOP_CLOVER_free_fermion(&c_psi);
                if (status) {
                    if (relaxed_p)
                        gstatus = 1;
                    else
                        return luaL_error(L, QOP_CLOVER_error(c->state));
                }

                env.s = rhs_n;
                QOP_CLOVER_export_fermion(q_CL_P_writer_scaled, &env, c_eta);
                lua_pushnumber(L, out_eps);
                lua_rawseti(L, -3, env.d + 1);
                lua_pushnumber(L, out_iters);
                lua_rawseti(L, -2, env.d + 1);
            }
            lua_rawseti(L, -3, env.c + 1);
            lua_rawseti(L, -3, env.c + 1);
        }
        QDP_D3_reset_P(psi->ptr);
        QDP_D3_reset_P(eta->ptr);
        QOP_CLOVER_free_fermion(&c_eta);

        if (gstatus)
            return 0;
        else
            return 3;
    }
#endif /* XXX solver */
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
        QOP_MDWF_free_deflator(&d->deflator);
    d->deflator = 0;

    return 0;
}

static struct luaL_Reg mtMDWFDeflatorState[] = {
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
    qlua_createLatticeTable(L, Sidx, mtMDWFDeflatorState, qMDWFDeflatorState,
                            MDWFDeflatorStateName);
    lua_setmetatable(L, -2);

    return d;
}

static mDeflatorState*
q_checkDeflatorState(lua_State *L, int idx, mLattice *S, int live)
{
    mDeflatorState *d = qlua_checkLatticeType(L, idx, qMDWFDeflatorState,
                                              MDWFDeflatorStateName);

    if (S) {
        mLattice *S1 = qlua_ObjLattice(L, idx);
        if (S1->id != S->id)
            luaL_error(L, "%s on a wrong lattice", MDWFDeflatorStateName);
        lua_pop(L, 1);
    }

    if (live && (d->deflator == 0))
        luaL_error(L, "Using closed MDWF.DeflatorState");

    return d;
}

/***** delfator interface */
static mMDWF *
q_Deflator_get_MDWF(lua_State *L, int idx, mLattice *S, int live)
{
    mMDWF *c;
    qlua_checktable(L, idx, "");

    lua_rawgeti(L, idx, 1);
    c = qlua_checkMDWF(L, -1, S, live);

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

    QOP_MDWF_free_deflator(&d->deflator);

    return 0;
}

static int
q_DF_reset(lua_State *L)
{
    mDeflatorState *d = q_Deflator_get_State(L, 1, NULL, 1);

    QOP_MDWF_deflator_reset(d->deflator);

    return 0;
}

static int
q_DF_stop(lua_State *L)
{
    mDeflatorState *d = q_Deflator_get_State(L, 1, NULL, 1);

    QOP_MDWF_deflator_stop(d->deflator);

    return 0;
}

static int
q_DF_resume(lua_State *L)
{
    mDeflatorState *d = q_Deflator_get_State(L, 1, NULL, 1);

    QOP_MDWF_deflator_resume(d->deflator);

    return 0;
}

static int
q_DF_eigenvalues(lua_State *L)
{
    mDeflatorState *d = q_Deflator_get_State(L, 1, NULL, 1);
    mVecReal *v = qlua_newVecReal(L, d->nev);
    double t[d->nev];
    int status = QOP_MDWF_deflator_eigen(t, d->deflator);

    if (status == 0) {
        int i;
        for (i = 0; i < d->nev; i++)
            v->val[i] = t[i];
    }
    if (status == 0)
        return 1;
    else
        return 0;
}

static int
q_DF_deflated_mixed_solver(lua_State *L,
                           struct QOP_MDWF_Fermion *solution,
                           int *out_iters,
                           double *out_epsilon,
                           const struct QOP_MDWF_Fermion *rhs,
                           int log_level)
{
    mMDWF          *c = qlua_checkMDWF(L, lua_upvalueindex(2), NULL, 1);
    mLattice       *S = qlua_ObjLattice(L, lua_upvalueindex(2));
    mDeflatorState *d = q_checkDeflatorState(L, lua_upvalueindex(3), S, 1);
    double      f_eps = luaL_checknumber(L, lua_upvalueindex(4));
    int   inner_iters = luaL_checkint(L, lua_upvalueindex(5));
    double        eps = luaL_checknumber(L, lua_upvalueindex(6));
    int     max_iters = luaL_checkint(L, lua_upvalueindex(7));

    lua_pop(L, 1);
    return QOP_MDWF_deflated_mixed_D_CG(solution, out_iters, out_epsilon,
                                        c->params, rhs, c->gauge, rhs,
                                        d->deflator,
                                        inner_iters, f_eps,
                                        max_iters, eps,
                                        log_level);
}

static MDWFSolver deflated_mixed_solver = {
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
    q_Deflator_get_MDWF(L, 1, NULL, 1);               /* clo[2]: MDWF */
    q_Deflator_get_State(L, 1, NULL, 1);              /* clo[3]: Deflator */
    lua_pushnumber(L, inner_eps);                     /* clo[4]: inner_eps */
    lua_pushnumber(L, inner_iter);                    /* clo[5]: inner_iter */
    lua_pushnumber(L, eps);                           /* clo[6]: epsilon */
    lua_pushnumber(L, max_iter);                      /* clo[7]: max_iter */
    lua_pushcclosure(L, q_dirac_solver, 7);

    return 1;
}

static struct luaL_Reg mtMDWFDeflator[] = {
    { "__tostring",         q_DF_fmt               },
    { "__newindex",         qlua_nowrite           },
    { "mixed_solver",       q_DF_make_mixed_solver },
    { "close",              q_DF_close             },
    { "reset",              q_DF_reset             },
    { "stop",               q_DF_stop              },
    { "resume",             q_DF_resume            },
    { "eigenvalues",        q_DF_eigenvalues       },
    { NULL,                 NULL                   }
};

static int
q_DW_make_deflator(lua_State *L)
{
    mMDWF *c = qlua_checkMDWF(L, 1, NULL, 1);
    qlua_ObjLattice(L, 1);
    int Sidx = lua_gettop(L);
    int vmax = luaL_checkint(L, 2);
    int nev = luaL_checkint(L, 3);
    double eps = luaL_checknumber(L, 4);
    int umax = luaL_checkint(L, 5);

    if ((nev <= 0) || (2 * nev >= vmax))
        return luaL_error(L, "bad eigenspace size");

    if ((c->state == 0) || (c->gauge == 0))
        return luaL_error(L, "closed MDWF used");

    lua_createtable(L, 2, 0);
    lua_pushvalue(L, 1);
    lua_rawseti(L, -2, 1);
    mDeflatorState *d = q_newDeflatorState(L, Sidx);
    d->nev = nev;
    d->vmax = vmax;
    d->umax = umax;

    CALL_QDP(L);
    if (QOP_MDWF_create_deflator(&d->deflator, c->state,
                                   vmax, nev, eps, umax))
        return luaL_error(L, "CLOVER_create_deflator() failed");

    lua_rawseti(L, -2, 2);
    qlua_createLatticeTable(L, Sidx, mtMDWFDeflator, qMDWFDeflator,
                            MDWFDeflatorName);
    lua_setmetatable(L, -2);

    return 1;
}

/***** clover interface */
static int
q_DW_fmt(lua_State *L)
{
    char fmt[72];
    mMDWF *c = qlua_checkMDWF(L, 1, NULL, 0);

    if (c->state)
        sprintf(fmt, "MDWF[%s]", c->name);
    else
        sprintf(fmt, "MDWF(closed)");

    lua_pushstring(L, fmt);

    return 1;
}

static int
q_DW_gc(lua_State *L)
{
    mMDWF *c = qlua_checkMDWF(L, 1, NULL, 0);

    if (c->gauge)
        QOP_MDWF_free_gauge(&c->gauge);
    c->gauge = 0;

    if (c->params)
        QOP_MDWF_free_parameters(&c->params);
    c->params = 0;
    
    if (c->state)
        QOP_MDWF_fini(&c->state);
    c->state = 0;
    
    return 0;
}

static int
q_DW_close(lua_State *L)
{
    mMDWF *c = qlua_checkMDWF(L, 1, NULL, 1);

    if (c->gauge)
        QOP_MDWF_free_gauge(&c->gauge);
    c->gauge = 0;

    if (c->params)
        QOP_MDWF_free_parameters(&c->params);
    c->params = 0;

    if (c->state)
        QOP_MDWF_fini(&c->state);
    c->state = 0;
    
    return 0;
}

static double
DW_5_reader(const int p[], int c, int d, int re_im, void *e)
{
    DW_5_env *env = e;
    int i = QDP_index_L(env->lat, p);
    QLA_D3_DiracFermion **f = env->f;
    QLA_D_Real xx;

    if (re_im == 0) {
        QLA_r_eq_Re_c(xx, QLA_elem_D(f[p[QOP_MDWF_DIM]][i], c, d));
    } else {
        QLA_r_eq_Im_c(xx, QLA_elem_D(f[p[QOP_MDWF_DIM]][i], c, d));
    }

    return xx;
}

static void
DW_5_writer(const int p[], int c, int d, int re_im, double v, void *e)
{
    DW_5_env *env = e;
    int i = QDP_index_L(env->lat, p);
    QLA_D3_DiracFermion **f = env->f;
 
    if (re_im == 0) {
        QLA_real(QLA_elem_D(f[p[QOP_MDWF_DIM]][i], c, d)) = v;
    } else {
        QLA_imag(QLA_elem_D(f[p[QOP_MDWF_DIM]][i], c, d)) = v;
    }
}

static int
q_DW_operator(lua_State *L,
              const char *name,
              int (*op)(struct QOP_D3_MDWF_Fermion *result,
                        const struct QOP_MDWF_Parameters *params,
                        const struct QOP_D3_MDWF_Gauge *gauge,
                        const struct QOP_D3_MDWF_Fermion *source))
{
    mMDWF *c = qlua_checkMDWF(L, 1, NULL, 1);
    mLattice *S = qlua_ObjLattice(L, 1);
    int Sidx = lua_gettop(L);

    switch (qlua_qtype(L, 2)) {
    case qTable: {
        mLatDirFerm3 *psi[c->Ls];
        mLatDirFerm3 *eta[c->Ls];
        struct QOP_MDWF_Fermion *c_psi;
        struct QOP_MDWF_Fermion *c_eta;
        QLA_D3_DiracFermion *e_psi[c->Ls];
        QLA_D3_DiracFermion *e_eta[c->Ls];
        DW_5_env env;
        int i;

        CALL_QDP(L);
        for (i = 0; i < c->Ls; i++) {
            lua_pushnumber(L, i + 1); /* [sic] lua indexing */
            lua_gettable(L, 2);
            psi[i] = qlua_checkLatDirFerm3(L, -1, S, 3);
            e_psi[i] = QDP_D3_expose_D(psi[i]->ptr);
            lua_pop(L, 1);
        }

        env.lat = S->lat;
        env.f = e_psi;
        if (QOP_MDWF_import_fermion(&c_psi, c->state, DW_5_reader, &env))
            return luaL_error(L, "MDWF_import_fermion() failed");

        for (i = 0; i < c->Ls; i++)
            QDP_D3_reset_D(psi[i]->ptr);

        if (QOP_MDWF_allocate_fermion(&c_eta, c->state))
            return luaL_error(L, "MDWF_create_fermion() failed");
        
        (*op)(c_eta, c->params, c->gauge, c_psi);

        QOP_MDWF_free_fermion(&c_psi);

        lua_createtable(L, c->Ls, 0);
        for (i = 0; i < c->Ls; i++) {
            eta[i] = qlua_newLatDirFerm3(L, Sidx, 3);
            e_eta[i] = QDP_D3_expose_D(eta[i]->ptr);
            lua_rawseti(L, -2, i + 1);
        }

        env.lat = S->lat;
        env.f = e_eta;
        if (QOP_MDWF_export_fermion(DW_5_writer, &env, c_eta))
            return luaL_error(L, "MDWF_export_fermion() failed");

        for (i = 0; i < c->Ls; i++) {
            QDP_D3_reset_D(eta[i]->ptr);
        }

        QOP_MDWF_free_fermion(&c_eta);

        return 1;
    }
    default:
        break;
    }
    return luaL_error(L, "bad arguments in Clover:%s", name);
}

static int
q_DW_D(lua_State *L)
{
    return q_DW_operator(L, "D", QOP_D3_MDWF_DDW_operator);
}

static int
q_DW_Dx(lua_State *L)
{
    return q_DW_operator(L, "Dx", QOP_D3_MDWF_DDW_operator_conjugated);
}

/* the standard clover solver */
static int
q_DW_std_solver(lua_State *L,
                struct QOP_MDWF_Fermion *solution,
                int *out_iters,
                double *out_epsilon,
                const struct QOP_MDWF_Fermion *rhs,
                int log_level)
{
    mMDWF *c = qlua_checkMDWF(L, lua_upvalueindex(2), NULL, 1);
    double eps = luaL_checknumber(L, lua_upvalueindex(3));
    int max_iters = luaL_checkint(L, lua_upvalueindex(4));
    
    return QOP_MDWF_DDW_CG(solution, out_iters, out_epsilon,
                           c->params, rhs, c->gauge, rhs, max_iters, eps,
                           log_level);
}

static MDWFSolver std_solver = { q_DW_std_solver, "CG" };

static int
q_DW_make_solver(lua_State *L)
{
    qlua_checkMDWF(L, 1, NULL, 1);   /* mClover */
    luaL_checknumber(L, 2);            /* double epsilon */
    luaL_checkint(L, 3);               /* int max_iter */

    lua_pushlightuserdata(L, &std_solver);  /* cl[1]: solver */
    lua_pushvalue(L, 1);                    /* cl[2]: mClover */
    lua_pushvalue(L, 2);                    /* cl[3]: epsilon */
    lua_pushvalue(L, 3);                    /* cl[4]: max_iters */
    lua_pushcclosure(L, q_dirac_solver, 4);

    return 1;
}

/* the mixed clover solver */
static int
q_DW_mixed_solver(lua_State *L,
                  struct QOP_MDWF_Fermion *solution,
                  int *out_iters,
                  double *out_epsilon,
                  const struct QOP_MDWF_Fermion *rhs,
                  int log_level)
{
    mMDWF *c = qlua_checkMDWF(L, lua_upvalueindex(2), NULL, 1);
    double f_eps    = luaL_checknumber(L, lua_upvalueindex(3));
    int inner_iters = luaL_checkint(L, lua_upvalueindex(4));
    double eps      = luaL_checknumber(L, lua_upvalueindex(5));
    int max_iters   = luaL_checkint(L, lua_upvalueindex(6));
    
    return QOP_MDWF_mixed_DDW_CG(solution, out_iters, out_epsilon,
                                 c->params, rhs, c->gauge, rhs,
                                 inner_iters, f_eps,
                                 max_iters, eps,
                                 log_level);
}

static MDWFSolver mixed_solver = { q_DW_mixed_solver, "mixedCG" };

static int
q_DW_make_mixed_solver(lua_State *L)
{
    qlua_checkMDWF(L, 1, NULL, 1);   /* mClover */
    luaL_checknumber(L, 2);            /* double inner_epsilon */
    luaL_checkint(L, 3);               /* int inner_iter */
    luaL_checknumber(L, 4);            /* double epsilon */
    luaL_checkint(L, 5);               /* int max_iter */

    lua_pushlightuserdata(L, &mixed_solver);  /* cl[1]: solver */
    lua_pushvalue(L, 1);                      /* cl[2]: mMDWF */
    lua_pushvalue(L, 2);                      /* cl[3]: inner_epsilon */
    lua_pushvalue(L, 3);                      /* cl[4]: inner_iters */
    lua_pushvalue(L, 4);                      /* cl[5]: epsilon */
    lua_pushvalue(L, 5);                      /* cl[6]: max_iters */
    lua_pushcclosure(L, q_dirac_solver, 6);

    return 1;
}

typedef struct {
    QDP_Lattice *lat;
    int lattice[QOP_MDWF_DIM];
    QLA_D_Complex bf[QOP_MDWF_DIM];
    QLA_D3_ColorMatrix *uf[QOP_MDWF_DIM];
} QCArgs;

static double
q_DW_u_reader(int d, const int p[], int a, int b, int re_im, void *env)
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

/* MDWF constructors helper:
 *    Lua(U[4], bc[4], Ls, ...) => mMDWF
 */
static mMDWF *
q_mdwf(lua_State *L)
{
    int i;
    int Ls = luaL_checkint(L, 3);

    luaL_checktype(L, 1, LUA_TTABLE);
    lua_pushnumber(L, 1);
    lua_gettable(L, 1);
    qlua_checkLatColMat3(L, -1, NULL, 3);
    mLattice *S = qlua_ObjLattice(L, -1);
    int Sidx = lua_gettop(L);
    mMDWF *c = qlua_newMDWF(L, Sidx);

    if (Ls < 1)
        luaL_error(L, "Bad value of Ls");
    if (S->rank != QOP_MDWF_DIM)
        luaL_error(L, "MDWF is not implemented for #L=%d", S->rank);
    if (QDP_Ns != QOP_MDWF_FERMION_DIM)
        luaL_error(L, "MDWF does not support Ns=%d", QDP_Ns);

    c->Ls = Ls;
    QCArgs args;
    luaL_checktype(L, 2, LUA_TTABLE);
    for (i = 0; i < QOP_MDWF_DIM; i++) {
        lua_pushnumber(L, i + 1);
        lua_gettable(L, 2);
        switch (qlua_qtype(L, -1)) {
        case qReal:
            QLA_c_eq_r_plus_ir(args.bf[i], lua_tonumber(L, -1), 0);
            break;
        case qComplex:
            QLA_c_eq_c(args.bf[i], *qlua_checkComplex(L, -1));
            break;
        default:
            luaL_error(L, "bad MDWF boundary condition type");
        }
        lua_pop(L, 1);
    }
    
    QDP_D3_ColorMatrix *UF[QOP_MDWF_DIM];

    luaL_checktype(L, 1, LUA_TTABLE);
    CALL_QDP(L);

    /* extract U from the arguments */
    for (i = 0; i < QOP_MDWF_DIM; i++) {
        UF[i] = QDP_D3_create_M_L(S->lat);
        lua_pushnumber(L, i + 1); /* [sic] lua indexing */
        lua_gettable(L, 1);
        /* avoid aliased Us in arg[1] */
        QDP_D3_M_eq_M(UF[i], qlua_checkLatColMat3(L, -1, S, 3)->ptr, S->all);
        args.uf[i] = QDP_D3_expose_M(UF[i]);
        lua_pop(L, 1);
    }
    args.lat = S->lat;
    QDP_latsize_L(S->lat, args.lattice);

    struct QOP_MDWF_Config cc;
    cc.self = S->node;
    cc.master_p = QMP_is_primary_node();
    cc.rank = S->rank;
    cc.lat = S->dim;
    cc.ls = Ls;
    cc.net = S->net;
    cc.neighbor_up = S->neighbor_up;
    cc.neighbor_down = S->neighbor_down;
    cc.sublattice = qlua_sublattice;
    cc.env = S;
    if (QOP_MDWF_init(&c->state, &cc))
        luaL_error(L, "MDWF_init() failed");
    
    if (QOP_MDWF_import_gauge(&c->gauge, c->state, q_DW_u_reader, &args))
        luaL_error(L, "MDWF_import_gauge() failed");

    for (i = 0; i < QOP_MDWF_DIM; i++) {
        QDP_D3_reset_M(UF[i]);
        QDP_D3_destroy_M(UF[i]);
    }
    return c;
}

/*
 * qcd.MDWF.generic(U[4],          -- [1] Gauge, gauge field
 *                  ubc[4],        -- [2] double[4], gauge boundary conditions
 *                  Ls,            -- [3] int, flavor dimension size
 *                  M5,            -- [4] double
 *                  mf,            -- [5] double
 *                  b5[Ls],        -- [6] double[Ls]
 *                  c5[Ls])        -- [7] double[Ls]
 */
static int
q_mdwf_generic(lua_State *L)
{
    mMDWF *M = q_mdwf(L);
    int Ls = luaL_checkint(L, 3);
    double M5 = luaL_checknumber(L, 4);
    double mf = luaL_checknumber(L, 5);
    double b5[Ls]; /* [6] */
    double c5[Ls]; /* [7] */
    int i;

    for (i = 0; i < Ls; i++) {
        lua_pushnumber(L, i + 1); /* [sic] lua indexing */
        lua_gettable(L, 6);
        b5[i] = luaL_checknumber(L, -1);
        lua_pushnumber(L, i + 1); /* [sic] lua indexing */
        lua_gettable(L, 7);
        c5[i] = luaL_checknumber(L, -1);
        lua_pop(L, 2);
    }

    M->name = "generic";
    if (QOP_MDWF_set_generic(&M->params, M->state, b5, c5, M5, mf))
        return luaL_error(L, "Not enough space");

    return 1;
}

/*
 * qcd.MDWF.Moebius(U[4],          -- [1] Gauge, gauge field
 *                  ubc[4],        -- [2] double[4], gauge boundary conditions
 *                  Ls,            -- [3] int, flavor dimension size
 *                  M5,            -- [4] double
 *                  mf,            -- [5] double
 *                  b5[Ls],        -- [6] double[Ls]
 *                  kappa)         -- [7] double
 */
static int
q_mdwf_Moebius(lua_State *L)
{
    mMDWF *M = q_mdwf(L);
    int Ls = luaL_checkint(L, 3);
    double M5 = luaL_checknumber(L, 4);
    double mf = luaL_checknumber(L, 5);
    double b5[Ls]; /* [6] */
    double kappa = luaL_checknumber(L, 7);
    int i;

    for (i = 0; i < Ls; i++) {
        lua_pushnumber(L, i + 1); /* [sic] lua indexing */
        lua_gettable(L, 6);
        b5[i] = luaL_checknumber(L, -1);
        lua_pop(L, 1);
    }

    M->name = "Moebius";
    if (QOP_MDWF_set_Moebius(&M->params, M->state, b5, kappa, M5, mf))
        return luaL_error(L, "Not enough space");

    return 1;
}

/*
 * qcd.MDWF.Shamir(U[4],           -- [1] Gauge, gauge field
 *                  ubc[4],        -- [2] double[4], gauge boundary conditions
 *                  Ls,            -- [3] int, flavor dimension size
 *                  M5,            -- [4] double
 *                  mf,            -- [5] double
 *                  a5)            -- [6] double
 */
static int
q_mdwf_Shamir(lua_State *L)
{
    mMDWF *M = q_mdwf(L);
    double M5 = luaL_checknumber(L, 4);
    double mf = luaL_checknumber(L, 5);
    double a5 = luaL_checknumber(L, 6);

    M->name = "Shamir";
    if (QOP_MDWF_set_Shamir(&M->params, M->state, a5, M5, mf))
        return luaL_error(L, "Not enough space");

    return 1;
}

/*
 * qcd.MDWF.Borichi(U[4],          -- [1] Gauge, gauge field
 *                  ubc[4],        -- [2] double[4], gauge boundary conditions
 *                  Ls,            -- [3] int, flavor dimension size
 *                  M5,            -- [4] double
 *                  mf,            -- [5] double
 *                  a5)            -- [6] double
 */
static int
q_mdwf_Borichi(lua_State *L)
{
    mMDWF *M = q_mdwf(L);
    double M5 = luaL_checknumber(L, 4);
    double mf = luaL_checknumber(L, 5);
    double a5 = luaL_checknumber(L, 6);

    M->name = "Borichi";
    if (QOP_MDWF_set_Borichi(&M->params, M->state, a5, M5, mf))
        return luaL_error(L, "Not enough space");

    return 1;
}

/*
 * qcd.MDWF.Chiu(U[4],          -- [1] Gauge, gauge field
 *               ubc[4],        -- [2] double[4], gauge boundary conditions
 *               Ls,            -- [3] int, flavor dimension size
 *               M5,            -- [4] double
 *               mf,            -- [5] double
 *               a5[Ls])        -- [6] double[Ls]
 */
static int
q_mdwf_Chiu(lua_State *L)
{
    mMDWF *M = q_mdwf(L);
    int Ls = luaL_checkint(L, 3);
    double M5 = luaL_checknumber(L, 4);
    double mf = luaL_checknumber(L, 5);
    double a5[Ls]; /* [6] */
    int i;

    for (i = 0; i < Ls; i++) {
        lua_pushnumber(L, i + 1); /* [sic] lua indexing */
        lua_gettable(L, 6);
        a5[i] = luaL_checknumber(L, -1);
        lua_pop(L, 1);
    }

    M->name = "Chiu";
    if (QOP_MDWF_set_Chiu(&M->params, M->state, a5, M5, mf))
        return luaL_error(L, "Not enough space");

    return 1;
}

static struct luaL_Reg mtMDWF[] = {
    { "__tostring",   q_DW_fmt                },
    { "__gc",         q_DW_gc                 },
    { "close",        q_DW_close              },
    { "D",            q_DW_D                  },
    { "Dx",           q_DW_Dx                 },
    { "solver",       q_DW_make_solver        },
    { "mixed_solver", q_DW_make_mixed_solver  },
    { "eig_deflator", q_DW_make_deflator      },
    { NULL,           NULL                    }
};

static mMDWF *
qlua_newMDWF(lua_State *L, int Sidx)
{
    mMDWF *c = lua_newuserdata(L, sizeof (mMDWF));

    c->state = 0;
    c->gauge = 0;
    c->params = 0;
    qlua_createLatticeTable(L, Sidx, mtMDWF, qMDWF, MDWFName);
    lua_setmetatable(L, -2);

    return c;
}

static mMDWF *
qlua_checkMDWF(lua_State *L, int idx, mLattice *S, int live)
{
    mMDWF *c = qlua_checkLatticeType(L, idx, qMDWF, MDWFName);

    if (S) {
        mLattice *S1 = qlua_ObjLattice(L, idx);
        if (S1->id != S->id)
            luaL_error(L, "%s on a wrong lattice", MDWFName);
        lua_pop(L, 1);
    }

    if (live && (c->state == 0 || c->gauge == 0))
        luaL_error(L, "using closed qcd.MDWF");

    return c;
}

static struct luaL_Reg fMDWF[] = {
    { "generic",       q_mdwf_generic },
    { "Moebius",       q_mdwf_Moebius },
    { "Shamir",        q_mdwf_Shamir  },
    { "Borichi",       q_mdwf_Borichi },
    { "Chiu",          q_mdwf_Chiu    },
    { NULL,            NULL           }
};

int
init_mdwf(lua_State *L)
{
    lua_getglobal(L, qcdlib);
    lua_newtable(L);
    luaL_register(L, NULL, fMDWF);
    lua_setfield(L, -2, mdwf_name);
    lua_pop(L, 1);
    return 0;
}

int
fini_mdwf(lua_State *L)
{
    return 0;
}
