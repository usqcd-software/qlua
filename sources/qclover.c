#include <qlua.h>                                                    /* DEPS */
#include <qclover.h>                                                 /* DEPS */
#include <qcomplex.h>                                                /* DEPS */
#include <lattice.h>                                                 /* DEPS */
#include <latcolmat.h>                                               /* DEPS */
#include <latdirferm.h>                                              /* DEPS */
#include <latdirprop.h>                                              /* DEPS */
#define QOP_CLOVER_DEFAULT_PRECISION QDP_Precision
#include <qop-clover.h>
#include <qmp.h>
#include <math.h>

/* NB: Clover operator does not agrees with BMW conventions */

static char mtnClover[] = "qcd.clover";

typedef struct {
    struct QOP_CLOVER_State *state;
    struct QOP_CLOVER_Gauge *gauge;
    double kappa, c_sw;
} mClover;

static mClover *
qlua_newClover(lua_State *L)
{
    mClover *c = lua_newuserdata(L, sizeof (mClover));

    c->state = 0;
    c->gauge = 0;
    luaL_getmetatable(L, mtnClover);
    lua_setmetatable(L, -2);

    return c;
}

static mClover *
qlua_checkClover(lua_State *L, int idx, int live)
{
    void *v = luaL_checkudata(L, idx, mtnClover);
    mClover *c = (mClover *)v;


    luaL_argcheck(L, v != 0, idx, "qcd.Clover expected");

    if (live && (c->state == 0 || c->gauge == 0))
        luaL_error(L, "using closed qcd.Clover");

    return v;
}

static int
q_CL_fmt(lua_State *L)
{
    char fmt[72];
    mClover *c = qlua_checkClover(L, 1, 0);

    if (c->state)
        sprintf(fmt, "Clover[%g,%g]", c->kappa, c->c_sw);
    else
        sprintf(fmt, "Clover(closed)");

    lua_pushstring(L, fmt);

    return 1;
}

static int
q_CL_gc(lua_State *L)
{
    mClover *c = qlua_checkClover(L, 1, 0);

    if (c->gauge) {
        QOP_CLOVER_free_gauge(&c->gauge);
        c->gauge = 0;
    }
    if (c->state) {
        QOP_CLOVER_fini(&c->state);
        c->state = 0;
    }
    
    return 0;
}

static int
q_CL_close(lua_State *L)
{
    mClover *c = qlua_checkClover(L, 1, 1);

    if (c->gauge) {
        QOP_CLOVER_free_gauge(&c->gauge);
        c->gauge = 0;
    }
    if (c->state) {
        QOP_CLOVER_fini(&c->state);
        c->state = 0;
    }
    
    return 0;
}

static double
q_CL_D_reader(const int p[], int c, int d, int re_im, void *env)
{
    int i = QDP_index(p);
    QLA_DiracFermion *f = env;
    QLA_Real xx;

    if (re_im == 0) {
        QLA_r_eq_Re_c(xx, QLA_elem_D(f[i], c, d));
    } else {
        QLA_r_eq_Im_c(xx, QLA_elem_D(f[i], c, d));
    }

    return xx;
}

static void
q_CL_D_writer(const int p[], int c, int d, int re_im, double v, void *env)
{
    int i = QDP_index(p);
    QLA_DiracFermion *f = env;
 
    if (re_im == 0) {
        QLA_real(QLA_elem_D(f[i], c, d)) = v;
    } else {
        QLA_imag(QLA_elem_D(f[i], c, d)) = v;
    }
}

struct CL_D_env {
    QLA_DiracFermion *f;
    double s;
};

static double
q_CL_D_reader_scaled(const int p[], int c, int d, int re_im, void *e)
{
    int i = QDP_index(p);
    struct CL_D_env *env = (struct CL_D_env *)e;
    QLA_DiracFermion *f = env->f;
    double s = env->s;
    QLA_Real xx;

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
    int i = QDP_index(p);
    struct CL_D_env *env = (struct CL_D_env *)e;
    QLA_DiracFermion *f = env->f;
    double s = env->s;
 
    v = v * s;
    if (re_im == 0) {
        QLA_real(QLA_elem_D(f[i], c, d)) = v;
    } else {
        QLA_imag(QLA_elem_D(f[i], c, d)) = v;
    }
}

typedef struct {
    int c, d;
    QLA_DiracPropagator *in;
    QLA_DiracPropagator *out;
    double s;
} qCL_P_env;

static double
q_CL_P_reader_scaled(const int p[], int c, int d, int re_im, void *env)
{
    int i = QDP_index(p);
    qCL_P_env *e = env;
    QLA_Real xx;

    if (re_im == 0) {
        QLA_r_eq_Re_c(xx, QLA_elem_P(e->in[i], c, d, e->c, e->d));
    } else {
        QLA_r_eq_Im_c(xx, QLA_elem_P(e->in[i], c, d, e->c, e->d));
    }

    return xx * e->s;
}

static void
q_CL_P_writer_scaled(const int p[], int c, int d, int re_im, double v, void *env)
{
    int i = QDP_index(p);
    qCL_P_env *e = env;

    v = v * e->s;
    if (re_im == 0) {
        QLA_real(QLA_elem_P(e->out[i], c, d, e->c, e->d)) = v;
    } else {
        QLA_imag(QLA_elem_P(e->out[i], c, d, e->c, e->d)) = v;
    }
}

static double
q_CL_P_reader(const int p[], int c, int d, int re_im, void *env)
{
    int i = QDP_index(p);
    qCL_P_env *e = env;
    QLA_Real xx;

    if (re_im == 0) {
        QLA_r_eq_Re_c(xx, QLA_elem_P(e->in[i], c, d, e->c, e->d));
    } else {
        QLA_r_eq_Im_c(xx, QLA_elem_P(e->in[i], c, d, e->c, e->d));
    }

    return xx;
}

static void
q_CL_P_writer(const int p[], int c, int d, int re_im, double v, void *env)
{
    int i = QDP_index(p);
    qCL_P_env *e = env;

    if (re_im == 0) {
        QLA_real(QLA_elem_P(e->out[i], c, d, e->c, e->d)) = v;
    } else {
        QLA_imag(QLA_elem_P(e->out[i], c, d, e->c, e->d)) = v;
    }
}
static int
q_CL_D(lua_State *L)
{
    mClover *c = qlua_checkClover(L, 1, 1);

    switch (qlua_gettype(L, 2)) {
    case qLatDirFerm: {
        mLatDirFerm *psi = qlua_checkLatDirFerm(L, 2);
        mLatDirFerm *eta = qlua_newLatDirFerm(L);
        struct QOP_CLOVER_Fermion *c_psi;
        struct QOP_CLOVER_Fermion *c_eta;
        QLA_DiracFermion *e_psi;
        QLA_DiracFermion *e_eta;

        CALL_QDP(L);
        e_psi = QDP_expose_D(psi->ptr);
        if (QOP_CLOVER_import_fermion(&c_psi, c->state, q_CL_D_reader, e_psi))
            return luaL_error(L, "CLOVER_import_fermion() failed");
        QDP_reset_D(psi->ptr);

        if (QOP_CLOVER_allocate_fermion(&c_eta, c->state))
            return luaL_error(L, "CLOVER_allocate_fermion() failed");

        QOP_CLOVER_D_operator(c_eta, c->gauge, c_psi);
        
        e_eta = QDP_expose_D(eta->ptr);
        QOP_CLOVER_export_fermion(q_CL_D_writer, e_eta, c_eta);
        QDP_reset_D(eta->ptr);
        
        QOP_CLOVER_free_fermion(&c_eta);
        QOP_CLOVER_free_fermion(&c_psi);

        return 1;
    }
    case qLatDirProp: {
        mLatDirProp *psi = qlua_checkLatDirProp(L, 2);
        mLatDirProp *eta = qlua_newLatDirProp(L);
        struct QOP_CLOVER_Fermion *c_psi;
        struct QOP_CLOVER_Fermion *c_eta;
        qCL_P_env env;

        CALL_QDP(L);
        if (QOP_CLOVER_allocate_fermion(&c_eta, c->state))
            return luaL_error(L, "CLOVER_allocate_fermion() failed");

        env.in = QDP_expose_P(psi->ptr);
        env.out = QDP_expose_P(eta->ptr);
        for (env.c = 0; env.c < QOP_CLOVER_COLORS; env.c++) {
            for (env.d = 0; env.d < QOP_CLOVER_FERMION_DIM; env.d++) {
                if (QOP_CLOVER_import_fermion(&c_psi, c->state,
                                              q_CL_P_reader, &env))
                    return luaL_error(L, "CLOVER_import_fermion failed");
                QOP_CLOVER_D_operator(c_eta, c->gauge, c_psi);
                QOP_CLOVER_free_fermion(&c_psi);
                QOP_CLOVER_export_fermion(q_CL_P_writer, &env, c_eta);
            }
        }
        QOP_CLOVER_free_fermion(&c_eta);
        QDP_reset_P(psi->ptr);
        QDP_reset_P(eta->ptr);

        return 1;
    }
    }
    return luaL_error(L, "bad arguments in Clover:D");
}

static int
q_CL_Dx(lua_State *L)
{
    mClover *c = qlua_checkClover(L, 1, 1);

    switch (qlua_gettype(L, 2)) {
    case qLatDirFerm: {
        mLatDirFerm *psi = qlua_checkLatDirFerm(L, 2);
        mLatDirFerm *eta = qlua_newLatDirFerm(L);
        struct QOP_CLOVER_Fermion *c_psi;
        struct QOP_CLOVER_Fermion *c_eta;
        QLA_DiracFermion *e_psi;
        QLA_DiracFermion *e_eta;

        CALL_QDP(L);
        e_psi = QDP_expose_D(psi->ptr);
        if (QOP_CLOVER_import_fermion(&c_psi, c->state, q_CL_D_reader, e_psi))
            return luaL_error(L, "CLOVER_import_fermion() failed");
        QDP_reset_D(psi->ptr);

        if (QOP_CLOVER_allocate_fermion(&c_eta, c->state))
            return luaL_error(L, "CLOVER_allocate_fermion() failed");

        QOP_CLOVER_D_operator_conjugated(c_eta, c->gauge, c_psi);
        e_eta = QDP_expose_D(eta->ptr);
        QOP_CLOVER_export_fermion(q_CL_D_writer, e_eta, c_eta);
        QDP_reset_D(eta->ptr);
     
        QOP_CLOVER_free_fermion(&c_eta);
        QOP_CLOVER_free_fermion(&c_psi);

        return 1;
    }
    case qLatDirProp: {
        mLatDirProp *psi = qlua_checkLatDirProp(L, 2);
        mLatDirProp *eta = qlua_newLatDirProp(L);
        struct QOP_CLOVER_Fermion *c_psi;
        struct QOP_CLOVER_Fermion *c_eta;
        qCL_P_env env;

        CALL_QDP(L);
        if (QOP_CLOVER_allocate_fermion(&c_eta, c->state))
            return luaL_error(L, "CLOVER_allocate_fermion() failed");

        env.in = QDP_expose_P(psi->ptr);
        env.out = QDP_expose_P(eta->ptr);
        for (env.c = 0; env.c < QOP_CLOVER_COLORS; env.c++) {
            for (env.d = 0; env.d < QOP_CLOVER_FERMION_DIM; env.d++) {
                if (QOP_CLOVER_import_fermion(&c_psi, c->state,
                                              q_CL_P_reader, &env))
                    return luaL_error(L, "CLOVER_import_fermion failed");
                QOP_CLOVER_D_operator_conjugated(c_eta, c->gauge, c_psi);
                QOP_CLOVER_free_fermion(&c_psi);
                QOP_CLOVER_export_fermion(q_CL_P_writer, &env, c_eta);
            }
        }
        QOP_CLOVER_free_fermion(&c_eta);
        QDP_reset_P(psi->ptr);
        QDP_reset_P(eta->ptr);

        return 1;
    }
    }
    return luaL_error(L, "bad arguments in Clover:Dx");
}

/* the solver */
static int
q_CL_solve(lua_State *L)
{
    mClover *c = qlua_checkClover(L, lua_upvalueindex(1), 1);
    double eps = luaL_checknumber(L, lua_upvalueindex(2));
    int max_iters = luaL_checkint(L, lua_upvalueindex(3));
    long long fl1;
    double t1;
    double out_eps;
    int out_iters;

    switch (qlua_gettype(L, 1)) {
    case qLatDirFerm: {
        mLatDirFerm *psi = qlua_checkLatDirFerm(L, 1);
        mLatDirFerm *eta = qlua_newLatDirFerm(L);
        struct QOP_CLOVER_Fermion *c_psi;
        struct QOP_CLOVER_Fermion *c_eta;
        struct CL_D_env env;
        QLA_Real rhs_norm2 = 0;
        double rhs_n;
        int status;

        CALL_QDP(L);
        QDP_r_eq_norm2_D(&rhs_norm2, psi->ptr, QDP_all);
        if (rhs_norm2 == 0) {
            QDP_D_eq_zero(eta->ptr, QDP_all);
            lua_pushnumber(L, 0.0);
            lua_pushnumber(L, 0);
            lua_pushnumber(L, 0);
            lua_pushnumber(L, 0);
            return 5;
        }
        rhs_n = sqrt(rhs_norm2);
        env.f = QDP_expose_D(psi->ptr);
        env.s = 1 / rhs_n;
        if (QOP_CLOVER_import_fermion(&c_psi, c->state, q_CL_D_reader_scaled, &env))
            return luaL_error(L, "CLOVER_import_fermion() failed");
        QDP_reset_D(psi->ptr);

        if (QOP_CLOVER_allocate_fermion(&c_eta, c->state))
            return luaL_error(L, "CLOVER_allocate_fermion() failed");

        status = QOP_CLOVER_D_CG(c_eta, &out_iters, &out_eps,
                                 c_psi, c->gauge, c_psi, max_iters, eps,
                                 0 /* QOP_CLOVER_FINAL_DIRAC_RESIDUAL */);
        QOP_CLOVER_performance(&t1, &fl1, NULL, NULL, c->state);
        if (qlua_primary_node)
            printf("CLOVER CG solver: status = %d,"
                   " eps = %.4e, iters = %d, time = %.3f sec,"
                   " perf = %.2f MFlops/sec\n",
                   status, out_eps, out_iters, t1, fl1 * 1e-6 / t1);

        if (status)
            return luaL_error(L, QOP_CLOVER_error(c->state));

        env.f = QDP_expose_D(eta->ptr);
        env.s = rhs_n;
        QOP_CLOVER_export_fermion(q_CL_D_writer_scaled, &env, c_eta);
        QDP_reset_D(eta->ptr);
        
        QOP_CLOVER_free_fermion(&c_eta);
        QOP_CLOVER_free_fermion(&c_psi);

        /* eta is on the stack already */
        lua_pushnumber(L, out_eps);
        lua_pushnumber(L, out_iters);
        lua_pushnumber(L, t1);
        lua_pushnumber(L, (double)fl1);

        return 5;
    }
    case qLatDirProp: {
        mLatDirProp *psi = qlua_checkLatDirProp(L, 1);
        mLatDirProp *eta = qlua_newLatDirProp(L);
        struct QOP_CLOVER_Fermion *c_psi;
        struct QOP_CLOVER_Fermion *c_eta;
        int status;
        qCL_P_env env;
        QLA_Real rhs_norm2 = 0;
        double rhs_n;

        CALL_QDP(L);
        if (QOP_CLOVER_allocate_fermion(&c_eta, c->state))
            return luaL_error(L, "CLOVER_allocate_fermion() failed");

        QDP_r_eq_norm2_P(&rhs_norm2, psi->ptr, QDP_all);
        if (rhs_norm2 == 0) {
            return 1;
        }
        rhs_n = sqrt(rhs_norm2);
        env.in = QDP_expose_P(psi->ptr);
        env.out = QDP_expose_P(eta->ptr);
        lua_createtable(L, QOP_CLOVER_COLORS, 0);  /* eps */
        lua_createtable(L, QOP_CLOVER_COLORS, 0);  /* iters */
        for (env.c = 0; env.c < QOP_CLOVER_COLORS; env.c++) {
            lua_createtable(L, QOP_CLOVER_FERMION_DIM, 0); /* eps.c */
            lua_createtable(L, QOP_CLOVER_FERMION_DIM, 0); /* iters.c */

            for (env.d = 0; env.d < QOP_CLOVER_FERMION_DIM; env.d++) {
                env.s = 1 / rhs_n;
                if (QOP_CLOVER_import_fermion(&c_psi, c->state,
                                              q_CL_P_reader_scaled, &env))
                    return luaL_error(L, "CLOVER_import_fermion() failed");
                status = QOP_CLOVER_D_CG(c_eta, &out_iters, &out_eps,
                                         c_psi, c->gauge, c_psi, max_iters, eps,
                                         0);
                /* QOP_CLOVER_FINAL_DIRAC_RESIDUAL); */

                QOP_CLOVER_performance(&t1, &fl1, NULL, NULL, c->state);
                
                if (qlua_primary_node)
                    printf("CLOVER CG solver: status = %d, c = %d, d = %d,"
                           " eps = %.4e, iters = %d, time = %.3f sec,"
                           " perf = %.2f MFlops/sec\n",
                           status, env.c, env.d, out_eps, out_iters, t1,
                           fl1 * 1e-6 / t1);
                QOP_CLOVER_free_fermion(&c_psi);
                if (status)
                    return luaL_error(L, QOP_CLOVER_error(c->state));

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
        QDP_reset_P(psi->ptr);
        QDP_reset_P(eta->ptr);
        QOP_CLOVER_free_fermion(&c_eta);

        return 3;
    }
    }
    return luaL_error(L, "bad argument to CLOVER solver");
}

static int
q_CL_make_solver(lua_State *L)
{

    qlua_checkClover(L, 1, 1);   /* mClover *c */
    luaL_checknumber(L, 2);   /* double epsilon */
    luaL_checkint(L, 3);      /* int max_iter */
    lua_pushcclosure(L, q_CL_solve, 3);

    return 1;
}

/* the solver */
static int
q_CL_mixed_solve(lua_State *L)
{
    mClover *c      = qlua_checkClover(L, lua_upvalueindex(1), 1);
    double eps      = luaL_checknumber(L, lua_upvalueindex(2));
    int inner_iters = luaL_checkint(L, lua_upvalueindex(3));
    int max_iters   = luaL_checkint(L, lua_upvalueindex(4));
    long long fl1;
    double t1;
    double out_eps;
    int out_iters;

    switch (qlua_gettype(L, 1)) {
    case qLatDirFerm: {
        mLatDirFerm *psi = qlua_checkLatDirFerm(L, 1);
        mLatDirFerm *eta = qlua_newLatDirFerm(L);
        struct QOP_CLOVER_Fermion *c_psi;
        struct QOP_CLOVER_Fermion *c_eta;
        struct CL_D_env env;
        QLA_Real rhs_norm2;
        double rhs_n;
        int status;

        CALL_QDP(L);
        QDP_r_eq_norm2_D(&rhs_norm2, psi->ptr, QDP_all);
        if (rhs_norm2 == 0) {
            lua_pushnumber(L, 0);
            lua_pushnumber(L, 0);
            lua_pushnumber(L, 0);
            lua_pushnumber(L, 0);
            return 5;
        }
        rhs_n = sqrt(rhs_norm2);
        env.f = QDP_expose_D(psi->ptr);
        env.s = 1 / rhs_n;
        if (QOP_CLOVER_import_fermion(&c_psi, c->state, q_CL_D_reader_scaled, &env))
            return luaL_error(L, "CLOVER_import_fermion() failed");
        QDP_reset_D(psi->ptr);

        if (QOP_CLOVER_allocate_fermion(&c_eta, c->state))
            return luaL_error(L, "CLOVER_allocate_fermion() failed");

        status = QOP_CLOVER_mixed_D_CG(c_eta, &out_iters, &out_eps,
                                       c_psi, c->gauge, c_psi,
                                       inner_iters, max_iters, eps,
                                       0);
        /* QOP_CLOVER_FINAL_DIRAC_RESIDUAL); */
        QOP_CLOVER_performance(&t1, &fl1, NULL, NULL, c->state);
        if (qlua_primary_node)
            printf("CLOVER mCG solver: status = %d,"
                   " eps = %.4e, iters = %d, time = %.3f sec,"
                   " perf = %.2f MFlops/sec\n",
                   status, out_eps, out_iters, t1, fl1 * 1e-6 / t1);

        if (status)
            return luaL_error(L, QOP_CLOVER_error(c->state));

        env.f = QDP_expose_D(eta->ptr);
        env.s = rhs_n;
        QOP_CLOVER_export_fermion(q_CL_D_writer_scaled, &env, c_eta);
        QDP_reset_D(eta->ptr);
        
        QOP_CLOVER_free_fermion(&c_eta);
        QOP_CLOVER_free_fermion(&c_psi);

        /* eta is on the stack already */
        lua_pushnumber(L, out_eps);
        lua_pushnumber(L, out_iters);
        lua_pushnumber(L, t1);
        lua_pushnumber(L, (double)fl1);

        return 5;
    }
    case qLatDirProp: {
        mLatDirProp *psi = qlua_checkLatDirProp(L, 1);
        mLatDirProp *eta = qlua_newLatDirProp(L);
        struct QOP_CLOVER_Fermion *c_psi;
        struct QOP_CLOVER_Fermion *c_eta;
        QLA_Real rhs_norm2;
        double rhs_n;
        int status;
        qCL_P_env env;

        CALL_QDP(L);
        if (QOP_CLOVER_allocate_fermion(&c_eta, c->state))
            return luaL_error(L, "CLOVER_allocate_fermion() failed");

        QDP_r_eq_norm2_P(&rhs_norm2, psi->ptr, QDP_all);
        if (rhs_norm2 == 0) {
            return 1;
        }
        rhs_n = sqrt(rhs_norm2);
        env.in = QDP_expose_P(psi->ptr);
        env.out = QDP_expose_P(eta->ptr);
        lua_createtable(L, QOP_CLOVER_COLORS, 0);  /* eps */
        lua_createtable(L, QOP_CLOVER_COLORS, 0);  /* iters */
        for (env.c = 0; env.c < QOP_CLOVER_COLORS; env.c++) {
            lua_createtable(L, QOP_CLOVER_FERMION_DIM, 0); /* eps.c */
            lua_createtable(L, QOP_CLOVER_FERMION_DIM, 0); /* iters.c */

            for (env.d = 0; env.d < QOP_CLOVER_FERMION_DIM; env.d++) {
                env.s = 1 / rhs_n;
                if (QOP_CLOVER_import_fermion(&c_psi, c->state,
                                              q_CL_P_reader_scaled, &env))
                    return luaL_error(L, "CLOVER_import_fermion() failed");
                status = QOP_CLOVER_mixed_D_CG(c_eta, &out_iters, &out_eps,
                                               c_psi, c->gauge, c_psi,
                                               inner_iters, max_iters, eps,
                                               0);
                /* QOP_CLOVER_FINAL_DIRAC_RESIDUAL); */
                QOP_CLOVER_performance(&t1, &fl1, NULL, NULL, c->state);
                
                if (qlua_primary_node)
                    printf("CLOVER mCG solver: status = %d, c = %d, d = %d,"
                           " eps = %.4e, iters = %d, time = %.3f sec,"
                           " perf = %.2f MFlops/sec\n",
                           status, env.c, env.d, out_eps, out_iters, t1,
                           fl1 * 1e-6 / t1);
                QOP_CLOVER_free_fermion(&c_psi);
                if (status)
                    return luaL_error(L, QOP_CLOVER_error(c->state));

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
        QDP_reset_P(psi->ptr);
        QDP_reset_P(eta->ptr);
        QOP_CLOVER_free_fermion(&c_eta);

        return 3;
    }
    }
    return luaL_error(L, "bad argument to CLOVER solver");
}

static int
q_CL_make_mixed_solver(lua_State *L)
{

    qlua_checkClover(L, 1, 1);   /* mClover *c */
    luaL_checknumber(L, 2);   /* double epsilon */
    luaL_checkint(L, 3);      /* int inner_iter */
    luaL_checkint(L, 4);      /* int max_iter */
    lua_pushcclosure(L, q_CL_mixed_solve, 4);

    return 1;
}

#define Nu  QOP_CLOVER_DIM
#define Nf  ((QOP_CLOVER_DIM * (QOP_CLOVER_DIM - 1)) / 2)
#define Nt  (Nu + Nf)
#define Nz  (Nu + Nf + 6)

typedef struct {
    int lattice[QOP_CLOVER_DIM];
    int network[QOP_CLOVER_DIM];
    QLA_Complex bf[QOP_CLOVER_DIM];
    QLA_ColorMatrix *uf[Nu + Nf];
} QCArgs;

static void
q_clover_sublattice(int lo[], int hi[], const int node[], void *env)
{
    QCArgs *args = env;
    int i;

    for (i = 0; i < QOP_CLOVER_DIM; i++) {
        lo[i] = (args->lattice[i] * node[i]) / args->network[i];
        hi[i] = (args->lattice[i] * (node[i] + 1)) / args->network[i];
    }
}

static void
get_vector(int v[], int def, int dim, const int d[])
{
    int i;

    for (i = 0; i < dim && i < QOP_CLOVER_DIM; i++)
        v[i] = d[i];
    for (;i < QOP_CLOVER_DIM; i++)
        v[i] = def;
}

static double
q_CL_u_reader(int d, const int p[], int a, int b, int re_im, void *env)
{
    QLA_Complex z;
    QCArgs *args = env;
    int i = QDP_index(p);

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
    QLA_Real xx;
    QCArgs *args = env;
    int i = QDP_index(p);
    int d, xm, xn;

    for (d = 0, xm = 0; xm < QOP_CLOVER_DIM; xm++) {
        for (xn = xm + 1; xn < QOP_CLOVER_DIM; xn++, d++) {
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

static int
q_clover(lua_State *L)
{
    mClover *c = qlua_newClover(L);
    double kappa = luaL_checknumber(L, 2);
    double c_sw = luaL_checknumber(L, 3);
    QDP_ColorMatrix *UF[Nz];
    int mu, nu, i;
    QCArgs args;
    int node[QOP_CLOVER_DIM];

    luaL_checktype(L, 4, LUA_TTABLE);
    for (i = 0; i < QOP_CLOVER_DIM; i++) {
        lua_pushnumber(L, i + 1);
        lua_gettable(L, 4);
        switch (qlua_gettype(L, -1)) {
        case qReal:
            QLA_c_eq_r_plus_ir(args.bf[i], luaL_checknumber(L, -1), 0);
            break;
        case qComplex:
            QLA_c_eq_c(args.bf[i], *qlua_checkComplex(L, -1));
            break;
        default:
            luaL_error(L, "bad clover boundary condition type");
        }
        lua_pop(L, 1);
    }

    luaL_checktype(L, 1, LUA_TTABLE);

    if (qRank != QOP_CLOVER_DIM)
        return luaL_error(L, "clover is not implemented for #L=%d", qRank);
    if (QDP_Nc != QOP_CLOVER_COLORS)
        return luaL_error(L, "clover does not support Nc=%d", QDP_Nc);
    if (QDP_Ns != QOP_CLOVER_FERMION_DIM)
        return luaL_error(L, "clover does not support Ns=%d", QDP_Ns);

    c->kappa = kappa;
    c->c_sw = c_sw;

    CALL_QDP(L);

    /* create a temporary U, F, and temp M */
    for (i = 0; i < Nz; i++)
        UF[i] = QDP_create_M();

    /* extract U from the arguments */
    for (i = 0; i < QOP_CLOVER_DIM; i++) {
        lua_pushnumber(L, i + 1); /* [sic] lua indexing */
        lua_gettable(L, 1);
        /* avoid aliased Us in arg[1] */
        QDP_M_eq_M(UF[i], qlua_checkLatColMat(L, -1)->ptr, QDP_all);
        lua_pop(L, 1);
    }

    CALL_QDP(L);
    /* compute 8i*F[mu,nu] in UF[Nf...] */
    for (i = 0, mu = 0; mu < QOP_CLOVER_DIM; mu++) {
        for (nu = mu + 1; nu < QOP_CLOVER_DIM; nu++, i++) {
            /* clover in [mu, nu] --> UF[Nu + i] */
            QDP_M_eq_sM(UF[Nt], UF[nu], QDP_neighbor[mu], QDP_forward,
                        QDP_all);
            QDP_M_eq_Ma_times_M(UF[Nt+1], UF[nu], UF[mu], QDP_all);
            QDP_M_eq_M_times_M(UF[Nt+2], UF[Nt+1], UF[Nt], QDP_all);
            QDP_M_eq_sM(UF[Nt+3], UF[Nt+2], QDP_neighbor[nu], QDP_backward,
                        QDP_all);
            QDP_M_eq_M_times_Ma(UF[Nt+4], UF[Nt+3], UF[mu], QDP_all);
            QDP_M_eq_sM(UF[Nt+1], UF[mu], QDP_neighbor[nu], QDP_forward,
                        QDP_all);
            QDP_M_eq_Ma_times_M(UF[Nt+5], UF[mu], UF[Nt+3], QDP_all);
            QDP_M_eq_M_times_Ma(UF[Nt+2], UF[Nt], UF[Nt+1], QDP_all);
            QDP_M_eq_M_times_Ma(UF[Nt+3], UF[Nt+2], UF[nu], QDP_all);
            QDP_M_peq_M_times_M(UF[Nt+4], UF[mu], UF[Nt+3], QDP_all);
            QDP_M_peq_M_times_M(UF[Nt+5], UF[Nt+3], UF[mu], QDP_all);
            QDP_M_eq_sM(UF[Nt+2], UF[Nt+5], QDP_neighbor[mu], QDP_backward,
                        QDP_all);
            QDP_M_peq_M(UF[Nt+4], UF[Nt+2], QDP_all);
            QDP_M_eq_M(UF[Nu+i], UF[Nt+4], QDP_all);
            QDP_M_meq_Ma(UF[Nu+i], UF[Nt+4], QDP_all);
        }
    }

    /* create the clover state */
    get_vector(args.network, 1, QMP_get_logical_number_of_dimensions(),
               QMP_get_logical_dimensions());
    get_vector(node, 0, QMP_get_logical_number_of_dimensions(),
               QMP_get_logical_coordinates());
    QDP_latsize(args.lattice);
    if (QOP_CLOVER_init(&c->state, args.lattice, args.network, node,
                        QMP_is_primary_node(), q_clover_sublattice, &args))
        return luaL_error(L, "CLOVER_init() failed");
    
    /* import the gauge field */
    for (i = 0; i < Nt; i++) {
        /* NB: QDP requires all UF to be distinct */ 
        args.uf[i] = QDP_expose_M(UF[i]);
    }

    if (QOP_CLOVER_import_gauge(&c->gauge, c->state, kappa, c_sw,
                                q_CL_u_reader, q_CL_f_reader, &args)) {
        return luaL_error(L, "CLOVER_import_gauge() failed");
    }

    for (i = 0; i < Nt; i++)
        QDP_reset_M(UF[i]);

    /* clean up temporaries */
    for (i = 0; i < Nz; i++)
        QDP_destroy_M(UF[i]);

    return 1;
}

static struct luaL_Reg mtClover[] = {
    { "__tostring",   q_CL_fmt },
    { "__gc",         q_CL_gc },
    { "close",        q_CL_close },
    { "D",            q_CL_D },
    { "Dx",           q_CL_Dx },
    { "solver",       q_CL_make_solver },
    { "mixed_solver", q_CL_make_mixed_solver },
    { NULL, NULL }
};

static struct luaL_Reg fClover[] = {
    { "Clover",       q_clover },
    { NULL,           NULL }
};

int
init_clover(lua_State *L)
{
    luaL_register(L, qcdlib, fClover);
    qlua_metatable(L, mtnClover, mtClover);

    return 0;
}

int
fini_clover(lua_State *L)
{
    return 0;
}
