#include "modules.h"                                                 /* DEPS */
#include "qlua.h"                                                    /* DEPS */
#include "qmdwf.h"                                                   /* DEPS */
#include "lattice.h"                                                 /* DEPS */ 
#include "latdirferm.h"                                              /* DEPS */
#include "qop-mdwf3.h"


/* use aliases to avoid FORTRAN symbol renaming */

typedef struct {
    struct QOP_MDWF_State       *mdwf_state;
    struct QOP_MDWF_Parameters  *mdwf_param;
    struct QOP_F3_MDWF_Gauge    *mdwf_gauge;
    QOP_F3_MDWF_HalfFermion     *x, y; /* must be allocated before 
                                          calling mdwf_eoprec_op */

    /* polynomial acc parameters */
    /* TODO switch between Chebyshev/etc polynomials? */
    int                         poly_n;
    double                      poly_a, poly_b;
} op_MDWF_F3_eoprec_MdagM_arg_s;

void
op_MDWF_F3_eoprec_MdagM_op(
        int loc_dim,
        float complex *x, 
        float complex *y, 
        void *op_arg) /* x<-op(y) */
{
    op_MDWF_F3_eoprec_MdagM_arg_s *a = (op_MDWF_F3_eoprec_MdagM_arg_s *)op_arg;
    assert(2 * loc_dim == QOP_MDWF_half_fermion_size(a->mdwf_state));

    QOP_F3_MDWF_half_fermion_from_blas(a->y, (float *)y, 2 * loc_dim);
    if (0 < a->poly_n) {
        /* TODO polynomial acc */
    } else {
        QOP_F3_MDWF_M_operator(
                a->x, a->mdwf_param, a->mdwf_gauge, a->y);
        QOP_F3_MDWF_M_operator_conjugated(
                a->y, a->mdwf_param, a->mdwf_gauge, a->x);
    }
    QOP_F3_MDWF_blas_from_half_fermion((float *)x, 2 * loc_dim, a->y);
}
extern int 
pcnaupd(int             *COMM,
        int             *IDO, 
        char            *BMAT, 
        int             *N, 
        char            *WHICH, 
        int             *NEV, 
        float           *TOL, 
        float complex   *RESID, 
        int             *NCV, 
        float complex   *V, 
        int             *LDV, 
        int             *IPARAM, 
        int             *IPNTR, 
        float complex   *WORKD,
        float complex   *WORKL, 
        int             *LWORKL, 
        float           *RWORK, 
        int             *INFO, 
        short           bmat_len, 
        short           which_len);
extern int 
pcneupd(int             *COMM, 
        int             *RVEC, 
        char            *HOWMNY,
        int             *SELECT, 
        float complex   *D, 
        float complex   *Z, 
        int             *LDZ, 
        float complex   *SIGMA, 
        float complex   *WORKEV, 
        char            *BMAT, 
        int             *N, 
        char            *WHICH, 
        int             *NEV, 
        float           *TOL, 
        float complex   *RESID, 
        int             *NCV, 
        float complex   *V, 
        int             *LDV, 
        int             *IPARAM, 
        int             *IPNTR, 
        float complex   *WORKD, 
        float complex   *WORKL,
        int             *LWORKL, 
        float           *RWORK, 
        int             *INFO, 
        short           howmny_len,
        short           bmat_len, 
        short           which_len);

#define PARPACK_CNAUPD  pcnaupd
#define PARPACK_CNEUPD  pcneupd

int 
lanczos_internal_float(
        lua_State *L,
        MPI_Comm *mpi_comm,
        void (*op)(int loc_dim, 
                   float complex *x,
                   float complex *y,
                   void *op_arg),   /* x<-Op(y) */
        void *op_arg,
        int loc_dim, 
        int nev,
        int ncv,
        int max_iter,
        float tol,
        float complex **evec,   /* return buffer for evectors, [nev, n] */
        float complex **eval,   /* return buffer for evalues,  [nev] */
        int *n_iters,           /* return the iteration count */
        int *nconv              /* number of converged evals/evecs */
        )
{
    MPI_Fint mpi_comm_f = MPI_Comm_c2f(mpi_comm);
    CALL_QDP(L);

    /* all FORTRAN communication uses underscored variables */
    int ido_info_[2];
    int iparam_[11];
    int ipntr_[14];
    int n_      = loc_dim,
        nev_    = nev,
        ncv_    = ncv,
        ldv_    = loc_dim;
    float complex *resid_      = malloc(sizeof(float complex) * n_);
    memset(resid_, 0, sizeof(float complex) * n_);
    float complex *w_v_        = malloc(sizeof(float complex) * ldv_ * ncv_);
    float complex *w_workd_    = malloc(sizeof(float complex) * n_ * 3);
    memset(w_workd_, 0, sizeof(float complex) * n_ * 3);
    int lworkl_ = (3 * ncv_ * ncv_ + 5 * ncv_) * 2;
    float complex *w_workl_    = malloc(sizeof(float complex) * lworkl_);
    float *w_rwork_            = malloc(sizeof(float) *ncv_);
    /* __neupd-only workspace */
    float complex *w_d_ = malloc(sizeof(float complex) * (nev_ + 1));
    int rvec_ = 1;
    float complex sigma_ = 0;
    float complex *w_workev_ = malloc(sizeof(float complex) * 2 * ncv_);
    int *select_ = malloc(sizeof(int) * ncv_);
#define lanczosC_free_workspace do {\
    free(resid_);\
    free(w_v_);\
    free(w_workd_);\
    free(w_workl_);\
    free(w_rwork_);\
    free(w_d_);\
    free(w_workev_);\
    free(select_);\
} while(0)

    if(NULL == resid_ ||
           NULL == w_v_ ||
           NULL == w_workd_ ||
           NULL == w_workl_ ||
           NULL == w_rwork_ ||
           NULL == w_d_ ||
           NULL == w_workev_ ||
           NULL == select_) {
        return luaL_error(L, "cannot allocate workspace for CN*UPD");
        
    /* cnaupd cycle */
    ido_        = 0;
    info_       = 0;
    iparam_[0]  = 1;
    iparam_[2]  = max_iter;
    iparam_[3]  = 1;
    iparam_[6]  = 1;
    int iter_cnt= 0;
    do {
        PARPACK_CNAUPD(&mpi_comm_f, &ido_, "I", &n_, "SR",
                       &nev_, &tol, resid_, &ncv_, w_v_,
                       &ldv_, iparam_, ipntr_, w_workd_, w_workl_,
                       &lworkl_, w_rwork_, &info_, 1, 2);
        if (info_ < 0 || 1 < info_) {
            lanczosC_free_workspace;
            return luaL_error(L, "CNAUPD returned INFO=%d", info_);
        }
        if (99 == ido_ || 1 == info_)
            break;
        if (-1 == ido_ || 1 == ido_)
            op(n_, w_workd_ + ipntr_[1] - 1, w_workd_ + ipntr_[0] - 1,
               op_arg);
        else {
            lanczosC_free_workspace;
            return luaL_error(L, "CNAUPD returned IDO=%d", ido_);
        }

        iter_cnt++;
    } while (99 != ido_ && iter_cnt < max_iter);

    if (0 == info_) {
        assert(iparam[4] == nev);
        *nconv = iparam[4];
    } else if (1 == info_) {
        *nconv = iparam[4];
    } else {
        lanczosC_free_workspace;
        return luaL_error(L, "CNAUPD returned INFO=%d", info_);
    }

    /* for howmny="P", no additional space is required */
    PARPACK_CNEUPD(&mpi_comm_f, &rvec_, "P", select_, w_d_,
                   w_v_, &ldv_, &sigma_, w_workev_, "I", &n_, "SR",
                   &nev_, &tol_, resid_, &ncv_, w_v_,
                   &ldv_, iparam_, ipntr_, w_workd_, w_workl_,
                   &lworkl_, w_rwork_, &info_, 1, 1, 2);

    if (0 != info_) {
        lanczosC_free_workspace;
        return luaL_error(L, "CNEUPD returned INFO=%d", info_);
    }
    
    /* cleanup */
    free(resid_);
    free(w_workd_);
    free(w_workl_);
    free(w_rwork_);
    free(w_workev_);
    free(select_);

    /* free extra memory in the eigenspace buffer */
    assert(ldv_ == loc_dim);
    *evec = realloc(w_v_, sizeof(float complex) * ldv_ * nev_);

    /* copy real part of eigenvalues */
    *eval = malloc(sizeof(float) * nev_);
    for (i = 0 ; i < nev_ ; i++)
        (*eval)[i] = creal(w_d_[i]);
    free(w_d_);

    *n_iters    = iter_cnt;

    return 0;
#undef lanczosC_free_workspace
}


/* TODO move to qmdwf.c ? */
/* MDWF:eig_deflator_lanczos(nev, ncv, max_iter, tol, [param])
   param = {
     [cheb_accel]={n, a, b}, 
     [precision] = "single"|"double",
     ... } */

static int 
q_DW_make_deflator_lanczos(lua_State *L)
{
    mMDWF *c = qlua_checkMDWF(L, 1, NULL, 1);
    /* TODO parse nev, ncv, max_iter, tol */
    int nev, ncv, max_iter;
    double tol;
    nev     = lua_checkint(L, 2);
    ncv     = lua_checkint(L, 3);
    max_iter= lua_checkint(L, 4);
    tol     = luaL_checknumber(L, 5);
    int cheb_n = -1;
    double cheb_a = 0,
           cheb_b = 0;
    int eigcg_vmax  = 0,
        eigcg_umax  = 0,
        eigcg_nev   = 0;
    double eigcg_eps= 0.;

    if (qlua_checkopt_table(L, 6)) {
        if (qlua_tabpushopt_key(L, 6, "cheb_accel")) {
            cheb_n  = qlua_tabidx_int(L, -1, 1);
            cheb_a  = qlua_tabidx_double(L, -1, 2);
            cheb_b  = qlua_tabidx_double(L, -1, 3);
            lua_pop(L, 1);
        }
        if (qlua_tabpushopt_key(L, 6, "eigcg")) {
            eigcg_vmax  = qlua_tabidx_int(L, -1, 1);
            eigcg_nev   = qlua_tabidx_int(L, -1, 2);
            eigcg_eps   = qlua_tabidx_double(L, -1, 3);
            eigcg_umax  = qlua_tabidx_int(L, -1, 4);
            lua_pop(L, 1);
        }
    }

    /* single precision for lanczos */
    op_MDWF_F3_eoprec_MdagM_arg_s op_arg;
    op_arg.mdwf_state   = c->state;
    op_arg.mdwf_params  = c->params;
    op_arg.mdwf_gauge   = c->gauge;
    if (QOP_F3_MDWF_allocate_half_fermion(&op_arg.x, c->state))
        return luaL_error(L, "cannot allocate HalfFermion");
    if (QOP_F3_MDWF_allocate_half_fermion(&op_arg.y, c->state))
        return luaL_error(L, "cannot allocate HalfFermion");

    if (0 < cheb_n) {
        op_arg.poly_n   = cheb_n;
        op_arg.poly_a   = cheb_a;
        op_arg.poly_b   = cheb_b;
    } else
        op_arg.poly_n = -1;
    
    /* TODO get MPI_Comm; use WORLD for now? */
    MPI_Comm mpi_comm = ;
    
    int loc_dim = QOP_MDWF_half_fermion_size(c->state) / 2;
    assert(0 == QOP_MDWF_half_fermion_size(c->state) % 2);

    float complex *evec,
                  *eval;
    int n_iters, nconv;
    int status = lanczos_internal_float(
            L, &mpi_comm, op_MDWF_F3_eoprec_MdagM_op, &op_arg, 
            loc_dim, nev, ncv, max_iter, 
            tol, &evec, &eval, &n_iters
            &nconv);
    if (0 != status)
        return luaL_error(L, "lanczos_internal_float returned %d", status);

    /* TODO create empty deflator */
    /* mDeflatorState *d = q_newDeflatorState(L, Sidx); */

    for (int i = 0 ; i < nconv ; i++) {
        QOP_F3_MDWF_half_fermion_from_blas(&op_arg.x, 
                (float *)(evec + i * loc_dim), 2 * loc_dim);
        QOP_F3_MDWF_deflator_add_vector(d->deflator, &op_arg.x);
    }

    /* cleanup */
    QOP_F3_MDWF_allocate_half_fermion(&op_arg.x);
    QOP_F3_MDWF_allocate_half_fermion(&op_arg.y);
    free(evec);
    free(eval);
    
    /* TODO return something Lua-meaningful */
}


#if 0

/* TODO put double version into a separate file */
/* TODO 'double' needs revision according to changes in 'single' */
/* use aliases to avoid FORTRAN symbol renaming */
#define PARPACK_ZNAUPD  pznaupd
#define PARPACK_ZNEUPD  pzneupd
/* ARPACK function prototypes */
extern int
pznaupd(int *COMM,
        int  *IDO,
        char *BMAT,
        int  *N,
        char *WHICH,
        int  *NEV,
        double *TOL,
        double complex *RESID,  /* [N] */
        int *NCV,
        double complex *V, /* (N,NCV) */
        int *LDV,
        int *IPARAM, /* [11] */
        int *IPNTR, /* [14] */
        double complex *WORKD, /* [3*N] */
        double complex *WORKL, /* [LWORKL] */
        int *LWORKL, /* >=3*NCV**2 + 5*NCV */
        double *RWORK, /* [NCV] */
        int *INFO,
        int bmat_len,
        int which_len);
extern int
pzneupd(int *COMM,
        int *RVEC,      /* compute Ritz vectors? */
        char *HOWMNY,   /* A/P/S */
        int *SELECT,    /* [NCV] */
        double complex *D, /* (NEV+1) */
        double complex *Z, /* (N, NEV) ; may be equal to V */
        int *LDZ,
        double complex *SIGMA,
        double complex *WORKEV, /* (2*NCV) */
        char *BMAT,     /*+ */
        int *N,         /*+ */
        char *WHICH,    /*+ */
        int *NEV,       /*+ */
        double *TOL,    /*+ */
        double complex *RESID,  /*+ [N] */
        int *NCV, /*+ */
        double complex *V, /*+ (N, NCV) */
        int *LDV, /*+ */
        int *IPARAM, /*+ */
        int *IPNTR, /*+ */
        double complex *WORKD, /*+ */
        double complex *WORKL, /*+ */
        int *LWORKL, /*+ */
        double *RWORK, /*+ */
        int *INFO, /*+ */
        int howmny_len,
        int bmat_len,
        int which_len);


/* complex ILAM, double precision ; 
   the operator is assumed to be Hermitian: modifications may be required 
   for a general operator */
int 
lanczos_internal_double(
        lua_State *L,
        MPI_Comm *mpi_comm,
        void (*op)(int loc_dim, 
                   double complex *x,
                   double complex *y,
                   void *op_arg),   /* x<-Op(y) */
        void *op_arg,
        int loc_dim, 
        int nev,
        int ncv,
        int max_iter,
        double tol,
        double complex **evec,  /* return buffer for evectors, [nev, n] */
        double complex **eval,  /* return buffer for evalues,  [nev] */
        int *n_iters           /* return the iteration count */
        )
{
    MPI_Fint mpi_comm_f = MPI_Comm_c2f(mpi_comm);
    CALL_QDP(L);

    /* TODO alloc workspaces */
    /* all FORTRAN communication uses underscored variables */
    int ido_info_[2];
    int iparam_[11];
    int ipntr_[14];
    int n_      = loc_dim,
        nev_    = nev,
        ncv_    = ncv,
        ldv_    = loc_dim;
    double complex *resid_      = malloc(sizeof(double complex) * n_);
    memset(resid_, 0, sizeof(double complex) * n_);
    double complex *w_v_        = malloc(sizeof(double complex) * ldv_ * ncv_);
    double complex *w_workd_    = malloc(sizeof(double complex) * n_ * 3);
    memset(w_workd_, 0, sizeof(double complex) * n_ * 3);
    int lworkl_ = (3 * ncv_ * ncv_ + 5 * ncv_) * 2;
    double complex *w_workl_    = malloc(sizeof(double complex) * lworkl_);
    double *w_rwork_            = malloc(sizeof(double) *ncv_);
    /* __neupd-only workspace */
    double complex *w_d_ = malloc(sizeof(double complex) * (nev_ + 1));
    int rvec_ = 1;
    double complex sigma_ = 0;
    double complex *w_workev_ = malloc(sizeof(double complex) * 2 * ncv_);
    int *select_ = malloc(sizeof(int) * ncv_);
#define lanczosZ_free_workspace do {\
    free(resid_);\
    free(w_v_);\
    free(w_workd_);\
    free(w_workl_);\
    free(w_rwork_);\
    free(w_d_);\
    free(w_workev_);\
    free(select_);\
} while(0)

    if(NULL == resid_ ||
           NULL == w_v_ ||
           NULL == w_workd_ ||
           NULL == w_workl_ ||
           NULL == w_rwork_ ||
           NULL == w_d_ ||
           NULL == w_workev_ ||
           NULL == select_) {
        return luaL_error(L, "cannot allocate workspace for ZN*UPD");
        
    /* znaupd cycle */
    ido_        = 0;
    info_       = 0;
    iparam_[0]  = 1;
    iparam_[2]  = max_iter;
    iparam_[3]  = 1;
    iparam_[6]  = 1;
    int iter_cnt= 0;
    do {
        PARPACK_ZNAUPD(&mpi_comm_f, &ido_, "I", &n_, "SR",
                       &nev_, &tol, resid_, &ncv_, w_v_,
                       &ldv_, iparam_, ipntr_, w_workd_, w_workl_,
                       &lworkl_, w_rwork_, &info_, 1, 2);
        if (info_ < 0 || 1 < info_) {
            /* TODO set return message/status */
            lanczosZ_free_workspace;
            return luaL_error(L, "ZNAUPD returned INFO=%d", info_);
        }
        if (99 == ido_ || 1 == info_)
            break;
        if (-1 == ido_ || 1 == ido_)
            op(n_, w_workd_ + ipntr_[1] - 1, w_workd_ + ipntr_[0] - 1,
               op_arg);
        else {
            /* TODO set return message/status */
            lanczosZ_free_workspace;
            return luaL_error(L, "ZNAUPD returned IDO=%d", ido_);
        }

        iter_cnt++;
    } while (99 != ido_ && 0 == info_);
    /* for howmny="P", no additional space is required */
    PARPACK_ZNEUPD(&mpi_comm_f, &rvec_, "P", select_, w_d_,
                   w_v_, &ldv_, &sigma_, w_workev_, "I", &n_, "SR",
                   &nev_, &tol_, resid_, &ncv_, w_v_,
                   &ldv_, iparam_, ipntr_, w_workd_, w_workl_,
                   &lworkl_, w_rwork_, &info_, 1, 1, 2);

    if (0 != info_) {
        lanczosZ_free_workspace;
        return luaL_error(L, "ZNEUPD returned INFO=%d", info_);
    }
    
    /* cleanup */
    free(resid_);
    free(w_workd_);
    free(w_workl_);
    free(w_rwork_);
    free(w_workev_);
    free(select_);

    /* free extra memory in the eigenspace buffer */
    assert(ldv_ == loc_dim);
    *evec = realloc(w_v_, sizeof(double complex) * ldv_ * nev_);

    /* copy real part of eigenvalues */
    *eval = malloc(sizeof(double) * nev_);
    for (i = 0 ; i < nev_ ; i++)
        (*eval)[i] = creal(w_d_[i]);
    free(w_d_);

    *n_iters    = iter_cnt;

    return 0;
#undef lanczosZ_free_workspace
}


typedef struct {
    struct QOP_MDWF_State       *mdwf_state;
    struct QOP_MDWF_Parameters  *mdwf_param;
    struct QOP_D3_MDWF_Gauge    *mdwf_gauge;
    QOP_D3_MDWF_HalfFermion     *x, y; /* must be allocated before 
                                          calling mdwf_eoprec_op */
} op_MDWF_D3_eoprec_MdagM_arg_s;

void
op_MDWF_D3_eoprec_MdagM_op(
        int loc_dim,
        double complex *x, 
        double complex *y, 
        void *op_arg) /* x<-op(y) */
{
    op_MDWF_D3_eoprec_MdagM_arg_s *a = (op_MDWF_D3_eoprec_MdagM_arg_s *)op_arg;
    assert(2 * loc_dim == QOP_MDWF_half_fermion_size(a->mdwf_state));

    QOP_D3_MDWF_half_fermion_from_blas(a->y, (double *)y, 2 * loc_dim);
    QOP_D3_MDWF_M_operator(
            a->x, a->mdwf_param, a->mdwf_gauge, a->y);
    QOP_D3_MDWF_M_operator_conjugated(
            a->y, a->mdwf_param, a->mdwf_gauge, a->x);
    QOP_D3_MDWF_blas_from_half_fermion((double *)x, 2 * loc_dim, a->y);
}


#endif
