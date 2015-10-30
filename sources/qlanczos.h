#ifndef QLANCZOS_H_b0StbRGSFlkiIFcndZHb
#define QLANCZOS_H_b0StbRGSFlkiIFcndZHb

#include "mpi.h"
#include <complex.h>


/* driver to lanczos_internal_float:
   return allocated memory with e.vectors and e.values in `eval' ,`evec' */
int 
lanczos_float(
        lua_State *L,
        MPI_Comm mpi_comm,
        void (*op)(int loc_dim, 
                   float complex *x,
                   float complex *y,
                   void *op_arg),   /* x<-Op(y) */
        void *op_arg,
        const char *lanczos_which,  /* ARPACK which="{S,L}{R,I,M}" */
        int loc_dim, 
        int nev,
        int ncv,
        int max_iter,
        float tol,
        float complex *v0,      /* if not NULL, contains initial vector */
        float complex **eval,   /* return buffer for evalues,  [nev] */
        float complex **evec,   /* return buffer for evectors, [ld_evec, ncol_evec] */
        int *n_iters,           /* return the iteration count */
        int *nconv,             /* number of converged evals/evecs */
        const char *arpack_logfile/* file for ARPACK log output, if not NULL */
        );

/* driver to `lanczos_internal_float':
   simply pass buffers `eval', `evec' - for inplace deflator creation */
int lanczos_inplace_float(
        lua_State *L,
        MPI_Comm mpi_comm,
        void (*op)(int loc_dim, 
                   float complex *x,
                   float complex *y,
                   void *op_arg),   /* x<-Op(y) */
        void *op_arg,
        const char *lanczos_which,  /* ARPACK which="{S,L}{R,I,M}" */
        int loc_dim, 
        int nev,
        int ncv,
        int max_iter,
        float tol,
        float complex *v0,      /* if not NULL, contains initial vector */
        float complex *eval,    /* (allocated) return buffer for evalues,  [nev] */
        float complex *evec,    /* (allocated) return buffer for evectors, [ld_evec, ncol_evec] */
        int ld_evec,            /* BLAS leading dimension, >=loc_dim */
        int ncol_evec,          /* number of columns in evec, >=ncv */
        int *n_iters,           /* return the iteration count */
        int *nconv,             /* number of converged evals/evecs */
        const char *arpack_logfile/* file for ARPACK log output, if not NULL */
        );

int lanczos_internal_float(
        lua_State *L,
        MPI_Comm mpi_comm,
        void (*op)(int loc_dim, 
                   float complex *x,
                   float complex *y,
                   void *op_arg),   /* x<-Op(y) */
        void *op_arg,
        const char *lanczos_which,  /* ARPACK which="{S,L}{R,I,M}" */
        int loc_dim, 
        int nev,
        int ncv,
        int max_iter,
        float tol,
        float complex *v0,      /* if not NULL, contains initial vector */
        float complex *eval,    /* return buffer for evalues,  [nev + 1] */
        float complex *evec,    /* workspace for evectors (BLAS-matrix), [ld_evec, >=ncv] */
        int ld_evec,            /* BLAS leading dimension */
        int *n_iters,           /* return the iteration count */
        int *nconv,             /* number of converged evals/evecs */
        const char *arpack_logfile/* file for ARPACK log output, if not NULL */
        );

#endif/*QLANCZOS_H_b0StbRGSFlkiIFcndZHb*/
