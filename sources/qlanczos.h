#ifndef QLANCZOS_H_b0StbRGSFlkiIFcndZHb
#define QLANCZOS_H_b0StbRGSFlkiIFcndZHb

#include "mpi.h"
#include <complex.h>
int
lanczos_internal_float(
        lua_State *L,
        MPI_Comm mpi_comm,
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
        );

#endif/*QLANCZOS_H_b0StbRGSFlkiIFcndZHb*/
