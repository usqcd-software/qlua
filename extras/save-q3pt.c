#include <math.h>
#include <string.h>
#include <assert.h>

#include <hdf5.h>

#include <complex.h>

#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_matrix.h>

#include "qmp.h"

#include "modules.h"                                                /* DEPS */
#include "qlua.h"                                                   /* DEPS */
#include "lattice.h"                                                /* DEPS */
#include "extras.h"                                                 /* DEPS */
#include "latcolvec.h"                                              /* DEPS */
#include "latdirferm.h"                                             /* DEPS */

#include "h5common.h"
#include "qparam.h"
#include "laph_common.h"

/* definitions of gamma matrices */
#include "gamma_dgr_defs.h"
#include "gamma_defs.h"

#define NSPIN   4
#define NCOLOR  3
#define LDIM    4




/* TODO make GammaX, GammaTX 
#define cxP1(z) (z)
#define cxM1(z) (gsl_complex_rect(-real, -imag))
#define cxPI(z) (gsl_complex_rect(-imag, +real))
#define cxMI(z) (gsl_complex_rect(+imag, -real))
 
 */

typedef struct {
    h5output    h5o;

    int lt;
    int n_vec;

    int t_src, t_snk;       /* source and sink location: for Attr checking */
    int n_t12op;            /* number of op.insertion tslices */
    int n_t12op_max;        /* max. number of op.insertion tslices */
    int *i_t12op;           /* [n_t12op] index in t12op map */
    int *t_op;              /* [n_t12op] op.insertion tslices */

    int n_op;

    int n_qmom;             /* number of mom insertions */
    int *qmom;              /* mom insertions [n_qmom][LDIM-1] */

    int         buf_rank;
    hsize_t     buf_dims[3];
    hid_t       buf_dspace;
} q3pt_h5output;

/* TODO get right params for i_t12op */
const char *
q3pt_h5_open_write(lua_State *L,
        q3pt_h5output *q3pt_h5o,
        const char *h5file, const char *h5path,
        int lt, int n_vec,
        int t_src, int t_snk,
        int n_t12op, int n_t12op_max, const int *i_t12op, const int *t_op,
        int n_op,
        int n_qmom, const int *qmom)
/* store q3pt correlator:
   [i_t12op, n1, n2, s1, s2, i_op, i_qmom, re/im]
    
 */
{
    if (0 != QDP_this_node)
        return NULL;

    const char *err_str = NULL;
    /* TODO put right dimensions */
    hsize_t dims[8] = { n_t12op_max, n_vec, n_vec, NSPIN, NSPIN, n_op, n_qmom, 2};
    err_str = h5_open_write(L, &(q3pt_h5o->h5o), h5file, h5path, 8, dims);
    if (NULL != err_str)
        return err_str;

    q3pt_h5o->lt            = lt;
    q3pt_h5o->n_vec         = n_vec;
    q3pt_h5o->t_src         = t_src;
    q3pt_h5o->t_snk         = t_snk;
    q3pt_h5o->n_op          = n_op;

    q3pt_h5o->n_t12op       = n_t12op;
    q3pt_h5o->n_t12op_max   = n_t12op_max;
    
    q3pt_h5o->i_t12op       = qlua_malloc(L, sizeof(q3pt_h5o->i_t12op[0]) * n_t12op);
    assert(NULL != q3pt_h5o->i_t12op);
    q3pt_h5o->t_op          = qlua_malloc(L, sizeof(q3pt_h5o->t_op[0]) * n_t12op);
    assert(NULL != q3pt_h5o->t_op);
    for (int i = 0; i < n_t12op ; i++) {
        q3pt_h5o->i_t12op[i]    = i_t12op[i];
        q3pt_h5o->t_op[i]       = t_op[i];
    }

    q3pt_h5o->n_qmom        = n_qmom;
    q3pt_h5o->qmom          = qlua_malloc(L, sizeof(q3pt_h5o->qmom[0]) * n_qmom);
    assert(NULL != q3pt_h5o->qmom);
    for (int i = 0; i < (LDIM - 1) * n_qmom ; i++)
        q3pt_h5o->qmom[i] = qmom[i];

    q3pt_h5o->buf_rank      = 3;
    q3pt_h5o->buf_dims[1]   = n_op;
    q3pt_h5o->buf_dims[2]   = n_qmom; 
    q3pt_h5o->buf_dims[3]   = 2;

    if (0 > (q3pt_h5o->buf_dspace =
                H5Screate_simple(q3pt_h5o->buf_rank, q3pt_h5o->buf_dims, NULL))) {
        err_str = "cannot create memory dataspace";
        goto clearerr_0;
    }

    return NULL;

clearerr_0:
    qlua_free(L, q3pt_h5o->qmom);
    qlua_free(L, q3pt_h5o->t_op);
    qlua_free(L, q3pt_h5o->i_t12op);

    return err_str;
}

/*static */ const char *
q3pt_h5_close(lua_State *L, q3pt_h5output *q3pt_h5o)
{
    if (0 != QDP_this_node)
        return NULL;
    if (NULL != q3pt_h5o->i_t12op) 
        qlua_free(L, q3pt_h5o->i_t12op);
    if (NULL != q3pt_h5o->t_op)
        qlua_free(L, q3pt_h5o->t_op);
    if (NULL != q3pt_h5o->qmom)
        qlua_free(L, q3pt_h5o->qmom);

    if (H5Sclose(q3pt_h5o->buf_dspace))
        return "cannot close HDF5 membuf";
    
    return h5_close(L, &(q3pt_h5o->h5o));
}

/*static*/ const char *
q3pt_h5_write(q3pt_h5output *q3pt_h5o, int i_t12op,
        int i_vec_1, int i_vec_2, int s1, int s2,
        const gsl_complex *buf)
/* save [n_op, n_qmom, re/im] 
   datafile: [i_t12op, n1, n2, s1, s2, i_op, i_qmom, re/im]
 */
{
    hsize_t dset_off[8] = { q3pt_h5o->i_t12op[i_t12op],
                            i_vec_1,
                            i_vec_2,
                            s1,
                            s2,
                            0,
                            0,
                            0 };
    hsize_t dset_cnt[8] = { 1,
                            1,
                            1,
                            1,
                            1,
                            q3pt_h5o->n_op,
                            q3pt_h5o->n_qmom,
                            2 };
    if (H5Sselect_hyperslab(q3pt_h5o->h5o.dspace, H5S_SELECT_SET,
                    dset_off, NULL, dset_cnt, NULL))
        return "cannot select hyperslab";
    assert(H5Sget_select_npoints(q3pt_h5o->h5o.dspace) ==
            H5Sget_select_npoints(q3pt_h5o->buf_dspace));
    if (H5Dwrite(q3pt_h5o->h5o.dset, H5T_NATIVE_DOUBLE, 
                    q3pt_h5o->buf_dspace, q3pt_h5o->h5o.dspace,
                    H5P_DEFAULT, buf))
        return "cannot write to dataset";
    return NULL;
}

/* TODO get right params for i_t12op */
const char *
save_q3pt_0deriv(lua_State *L,
        mLattice *S,
        const char *h5_file,
        const char *h5_path,
        int n_vec_max,
        int t_src, int n_src, const int *vec_src, const int *spin_src, 
        QDP_D3_DiracFermion **sol_src,
        int t_snk, int n_snk, const int *vec_snk, const int *spin_snk,
        QDP_D3_DiracFermion **sol_snk,
        int n_qmom, const int *qmom_list,   /* [n_qmom][LDIM-1] */
        int t_axis, const int *c0,
        int n_t12op, int n_t12op_max, const int *i_t12op, const int *t_op)
/* calc correlator C3(tsrc,tsnk;t) 
        = < q^{vec_snk, spin_snk}(tsnk) 
            \sum_z3 (e^{i*q.z3}[q Gamma_G \bar q]_{z3,t}
            \bar{q}^{vec_src, spin_src}(tsrc) >

  t_op      [n_t12op] timeslice: operator location (abs)
  {}
  IN:
  sol_src[Jsrc={vec_src, spin_src}][iqdp][c,jspin]
  sol_snk[Jsnk={vec_snk, spin_snk}][iqdp][c,ispin]
  OUT:
  cycle over Jsrc, Jsnk, t [G]
 */
{
    const char *err_str = NULL;

    if (LDIM != S->rank
            || t_axis < 0 || LDIM <= t_axis)
        return "not implemented for this dim or t-axis";
        
    int latsize[LDIM];
    QDP_latsize_L(S->lat, latsize);
    int lt = latsize[t_axis];
        
    /* determine subgrid parameters */
    int node_latsize[LDIM];
    int node_x0[LDIM];                      
    calc_subgrid(S->lat, node_latsize, node_x0);
        
    int vol3_dim[LDIM-1];
    int vol3_axis[LDIM-1];
    int node_lt = node_latsize[t_axis];
    int node_t0 = node_x0[t_axis];
    int vol3_local  = 1,
        vol4_local  = 1;
    for (int i=0, k=0 ; i < LDIM ; i++) {
        vol4_local  *= node_latsize[i];
        if (i == t_axis)
            continue;
        vol3_axis[k]= i;
        vol3_dim[k] = node_latsize[i];
        k++;
        vol3_local  *= node_latsize[i];
    }       
    assert(QDP_sites_on_node_L(S->lat) == vol4_local);
#if LDIM == 4  
    assert(LDIM == S->rank);
#define i_vol3(c) ((c)[vol3_axis[0]] + vol3_dim[0]*(\
                   (c)[vol3_axis[1]] + vol3_dim[1]*(\
                   (c)[vol3_axis[2]])))
#define i_X_vol3(X,c) (i_vol3(c) + vol3_local*(X))
#define i_vol4(c) i_X_vol3(((c)[t_axis] - node_t0), c)
#define mom(i_qmom, k)   ((qmom_list)[(i_qmom)*(LDIM-1) + k])
#else
#error "LDIM!=4 is not consistent"
#endif

    /* make map i_qdp -> tsnk, i_x3 */
    int *tx2iqdp = NULL;
    tx2iqdp     = qlua_malloc(L, sizeof(tx2iqdp[0]) * vol4_local);
    assert(NULL != tx2iqdp);
#define Tx2iqdp(t,ix) ((tx2iqdp)[ix + (vol3_local)*((t) - (node_t0))])
#define have_timeslice(t)   (( (node_t0) <= (t) ) && ( (t) < (node_t0) + (node_lt) ))
    for (int iqdp = 0; iqdp < vol4_local ; iqdp++) {
        int coord[LDIM];
        QDP_get_coords_L(S->lat, coord, QDP_this_node, iqdp);
        /* FIXME change order of ix3 for fewer cache misses 
            when converting data [iqdp] <-> [t, i_x3] ? */
        int t = coord[t_axis],
            ix = i_vol3(coord);
        Tx2iqdp(t, ix) = iqdp;
    }
    
    /* exp(i\vec{q}.\vec{x} [i_qmom, i_x3] */
    gsl_complex *exp_iqx = NULL;
    exp_iqx = qlua_malloc(L, sizeof(exp_iqx[0]) * n_qmom * vol3_local);
    assert(NULL != exp_iqx);
    gsl_matrix_complex_view  gsl_exp_iqx = gsl_matrix_complex_view_array(
            (double *)exp_iqx, n_qmom, vol3_local);

    for (int i_qdp = 0 ; i_qdp < vol4_local ; i_qdp++) {
        int coord[LDIM];
        QDP_get_coords_L(S->lat, coord, QDP_this_node, i_qdp);
        if (coord[t_axis] != node_t0)
            continue;
        for (int i_qmom = 0; i_qmom < n_qmom ; i_qmom++) {
            double phase = 0.;
            for (int k = 0; k < LDIM - 1; k++)
                phase += mom(i_qmom, k) * (coord[vol3_axis[k]] - c0[vol3_axis[k]])
                            / (double)latsize[k];
            exp_iqx[i_X_vol3(i_qmom, coord)] = gsl_complex_polar(1., phase*2*M_PI);
        }
    }

    /* HDF5 output */
    q3pt_h5output q3pt_h5o;
    if (NULL != (err_str = q3pt_h5_open_write(L, &q3pt_h5o,
                    h5_file, h5_path, lt, n_vec_max,
                    t_src, t_snk, n_t12op, n_t12op_max, i_t12op, t_op,
                    NSPIN*NSPIN, n_qmom, qmom_list))) {
        goto clearerr_0;
    }

    /* expose src */
    QLA_D3_DiracFermion **qla_sol_src = NULL;
    qla_sol_src = qlua_malloc(L, sizeof(qla_sol_src[0]) *n_src);
    assert(NULL != qla_sol_src);
    for (int i = 0 ; i < n_src ; i++)
        qla_sol_src[i] = QDP_D3_expose_D(sol_src[i]);
    
    /* expose snk */
    QLA_D3_DiracFermion **qla_sol_snk = NULL;
    qla_sol_snk = qlua_malloc(L, sizeof(qla_sol_snk[0]) *n_snk);
    assert(NULL != qla_sol_snk);
    for (int i = 0 ; i < n_snk ; i++)
        qla_sol_snk[i] = QDP_D3_expose_D(sol_snk[i]);

    /* q3pt contraction [i_x3, iGamma] */
    gsl_complex *q3pt_matr = NULL;
    q3pt_matr = qlua_malloc(L, sizeof(q3pt_matr[0]) * vol3_local * NSPIN * NSPIN);
    assert(NULL != q3pt_matr);
    gsl_matrix_complex_view gsl_q3pt_matr = gsl_matrix_complex_view_array(
            (double *)q3pt_matr, vol3_local, NSPIN * NSPIN);
#define Q3PTmatr(ix3, iG) (q3pt_matr[(iG) + (NSPIN)*(NSPIN)*(ix3)])

    /* q3pt_qmom [i_qmom, iGamma] */
    gsl_complex *q3pt_mom_matr = NULL;
    q3pt_mom_matr = qlua_malloc(L, sizeof(q3pt_mom_matr[0]) * n_qmom * NSPIN * NSPIN);
    assert(NULL != q3pt_mom_matr);
    gsl_matrix_complex_view gsl_q3pt_mom = gsl_matrix_complex_view_array(
            (double *)q3pt_mom_matr, n_qmom, NSPIN * NSPIN);

//    const complex double GTx_mul[NSPIN * NSPIN][NSPIN] = GammaT_dot_mul;
//    const int GTx_ind[NSPIN * NSPIN][NSPIN] = GammaT_dot_ind;

    for (int i_t = 0 ; i_t < n_t12op ; i_t++)
        for (int i_src = 0; i_src < n_src ; i_src++) {
//            int v_src = vec_src[i_src], 
//                s_src = spin_src[i_src];
            QLA_D3_DiracFermion *d_src = qla_sol_src[i_src];

            for (int i_snk = 0; i_snk < n_snk ; i_snk++) {
                QLA_D3_DiracFermion *d_snk = qla_sol_snk[i_snk];
//                int v_snk = vec_snk[i_snk], s_snk = spin_snk[i_snk];
//                complex snk_GT15_mul = GTx_mul[15][s_snk];
//                int     snk_GT15_idx = GTx_ind[15][s_snk];

                if (have_timeslice(t_op[i_t])) {
                    for (int i_x3 = 0 ; i_x3 < vol3_local ; i_x3++) {
                        int i_qdp = Tx2iqdp(t_op[i_t], i_x3);
                        complex double z[NSPIN][NSPIN];
                        for (int is = 0; is < NSPIN ; is++) {
                            for (int js = 0; js < NSPIN ; js++) {
                                complex double x = 0;
                                for (int kc = 0 ; kc < NCOLOR ; kc++)
                                    x   +=  qla2complex(QLA_elem_D(d_src[i_qdp], kc, is)) *
                                            conj(qla2complex(QLA_elem_D(d_snk[i_qdp], kc, js)));
                                z[is][js] = x;
                            }
                        }
#define Q3tr(iG)  TrGamma_dot(iG, \
                z[0][0], z[0][1], z[0][2], z[0][3], \
                z[1][0], z[1][1], z[1][2], z[1][3], \
                z[2][0], z[2][1], z[2][2], z[2][3], \
                z[3][0], z[3][1], z[3][2], z[3][3])

                        /* replace Gamma with (gamma_5.Gamma) for sol_snk conjugation */
                        Q3PTmatr(i_x3, 0) = complex2gsl( Q3tr(15));
                        Q3PTmatr(i_x3, 1) = complex2gsl(-Q3tr(14));
                        Q3PTmatr(i_x3, 2) = complex2gsl( Q3tr(13));
                        Q3PTmatr(i_x3, 3) = complex2gsl(-Q3tr(12));
                                                                  
                        Q3PTmatr(i_x3, 4) = complex2gsl(-Q3tr(11));
                        Q3PTmatr(i_x3, 5) = complex2gsl( Q3tr(10));
                        Q3PTmatr(i_x3, 6) = complex2gsl(-Q3tr( 9));
                        Q3PTmatr(i_x3, 7) = complex2gsl( Q3tr( 8));
                                                                  
                        Q3PTmatr(i_x3, 8) = complex2gsl( Q3tr( 7));
                        Q3PTmatr(i_x3, 9) = complex2gsl(-Q3tr( 6));
                        Q3PTmatr(i_x3,10) = complex2gsl( Q3tr( 5));
                        Q3PTmatr(i_x3,11) = complex2gsl(-Q3tr( 4));
                                                                  
                        Q3PTmatr(i_x3,12) = complex2gsl(-Q3tr( 3));
                        Q3PTmatr(i_x3,13) = complex2gsl( Q3tr( 2));
                        Q3PTmatr(i_x3,14) = complex2gsl(-Q3tr( 1));
                        Q3PTmatr(i_x3,15) = complex2gsl( Q3tr( 0));
#undef Q3tr
                    }
                    /* Fourier transf */
                    gsl_blas_zgemm(CblasNoTrans, CblasNoTrans, gsl_complex_rect(1, 0),  
                            &gsl_exp_iqx.matrix, &gsl_q3pt_matr.matrix, 
                            gsl_complex_rect(0,0), &gsl_q3pt_mom.matrix);
                } else {
                    memset((void *)(q3pt_mom_matr), 0, 
                           sizeof(q3pt_mom_matr[0]) * n_qmom * NSPIN * NSPIN);
                }

                QMP_sum_double_array((double *)q3pt_mom_matr, 2 * n_qmom * NSPIN * NSPIN);

                /* save */
                /* FIXME multiply sol_snk with gamma5 from left */
                if (NULL != (err_str = (q3pt_h5_write(&q3pt_h5o, i_t12op[i_t], 
                                    vec_src[i_src], vec_snk[i_snk], 
                                    spin_src[i_src], spin_snk[i_snk], 
                                    q3pt_mom_matr)))) {
                    goto clearerr_2;
                }
            }
        }

clearerr_2:
    if (NULL != (err_str = q3pt_h5_close(L, &q3pt_h5o))) {
        goto clearerr_1;
    }

    /* cleanup */               
clearerr_1:

    for (int i = 0 ; i < n_src ; i++)
        QDP_D3_reset_D(sol_src[i]);
    qlua_free(L, qla_sol_src);
    for (int i = 0 ; i < n_snk ; i++)
        QDP_D3_reset_D(sol_snk[i]);
    qlua_free(L, qla_sol_snk);
    
    qlua_free(L, q3pt_mom_matr);
    qlua_free(L, q3pt_matr);

clearerr_0:
    qlua_free(L, tx2iqdp);
    qlua_free(L, exp_iqx);

    return err_str;

#undef i_vol3
#undef Q3tr
#undef Q3PTmatr
#undef Tx2iqdp
#undef have_timeslice
#undef mom
#undef i_vol4
#undef i_X_vol3
}

