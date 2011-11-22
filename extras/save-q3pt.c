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

#define DEBUG_LEVEL 1   
#include "debug.h"


#define Q3PT_INDEX_STRMAX   16
const char *q3pt_index_order[] = {
    "i_t12op",      /* [ 0] */
    "vec1",         /* [ 1] */
    "vec2",         /* [ 2] */
    "spin1",        /* [ 3] */
    "spin2",        /* [ 4] */
    "i_op",         /* [ 5] */
    "i_qmom",       /* [ 6] */
    "re_im"         /* [ 7] */ };


typedef struct {
    h5output    h5o;

    int latsize[LDIM];
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



/* check all attributes 
    IN:
    OUT:
        p_x_name    name of the problem attribute
    result:
        OK:
            0
        Fail:
            -1      cannot open+read or creat+write an attribute
            -2      wrong data type of an attribute
            -3      wrong data space of an attribute
            -4      wrong data of an attribute
 
*/
/*static */int 
q3pt_check_meta(lua_State *L,
        const char **p_x_name, q3pt_h5output *q3pt_h5o)
{
    int x_status = 0;
    const char *x_name = NULL;
    {   const char *a_name = "latsize";
        hsize_t ls_dims[] = {4};
        hid_t a_space = H5Screate_simple(1, ls_dims, NULL);
        x_status = h5_check_attr_data(L, NULL, q3pt_h5o->h5o.dset, 
                    a_name, H5T_STD_I32BE, a_space, 
                    q3pt_h5o->latsize, H5T_NATIVE_INT);
        H5Sclose(a_space);
        if (x_status < 0) {
            x_name = a_name;
            goto clearerr_0;
        }
    }
    {   const char *a_name = "nvec";
        hid_t a_space = H5Screate(H5S_SCALAR);
        x_status = h5_check_attr_data(L, NULL, q3pt_h5o->h5o.dset, 
                    a_name, H5T_STD_I32BE, a_space, 
                    &(q3pt_h5o->n_vec), H5T_NATIVE_INT);
        H5Sclose(a_space);
        if (x_status < 0) {
            x_name = a_name;
            goto clearerr_0;
        }
    }
    {   const char *a_name = "qmom";
        hsize_t ls_dims[] = { q3pt_h5o->n_qmom, LDIM - 1 };
        hid_t a_space = H5Screate_simple(2, ls_dims, NULL);
        x_status = h5_check_attr_data(L, NULL, q3pt_h5o->h5o.dset, 
                    a_name, H5T_STD_I32BE, a_space, 
                    q3pt_h5o->qmom, H5T_NATIVE_INT);
        H5Sclose(a_space);
        if (x_status < 0) {
            x_name = a_name;
            goto clearerr_0;
        }
    }
    {   const char *a_name = "index_order";
        hsize_t ls_dims[] = { sizeof(q3pt_index_order) / sizeof(q3pt_index_order[0]),
                              Q3PT_INDEX_STRMAX };
        hid_t a_space = H5Screate_simple(2, ls_dims, NULL);
        char *index_order_data = h5_make_str_data(L, ls_dims, q3pt_index_order);
        x_status = h5_check_attr_data(L, NULL, q3pt_h5o->h5o.dset, 
                    a_name, H5T_STD_I8BE, a_space, 
                    index_order_data, H5T_NATIVE_CHAR);
        H5Sclose(a_space);
        if (x_status < 0) {
            x_name = a_name;
            qlua_free(L, index_order_data);
            goto clearerr_0;
        }
        qlua_free(L, index_order_data);
    }
    {   const char *a_name = "t12op";
        hsize_t ls_dims[] = { q3pt_h5o->n_t12op_max, 3 };
        hid_t a_space = H5Screate_simple(2, ls_dims, NULL);
        hid_t a_id;
        x_status = h5_check_attr(L, &a_id, q3pt_h5o->h5o.dset,
                    a_name, H5T_STD_I32BE, a_space);
        H5Sclose(a_space);
        if (x_status < 0) {
            x_name = a_name;
            goto clearerr_0;
        }
        int *buf_r = qlua_malloc(L, sizeof(int) * ls_dims[0] * ls_dims[1]);
        assert(NULL != buf_r);
        if (H5Aread(a_id, H5T_STD_I32BE, buf_r)) {
            x_status = -1;
            x_name = a_name;
            qlua_free(L, buf_r);
            H5Aclose(a_id);
            goto clearerr_0;
        }
        /* set default to invalid to indicate uninitialized entries */
        if (1 == x_status) {
            for (int k = 0; k < ls_dims[0]; k++)
                buf_r[3*k] = buf_r[3*k+1] = buf_r[3*k+2] = -1;
        } 
        /* initialize or check */
        for (int i = 0 ; i < q3pt_h5o->n_t12op; i++) {
            int k   = q3pt_h5o->i_t12op[i];
            int t_op    = q3pt_h5o->t_op[i];
            if (-1 == buf_r[3*k] && -1 == buf_r[3*k+1] && -1 == buf_r[3*k+2]) {
                buf_r[3*k]      = q3pt_h5o->t_src;
                buf_r[3*k+1]    = q3pt_h5o->t_snk;
                buf_r[3*k+2]    = t_op;
            } else if (q3pt_h5o->t_src != buf_r[3*k] 
                        || q3pt_h5o->t_snk != buf_r[3*k+1] 
                        || t_op != buf_r[3*k+2]) {
                x_status = -4;
                x_name = a_name;
                qlua_free(L, buf_r);
                H5Aclose(a_id);
                goto clearerr_0;
            }
        }
        if (H5Awrite(a_id, H5T_STD_I32BE, buf_r)) {
            x_status = -1;
            x_name = a_name;
            qlua_free(L, buf_r);
            H5Aclose(a_id);
            goto clearerr_0;
        }
        qlua_free(L, buf_r);
        H5Aclose(a_id);
    }
    return 0;

clearerr_0:
    if (NULL != p_x_name)
        *p_x_name = x_name;
    return x_status;
}

/*static*/ const char *
q3pt_h5_open_write(lua_State *L,
        q3pt_h5output *q3pt_h5o,
        const char *h5file, const char *h5path,
        const int *latsize, int n_vec,
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

    for (int d = 0 ; d < LDIM ; d++)
        q3pt_h5o->latsize[d] = latsize[d];
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
    q3pt_h5o->qmom          = qlua_malloc(L, sizeof(q3pt_h5o->qmom[0]) 
                                                * (LDIM - 1) * n_qmom);
    assert(NULL != q3pt_h5o->qmom);
    for (int i = 0; i < (LDIM - 1) * n_qmom ; i++)
        q3pt_h5o->qmom[i] = qmom[i];

    q3pt_h5o->buf_rank      = 3;
    q3pt_h5o->buf_dims[0]   = n_op;
    q3pt_h5o->buf_dims[1]   = n_qmom; 
    q3pt_h5o->buf_dims[2]   = 2;

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
/*static */const char *
save_q3pt_0deriv_selectspin(lua_State *L,
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
//    int lt = latsize[t_axis];
        
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
#if 1
    if (NULL != (err_str = q3pt_h5_open_write(L, &q3pt_h5o,
                    h5_file, h5_path, latsize, n_vec_max,
                    t_src, t_snk, n_t12op, n_t12op_max, i_t12op, t_op,
                    NSPIN*NSPIN, n_qmom, qmom_list))) {
        goto clearerr_0;
    }
#endif
    const char *x_attr;
    if (q3pt_check_meta(L, &x_attr, &q3pt_h5o)) {
        err_str = "cannot update meta info";
        luaL_error(L, err_str);
        goto clearerr_1;
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
    gsl_vector_complex_view gsl_q3pt_mom_vec = gsl_vector_complex_view_array(
            (double *)q3pt_mom_matr, n_qmom * NSPIN * NSPIN);


    for (int i_t = 0 ; i_t < n_t12op ; i_t++)
        for (int i_src = 0; i_src < n_src ; i_src++) {
            QLA_D3_DiracFermion *d_src = qla_sol_src[i_src];

            for (int i_snk = 0; i_snk < n_snk ; i_snk++) {
                QLA_D3_DiracFermion *d_snk = qla_sol_snk[i_snk];

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

                        /* replace Gamma -> (gamma_5.Gamma) in the contraction
                           because of sol_snk conjugation */
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
                    /* Fourier transform */
                    gsl_blas_zgemm(CblasNoTrans, CblasNoTrans, gsl_complex_rect(1, 0),
                            &gsl_exp_iqx.matrix, &gsl_q3pt_matr.matrix, 
                            gsl_complex_rect(0,0), &gsl_q3pt_mom.matrix);
                } else {
                    gsl_vector_complex_set_all(&gsl_q3pt_mom_vec.vector, 
                            gsl_complex_rect(0,0));
                }

                QMP_sum_double_array((double *)q3pt_mom_matr, 2 * n_qmom * NSPIN * NSPIN);

                /* multiply sol_snk with gamma5 from left: 
                   use GammaN_{ij} = A^N_i \delta_{I^N_i, j} 
                                 === A^{N(T)}_j \delta_{i,I^{N(T)}_j}
                 */
                int spin_snk_orig = spin_snk[i_snk];
                const int GTx_ind[NSPIN * NSPIN][NSPIN] = GammaT_dot_ind;
                int spin_snk_conj = GTx_ind[NSPIN * NSPIN - 1][spin_snk_orig];
                const complex double GTx_mul[NSPIN * NSPIN][NSPIN] = GammaT_dot_mul;
                gsl_complex mul_snk = complex2gsl(GTx_mul[NSPIN * NSPIN - 1][spin_snk_orig]);
                gsl_blas_zscal(mul_snk, &gsl_q3pt_mom_vec.vector);

                /* save */
                /* FIXME change order [i_qmom, iGamma] -> [iGamma, i_qmom] */
                if (NULL != (err_str = (q3pt_h5_write(&q3pt_h5o, i_t12op[i_t], 
                                    vec_src[i_src], vec_snk[i_snk], 
                                    spin_src[i_src], spin_snk_conj, 
                                    q3pt_mom_matr)))) {
                    goto clearerr_2;
                }
            }
        }

    /* cleanup */
    for (int i = 0 ; i < n_src ; i++)
        QDP_D3_reset_D(sol_src[i]);
    qlua_free(L, qla_sol_src);
    for (int i = 0 ; i < n_snk ; i++)
        QDP_D3_reset_D(sol_snk[i]);
    qlua_free(L, qla_sol_snk);
    
    qlua_free(L, q3pt_mom_matr);
    qlua_free(L, q3pt_matr);
    if (NULL != (err_str = q3pt_h5_close(L, &q3pt_h5o))) {
        goto clearerr_0;
    }
    qlua_free(L, tx2iqdp);
    qlua_free(L, exp_iqx);
    DEBUG_LEAVE;
    return NULL;


clearerr_2:
    for (int i = 0 ; i < n_src ; i++)
        QDP_D3_reset_D(sol_src[i]);
    qlua_free(L, qla_sol_src);
    for (int i = 0 ; i < n_snk ; i++)
        QDP_D3_reset_D(sol_snk[i]);
    qlua_free(L, qla_sol_snk);
    
    qlua_free(L, q3pt_mom_matr);
    qlua_free(L, q3pt_matr);

clearerr_1:
    q3pt_h5_close(L, &q3pt_h5o);

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



/* Qlua wrapper
    qcd.save_q3pt_selectspin(h5_file, h5_path,
                             max_vec, src_list, snk_list,
                             t_axis, c0, qmom_list
                             max_t12op, t12op_list, i_t12op)
    h5_file     HDF5 file name
    h5_path     keypath
    max_vec     maximum number of laph-vectors
    src_list    {{tsrc, jvec, jspin, sol}}  Note: tsrc must be the same for all
    snk_list    {{tsnk, jvec, jspin, sol}}  Note: tsnk must be the same for all
    t_axis      time axis
    c0          initial coordinates for FT
    qmom_list   insertion momenta
    max_t12op   maximal number of (t_src, t_snk, t_op) combinations
    t_op_list   list of t_op
    i_t12op     a starting number OR a list [len(t_op_list)] of indices
*/
int
q_save_q3pt_0deriv_selectspin(lua_State *L)
{
    int argc = lua_gettop(L);
    if (11 != argc) {
        luaL_error(L, "expect 11 arguments");
        return 1;
    }

    const char *h5_file     = luaL_checkstring(L, 1);
    const char *h5_path     = luaL_checkstring(L, 2);
    
    int n_src, *t_src, *vec_src, *spin_src,
        n_snk, *t_snk, *vec_snk, *spin_snk;
    QDP_D3_DiracFermion **sol_src, **sol_snk;
    mLattice *S;
    
    int max_vec             = luaL_checkint(L, 3);

    if (qlua_check_laph_sol_list(L, 4, NULL, &n_src, &t_src, &vec_src, 
                &spin_src, &sol_src, &S)) {
        luaL_error(L, "list of {tsrc, jvec, jspin, sol} objects expected in #4");
        goto clearerr_0;
    }
    for (int i = 1; i < n_src ; i++) {
        if (t_src[i] != t_src[0]) {
            luaL_error(L, "all vectors in src_list must have the same t_src");
            goto clearerr_1;
        }
    }

    if (qlua_check_laph_sol_list(L, 5, S, &n_snk, &t_snk, &vec_snk,
                &spin_snk, &sol_snk, NULL)) {
        luaL_error(L, "list of {tsrc, jvec, jspin, sol} objects expected in #5");
        goto clearerr_1;
    }
    for (int i = 1; i < n_snk ; i++) {
        if (t_snk[i] != t_snk[0]) {
            luaL_error(L, "all vectors in snk_list must have the same t_src");
            goto clearerr_2;
        }
    }
    
    int t_axis              = luaL_checkint(L, 6);

    int *c0 = qlua_check_array1d_int(L, 7, LDIM, NULL);
    if (NULL == c0) {
        luaL_error(L, "initial FT coordinate {x0,y0,z0,t0} expected in #7");
        goto clearerr_2;
    }

    int n_qmom;
    int *qmom_list = qlua_check_array2d_int(L, 8, -1, LDIM - 1, &n_qmom, NULL);
    if (NULL == qmom_list) {
        luaL_error(L, "list of momenta {{qx,qy,qz}...} expected in #8");
        goto clearerr_3;
    }

    int max_t12op           = luaL_checkint(L, 9);

    int n_t12op;
    int *t_op = qlua_check_array1d_int(L, 10, -1, &n_t12op);
    if (NULL == t_op) {
        luaL_error(L, "list of {t_op,...} expected in #10");
        goto clearerr_4;
    }

    int *i_t12op = NULL;
    if (LUA_TTABLE == lua_type(L, 11)) { 
        /* initialize i_t12op from a list */
        i_t12op = qlua_check_array1d_int(L, 11, n_t12op, NULL);
        if (NULL == i_t12op) {
            luaL_error(L, "list {i_t12op,...} must match list {t_op,...} in #11");
            goto clearerr_5;
        }
        for (int i = 0; i < n_t12op ; i++) {
            if (max_t12op <= i_t12op[i]) {
                luaL_error(L, "i_t12op=%d exceeds max_t12op=%d in #11[%d]", 
                           i_t12op[i], max_t12op, i);
                goto clearerr_6;
            }
        }
    } else if (LUA_TNUMBER == lua_type(L, 11)) {    
        /* initialize i_t12op[] with [x, x+n) */
        int i_t12op_0 = luaL_checkint(L, 11);
        if (i_t12op_0 < 0 || max_t12op < i_t12op_0 + n_t12op) {
            luaL_error(L, "incorrect range [%d, %d) for i_t12op", 
                    i_t12op_0, i_t12op + n_t12op);
            goto clearerr_5;
        }
        i_t12op = qlua_malloc(L, sizeof(i_t12op[0]) * n_t12op);
        assert(NULL != i_t12op);
        for (int i = 0 ; i < n_t12op ; i++)
            i_t12op[i] = i_t12op_0 + i;
    }

#if 1
    /* print all parameters */
    DEBUG_L1("H5='%s'['%s']\n", h5_file, h5_path);
    DEBUG_L1("max_vec=%d  max_t12op=%d  t_axis=%d  c0={%d,%d,%d,%d}\n",
            max_vec, max_t12op, t_axis, c0[0], c0[1], c0[2], c0[3]);

    DEBUG_L1("qmom[%d] = { ", n_qmom);
    for (int i = 0; i < n_qmom ; i++)
        DEBUG_L1("{%d,%d,%d} ", qmom_list[(LDIM-1)*i], qmom_list[(LDIM-1)*i+1], 
                qmom_list[(LDIM-1)*i+2]);
    DEBUG_L1("}\n");

    DEBUG_L1("src[%d] = { ", n_src);
    for (int i = 0 ; i < n_src ; i++)
        DEBUG_L1("{%d,%d,%d,%p} ", t_src[i], vec_src[i], spin_src[i], sol_src[i]);
    DEBUG_L1("}\n");

    DEBUG_L1("snk[%d] = { ", n_snk);
    for (int i = 0 ; i < n_snk ; i++)
        DEBUG_L1("{%d,%d,%d,%p} ", t_snk[i], vec_snk[i], spin_snk[i], sol_snk[i]);
    DEBUG_L1("}\n");

    DEBUG_L1("t12op[%d] = { ", n_t12op);
    for (int i = 0 ; i < n_t12op ; i++)
        DEBUG_L1("{%d->%d} ", i_t12op[i], t_op[i]);
    DEBUG_L1("}\n");
#endif
    const char *status = save_q3pt_0deriv_selectspin(L, S, h5_file, h5_path,
                            max_vec, 
                            t_src[0], n_src, vec_src, spin_src, sol_src,
                            t_snk[0], n_snk, vec_snk, spin_snk, sol_snk,
                            n_qmom, qmom_list, t_axis, c0,
                            n_t12op, max_t12op, i_t12op, t_op);
    if (NULL != status) {
        luaL_error(L, status);
        goto clearerr_6;
    }


    /* cleanup */
    qlua_free(L, i_t12op);
    qlua_free(L, t_op);
    qlua_free(L, qmom_list);
    qlua_free(L, c0);

    qlua_free(L, t_snk);
    qlua_free(L, vec_snk);
    qlua_free(L, spin_snk);
    qlua_free(L, sol_snk);

    qlua_free(L, t_src);
    qlua_free(L, vec_src);
    qlua_free(L, spin_src);
    qlua_free(L, sol_src);
    
    DEBUG_LEAVE;
    return 0;    


clearerr_6:     qlua_free(L, i_t12op);
clearerr_5:     qlua_free(L, t_op);
clearerr_4:     qlua_free(L, qmom_list);
clearerr_3:     qlua_free(L, c0);
clearerr_2:
    qlua_free(L, t_snk);
    qlua_free(L, vec_snk);
    qlua_free(L, spin_snk);
    qlua_free(L, sol_snk);
clearerr_1:
    qlua_free(L, t_src);
    qlua_free(L, vec_src);
    qlua_free(L, spin_src);
    qlua_free(L, sol_src);
clearerr_0:     return 1;
}
