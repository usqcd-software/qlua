#include <math.h>
#include <string.h>
#include <assert.h>

#include <hdf5.h>

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

#define NSPIN   4
#define LDIM    4

typedef struct {
    h5output    h5o;

    int lt;
    int n_vec;

    int t_src;
    int j_vec_prop;
    int j_spin_prop;

    int         buf_rank;
    hsize_t     buf_dims[4];
    hid_t       buf_dspace;
} q2pt_h5output;

/*static*/ const char *
q2pt_h5_open_write(lua_State *L,
            q2pt_h5output *q2pt_h5o,
            const char *h5file, const char *h5path,
            int lt, int n_vec, 
            int t_src, int j_vec_prop, int j_spin_prop)
{
    if (0 != QDP_this_node)
        return NULL;

    const char *err_str = NULL;

    hsize_t dims[7] = { lt, lt, n_vec, n_vec, NSPIN, NSPIN, 2};
    err_str = h5_open_write(L, &(q2pt_h5o->h5o), h5file, h5path, 7, dims);
    if (NULL != err_str)
        return err_str;

    q2pt_h5o->lt            = lt;
    q2pt_h5o->n_vec         = n_vec;
    q2pt_h5o->t_src         = t_src;
    q2pt_h5o->j_vec_prop    = j_vec_prop;
    q2pt_h5o->j_spin_prop   = j_spin_prop;

    q2pt_h5o->buf_rank = 4;
    q2pt_h5o->buf_dims[0]   = lt;
    q2pt_h5o->buf_dims[1]   = n_vec;
    q2pt_h5o->buf_dims[2]   = NSPIN;
    q2pt_h5o->buf_dims[3]   = 2;

    if (0 > (q2pt_h5o->buf_dspace = 
                H5Screate_simple(q2pt_h5o->buf_rank, q2pt_h5o->buf_dims, NULL)))
        return "cannot create memory dataspace";

    return NULL;
}


/* close q2pt_h5 structures
   XXX h5o is NOT deallocated
 */
/*static*/ const char *
q2pt_h5_close(lua_State *L, q2pt_h5output *q2pt_h5o)
{
    if (0 != QDP_this_node)
        return NULL;
    if (H5Sclose(q2pt_h5o->buf_dspace))
        return "cannot close HDF5 membuf";
    q2pt_h5o->buf_dspace = -1;

    return h5_close(L, &(q2pt_h5o->h5o));
}

/* write a portion of data 
   buf      complex array [t1][n1][s1] of 
            q2pt[t1][t2][n1][n2][s1][s2][real/imag?0:1]
 */
/*static*/ const char *
q2pt_h5_write(q2pt_h5output *q2pt_h5o, 
              const gsl_complex *buf)
{
    if (0 != QDP_this_node)
        return NULL;
    hsize_t dset_off[7] = { 0, 
                            q2pt_h5o->t_src, 
                            0,
                            q2pt_h5o->j_vec_prop,
                            0,
                            q2pt_h5o->j_spin_prop,
                            0 };
    hsize_t dset_cnt[7] = { q2pt_h5o->lt,
                            1,
                            q2pt_h5o->n_vec,
                            1,
                            NSPIN,
                            1,
                            2 };
    if (H5Sselect_hyperslab(q2pt_h5o->h5o.dspace, H5S_SELECT_SET, 
                dset_off, NULL, dset_cnt, NULL))
        return "cannot select hyperslab";
    assert(H5Sget_select_npoints(q2pt_h5o->h5o.dspace) == 
            H5Sget_select_npoints(q2pt_h5o->buf_dspace));
    if (H5Dwrite(q2pt_h5o->h5o.dset, H5T_NATIVE_DOUBLE, 
                 q2pt_h5o->buf_dspace, q2pt_h5o->h5o.dspace,
                 H5P_DEFAULT, buf))
        return "cannot write to dataset";
    return NULL;
}


/* main:
   cycle through propagator(s) components and save its projection(s) on 
   LapH eigenmodes for each timeslice
 */
/*static*/ const char *
save_q2pt(lua_State *L, 
        mLattice *S,
        const char *h5_file,
        const char *h5_path,
        int t_src, int j_vec_prop, int j_spin_prop, QDP_D3_DiracFermion *sol,
        int n_vec, QDP_D3_ColorVector **v,
        int t_axis)
{
    const char *err_str = NULL;

    if (LDIM != S->rank
            || t_axis < 0 || LDIM <= t_axis) {
        /* FIXME test for other t-axis and remove the check */
        return "not implemented for this dim or t-axis";
    }
    int latsize[LDIM];
    QDP_latsize_L(S->lat, latsize);
    int lt = latsize[t_axis];

    int vol4_local = QDP_sites_on_node_L(S->lat);

    /* HDF5 output */
    q2pt_h5output q2pt_h5o;
    if (NULL != (err_str = q2pt_h5_open_write(L, &q2pt_h5o, 
                    h5_file, h5_path, lt, n_vec, 
                    t_src, j_vec_prop, j_spin_prop))) {
        goto clearerr_0;
    }

    /* expose fields */
    QLA_D3_ColorVector **qla_v = NULL;
    qla_v = qlua_malloc(L, sizeof(qla_v[0]) * n_vec);
    for (int i = 0 ; i < n_vec ; i++)
        qla_v[i] = QDP_D3_expose_V(v[i]);

    QLA_D3_DiracFermion *qla_sol;
    qla_sol = QDP_D3_expose_D(sol);

    /* buffer for q2pt */
    gsl_complex *q2pt_buf = NULL;
    q2pt_buf = qlua_malloc(L, sizeof(q2pt_buf[0]) * lt * n_vec * NSPIN);

#if LDIM == 4  
    assert(LDIM == S->rank);
#define q2pt(t, n, s) ((q2pt_buf)[(s) + (NSPIN) * ((n) + (n_vec) * (t))])
#else
#error "LDIM!=4 is not consistent"
#endif

    memset(q2pt_buf, 0, sizeof(q2pt_buf[0]) * lt * n_vec * NSPIN);
    for (int i_qdp = 0 ; i_qdp < vol4_local ; i_qdp++) {
        int coord[LDIM];
        QDP_get_coords_L(S->lat, coord, QDP_this_node, i_qdp);
        int t = coord[t_axis];

        for (int i_v = 0 ; i_v < n_vec ; i_v++) {
            for (int s = 0 ; s < NSPIN ; s++) {
                /* calculate product */
                QLA_Complex res;
                QLA_c_eq_c_times_ca(res, QLA_elem_D(qla_sol[i_qdp], 0, s), 
                        QLA_elem_V(qla_v[i_v][i_qdp], 0));
                QLA_c_peq_c_times_ca(res, QLA_elem_D(qla_sol[i_qdp], 1, s), 
                        QLA_elem_V(qla_v[i_v][i_qdp], 1));
                QLA_c_peq_c_times_ca(res, QLA_elem_D(qla_sol[i_qdp], 2, s),
                        QLA_elem_V(qla_v[i_v][i_qdp], 2));
                q2pt(t, i_v, s) = gsl_complex_add(q2pt(t, i_v, s),
                                            gsl_complex_rect(QLA_real(res), QLA_imag(res)));
            }
        }
    }
    QMP_sum_double_array((double *)q2pt_buf, 2 * lt * n_vec * NSPIN);
    if (NULL != (err_str = q2pt_h5_write(&q2pt_h5o, q2pt_buf)))
        goto clearerr_1;
    if (NULL != (err_str = q2pt_h5_close(L, &(q2pt_h5o))))
        goto clearerr_1;

clearerr_1:
    /* cleanup */
    for (int i = 0 ; i < n_vec ; i++)
        QDP_D3_reset_V(v[i]);
    qlua_free(L, qla_v);
    QDP_D3_reset_D(sol);
    
    qlua_free(L, q2pt_buf);

clearerr_0:

    return err_str;
#undef q2pt
}

/* QLua wrapper 
   qcd.save_squark_wf(h5_file, h5_path, 
                      t_src, j_vec, j_spin, sol,
                      vec_list, t_axis)
    h5_file     HDF5 file name
    h5_path     keypath
    t_src, j_vec, j_spin    
                source time, evec and spin
    vec_list    LapH vectors
    t_axis

 */
int         
q_save_q2pt(lua_State *L)
{
    int argc = lua_gettop(L);
    if (8 != argc) {
        luaL_error(L, "expect 8 arguments");
        return 1;
    }

    const char *h5_file     = luaL_checkstring(L, 1);
    const char *h5_path     = luaL_checkstring(L, 2);

    int t_src               = luaL_checkint(L, 3);
    int j_vec               = luaL_checkint(L, 4);
    int j_spin              = luaL_checkint(L, 5);

    mLattice *S             = qlua_ObjLattice(L, 6);
    QDP_D3_DiracFermion *sol= qlua_checkLatDirFerm3(L, 6, S, 3)->ptr;
    int n_vec;
    QDP_D3_ColorVector **v  = qlua_check_latcolvec_table(L, 7, S, -1, NULL, &n_vec);
    int t_axis              = luaL_checkint(L, 8);

    const char *status = save_q2pt(L, S, h5_file, h5_path,
                            t_src, j_vec, j_spin, sol,
                            n_vec, v,
                            t_axis);

    qlua_free(L, v);
    
    if (NULL != status)
        luaL_error(L, status);

    return 0;
}

#undef LDIM
#undef NSPIN
