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

#include "h5common.h"
#include "qparam.h"
#include "laph_common.h"

#define V123_INDEX_STRMAX   16
const char *v123_index_order[] = {
    "vec1",
    "vec2",
    "vec3",
    "t",
    "i_mom",
    "re_im"
};

static void
cross_VV(QLA_D3_ColorVector *r, 
         QLA_D3_ColorVector *c1, QLA_D3_ColorVector *c2)
{
    QLA_c_eq_c_times_c(QLA_D3_elem_V(*r, 0), QLA_D3_elem_V(*c1, 1), QLA_D3_elem_V(*c2, 2));
    QLA_c_meq_c_times_c(QLA_D3_elem_V(*r, 0), QLA_D3_elem_V(*c1, 2), QLA_D3_elem_V(*c2, 1));

    QLA_c_eq_c_times_c(QLA_D3_elem_V(*r, 1), QLA_D3_elem_V(*c1, 2), QLA_D3_elem_V(*c2, 0));
    QLA_c_meq_c_times_c(QLA_D3_elem_V(*r, 1), QLA_D3_elem_V(*c1, 0), QLA_D3_elem_V(*c2, 2));

    QLA_c_eq_c_times_c(QLA_D3_elem_V(*r, 2), QLA_D3_elem_V(*c1, 0), QLA_D3_elem_V(*c2, 1));
    QLA_c_meq_c_times_c(QLA_D3_elem_V(*r, 2), QLA_D3_elem_V(*c1, 1), QLA_D3_elem_V(*c2, 0));
}

static gsl_complex
cmult_VV(QLA_D3_ColorVector c1, QLA_D3_ColorVector c2)
{
    QLA_D_Complex c;
    QLA_c_eq_c_times_c(c, QLA_D3_elem_V(c1, 0), QLA_D3_elem_V(c2, 0));
    QLA_c_peq_c_times_c(c, QLA_D3_elem_V(c1, 1), QLA_D3_elem_V(c2, 1));
    QLA_c_peq_c_times_c(c, QLA_D3_elem_V(c1, 2), QLA_D3_elem_V(c2, 2));
    return gsl_complex_rect(QLA_real(c), QLA_imag(c));
}

typedef struct {
    h5output    h5o;

    int         latsize[LDIM];
    int         t_axis;

    int         n_mom;
    int         *ft_mom;
    int         ft_x0[LDIM];

    int         n_v1, n_v2, n_v3;
    
    int         buf_rank;
    hsize_t     buf_dims[3];
    hid_t       buf_dspace;
} barpw_h5output;


static int
bpw_h5_check_meta(lua_State *L,
        const char **p_x_name, barpw_h5output *bpw_h5o)
{
    if (! is_masternode())
        return 0;

    hid_t d_id = bpw_h5o->h5o.dset;
    int x_status = 0;
    const char *x_name = NULL;
    int n_vec_a[] = { bpw_h5o->n_v1, bpw_h5o->n_v2, bpw_h5o->n_v3 };

    if ((x_status = h5_check_attr_str_list(L, NULL, d_id, "index_order", 
                        sizeof(v123_index_order) / sizeof(v123_index_order[0]),
                        V123_INDEX_STRMAX, v123_index_order)) < 0) {
        x_name = "index_order";
        goto clearerr_0;
    }
    if ((x_status = h5_check_attr_array1d_int(L, NULL, d_id,
                            "latsize", LDIM, bpw_h5o->latsize)) < 0) {
        x_name = "latsize";
        goto clearerr_0;
    }
    if ((x_status = h5_check_attr_int(L, NULL, d_id,
                        "t_axis", &bpw_h5o->t_axis)) < 0) {
        x_name = "t_axis";
        goto clearerr_0;
    }
    
    if ((x_status = h5_check_attr_array1d_int(L, NULL, d_id,
                        "nvec123", 3, n_vec_a)) < 0) {
        x_name = "nvec123";
        goto clearerr_0;
    }
    if ((x_status = h5_check_attr_array2d_int(L, NULL, d_id,
                        "ft_mom", bpw_h5o->n_mom, LDIM - 1, bpw_h5o->ft_mom)) < 0) {
        x_name = "ft_mom";
        goto clearerr_0;
    }
    if ((x_status = h5_check_attr_array1d_int(L, NULL, d_id,
                            "ft_x0", LDIM, bpw_h5o->ft_x0)) < 0) {
        x_name = "ft_x0";
        goto clearerr_0;
    }

    return 0;

clearerr_0:

    if (NULL != p_x_name)
        *p_x_name = x_name;
    return x_status;
}

static const char *
bpw_h5_open_write(lua_State *L, 
        barpw_h5output *bpw_h5o, 
        const char *h5file, const char *h5path,
        const int latsize[], int t_axis, 
        int n_mom, const int ft_mom[], const int ft_x0[],
        int n_v1, int n_v2, int n_v3)
{
    if (! is_masternode())
        return NULL;

    const char *err_str = NULL;
    int lt = latsize[t_axis];
    hsize_t dims[6] = {n_v1, n_v2, n_v3, lt, n_mom, 2};
    err_str = h5_open_write(L, &(bpw_h5o->h5o), h5file, h5path, 6, dims);
    if (NULL != err_str)
        return err_str;

    for (int d = 0 ; d < LDIM ; d++)
        bpw_h5o->latsize[d] = latsize[d];
    bpw_h5o->t_axis = t_axis;
    bpw_h5o->n_v1   = n_v1;
    bpw_h5o->n_v2   = n_v2;
    bpw_h5o->n_v3   = n_v3;

    bpw_h5o->n_mom  = n_mom;
    bpw_h5o->ft_mom = qlua_malloc(L, sizeof(bpw_h5o->ft_mom[0])
                                * (LDIM - 1) * n_mom);
    assert(NULL != bpw_h5o->ft_mom);
    for (int i = 0; i < (LDIM - 1) * n_mom ; i++)
        bpw_h5o->ft_mom[i] = ft_mom[i];
    for (int i = 0 ; i < LDIM ; i++)
        bpw_h5o->ft_x0[i] = ft_x0[i];

    bpw_h5o->buf_rank       = 3;
    bpw_h5o->buf_dims[0]    = lt;
    bpw_h5o->buf_dims[1]    = n_mom;
    bpw_h5o->buf_dims[2]    = 2;
    if (0 > (bpw_h5o->buf_dspace = 
                H5Screate_simple(bpw_h5o->buf_rank, bpw_h5o->buf_dims, NULL)))
        return "cannot create memory dataspace";

    return NULL;
}

static const char *
bpw_h5_close(lua_State *L, barpw_h5output *bpw_h5o)
{
    if (! is_masternode())
        return NULL;

    if (NULL != bpw_h5o->ft_mom)
        qlua_free(L, bpw_h5o->ft_mom);

    if (H5Sclose(bpw_h5o->buf_dspace))
        return "cannot close HDF5 membuf";

    return h5_close(L, &(bpw_h5o->h5o));
}


static const char *
bpw_h5_write(barpw_h5output *bpw_h5o, 
        const gsl_complex *v123_ft, 
        int i_v1, int i_v2, int i_v3)
{
    if (! is_masternode())
        return NULL;

    int lt = bpw_h5o->latsize[bpw_h5o->t_axis];
    hsize_t dset_off[6] = { i_v1, i_v2, i_v3,  0,              0, 0 };
    hsize_t dset_cnt[6] = {    1,    1,    1, lt, bpw_h5o->n_mom, 2 };
    if (0 != H5Sselect_hyperslab(bpw_h5o->h5o.dspace, H5S_SELECT_SET, 
                dset_off, NULL, dset_cnt, NULL))
        return "cannot select hyperslab";
    assert(H5Sget_select_npoints(bpw_h5o->h5o.dspace) == 
            H5Sget_select_npoints(bpw_h5o->buf_dspace));
    if (0 != H5Dwrite(bpw_h5o->h5o.dset, H5T_NATIVE_DOUBLE, 
                bpw_h5o->buf_dspace, bpw_h5o->h5o.dspace,
                H5P_DEFAULT, v123_ft))
        return "cannot write to dataset";

    return NULL;
}

static const char *
save_laph_wf_baryon_pwave(lua_State *L,
        mLattice *S,
        const char *h5_file,
        const char *h5_path,
        int n_v1, QDP_D3_ColorVector **v1,
        int n_v2, QDP_D3_ColorVector **v2,
        int n_v3, QDP_D3_ColorVector **v3,
        int n_mom, const int *ft_mom,     /* [n_mom][LDIM-1] */
        int t_axis,                         /* 0-based */
        const int *ft_x0                       /* ref.point for momentum phase */
        )
/* FIXME economize: if v1==v2==v3, copy (not recalculate) 5 of 6=3! times */
/* TODO save list of momenta in attributes? */
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

    /* determine subgrid parameters */
    int node_latsize[LDIM];
    int node_x0[LDIM];
    calc_subgrid(S->lat, node_latsize, node_x0);

    int vol3_dim[LDIM-1];
    int vol3_axis[LDIM-1];
    int node_lt = node_latsize[t_axis];
    int node_t0 = node_x0[t_axis];
    int vol3_local  = 1,
        vol4_local  = 1,
        vol3        = 1;
    for (int i=0, k=0 ; i < LDIM ; i++) {
        vol4_local  *= node_latsize[i];
        if (i == t_axis)
            continue;
        vol3_axis[k]= i;
        vol3_dim[k] = node_latsize[i];
        k++;
        vol3_local  *= node_latsize[i];
        vol3        *= latsize[i];
    }
    assert(QDP_sites_on_node_L(S->lat) == vol4_local);
    
    /* HDF5 output */
    barpw_h5output bpw_h5o;
    if (NULL != (err_str = bpw_h5_open_write(L, &bpw_h5o, 
                    h5_file, h5_path, latsize, t_axis, 
                    n_mom, ft_mom, ft_x0, n_v1, n_v2, n_v3))) {
        goto clearerr_0;
    }
    const char *x_attr;
    if (bpw_h5_check_meta(L, &x_attr, &bpw_h5o)) {
        err_str = "cannot update meta info";
        luaL_error(L, err_str);
        bpw_h5_close(L, &bpw_h5o);
        goto clearerr_0;
    }

    /* expose fields */
    /* FIXME handle case when fields are the same QDP objects */
    QLA_D3_ColorVector **qla_v[3];
    qla_v[0] = qlua_malloc(L, sizeof(qla_v[0][0]) *n_v1);
    for (int i = 0 ; i < n_v1 ; i++)
        qla_v[0][i] = QDP_D3_expose_V(v1[i]);
    qla_v[1] = qlua_malloc(L, sizeof(QLA_D3_ColorVector *) *n_v2);
    for (int i = 0 ; i < n_v2 ; i++)
        qla_v[1][i] = QDP_D3_expose_V(v2[i]);
    qla_v[2] = qlua_malloc(L, sizeof(QLA_D3_ColorVector *) *n_v3);
    for (int i = 0 ; i < n_v3 ; i++)
        qla_v[2][i] = QDP_D3_expose_V(v3[i]);
    /* TODO check for expose errors */

#if LDIM == 4  
    assert(LDIM == S->rank);
#define i_vol3(c) ((c)[vol3_axis[0]] - (node_x0)[vol3_axis[0]] + vol3_dim[0]*(\
                   (c)[vol3_axis[1]] - (node_x0)[vol3_axis[1]] + vol3_dim[1]*(\
                   (c)[vol3_axis[2]] - (node_x0)[vol3_axis[2]])))
#define i_X_vol3(X,c) (i_vol3(c) + vol3_local * (X))
#define i_vol4(c) i_X_vol3(((c)[t_axis] - node_t0), c)
#define mom(i_mom, k)   ((ft_mom)[(i_mom)*(LDIM-1) + k])
#else
#error "LDIM!=4 is not consistent"
#endif

    /* map i_qdp -> i_vol4 */
    int *i_vol4_qdp = NULL;
    i_vol4_qdp  = qlua_malloc(L, sizeof(i_vol4_qdp[0]) * vol4_local);
    assert(NULL != i_vol4_qdp);

    /* exp(-i\vec{p}*\vec{x})  [i_mom][i_vol3] */
    gsl_complex *exp_mipx = NULL;
    exp_mipx = qlua_malloc(L, sizeof(exp_mipx[0]) * n_mom * vol3_local);
    assert(NULL != exp_mipx);
    gsl_matrix_complex_view  gsl_exp_mipx = gsl_matrix_complex_view_array(
            (double *)exp_mipx, n_mom, vol3_local);

    for (int i_qdp = 0 ; i_qdp < vol4_local ; i_qdp++) {
        int coord[LDIM];
        QDP_get_coords_L(S->lat, coord, QDP_this_node, i_qdp);
        int i_lin = i_vol4(coord);
        i_vol4_qdp[i_qdp] = i_lin;
        if (coord[t_axis] == node_t0) {
            /* FT matrix ; XXX note the (1/V_3) factor */
            for (int i_mom = 0; i_mom < n_mom ; i_mom++) {
                double phase = 0.;
                for (int k = 0; k < LDIM-1; k++)
                    phase += mom(i_mom, k) * (coord[vol3_axis[k]] 
                            - ft_x0[vol3_axis[k]]) / (double)latsize[k];
                /* FIXME address matrix elements with gsl func/macro */ 
                exp_mipx[i_X_vol3(i_mom, coord)] = 
                            gsl_complex_polar(1. / vol3, -phase*2*M_PI);
            }
        }
    }
    
    
    QLA_D3_ColorVector *qla_v12 = NULL;
    qla_v12 = qlua_malloc(L, sizeof(qla_v12[0]) * vol4_local);

    gsl_complex *v123 = NULL;
    /* v123 [t - node_t0][i_vol3] */
    v123    = qlua_malloc(L, sizeof(v123[0]) * vol4_local);
    assert(NULL != v123);
    gsl_matrix_complex_view  gsl_v123 = gsl_matrix_complex_view_array(
                    (double *)v123, node_lt, vol3_local);

    /* Fourier-transformed v123 [t] [i_mom]
       XXX full time extent for global reduction */
    gsl_complex *v123_ft = NULL;
    v123_ft = qlua_malloc(L, sizeof(v123_ft[0]) * lt * n_mom);
    assert(NULL != v123_ft);
    gsl_matrix_complex_view  gsl_v123_ft = gsl_matrix_complex_view_array(
                    (double *)(v123_ft + node_t0 * n_mom), node_lt, n_mom);

    /* cycle over v1, v2, v3: calc v123-ft and save */
    for (int i = 0 ; i < n_v1 ; i++)
        for(int j = 0 ; j < n_v2 ; j++) {
            for (int i_site = 0; i_site < vol4_local; i_site++)
                cross_VV(qla_v12 + i_site, 
                        qla_v[0][i] + i_site, qla_v[1][j] + i_site);
            for (int k = 0; k < n_v3; k++) {
                for (int i_site = 0; i_site < vol4_local; i_site++)
                    v123[i_vol4_qdp[i_site]] = cmult_VV(qla_v12[i_site], 
                                                    qla_v[2][k][i_site]);
            
                memset(v123_ft, 0, sizeof(v123_ft[0]) * lt * n_mom);
                /* do FT */
                gsl_blas_zgemm(CblasNoTrans, CblasTrans, gsl_complex_rect(1,0), 
                        &gsl_v123.matrix, &gsl_exp_mipx.matrix,
                        gsl_complex_rect(0,0), &gsl_v123_ft.matrix);
            
                QMP_sum_double_array((double *)v123_ft, 2 * lt * n_mom);

                if (NULL != (err_str = bpw_h5_write(&bpw_h5o, v123_ft, i, j, k))) {
                    goto clearerr_2;
                }
            }
        }

clearerr_2:
    if (NULL != (err_str = bpw_h5_close(L, &bpw_h5o))) {
        goto clearerr_1;
    }

    /* cleanup */
clearerr_1:
    for (int i = 0 ; i < n_v1 ; i++)
        QDP_D3_reset_V(v1[i]);
    for (int i = 0 ; i < n_v2 ; i++)
        QDP_D3_reset_V(v2[i]);
    for (int i = 0 ; i < n_v3 ; i++)
        QDP_D3_reset_V(v3[i]);
    qlua_free(L, qla_v[0]);
    qlua_free(L, qla_v[1]);
    qlua_free(L, qla_v[2]);

    qlua_free(L, i_vol4_qdp);
    qlua_free(L, exp_mipx);
    qlua_free(L, qla_v12);
    qlua_free(L, v123);
    qlua_free(L, v123_ft);

clearerr_0:
    return err_str;

#undef i_vol3
#undef i_X_vol3
#undef i_vol4
#undef mom
}



/* 
   qcd.save_laph_wf_baryon_pwave(h5_file, h5_path, v1, v2, v3, ft_mom, t_axis, ft_x0)
   h5_file      HDF5 file name
   h5_path      keypath
   v1, v2, v3   LapH evectors
   ft_mom
   t_axis
   ft_x0
 */
int
q_save_laph_wf_baryon_pwave(lua_State *L)
{
    /* check & parse parameters */
    int argc = lua_gettop(L);
    if (8 != argc) {
        luaL_error(L, "expect 8 arguments");
        return 1;
    }

    const char *h5_file     = luaL_checkstring(L, 1);
    const char *h5_path     = luaL_checkstring(L, 2);
    int n_v1, n_v2, n_v3;
    mLattice *S = NULL;
    QDP_D3_ColorVector **v1   = qlua_check_latcolvec_table(L, 3, NULL, -1, &S, &n_v1);
    if (NULL == v1) {
        luaL_error(L, "a list of latcolvec expected in #3");
        goto clearerr_0;
    }

    QDP_D3_ColorVector **v2   = qlua_check_latcolvec_table(L, 4, S, -1, NULL, &n_v2);
    if (NULL == v2) {
        luaL_error(L, "a list of latcolvec expected in #4");
        goto clearerr_1;
    }

    QDP_D3_ColorVector **v3   = qlua_check_latcolvec_table(L, 5, S, -1, NULL, &n_v3);
    if (NULL == v3) {
        luaL_error(L, "a list of latcolvec expected in #5");
        goto clearerr_2;
    }

    int n_mom = 0;
    int *ft_mom           = qlua_check_array2d_int(L, 6, -1, LDIM - 1, &n_mom, NULL);
    if (NULL == ft_mom) {
        luaL_error(L, "list of momenta {{px,py,pz}...} expected in #6");
        goto clearerr_3;
    }

    int t_axis              = luaL_checkint(L, 7);

    int *ft_x0              = qlua_check_array1d_int(L, 8, LDIM, NULL);
    if (NULL == ft_x0) {        
        luaL_error(L, "initial FT coordinate {x0,y0,z0,t0} expected in #8");
        goto clearerr_4;
    }

    CALL_QDP(L);
    const char *status = save_laph_wf_baryon_pwave(L, S, h5_file, h5_path,
            n_v1, v1, n_v2, v2, n_v3, v3,
            n_mom, ft_mom,
            t_axis, ft_x0);
    if (NULL != status) {
        luaL_error(L, status);
        goto clearerr_5;
    }


    /* cleanup */
    qlua_free(L, v1);
    qlua_free(L, v2);
    qlua_free(L, v3);
    qlua_free(L, ft_mom);
    qlua_free(L, ft_x0);

    return 0;

clearerr_5:     qlua_free(L, ft_x0);
clearerr_4:     qlua_free(L, ft_mom);
clearerr_3:     qlua_free(L, v3);
clearerr_2:     qlua_free(L, v2);
clearerr_1:     qlua_free(L, v1);
clearerr_0:     return 1;
}


#undef LDIM
