#include <math.h>
#include <complex.h>
#include <string.h>
#include <assert.h>
#include <sys/stat.h>

#include <hdf5.h>

#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_matrix.h>

#include "qmp.h"
#include "lhpc-aff.h"

#include "modules.h"                                                /* DEPS */
#include "qlua.h"                                                   /* DEPS */
#include "lattice.h"                                                /* DEPS */
#include "qlua.h"
#include "extras.h"                                                 /* DEPS */
#include "latcolvec.h"                                              /* DEPS */


#define LDIM 4      /* number of lattice dimensions */

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

/* calculate subgrid dimensions and "lowest corner" coordinate for `lat'
   TODO check that the subgrid is rectilinear 
 */
static int 
calc_subgrid(QDP_Lattice *lat, int dims[], int cmin[])
{
    assert(LDIM == QDP_ndim_L(lat));
    int n_site = QDP_sites_on_node_L(lat);
    assert(0 < n_site);

    QDP_get_coords_L(lat, cmin, QDP_this_node, 0);
    int cmax[LDIM];
    for (int k = 0 ; k < LDIM ; k++)
        cmax[k] = cmin[k];
    
    int coord[LDIM];
    for (int i_site = 1; i_site < n_site ; i_site++) {
        QDP_get_coords_L(lat, coord, QDP_this_node, i_site);
        for (int k = 0 ; k < LDIM ; k++) {
            if (coord[k] < cmin[k])
                cmin[k] = coord[k];
            if (cmax[k] < coord[k])
                cmax[k] = coord[k];
        }
    }
    
    int vol = 1;
    for (int k = 0; k < LDIM; k++) {
        dims[k] = 1 + cmax[k] - cmin[k];
        assert(0 < dims[k]);
        vol *= dims[k];
    }
    assert(vol == n_site);

    return 0;
}

typedef struct {
    int         lt,
                n_mom;
    int         n_v1, 
                n_v2, 
                n_v3;
    
    hid_t       h5f;
    hid_t       dset;
    
    int         rank;
    hsize_t     dims[6];
    hid_t       dspace;

    int         buf_rank;
    hsize_t     buf_dims[3];
    hid_t       buf_dspace;
} h5output;

static int 
is_hsizearray_equal(int n, hsize_t *a, hsize_t *b)
{
    for (; n-- ; a++, b++)
        if (*a != *b)
            return 0;
    return 1;
}

static const char *
h5_open_v123ft_write(h5output *h5o, 
        const char *h5file, const char *h5path,
        int lt, int n_mom, 
        int n_v1, int n_v2, int n_v3)
{
    H5E_auto_t old_ehandler;
    void *old_client_data;
    hid_t error_stack = H5Eget_current_stack();

    if (0 != QDP_this_node)
        return NULL;

    assert(NULL != h5o &&
            NULL != h5file &&
            NULL != h5path);
    assert(0 < lt &&
            0 < n_mom &&
            0 < n_v1 &&
            0 < n_v2 &&
            0 < n_v3);

    /* try opening the datafile */
    H5Eget_auto(error_stack, &old_ehandler, &old_client_data);
    H5Eset_auto(error_stack, NULL, NULL);
    h5o->h5f = H5Fopen(h5file, H5F_ACC_RDWR, H5P_DEFAULT);
    H5Eset_auto(error_stack, old_ehandler, old_client_data);

    if (h5o->h5f < 0) {    /* file does not exist */
        H5Eclear(error_stack);
        if (0 > (h5o->h5f = H5Fcreate(h5file, H5F_ACC_EXCL,
                            H5P_DEFAULT, H5P_DEFAULT))) 
            return "cannot create HDF5 file";
    } 

    h5o->rank       = 6;
    h5o->dims[0]    = h5o->n_v1     = n_v1;
    h5o->dims[1]    = h5o->n_v2     = n_v2;
    h5o->dims[2]    = h5o->n_v3     = n_v3;
    h5o->dims[3]    = h5o->lt       = lt;
    h5o->dims[4]    = h5o->n_mom    = n_mom;
    h5o->dims[5]    = 2; /* real/imag */
    if (0 > (h5o->dspace = H5Screate_simple(h5o->rank, h5o->dims, NULL)))
        return "cannot create dataspace" ;
    
    /* try opening the dataset */

    H5Eget_auto(error_stack, &old_ehandler, &old_client_data);
    H5Eset_auto(error_stack, NULL, NULL);
    h5o->dset = H5Dopen(h5o->h5f, h5path, H5P_DEFAULT);
    H5Eset_auto(error_stack, old_ehandler, old_client_data);

    if (0 < h5o->dset) { /* dataset exists */
        hid_t dspace_orig;
        if (0 > (dspace_orig = H5Dget_space(h5o->dset)))
            return "cannot get dataspace from dataset";
        int rank_orig;
        rank_orig = H5Sget_simple_extent_ndims(dspace_orig);
        if (6 != rank_orig)
            return "dspace rank is not consistent";
        hsize_t dims_orig[6];
        H5Sget_simple_extent_dims(dspace_orig, dims_orig, NULL);
        if (!is_hsizearray_equal(rank_orig, dims_orig, h5o->dims))
            return "dspace dims are not consistent";
    } else { /* dataset does not exist */
        H5Eclear(error_stack);
        if (0 > (h5o->dset = H5Dcreate(h5o->h5f, h5path, H5T_IEEE_F64BE, 
                        h5o->dspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT)))
            return "cannot create dataset" ;
    }

    h5o->buf_rank       = 3;
    h5o->buf_dims[0]    = lt;
    h5o->buf_dims[1]    = n_mom;
    h5o->buf_dims[2]    = 2;
    if (0 > (h5o->buf_dspace = 
                H5Screate_simple(h5o->buf_rank, h5o->buf_dims, NULL)))
        return "cannot create memory dataspace";

    return NULL;
}

static const char *
h5_close_v123ft(h5output *h5o)
{
    if (0 != QDP_this_node)
        return NULL;

    if (H5Sclose(h5o->dspace) 
            || H5Sclose(h5o->buf_dspace)
            || H5Dclose(h5o->dset)
            || H5Fclose(h5o->h5f))
        return "cannot close HDF5 file, dataspace or dataset";
    else
        return NULL;
}


static const char *
h5_save_v123ft(h5output *h5o, const gsl_complex *v123_ft, 
        int i_v1, int i_v2, int i_v3)
{
    if (0 != QDP_this_node)
        return NULL;

    hsize_t dset_off[6] = { i_v1, i_v2, i_v3,       0,          0, 0 };
    hsize_t dset_cnt[6] = {    1,    1,    1, h5o->lt, h5o->n_mom, 2 };
    if (0 != H5Sselect_hyperslab(h5o->dspace, H5S_SELECT_SET, 
                dset_off, NULL, dset_cnt, NULL))
        return "cannot select hyperslab";
    assert(H5Sget_select_npoints(h5o->dspace) == 
            H5Sget_select_npoints(h5o->buf_dspace));
    if (0 != H5Dwrite(h5o->dset, H5T_NATIVE_DOUBLE, h5o->buf_dspace, h5o->dspace,
                H5P_DEFAULT, v123_ft))
        return "cannot write to dataset";

    return 0;
}

static const char *
save_squark_wf(lua_State *L,
        mLattice *S,
        const char *h5_file,
        const char *h5_path,
        int n_v1, QDP_D3_ColorVector **v1,
        int n_v2, QDP_D3_ColorVector **v2,
        int n_v3, QDP_D3_ColorVector **v3,
        int n_mom, const int *mom_list,     /* [n_mom][LDIM-1] */
        int t_axis,                         /* 0-based */
        int *c0                             /* ref.point for momentum phase */
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
    
    /* HDF5 output */
    h5output h5o;
    if (NULL != (err_str = h5_open_v123ft_write(&h5o, h5_file, h5_path, 
                lt, n_mom, n_v1, n_v2, n_v3))) {
        goto clearerr_0;
    }

    /* expose fields */
    QLA_D3_ColorVector **qla_v[3];
    qla_v[0] = qlua_malloc(L, sizeof(QLA_D3_ColorVector *) *n_v1);
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
#define i_vol3(c) ((c)[vol3_axis[0]] + vol3_dim[0] * (\
                   (c)[vol3_axis[1]] + vol3_dim[1] * (\
                   (c)[vol3_axis[2]])))
#define i_X_vol3(X,c) (i_vol3(c) + vol3_local * (X))
#define i_vol4(c) i_X_vol3((c)[t_axis], c)
#define mom(i_mom, k)   ((mom_list)[(i_mom)*(LDIM-1) + k])
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
        for (int i_mom = 0; i_mom < n_mom ; i_mom++) {
            double phase = 0.;
            for (int k = 0; k < LDIM-1; k++)
                phase += mom(i_mom, k) * (coord[vol3_axis[k]] - c0[vol3_axis[k]]) 
                            / (double)latsize[k];
            exp_mipx[i_X_vol3(i_mom, coord)] = gsl_complex_polar(1., -phase*2*M_PI);
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
                gsl_blas_zgemm(CblasNoTrans, CblasTrans, gsl_complex_rect(1,0), 
                        &gsl_v123.matrix, &gsl_exp_mipx.matrix,
                        gsl_complex_rect(0,0), &gsl_v123_ft.matrix);
            
                QMP_sum_double_array((double *)v123_ft, 2 * lt * n_mom);

                if (NULL != (err_str = h5_save_v123ft(&h5o, v123_ft, i, j, k))) {
                    goto clearerr_2;
                }
            }
        }

clearerr_2:
    if (NULL != (err_str = h5_close_v123ft(&h5o))) {
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


int *
qlua_check_array1d_int(lua_State *L, int idx, 
        int need_dim1, 
        int *have_dim1)
{
    if (LUA_TTABLE != lua_type(L, idx))
        luaL_error(L, "array1d(int) expected");

    int d = lua_objlen(L, idx);
    if (0 <= need_dim1 && d != need_dim1) {
        luaL_error(L, "array1d(int)[%d] expected", need_dim1);
        return NULL;
    }
    if (NULL != have_dim1)
        *have_dim1 = d;

    int *res = qlua_malloc(L, d * sizeof(int));
    for (int i = 0; i < d; i++) {
        lua_pushnumber(L, i + 1);
        lua_gettable(L, idx);
        if (lua_type(L, -1) != LUA_TNUMBER) {
            qlua_free(L, res);
            luaL_error(L, "non-int encountered in array1d(int)");
            return NULL;
        }
        res[i] = qlua_checkint(L, -1, "array element %d", i + 1);
        lua_pop(L, 1);
    }
    return res;
}


int *
qlua_check_array2d_int(lua_State *L, int idx, 
        int need_dim1, int need_dim2,
        int *have_dim1, int *have_dim2)
{
    if (LUA_TTABLE != lua_type(L, idx))
        luaL_error(L, "array2d(int) expected");

    int d1 = lua_objlen(L, idx);
    if (0 <= need_dim1 && d1 != need_dim1) {
        luaL_error(L, "array2d(int)[%d,..] expected", need_dim1, need_dim2);
        return NULL;
    }
    if (NULL != have_dim1)
        *have_dim1 = d1;
    
    lua_pushnumber(L, 1);
    lua_gettable(L, idx);
    int d2;
    int *res0   = qlua_check_array1d_int(L, lua_gettop(L), need_dim2, &d2);
    lua_pop(L, 1);
    int *res    = qlua_malloc(L, d1 * d2 * sizeof(int));
    memcpy(res, res0, d2 * sizeof(int));
    qlua_free(L, res0);

    if (NULL != have_dim2)
        *have_dim2 = d2;

    for (int i = 1; i < d1; i++) {
        lua_pushnumber(L, i + 1);
        lua_gettable(L, idx);
        int *res0 = qlua_check_array1d_int(L, lua_gettop(L), d2, NULL);
        lua_pop(L, 1);
        memcpy(res + i * d2, res0, d2 * sizeof(int));
        qlua_free(L, res0);
    }
    lua_pop(L, 1);

    return res;
}

QDP_D3_ColorVector **
qlua_check_latcolvec_table(lua_State *L, int idx, 
        mLattice *S, int need_dim1, 
        mLattice **have_S, int *have_dim1)
{
    if (LUA_TTABLE != lua_type(L, idx))
        luaL_error(L, "array1d(LatColVec) expected");
    int d = lua_objlen(L, idx);
    if (0 <= need_dim1 && d != need_dim1) {
        luaL_error(L, "array1d(LatColVec)[%d] expected", need_dim1);
        return NULL;
    }
    if (NULL != have_dim1)
        *have_dim1 = d;

    QDP_D3_ColorVector **res = qlua_malloc(L, d * sizeof(QDP_D3_ColorVector *));
    for (int i = 0; i < d; i++) {
        lua_pushnumber(L, i + 1);
        lua_gettable(L, idx);
        mLatColVec3 *ch = qlua_checkLatColVec3(L, lua_gettop(L), S, 3);
        assert(NULL != ch);
        res[i] = ch->ptr;
        mLattice *S1 = qlua_ObjLattice(L, lua_gettop(L));
        if (NULL == S && 0 == i)
            S = S1;
        else {
            if (S1->id != S->id) {
                luaL_error(L, "expect array of same lattice fields");
                return NULL;
            }
        }
        lua_pop(L, 1);
    }
    if (NULL != have_S)
        *have_S = S;

    return res;
}

/* 
   qcd.save_squark_wf(h5_file, h5_path, v1, v2, v3, mom_list, t_axis, c0)
   h5_file      HDF5 file name
   h5_path      keypath
   v1, v2, v3   LapH evectors
   mom_list
   t_axis
   c0
 */
int
q_save_squark_wf(lua_State *L)
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
    QDP_D3_ColorVector **v2   = qlua_check_latcolvec_table(L, 4, S, -1, NULL, &n_v2);
    QDP_D3_ColorVector **v3   = qlua_check_latcolvec_table(L, 5, S, -1, NULL, &n_v3);
    int n_mom = 0;
    int *mom_list           = qlua_check_array2d_int(L, 6, -1, LDIM - 1, &n_mom, NULL);
    int t_axis              = luaL_checkint(L, 7);
    int *c0                 = qlua_check_array1d_int(L, 8, LDIM, NULL);

    /* TODO call save_quark_wf */
    const char *status = save_squark_wf(L, S, h5_file, h5_path,
            n_v1, v1, n_v2, v2, n_v3, v3,
            n_mom, mom_list,
            t_axis, c0);

    /* cleanup */
    qlua_free(L, v1);
    qlua_free(L, v2);
    qlua_free(L, v3);
    qlua_free(L, mom_list);
    qlua_free(L, c0);

    if (NULL != status)
        luaL_error(L, status);

    return 0;
}


#undef LDIM
