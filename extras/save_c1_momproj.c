#include <assert.h>
#include <math.h>
#include <complex.h>
#include <string.h>

#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_matrix.h>

#include "qmp.h"
#include "lhpc-aff.h"

#include "modules.h"                                                /* DEPS */
#include "qlua.h"                                                   /* DEPS */
#include "lattice.h"                                                /* DEPS */
#include "extras.h"                                                 /* DEPS */
#include "extras_common.h"                                          /* DEPS */
#include "latcomplex.h"
#include "qparam.h"


static double complex
calc_exp_iphase_space(int ndim, int t_axis,
                const int coord[], const int c0[], 
                const int latsize[], const int mom[])
/* momentum phases for (D-1) phase transform 
    ndim        == lat. dimension D (sizes of coord, c0, latsize)
    t_axis      skipped dimension
    c0,coord    (initial) coordinate D-vector; {t_axis} component is present but ignored
    latsize     total lattice size
    mom         momentum (D-1)-vector; {t_axis} component is skipped
 */
{
    double phase = 0.;
    for (int i = 0; i < ndim; i++) {
        if (i < t_axis)
            phase += mom[i] * (double)(coord[i] - c0[i]) / latsize[i];
        else if (t_axis < i)
            phase += mom[i - 1] * (double)(coord[i] - c0[i]) / latsize[i];
    }
    phase *= 2 * M_PI;
    return cos(phase) + I * sin(phase) ;
}


static const char *
save_c1_momproj(lua_State *L, 
        mLattice *S, 
        mAffWriter *aff_w, 
        const char *aff_kpath, 
        int n_fields, 
        const char **strkey_list, 
        QDP_D_Complex **field_list, 
        int n_mom, 
        const int *mom_list, 
        const int *csrc, 
        int tlen,
        int t_axis, 
        double bc_t,
        int ft_sign)
/*  sum over sites with phase exp(i\vec p\vec x)
    XXX note the sign : analog of save_hadspec has to flip the sign of momenta!
        the sign is chosen to agree with save_bb
    strkey_list,
    field_list  [n_fields] arrays of keys,LatComplex to save
    csrc        [ndim] initial coordinates
    mom_list   array[n_mom, (ndim-1)] of momentum components
    t_axis      direction along which Fourier transf is not done
    tlen        number of time slices to save, starting from csrc[t_axis]
 */
{
    const char *errstr = NULL;
    int ndim = S->rank;
    assert(ndim <= MAX_LDIM);
    assert(t_axis < ndim);
    assert(1 == ft_sign || -1 == ft_sign);

    QDP_Lattice *lat = S->lat;
    int latsize[MAX_LDIM];
    QDP_latsize_L(lat, latsize);
    for (int i = 0 ; i < S->rank; i++) {
        if (csrc[i] < 0 || latsize[i] <= csrc[i])
            return "incorrect source coordinates";
    }

    if (4 != ndim || 
            4 != QDP_Ns ||
            3 != t_axis)
        return "not implemented for this dim, spin, color, or t-axis";

    if (NULL == aff_w ||
            NULL == aff_kpath || 
            NULL == mom_list ||
            n_mom < 0) 
        return "incorrect pointer parameters";
    
    if (n_mom <= 0)
        return NULL;       /* relax */
    
    int lt = latsize[t_axis],
        tsrc = csrc[t_axis];
    if (tlen <= 0 || lt < tlen) 
        tlen = lt;       /* default set by passing tlen<=0 or MAX_INT */

    int locsize[MAX_LDIM];
    int loc_c0[MAX_LDIM];
    calc_subgrid(lat, locsize, loc_c0);

    int loc_lt  = locsize[t_axis],
        loc_t0  = loc_c0[t_axis];
//        i_t0    = (lt + loc_t0 - tsrc) % lt;
    int *x3_idx = NULL;
    x3_idx      = locindex_space(lat, QDP_this_node, t_axis, 
            locsize, loc_c0);
    if (NULL == x3_idx) 
        return "not enough memory";

    int x3_size = 1;
    for (int d = 0 ; d < ndim ; d++)
        if (d != t_axis) x3_size *= locsize[d];
    int n_sites = QDP_sites_on_node_L(S->lat);

    /* track malloc errors of these pointers */
    double *c2lat_arr_r = NULL,
           *ph_arr_r    = NULL,
           *c2_arr_r    = NULL;
    QLA_D_Complex **exp_c2 = NULL;

    /* original lattice fields */
    /* XXX full lattice size is allocated ; can be reduced lt->loc_lt */
    /* XXX full lattice size is ZGEMM'ed at every node; can be reduced */
    int c2lat_arr_size  = 2 * lt * n_fields * x3_size * sizeof(c2lat_arr_r[0]);
    c2lat_arr_r = qlua_malloc(L, c2lat_arr_size);
    /* i_t time index counting from csrc[t_axis] */
    /* BLAS matr [{i_t,i_c2}, i_x3] ; i_t is slowest changing because local vol has subset of i_t */
#define c2lat_arr(i_c2, i_t, i_x3)   (((double complex *)c2lat_arr_r)[(i_x3) + (x3_size)*((i_c2) + (n_fields)*(i_t))])
    gsl_matrix_complex_view c2lat_gslmatr = gsl_matrix_complex_view_array(
                (double *)&c2lat_arr(0, loc_t0, 0), loc_lt * n_fields, x3_size);

    /* planewave phases */
    int ph_arr_size = 2 * n_mom * x3_size * sizeof(ph_arr_r[0]);
    ph_arr_r= qlua_malloc(L, ph_arr_size);
    /* BLAS matr [i_mom, i_x3] */
#define ph_arr(i_mom, i_x3) (((double complex *)ph_arr_r)[(i_x3) + (x3_size)*(i_mom)])
    gsl_matrix_complex_view ph_gslmatr = gsl_matrix_complex_view_array(
                ph_arr_r, n_mom, x3_size);
#define get_mom(i_mom) ((mom_list) + (ndim-1)*(i_mom))

    /* projected fields */
    int c2_arr_size = lt * 2 * n_mom * n_fields * sizeof(c2_arr_r[0]);
    c2_arr_r= qlua_malloc(L, c2_arr_size);
    /* BLAS matr [{i_t,i_c2}, i_mom] */
#define c2_arr(i_c2, i_t, i_mom) (((double complex *)c2_arr_r)[(i_mom) + (n_mom)*((i_c2) + (n_fields)*(i_t))])
    gsl_matrix_complex_view c2_gslmatr = gsl_matrix_complex_view_array(
                (double *)&c2_arr(0, loc_t0, 0), loc_lt * n_fields, n_mom);

    /* expose fields */
    exp_c2 = qlua_malloc(L, n_fields * sizeof(exp_c2[0]));
    if (    NULL == c2lat_arr_r || 
            NULL == ph_arr_r ||
            NULL == c2_arr_r ||
            NULL == exp_c2) {
        errstr = "not enough memory";
        goto clearerr_0_0;
    }

    /* init c2lat and ph */
    memset(c2lat_arr_r, 0, c2lat_arr_size);
    memset(ph_arr_r, 0, ph_arr_size);
    for (int i_c2 = 0 ; i_c2 < n_fields ; i_c2++)
        exp_c2[i_c2] = QDP_D_expose_C(field_list[i_c2]);

#pragma omp parallel for          /* XXX probably not worth omp overhead? */
    for (int i_site = 0 ; i_site < n_sites ; i_site++) {
        int i_x3 = x3_idx[i_site];
        int coord[MAX_LDIM];
        QDP_get_coords_L(lat, coord, QDP_this_node, i_site);
        /* init Fourier phase */
        for (int i_mom = 0 ; i_mom < n_mom ; i_mom++)
            ph_arr(i_mom, i_x3) = calc_exp_iphase_space(ndim, t_axis, 
                        coord, csrc, latsize, get_mom(i_mom));
        /* init orig lattice fields, including time shift and time BC */
        int tcoord = coord[t_axis];
        int bc_factor = (tsrc <= tcoord ? 1 : bc_t);
        /* TODO select only relevant timeslices? */
        for (int i_c2 = 0 ; i_c2 < n_fields ; i_c2++) {
            QLA_D_Complex qla_c = exp_c2[i_c2][i_site];
            c2lat_arr(i_c2, tcoord, i_x3) = bc_factor * (
                    QLA_real(qla_c) + I*QLA_imag(qla_c));
        }
    }
    for (int i_c2 = 0 ; i_c2 < n_fields ; i_c2++)
        QDP_D_reset_C(field_list[i_c2]);

    /* project fields */
    memset(c2_arr_r, 0, c2_arr_size);
    if (1 == ft_sign)
        gsl_blas_zgemm(CblasNoTrans, CblasTrans, 
                gsl_complex_rect(1,0), &c2lat_gslmatr.matrix, &ph_gslmatr.matrix, 
                gsl_complex_rect(0,0), &c2_gslmatr.matrix);
    else if (-1 == ft_sign)
        gsl_blas_zgemm(CblasNoTrans, CblasConjTrans, 
                gsl_complex_rect(1,0), &c2lat_gslmatr.matrix, &ph_gslmatr.matrix, 
                gsl_complex_rect(0,0), &c2_gslmatr.matrix);
    
    /* global sum */
    if (QMP_sum_double_array(c2_arr_r, c2_arr_size / sizeof(double))) {
        errstr = "QMP_sum_double_array error";
        goto clearerr_0_0;
    }
    
    /* save to AFF
       TODO add hdf5 as an option, perhaps splitting this function
       into computing and writing parts */
    if (aff_w->master) {
        double complex *cplx_buf = NULL;
        char *momkey_list = NULL;

        /* check that path can be created */
        struct AffNode_s *aff_top = NULL;
        aff_top = aff_writer_mkpath(aff_w->ptr, aff_w->dir, aff_kpath);
        if (NULL == aff_top) {
            errstr = aff_writer_errstr(aff_w->ptr);
            goto clearerr_1_0;
        }

        cplx_buf = qlua_malloc(L, tlen * sizeof(double complex));
        int momkey_len  = 6*MAX_LDIM;
        momkey_list = qlua_malloc(L, momkey_len * n_mom);
#define get_momkey(i_mom) ((momkey_list) + (i_mom)*(momkey_len))
        if (NULL == cplx_buf || NULL == momkey_list) {
            errstr = "not enough memory" ;
            goto clearerr_1_1;
        }
        /* generate momkey strings */
        for (int i_mom = 0 ; i_mom < n_mom ; i_mom++) {
            const int *mom = get_mom(i_mom);
            /* FIXME support arbitrary momenta: need to change momkey format */
            assert(4 == ndim);          
            snprintf(get_momkey(i_mom), momkey_len,
                    "PX%d_PY%d_PZ%d", mom[0], mom[1], mom[2]);
        }
        /* save data */
        for (int i_c2 = 0; i_c2 < n_fields; i_c2++) {
            struct AffNode_s *aff_c2 = aff_writer_mkpath(
                        aff_w->ptr, aff_top, strkey_list[i_c2]);
            if (NULL == aff_c2) {
                errstr = aff_writer_errstr(aff_w->ptr);
                goto clearerr_1_1;
            }

            for (int i_mom = 0; i_mom < n_mom; i_mom++) {
                for (int i_t = 0 ; i_t < tlen ; i_t++)
                    cplx_buf[i_t] = c2_arr(i_c2, (tsrc + i_t)%lt, i_mom);
                
                struct AffNode_s *node = aff_writer_mkdir(aff_w->ptr, aff_c2, get_momkey(i_mom));
                if (NULL == node) {
                    errstr = aff_writer_errstr(aff_w->ptr);
                    goto clearerr_1_1;
                }
                if (aff_node_put_complex(aff_w->ptr, node, cplx_buf, tlen)) {
                    errstr = aff_writer_errstr(aff_w->ptr);
                    goto clearerr_1_1;
                }
            }
        }
clearerr_1_1:
        qlua_free(L, momkey_list);
        qlua_free(L, cplx_buf);
clearerr_1_0:
        ;
    }

clearerr_0_0:
    if (NULL != c2lat_arr_r)    qlua_free(L, c2lat_arr_r);
    if (NULL != ph_arr_r)       qlua_free(L, ph_arr_r);
//clearerr_0_0:
    if (NULL != c2_arr_r)       qlua_free(L, c2_arr_r);
//clearerr_0_0:
    if (NULL != exp_c2)         qlua_free(L, exp_c2);
    if (NULL != x3_idx)         free(x3_idx);

#undef c2lat_arr
#undef ph_arr
#undef c2_arr
#undef get_mom
#undef get_momkey

    return errstr;
}


/* save momentum projections of latcomplex to aff_kpath_prefix/relkpath/PX_PY_PZ
   
   save_momproj(
       aff_writer,                  -- 1 TODO make hdf5 version
       aff_kpath_prefix,            -- 2 
       {relkpath->latcomplex},      -- 3
       csrc,                        -- 4
       mom_list,                    -- 5
       t_axis,                      -- 6
       bc_t,                        -- 7
       {                            -- 8 options
         tlen=-1, 
         ft_sign=1 }
       )
 */
int
q_save_momproj(lua_State *L)
{
    int argc = lua_gettop(L);
    if (argc < 7) {
        luaL_error(L, "expect 7 arguments");
        return 1;
    }

    mAffWriter *aff_w = qlua_checkAffWriter(L, 1);
    const char *aff_kpath = luaL_checkstring(L, 2);

    int n_fields    = -1;
    const char **strkey_list = NULL;
    QDP_Complex **field_list = NULL;
    mLattice *S     = NULL; 
    qlua_check_latcomplex_strmap(L, 3, &S, &n_fields, &strkey_list, &field_list);

    int *csrc       = qlua_check_array1d_int(L, 4, S->rank, NULL);
    if (csrc == NULL)
        return luaL_error(L, "bad value for coord_src");
    int n_mom       = -1;
    int *mom_list   = qlua_check_array2d_int(L, 5, -1, S->rank - 1, &n_mom, NULL);
    int t_axis      = luaL_checkint(L, 6);
    double bc_t     = luaL_checknumber(L, 7);

    int tlen    = -1,
        ft_sign =  1;
    if (qlua_checkopt_paramtable(L, 8)) {
        tlen    = qlua_tabkey_intopt(L, 8, "tlen",  -1);
        ft_sign = qlua_tabkey_intopt(L, 8, "ft_sign",   1);
    }                                                           /* [-0,+0,-] */
    

//    printf("csrc={");for(int i=0; i<4; i++) printf("%d", csrc[i]);printf("}\n");
//    for(int j=0; j<n_mom; j++) { printf("mom[%d]={", j);for(int i=0; i<3; i++) printf("%d", mom_list[3*j+i]);printf("}\n"); }
//    printf("t_axis=%d\n", t_axis);
//    printf("bc_baryon=%f\n", bc_baryon);
//    printf("aff_kpath=%s\n", aff_kpath);

    qlua_Aff_enter(L);
    CALL_QDP(L);
    
    const char *status = NULL;
    status = save_c1_momproj(
            L, S, aff_w, aff_kpath, 
            n_fields, strkey_list, field_list, 
            n_mom, mom_list, 
            csrc, tlen, t_axis, bc_t, ft_sign);
    qlua_Aff_leave();
    
    qlua_free(L, csrc);
    qlua_free(L, mom_list);
    qlua_free(L, strkey_list);
    qlua_free(L, field_list);

    if (status)
        luaL_error(L, status);

    return 0;
}

