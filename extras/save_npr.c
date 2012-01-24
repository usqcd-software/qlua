#include <math.h>
#include <complex.h>
#include <string.h>
#include <assert.h>

#include "qmp.h"
#include "lhpc-aff.h"

#include "modules.h"                                                /* DEPS */
#include "qlua.h"                                                   /* DEPS */
#include "lattice.h"                                                /* DEPS */
#include "extras.h"                                                 /* DEPS */
#include "latdirprop.h"
#include "gamma_defs.h"
#include "extras_common.h"
#include "qparam.h"

#define NSPIN   4
#define NCOLOR  3
#define LDIM    4

#if 0
static gsl_complex *
gsl_make_exp_ipx(lua_State *L, QDP_Lattice *lat,
            const int ft_mom[], const int ft_x0[])
{
    int n_site = QDP_sites_on_node_L(lat);
    gsl_complex *res = qlua_malloc(L, sizeof(gsl_complex) * n_site);
    int *latsize=qlua_malloc(L, sizeof(QDP_Int *) * QDP_ndim_L(lat));
    int *coord = qlua_malloc(L, sizeof(QDP_Int *) * QDP_ndim_L(lat));
    for (int i_site = 0; i_site < n_site; i_site++) {
        double phase = 0.;
        QDP_get_coords_L(lat, coord, QDP_this_node, i_site);
        for (int d = 0 ; d < QDP_ndim() ; d++) 
            phase += 2 * M_PI * (coord[d] - ft_x0[d]) * (double)ft_mom[d] / latsize[d];
        res[i_site] = gsl_complex_polar(1., phase);
    }
    qlua_free(L, coord);
    qlua_free(L, latsize);
    return res;
}
#endif

static const char *
save_npr_2qvertex(lua_State *L, mAffWriter *aff_w, const char *aff_kpath,
        QDP_D3_DiracPropagator *frw_prop, QDP_D3_DiracPropagator *bkw_prop,
        const int ft_mom[], const int ft_x0[])
{
    const char *err_str = NULL;

    QDP_Lattice *lat = QDP_get_lattice_P(frw_prop);
    int vol4 = QDP_volume_L(lat);
#define NGAMMA  (NSPIN*NSPIN)
#define iGamma_5 (NGAMMA-1)

    double complex *gsl_qbarq = qlua_malloc(L, sizeof(double complex) * 
                    NSPIN*NSPIN * NCOLOR*NCOLOR);
    /* order in AFF file 
       [ic][is]qbar . Gamma[i_Gamma] . q[jc,js] */ 
#define QbarQ(ic, is, jc, js) (gsl_qbarq[\
                    (js) + NSPIN*((jc) + NCOLOR *(\
                    (is) + NSPIN*(ic)))])

    QDP_D_Complex *exp_ipx = QDP_D_create_C_L(lat);
    calc_D_C_eq_planewave(exp_ipx, ft_mom, ft_x0); /* sic! will be conj'ed later */

    QDP_D3_DiracPropagator *dirprop_aux1 = QDP_D3_create_P_L(lat);
    QDP_D3_DiracPropagator *dirprop_aux2 = QDP_D3_create_P_L(lat);
    QDP_D3_DiracPropagator *dirprop_aux3 = QDP_D3_create_P_L(lat);

    QDP_D3_P_eq_P_times_gamma(dirprop_aux1, bkw_prop, iGamma_5, QDP_all_L(lat));
    QDP_D3_P_eq_gamma_times_P(dirprop_aux2, dirprop_aux1, iGamma_5, QDP_all_L(lat));
    QDP_D3_P_eq_C_times_P(dirprop_aux1, exp_ipx, dirprop_aux2, QDP_all_L(lat));

    struct AffNode_s *aff_top = NULL;
    if (aff_w->master) {
        aff_top = aff_writer_mkpath(aff_w->ptr, aff_w->dir, aff_kpath);
        if (NULL == aff_top) {
            err_str = aff_writer_errstr(aff_w->ptr);
            goto clearerr_1;
        }
    }
    char buf[1000];

    for (int iG = 0 ; iG < NGAMMA ; iG++) {
        QDP_P_eq_gamma_times_P(dirprop_aux2, frw_prop, iG, QDP_all_L(lat));
        QDP_P_eq_Pa_times_P(dirprop_aux3, dirprop_aux1, dirprop_aux2, QDP_all_L(lat));
        QLA_DiracPropagator qla_qbarq;
        QDP_p_eq_sum_P(&qla_qbarq, dirprop_aux3, QDP_all_L(lat));
        for (int ic = 0 ; ic < NCOLOR ; ic++ )  for (int jc = 0 ; jc < NCOLOR ; jc++) 
            for (int is = 0 ; is < NSPIN ; is++)  for (int js = 0 ; js < NSPIN ; js++)
                QbarQ(ic, is, jc, js) = qla2complex(QLA_elem_P(qla_qbarq, ic, is, jc, js)) / vol4;
        /* save to AFF */
        if (aff_w->master) {
            snprintf(buf, sizeof(buf), "g%d", iG);
            struct AffNode_s *node = aff_writer_mkpath(aff_w->ptr, aff_top, buf);
            if (NULL == node) {
                err_str = aff_writer_errstr(aff_w->ptr);
                goto clearerr_1;
            }
            if (aff_node_put_complex(aff_w->ptr, node, gsl_qbarq, 
                        NSPIN*NSPIN * NCOLOR*NCOLOR)) {
                err_str = aff_writer_errstr(aff_w->ptr);
                goto clearerr_1;
            }
        }
    }


    qlua_free(L, gsl_qbarq);
    QDP_D_destroy_C(exp_ipx);
    QDP_D3_destroy_P(dirprop_aux1);
    QDP_D3_destroy_P(dirprop_aux2);
    QDP_D3_destroy_P(dirprop_aux3);
    return NULL;

clearerr_1:
    qlua_free(L, gsl_qbarq);
    QDP_D_destroy_C(exp_ipx);
    QDP_D3_destroy_P(dirprop_aux1);
    QDP_D3_destroy_P(dirprop_aux2);
    QDP_D3_destroy_P(dirprop_aux3);
    return err_str;
}

/* semantics 
   save_npr_2qvertex(aff_writer, aff_kpath, frw_prop, bkw_prop, ft_mom, ft_x0) 
   * aff_writer,aff_kpath   
            destination
   * frw_prop,bkw_prop      
            propagators
   * ft_mom,ft_x0:          
            momentum projection: sum_x vertex(x) * exp(-i*ft_mom.(x-ft_x0))
*/
int 
q_save_npr_2qvertex(lua_State *L)
{
    int argc = lua_gettop(L); 
    if (6 != argc) {
        luaL_error(L, "expect 5 arguments");
        return 1;
    }

    mAffWriter *aff_w = qlua_checkAffWriter(L, 1);
    const char *key_path = luaL_checkstring(L, 2);
    
    mLatDirProp3 *F = qlua_checkLatDirProp3(L, 3, NULL, 3);
    mLattice *S     = qlua_ObjLattice(L, 3);
    mLatDirProp3 *B = qlua_checkLatDirProp3(L, 4, S, 3);

    int *ft_mom     = qlua_checkintarray(L, 5, S->rank, NULL);
    int *ft_x0      = qlua_checkintarray(L, 6, S->rank, NULL);

    qlua_Aff_enter(L);
    CALL_QDP(L);
    const char *status = save_npr_2qvertex(L, aff_w, key_path,
            F->ptr, B->ptr, ft_mom, ft_x0);
    qlua_Aff_leave();

    if (NULL != status) {
        luaL_error(L, status);
        goto clearerr_1;
    }

    qlua_free(L, ft_mom);
    qlua_free(L, ft_x0);

    return 0;

clearerr_1:

    qlua_free(L, ft_mom);
    qlua_free(L, ft_x0);

    return 1;
}


static const char *
save_npr_prop(lua_State *L, mAffWriter *aff_w, const char *aff_kpath,
        QDP_D3_DiracPropagator *prop, const int ft_mom[], const int ft_x0[])
{
    const char *err_str = NULL;

    QDP_Lattice *lat = QDP_get_lattice_P(prop);
    int vol4 = QDP_volume_L(lat);
    
    double complex *gsl_prop = qlua_malloc(L, sizeof(double complex) * 
                    NSPIN*NSPIN * NCOLOR*NCOLOR);
    /* order in AFF file [ic, is, jc, js] */ 
#define buf_prop(ic, is, jc, js) (gsl_prop[\
                    (js) + NSPIN*((jc) + NCOLOR *(\
                    (is) + NSPIN*(ic)))])

    QDP_D_Complex *exp_mipx = QDP_D_create_C_L(lat);
    int *ft_m_mom = qlua_malloc(L, sizeof(int) * LDIM);
    for (int i = 0 ; i < LDIM ; i++) 
        ft_m_mom[i] = -ft_mom[i];
    calc_D_C_eq_planewave(exp_mipx, ft_m_mom, ft_x0);
    qlua_free(L, ft_m_mom);

    QDP_D3_DiracPropagator *dirprop_aux = QDP_D3_create_P_L(lat);
    QDP_D3_P_eq_C_times_P(dirprop_aux, exp_mipx, prop, QDP_all_L(lat));
    QLA_D3_DiracPropagator qla_prop;
    QDP_D3_p_eq_sum_P(&qla_prop, dirprop_aux, QDP_all_L(lat));
    for (int ic = 0 ; ic < NCOLOR ; ic++ ) for (int jc = 0 ; jc < NCOLOR ; jc++) 
        for (int is = 0 ; is < NSPIN ; is++) for (int js = 0 ; js < NSPIN ; js++)
            buf_prop(ic, is, jc, js) = 
                qla2complex(QLA_elem_P(qla_prop, ic, is, jc, js)) / vol4;

    /* save to AFF */
    if (aff_w->master) {
        struct AffNode_s *node = aff_writer_mkpath(aff_w->ptr, aff_w->dir, aff_kpath);
        if (NULL == node) {
            err_str = aff_writer_errstr(aff_w->ptr);
            goto clearerr_1;
        }
        if (aff_node_put_complex(aff_w->ptr, node, gsl_prop, 
                    NSPIN*NSPIN * NCOLOR*NCOLOR)) {
            err_str = aff_writer_errstr(aff_w->ptr);
            goto clearerr_1;
        }
    }
    
    qlua_free(L, gsl_prop);
    QDP_D_destroy_C(exp_mipx);
    QDP_D3_destroy_P(dirprop_aux);
    return NULL;

clearerr_1:
    qlua_free(L, gsl_prop);
    QDP_D_destroy_C(exp_mipx);
    QDP_D3_destroy_P(dirprop_aux);
    return err_str;
}
/* semantics 
   save_npr_prop(aff_writer, aff_kpath, prop, ft_mom, ft_x0) 
   * aff_writer,aff_kpath   
            destination
   * prop   
            propagators
   * ft_mom,ft_x0:          
            momentum projection: sum_x vertex(x) * exp(-i*ft_mom.(x-ft_x0))
*/
int 
q_save_npr_prop(lua_State *L)
{
    int argc = lua_gettop(L); 
    if (5 != argc) {
        luaL_error(L, "expect 5 arguments");
        return 1;
    }

    mAffWriter *aff_w = qlua_checkAffWriter(L, 1);
    const char *key_path = luaL_checkstring(L, 2);
    
    mLatDirProp3 *F = qlua_checkLatDirProp3(L, 3, NULL, 3);

    int *ft_mom     = qlua_check_array1d_int(L, 4, LDIM, NULL);
    int *ft_x0      = qlua_check_array1d_int(L, 5, LDIM, NULL);

    qlua_Aff_enter(L);
    CALL_QDP(L);
    const char *status = NULL;
    status = save_npr_prop(L, aff_w, key_path, F->ptr, ft_mom, ft_x0);
    qlua_Aff_leave();

    if (NULL != status) {
        luaL_error(L, status);
        goto clearerr_1;
    }

    qlua_free(L, ft_mom);
    qlua_free(L, ft_x0);

    return 0;

clearerr_1:

    qlua_free(L, ft_mom);
    qlua_free(L, ft_x0);

    return 1;
}
