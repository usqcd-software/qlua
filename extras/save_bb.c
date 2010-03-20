#include <math.h>
#include <complex.h>
#include <string.h>

#include "qmp.h"
#include "lhpc-aff.h"

#include "qlua.h"                                                   /* DEPS */
#include "aff_io.h"                                                 /* DEPS */
#include "extras.h"                                                 /* DEPS */

static double complex
calc_exp_iphase(const int coord[], const int c0[], 
                const int latsize[], const int mom[])
{
    double phase = 0.;
    int i;
    for (i = 0; i < 4; i++)
        phase += mom[i] * (double)(coord[i] - c0[i]) / latsize[i];
    phase *= 2 * M_PI;
    return cos(phase) + I * sin(phase) ;
}




/* optimized version of building blocks 
   compute tr(B^+ \Gamma_n F) [n=0..15] and projects on n_qext momenta
   select time interval [tsrc: tsnk] and does time reversal if time_rev==1
   save results to aff_w[aff_kpath . 'g%d/qx%d_qy%d_qz%d']
   Parameters:
    csrc = { xsrc, ysrc, zsrc, tsrc }
    tsnk
    qext[4 * i_qext + dir]  ext.mom components
    time_rev ==0 for proton_3, ==1 for proton_negpar_3
    bc_baryon_t =+/-1 boundary condition for baryon 2pt[sic!] function; =bc_quark^3
 */
const char *
save_bb(lua_State *L,
        mAffWriter *aff_w,
        const char *aff_kpath,
        QDP_DiracPropagator *F,
        QDP_DiracPropagator *B,
        const int *csrc,             /* [qRank] */
        int tsnk,
        int n_mom,
        const int *mom,             /* [n_mom][qRank] */
        int time_rev,                /* 1 to reverse, 0 to not */
        int t_axis,                  /* 0-based */
        double bc_baryon_t)
{
#define trace_with_Gamma(dst, s4x4) do { \
    (dst)[ 0]   = +  (s4x4)[0][0] +  (s4x4)[1][1] +  (s4x4)[2][2] +  (s4x4)[3][3];\
    (dst)[ 1]   = +I*(s4x4)[3][0] +I*(s4x4)[2][1] -I*(s4x4)[1][2] -I*(s4x4)[0][3];\
    (dst)[ 2]   = -  (s4x4)[3][0] +  (s4x4)[2][1] +  (s4x4)[1][2] -  (s4x4)[0][3];\
    (dst)[ 3]   = -I*(s4x4)[0][0] +I*(s4x4)[1][1] -I*(s4x4)[2][2] +I*(s4x4)[3][3];\
    (dst)[ 4]   = +I*(s4x4)[2][0] -I*(s4x4)[3][1] -I*(s4x4)[0][2] +I*(s4x4)[1][3];\
    (dst)[ 5]   = -  (s4x4)[1][0] +  (s4x4)[0][1] -  (s4x4)[3][2] +  (s4x4)[2][3];\
    (dst)[ 6]   = -I*(s4x4)[1][0] -I*(s4x4)[0][1] -I*(s4x4)[3][2] -I*(s4x4)[2][3];\
    (dst)[ 7]   = +  (s4x4)[2][0] +  (s4x4)[3][1] -  (s4x4)[0][2] -  (s4x4)[1][3];\
    (dst)[ 8]   = +  (s4x4)[2][0] +  (s4x4)[3][1] +  (s4x4)[0][2] +  (s4x4)[1][3];\
    (dst)[ 9]   = +I*(s4x4)[1][0] +I*(s4x4)[0][1] -I*(s4x4)[3][2] -I*(s4x4)[2][3];\
    (dst)[10]   = -  (s4x4)[1][0] +  (s4x4)[0][1] +  (s4x4)[3][2] -  (s4x4)[2][3];\
    (dst)[11]   = -I*(s4x4)[2][0] +I*(s4x4)[3][1] -I*(s4x4)[0][2] +I*(s4x4)[1][3];\
    (dst)[12]   = +I*(s4x4)[0][0] -I*(s4x4)[1][1] -I*(s4x4)[2][2] +I*(s4x4)[3][3];\
    (dst)[13]   = -  (s4x4)[3][0] +  (s4x4)[2][1] -  (s4x4)[1][2] +  (s4x4)[0][3];\
    (dst)[14]   = -I*(s4x4)[3][0] -I*(s4x4)[2][1] -I*(s4x4)[1][2] -I*(s4x4)[0][3];\
    (dst)[15]   = +  (s4x4)[0][0] +  (s4x4)[1][1] -  (s4x4)[2][2] -  (s4x4)[3][3];\
    } while(0)

#define get_mom(mom_list, i_mom) ((mom_list) + 4*(i_mom))
    if (4 != QDP_ndim() || 
            4 != QDP_Ns ||
            3 != QDP_Nc ||
            3 != t_axis) {
        return "not implemented for this dim, spin, color, or t-axis";
    }

    int latsize[4];
    QDP_latsize(latsize);
    if (NULL == aff_w ||
            NULL == aff_kpath || 
            NULL == mom ||
            n_mom < 0) {
        return "incorrect pointer parameters";
    }
    int i;
    for (i = 0 ; i < 4; i++) {
        if (csrc[i] < 0 || latsize[i] <= csrc[i]) {
            return "incorrect source coordinates";
        }
    }
    if (tsnk < 0 || latsize[t_axis] <= tsnk) {
        return "incorrect sink t-coordinate";
    }
    
    if (n_mom <= 0)
        return NULL;       /* relax */

    int src_snk_dt = -1;
    int lt = latsize[t_axis];
    if (!time_rev) {
        src_snk_dt = (lt + tsnk - csrc[t_axis]) % lt;
    } else {
        src_snk_dt = (lt + csrc[t_axis] - tsnk) % lt;
    }
           
    int bb_arr_size = 16 * n_mom * (src_snk_dt + 1) * 2 * sizeof(double);
    double *bb_arr = qlua_malloc(L, bb_arr_size);
    if (NULL == bb_arr) {
        return "not enough memory";
    }

    memset(bb_arr, 0, bb_arr_size);
#define bb_real(i_gamma, i_mom) ((bb_arr) + (src_snk_dt + 1) * (0 + 2 * ((i_mom) + n_mom * (i_gamma))))
#define bb_imag(i_gamma, i_mom) ((bb_arr) + (src_snk_dt + 1) * (1 + 2 * ((i_mom) + n_mom * (i_gamma))))

    double complex *exp_iphase = qlua_malloc(L, n_mom * sizeof(double complex));
    
    int coord[4];
    double complex trc_FBd[4][4];
    QLA_DiracPropagator *F_exp = QDP_expose_P(F);
    QLA_DiracPropagator *B_exp = QDP_expose_P(B);

    int i_site;
    for (i_site = 0; i_site < QDP_sites_on_node; i_site++) {
        QDP_get_coords(coord, QDP_this_node, i_site);
        
        int t = -1;
        if (!time_rev) {
            t = (lt + coord[t_axis] - csrc[t_axis]) % lt;
        } else {
            t = (lt + csrc[t_axis] - coord[t_axis]) % lt;
        }
        if (src_snk_dt < t)
            continue;

        /* precalc phases for inner contraction loop */
        int i_mom;
        for (i_mom = 0 ; i_mom < n_mom ; i_mom++) {
            exp_iphase[i_mom] = calc_exp_iphase(coord, csrc, 
                    latsize, get_mom(mom, i_mom));
        }

        /* compute trace_{color} [ F * B^\dag] 
           [is,js]  = sum_{ic,jc,ks} F[ic,is; jc,ks] * (B[ic,js; jc,ks])^*
           is,js,ks - spin, ic,jc - color
         */
        int is, js, ks, ic, jc;
        for (is = 0; is < 4; is++) {
            for (js = 0; js < 4; js++) {
                QLA_Complex sum;
                QLA_c_eq_r(sum, 0);
                for (ks = 0; ks < 4; ks++) {
                    for (ic = 0; ic < 3 ; ic++)
                        for (jc = 0; jc < 3 ; jc++)
                            QLA_c_peq_c_times_ca(sum, 
                                    QLA_elem_P(F_exp[i_site], ic,is, jc,ks),
                                    QLA_elem_P(B_exp[i_site], ic,js, jc,ks));
                }
                trc_FBd[is][js] = QLA_real(sum) + I*QLA_imag(sum);
            }
        }

        /* cycle over Gamma */
        int gn;
        double complex tr_G_F_Bd[16];
        trace_with_Gamma(tr_G_F_Bd, trc_FBd);
        for (gn = 0; gn < 16 ; gn++) {
            /* mult. by phase and add to timeslice sum */
            for (i_mom = 0; i_mom < n_mom; i_mom++) {
                /* FIXME do Fourier transf with optimized matrix-dot
                   exp_ipx[\vec q, \vec x] . bb[\vec x, {Gamma,tau}] 
                   #\vec q ~ 147,  #\vec x ~ 128, #{Gamma,tau} ~ 96
                 */
                double complex aux = exp_iphase[i_mom] * tr_G_F_Bd[gn];
                bb_real(gn, i_mom)[t] += creal(aux);
                bb_imag(gn, i_mom)[t] += cimag(aux);
            }
        }
    }
    
    qlua_free(L, exp_iphase);

    /* global sum */
    if (QMP_sum_double_array(bb_arr, bb_arr_size / sizeof(double))) {
        qlua_free(L, bb_arr);
        return "QMP_sum_double_array error";
    }
    
    /* save to AFF */
    if (aff_w->master) {
        struct AffNode_s *aff_top = aff_writer_mkpath(aff_w->ptr, 
                                            aff_w->dir, aff_kpath);
        if (NULL == aff_top) {
            qlua_free(L, bb_arr);
            return aff_writer_errstr(aff_w->ptr);
        }

        double complex *cplx_buf = qlua_malloc(L, (src_snk_dt + 1) * sizeof(double complex));
        if (NULL == cplx_buf) {
            qlua_free(L, bb_arr);
            return "not enough memory";
        }
        int gn, i_mom, t;

        /* generate qext keys */
#define MOM_STR_STRIDE   32
        char *mom_str_ = qlua_malloc(L, n_mom * MOM_STR_STRIDE * sizeof(char));
        if (NULL == mom_str_) {
            qlua_free(L, bb_arr);
            qlua_free(L, cplx_buf);
            return "not enough memory";
        }
#define MOM_STR(i_mom) ((mom_str_) + (i_mom) * (MOM_STR_STRIDE))
        for (i_mom = 0; i_mom < n_mom; i_mom++)
            snprintf(MOM_STR(i_mom), MOM_STR_STRIDE,
                     "qx%d_qy%d_qz%d", get_mom(mom, i_mom)[0], 
                     get_mom(mom, i_mom)[1], get_mom(mom, i_mom)[2]);

        /* loop over gamma and qext */
        for (gn = 0; gn < 16; gn++) {
            char buf[16];
            snprintf(buf, sizeof(buf), "g%d", gn);
            struct AffNode_s *node_gn = aff_writer_mkdir(aff_w->ptr, aff_top, buf);
            if (NULL == node_gn) {
                qlua_free(L, bb_arr);
                qlua_free(L, cplx_buf);
                qlua_free(L, mom_str_);
                return aff_writer_errstr(aff_w->ptr);
            }
            for (i_mom = 0; i_mom < n_mom; i_mom++) {
                /* copy & mult by bc, if necessary */
                /* FIXME avoid copying: transform in-place and send to AFF */
                const double *bb_re_cur = bb_real(gn, i_mom),
                             *bb_im_cur = bb_imag(gn, i_mom);
                if (!time_rev) {    /* no bc */
                    for (t = 0 ; t <= src_snk_dt; t++)
                        cplx_buf[t] = bb_re_cur[t] + I*bb_im_cur[t];
                } else {
                    if (gn < 8) {
                        for (t = 0 ; t <= src_snk_dt; t++)
                            cplx_buf[t] = bc_baryon_t * (bb_re_cur[t] + I*bb_im_cur[t]);
                    } else {
                        for (t = 0 ; t <= src_snk_dt; t++)
                            cplx_buf[t] = -bc_baryon_t * (bb_re_cur[t] + I*bb_im_cur[t]);
                    }
                }
                /* write to AFF */
                struct AffNode_s *node_gn_q = aff_writer_mkdir(aff_w->ptr, 
                                    node_gn, MOM_STR(i_mom));
                if (NULL == node_gn_q) {
                    qlua_free(L, bb_arr);
                    qlua_free(L, cplx_buf);
                    qlua_free(L, mom_str_);
                    return aff_writer_errstr(aff_w->ptr);
                }
                if (aff_node_put_complex(aff_w->ptr, node_gn_q, cplx_buf, src_snk_dt + 1)) {
                    qlua_free(L, bb_arr);
                    qlua_free(L, cplx_buf);
                    qlua_free(L, mom_str_);
                    return aff_writer_errstr(aff_w->ptr);
                }
            }
        }
        qlua_free(L, cplx_buf);
        qlua_free(L, mom_str_);
    }

#undef bb_real
#undef bb_imag
#undef get_mom
    qlua_free(L, bb_arr);   
    QDP_reset_P(F);
    QDP_reset_P(B);
    return 0;
}
