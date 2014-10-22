#include <assert.h>
#include <string.h>
#include <complex.h>

#include <stdio.h> /**/

#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>

#include "qmp.h"

#include "modules.h"                                                /* DEPS */
#include "qlua.h"                                                   /* DEPS */
#include "lattice.h"                                                /* DEPS */
#include "aff_io.h"                                                 /* DEPS */
#include "extras.h"                                                 /* DEPS */
#include "qlayout.h"                                                /* DEPS */


/* lattice (L->L) Fourier Transform of a Complex 
    y_p <- \sum_a(e^{i p*a} x_a) : only along 'dir'
    lat_c_y, lat_c_x : lat_c_y <- FT(lat_c_x)
            QLA_D_Complex arrays [site_ind*rec_len + rec_ind] 
            where site_ind is QDP site index within "exposed" fields and
            0<=rec_ind<rec_len is internal "field" index (e.g. 1 for Complex, 
            12 for DirFerm, 144 for DirProp, etc)
    rec_len:   field var size (# of complex numbers)
    ft_sign:    sign of Fourier transform to do
    ft_dir:     0<=ft_dir<NDIM: direction of ft transform
 */
const char *
extra_lat_fourier_qla_onedim(lua_State *L,
        mLattice *S, QLA_D_Complex *qla_y, QLA_D_Complex *qla_x, 
        int rec_len, int ft_sign, int ft_dir)
{
    const char *status = NULL;
    int ndim = QDP_ndim_L(S->lat);

    assert(0 <= ft_dir && ft_dir < ndim);
    /* collect all subvolumes along ft_dir by near-neighbor comm */
    /* assume that all subvol.dims are equal on all nodes */
    int subvol_lo[QLUA_MAX_LATTICE_RANK],
        subvol_hi[QLUA_MAX_LATTICE_RANK],
        subvol_ls[QLUA_MAX_LATTICE_RANK],
        c_ind_mul[QLUA_MAX_LATTICE_RANK],
        c0_bkw[QLUA_MAX_LATTICE_RANK],
        c0_frw[QLUA_MAX_LATTICE_RANK],
        ls[QLUA_MAX_LATTICE_RANK];
    
    /* FIXME init geometry using qlua or QDP/QMP functions? */
    /* get subvol size and create indices for "navigating" buffers */
    QDP_latsize_L(S->lat, ls);
    qlua_sublattice(subvol_lo, subvol_hi, QDP_this_node, (void *)S);

    /* init geometry and indices */
    int dir_mul, subvol, ft_len, c0_ft;
    subvol  = 1;
    dir_mul = 1;
    for (int mu = 0; mu < ndim ; mu++) {
        subvol_ls[mu] = subvol_hi[mu] - subvol_lo[mu];
        subvol  *= subvol_ls[mu];
        if (mu != ft_dir) {
            c_ind_mul[mu] = dir_mul;
            dir_mul *= subvol_ls[mu];
            
            c0_bkw[mu] = subvol_lo[mu];
            c0_frw[mu] = subvol_lo[mu];
        }
        else {
            c_ind_mul[mu] = 0;

            c0_bkw[mu] = (subvol_lo[mu] + ls[mu] - subvol_ls[mu]) % ls[mu];
            c0_frw[mu] = (subvol_lo[mu] + subvol_ls[mu]) % ls[mu];
        }
    }
    int node_bkw = QDP_node_number_L(S->lat, c0_bkw),
        node_frw = QDP_node_number_L(S->lat, c0_frw);

    ft_len  = subvol_ls[ft_dir];
    c0_ft   = subvol_lo[ft_dir];
    assert(ft_len * dir_mul == subvol);
    assert(QDP_sites_on_node_L(S->lat) == subvol);
    
    assert(0 == (ls[ft_dir] % ft_len)); /* latsize must be divisible by the mesh */
    int mesh_ft = ls[ft_dir] / ft_len;
    
    /* verify that all subvolumes are of the same size;
       XXX how does QMP handle mismatch in send/recv buf size? */
    for (int mu = 0 ; mu < ndim ; mu++) {
        /* QMP has no functions for global min/max for int or arrays */
        float subvol_ls_glmax = subvol_ls[mu];
        QMP_max_float(&subvol_ls_glmax);
        assert((int)subvol_ls_glmax == subvol_ls[mu]);

        float subvol_ls_glmin = subvol_ls[mu];
        QMP_min_float(&subvol_ls_glmin);
        assert((int)subvol_ls_glmin == subvol_ls[mu]);
    }

    /* init buffers for fourier transf result res_ft, mat_ft, comm */
    /* XXX send/recv buffers are of the same size */
    int ortho_len   = dir_mul * rec_len; /* length "orthogonal" to ft dir */
    int msgmem_size = ortho_len * ft_len * sizeof(double complex);
    /* FIXME switch (double complex) -> GSL_COMPLEX */
    double complex 
        *mat_ft = qlua_malloc(L, ft_len * ft_len * sizeof(double complex)),
        *res_ft = qlua_malloc(L, msgmem_size);
    QMP_mem_t 
        *qmp_send_b = QMP_allocate_memory(msgmem_size),
        *qmp_recv_b = QMP_allocate_memory(msgmem_size);
    if (NULL == res_ft 
            || NULL == qmp_send_b
            || NULL == qmp_recv_b
            || NULL == mat_ft) {
        status = "not enough memory";
        goto clearerr_0;
    }
    double complex
        *send_b = QMP_get_memory_pointer(qmp_send_b),
        *recv_b = QMP_get_memory_pointer(qmp_recv_b);
    assert(NULL != send_b);
    assert(NULL != recv_b);
    
#if 0 
    /**/printf("ortho_len=%d ft_len=%d rec_len=%d\n", ortho_len, ft_len, rec_len);
    /**/for (int mu = 0; mu < ndim ; mu++) {
    /**/    printf("%4d\t%7d\t%7d\t%7d\t%7d\t%7d\n", 
    /**/           mu, ls[mu], 
    /**/           subvol_lo[mu], subvol_hi[mu], 
    /**/           subvol_ls[mu], c_ind_mul[mu]);
    /**/}
#endif

    /* init fourier transform matrix 
       mat_ft[m,n] = exp(2*pi*i *(kLo[ft_dir]+m) *(xLo[ft_dir]+n)/L[ft_dir]) */
    for (int m = 0 ; m < ft_len ; m++)
        for (int n = 0 ; n < ft_len ; n++) {
            double ph = 2 * M_PI * (c0_ft + m) * (c0_ft + n) /(double)(ls[ft_dir]);
            mat_ft[m * ft_len + n] = cos(ph) + I*sin(ph)*ft_sign;
        }
    /* init gsl matrices to use gsl BLAS */
    gsl_matrix_complex_view gsl_res_ft = gsl_matrix_complex_view_array(
                (double *)res_ft, ortho_len, ft_len);
    gsl_matrix_complex_set_zero(&gsl_res_ft.matrix);
    gsl_matrix_complex_view gsl_mat_ft = gsl_matrix_complex_view_array(
                (double *)mat_ft, ft_len, ft_len);
    gsl_matrix_complex_view gsl_send_b = gsl_matrix_complex_view_array(
                (double *)send_b, ortho_len, ft_len);
    

    /* transform local orig.field, x -> sendbuf 
            src_co[ cInd(x[0],...,^x[ft_dir],...,x[DIM-1]), n] 
                = x[x[0], ..., xLo[ft_dir]+n, ..., x[DIM-1]] */
    for (int ind = 0 ; ind < subvol ; ind++) {
        int c[QLUA_MAX_LATTICE_RANK];
        QDP_get_coords_L(S->lat, c, QDP_this_node, ind);
        int c_ft  = c[ft_dir] - c0_ft;
        int c_ind = 0;
        for (int mu = 0 ; mu < ndim ; mu++)
            c_ind += c_ind_mul[mu] * (c[mu] - subvol_lo[mu]);
        for (int i = 0 ; i < rec_len ; i++) {
            QLA_D_Complex pz = qla_x[ind * rec_len + i]; /* XXX ind != c_ind */
            send_b[(c_ind * rec_len + i) * ft_len + c_ft] = QLA_real(pz) + I*QLA_imag(pz);
        }
    }
    
    /* iterate over i=0..MESH[ft_dir]-1 : 
       1) if not the last iter, start ASYNC send(sendbuf, neighbor_up[ft_dir]) 
            & receive(recvbuf, neighbor_down[ft_dir]) 
       2) add res_ft += dot(src_co, mat_ft^T)
       3) if not the last iter, 
          * scale mat_ft[m,n] *= exp(2*pi*i *(kLo[ft_dir]+m) *subvol_ls[ft_dir] / L[ft_dir])
          * wait for recv&send to finish; 
          * copy recv->send
     */ 

    QMP_msgmem_t qmp_msg_sendrecv[2];
    qmp_msg_sendrecv[0] = QMP_declare_msgmem(send_b, msgmem_size);
    qmp_msg_sendrecv[1] = QMP_declare_msgmem(recv_b, msgmem_size);
    QMP_msghandle_t qmp_msg[2];
    QMP_status_t qmp_st;
    for (int imesh = 0 ; imesh < mesh_ft ; imesh++) {
        if (imesh < mesh_ft - 1) { /* not the last iter */
            /*qmp_msg[0] = QMP_declare_send_relative(qmp_msg_sendrecv[0], ft_dir, +1, 0);*/
            /*qmp_msg[1] = QMP_declare_receive_relative(qmp_msg_sendrecv[1], ft_dir, -1, 0);*/
            qmp_msg[0] = QMP_declare_send_to(qmp_msg_sendrecv[0], node_frw, 0);
            qmp_msg[1] = QMP_declare_receive_from(qmp_msg_sendrecv[1], node_bkw, 0);
            if (QMP_SUCCESS != (qmp_st = QMP_start(qmp_msg[0]))) {
                status = QMP_error_string(qmp_st);
                goto clearerr_1;
            }
            if (QMP_SUCCESS != (qmp_st = QMP_start(qmp_msg[1]))) {
                status = QMP_error_string(qmp_st);
                goto clearerr_1;
            }
        }

        /* add FT(send_b) to res_ft */
        gsl_complex z1;
        GSL_SET_COMPLEX(&z1, 1.0, 0.0);
        gsl_blas_zgemm(CblasNoTrans, CblasTrans, 
                       z1, &gsl_send_b.matrix, &gsl_mat_ft.matrix,
                       z1, &gsl_res_ft.matrix);

        if (imesh < mesh_ft - 1) { /* not the last iter */
            /* multiply the ft matrix */
            for (int m = 0 ; m < ft_len ; m++) {
                gsl_vector_complex_view gsl_mat_ft_m = gsl_matrix_complex_row(
                            &gsl_mat_ft.matrix, m);
                gsl_complex z = gsl_complex_polar(1.,
                            2 * M_PI * (c0_ft + m) * ft_len / (double)(ls[ft_dir]) 
                            * ft_sign);
                gsl_vector_complex_scale(&gsl_mat_ft_m.vector, z);
            }


            /* finalize comm */
            if (QMP_SUCCESS != (qmp_st = QMP_wait_all(qmp_msg, 2))) {
                status = QMP_error_string(qmp_st);
                goto clearerr_1;
            }
            QMP_free_msghandle(qmp_msg[0]);
            QMP_free_msghandle(qmp_msg[1]);
            /* copy send_b <- recv_b for next iter */
            memcpy(send_b, recv_b, msgmem_size);
        }
    }

    /* transform local FT'ed field, res_ft -> y 
            src_co[ cInd(x[0],...,^x[ft_dir],...,x[DIM-1]), n] 
                = x[x[0], ..., xLo[ft_dir]+n, ..., x[DIM-1]] */
    for (int ind = 0 ; ind < subvol ; ind++) {
        int c[QLUA_MAX_LATTICE_RANK];
        QDP_get_coords_L(S->lat, c, QDP_this_node, ind);
        int c_ft  = c[ft_dir] - c0_ft;
        int c_ind = 0;
        for (int mu = 0 ; mu < ndim ; mu++)
            c_ind += c_ind_mul[mu] * (c[mu] - subvol_lo[mu]);
        for (int i = 0 ; i < rec_len ; i++) {
            double complex z = res_ft[(c_ind * rec_len + i) * ft_len + c_ft];
            QLA_c_eq_r_plus_ir(qla_y[ind * rec_len + i], creal(z), cimag(z));
        }
    }
    

    
    /* TODO free all buffers */
    QMP_free_msgmem(qmp_msg_sendrecv[0]);
    QMP_free_msgmem(qmp_msg_sendrecv[1]);
    QMP_free_memory(qmp_send_b);
    QMP_free_memory(qmp_recv_b);
    qlua_free(L, res_ft);
    qlua_free(L, mat_ft);

    return NULL;

clearerr_1:
    QMP_free_msghandle(qmp_msg[0]);
    QMP_free_msghandle(qmp_msg[1]);
    QMP_free_msgmem(qmp_msg_sendrecv[0]);
    QMP_free_msgmem(qmp_msg_sendrecv[1]);

clearerr_0:
    if (NULL != qmp_send_b) QMP_free_memory(qmp_send_b);
    if (NULL != qmp_recv_b) QMP_free_memory(qmp_recv_b);
    if (NULL != res_ft) qlua_free(L, res_ft);
    if (NULL != mat_ft) qlua_free(L, mat_ft);
    return status;

}

    
const char *
extra_lat_fourier_qla_full(lua_State *L,
        mLattice *S, QLA_D_Complex *qla_y, QLA_D_Complex *qla_x, 
        int rec_len, int ft_sign)
{
    const char *status = NULL;
    int ndim = QDP_ndim_L(S->lat);
    int subvol = QDP_sites_on_node_L(S->lat);
    int buf_size = subvol * rec_len * sizeof(QLA_D_Complex);
    QLA_D_Complex *qla_aux = qlua_malloc(L, buf_size);
    if (NULL == qla_aux) {
        status = "cannot allocate memory";
        goto clearerr_0;
    }
    QLA_D_Complex *p1, *p2; 
    /* before 1 step : p1 <- lat_c_x ; after last : result in p1=lat_c_y */
    if (0 == ndim % 2) {
        p1 = qla_y;
        p2 = qla_aux ; 
    } else {
        p1 = qla_aux;
        p2 = qla_y ; 
    }
    memcpy(p1, qla_x, buf_size);

    for (int mu = 0 ; mu < ndim ; mu++) {
        if (NULL != (status = extra_lat_fourier_qla_onedim(L, S, p2, p1,
                    rec_len, ft_sign, mu)))
            goto clearerr_1;
        QLA_D_Complex *p_aux = p1 ; p1 = p2 ; p2 = p_aux;
    }
    
    qlua_free(L, qla_aux);
    return NULL;

clearerr_1:
    qlua_free(L, qla_aux);
clearerr_0:
    return status;
}

const char *extra_lat_fourier_qla(lua_State *L,
                mLattice *S, QLA_D_Complex *qla_y, QLA_D_Complex *qla_x, 
                int rec_len, int ft_sign, int ft_dir)
{
    int ndim = QDP_ndim_L(S->lat);
    if (ft_dir < 0)
        return extra_lat_fourier_qla_full(L, S, qla_y, qla_x, 
                rec_len, ft_sign);
    else if (ft_dir < ndim)
        return extra_lat_fourier_qla_onedim(L, S, qla_y, qla_x, 
                rec_len, ft_sign, ft_dir);
    else 
        return "bad ft_dir";
}
