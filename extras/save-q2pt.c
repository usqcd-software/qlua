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
#include "laph_common.h"


typedef struct {
    h5output    h5o;

    int lt;
    int n_vec;

    int         buf_rank;
    hsize_t     buf_dims[4];
    hid_t       buf_dspace;
} q2pt_h5output;

/*static*/ const char *
q2pt_h5_open_write(lua_State *L,
            q2pt_h5output *q2pt_h5o,
            const char *h5file, const char *h5path,
            int lt, int n_vec)
{
    if (! is_masternode())
        return NULL;

    const char *err_str = NULL;

    hsize_t dims[7] = { lt, lt, n_vec, n_vec, NSPIN, NSPIN, 2};
    err_str = h5_open_write(L, &(q2pt_h5o->h5o), h5file, h5path, 7, dims);
    if (NULL != err_str)
        return err_str;

    q2pt_h5o->lt            = lt;
    q2pt_h5o->n_vec         = n_vec;

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
    if (! is_masternode())
        return NULL;

    if (H5Sclose(q2pt_h5o->buf_dspace))
        return "cannot close HDF5 membuf";
    q2pt_h5o->buf_dspace = -1;

    return h5_close(L, &(q2pt_h5o->h5o));
}

/* write a portion of data 
   buf      complex array [t1][n1][s1] of 
            q2pt[t1][t2][n1][n2][s1][s2][real/imag?0:1]
   assume t1, n1, s1 = sink
          t2, n2, s2 = source
         
 */
/*static*/ const char *
q2pt_h5_write(q2pt_h5output *q2pt_h5o, 
              int t_src, int j_vec_prop, int j_spin_prop,
              const gsl_complex *buf)
{
    if (! is_masternode())
        return NULL;

    hsize_t dset_off[7] = { 0, 
                            t_src, 
                            0,
                            j_vec_prop,
                            0,
                            j_spin_prop,
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
static const char *
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
            || t_axis < 0 || LDIM <= t_axis) 
        return "not implemented for this dim or t-axis";
    
    int latsize[LDIM];
    QDP_latsize_L(S->lat, latsize);
    int lt = latsize[t_axis];

    int vol4_local = QDP_sites_on_node_L(S->lat);

    /* HDF5 output */
    q2pt_h5output q2pt_h5o;
    if (NULL != (err_str = q2pt_h5_open_write(L, &q2pt_h5o, 
                    h5_file, h5_path, lt, n_vec))) {
        goto clearerr_0;
    }

    /* expose fields */
    QLA_D3_ColorVector **qla_v = NULL;
    qla_v = qlua_malloc(L, sizeof(qla_v[0]) * n_vec);
    assert(qla_v != NULL);
    for (int i = 0 ; i < n_vec ; i++)
        qla_v[i] = QDP_D3_expose_V(v[i]);

    QLA_D3_DiracFermion *qla_sol;
    qla_sol = QDP_D3_expose_D(sol);

    /* buffer for q2pt */
    gsl_complex *q2pt_buf = NULL;
    q2pt_buf = qlua_malloc(L, sizeof(q2pt_buf[0]) * lt * n_vec * NSPIN);
    assert(NULL != q2pt_buf);

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
                for (int i_c = 1 ; i_c < NCOLOR ; i_c++) 
                    QLA_c_peq_c_times_ca(res, QLA_elem_D(qla_sol[i_qdp], i_c, s), 
                            QLA_elem_V(qla_v[i_v][i_qdp], i_c));
                q2pt(t, i_v, s) = gsl_complex_add(q2pt(t, i_v, s),
                                            gsl_complex_rect(QLA_real(res), QLA_imag(res)));
            }
        }
    }
    QMP_sum_double_array((double *)q2pt_buf, 2 * lt * n_vec * NSPIN);
    if (NULL != (err_str = q2pt_h5_write(&q2pt_h5o, 
                    t_src, j_vec_prop, j_spin_prop, 
                    q2pt_buf)))
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
   qcd.save_q2pt(h5_file, h5_path, 
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
    QDP_D3_DiracFermion *sol= qlua_checkLatDirFerm3(L, 6, S, NCOLOR)->ptr;
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

#if 1
    
/*static*/ const char *
save_q2pt_list(lua_State *L, 
        mLattice *S,
        const char *h5_file,
        const char *h5_path,
        int n_sol, const int *tsrc, const int *jvec, const int *jspin, 
        QDP_D3_DiracFermion **sol,
        int n_vec, QDP_D3_ColorVector **v,
        int t_axis)
/* calc correlator  < q^{ivec,ispin}(tsnk) \bar{q}^{jvec,jspin}(tsrc) >
   ivec,jvec        sink, source LapH vector index
   tsnk,tsrc        sink, source time slices
   ispin,jspin      sink, source spin index   
   x3               site spatial index
   c                color

   IN:
   v[ivec][iqdp][c] 
        --> [for all tsnk] v_matr[ivec, {x3, c}]
   sol[Jsol={tsrc, jvec, jspin}][iqdp][c,ispin]  
        --> [for all tsnk] sol_matr[Jsol, ispin, {x3, c}]

   for each tsnk: calc matrix product over {x3, c}:
     Q[Jsol, ispin, ivec] := 
            v_matr^\dag[ivec, {x3, c}] . sol_matr[tsnk][Jsol, ispin, {x3, c}]
        
   OUT to hdf5:
   cycle over Jsol [tsnk, ivec, ispin]
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
#define i_vol3(c) ((c)[vol3_axis[0]] - (node_x0)[vol3_axis[0]] + vol3_dim[0]*(\
                   (c)[vol3_axis[1]] - (node_x0)[vol3_axis[1]] + vol3_dim[1]*(\
                   (c)[vol3_axis[2]] - (node_x0)[vol3_axis[2]])))
#define i_X_vol3(X,c) (i_vol3(c) + vol3_local*(X))
#define i_vol4(c) i_X_vol3(((c)[t_axis] - node_t0), c)
#else
#error "LDIM!=4 is not consistent"
#endif

    /* make map i_qdp -> tsnk, i_x3 */
    int *tx2iqdp = NULL;
    tx2iqdp     = qlua_malloc(L, sizeof(tx2iqdp[0]) * vol4_local);
    assert(NULL != tx2iqdp);
#define Tx2iqdp(t,ix) ((tx2iqdp)[ix + (vol3_local)*((t) - (node_t0))])
    for (int iqdp = 0; iqdp < vol4_local ; iqdp++) {
        int coord[LDIM];
        QDP_get_coords_L(S->lat, coord, QDP_this_node, iqdp);
        /* FIXME change order of ix3 for fewer cache misses 
            when converting data [iqdp] <-> [t, i_x3] ? */
        int t = coord[t_axis], 
            ix = i_vol3(coord);
        Tx2iqdp(t, ix) = iqdp;
    }   

    /* HDF5 output */
    q2pt_h5output q2pt_h5o;
    if (NULL != (err_str = q2pt_h5_open_write(L, &q2pt_h5o, 
                    h5_file, h5_path, lt, n_vec))) {
        goto clearerr_0;
    }

    /* v[ivec][iqdp][c] 
        --> [for all tsnk] v_matr[ivec, {x3, c}] */
    gsl_complex *v_matr = NULL;
    v_matr      = qlua_malloc(L, 
            sizeof(v_matr[0]) * n_vec * vol3_local * NCOLOR);
    assert(NULL != v_matr);
#define Vmatr(ivec, ix3, c) ((v_matr)[(c) + (NCOLOR)*(\
                                      (ix3) + (vol3_local)*\
                                      (ivec))])
    gsl_matrix_complex_view gsl_v_matr = gsl_matrix_complex_view_array(
            (double *)v_matr, n_vec, vol3_local * NCOLOR);

    /* sol[Jsol={tsrc, jvec, jspin}][iqdp][c,ispin]  
         -->[for all tsnk] sol_matr[Jsol, ispin, {x3, c}] */
    gsl_complex *sol_matr = NULL;
    sol_matr    = qlua_malloc(L, 
            sizeof(sol_matr[0]) * n_sol * vol3_local * NCOLOR * NSPIN);
    assert(NULL != sol_matr);
#define Smatr(Jsol, ispin, ix3, c) ((sol_matr)[(c) + (NCOLOR)*(\
                                               (ix3) + (vol3_local)*(\
                                               (ispin) + (NSPIN)*\
                                               (Jsol)))])
    gsl_matrix_complex_view gsl_sol_matr = gsl_matrix_complex_view_array(
            (double *)sol_matr, n_sol * NSPIN, vol3_local * NCOLOR);

    /* prod_matr[Jsol, ispin, ivec] -- full time extent for global sum */
    gsl_complex *prod_matr = NULL;
    prod_matr   = qlua_malloc(L, 
            sizeof(prod_matr[0]) * n_sol * NSPIN * n_vec);
    assert(NULL != prod_matr);
#define Pmatr(Jsol, ispin, ivec) ((prod_matr)[(ivec) + (n_vec)*(\
                                              (ispin) + (NSPIN)*\
                                              (Jsol))])
    gsl_matrix_complex_view gsl_prod_matr = gsl_matrix_complex_view_array(
            (double *)prod_matr, n_sol * NSPIN, n_vec);

    /* output buf */
    int buf_out_len = n_sol * lt * n_vec * NSPIN;
    gsl_complex *buf_out = NULL;
    buf_out     = qlua_malloc(L, sizeof(prod_matr[0]) * buf_out_len);
    assert(NULL != buf_out);
    memset(buf_out, 0, sizeof(prod_matr[0]) * buf_out_len);
#define Obuf(Jsol, tsnk, ivec, ispin) ((buf_out)[(ispin) + (NSPIN)*(\
                                                 (ivec) + (n_vec)*(\
                                                 (tsnk) + (lt)*\
                                                 (Jsol)))])

    /* expose fields */
    QLA_D3_ColorVector **qla_v = NULL;
    qla_v = qlua_malloc(L, sizeof(qla_v[0]) *n_vec);
    assert(NULL != qla_v);
    for (int i = 0 ; i < n_vec ; i++)
        qla_v[i] = QDP_D3_expose_V(v[i]);

    QLA_D3_DiracFermion **qla_sol = NULL;
    qla_sol = qlua_malloc(L, sizeof(qla_sol[0]) *n_sol);
    assert(NULL != qla_sol);
    for (int i = 0 ; i < n_sol ; i++)
        qla_sol[i] = QDP_D3_expose_D(sol[i]);
    
    /* cycle through tsnk */
    for (int tsnk = node_t0; tsnk < node_t0 + node_lt ; tsnk++) {
        /* copy timeslice fields */
        for (int ix3 = 0 ; ix3 < vol3_local ; ix3++) {
            for (int ivec = 0 ; ivec < n_vec ; ivec++) 
                for (int c = 0 ; c < NCOLOR ; c++) {
                    QLA_D_Complex x = QLA_D3_elem_V(
                            qla_v[ivec][Tx2iqdp(tsnk, ix3)], c);
                    Vmatr(ivec, ix3, c) = gsl_complex_rect(
                            QLA_real(x), QLA_imag(x));
                }
            for (int isol = 0 ; isol < n_sol ; isol ++)
                for (int c = 0 ; c < NCOLOR ; c++)
                    for (int s = 0 ; s < NSPIN ; s++) {
                        QLA_D_Complex x = QLA_D3_elem_D(
                                qla_sol[isol][Tx2iqdp(tsnk, ix3)], c, s);
                        Smatr(isol, s, ix3, c) = gsl_complex_rect(
                                QLA_real(x), QLA_imag(x));
                    }
        }
        /* The Heart */
        gsl_blas_zgemm(CblasNoTrans, CblasConjTrans, gsl_complex_rect(1,0),
                       &gsl_sol_matr.matrix, &gsl_v_matr.matrix,
                       gsl_complex_rect(0,0), &gsl_prod_matr.matrix);
        /* transpose for output */
        for (int isol = 0 ; isol < n_sol ; isol++)
            for (int ispin = 0 ; ispin < NSPIN ; ispin++)
                for (int ivec = 0 ; ivec < n_vec ; ivec++) {
                    Obuf(isol, tsnk, ivec, ispin) = Pmatr(isol, ispin, ivec);
                }
    }

    QMP_sum_double_array((double *)buf_out, 2 * buf_out_len);
    for (int isol = 0 ; isol < n_sol ; isol++) {
        if (NULL != (err_str = q2pt_h5_write(&q2pt_h5o,
                        tsrc[isol], jvec[isol], jspin[isol],
                        &Obuf(isol, 0, 0, 0))))
            goto clearerr_1;
    }
        

    /* cleanup, reset fields */
clearerr_1:
    qlua_free(L, v_matr);
    qlua_free(L, sol_matr);
    qlua_free(L, prod_matr);
    qlua_free(L, buf_out);

    for (int ivec = 0 ; ivec < n_vec ; ivec++)
        QDP_D3_reset_V(v[ivec]);
    qlua_free(L, qla_v);

    for (int isol = 0 ; isol < n_sol ; isol++)
        QDP_D3_reset_D(sol[isol]);
    qlua_free(L, qla_sol);

clearerr_0:
    qlua_free(L, tx2iqdp);


    return err_str;

#undef Vmatr
#undef Smatr
#undef Pmatr
#undef Obuf

#undef Tx2iqdp

#undef i_vol3
#undef i_X_vol3
#undef i_vol4
}

/* QLua wrapper 
   qcd.save_q2pt_list(h5_file, h5_path, 
                      sol_list,
                      vec_list, t_axis)
    h5_file     HDF5 file name
    h5_path     keypath
    sol_list    {{tsrc, jvec, jspin, sol} ... }
        tsrc, jvec, jspin    
                source time, evec and spin
        sol     solution vector (DiracFermion)
    vec_list    list of LapH vectors (ColorVector)
    t_axis

 */
/*DEBUG*/void
print_sol_list(int n_sol, int *tsrc, int *jvec, int *jspin, QDP_D3_DiracFermion **sol)
{
    for (int i = 0 ; i < n_sol; i++)
        printf("%d\t%d\t%d\t%d\t%p\n", i, tsrc[i], jvec[i], jspin[i], sol[i]);
}


int         
q_save_q2pt_list(lua_State *L)
{
    int argc = lua_gettop(L);
    if (5 != argc) {
        luaL_error(L, "expect 5 arguments");
        return 1;
    }

    const char *h5_file     = luaL_checkstring(L, 1);
    const char *h5_path     = luaL_checkstring(L, 2);

//    int t_src               = luaL_checkint(L, 3);
//    int j_vec               = luaL_checkint(L, 4);
//    int j_spin              = luaL_checkint(L, 5);
    
    int n_vec;
    mLattice *S = NULL;
    QDP_D3_ColorVector **v  = qlua_check_latcolvec_table(L, 4, NULL, -1, &S, &n_vec);
    if (NULL == v) {
        luaL_error(L, "a list of latcolvec expected in #5");
        goto clearerr_0;
    }

    
    int t_axis              = luaL_checkint(L, 5);

    int sol_list_idx = 3;
    int n_sol;
    int *tsrc   = NULL, *jvec   = NULL, *jspin  = NULL;
    QDP_D3_DiracFermion **sol = NULL;
#if 0
    if (LUA_TTABLE != lua_type(L, sol_list_idx)) {
        luaL_error(L, "list of {tsrc, jvec, jspin, sol} objects expected");
        return 1;
    }
    n_sol = lua_objlen(L, sol_list_idx);
    tsrc    = qlua_malloc(L, sizeof(tsrc[0]) * n_sol);
    assert(NULL != tsrc);
    jvec    = qlua_malloc(L, sizeof(jvec[0]) * n_sol);
    assert(NULL != jvec);
    jspin   = qlua_malloc(L, sizeof(jspin[0]) * n_sol);
    assert(NULL != jspin);
    sol     = qlua_malloc(L, sizeof(sol[0]) * n_sol);
    assert(NULL != sol);

    for (int isol = 0 ; isol < n_sol ; isol++) {
        lua_pushnumber(L, isol + 1);
        lua_gettable(L, sol_list_idx);
        int sol_idx = lua_gettop(L);

        if (LUA_TTABLE != lua_type(L, sol_idx) || 4 != lua_objlen(L, sol_idx)){
            qlua_free(L, tsrc);
            qlua_free(L, jvec);
            qlua_free(L, jspin);
            qlua_free(L, sol);
            qlua_free(L, v);
            luaL_error(L, "list of {tsrc, jvec, jspin, sol} objects expected");
            return 1;
        }
        
        lua_pushnumber(L, 1); 
        lua_gettable(L, sol_idx);
        tsrc[isol]  = luaL_checkint(L, -1);
        lua_pop(L, 1);

        lua_pushnumber(L, 2); 
        lua_gettable(L, sol_idx);
        jvec[isol]  = luaL_checkint(L, -1);
        lua_pop(L, 1);

        lua_pushnumber(L, 3); 
        lua_gettable(L, sol_idx);
        jspin[isol]  = luaL_checkint(L, -1);
        lua_pop(L, 1);

        lua_pushnumber(L, 4); 
        lua_gettable(L, sol_idx);
        sol[isol]  = qlua_checkLatDirFerm3(L, -1, S, NCOLOR)->ptr;
        lua_pop(L, 1);

        lua_pop(L, 1);
    }
#else
    if (qlua_check_laph_sol_list(L, sol_list_idx, S, // <<< error is in thsi function!@@
                &n_sol, &tsrc, &jvec, &jspin, &sol, NULL)) {
        luaL_error(L, "list of {tsrc, jvec, jspin, sol} objects expected");
        goto clearerr_1;
    }
#endif


    const char *status = save_q2pt_list(L, S, h5_file, h5_path,
                            n_sol, tsrc, jvec, jspin, sol,
                            n_vec, v,
                            t_axis); 
    if (NULL != status) {
        luaL_error(L, status);
        goto clearerr_2;
    }

    qlua_free(L, tsrc);
    qlua_free(L, jvec);
    qlua_free(L, jspin);
    qlua_free(L, sol);
    qlua_free(L, v);

    return 0;

clearerr_2:
    qlua_free(L, tsrc);
    qlua_free(L, jvec);
    qlua_free(L, jspin);
    qlua_free(L, sol);
clearerr_1:     qlua_free(L, v);
clearerr_0:     return 1;
}

#endif



#undef LDIM
#undef NSPIN
#undef NCOLOR
