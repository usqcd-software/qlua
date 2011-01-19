#define QOP_Precision 2
#include "modules.h"                                                 /* DEPS */
#include "qlua.h"                                                    /* DEPS */
#include "qclover.h"                                                 /* DEPS */
#include "qcomplex.h"                                                /* DEPS */
#include "qvector.h"                                                 /* DEPS */
#ifdef HAS_GSL
#include "qmatrix.h"                                                 /* DEPS */
#include "linalg_gsl.h"                                              /* DEPS */
#endif
#include "lattice.h"                                                 /* DEPS */
#include "qlayout.h"                                                 /* DEPS */
#include "latcolvec.h"                                               /* DEPS */
#include "latcolmat.h"                                               /* DEPS */
#include "latdirferm.h"                                              /* DEPS */
#include "latdirprop.h"                                              /* DEPS */
#include "qop_qdp.h"
#include "qdp_f3.h"
#include "qdp_df3.h"

#define printf0 if(QDP_this_node==0) printf
#define TRACE printf0("%s\t%s\t%i\n", __FILE__, __func__, __LINE__)
#if USE_Nc3

static void
Init(QDP_Lattice *lat)
{
  static int inited=0;

  if(!QOP_is_initialized()) {
    QOP_layout_t qoplayout;
    qoplayout.latdim = QDP_ndim_L(lat);
    qoplayout.latsize = (int *) malloc(QDP_ndim_L(lat)*sizeof(int));
    qoplayout.machdim = -1;
    QDP_latsize_L(lat, qoplayout.latsize);
    QDP_set_default_lattice(lat);
    printf0("begin QOP init... ");
    QOP_init(&qoplayout);
    printf0("done\n");
  }
  inited = 1;
}

static void
add_default(lua_State *L, const char *key, double value)
{
  lua_pushnumber(L, value);
  lua_setfield(L, -2, key);
}

static double
set_default(lua_State *L, const char *key, double def)
{
  double v = def;

  lua_getfield(L, -1, key);
  if (qlua_qtype(L, -1) == qReal) {
    v = luaL_checknumber(L, -1);
  }
  lua_pop(L, 1);
  return v;
}

#include "qqopqdp-asqtad.c"
#include "qqopqdp-clover.c"
#include "qqopqdp-hisq.c"

#ifndef HAS_GSL

static int
q_orthonormalize(lua_State *L)
{
  return luaL_error(L, "need GSL for Orthonormalize");
}

static int
q_rr(lua_State *L)
{
  return luaL_error(L, "need GSL for RayleighRitz");
}

#else

void
do_orthonormalize(gsl_matrix_complex *v, int nv,
		  gsl_matrix_complex *v2, int nv2, int ne)
{
  gsl_complex zero = gsl_complex_rect(0,0);
  gsl_complex one = gsl_complex_rect(1,0);

  if(nv2>0) { // orthogonalize v against v2
    // vv[1:nv2,1:nv2] = v2[1:nv2,:]*v2'[:,1:nv2]
    // mv[1:nv,1:nv2] = v[1:nv,:]*v2'[:,1:nv2]
    gsl_matrix_complex *vv = gsl_matrix_complex_alloc(nv2,nv2);
    gsl_matrix_complex *mv = gsl_matrix_complex_alloc(nv,nv2);
    gsl_blas_zgemm(CblasNoTrans, CblasConjTrans, one, v2, v2, zero, vv);
    gsl_blas_zgemm(CblasNoTrans, CblasConjTrans, one, v, v2, zero, mv);
    QMP_sum_double_array(vv->data, 2*nv2*nv2);
    QMP_sum_double_array(mv->data, 2*nv*nv2);

    // s = -mv*(vv^-1)
    int s;
    gsl_permutation *p = gsl_permutation_alloc(nv2);
    gsl_linalg_complex_LU_decomp(vv, p, &s);
    // mv = mv * U^-1
    gsl_blas_ztrsm(CblasRight,CblasUpper,CblasNoTrans,CblasNonUnit,one,vv,mv);
    // mv = mv * L^-1
    gsl_blas_ztrsm(CblasRight,CblasLower,CblasNoTrans,CblasUnit,one,vv,mv);
    // mv = -mv * p;
    gsl_permutation *q = gsl_permutation_alloc(nv2);
    gsl_permutation_linear_to_canonical(q, p);
    for(int i=0; i<nv; i++) {
      gsl_complex z = zero;
      s = nv2;
      for(int j=0; j<nv2; j++) {
	int so = s;
	s = gsl_permutation_get(q, j);
	if(s>so) { // in cycle
	  gsl_complex t;
	  t = gsl_complex_negative(gsl_matrix_complex_get(mv, i, s));
	  gsl_matrix_complex_set(mv, i, so, t);
	} else { // new cycle
	  if(so<nv2) gsl_matrix_complex_set(mv, i, so, z);
	  z = gsl_complex_negative(gsl_matrix_complex_get(mv, i, s));
	}
      }
      gsl_matrix_complex_set(mv, i, s, z);
    }
    gsl_permutation_free(q);
    gsl_permutation_free(p);

    // v[1:nv,:] = v[1:nv,:] + c[1:nv,1:nv2]*v2[1:nv2,:]
    gsl_blas_zgemm(CblasNoTrans, CblasNoTrans, one, mv, v2, one, v);
  }

  gsl_matrix_complex *mv0 = gsl_matrix_complex_alloc(nv,1);
  for(int i=0; i<nv; i++) {
    // vv[nv-i,1] = v[j>=i,:]*v'[:,i]
    gsl_matrix_complex_view mj = gsl_matrix_complex_submatrix(v,i,0,nv-i,ne);
    gsl_matrix_complex_view mi = gsl_matrix_complex_submatrix(v,i,0,1,ne);
    gsl_matrix_complex_view mv = gsl_matrix_complex_submatrix(mv0,0,0,nv-i,1);
    gsl_blas_zgemm(CblasNoTrans, CblasConjTrans, one, &mj.matrix, &mi.matrix,
		   zero, &mv.matrix);
    QMP_sum_double_array(mv0->data, 2*(nv-i));

    // normalize v[i]
    gsl_vector_complex_view vi = gsl_matrix_complex_row(v, i);
    double s = 1.0/sqrt(GSL_REAL(gsl_matrix_complex_get(mv0,0,0)));
    gsl_blas_zdscal(s, &vi.vector);

    if(i<nv-1) {
      gsl_vector_complex_view vv = gsl_matrix_complex_subcolumn(mv0,0,1,nv-1-i);
      gsl_blas_zdscal(-s, &vv.vector);
      // v[j>i,:] = v[j,:] + c[j,1]*v[i,:]
      mj = gsl_matrix_complex_submatrix(v,i+1,0,nv-1-i,ne);
      gsl_blas_zgeru(one, &vv.vector, &vi.vector, &mj.matrix);
    }
  }

  gsl_matrix_complex_free(mv0);
}

void
do_rr(gsl_matrix_complex *v, gsl_matrix_complex *Av, int nv, int ne, int no)
{
  gsl_matrix_complex *vt = gsl_matrix_complex_alloc(nv, ne);
  gsl_matrix_complex *Avt = gsl_matrix_complex_alloc(nv, no);
  gsl_matrix_complex *vv = gsl_matrix_complex_alloc(nv, nv);
  gsl_matrix_complex *vAAv = gsl_matrix_complex_alloc(nv, nv);
  gsl_matrix_complex *evec = gsl_matrix_complex_alloc(nv, nv);
  gsl_vector *eval = gsl_vector_alloc(nv);

#if 0
  // orthogonalize
  gsl_blas_zherk(CblasLower, CblasNoTrans, 1, v, 0, vv);
  QMP_sum_double_array(vv->data, 2*nv*nv);
  {
    gsl_eigen_hermv_workspace *work = gsl_eigen_hermv_alloc(nv);
    gsl_eigen_hermv(vv, eval, evec, work);
    gsl_eigen_hermv_free(work);
  }
  gsl_blas_zgemm(CblasConjTrans, CblasNoTrans, gsl_complex_rect(1,0),
		 evec, v, gsl_complex_rect(0,0), vt);
  gsl_matrix_complex_memcpy(v, vt);
  gsl_blas_zgemm(CblasConjTrans, CblasNoTrans, gsl_complex_rect(1,0),
		 evec, Av, gsl_complex_rect(0,0), Avt);
  gsl_matrix_complex_memcpy(Av, Avt);
#endif

  // get matrices
  gsl_blas_zherk(CblasLower, CblasNoTrans, 1, v, 0, vv);
  QMP_sum_double_array(vv->data, 2*nv*nv);
#if 0
  for(int i=0; i<nv; i++) {
    for(int j=0; j<i; j++) {
      gsl_complex z = gsl_matrix_complex_get(vv,j,i);
      gsl_matrix_complex_set(vv,i,j,gsl_complex_conjugate(z));
    }
  }
#endif
  gsl_blas_zherk(CblasLower, CblasNoTrans, 1, Av, 0, vAAv);
  QMP_sum_double_array(vAAv->data, 2*nv*nv);
#if 0
  for(int i=0; i<nv; i++) {
    for(int j=0; j<i; j++) {
      gsl_complex z = gsl_matrix_complex_get(vAAv,j,i);
      gsl_matrix_complex_set(vAAv,i,j,gsl_complex_conjugate(z));
    }
  }
#endif

#if 0
  for(int i=0; i<nv; i++) {
    for(int j=0; j<nv; j++) {
      gsl_complex z = gsl_matrix_complex_get(vv,i,j);
      printf0("vv[%i,%i] = %12g %12g\n", i, j, GSL_REAL(z), GSL_IMAG(z));
    }
  }
  for(int i=0; i<nv; i++) {
    for(int j=0; j<nv; j++) {
      gsl_complex z = gsl_matrix_complex_get(vAAv,i,j);
      printf0("vAAv[%i,%i] = %12g %12g\n", i, j, GSL_REAL(z), GSL_IMAG(z));
    }
  }
#endif

  // generalized eigenvalues
  {
    gsl_eigen_genhermv_workspace *work = gsl_eigen_genhermv_alloc(nv);
    gsl_eigen_genhermv(vAAv, vv, eval, evec, work);
    gsl_eigen_genhermv_free(work);
  }

#if 0
  for(int i=0; i<nv; i++) {
    printf0("ev[%i] = %g\n", i, gsl_vector_get(eval,i));
  }
  for(int i=0; i<nv; i++) {
    for(int j=0; j<nv; j++) {
      gsl_complex z = gsl_matrix_complex_get(evec,i,j);
      printf0("evec[%i,%i] = %12g %12g\n", i, j, GSL_REAL(z), GSL_IMAG(z));
    }
  }
#endif

  // rotate vectors
  gsl_blas_zgemm(CblasConjTrans, CblasNoTrans, gsl_complex_rect(1,0), evec, v, gsl_complex_rect(0,0), vt);
  gsl_matrix_complex_memcpy(v, vt);

  gsl_blas_zgemm(CblasConjTrans, CblasNoTrans, gsl_complex_rect(1,0), evec, Av, gsl_complex_rect(0,0), Avt);
  gsl_matrix_complex_memcpy(Av, Avt);

  gsl_vector_free(eval);
  gsl_matrix_complex_free(evec);
  gsl_matrix_complex_free(vAAv);
  gsl_matrix_complex_free(vv);
  gsl_matrix_complex_free(Avt);
  gsl_matrix_complex_free(vt);
}

static int
q_orthonormalize(lua_State *L)
{
  luaL_checktype(L, 1, LUA_TTABLE);
  lua_pushnumber(L, 1);
  lua_gettable(L, 1);
  mLattice *S = qlua_ObjLattice(L, -1);
  lua_pop(L, 2);

  int nv = -1, nv2 = -1;
  mLatColVec3 **v = qlua_checkTableLatColVec3(L, 1, S, 3, &nv);
  mLatColVec3 **v2 = NULL;
  if(lua_gettop(L)>1) {
    v2 = qlua_checkTableLatColVec3(L, 2, S, 3, &nv2);
  }
  //printf0("nv2 %i ptr %p\n", nv2, v2);

  int ne = 3*QDP_subset_len(S->even);

  gsl_matrix_complex *gv = gsl_matrix_complex_alloc(nv, ne);
  gsl_matrix_complex *gv2 = NULL;
  if(nv2>0) gv2 = gsl_matrix_complex_alloc(nv2, ne);
  for(int i=0; i<nv; i++) {
    QDP_D3_extract_packed_V((QLA_D3_ColorVector *)
			    gsl_matrix_complex_ptr(gv, i, 0),
			    v[i]->ptr, S->even);
  }
  for(int i=0; i<nv2; i++) {
    QDP_D3_extract_packed_V((QLA_D3_ColorVector *)
			    gsl_matrix_complex_ptr(gv2, i, 0),
			    v2[i]->ptr, S->even);
  }

  //TRACE;
  do_orthonormalize(gv, nv, gv2, nv2, ne);
  //TRACE;

  for(int i=0; i<nv; i++) {
    QDP_D3_insert_packed_V(v[i]->ptr, (QLA_D3_ColorVector *)
			   gsl_matrix_complex_ptr(gv, i, 0), S->even);
  }

  if(nv2>0) gsl_matrix_complex_free(gv2);
  gsl_matrix_complex_free(gv);
  if(v) qlua_free(L, v2);
  qlua_free(L, v);

  return 0;
}

static int
q_rr(lua_State *L)
{
  luaL_checktype(L, 1, LUA_TTABLE);
  lua_pushnumber(L, 1);
  lua_gettable(L, 1);
  mLattice *S = qlua_ObjLattice(L, -1);
  lua_pop(L, 1);

  int nv = -1;
  mLatColVec3 **v = qlua_checkTableLatColVec3(L, 1, S, 3, &nv);
  mLatColVec3 **Av = qlua_checkTableLatColVec3(L, 2, S, 3, &nv);

  int ne = 3*QDP_subset_len(S->even);
  int no = 3*QDP_subset_len(S->odd);

  gsl_matrix_complex *gv = gsl_matrix_complex_alloc(nv, ne);
  gsl_matrix_complex *gAv = gsl_matrix_complex_alloc(nv, no);
  for(int i=0; i<nv; i++) {
    QDP_D3_extract_packed_V((QLA_D3_ColorVector *)
			    gsl_matrix_complex_ptr(gv, i, 0),
			    v[i]->ptr, S->even);
  }
  for(int i=0; i<nv; i++) {
    QDP_D3_extract_packed_V((QLA_D3_ColorVector *)
			    gsl_matrix_complex_ptr(gAv, i, 0),
			    Av[i]->ptr, S->odd);
  }

  //TRACE;
  do_rr(gv, gAv, nv, ne, no);
  //TRACE;

  for(int i=0; i<nv; i++) {
    QDP_D3_insert_packed_V(v[i]->ptr, (QLA_D3_ColorVector *)
			   gsl_matrix_complex_ptr(gv, i, 0), S->even);
  }
  for(int i=0; i<nv; i++) {
    QDP_D3_insert_packed_V(Av[i]->ptr, (QLA_D3_ColorVector *)
			   gsl_matrix_complex_ptr(gAv, i, 0), S->odd);
  }

  gsl_matrix_complex_free(gv);
  gsl_matrix_complex_free(gAv);
  qlua_free(L, Av);
  qlua_free(L, v);

  return 0;
}

#endif

static struct luaL_Reg fQOPQDP[] = {
    { "Asqtad",         q_asqtad },
    { "Clover",         q_clover },
    { "Hisq",           q_hisq },
    { "Orthonormalize", q_orthonormalize },
    { "RayleighRitz",   q_rr },
    { NULL,           NULL }
};

int
init_qopqdp(lua_State *L)
{
  lua_getglobal(L, qcdlib);
  lua_newtable(L);
  luaL_register(L, NULL, fQOPQDP);
  lua_setfield(L, -2, "qopqdp");
  lua_pop(L, 1);
  return 0;
}
#else /* USE_Nc3 */
int
init_qopqdp(lua_State *L)
{
  return 0;
}
#endif /* USE_Nc3 */

int
fini_qopqdp(lua_State *L)
{
  return 0;
}
