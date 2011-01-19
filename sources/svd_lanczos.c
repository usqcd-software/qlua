#include "qlua.h"                                                    /* DEPS */
#include "linalg_gsl.h"                                              /* DEPS */
#include <qmp.h>
#include <float.h>
#include <math.h>
//#include <complex.h>

extern int gsl_linalg_SV_decomp_bidiag (gsl_vector * S, gsl_vector * E);

static int nev;

#undef TRACE
#define TRACE
//#define TRACE printf0("%s %s %i\n", __FILE__, __func__, __LINE__)
//#define printf0
#define printf0 if(QDP_this_node==0) printf

#if 0
static void
svd_bi(dvec *e, dmat *m, dmat *ma, dvec *a, dvec *b, int k, int n, int na)
{
  TRACE;
  gsl_matrix *u = gsl_matrix_calloc(k,k);
  gsl_matrix *v = gsl_matrix_alloc(k,k);
  dvec *sv;
  dvec_alloc(&sv, k);

  for(int i=0; i<k-1; i++) {
    //printf0("%g\t%g\n", dvec_get(a,i), dvec_get(b,i));
    gsl_matrix_set(u, i, i+1, dvec_get(b,i));
    gsl_matrix_set(u, i, i, dvec_get(a,i));
  }
  //printf0("%g\n", dvec_get(a,k-1));
  gsl_matrix_set(u, k-1, k-1, dvec_get(a,k-1));

  dvec *work; dvec_alloc(&work, k);
  gsl_linalg_SV_decomp(u, v, sv, work);
  dvec_free(work);
  //gsl_linalg_SV_decomp_jacobi(u, v, sv);

  for(int i=0; i<k; i++) {
    //printf0("sv[%i] = %g\n", i, dvec_get(sv,i));
    dvec_set(e,i, dvec_get(sv,k-1-i));
  }
  dvec_free(sv);

#if 0
  printf0("sv[%i] %g ", 0, dvec_get(e, 0));
  printf0("sv[%i] %g ", 1, dvec_get(e, 1));
  printf0("sv[%i] %g ", nev-1, dvec_get(e, nev-1));
  printf0("sv[%i] %g\n", k-1, dvec_get(e, k-1));
#endif

  for(int i=0; i<n; i++) {
    gsl_vector_view vv;
    vv = gsl_matrix_column(v, k-1-i);
    gsl_matrix_set_col(m, i, &(vv.vector));
    // m(:,i) = v(:,k-1-i);
  }
  gsl_matrix_free(v);

  for(int i=0; i<na; i++) {
    gsl_vector_view vv;
    vv = gsl_matrix_column(u, k-1-i);
    gsl_matrix_set_col(ma, i, &(vv.vector));
    // ma(:,i) = u(:,k-1-i);
  }
  gsl_matrix_free(u);

#if 0
  if(n) {
    for(int i=0; i<k; i++) {
      for(int j=0; j<n; j++) {
	printf0(" %g", gsl_matrix_get(m, i, j));
      }
      printf0("\n");
    }
  }
#endif

  TRACE;
}
#else
#define A(i) dvec_get(a,i)
#define B(i) dvec_get(b,i)
#define V(i) dvec_get(vv,i)
#define U(i) dvec_get(uv,i)
static void
get_svec(dvec *vv, dvec *uv, dvec *a, dvec *b, double ev)
{
  double u, v, us2, vs2;
  int i, n;
  n = vec_size(vv);

  v = 1;
  u = ev/A(0);
  dvec_set(vv, 0, v);
  dvec_set(uv, 0, u);
  if(n>1) {
    vs2 = v*v;
    us2 = u*u;
    for(i=1; i<n; i++) {
      v = (ev*u-A(i-1)*v)/B(i-1);
      u = (ev*v-B(i-1)*u)/A(i);
      dvec_set(vv, i, v);
      dvec_set(uv, i, u);
      vs2 += v*v;
      us2 += u*u;
    }
    vs2 = 1/sqrt(vs2);
    us2 = 1/sqrt(us2);
    //printf0("%g %g\n", vs2, us2);
    for(i=0; i<n; i++) {
      dvec_set(vv, i, vs2*V(i));
      dvec_set(uv, i, us2*U(i));
    }
  }
}
#undef A
#undef B
#undef V
#undef U

static void
svd_bi(dvec *e, dmat *m, dmat *ma, dvec *a, dvec *b, int k, int n, int na)
{
  TRACE;
  gsl_vector_view d;
  d = gsl_vector_subvector(e, 0, k);
  for(int i=0; i<k; i++) dvec_set(&d.vector, i, dvec_get(a,i));
  dvec *t;
  dvec_alloc(&t, k-1);
  for(int i=0; i<k-1; i++) dvec_set(t, i, dvec_get(b,i));

  gsl_linalg_SV_decomp_bidiag(&d.vector, t);

  dvec_free(t);

#if 0
  printf0("sv[%i] %g ", 0, dvec_get(e, 0));
  printf0("sv[%i] %g ", 1, dvec_get(e, 1));
  printf0("sv[%i] %g ", nev-1, dvec_get(e, nev-1));
  printf0("sv[%i] %g\n", k-1, dvec_get(e, k-1));
#endif

  int nn = n<na ? na : n;
  if(nn>0) {
    //printf("%i\n", nn);
    gsl_vector_view vv, vva;
    dvec *v, *va;
    for(int i=0; i<nn; i++) {
      if(i<n) { vv = gsl_matrix_column(m, i); v = &vv.vector; }
      else if(i==n) dvec_alloc(&v, k);
      if(i<na) { vva = gsl_matrix_column(ma, i); va = &vva.vector; }
      else if(i==na) dvec_alloc(&va, k);
      get_svec(v, va, a, b, dvec_get(e,i));
    }
    if(nn>n) dvec_free(v);
    if(nn>na) dvec_free(va);
  }

#if 0
  if(n) {
    for(int i=0; i<k; i++) {
      for(int j=0; j<n; j++) {
	printf0(" %g", gsl_matrix_get(m, i, j));
      }
      printf0("\n");
    }
  }
#endif

  TRACE;
}
#endif

void
svd_lanczos(linopnh_t linop,
	    void *args,
	    zvec *in,
	    mVecReal *sv,
	    int nv,
	    zvec *qv[],
	    int nva,
	    zvec *qva[],
	    int na,
	    double rsq,
	    int kmax)
{
  double alpha, beta;
  zvec *r, *p, *u, *v;
  int i, k, kcheck;
  double dtime1, dtime2, dtime3;
  int n = in->size;

  nev = nv;
  kcheck = 1;

  TRACE;
  //printf0("kmax = %i\n", kmax);

  dvec *a, *b, *e;
  dvec_alloc(&a, kmax);
  dvec_alloc(&b, kmax);
  dvec_alloc(&e, kmax);

  TRACE;

  //printf0("n %i  na %i\n", n, na);
  zvec_alloc(&r, na);
  zvec_alloc(&u, na);
  zvec_alloc(&p, n);
  zvec_alloc(&v, n);

  dtime1 = -QMP_time();

  beta = sqrt(norm2_zv(in));
  //printf0("beta = %g\n", beta);
  zv_eq_r_times_v(v, 1/beta, in);
  zv_eq_v(p, v);
  beta = 1.;
  zv_eq_zero(u);

  k = 0;
  do {

    TRACE;
    zv_eq_r_times_v(v, 1/beta, p);
    linop(r, v, 0, args);
    zv_meq_r_times_v(r, beta, u);
    alpha = sqrt(norm2_zv(r));
    dvec_set(a, k, alpha);
    k++;
    TRACE;

#if 0
    //check singular values
    if(k>=kcheck || k>=kmax) {
      svd_bi(e, NULL, NULL, a, b, k, 0, 0);
      printf0("%i ", k);
      printf0(" sv[%i] %-9g", 0, dvec_get(e,0));
      printf0(" sv[%i] %-9g", 1, dvec_get(e,1));
      printf0(" sv[%i] %-9g", nv-1, dvec_get(e,nv-1));
      printf0(" sv[%i] %-9g\n", k-1, dvec_get(e,k-1));
      kcheck = 1 + 1.5*kcheck;
      if(k>=kmax) break;
    }
#else
    if(k>=kmax) break;
#endif

    TRACE;
    zv_eq_r_times_v(u, 1/alpha, r);
    linop(p, u, 1, args);
    zv_meq_r_times_v(p, alpha, v);
    beta = sqrt(norm2_zv(p));
    dvec_set(b, k-1, beta);
    TRACE;
  } while(1);

  dtime1 += QMP_time();
  printf0("%s %g secs\n", __func__, dtime1);
  dtime2 = -QMP_time();

  for(i=k; i<kmax; i++) {
    sv->val[i] = FLT_MAX;
  }
  kmax = k;
  if(nv>kmax) nv = kmax;
  if(nva>kmax) nva = kmax;

  dmat *vr, *ur;
  dmat_alloc(&vr, kmax, nv);
  dmat_alloc(&ur, kmax, nva);
  svd_bi(e, vr, ur, a, b, kmax, nv, nva);
  printf0("%i ", kmax);
  printf0(" sv[%i] %-9g", 0, dvec_get(e,0));
  printf0(" sv[%i] %-9g", 1, dvec_get(e,1));
  if(nv>0) printf0(" sv[%i] %-9g", nv-1, dvec_get(e,nv-1));
  printf0(" sv[%i] %-9g\n", k-1, dvec_get(e,kmax-1));
  for(i=0; i<kmax; i++) {
    sv->val[i] = dvec_get(e, i);
  }

  TRACE;
  dvec_free(a);
  dvec_free(b);
  dvec_free(e);
  TRACE;

  dtime2 += QMP_time();
  printf0("%s %g secs %g\n", __func__, dtime2, dtime1+dtime2);
  dtime3 = -QMP_time();

  for(i=0; i<nv; i++) {
    zv_eq_zero(qv[i]);
  }
  for(i=0; i<nva; i++) {
    zv_eq_zero(qva[i]);
  }

  beta = sqrt(norm2_zv(in));
  zv_eq_r_times_v(v, 1/beta, in);
  zv_eq_v(p, v);
  beta = 1.;
  zv_eq_zero(u);

  k = 0;
  do {

    //{ gsl_complex gz; gsl_blas_zdotc(v,p,&gz); printf0("v.p = %g\n", GSL_REAL(gz)); }
    zv_eq_r_times_v(v, 1/beta, p);

    //{ double nrm = norm2_zv(v); printf0("|v| = %g\n", nrm); }
    for(i=0; i<nv; i++) {
      zv_peq_r_times_v(qv[i], dmat_get(vr,k,i), v);
      //{ double nrm = norm2_zv(qv[i]); printf0("|qv[%i]| = %g\n", i, nrm); }
    }

    linop(r, v, 0, args);
    zv_meq_r_times_v(r, beta, u);
    alpha = sqrt(norm2_zv(r));
    zv_eq_r_times_v(u, 1/alpha, r);

    for(i=0; i<nva; i++) {
      zv_peq_r_times_v(qva[i], dmat_get(ur,k,i), u);
    }
    k++;
    if(k>=kmax) break;

    linop(p, u, 1, args);
    zv_meq_r_times_v(p, alpha, v);
    beta = sqrt(norm2_zv(p));
  } while(1);

  TRACE;
  zvec_free(r);
  zvec_free(p);
  zvec_free(u);
  zvec_free(v);
  TRACE;
  dmat_free(vr);
  dmat_free(ur);
  TRACE;

  dtime3 += QMP_time();
  printf0("%s %g secs %g\n", __func__, dtime3, dtime1+dtime2+dtime3);
}
