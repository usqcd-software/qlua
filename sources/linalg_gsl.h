#include <qvector.h>
#include <qmatrix.h>
#include <gsl/gsl_complex_math.h>

typedef gsl_vector dvec;
typedef gsl_vector_complex zvec;
typedef gsl_matrix dmat;

#define QLUA_INLINE inline

#define r2c(r) gsl_complex_rect(r,0)
#define vec_size(x) ((x)->size)

#define dvec_alloc(x,n) *(x) = gsl_vector_alloc(n)
#define dvec_free(x) gsl_vector_free(x)
#define dvec_get(x,i) gsl_vector_get(x,i)
#define dvec_set(x,i,v) gsl_vector_set(x,i,v)
#define dv_eq_zero(x) gsl_vector_set_zero(x)
#define dsubvec(y,x,i,n) { static gsl_vector_view _vv; _vv = gsl_vector_subvector(x,i,n); *(y) = &(_vv.vector); }
#define dv_eq_v(y,x) gsl_blas_dcopy(x,y)

#define zvec_alloc(x,n) *(x) = gsl_vector_complex_alloc(n)
#define zvec_free(x) gsl_vector_complex_free(x)
#define zvec_get(x,i) gsl_vector_complex_get(x,i)
#define zvec_set(x,i,v) gsl_vector_complex_set(x,i,v)
#define zv_eq_zero(x) gsl_vector_complex_set_zero(x)
#define zv_eq_v(y,x) gsl_blas_zcopy(x,y)
#define zv_eq_r_times_v(y,a,x) { zv_eq_v(y,x); zv_teq_r(y,a); }
#define zv_peq_r_times_v(y,a,x) gsl_blas_zaxpy(r2c(a), x, y)
#define zv_meq_r_times_v(y,a,x) gsl_blas_zaxpy(r2c(-(a)), x, y)
#define zv_teq_r(x,a) gsl_blas_zdscal(a,x)

#define dmat_alloc(x,nr,nc) { *(x) = gsl_matrix_alloc(nr,nc); }
#define dmat_free(x) gsl_matrix_free(x)
#define dmat_get(x,i,j) gsl_matrix_get(x,i,j)
#define dmat_set(x,i,j,v) gsl_matrix_set(x,i,j,v)


QLUA_INLINE double
norm2_zv(zvec *x)
{
  double lnrm2;
  lnrm2 = gsl_blas_dznrm2(x);
  lnrm2 *= lnrm2;
  QMP_sum_double(&lnrm2);
  return lnrm2;
}


typedef void (*linopnh_t)(zvec *out, zvec *in, int adj, void *args);

extern void svd_lanczos(linopnh_t linop, void *args, zvec *in, mVecReal *sv,
			int nv, zvec *qv[], int nva, zvec *qva[],
			int na, double rsq, int kmax);
