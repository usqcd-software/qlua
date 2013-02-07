#ifndef MARK_FE647744_2081_4E52_9AE1_899F1FA89AE5
#define MARK_FE647744_2081_4E52_9AE1_899F1FA89AE5

#define GSL_RANGE_CHECK_OFF
#define HAVE_INLINE
#include <gsl/gsl_errno.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_matrix_double.h>
#include <gsl/gsl_vector_double.h>
#include <gsl/gsl_permutation.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_complex_math.h>

typedef struct {
    int l_size;
    int r_size;
    gsl_matrix *m;
} mMatReal;

typedef struct {
    int l_size;
    int r_size;
    gsl_matrix_complex *m;
} mMatComplex;

extern const char mtnMatReal[];
extern const char mtnMatComplex[];

int init_matrix(lua_State *L);
int fini_matrix(lua_State *L);

mMatReal *qlua_checkMatReal(lua_State *L, int idx);
mMatReal *qlua_newMatReal(lua_State *L, int sl, int sr);
mMatComplex *qlua_checkMatComplex(lua_State *L, int idx);
mMatComplex *qlua_newMatComplex(lua_State *L, int sl, int sr);

#endif /* !defined(MARK_FE647744_2081_4E52_9AE1_899F1FA89AE5) */
