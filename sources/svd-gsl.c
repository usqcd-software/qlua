// modified for bidiagonal SVD

/* linalg/svd.c
 * 
 * Copyright (C) 1996, 1997, 1998, 1999, 2000, 2004, 2007, 2010 Gerard Jungman, Brian Gough
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

//#include <config.h>
#include <stdlib.h>
#include <string.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>

#include <gsl/gsl_linalg.h>

#include "givens-gsl.c"
#include "svdstep-gsl.c"

/* Factorise a general M x N matrix A into,
 *
 *   A = U D V^T
 *
 * where U is a column-orthogonal M x N matrix (U^T U = I), 
 * D is a diagonal N x N matrix, 
 * and V is an N x N orthogonal matrix (V^T V = V V^T = I)
 *
 * U is stored in the original matrix A, which has the same size
 *
 * V is stored as a separate matrix (not V^T). You must take the
 * transpose to form the product above.
 *
 * The diagonal matrix D is stored in the vector S,  D_ii = S_i
 */

int
gsl_linalg_SV_decomp_bidiag (gsl_vector * S, gsl_vector * E)
{
  size_t a, b, i, j, iter;

  const size_t N = S->size;

  /* Handle the case of N = 1 */
  if (N == 1) {
    return GSL_SUCCESS;
  }

  if (E->size != N-1) {
    GSL_ERROR ("length of E must be one less than length of D",
	       GSL_EBADLEN);
  }

  /* apply reduction steps to B=(S,Sd) */

  chop_small_elements (S, E);

  /* Progressively reduce the matrix until it is diagonal */
  b = N - 1;
  iter = 0;
  while (b > 0) {
    double fbm1 = gsl_vector_get (E, b - 1);

    if (fbm1 == 0.0 || gsl_isnan (fbm1)) {
      b--;
      continue;
    }

    /* Find the largest unreduced block (a,b) starting from b
       and working backwards */

    a = b - 1;
    while (a > 0) {
      double fam1 = gsl_vector_get (E, a - 1);
      if (fam1 == 0.0 || gsl_isnan (fam1)) break;
      a--;
    }

    iter++;

    if (iter > 100 * N) {
      GSL_ERROR("SVD decomposition failed to converge", GSL_EMAXITER);
    }

    {
      const size_t n_block = b - a + 1;
      gsl_vector_view S_block = gsl_vector_subvector (S, a, n_block);
      gsl_vector_view f_block = gsl_vector_subvector (E, a, n_block - 1);

      int rescale = 0;
      double scale = 1; 
      double norm = 0;

      /* Find the maximum absolute values of the diagonal and subdiagonal */
      for (i = 0; i < n_block; i++) {
	double s_i = gsl_vector_get (&S_block.vector, i);
	double a = fabs(s_i);
	if (a > norm) norm = a;
      }
      for (i = 0; i < n_block - 1; i++) {
	double f_i = gsl_vector_get (&f_block.vector, i);
	double a = fabs(f_i);
	if (a > norm) norm = a;
      }

      /* Temporarily scale the submatrix if necessary */
      if (norm > GSL_SQRT_DBL_MAX) {
	scale = (norm / GSL_SQRT_DBL_MAX);
	rescale = 1;
      }
      else if (norm < GSL_SQRT_DBL_MIN && norm > 0) {
	scale = (norm / GSL_SQRT_DBL_MIN);
	rescale = 1;
      }
      if (rescale) {
	gsl_blas_dscal(1.0 / scale, &S_block.vector);
	gsl_blas_dscal(1.0 / scale, &f_block.vector);
      }

      /* Perform the implicit QR step */
      qrstep (&S_block.vector, &f_block.vector);

      /* remove any small off-diagonal elements */
      chop_small_elements (&S_block.vector, &f_block.vector);

      /* Undo the scaling if needed */
      if (rescale) {
	gsl_blas_dscal(scale, &S_block.vector);
	gsl_blas_dscal(scale, &f_block.vector);
      }
    }

  }

  /* Make singular values positive by reflections if necessary */
  for (j = 0; j < N; j++) {
    double Sj = gsl_vector_get (S, j);
    if (Sj < 0.0) {
      gsl_vector_set (S, j, -Sj);
    }
  }

  /* Sort singular values into increasing order */
  for (i = 0; i < N; i++) {
    double S_max = gsl_vector_get (S, i);
    size_t i_max = i;

    for (j = i + 1; j < N; j++) {
      double Sj = gsl_vector_get (S, j);
      if (Sj < S_max) {
	S_max = Sj;
	i_max = j;
      }
    }
    if (i_max != i) {
      /* swap eigenvalues */
      gsl_vector_swap_elements (S, i, i_max);
    }
  }

  return GSL_SUCCESS;
}
