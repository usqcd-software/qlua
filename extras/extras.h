#ifndef MARK_E390D5B2_478D_48BF_B753_64072003402B
#define MARK_E390D5B2_478D_48BF_B753_64072003402B

#ifdef HAS_GSL
#define GSL_RANGE_CHECK_OFF
#define HAVE_INLINE
#include <gsl/gsl_blas.h>
#include <gsl/gsl_matrix_double.h>

#include "aff_io.h"                                                 /* DEPS */

void
baryon_duu(mLattice *S,
		   QDP_D_Complex *B,
		   QDP_D3_DiracPropagator *Pd,
		   QDP_D3_DiracPropagator *Pu11,
		   QDP_D3_DiracPropagator *Pu12,
		   QDP_D3_DiracPropagator *Pu21,
		   QDP_D3_DiracPropagator *Pu22,
		   gsl_matrix_complex *Sf,
		   gsl_matrix_complex *Si_bar,
		   gsl_matrix_complex *RTR);
#endif /* defined HAS_GSL */

const char *
gen_laplace_V(lua_State *L,
              mLattice *S,
              QDP_D3_ColorVector *result,
              double a,
              double b,
              QDP_D3_ColorMatrix **gauge,
              QDP_D3_ColorVector *source,
              int skip_axis);


const char *
gen_laplace_M(lua_State *L,
              mLattice *S,
              QDP_D3_ColorMatrix *result,
              double a,
              double b,
              QDP_D3_ColorMatrix **gauge,
              QDP_D3_ColorMatrix *source,
              int skip_axis);

const char *
gen_laplace_D(lua_State *L,
              mLattice *S,
              QDP_D3_DiracFermion *result,
              double a,
              double b,
              QDP_D3_ColorMatrix **gauge,
              QDP_D3_DiracFermion *source,
              int skip_axis);

const char *
gen_laplace_P(lua_State *L,
              mLattice *S,
              QDP_D3_DiracPropagator *result,
              double a,
              double b,
              QDP_D3_ColorMatrix **gauge,
              QDP_D3_DiracPropagator *source,
              int skip_axis);

int q_save_bb(lua_State *L);


int q_save_laph_wf_baryon_pwave(lua_State *L);
int q_save_q2pt_list(lua_State *L);
int q_save_q2pt(lua_State *L);
int q_save_q3pt_0deriv_selectspin(lua_State *L);
int q_save_npr_prop(lua_State *L);
int q_save_npr_2qvertex(lua_State *L);

const char *
extra_lat_fourier_qla(lua_State *L,
                      mLattice *S, 
                      QLA_D_Complex *qla_y, 
                      QLA_D_Complex *qla_x,
                      int rec_len, 
                      int ft_sign, 
                      int ft_dir);

const char *
extra_lat_fourier_qla_onedim(lua_State *L,
                             mLattice *S, 
                             QLA_D_Complex *qla_y, 
                             QLA_D_Complex *qla_x,
                             int rec_len, 
                             int ft_sign, 
                             int ft_dir);

const char *
extra_lat_fourier_qla_full(lua_State *L,
                           mLattice *S, 
                           QLA_D_Complex *qla_y, 
                           QLA_D_Complex *qla_x,
                           int rec_len, 
                           int ft_sign);

int init_extras(lua_State *L);
int q_test_qopt(lua_State *L);
void fini_extras(void);

#endif /* !defined(MARK_E390D5B2_478D_48BF_B753_64072003402B) */
