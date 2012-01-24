#ifndef GAMMA_DEFS_H_jyn2akiqoKxfqtkLzwc1
#define GAMMA_DEFS_H_jyn2akiqoKxfqtkLzwc1

#define Gamma_dot_vec(i_G, i, c0,c1,c2,c3)    (Gamma ## i_G ##  _dot_vec_##i(c0,c1,c2,c3))
#define GammaT_dot_vec(i_G, i, c0,c1,c2,c3)   (Gamma ## i_G ## T_dot_vec_##i(c0,c1,c2,c3))
/* calculate Tr(Gamma[i_G] . C) */
#define TrGamma_dot(i_G, c00,c01,c02,c03, c10,c11,c12,c13, c20,c21,c22,c23, c30,c31,c32,c33) \
      (Gamma_dot_vec(i_G,0, c00,c10,c20,c30) \
     + Gamma_dot_vec(i_G,1, c01,c11,c21,c31) \
     + Gamma_dot_vec(i_G,2, c02,c12,c22,c32) \
     + Gamma_dot_vec(i_G,3, c03,c13,c23,c33))
/* calculate Tr(Gamma[i_G]^T . C) */
#define TrGammaT_dot(i_G, c00,c01,c02,c03, c10,c11,c12,c13, c20,c21,c22,c23, c30,c31,c32,c33) \
      (GammaT_dot_vec(i_G,0, c00,c10,c20,c30) \
     + GammaT_dot_vec(i_G,1, c01,c11,c21,c31) \
     + GammaT_dot_vec(i_G,2, c02,c12,c22,c32) \
     + GammaT_dot_vec(i_G,3, c03,c13,c23,c33))

#define complex2gsl(c) (gsl_complex_rect(creal(c), cimag(c)))
#define gsl2complex(c) (GSL_REAL(c) + GSL_IMAG(c)*I)

#define qla2complex(c) (QLA_real(c) + QLA_imag(c)*I)
#define qla2gsl(c)      (gsl_complex_rect(QLA_real(c), QLA_imag(c)))

/* tests:
    Gamma_dot_vec(7, 2, c0,c1,c2,c3)
    GammaT_dot_vec(7, 2, c0,c1,c2,c3)
    TrGamma_dot(7, c00,c01,c02,c03, c10,c11,c12,c13, c20,c21,c22,c23, c30,c31,c32,c33)
    TrGammaT_dot(7, c00,c01,c02,c03, c10,c11,c12,c13, c20,c21,c22,c23, c30,c31,c32,c33)
  run as
  $ gcc -E gamma_dgr.h
 */

#endif/*GAMMA_DEFS_H_jyn2akiqoKxfqtkLzwc1*/
