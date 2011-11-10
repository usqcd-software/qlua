#ifndef GAMMA_DGR_GEN_H_zx7Jzuibswaxfwcmtg7O
#define GAMMA_DGR_GEN_H_zx7Jzuibswaxfwcmtg7O
/* DeGrand-Rossi representation of gamma matrices */

/* G = [1           g1          g2          g1 g2           
        g3          g1 g3       g2 g3       g1 g2 g3 
        g4          g1 g4       g2 g4       g1 g2 g4 
        g3 g4       g1 g3 g4    g2 g3 g4    g1 g2 g3 g4 ]
 */


/* representation of (Gamma.X) as permutation:
   (Gamma.X)[i] = Gamma_dot_mul[i] * X[Gamma_dot_ind[i]]
   or, equivalently,
   Gamma[i,j] = Gamma_dot_mul[i] * \delta(Gamma_dot_ind[i]==j)
 */
#define Gamma_dot_ind   { \
        { 0, 1, 2, 3 }, { 3, 2, 1, 0 }, { 3, 2, 1, 0 }, { 0, 1, 2, 3 }, \
        { 2, 3, 0, 1 }, { 1, 0, 3, 2 }, { 1, 0, 3, 2 }, { 2, 3, 0, 1 }, \
        { 2, 3, 0, 1 }, { 1, 0, 3, 2 }, { 1, 0, 3, 2 }, { 2, 3, 0, 1 }, \
        { 0, 1, 2, 3 }, { 3, 2, 1, 0 }, { 3, 2, 1, 0 }, { 0, 1, 2, 3 }  }
#define Gamma_dot_mul   { \
        { 1, 1, 1, 1 }, { I, I,-I,-I }, {-1, 1, 1,-1 }, {-I, I,-I, I }, \
        { I,-I,-I, I }, {-1, 1,-1, 1 }, {-I,-I,-I,-I }, { 1, 1,-1,-1 }, \
        { 1, 1, 1, 1 }, { I, I,-I,-I }, {-1, 1, 1,-1 }, {-I, I,-I, I }, \
        { I,-I,-I, I }, {-1, 1,-1, 1 }, {-I,-I,-I,-I }, { 1, 1,-1,-1 }  }

/* representation of (Gamma^T.X) as permutation:
   (Gamma^T.X)[i] = GammaT_dot_mul[i] * X[GammaT_dot_ind[i]]
   or, equivalently,
   Gamma^T[i,j] = GammaT_dot_mul[i] * \delta(GammaT_dot_ind[i]==j)
 */
#define GammaT_dot_ind  { \
        { 0, 1, 2, 3 }, { 3, 2, 1, 0 }, { 3, 2, 1, 0 }, { 0, 1, 2, 3 }, \
        { 2, 3, 0, 1 }, { 1, 0, 3, 2 }, { 1, 0, 3, 2 }, { 2, 3, 0, 1 }, \
        { 2, 3, 0, 1 }, { 1, 0, 3, 2 }, { 1, 0, 3, 2 }, { 2, 3, 0, 1 }, \
        { 0, 1, 2, 3 }, { 3, 2, 1, 0 }, { 3, 2, 1, 0 }, { 0, 1, 2, 3 }  }
#define GammaT_dot_mul  { \
        { 1, 1, 1, 1 }, {-I,-I, I, I }, {-1, 1, 1,-1 }, {-I, I,-I, I }, \
        {-I, I, I,-I }, { 1,-1, 1,-1 }, {-I,-I,-I,-I }, {-1,-1, 1, 1 }, \
        { 1, 1, 1, 1 }, { I, I,-I,-I }, { 1,-1,-1, 1 }, {-I, I,-I, I }, \
        { I,-I,-I, I }, { 1,-1, 1,-1 }, {-I,-I,-I,-I }, { 1, 1,-1,-1 }  }


/* #define representation (GammaN.C)[i] and (GammaN^T.C), C=(c0,c1,c2,c3) */
#define Gamma0_dot_vec_0(c0,c1,c2,c3) ( 1*c0)
#define Gamma0T_dot_vec_0(c0,c1,c2,c3) ( 1*c0)
#define Gamma0_dot_vec_1(c0,c1,c2,c3) ( 1*c1)
#define Gamma0T_dot_vec_1(c0,c1,c2,c3) ( 1*c1)
#define Gamma0_dot_vec_2(c0,c1,c2,c3) ( 1*c2)
#define Gamma0T_dot_vec_2(c0,c1,c2,c3) ( 1*c2)
#define Gamma0_dot_vec_3(c0,c1,c2,c3) ( 1*c3)
#define Gamma0T_dot_vec_3(c0,c1,c2,c3) ( 1*c3)
#define Gamma1_dot_vec_0(c0,c1,c2,c3) ( I*c3)
#define Gamma1T_dot_vec_0(c0,c1,c2,c3) (-I*c3)
#define Gamma1_dot_vec_1(c0,c1,c2,c3) ( I*c2)
#define Gamma1T_dot_vec_1(c0,c1,c2,c3) (-I*c2)
#define Gamma1_dot_vec_2(c0,c1,c2,c3) (-I*c1)
#define Gamma1T_dot_vec_2(c0,c1,c2,c3) ( I*c1)
#define Gamma1_dot_vec_3(c0,c1,c2,c3) (-I*c0)
#define Gamma1T_dot_vec_3(c0,c1,c2,c3) ( I*c0)
#define Gamma2_dot_vec_0(c0,c1,c2,c3) (-1*c3)
#define Gamma2T_dot_vec_0(c0,c1,c2,c3) (-1*c3)
#define Gamma2_dot_vec_1(c0,c1,c2,c3) ( 1*c2)
#define Gamma2T_dot_vec_1(c0,c1,c2,c3) ( 1*c2)
#define Gamma2_dot_vec_2(c0,c1,c2,c3) ( 1*c1)
#define Gamma2T_dot_vec_2(c0,c1,c2,c3) ( 1*c1)
#define Gamma2_dot_vec_3(c0,c1,c2,c3) (-1*c0)
#define Gamma2T_dot_vec_3(c0,c1,c2,c3) (-1*c0)
#define Gamma3_dot_vec_0(c0,c1,c2,c3) (-I*c0)
#define Gamma3T_dot_vec_0(c0,c1,c2,c3) (-I*c0)
#define Gamma3_dot_vec_1(c0,c1,c2,c3) ( I*c1)
#define Gamma3T_dot_vec_1(c0,c1,c2,c3) ( I*c1)
#define Gamma3_dot_vec_2(c0,c1,c2,c3) (-I*c2)
#define Gamma3T_dot_vec_2(c0,c1,c2,c3) (-I*c2)
#define Gamma3_dot_vec_3(c0,c1,c2,c3) ( I*c3)
#define Gamma3T_dot_vec_3(c0,c1,c2,c3) ( I*c3)
#define Gamma4_dot_vec_0(c0,c1,c2,c3) ( I*c2)
#define Gamma4T_dot_vec_0(c0,c1,c2,c3) (-I*c2)
#define Gamma4_dot_vec_1(c0,c1,c2,c3) (-I*c3)
#define Gamma4T_dot_vec_1(c0,c1,c2,c3) ( I*c3)
#define Gamma4_dot_vec_2(c0,c1,c2,c3) (-I*c0)
#define Gamma4T_dot_vec_2(c0,c1,c2,c3) ( I*c0)
#define Gamma4_dot_vec_3(c0,c1,c2,c3) ( I*c1)
#define Gamma4T_dot_vec_3(c0,c1,c2,c3) (-I*c1)
#define Gamma5_dot_vec_0(c0,c1,c2,c3) (-1*c1)
#define Gamma5T_dot_vec_0(c0,c1,c2,c3) ( 1*c1)
#define Gamma5_dot_vec_1(c0,c1,c2,c3) ( 1*c0)
#define Gamma5T_dot_vec_1(c0,c1,c2,c3) (-1*c0)
#define Gamma5_dot_vec_2(c0,c1,c2,c3) (-1*c3)
#define Gamma5T_dot_vec_2(c0,c1,c2,c3) ( 1*c3)
#define Gamma5_dot_vec_3(c0,c1,c2,c3) ( 1*c2)
#define Gamma5T_dot_vec_3(c0,c1,c2,c3) (-1*c2)
#define Gamma6_dot_vec_0(c0,c1,c2,c3) (-I*c1)
#define Gamma6T_dot_vec_0(c0,c1,c2,c3) (-I*c1)
#define Gamma6_dot_vec_1(c0,c1,c2,c3) (-I*c0)
#define Gamma6T_dot_vec_1(c0,c1,c2,c3) (-I*c0)
#define Gamma6_dot_vec_2(c0,c1,c2,c3) (-I*c3)
#define Gamma6T_dot_vec_2(c0,c1,c2,c3) (-I*c3)
#define Gamma6_dot_vec_3(c0,c1,c2,c3) (-I*c2)
#define Gamma6T_dot_vec_3(c0,c1,c2,c3) (-I*c2)
#define Gamma7_dot_vec_0(c0,c1,c2,c3) ( 1*c2)
#define Gamma7T_dot_vec_0(c0,c1,c2,c3) (-1*c2)
#define Gamma7_dot_vec_1(c0,c1,c2,c3) ( 1*c3)
#define Gamma7T_dot_vec_1(c0,c1,c2,c3) (-1*c3)
#define Gamma7_dot_vec_2(c0,c1,c2,c3) (-1*c0)
#define Gamma7T_dot_vec_2(c0,c1,c2,c3) ( 1*c0)
#define Gamma7_dot_vec_3(c0,c1,c2,c3) (-1*c1)
#define Gamma7T_dot_vec_3(c0,c1,c2,c3) ( 1*c1)
#define Gamma8_dot_vec_0(c0,c1,c2,c3) ( 1*c2)
#define Gamma8T_dot_vec_0(c0,c1,c2,c3) ( 1*c2)
#define Gamma8_dot_vec_1(c0,c1,c2,c3) ( 1*c3)
#define Gamma8T_dot_vec_1(c0,c1,c2,c3) ( 1*c3)
#define Gamma8_dot_vec_2(c0,c1,c2,c3) ( 1*c0)
#define Gamma8T_dot_vec_2(c0,c1,c2,c3) ( 1*c0)
#define Gamma8_dot_vec_3(c0,c1,c2,c3) ( 1*c1)
#define Gamma8T_dot_vec_3(c0,c1,c2,c3) ( 1*c1)
#define Gamma9_dot_vec_0(c0,c1,c2,c3) ( I*c1)
#define Gamma9T_dot_vec_0(c0,c1,c2,c3) ( I*c1)
#define Gamma9_dot_vec_1(c0,c1,c2,c3) ( I*c0)
#define Gamma9T_dot_vec_1(c0,c1,c2,c3) ( I*c0)
#define Gamma9_dot_vec_2(c0,c1,c2,c3) (-I*c3)
#define Gamma9T_dot_vec_2(c0,c1,c2,c3) (-I*c3)
#define Gamma9_dot_vec_3(c0,c1,c2,c3) (-I*c2)
#define Gamma9T_dot_vec_3(c0,c1,c2,c3) (-I*c2)
#define Gamma10_dot_vec_0(c0,c1,c2,c3) (-1*c1)
#define Gamma10T_dot_vec_0(c0,c1,c2,c3) ( 1*c1)
#define Gamma10_dot_vec_1(c0,c1,c2,c3) ( 1*c0)
#define Gamma10T_dot_vec_1(c0,c1,c2,c3) (-1*c0)
#define Gamma10_dot_vec_2(c0,c1,c2,c3) ( 1*c3)
#define Gamma10T_dot_vec_2(c0,c1,c2,c3) (-1*c3)
#define Gamma10_dot_vec_3(c0,c1,c2,c3) (-1*c2)
#define Gamma10T_dot_vec_3(c0,c1,c2,c3) ( 1*c2)
#define Gamma11_dot_vec_0(c0,c1,c2,c3) (-I*c2)
#define Gamma11T_dot_vec_0(c0,c1,c2,c3) (-I*c2)
#define Gamma11_dot_vec_1(c0,c1,c2,c3) ( I*c3)
#define Gamma11T_dot_vec_1(c0,c1,c2,c3) ( I*c3)
#define Gamma11_dot_vec_2(c0,c1,c2,c3) (-I*c0)
#define Gamma11T_dot_vec_2(c0,c1,c2,c3) (-I*c0)
#define Gamma11_dot_vec_3(c0,c1,c2,c3) ( I*c1)
#define Gamma11T_dot_vec_3(c0,c1,c2,c3) ( I*c1)
#define Gamma12_dot_vec_0(c0,c1,c2,c3) ( I*c0)
#define Gamma12T_dot_vec_0(c0,c1,c2,c3) ( I*c0)
#define Gamma12_dot_vec_1(c0,c1,c2,c3) (-I*c1)
#define Gamma12T_dot_vec_1(c0,c1,c2,c3) (-I*c1)
#define Gamma12_dot_vec_2(c0,c1,c2,c3) (-I*c2)
#define Gamma12T_dot_vec_2(c0,c1,c2,c3) (-I*c2)
#define Gamma12_dot_vec_3(c0,c1,c2,c3) ( I*c3)
#define Gamma12T_dot_vec_3(c0,c1,c2,c3) ( I*c3)
#define Gamma13_dot_vec_0(c0,c1,c2,c3) (-1*c3)
#define Gamma13T_dot_vec_0(c0,c1,c2,c3) ( 1*c3)
#define Gamma13_dot_vec_1(c0,c1,c2,c3) ( 1*c2)
#define Gamma13T_dot_vec_1(c0,c1,c2,c3) (-1*c2)
#define Gamma13_dot_vec_2(c0,c1,c2,c3) (-1*c1)
#define Gamma13T_dot_vec_2(c0,c1,c2,c3) ( 1*c1)
#define Gamma13_dot_vec_3(c0,c1,c2,c3) ( 1*c0)
#define Gamma13T_dot_vec_3(c0,c1,c2,c3) (-1*c0)
#define Gamma14_dot_vec_0(c0,c1,c2,c3) (-I*c3)
#define Gamma14T_dot_vec_0(c0,c1,c2,c3) (-I*c3)
#define Gamma14_dot_vec_1(c0,c1,c2,c3) (-I*c2)
#define Gamma14T_dot_vec_1(c0,c1,c2,c3) (-I*c2)
#define Gamma14_dot_vec_2(c0,c1,c2,c3) (-I*c1)
#define Gamma14T_dot_vec_2(c0,c1,c2,c3) (-I*c1)
#define Gamma14_dot_vec_3(c0,c1,c2,c3) (-I*c0)
#define Gamma14T_dot_vec_3(c0,c1,c2,c3) (-I*c0)
#define Gamma15_dot_vec_0(c0,c1,c2,c3) ( 1*c0)
#define Gamma15T_dot_vec_0(c0,c1,c2,c3) ( 1*c0)
#define Gamma15_dot_vec_1(c0,c1,c2,c3) ( 1*c1)
#define Gamma15T_dot_vec_1(c0,c1,c2,c3) ( 1*c1)
#define Gamma15_dot_vec_2(c0,c1,c2,c3) (-1*c2)
#define Gamma15T_dot_vec_2(c0,c1,c2,c3) (-1*c2)
#define Gamma15_dot_vec_3(c0,c1,c2,c3) (-1*c3)
#define Gamma15T_dot_vec_3(c0,c1,c2,c3) (-1*c3)

#endif/*GAMMA_DGR_GEN_H_zx7Jzuibswaxfwcmtg7O*/
