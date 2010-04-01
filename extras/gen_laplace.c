#include "qlua.h"                                                   /* DEPS */
#include "lattice.h"                                                /* DEPS */

#include "aff_io.h"                                                 /* DEPS */
#include "extras.h"                                                 /* DEPS */

#define gen_laplace_ID              gen_laplace_P
#define QDP_Type                    QDP_D3_DiracPropagator
#define QDP_create_ID               QDP_D3_create_P_L
#define QDP_destroy_ID              QDP_D3_destroy_P
#define QDP_ID_eq_zero              QDP_D3_P_eq_zero
#define QDP_ID_eq_r_times_ID        QDP_D3_P_eq_r_times_P
#define QDP_ID_peq_M_times_sID      QDP_D3_P_peq_M_times_sP
#define QDP_ID_peq_sMa_times_sID    QDP_D3_P_peq_sMa_times_sP
#define QDP_ID_peq_r_times_ID       QDP_D3_P_peq_r_times_P
#include "gen_laplace_template.c"                       /* DEPS */

#define gen_laplace_ID              gen_laplace_V
#define QDP_Type                    QDP_D3_ColorVector
#define QDP_create_ID               QDP_D3_create_V_L
#define QDP_destroy_ID              QDP_D3_destroy_V
#define QDP_ID_eq_zero              QDP_D3_V_eq_zero
#define QDP_ID_eq_r_times_ID        QDP_D3_V_eq_r_times_V
#define QDP_ID_peq_M_times_sID      QDP_D3_V_peq_M_times_sV
#define QDP_ID_peq_sMa_times_sID    QDP_D3_V_peq_sMa_times_sV
#define QDP_ID_peq_r_times_ID       QDP_D3_V_peq_r_times_V
#include "gen_laplace_template.c"                       /* DEPS */

#define gen_laplace_ID              gen_laplace_M
#define QDP_Type                    QDP_D3_ColorMatrix
#define QDP_create_ID               QDP_D3_create_M_L
#define QDP_destroy_ID              QDP_D3_destroy_M
#define QDP_ID_eq_zero              QDP_D3_M_eq_zero
#define QDP_ID_eq_r_times_ID        QDP_D3_M_eq_r_times_M
#define QDP_ID_peq_M_times_sID      QDP_D3_M_peq_M_times_sM
#define QDP_ID_peq_sMa_times_sID    QDP_D3_M_peq_sMa_times_sM
#define QDP_ID_peq_r_times_ID       QDP_D3_M_peq_r_times_M
#include "gen_laplace_template.c"                       /* DEPS */

#define gen_laplace_ID              gen_laplace_D
#define QDP_Type                    QDP_D3_DiracFermion
#define QDP_create_ID               QDP_D3_create_D_L
#define QDP_destroy_ID              QDP_D3_destroy_D
#define QDP_ID_eq_zero              QDP_D3_D_eq_zero
#define QDP_ID_eq_r_times_ID        QDP_D3_D_eq_r_times_D
#define QDP_ID_peq_M_times_sID      QDP_D3_D_peq_M_times_sD
#define QDP_ID_peq_sMa_times_sID    QDP_D3_D_peq_sMa_times_sD
#define QDP_ID_peq_r_times_ID       QDP_D3_D_peq_r_times_D
#include "gen_laplace_template.c"                       /* DEPS */
