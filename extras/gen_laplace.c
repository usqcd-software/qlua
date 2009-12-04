#include <qlua.h>                                                   /* DEPS */
#include <lattice.h>                                                /* DEPS */

#include <aff_io.h>                                                 /* DEPS */
#include <extras.h>                                                 /* DEPS */

#define gen_laplace_ID              gen_laplace_P
#define QDP_Type                    QDP_DiracPropagator
#define QDP_create_ID               QDP_create_P
#define QDP_destroy_ID              QDP_destroy_P
#define QDP_ID_eq_zero              QDP_P_eq_zero
#define QDP_ID_eq_r_times_ID        QDP_P_eq_r_times_P
#define QDP_ID_peq_M_times_sID      QDP_P_peq_M_times_sP
#define QDP_ID_peq_sMa_times_sID    QDP_P_peq_sMa_times_sP
#define QDP_ID_peq_r_times_ID       QDP_P_peq_r_times_P
#include <gen_laplace_template.c>                       /* DEPS */

#define gen_laplace_ID              gen_laplace_V
#define QDP_Type                    QDP_ColorVector
#define QDP_create_ID               QDP_create_V
#define QDP_destroy_ID              QDP_destroy_V
#define QDP_ID_eq_zero              QDP_V_eq_zero
#define QDP_ID_eq_r_times_ID        QDP_V_eq_r_times_V
#define QDP_ID_peq_M_times_sID      QDP_V_peq_M_times_sV
#define QDP_ID_peq_sMa_times_sID    QDP_V_peq_sMa_times_sV
#define QDP_ID_peq_r_times_ID       QDP_V_peq_r_times_V
#include <gen_laplace_template.c>                       /* DEPS */

#define gen_laplace_ID              gen_laplace_M
#define QDP_Type                    QDP_ColorMatrix
#define QDP_create_ID               QDP_create_M
#define QDP_destroy_ID              QDP_destroy_M
#define QDP_ID_eq_zero              QDP_M_eq_zero
#define QDP_ID_eq_r_times_ID        QDP_M_eq_r_times_M
#define QDP_ID_peq_M_times_sID      QDP_M_peq_M_times_sM
#define QDP_ID_peq_sMa_times_sID    QDP_M_peq_sMa_times_sM
#define QDP_ID_peq_r_times_ID       QDP_M_peq_r_times_M
#include <gen_laplace_template.c>                       /* DEPS */

#define gen_laplace_ID              gen_laplace_D
#define QDP_Type                    QDP_DiracFermion
#define QDP_create_ID               QDP_create_D
#define QDP_destroy_ID              QDP_destroy_D
#define QDP_ID_eq_zero              QDP_D_eq_zero
#define QDP_ID_eq_r_times_ID        QDP_D_eq_r_times_D
#define QDP_ID_peq_M_times_sID      QDP_D_peq_M_times_sD
#define QDP_ID_peq_sMa_times_sID    QDP_D_peq_sMa_times_sD
#define QDP_ID_peq_r_times_ID       QDP_D_peq_r_times_D
#include <gen_laplace_template.c>                       /* DEPS */
