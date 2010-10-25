#include "modules.h"                                                /* DEPS */
#include "qlua.h"                                                   /* DEPS */
#include "lattice.h"                                                /* DEPS */
#include "qmatrix.h"                                                /* DEPS */
#include "latdirprop.h"                                             /* DEPS */
#include "aff_io.h"                                                 /* DEPS */
#include "extras.h"                                                 /* DEPS */

struct b_arg {
	QLA_D3_DiracPropagator *Pd;
	QLA_D3_DiracPropagator *Pu11;
	QLA_D3_DiracPropagator *Pu12;
	QLA_D3_DiracPropagator *Pu21;
	QLA_D3_DiracPropagator *Pu22;
	QLA_D_Complex Sf[QDP_Ns][QDP_Ns];
	QLA_D_Complex Si_bar[QDP_Ns][QDP_Ns];
	QLA_D_Complex RTR[QDP_Ns][QDP_Ns];
};

static const struct { unsigned char bL, bR, ac[4][4]; } eps_eps[] = {
	{0, 0, {{1, 1, 2, 2}, {2, 2, 1, 1}, {1, 2, 2, 1}, {2, 1, 1, 2}}},
	{0, 1, {{1, 2, 2, 0}, {2, 0, 1, 2}, {1, 0, 2, 2}, {2, 2, 1, 0}}},
	{0, 2, {{1, 0, 2, 1}, {2, 1, 1, 0}, {1, 1, 2, 0}, {2, 0, 1, 1}}},
	{1, 0, {{0, 2, 2, 1}, {2, 1, 0, 2}, {0, 1, 2, 2}, {2, 2, 0, 1}}},
	{1, 1, {{0, 0, 2, 2}, {2, 2, 0, 0}, {0, 2, 2, 0}, {2, 0, 0, 2}}},
	{1, 2, {{0, 1, 2, 0}, {2, 0, 0, 1}, {0, 0, 2, 1}, {2, 1, 0, 0}}},
	{2, 0, {{0, 1, 1, 2}, {1, 2, 0, 1}, {0, 2, 1, 1}, {1, 1, 0, 2}}},
	{2, 1, {{0, 2, 1, 0}, {1, 0, 0, 2}, {0, 0, 1, 2}, {1, 2, 0, 0}}},
	{2, 2, {{0, 0, 1, 1}, {1, 1, 0, 0}, {0, 1, 1, 0}, {1, 0, 0, 1}}}};

static void
left_mul(QLA_D3_DiracPropagator *V,
		 QLA_D_Complex M[QDP_Ns][QDP_Ns],
		 QLA_D3_DiracPropagator *X)
{
	int mu, nu, eta, a, b;
	QLA_D3_P_eq_zero(V);
	for (mu = 0; mu < QDP_Ns; mu++) {
		for (nu = 0; nu < QDP_Ns; nu++) {
			QLA_D_Complex *z = &M[mu][nu];
			if ((QLA_real(*z) == 0.0) && (QLA_imag(*z) == 0.0))
				continue;
			for (eta = 0; eta < QDP_Ns; eta++) {
				for (a = 0; a < 3; a++) {
					for (b = 0; b < 3; b++) {
						QLA_D_Complex *v = &QLA_D3_elem_P(*V, a, mu, b, eta);
						QLA_D_C_peq_C_times_C(v, z, &QLA_D3_elem_P(*X, a, nu, b, eta));
					}
				}
			}
		}
	}
}

static void
right_mul(QLA_D3_DiracPropagator *V,
		  QLA_D3_DiracPropagator *X,
		  QLA_D_Complex M[QDP_Ns][QDP_Ns])
{
	int mu, nu, eta, a, b;
	QLA_D3_P_eq_zero(V);
	for (mu = 0; mu < QDP_Ns; mu++) {
		for (nu = 0; nu < QDP_Ns; nu++) {
			QLA_D_Complex *z = &M[mu][nu];
			if ((QLA_real(*z) == 0.0) && (QLA_imag(*z) == 0.0))
				continue;
			for (eta = 0; eta < QDP_Ns; eta++) {
				for (a = 0; a < 3; a++) {
					for (b = 0; b < 3; b++) {
						QLA_D_Complex *v = &QLA_D3_elem_P(*V, a, eta, b, nu);
						QLA_D_C_peq_C_times_C(v, &QLA_D3_elem_P(*X, a, eta, b, mu), z);
					}
				}
			}
		}
	}
}

static void
spin_trace(QLA_D3_ColorMatrix *V,
		   QLA_D_Complex M[QDP_Ns][QDP_Ns],
		   QLA_D3_DiracPropagator *X)
{
	int mu, nu, a, b;
	QLA_D3_M_eq_zero(V);
	for (mu = 0; mu < QDP_Ns; mu++) {
		for (nu = 0; nu < QDP_Ns; nu++) {
			QLA_D_Complex *z = &M[mu][nu];
			if ((QLA_real(*z) == 0.0) && (QLA_imag(*z) == 0.0))
				continue;
			for (a = 0; a < 3; a++) {
				for (b = 0; b < 3; b++) {
					QLA_D_Complex *v = &QLA_D3_elem_M(*V, a, b);
					QLA_D_C_peq_C_times_C(v, z, &QLA_D3_elem_P(*X, a, nu, b, mu));
				}
			}
		}
	}
}

static void
compute_uu(QLA_D_Complex *uu,
		   int alpha,
		   int alpha1,
		   const unsigned char ac[4],
		   QLA_D3_DiracPropagator *U11,
		   QLA_D3_ColorMatrix *U22,
		   QLA_D3_DiracPropagator *U12,
		   QLA_D3_DiracPropagator *U21)
{
	int gamma;

	/* U11(alpha,alpha')(aa') * U22(cc')  */
	QLA_C_eq_C_times_C(uu,
					   &QLA_D3_elem_P(*U11, ac[0], alpha, ac[1], alpha1),
					   &QLA_D3_elem_M(*U22, ac[2], ac[3]));
	for (gamma = 0; gamma < QDP_Ns; gamma++) {
		/* U12(alpha,gamma)(aa') * U21(gamma,alpha')(cc')*/
		QLA_C_peq_C_times_C(uu,
							&QLA_D3_elem_P(*U12, ac[0], alpha, ac[1], gamma),
							&QLA_D3_elem_P(*U21, ac[2], gamma, ac[3], alpha1));
	}
}

static void
duu(QLA_D_Complex *B, int idx, void *env)
{
	struct b_arg *arg = env;
	QLA_D3_DiracPropagator SDS;
	QLA_D3_DiracPropagator *U11 = &arg->Pu11[idx];
	QLA_D3_DiracPropagator *U12 = &arg->Pu12[idx];
	QLA_D3_DiracPropagator U21;
	QLA_D3_ColorMatrix U22;
	QLA_D_Complex uu[4];
	int k;

	left_mul(&U21, arg->Sf, &arg->Pd[idx]);
	right_mul(&SDS, &U21, arg->Si_bar);
	left_mul(&U21, arg->RTR, &arg->Pu21[idx]);
	spin_trace(&U22, arg->RTR, &arg->Pu22[idx]);
	QLA_D_C_eq_zero(B);
	for (k = 0; k < 9; k++) {
		int bL = eps_eps[k].bL;
		int bR = eps_eps[k].bR;
		int alpha, alpha1;
		for (alpha = 0; alpha < QDP_Ns; alpha++) {
			for (alpha1 = 0; alpha1 < QDP_Ns; alpha1++) {
				int j;
				for (j = 0; j < 4; j++)
					compute_uu(&uu[j], alpha, alpha1, eps_eps[k].ac[j], U11, &U22, U12, &U21);
				QLA_D_C_peq_C(&uu[0], &uu[1]);
				QLA_D_C_peq_C(&uu[2], &uu[3]);
				QLA_D_C_meq_C(&uu[0], &uu[2]);
				QLA_D_C_peq_C_times_C(B, &QLA_D3_elem_P(SDS, bL, alpha, bR, alpha1), &uu[0]);
			}
		}
	}
}

static void
c_get(int i, int j, QLA_D_Complex V[QDP_Ns][QDP_Ns], gsl_matrix_complex *m)
{
	gsl_complex z = gsl_matrix_complex_get(m, i, j);
	QLA_real(V[i][j]) = GSL_REAL(z);
	QLA_imag(V[i][j]) = GSL_IMAG(z);
}

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
		   gsl_matrix_complex *RTR)
{
	struct b_arg arg;
	int i, j;

	arg.Pd = QDP_D3_expose_P(Pd);
	arg.Pu11 = QDP_D3_expose_P(Pu11);
	arg.Pu12 = QDP_D3_expose_P(Pu12);
	arg.Pu21 = QDP_D3_expose_P(Pu21);
	arg.Pu22 = QDP_D3_expose_P(Pu22);

	for (i = 0; i < QDP_Ns; i++) {
		for (j = 0; j < QDP_Ns; j++) {
			c_get(i, j, arg.Sf, Sf);
			c_get(i, j, arg.Si_bar, Si_bar);
			c_get(i, j, arg.RTR, RTR);
		}
	}
	QDP_D_C_eq_funcia(B, duu, &arg, *S->qss);

	QDP_D3_reset_P(Pd);
	QDP_D3_reset_P(Pu11);
	QDP_D3_reset_P(Pu12);
	QDP_D3_reset_P(Pu21);
	QDP_D3_reset_P(Pu22);
}
