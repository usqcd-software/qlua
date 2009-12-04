
const char *
gen_laplace_P(QDP_DiracPropagator *res, double a, double b, 
        QDP_ColorMatrix **u, QDP_DiracPropagator *x, int skip_axis)
{
    QLA_Real r1;
    if (skip_axis < 0 || QDP_ndim() <= skip_axis)
        r1 = a - 2 * b * QDP_ndim();
    else 
        r1 = a - 2 * b * (QDP_ndim() - 1);
    QDP_P_eq_r_times_P(res, &r1, x, *qCurrent);
    int i; 
    QDP_DiracPropagator *aux = QDP_create_P();
    QDP_P_eq_zero(aux, *qCurrent);
    for (i = 0; i < QDP_ndim(); i++) {
        if (i == skip_axis)
            continue;
        QDP_P_peq_M_times_sP(aux, u[i], x, QDP_neighbor[i], QDP_forward, 
                *qCurrent);
        QDP_P_peq_sMa_times_sP(aux, u[i], x, QDP_neighbor[i], QDP_backward,
                *qCurrent);
    }
    QLA_Real r2 = b;
    QDP_P_peq_r_times_P(res, &r2, aux, *qCurrent);
    QDP_destroy_P(aux);
    return NULL;
}
