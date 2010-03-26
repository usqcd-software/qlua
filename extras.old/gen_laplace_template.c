
const char *
gen_laplace_ID(lua_State *L,
        QDP_Type *res, double a, double b, 
        QDP_ColorMatrix **u, QDP_Type *x, int skip_axis)
{
    QLA_Real r1;
    if (skip_axis < 0 || qRank <= skip_axis)
        r1 = a - 2 * b * qRank;
    else 
        r1 = a - 2 * b * (qRank - 1);
    QDP_ID_eq_r_times_ID(res, &r1, x, *qCurrent);
    int i; 
    QDP_Type *aux = QDP_create_ID();
    QDP_ID_eq_zero(aux, *qCurrent);
    for (i = 0; i < qRank; i++) {
        if (i == skip_axis)
            continue;
        QDP_ID_peq_M_times_sID(aux, u[i], x, QDP_neighbor[i], QDP_forward, 
                *qCurrent);
        QDP_ID_peq_sMa_times_sID(aux, u[i], x, QDP_neighbor[i], QDP_backward,
                *qCurrent);
    }
    QLA_Real r2 = b;
    QDP_ID_peq_r_times_ID(res, &r2, aux, *qCurrent);
    QDP_destroy_ID(aux);
    return NULL;
}
#undef gen_laplace_ID
#undef QDP_Type
#undef QDP_create_ID
#undef QDP_destroy_ID
#undef QDP_ID_eq_zero
#undef QDP_ID_eq_r_times_ID
#undef QDP_ID_peq_M_times_sID
#undef QDP_ID_peq_sMa_times_sID
#undef QDP_ID_peq_r_times_ID
