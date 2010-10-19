
const char *
gen_laplace_ID(lua_State *L, mLattice *S,
               QDP_Type *res, double a, double b, 
               QDP_D3_ColorMatrix **u, QDP_Type *x, int skip_axis)
{
    QLA_Real r1;
    QDP_Shift *neighbor = QDP_neighbor_L(S->lat);

    if (skip_axis < 0 || S->rank <= skip_axis)
        r1 = a - 2 * b * S->rank;
    else 
        r1 = a - 2 * b * (S->rank - 1);
    QDP_ID_eq_r_times_ID(res, &r1, x, *S->qss);
    int i; 
    QDP_Type *aux = QDP_create_ID(S->lat);
    QDP_ID_eq_zero(aux, *S->qss);
    for (i = 0; i < S->rank; i++) {
        if (i == skip_axis)
            continue;
        QDP_ID_peq_M_times_sID(aux, u[i], x, neighbor[i], QDP_forward, 
                *S->qss);
        QDP_ID_peq_sMa_times_sID(aux, u[i], x, neighbor[i], QDP_backward,
                *S->qss);
    }
    QLA_Real r2 = b;
    QDP_ID_peq_r_times_ID(res, &r2, aux, *S->qss);
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
