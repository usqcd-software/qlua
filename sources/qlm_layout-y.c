#define QLM_INDEX_V(a,c)  ((c) + Qnc*(a))
static void
Qppc(fc_qla2int_V)(void *restrict dst_, int dst_idx, const void *restrict src_, int src_idx, int iarr, void *arg_)
{
    qlmData *d = (qlmData *)arg_;
    QT(_ColorVector) *qx0 = ((QT(_ColorVector) **)src_)[iarr] + src_idx;
    pcint *ix0 = (pcint *)dst_ + dst_idx * d->site_num_len;
    for (int c = 0 ; c < Qnc ; c++) {
        QT(_Complex) *qx = &QTc(_elem_V)(*qx0, c);
        pcint *ix = ix0 + QLM_INDEX_V(iarr, c);
        cplx_qla2int(*ix, *qx);
    }
}
static void
Qppc(fc_int2qla_V)(void *restrict dst_, int dst_idx, const void *restrict src_, int src_idx, int iarr, void *arg_)
{
    qlmData *d = (qlmData *)arg_;
    pcint *ix0 = (pcint *)src_ + src_idx * d->site_num_len;
    QT(_ColorVector) *qx0 = ((QT(_ColorVector) **)dst_)[iarr] + dst_idx;
    for (int c = 0 ; c < Qnc ; c++) {
        pcint *ix = ix0 + QLM_INDEX_V(iarr, c);
        QT(_Complex) *qx = &QTc(_elem_V)(*qx0, c);
        cplx_int2qla(*qx, *ix);
    }
}
#undef QLM_INDEX_V

#define QLM_INDEX_M(a,c,c2)  ((c) + Qnc*((c2) + Qnc*(a)))
static void
Qppc(fc_qla2int_M)(void *restrict dst_, int dst_idx, const void *restrict src_, int src_idx, int iarr, void *arg_)
{
    qlmData *d = (qlmData *)arg_;
    QT(_ColorMatrix) *qx0 = ((QT(_ColorMatrix) **)src_)[iarr] + src_idx;
    pcint *ix0 = (pcint *)dst_ + dst_idx * d->site_num_len;
    for (int c = 0 ; c < Qnc ; c++) for (int c2 = 0 ; c2 < Qnc ; c2++) {
        QT(_Complex) *qx = &QTc(_elem_M)(*qx0, c,c2);
        pcint *ix = ix0 + QLM_INDEX_M(iarr, c,c2);
        cplx_qla2int(*ix, *qx);
    }
}
static void
Qppc(fc_int2qla_M)(void *restrict dst_, int dst_idx, const void *restrict src_, int src_idx, int iarr, void *arg_)
{
    qlmData *d = (qlmData *)arg_;
    pcint *ix0 = (pcint *)src_ + src_idx * d->site_num_len;
    QT(_ColorMatrix) *qx0 = ((QT(_ColorMatrix) **)dst_)[iarr] + dst_idx;
    for (int c = 0 ; c < Qnc ; c++) for (int c2 = 0 ; c2 < Qnc ; c2++) {
        pcint *ix = ix0 + QLM_INDEX_M(iarr, c,c2);
        QT(_Complex) *qx = &QTc(_elem_M)(*qx0, c,c2);
        cplx_int2qla(*qx, *ix);
    }
}
#undef QLM_INDEX_M

#define QLM_INDEX_D(a,c,s)  ((s) + Qns*((c) + Qnc*(a)))
static void
Qppc(fc_qla2int_D)(void *restrict dst_, int dst_idx, const void *restrict src_, int src_idx, int iarr, void *arg_)
{
    qlmData *d = (qlmData *)arg_;
    QT(_DiracFermion) *qx0 = ((QT(_DiracFermion) **)src_)[iarr] + src_idx;
    pcint *ix0 = (pcint *)dst_ + dst_idx * d->site_num_len;
    for (int c = 0 ; c < Qnc ; c++) for (int s = 0 ; s < Qns ; s++) {
        QT(_Complex) *qx = &QTc(_elem_D)(*qx0, c, s);
        pcint *ix = ix0 + QLM_INDEX_D(iarr, c, s);
        cplx_qla2int(*ix, *qx);
    }
}
static void
Qppc(fc_int2qla_D)(void *restrict dst_, int dst_idx, const void *restrict src_, int src_idx, int iarr, void *arg_)
{
    qlmData *d = (qlmData *)arg_;
    pcint *ix0 = (pcint *)src_ + src_idx * d->site_num_len;
    QT(_DiracFermion) *qx0 = ((QT(_DiracFermion) **)dst_)[iarr] + dst_idx;
    for (int c = 0 ; c < Qnc ; c++) for (int s = 0 ; s < Qns ; s++) {
        pcint *ix = ix0 + QLM_INDEX_D(iarr, c, s);
        QT(_Complex) *qx = &QTc(_elem_D)(*qx0, c, s);
        cplx_int2qla(*qx, *ix);
    }
}
#undef QLM_INDEX_D

#define QLM_INDEX_P(a,c,s,c2,s2)  ((s) + Qns*((c) + Qnc*((s2) + Qns*((c2) + Qnc*(a)))))
static void
Qppc(fc_qla2int_P)(void *restrict dst_, int dst_idx, const void *restrict src_, int src_idx, int iarr, void *arg_)
{
    qlmData *d = (qlmData *)arg_;
    QT(_DiracPropagator) *qx0 = ((QT(_DiracPropagator) **)src_)[iarr] + src_idx;
    pcint *ix0 = (pcint *)dst_ + dst_idx * d->site_num_len;
    for (int c = 0 ; c < Qnc ; c++ ) for (int s = 0 ; s < Qns ; s++ ) 
    for (int c2= 0 ; c2< Qnc ; c2++) for (int s2= 0 ; s2< Qns ; s2++) {
        QT(_Complex) *qx = &QTc(_elem_P)(*qx0, c,s,c2,s2);
        pcint *ix = ix0 + QLM_INDEX_P(iarr, c,s,c2,s2);
        cplx_qla2int(*ix, *qx);
    }
}
static void
Qppc(fc_int2qla_P)(void *restrict dst_, int dst_idx, const void *restrict src_, int src_idx, int iarr, void *arg_)
{
    qlmData *d = (qlmData *)arg_;
    pcint *ix0 = (pcint *)src_ + src_idx * d->site_num_len;
    QT(_DiracPropagator) *qx0 = ((QT(_DiracPropagator) **)dst_)[iarr] + dst_idx;
    for (int c = 0 ; c < Qnc ; c++) for (int s = 0 ; s < Qns ; s++)
    for (int c2= 0 ; c2< Qnc ; c2++) for (int s2= 0 ; s2< Qns ; s2++) {
        pcint *ix = ix0 + QLM_INDEX_P(iarr, c,s,c2,s2);
        QT(_Complex) *qx = &QTc(_elem_P)(*qx0, c,s,c2,s2);
        cplx_int2qla(*qx, *ix);
    }
}
#undef QLM_INDEX_P

/*TODO V M P */

#undef Qnc
#undef Qns
#undef Qppc
#undef QTc
