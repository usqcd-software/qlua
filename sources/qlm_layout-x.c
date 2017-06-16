/* Let the ugliness begin! */
#define cplx_qla2int(int_x,qla_x) int_x = QLA_real(qla_x) + I*QLA_imag(qla_x)
#define cplx_int2qla(qla_x,int_x) QLA_c_eq_r_plus_ir((qla_x), creal(int_x), cimag(int_x))
#define pcint pint complex
static void
Qpp(fc_qla2int_R)(void *restrict dst_, int dst_idx, const void *restrict src_, int src_idx, int iarr, void *arg_)
{
    qlmData *d = (qlmData *)arg_;
    const QT(_Real) *qx = ((QT(_Real) **)src_)[iarr] + src_idx;
    pint *ix = ((pint *)dst_) + dst_idx * d->site_num_len + iarr;
    *ix = *qx;
}
static void
Qpp(fc_int2qla_R)(void *restrict dst_, int dst_idx, const void *restrict src_, int src_idx, int iarr, void *arg_)
{
    qlmData *d = (qlmData *)arg_;
    pint *ix = ((pint *)src_) + src_idx * d->site_num_len + iarr;
    QT(_Real) *qx = ((QT(_Real) **)dst_)[iarr] + dst_idx;
    *qx = *ix;
}

static void
Qpp(fc_qla2int_C)(void *restrict dst_, int dst_idx, const void *restrict src_, int src_idx, int iarr, void *arg_)
{
    qlmData *d = (qlmData *)arg_;
    QT(_Complex) *qx = ((QT(_Complex) **)src_)[iarr] + src_idx;
    pcint *ix = (pcint *)dst_ + dst_idx * d->site_num_len + iarr;
    cplx_qla2int(*ix, *qx);
}
static void
Qpp(fc_int2qla_C)(void *restrict dst_, int dst_idx, const void *restrict src_, int src_idx, int iarr, void *arg_)
{
    qlmData *d = (qlmData *)arg_;
    pcint *ix = (pcint *)src_ + src_idx * d->site_num_len;
    QT(_Complex) *qx = ((QT(_Complex) **)dst_)[iarr] + dst_idx;
    cplx_int2qla(*qx, *ix);
}
#undef QLM_INDEX_C

/* TODO C */

#define Qnc 2
#define Qns 4
#define Qppc(x) Qpp(x##2)
#define QTc(x) QT(2##x)
#include "qlm_layout-y.c"

#define Qnc 3
#define Qns 4
#define Qppc(x) Qpp(x##3)
#define QTc(x) QT(3##x)
#include "qlm_layout-y.c"

#undef Qpp
#undef pint
#undef QT

