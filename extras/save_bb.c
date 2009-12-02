#include <qlua.h>                                                   /* DEPS */
#include <aff_io.h>                                                 /* DEPS */
#include <extras.h>                                                 /* DEPS */

const char *
save_bb(mAffWriter *aff_w,
        const char *aff_kpath,
        QDP_DiracPropagator *F,
        QDP_DiracPropagator *B,
        const int *csrc,             /* [qRank] */
        int tsnk,
        int n_qext,
        const int *qext,             /* [n_qext][qRank] */
        int time_rev,                /* 1 to reverse, 0 to not */
        int t_axis,                  /* 0-based */
        double bc_baryon_t)
{
    /* XXX WRITE ME */
    return NULL;
}
