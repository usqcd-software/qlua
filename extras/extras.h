#ifndef MARK_E390D5B2_478D_48BF_B753_64072003402B
#define MARK_E390D5B2_478D_48BF_B753_64072003402B

const char *
save_bb(lua_State *L,
        mAffWriter *aff_w,
        const char *aff_kpath,
        QDP_DiracPropagator *F,
        QDP_DiracPropagator *B,
        const int *csrc,             /* [qRank] */
        int tsnk,
        int n_qext,
        const int *qext,             /* [n_qext][qRank] */
        int time_rev,
        int t_axis,                  /* 0-based */
        double bc_baryon_t);

int init_extras(lua_State *L);
int fini_extras(lua_State *L);

#endif /* !defined(MARK_E390D5B2_478D_48BF_B753_64072003402B) */
