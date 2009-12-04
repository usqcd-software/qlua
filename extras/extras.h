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

const char *
gen_laplace_V(lua_State *L,
              QDP_ColorVector *result,
              double a,
              double b,
              QDP_ColorMatrix **gauge,
              QDP_ColorVector *source,
              int skip_axis);


const char *
gen_laplace_M(lua_State *L,
              QDP_ColorMatrix *result,
              double a,
              double b,
              QDP_ColorMatrix **gauge,
              QDP_ColorMatrix *source,
              int skip_axis);

const char *
gen_laplace_D(lua_State *L,
              QDP_DiracFermion *result,
              double a,
              double b,
              QDP_ColorMatrix **gauge,
              QDP_DiracFermion *source,
              int skip_axis);

const char *
gen_laplace_P(lua_State *L,
              QDP_DiracPropagator *result,
              double a,
              double b,
              QDP_ColorMatrix **gauge,
              QDP_DiracPropagator *source,
              int skip_axis);

int init_extras(lua_State *L);
int fini_extras(lua_State *L);

#endif /* !defined(MARK_E390D5B2_478D_48BF_B753_64072003402B) */
