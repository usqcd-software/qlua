#include "modules.h"                                                 /* DEPS */
#include "qlua.h"                                                    /* DEPS */
#include "lattice.h"                                                 /* DEPS */
#include "qquda.h"                                                   /* DEPS */
#include "quda.h"
#include <string.h>

static const char qudalib[] = "quda";
static const char mtnGaugeParam[] = "quda.mtGaugeParam";
static const char mtnInvertParam[] = "quda.mtInvertParam";
static const char mtnEigParam[] = "quda.mtEigParam";

/**** QudaGaugeParam */
static QudaGaugeParam *
qq_checkGaugeParam(lua_State *L, int idx)
{
  void *v = luaL_checkudata(L, idx, mtnGaugeParam);
  luaL_argcheck(L, v != 0, idx, "quda.GaugeParam expected");
  return v;
}

static int
qq_gauge_param_fmt(lua_State *L)
{
  QudaGaugeParam *p = qq_checkGaugeParam(L, 1);
  lua_pushfstring(L, "quda.GaugeParam(%p)", p);
  return 1;
}

#define QQ_ENUM(t) static void qq_named_##t(lua_State *L, int idx, t *ptr, const char *name)

QQ_ENUM(QudaReconstructType)
{
  const char *val = qlua_tabkey_stringopt(L, idx, name, NULL);
  if (val == NULL) return;
#define CHECK_VALUE(str_val, t_val) if (strcmp(val, str_val) == 0) { *ptr = t_val; return; }
  CHECK_VALUE("RECONSTRUCT_NO", QUDA_RECONSTRUCT_NO);
  CHECK_VALUE("RECONSTRUCT_12", QUDA_RECONSTRUCT_12);
  CHECK_VALUE("RECONSTRUCT_8", QUDA_RECONSTRUCT_8);
  CHECK_VALUE("RECONSTRUCT_9", QUDA_RECONSTRUCT_9);
  CHECK_VALUE("RECONSTRUCT_13", QUDA_RECONSTRUCT_13);
  CHECK_VALUE("RECONSTRUCT_10", QUDA_RECONSTRUCT_10);
#undef CHECK_VALUE
  luaL_error(L, "unexpected value for quda.GaugeParam.%s", name);
}

QQ_ENUM(QudaPrecision)
{
  const char *val = qlua_tabkey_stringopt(L, idx, name, NULL);
  if (val == NULL) return;
#define CHECK_VALUE(str_val, t_val) if (strcmp(val, str_val) == 0) { *ptr = t_val; return; }
  CHECK_VALUE("HALF_PRECISION", QUDA_HALF_PRECISION);
  CHECK_VALUE("SINGLE_PRECISION", QUDA_SINGLE_PRECISION);
  CHECK_VALUE("DOUBLE_PRECISION", QUDA_DOUBLE_PRECISION);
#undef CHECK_VALUE
  luaL_error(L, "unexpected value for quda.GaugeParam.%s", name);
}


QQ_ENUM(QudaFieldLocation)
{
  const char *val = qlua_tabkey_stringopt(L, idx, name, NULL);
  if (val == NULL) return;
#define CHECK_VALUE(str_val, t_val) if (strcmp(val, str_val) == 0) { *ptr = t_val; return; }
  CHECK_VALUE("CPU_FIELD_LOCATION", QUDA_CPU_FIELD_LOCATION);
  CHECK_VALUE("CUDA_FIELD_LOCATION", QUDA_CUDA_FIELD_LOCATION);
#undef CHECK_VALUE
  luaL_error(L, "unexpected value for quda.GaugeParam.%s", name);
}

QQ_ENUM(QudaGaugeFieldOrder)
{
  const char *val = qlua_tabkey_stringopt(L, idx, name, NULL);
  if (val == NULL) return;
#define CHECK_VALUE(str_val, t_val) if (strcmp(val, str_val) == 0) { *ptr = t_val; return; }
  CHECK_VALUE("FLOAT_GAUGE_ORDER", QUDA_FLOAT_GAUGE_ORDER);
  CHECK_VALUE("FLOAT2_GAUGE_ORDER", QUDA_FLOAT2_GAUGE_ORDER);
  CHECK_VALUE("FLOAT4_GAUGE_ORDER", QUDA_FLOAT4_GAUGE_ORDER);
  CHECK_VALUE("QDP_GAUGE_ORDER,", QUDA_QDP_GAUGE_ORDER);
  CHECK_VALUE("QDPJIT_GAUGE_ORDER,", QUDA_QDPJIT_GAUGE_ORDER);
  CHECK_VALUE("CPS_WILSON_GAUGE_ORDER,", QUDA_CPS_WILSON_GAUGE_ORDER);
  CHECK_VALUE("MILC_GAUGE_ORDER,", QUDA_MILC_GAUGE_ORDER);
  CHECK_VALUE("BQCD_GAUGE_ORDER,", QUDA_BQCD_GAUGE_ORDER);
  CHECK_VALUE("TIFR_GAUGE_ORDER,", QUDA_TIFR_GAUGE_ORDER);
#undef CHECK_VALUE
  luaL_error(L, "unexpected value for quda.GaugeParam.%s", name);
}

QQ_ENUM(QudaTboundary)
{
  const char *val = qlua_tabkey_stringopt(L, idx, name, NULL);
  if (val == NULL) return;
#define CHECK_VALUE(str_val, t_val) if (strcmp(val, str_val) == 0) { *ptr = t_val; return; }
  CHECK_VALUE("ANTI_PERIODIC_T", QUDA_ANTI_PERIODIC_T);
  CHECK_VALUE("PERIODIC_T", QUDA_PERIODIC_T);
#undef CHECK_VALUE
  luaL_error(L, "unexpected value for quda.GaugeParam.%s", name);
}

QQ_ENUM(QudaGaugeFixed)
{
  const char *val = qlua_tabkey_stringopt(L, idx, name, NULL);
  if (val == NULL) return;
#define CHECK_VALUE(str_val, t_val) if (strcmp(val, str_val) == 0) { *ptr = t_val; return; }
  CHECK_VALUE("GAUGE_FIXED_NO", QUDA_GAUGE_FIXED_NO);
  CHECK_VALUE("GAUGE_FIXED_YES", QUDA_GAUGE_FIXED_YES);
#undef CHECK_VALUE
  luaL_error(L, "unexpected value for quda.GaugeParam.%s", name);
}

QQ_ENUM(QudaLinkType)
{
  const char *val = qlua_tabkey_stringopt(L, idx, name, NULL);
  if (val == NULL) return;
#define CHECK_VALUE(str_val, t_val) if (strcmp(val, str_val) == 0) { *ptr = t_val; return; }
  CHECK_VALUE("SU3_LINKS", QUDA_SU3_LINKS);
  CHECK_VALUE("GENERAL_LINKS", QUDA_GENERAL_LINKS);
  CHECK_VALUE("THREE_LINKS", QUDA_THREE_LINKS);
  CHECK_VALUE("MOMENTUM", QUDA_MOMENTUM);
  CHECK_VALUE("WILSON_LINKS", QUDA_WILSON_LINKS);
  CHECK_VALUE("ASQTAD_FAT_LINKS", QUDA_ASQTAD_FAT_LINKS);
  CHECK_VALUE("ASQTAD_LONG_LINKS", QUDA_ASQTAD_LONG_LINKS);
  CHECK_VALUE("ASQTAD_MOM_LINKS", QUDA_ASQTAD_MOM_LINKS);
  CHECK_VALUE("ASQTAD_GENERAL_LINKS", QUDA_ASQTAD_GENERAL_LINKS);
#undef CHECK_VALUE
  luaL_error(L, "unexpected value for quda.GaugeParam.%s", name);
}

QQ_ENUM(QudaStaggeredPhase)
{
  const char *val = qlua_tabkey_stringopt(L, idx, name, NULL);
  if (val == NULL) return;
#define CHECK_VALUE(str_val, t_val) if (strcmp(val, str_val) == 0) { *ptr = t_val; return; }
  CHECK_VALUE("MILC_STAGGERED_PHASE", QUDA_MILC_STAGGERED_PHASE);
  CHECK_VALUE("CPS_STAGGERED_PHASE", QUDA_CPS_STAGGERED_PHASE);
  CHECK_VALUE("TIFR_STAGGERED_PHASE", QUDA_TIFR_STAGGERED_PHASE);
#undef CHECK_VALUE
  luaL_error(L, "unexpected value for quda.GaugeParam.%s", name);
}
#undef QQ_ENUM

static int
qq_gauge_param_set(lua_State *L)
{
  QudaGaugeParam *p = qq_checkGaugeParam(L, 1);
#define GET_INT_VALUE(name)  p->name = qlua_tabkey_intopt(L, 2, #name, p->name)
#define GET_DOUBLE_VALUE(name)  p->name = qlua_tabkey_doubleopt(L, 2, #name, p->name)
#define GET_NAMED_VALUE(t,n) qq_named_##t(L, 2, &p->n, #n)
  GET_INT_VALUE(ga_pad);
  GET_INT_VALUE(site_ga_pad);
  GET_INT_VALUE(staple_pad);
  GET_INT_VALUE(llfat_ga_pad);
  GET_INT_VALUE(mom_ga_pad);
  GET_INT_VALUE(preserve_gauge);
  GET_INT_VALUE(staggered_phase_applied);
  GET_INT_VALUE(overlap);
  GET_INT_VALUE(overwrite_mom);
  GET_INT_VALUE(use_resident_gauge);
  GET_INT_VALUE(use_resident_mom);
  GET_INT_VALUE(make_resident_gauge);
  GET_INT_VALUE(make_resident_mom);
  GET_INT_VALUE(return_result_gauge);
  GET_INT_VALUE(return_result_mom);
  GET_DOUBLE_VALUE(anisotropy);
  GET_DOUBLE_VALUE(tadpole_coeff);
  GET_DOUBLE_VALUE(scale);
  GET_DOUBLE_VALUE(gaugeGiB);
  GET_DOUBLE_VALUE(i_mu);
  GET_NAMED_VALUE(QudaReconstructType, reconstruct);
  GET_NAMED_VALUE(QudaReconstructType, reconstruct_precondition);
  GET_NAMED_VALUE(QudaReconstructType, reconstruct_sloppy);
  GET_NAMED_VALUE(QudaPrecision, cpu_prec);
  GET_NAMED_VALUE(QudaPrecision, cuda_prec);
  GET_NAMED_VALUE(QudaPrecision, cuda_prec_precondition);
  GET_NAMED_VALUE(QudaPrecision, cuda_prec_sloppy);
  GET_NAMED_VALUE(QudaFieldLocation, location);
  GET_NAMED_VALUE(QudaGaugeFieldOrder, gauge_order);
  GET_NAMED_VALUE(QudaGaugeFixed, gauge_fix);
  GET_NAMED_VALUE(QudaLinkType, type);
  GET_NAMED_VALUE(QudaStaggeredPhase, staggered_phase_type);
  GET_NAMED_VALUE(QudaTboundary, t_boundary);
#undef GET_NAMED_VALUE
#undef GET_DOUBLE_VALUE
#undef GET_INT_VALUE
  if (qlua_tabkey_tableopt(L, 2, "X")) {
    int *xx = qlua_checkintarray(L, -1, 4, NULL);
    int i;
    for (i = 0; i < 4; i++)
      p->X[i] = xx[i];
    qlua_free(L, xx);
    lua_pop(L, 1);
  }
  return 0;
}

/* XXX */
static struct luaL_Reg mtGaugeParam[] = {
  { "__tostring", qq_gauge_param_fmt },
  { "set",        qq_gauge_param_set },
  { NULL, NULL }
};

static int
qq_gauge_param(lua_State *L)
{
  QudaGaugeParam *v = lua_newuserdata(L, sizeof (QudaGaugeParam));
  luaL_getmetatable(L, mtnGaugeParam);
  lua_setmetatable(L, -2);
  *v = newQudaGaugeParam();

  return 1;
}

static struct luaL_Reg mtInvertParam[] = {
  { NULL, NULL }
};

static int
qq_invert_param(lua_State *L)
{
  /* XXX */
  return 0;
}

static struct luaL_Reg mtEigParam[] = {
  { NULL, NULL }
};

static int
qq_eig_param(lua_State *L)
{
  /* XXX */
  return 0;
}

static struct luaL_Reg fquda[] = {
  {"GaugeParam",    qq_gauge_param },
  {"InvertParam",   qq_invert_param },
  {"EigParam",      qq_eig_param },
  { NULL,     NULL }
};

int
init_quda(lua_State *L)
{
  luaL_register(L, qudalib, fquda);
  qlua_metatable(L, mtnGaugeParam, mtGaugeParam, qQudaGaugeParam);
  qlua_metatable(L, mtnInvertParam, mtInvertParam, qQudaInvertParam);
  qlua_metatable(L, mtnEigParam, mtEigParam, qQudaEigParam);
  return 0;
}

void
fini_quda(void)
{
  /* free resources of GPU */
#if 0 /* XXX */
  if (in_quda_p)
    qudaFinalize();
  in_quda_p = 0;
#endif /* XXX */
}
