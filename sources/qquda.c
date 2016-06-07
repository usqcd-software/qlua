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

#define QQ_ENUM(t) static void qq_named_##t(lua_State *L, int idx, t *ptr, const char *name, const char *tp)

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
  luaL_error(L, "unexpected value for quda.%s.%s", tp, name);
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
  luaL_error(L, "unexpected value for quda.%s.%s", tp, name);
}


QQ_ENUM(QudaFieldLocation)
{
  const char *val = qlua_tabkey_stringopt(L, idx, name, NULL);
  if (val == NULL) return;
#define CHECK_VALUE(str_val, t_val) if (strcmp(val, str_val) == 0) { *ptr = t_val; return; }
  CHECK_VALUE("CPU_FIELD_LOCATION", QUDA_CPU_FIELD_LOCATION);
  CHECK_VALUE("CUDA_FIELD_LOCATION", QUDA_CUDA_FIELD_LOCATION);
#undef CHECK_VALUE
  luaL_error(L, "unexpected value for quda.%s.%s", tp, name);
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
  luaL_error(L, "unexpected value for quda.%s.%s", tp, name);
}

QQ_ENUM(QudaTboundary)
{
  const char *val = qlua_tabkey_stringopt(L, idx, name, NULL);
  if (val == NULL) return;
#define CHECK_VALUE(str_val, t_val) if (strcmp(val, str_val) == 0) { *ptr = t_val; return; }
  CHECK_VALUE("ANTI_PERIODIC_T", QUDA_ANTI_PERIODIC_T);
  CHECK_VALUE("PERIODIC_T", QUDA_PERIODIC_T);
#undef CHECK_VALUE
  luaL_error(L, "unexpected value for quda.%s.%s", tp, name);
}

QQ_ENUM(QudaGaugeFixed)
{
  const char *val = qlua_tabkey_stringopt(L, idx, name, NULL);
  if (val == NULL) return;
#define CHECK_VALUE(str_val, t_val) if (strcmp(val, str_val) == 0) { *ptr = t_val; return; }
  CHECK_VALUE("GAUGE_FIXED_NO", QUDA_GAUGE_FIXED_NO);
  CHECK_VALUE("GAUGE_FIXED_YES", QUDA_GAUGE_FIXED_YES);
#undef CHECK_VALUE
  luaL_error(L, "unexpected value for quda.%s.%s", tp, name);
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
  luaL_error(L, "unexpected value for quda.%s.%s", tp, name);
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
  luaL_error(L, "unexpected value for quda.%s.%s", tp, name);
}

QQ_ENUM(QudaInverterType)
{
  const char *val = qlua_tabkey_stringopt(L, idx, name, NULL);
  if (val == NULL) return;
#define CHECK_VALUE(str_val, t_val) if (strcmp(val, str_val) == 0) { *ptr = t_val; return; }
  CHECK_VALUE("CG_INVERTER", QUDA_CG_INVERTER);
  CHECK_VALUE("BICGSTAB_INVERTER", QUDA_BICGSTAB_INVERTER);
  CHECK_VALUE("GCR_INVERTER", QUDA_GCR_INVERTER);
  CHECK_VALUE("MR_INVERTER", QUDA_MR_INVERTER);
  CHECK_VALUE("MPBICGSTAB_INVERTER", QUDA_MPBICGSTAB_INVERTER);
  CHECK_VALUE("SD_INVERTER", QUDA_SD_INVERTER);
  CHECK_VALUE("XSD_INVERTER", QUDA_XSD_INVERTER);
  CHECK_VALUE("PCG_INVERTER", QUDA_PCG_INVERTER);
  CHECK_VALUE("MPCG_INVERTER", QUDA_MPCG_INVERTER);
  CHECK_VALUE("EIGCG_INVERTER", QUDA_EIGCG_INVERTER);
  CHECK_VALUE("INC_EIGCG_INVERTER", QUDA_INC_EIGCG_INVERTER);
  CHECK_VALUE("GMRESDR_INVERTER", QUDA_GMRESDR_INVERTER);
  CHECK_VALUE("GMRESDR_PROJ_INVERTER", QUDA_GMRESDR_PROJ_INVERTER);
  CHECK_VALUE("GMRESDR_SH_INVERTER", QUDA_GMRESDR_SH_INVERTER);
  CHECK_VALUE("FGMRESDR_INVERTER", QUDA_FGMRESDR_INVERTER);
#undef CHECK_VALUE
  luaL_error(L, "unexpected value for quda.%s.%s", tp, name);
}

QQ_ENUM(QudaDslashType)
{
  const char *val = qlua_tabkey_stringopt(L, idx, name, NULL);
  if (val == NULL) return;
#define CHECK_VALUE(str_val, t_val) if (strcmp(val, str_val) == 0) { *ptr = t_val; return; }
  CHECK_VALUE("WILSON_DSLASH", QUDA_WILSON_DSLASH);
  CHECK_VALUE("CLOVER_WILSON_DSLASH", QUDA_CLOVER_WILSON_DSLASH);
  CHECK_VALUE("DOMAIN_WALL_DSLASH", QUDA_DOMAIN_WALL_DSLASH);
  CHECK_VALUE("DOMAIN_WALL_4D_DSLASH", QUDA_DOMAIN_WALL_4D_DSLASH);
  CHECK_VALUE("MOBIUS_DWF_DSLASH", QUDA_MOBIUS_DWF_DSLASH);
  CHECK_VALUE("STAGGERED_DSLASH", QUDA_STAGGERED_DSLASH);
  CHECK_VALUE("ASQTAD_DSLASH", QUDA_ASQTAD_DSLASH);
  CHECK_VALUE("TWISTED_MASS_DSLASH", QUDA_TWISTED_MASS_DSLASH);
  CHECK_VALUE("TWISTED_CLOVER_DSLASH", QUDA_TWISTED_CLOVER_DSLASH);
#undef CHECK_VALUE
  luaL_error(L, "unexpected value for quda.%s.%s", tp, name);
}

QQ_ENUM(QudaGammaBasis)
{
  const char *val = qlua_tabkey_stringopt(L, idx, name, NULL);
  if (val == NULL) return;
#define CHECK_VALUE(str_val, t_val) if (strcmp(val, str_val) == 0) { *ptr = t_val; return; }
  CHECK_VALUE("DEGRAND_ROSSI_GAMMA_BASIS", QUDA_DEGRAND_ROSSI_GAMMA_BASIS);
  CHECK_VALUE("UKQCD_GAMMA_BASIS", QUDA_UKQCD_GAMMA_BASIS);
  CHECK_VALUE("CHIRAL_GAMMA_BASIS", QUDA_CHIRAL_GAMMA_BASIS);
#undef CHECK_VALUE
  luaL_error(L, "unexpected value for quda.%s.%s", tp, name);
}

QQ_ENUM(QudaMassNormalization)
{
  const char *val = qlua_tabkey_stringopt(L, idx, name, NULL);
  if (val == NULL) return;
#define CHECK_VALUE(str_val, t_val) if (strcmp(val, str_val) == 0) { *ptr = t_val; return; }
  CHECK_VALUE("KAPPA_NORMALIZATION", QUDA_KAPPA_NORMALIZATION);
  CHECK_VALUE("MASS_NORMALIZATION", QUDA_MASS_NORMALIZATION);
  CHECK_VALUE("ASYMMETRIC_MASS_NORMALIZATION", QUDA_ASYMMETRIC_MASS_NORMALIZATION);
#undef CHECK_VALUE
  luaL_error(L, "unexpected value for quda.%s.%s", tp, name);
}

QQ_ENUM(QudaSolverNormalization)
{
  const char *val = qlua_tabkey_stringopt(L, idx, name, NULL);
  if (val == NULL) return;
#define CHECK_VALUE(str_val, t_val) if (strcmp(val, str_val) == 0) { *ptr = t_val; return; }
  CHECK_VALUE("DEFAULT_NORMALIZATION", QUDA_DEFAULT_NORMALIZATION);
  CHECK_VALUE("SOURCE_NORMALIZATION", QUDA_SOURCE_NORMALIZATION);
#undef CHECK_VALUE
  luaL_error(L, "unexpected value for quda.%s.%s", tp, name);
}

QQ_ENUM(QudaCloverFieldOrder)
{
  const char *val = qlua_tabkey_stringopt(L, idx, name, NULL);
  if (val == NULL) return;
#define CHECK_VALUE(str_val, t_val) if (strcmp(val, str_val) == 0) { *ptr = t_val; return; }
  CHECK_VALUE("FLOAT_CLOVER_ORDER", QUDA_FLOAT_CLOVER_ORDER);
  CHECK_VALUE("FLOAT2_CLOVER_ORDER", QUDA_FLOAT2_CLOVER_ORDER);
  CHECK_VALUE("FLOAT4_CLOVER_ORDER", QUDA_FLOAT4_CLOVER_ORDER);
  CHECK_VALUE("PACKED_CLOVER_ORDER", QUDA_PACKED_CLOVER_ORDER);
  CHECK_VALUE("QDPJIT_CLOVER_ORDER", QUDA_QDPJIT_CLOVER_ORDER);
  CHECK_VALUE("BQCD_CLOVER_ORDER", QUDA_BQCD_CLOVER_ORDER);
#undef CHECK_VALUE
  luaL_error(L, "unexpected value for quda.%s.%s", tp, name);
}

QQ_ENUM(QudaVerbosity)
{
  const char *val = qlua_tabkey_stringopt(L, idx, name, NULL);
  if (val == NULL) return;
#define CHECK_VALUE(str_val, t_val) if (strcmp(val, str_val) == 0) { *ptr = t_val; return; }
  CHECK_VALUE("SILENT", QUDA_SILENT);
  CHECK_VALUE("SUMMARIZE", QUDA_SUMMARIZE);
  CHECK_VALUE("VERBOSE", QUDA_VERBOSE);
  CHECK_VALUE("DEBUG_VERBOSE", QUDA_DEBUG_VERBOSE);
#undef CHECK_VALUE
  luaL_error(L, "unexpected value for quda.%s.%s", tp, name);
}

QQ_ENUM(QudaTune)
{
  const char *val = qlua_tabkey_stringopt(L, idx, name, NULL);
  if (val == NULL) return;
#define CHECK_VALUE(str_val, t_val) if (strcmp(val, str_val) == 0) { *ptr = t_val; return; }
  CHECK_VALUE("TUNE_NO", QUDA_TUNE_NO);
  CHECK_VALUE("TUNE_YES", QUDA_TUNE_YES);
#undef CHECK_VALUE
  luaL_error(L, "unexpected value for quda.%s.%s", tp, name);
}

QQ_ENUM(QudaDagType)
{
  const char *val = qlua_tabkey_stringopt(L, idx, name, NULL);
  if (val == NULL) return;
#define CHECK_VALUE(str_val, t_val) if (strcmp(val, str_val) == 0) { *ptr = t_val; return; }
  CHECK_VALUE("DAG_NO", QUDA_DAG_NO);
  CHECK_VALUE("DAG_YES", QUDA_DAG_YES);
#undef CHECK_VALUE
  luaL_error(L, "unexpected value for quda.%s.%s", tp, name);
}
  
QQ_ENUM(QudaDiracFieldOrder)
{
  const char *val = qlua_tabkey_stringopt(L, idx, name, NULL);
  if (val == NULL) return;
#define CHECK_VALUE(str_val, t_val) if (strcmp(val, str_val) == 0) { *ptr = t_val; return; }
  CHECK_VALUE("INTERNAL_DIRAC_ORDER", QUDA_INTERNAL_DIRAC_ORDER);
  CHECK_VALUE("DIRAC_ORDER", QUDA_DIRAC_ORDER);
  CHECK_VALUE("QDP_DIRAC_ORDER", QUDA_QDP_DIRAC_ORDER);
  CHECK_VALUE("QDPJIT_DIRAC_ORDER", QUDA_QDPJIT_DIRAC_ORDER);
  CHECK_VALUE("CPS_WILSON_DIRAC_ORDER", QUDA_CPS_WILSON_DIRAC_ORDER);
  CHECK_VALUE("LEX_DIRAC_ORDER", QUDA_LEX_DIRAC_ORDER);
#undef CHECK_VALUE
  luaL_error(L, "unexpected value for quda.%s.%s", tp, name);
}

QQ_ENUM(QudaMatPCType)
{
  const char *val = qlua_tabkey_stringopt(L, idx, name, NULL);
  if (val == NULL) return;
#define CHECK_VALUE(str_val, t_val) if (strcmp(val, str_val) == 0) { *ptr = t_val; return; }
  CHECK_VALUE("MATPC_EVEN_EVEN", QUDA_MATPC_EVEN_EVEN);
  CHECK_VALUE("MATPC_ODD_ODD", QUDA_MATPC_ODD_ODD);
  CHECK_VALUE("MATPC_EVEN_EVEN_ASYMMETRIC", QUDA_MATPC_EVEN_EVEN_ASYMMETRIC);
  CHECK_VALUE("MATPC_ODD_ODD_ASYMMETRIC", QUDA_MATPC_ODD_ODD_ASYMMETRIC);
#undef CHECK_VALUE
  luaL_error(L, "unexpected value for quda.%s.%s", tp, name);
}

QQ_ENUM(QudaPreserveSource)
{
  const char *val = qlua_tabkey_stringopt(L, idx, name, NULL);
  if (val == NULL) return;
#define CHECK_VALUE(str_val, t_val) if (strcmp(val, str_val) == 0) { *ptr = t_val; return; }
  CHECK_VALUE("PRESERVE_SOURCE_NO", QUDA_PRESERVE_SOURCE_NO);
  CHECK_VALUE("PRESERVE_SOURCE_YES", QUDA_PRESERVE_SOURCE_YES);
#undef CHECK_VALUE
  luaL_error(L, "unexpected value for quda.%s.%s", tp, name);
}

QQ_ENUM(QudaResidualType)
{
  const char *val = qlua_tabkey_stringopt(L, idx, name, NULL);
  if (val == NULL) return;
#define CHECK_VALUE(str_val, t_val) if (strcmp(val, str_val) == 0) { *ptr = t_val; return; }
  CHECK_VALUE("L2_RELATIVE_RESIDUAL", QUDA_L2_RELATIVE_RESIDUAL);
  CHECK_VALUE("L2_ABSOLUTE_RESIDUAL", QUDA_L2_ABSOLUTE_RESIDUAL);
  CHECK_VALUE("HEAVY_QUARK_RESIDUAL", QUDA_HEAVY_QUARK_RESIDUAL);
#undef CHECK_VALUE
  luaL_error(L, "unexpected value for quda.%s.%s", tp, name);
}

QQ_ENUM(QudaSchwarzType)
{
  const char *val = qlua_tabkey_stringopt(L, idx, name, NULL);
  if (val == NULL) return;
#define CHECK_VALUE(str_val, t_val) if (strcmp(val, str_val) == 0) { *ptr = t_val; return; }
  CHECK_VALUE("ADDITIVE_SCHWARZ", QUDA_ADDITIVE_SCHWARZ);
  CHECK_VALUE("MULTIPLICATIVE_SCHWARZ", QUDA_MULTIPLICATIVE_SCHWARZ);
#undef CHECK_VALUE
  luaL_error(L, "unexpected value for quda.%s.%s", tp, name);
}

QQ_ENUM(QudaSolutionType)
{
  const char *val = qlua_tabkey_stringopt(L, idx, name, NULL);
  if (val == NULL) return;
#define CHECK_VALUE(str_val, t_val) if (strcmp(val, str_val) == 0) { *ptr = t_val; return; }
  CHECK_VALUE("MAT_SOLUTION", QUDA_MAT_SOLUTION);
  CHECK_VALUE("MATDAG_MAT_SOLUTION", QUDA_MATDAG_MAT_SOLUTION);
  CHECK_VALUE("MATPC_SOLUTION", QUDA_MATPC_SOLUTION);
  CHECK_VALUE("MATPC_DAG_SOLUTION", QUDA_MATPC_DAG_SOLUTION);
  CHECK_VALUE("MATPCDAG_MATPC_SOLUTION", QUDA_MATPCDAG_MATPC_SOLUTION);
  CHECK_VALUE("MATPCDAG_MATPC_SHIFT_SOLUTION", QUDA_MATPCDAG_MATPC_SHIFT_SOLUTION);
#undef CHECK_VALUE
  luaL_error(L, "unexpected value for quda.%s.%s", tp, name);
}

QQ_ENUM(QudaSolveType)
{
  const char *val = qlua_tabkey_stringopt(L, idx, name, NULL);
  if (val == NULL) return;
#define CHECK_VALUE(str_val, t_val) if (strcmp(val, str_val) == 0) { *ptr = t_val; return; }
  CHECK_VALUE("DIRECT_SOLVE", QUDA_DIRECT_SOLVE);
  CHECK_VALUE("NORMOP_SOLVE", QUDA_NORMOP_SOLVE);
  CHECK_VALUE("DIRECT_PC_SOLVE", QUDA_DIRECT_PC_SOLVE);
  CHECK_VALUE("NORMOP_PC_SOLVE", QUDA_NORMOP_PC_SOLVE);
  CHECK_VALUE("NORMERR_SOLVE", QUDA_NORMERR_SOLVE);
  CHECK_VALUE("NORMERR_PC_SOLVE", QUDA_NORMERR_PC_SOLVE);
  CHECK_VALUE("NORMEQ_SOLVE", QUDA_NORMEQ_SOLVE);
  CHECK_VALUE("NORMEQ_PC_SOLVE", QUDA_NORMEQ_PC_SOLVE);
#undef CHECK_VALUE
  luaL_error(L, "unexpected value for quda.%s.%s", tp, name);
}

QQ_ENUM(QudaTwistFlavorType)
{
  const char *val = qlua_tabkey_stringopt(L, idx, name, NULL);
  if (val == NULL) return;
#define CHECK_VALUE(str_val, t_val) if (strcmp(val, str_val) == 0) { *ptr = t_val; return; }
  CHECK_VALUE("TWIST_MINUS", QUDA_TWIST_MINUS);
  CHECK_VALUE("TWIST_PLUS", QUDA_TWIST_PLUS);
  CHECK_VALUE("TWIST_NONDEG_DOUBLET", QUDA_TWIST_NONDEG_DOUBLET);
  CHECK_VALUE("TWIST_DEG_DOUBLET", QUDA_TWIST_DEG_DOUBLET);
  CHECK_VALUE("TWIST_NO", QUDA_TWIST_NO);
#undef CHECK_VALUE
  luaL_error(L, "unexpected value for quda.%s.%s", tp, name);
}

QQ_ENUM(QudaUseInitGuess)
{
  const char *val = qlua_tabkey_stringopt(L, idx, name, NULL);
  if (val == NULL) return;
#define CHECK_VALUE(str_val, t_val) if (strcmp(val, str_val) == 0) { *ptr = t_val; return; }
  CHECK_VALUE("USE_INIT_GUESS_NO", QUDA_USE_INIT_GUESS_NO);
  CHECK_VALUE("USE_INIT_GUESS_YES", QUDA_USE_INIT_GUESS_YES);
#undef CHECK_VALUE
  luaL_error(L, "unexpected value for quda.%s.%s", tp, name);
}
#undef QQ_ENUM

static int
qq_gauge_param_set(lua_State *L)
{
  QudaGaugeParam *p = qq_checkGaugeParam(L, 1);
#define GET_INT_VALUE(name)  p->name = qlua_tabkey_intopt(L, 2, #name, p->name)
#define GET_DOUBLE_VALUE(name)  p->name = qlua_tabkey_doubleopt(L, 2, #name, p->name)
#define GET_NAMED_VALUE(t,n) qq_named_##t(L, 2, &p->n, #n, "GaugeParam")
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

static int
qq_gauge_param_print(lua_State *L)
{
  QudaGaugeParam *p = qq_checkGaugeParam(L, 1);
  printQudaGaugeParam(p);
  return 0;
}

static struct luaL_Reg mtGaugeParam[] = {
  { "__tostring", qq_gauge_param_fmt },
  { "set",        qq_gauge_param_set },
  { "print",      qq_gauge_param_print },
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

/**** QudaInvertParam */
static QudaInvertParam *
qq_checkInvertParam(lua_State *L, int idx)
{
  void *v = luaL_checkudata(L, idx, mtnInvertParam);
  luaL_argcheck(L, v != 0, idx, "quda.InvertParam expected");
  return v;
}

static int
qq_invert_param_fmt(lua_State *L)
{
  QudaInvertParam *p = qq_checkInvertParam(L, 1);
  lua_pushfstring(L, "quda.InvertParam(%p)", p);
  return 1;
}

static int
qq_invert_param_set(lua_State *L)
{
  QudaInvertParam *p = qq_checkInvertParam(L, 1);
#define GET_INT_VALUE(name)  p->name = qlua_tabkey_intopt(L, 2, #name, p->name)
#define GET_DOUBLE_VALUE(name)  p->name = qlua_tabkey_doubleopt(L, 2, #name, p->name)
#define GET_NAMED_VALUE(t,n) qq_named_##t(L, 2, &p->n, #n, "InvertParam")

  GET_NAMED_VALUE(QudaCloverFieldOrder, clover_order);
  GET_NAMED_VALUE(QudaDagType, dagger);
  GET_NAMED_VALUE(QudaDiracFieldOrder, dirac_order);
  GET_NAMED_VALUE(QudaDslashType, dslash_type);
  GET_NAMED_VALUE(QudaDslashType, dslash_type_precondition);
  GET_NAMED_VALUE(QudaFieldLocation, clover_location);
  GET_NAMED_VALUE(QudaFieldLocation, input_location);
  GET_NAMED_VALUE(QudaFieldLocation, output_location);
  GET_NAMED_VALUE(QudaGammaBasis, gamma_basis);
  GET_NAMED_VALUE(QudaInverterType, inv_type);
  GET_NAMED_VALUE(QudaInverterType, inv_type_precondition);
  GET_NAMED_VALUE(QudaMassNormalization, mass_normalization);
  GET_NAMED_VALUE(QudaMatPCType, matpc_type);
  GET_NAMED_VALUE(QudaPrecision, clover_cpu_prec);
  GET_NAMED_VALUE(QudaPrecision, clover_cuda_prec);
  GET_NAMED_VALUE(QudaPrecision, clover_cuda_prec_precondition);
  GET_NAMED_VALUE(QudaPrecision, clover_cuda_prec_sloppy);
  GET_NAMED_VALUE(QudaPrecision, cpu_prec);
  GET_NAMED_VALUE(QudaPrecision, cuda_prec);
  GET_NAMED_VALUE(QudaPrecision, cuda_prec_precondition);
  GET_NAMED_VALUE(QudaPrecision, cuda_prec_ritz);
  GET_NAMED_VALUE(QudaPrecision, cuda_prec_sloppy);
  GET_NAMED_VALUE(QudaPreserveSource, preserve_source);
  GET_NAMED_VALUE(QudaResidualType, residual_type);
  GET_NAMED_VALUE(QudaSchwarzType, schwarz_type);
  GET_NAMED_VALUE(QudaSolutionType, solution_type);
  GET_NAMED_VALUE(QudaSolveType, solve_type);
  GET_NAMED_VALUE(QudaSolverNormalization, solver_normalization);
  GET_NAMED_VALUE(QudaTune, tune);
  GET_NAMED_VALUE(QudaTwistFlavorType, twist_flavor);
  GET_NAMED_VALUE(QudaUseInitGuess, use_init_guess);
  GET_NAMED_VALUE(QudaVerbosity, verbosity);
  GET_NAMED_VALUE(QudaVerbosity, verbosity_precondition);
  GET_DOUBLE_VALUE(cg_iterref_tol);
  GET_DOUBLE_VALUE(cloverGiB);
  GET_DOUBLE_VALUE(clover_coeff);
  GET_DOUBLE_VALUE(eigenval_tol);
  GET_DOUBLE_VALUE(epsilon);
  GET_DOUBLE_VALUE(gflops);
  GET_DOUBLE_VALUE(inc_tol);
  GET_DOUBLE_VALUE(kappa);
  GET_DOUBLE_VALUE(mass);
  GET_DOUBLE_VALUE(mu);
  GET_DOUBLE_VALUE(omega);
  GET_DOUBLE_VALUE(reliable_delta);
  GET_DOUBLE_VALUE(secs);
  GET_DOUBLE_VALUE(spinorGiB);
  GET_DOUBLE_VALUE(tol);
  GET_DOUBLE_VALUE(tol_hq);
  GET_DOUBLE_VALUE(tol_precondition);
  GET_DOUBLE_VALUE(tol_restart);
  GET_DOUBLE_VALUE(true_res);
  GET_DOUBLE_VALUE(true_res_hq);
  GET_INT_VALUE(Nsteps);
  GET_INT_VALUE(cl_pad);
  GET_INT_VALUE(deflation_grid);
  GET_INT_VALUE(eigcg_max_restarts);
  GET_INT_VALUE(gcrNkrylov);
  GET_INT_VALUE(heavy_quark_check);
  GET_INT_VALUE(iter);
  GET_INT_VALUE(make_resident_solution);
  GET_INT_VALUE(max_res_increase);
  GET_INT_VALUE(max_res_increase_total);
  GET_INT_VALUE(max_restart_num);
  GET_INT_VALUE(max_search_dim);
  GET_INT_VALUE(maxiter);
  GET_INT_VALUE(maxiter_precondition);
  GET_INT_VALUE(nev);
  GET_INT_VALUE(num_offset);
  GET_INT_VALUE(overlap);
  GET_INT_VALUE(pipeline);
  GET_INT_VALUE(precondition_cycle);
  GET_INT_VALUE(rhs_idx);
  GET_INT_VALUE(sp_pad);
  GET_INT_VALUE(use_cg_updates);
  GET_INT_VALUE(use_reduced_vector_set);
  GET_INT_VALUE(use_resident_solution);
  GET_INT_VALUE(use_sloppy_partial_accumulator);
  /* XXX no clover trace log computation */
  // int compute_clover_trlog;
  // double trlogA[2];
  /* XXX no DWF stuff yet */
  // double b_5[QUDA_MAX_DWF_LS];
  //   double c_5[QUDA_MAX_DWF_LS];
  //   double m5;
  //   int Ls;
  /* XXX no multishift solvers for now */
  //   double iter_res_offset[QUDA_MAX_MULTI_SHIFT];
  //   double offset[QUDA_MAX_MULTI_SHIFT];
  //   double tol_hq_offset[QUDA_MAX_MULTI_SHIFT];
  //   double tol_offset[QUDA_MAX_MULTI_SHIFT];     
  //   double true_res_hq_offset[QUDA_MAX_MULTI_SHIFT]; 
  //   double true_res_offset[QUDA_MAX_MULTI_SHIFT]; 

#undef GET_NAMED_VALUE
#undef GET_DOUBLE_VALUE
#undef GET_INT_VALUE
  return 0;
}


static int
qq_invert_param_print(lua_State *L)
{
  QudaInvertParam *p = qq_checkInvertParam(L, 1);
  printQudaInvertParam(p);
  return 0;
}


static struct luaL_Reg mtInvertParam[] = {
  { "__tostring", qq_invert_param_fmt },
  { "set",        qq_invert_param_set },
  { "print",      qq_invert_param_print },
  { NULL, NULL }
};

static int
qq_invert_param(lua_State *L)
{
  QudaInvertParam *v = lua_newuserdata(L, sizeof (QudaInvertParam));
  luaL_getmetatable(L, mtnInvertParam);
  lua_setmetatable(L, -2);
  *v = newQudaInvertParam();

  return 1;
}

#if 0 /* XXX no QudaEig related stuff yet */
/***** QudaEigParam */
#if 0 /* XXX */
QudaEigType eig_type;
QudaInvertParam *invert_param; /* XXX this is troublesome */
QudaSolutionType  RitzMat_Convcheck;
QudaSolutionType  RitzMat_lanczos;
double *MatPoly_param; /* XXX this is troublesome too */
double Stp_residual;
double eigen_shift;
int NPoly;
int f_size;
int nk;
int np;
#endif /* XXX */

static struct luaL_Reg mtEigParam[] = {
  { NULL, NULL }
};

static int
qq_eig_param(lua_State *L)
{
  /* XXX */
  return 0;
}
#endif /* XXX Eig-related stuff off */

#if 0 /* XXXX Quda interface functions */
void initQudaMemory();
void endQuda();
void freeGaugeQuda();
void freeCloverQuda();
double qChargeCuda();
void openMagma();
void closeMagma();

void initQudaDevice(int device);
void initQuda(int device);

void plaqQuda(double plaq[3]);

void setVerbosityQuda(QudaVerbosity verbosity, const char prefix[], FILE *outfile);
void initCommsGridQuda(int nDim, const int *dims, QudaCommsMap func, void *fdata);
void loadGaugeQuda(void *h_gauge, QudaGaugeParam *param);
void saveGaugeQuda(void *h_gauge, QudaGaugeParam *param);
void loadCloverQuda(void *h_clover, void *h_clovinv, QudaInvertParam *inv_param);
void invertQuda(void *h_x, void *h_b, QudaInvertParam *param);
void performAPEnStep(unsigned int nSteps, double alpha);
void projectSU3Quda(void *gauge_h, double tol, QudaGaugeParam *param);
void destroyDeflationQuda(QudaInvertParam *param, const int *X, void *_h_u, double *inv_eigenvals);

void dslashQuda(void *h_out, void *h_in, QudaInvertParam *inv_param, QudaParity parity);
void dslashQuda_4dpc(void *h_out, void *h_in, QudaInvertParam *inv_param, QudaParity parity, int test_type);
void dslashQuda_mdwf(void *h_out, void *h_in, QudaInvertParam *inv_param, QudaParity parity, int test_type);
void cloverQuda(void *h_out, void *h_in, QudaInvertParam *inv_param, QudaParity *parity, int inverse);
void MatQuda(void *h_out, void *h_in, QudaInvertParam *inv_param);
void MatDagMatQuda(void *h_out, void *h_in, QudaInvertParam *inv_param);
void* createExtendedGaugeFieldQuda(void* gauge, int geometry, QudaGaugeParam* param);
void* createGaugeFieldQuda(void* gauge, int geometry, QudaGaugeParam* param);
void saveGaugeFieldQuda(void* outGauge, void* inGauge, QudaGaugeParam* param);
void  extendGaugeFieldQuda(void* outGauge, void* inGauge);
void destroyGaugeFieldQuda(void* gauge);
void createCloverQuda(QudaInvertParam* param);
void computeCloverTraceQuda(void* out, void* dummy, int mu, int nu, int dim[4]);

#endif /* XXXX Quda interface functions */

static struct luaL_Reg fquda[] = {
  {"GaugeParam",    qq_gauge_param },
  {"InvertParam",   qq_invert_param },
  // XXX {"EigParam",      qq_eig_param },
  /* XXX */
  { NULL,     NULL }
};

int
init_quda(lua_State *L)
{
  luaL_register(L, qudalib, fquda);
  qlua_metatable(L, mtnGaugeParam, mtGaugeParam, qQudaGaugeParam);
  qlua_metatable(L, mtnInvertParam, mtInvertParam, qQudaInvertParam);
  // XXX qlua_metatable(L, mtnEigParam, mtEigParam, qQudaEigParam);
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
