#include "modules.h"                                                 /* DEPS */
#include "qlua.h"                                                    /* DEPS */
#include "lattice.h"                                                 /* DEPS */
#include "qlayout.h"                                                 /* DEPS */
#include "latcolmat.h"                                               /* DEPS */
#include "latdirferm.h"                                              /* DEPS */
#include "qquda.h"                                                   /* DEPS */
#include "quda.h"
#include <string.h>

#define QUDA_REAL double
#define QUDA_Nc 3
#define QUDA_Ns 4
#define QUDA_DIM 4

static const char qudalib[] = "_quda";
static const char mtnGaugeParam[] = "_quda.mtGaugeParam";
static const char mtnInvertParam[] = "_quda.mtInvertParam";
static const char mtnEigParam[] = "_quda.mtEigParam";

/**** QudaGaugeParam */
static QudaGaugeParam *
qq_checkGaugeParam(lua_State *L, int idx)
{
  void *v = luaL_checkudata(L, idx, mtnGaugeParam);
  luaL_argcheck(L, v != 0, idx, "quda.GaugeParam expected");
  return v;
}

static int
qq_gp_fmt(lua_State *L)
{
  QudaGaugeParam *p = qq_checkGaugeParam(L, 1);
  lua_pushfstring(L, "quda.GaugeParam(%p)", p);
  return 1;
}

/* conversions between QUDA enums and Qlua strings */
#define QQ_ENUM(t) static t qlua2quda_##t(lua_State *L, const char *name, const char *tp, const char *val)
#define CHECK_VALUE(val, str_val) if (strcmp(val, #str_val) == 0) { return QUDA_##str_val; }
#define DUMMY_VALUE(str_val) QUDA_##str_val
#include "qq-enum.c"   /* DEPS */

#define QQ_ENUM(t) static const char * quda2qlua_##t(lua_State *L, const char *name, const char *tp, t val)
#define CHECK_VALUE(val, str_val) if (val == QUDA_##str_val) { return #str_val; }
#define DUMMY_VALUE(str_val) #str_val
#include "qq-enum.c"   /* DEPS */

static int
qq_gp_get(lua_State *L)
{
  QudaGaugeParam *p = qq_checkGaugeParam(L, 1);
  const char *fld = luaL_checkstring(L, 2);
  
#define GET_INT_VALUE(name) if (strcmp(fld, #name) == 0) { \
    lua_pushnumber(L, p->name); return 1; }
#define GET_DOUBLE_VALUE(name) if (strcmp(fld, #name) == 0) { \
    lua_pushnumber(L, p->name); return 1; }
#define GET_NAMED_VALUE(t,name)  if (strcmp(fld, #name) == 0) { \
    lua_pushstring(L, quda2qlua_##t(L, #name, "GaugeParam", p->name)); return 1; }

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
  if (strcmp(fld, "X") == 0) {
    int i;
    lua_createtable(L, QUDA_DIM, 0);
    for (i = 0; i < QUDA_DIM; i++) {
      lua_pushinteger(L, p->X[i]);
      lua_rawseti(L, -2, i+1); /* lua indexing */
    }
    return 1;
  }
  return qlua_lookup(L, 2, mtnGaugeParam);
}

static int
qq_gp_put(lua_State *L)
{
  QudaGaugeParam *p = qq_checkGaugeParam(L, 1);
  const char *fld = luaL_checkstring(L, 2);
  
#define PUT_INT_VALUE(name) if (strcmp(fld, #name) == 0) { \
    p->name = luaL_checkint(L, 3); return 0; }
#define PUT_DOUBLE_VALUE(name) if (strcmp(fld, #name) == 0) { \
    p->name = luaL_checknumber(L, 3); return 1; }
#define PUT_NAMED_VALUE(t,name)  if (strcmp(fld, #name) == 0) { \
    p->name = qlua2quda_##t(L, #name, "GaugeParam", luaL_checkstring(L, 3)); return 0; }

  PUT_INT_VALUE(ga_pad);
  PUT_INT_VALUE(site_ga_pad);
  PUT_INT_VALUE(staple_pad);
  PUT_INT_VALUE(llfat_ga_pad);
  PUT_INT_VALUE(mom_ga_pad);
  PUT_INT_VALUE(preserve_gauge);
  PUT_INT_VALUE(staggered_phase_applied);
  PUT_INT_VALUE(overlap);
  PUT_INT_VALUE(overwrite_mom);
  PUT_INT_VALUE(use_resident_gauge);
  PUT_INT_VALUE(use_resident_mom);
  PUT_INT_VALUE(make_resident_gauge);
  PUT_INT_VALUE(make_resident_mom);
  PUT_INT_VALUE(return_result_gauge);
  PUT_INT_VALUE(return_result_mom);
  PUT_DOUBLE_VALUE(anisotropy);
  PUT_DOUBLE_VALUE(tadpole_coeff);
  PUT_DOUBLE_VALUE(scale);
  PUT_DOUBLE_VALUE(gaugeGiB);
  PUT_DOUBLE_VALUE(i_mu);
  PUT_NAMED_VALUE(QudaReconstructType, reconstruct);
  PUT_NAMED_VALUE(QudaReconstructType, reconstruct_precondition);
  PUT_NAMED_VALUE(QudaReconstructType, reconstruct_sloppy);
  PUT_NAMED_VALUE(QudaPrecision, cpu_prec);
  PUT_NAMED_VALUE(QudaPrecision, cuda_prec);
  PUT_NAMED_VALUE(QudaPrecision, cuda_prec_precondition);
  PUT_NAMED_VALUE(QudaPrecision, cuda_prec_sloppy);
  PUT_NAMED_VALUE(QudaFieldLocation, location);
  PUT_NAMED_VALUE(QudaGaugeFieldOrder, gauge_order);
  PUT_NAMED_VALUE(QudaGaugeFixed, gauge_fix);
  PUT_NAMED_VALUE(QudaLinkType, type);
  PUT_NAMED_VALUE(QudaStaggeredPhase, staggered_phase_type);
  PUT_NAMED_VALUE(QudaTboundary, t_boundary);

#undef PUT_NAMED_VALUE
#undef PUT_DOUBLE_VALUE
#undef PUT_INT_VALUE
  if (strcmp(fld, "X") == 0) {
    int *xx = qlua_checkintarray(L, -1, QUDA_DIM, NULL);
    int i;
    for (i = 0; i < QUDA_DIM; i++)
      p->X[i] = xx[i];
    qlua_free(L, xx);
    lua_pop(L, 1);
    return 0;
  }    
  luaL_error(L, "invalid or unsettable GaugeParam element");
  return 0;
}

static int
qq_gp_print(lua_State *L)
{
  QudaGaugeParam *p = qq_checkGaugeParam(L, 1);

  printQudaGaugeParam(p);

  return 0;
}

static int
qq_gp_copy(lua_State *L)
{
  QudaGaugeParam *p = qq_checkGaugeParam(L, 1);
  QudaGaugeParam *v = lua_newuserdata(L, sizeof (QudaGaugeParam));

  luaL_getmetatable(L, mtnGaugeParam);
  lua_setmetatable(L, -2);
  *v = *p;

  return 1;
}

static struct luaL_Reg mtGaugeParam[] = {
  { "__tostring", qq_gp_fmt },
  { "__index",    qq_gp_get },
  { "__newindex", qq_gp_put },
  { "copy",       qq_gp_copy },
  { "print",      qq_gp_print },
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
qq_ip_fmt(lua_State *L)
{
  QudaInvertParam *p = qq_checkInvertParam(L, 1);
  lua_pushfstring(L, "quda.InvertParam(%p)", p);
  return 1;
}

#if 0 /* XXX */
static int
qq_ip_set(lua_State *L)
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
#undef GET_NAMED_VALUE
#undef GET_DOUBLE_VALUE
#undef GET_INT_VALUE
  return 0;
}
#endif /* XXX */

static int
qq_ip_get(lua_State *L)
{
  QudaInvertParam *p = qq_checkInvertParam(L, 1);
  const char *fld = luaL_checkstring(L, 2);
  
#define GET_INT_VALUE(name) if (strcmp(fld, #name) == 0) { \
    lua_pushnumber(L, p->name); return 1; }
#define GET_DOUBLE_VALUE(name) if (strcmp(fld, #name) == 0) { \
    lua_pushnumber(L, p->name); return 1; }
#define GET_NAMED_VALUE(t,name)  if (strcmp(fld, #name) == 0) { \
    lua_pushstring(L, quda2qlua_##t(L, #name, "InvertParam", p->name)); return 1; }

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

#undef GET_NAMED_VALUE
#undef GET_DOUBLE_VALUE
#undef GET_INT_VALUE
  return qlua_lookup(L, 2, mtnInvertParam);
}

static int
qq_ip_put(lua_State *L)
{
  QudaInvertParam *p = qq_checkInvertParam(L, 1);
  const char *fld = luaL_checkstring(L, 2);
  
#define PUT_INT_VALUE(name) if (strcmp(fld, #name) == 0) { \
    p->name = luaL_checkint(L, 3); return 0; }
#define PUT_DOUBLE_VALUE(name) if (strcmp(fld, #name) == 0) { \
    p->name = luaL_checknumber(L, 3); return 1; }
#define PUT_NAMED_VALUE(t,name)  if (strcmp(fld, #name) == 0) { \
    p->name = qlua2quda_##t(L, #name, "InvertParam", luaL_checkstring(L, 3)); return 0; }

  PUT_NAMED_VALUE(QudaCloverFieldOrder, clover_order);
  PUT_NAMED_VALUE(QudaDagType, dagger);
  PUT_NAMED_VALUE(QudaDiracFieldOrder, dirac_order);
  PUT_NAMED_VALUE(QudaDslashType, dslash_type);
  PUT_NAMED_VALUE(QudaDslashType, dslash_type_precondition);
  PUT_NAMED_VALUE(QudaFieldLocation, clover_location);
  PUT_NAMED_VALUE(QudaFieldLocation, input_location);
  PUT_NAMED_VALUE(QudaFieldLocation, output_location);
  PUT_NAMED_VALUE(QudaGammaBasis, gamma_basis);
  PUT_NAMED_VALUE(QudaInverterType, inv_type);
  PUT_NAMED_VALUE(QudaInverterType, inv_type_precondition);
  PUT_NAMED_VALUE(QudaMassNormalization, mass_normalization);
  PUT_NAMED_VALUE(QudaMatPCType, matpc_type);
  PUT_NAMED_VALUE(QudaPrecision, clover_cpu_prec);
  PUT_NAMED_VALUE(QudaPrecision, clover_cuda_prec);
  PUT_NAMED_VALUE(QudaPrecision, clover_cuda_prec_precondition);
  PUT_NAMED_VALUE(QudaPrecision, clover_cuda_prec_sloppy);
  PUT_NAMED_VALUE(QudaPrecision, cpu_prec);
  PUT_NAMED_VALUE(QudaPrecision, cuda_prec);
  PUT_NAMED_VALUE(QudaPrecision, cuda_prec_precondition);
  PUT_NAMED_VALUE(QudaPrecision, cuda_prec_ritz);
  PUT_NAMED_VALUE(QudaPrecision, cuda_prec_sloppy);
  PUT_NAMED_VALUE(QudaPreserveSource, preserve_source);
  PUT_NAMED_VALUE(QudaResidualType, residual_type);
  PUT_NAMED_VALUE(QudaSchwarzType, schwarz_type);
  PUT_NAMED_VALUE(QudaSolutionType, solution_type);
  PUT_NAMED_VALUE(QudaSolveType, solve_type);
  PUT_NAMED_VALUE(QudaSolverNormalization, solver_normalization);
  PUT_NAMED_VALUE(QudaTune, tune);
  PUT_NAMED_VALUE(QudaTwistFlavorType, twist_flavor);
  PUT_NAMED_VALUE(QudaUseInitGuess, use_init_guess);
  PUT_NAMED_VALUE(QudaVerbosity, verbosity);
  PUT_NAMED_VALUE(QudaVerbosity, verbosity_precondition);
  PUT_DOUBLE_VALUE(cg_iterref_tol);
  PUT_DOUBLE_VALUE(cloverGiB);
  PUT_DOUBLE_VALUE(clover_coeff);
  PUT_DOUBLE_VALUE(eigenval_tol);
  PUT_DOUBLE_VALUE(epsilon);
  PUT_DOUBLE_VALUE(gflops);
  PUT_DOUBLE_VALUE(inc_tol);
  PUT_DOUBLE_VALUE(kappa);
  PUT_DOUBLE_VALUE(mass);
  PUT_DOUBLE_VALUE(mu);
  PUT_DOUBLE_VALUE(omega);
  PUT_DOUBLE_VALUE(reliable_delta);
  PUT_DOUBLE_VALUE(secs);
  PUT_DOUBLE_VALUE(spinorGiB);
  PUT_DOUBLE_VALUE(tol);
  PUT_DOUBLE_VALUE(tol_hq);
  PUT_DOUBLE_VALUE(tol_precondition);
  PUT_DOUBLE_VALUE(tol_restart);
  PUT_DOUBLE_VALUE(true_res);
  PUT_DOUBLE_VALUE(true_res_hq);
  PUT_INT_VALUE(Nsteps);
  PUT_INT_VALUE(cl_pad);
  PUT_INT_VALUE(deflation_grid);
  PUT_INT_VALUE(eigcg_max_restarts);
  PUT_INT_VALUE(gcrNkrylov);
  PUT_INT_VALUE(heavy_quark_check);
  PUT_INT_VALUE(iter);
  PUT_INT_VALUE(make_resident_solution);
  PUT_INT_VALUE(max_res_increase);
  PUT_INT_VALUE(max_res_increase_total);
  PUT_INT_VALUE(max_restart_num);
  PUT_INT_VALUE(max_search_dim);
  PUT_INT_VALUE(maxiter);
  PUT_INT_VALUE(maxiter_precondition);
  PUT_INT_VALUE(nev);
  PUT_INT_VALUE(num_offset);
  PUT_INT_VALUE(overlap);
  PUT_INT_VALUE(pipeline);
  PUT_INT_VALUE(precondition_cycle);
  PUT_INT_VALUE(rhs_idx);
  PUT_INT_VALUE(sp_pad);
  PUT_INT_VALUE(use_cg_updates);
  PUT_INT_VALUE(use_reduced_vector_set);
  PUT_INT_VALUE(use_resident_solution);
  PUT_INT_VALUE(use_sloppy_partial_accumulator);

#undef PUT_NAMED_VALUE
#undef PUT_DOUBLE_VALUE
#undef PUT_INT_VALUE

  luaL_error(L, "invalid or unsettable InvertParam element");
  return 0;
}

static int
qq_ip_print(lua_State *L)
{
  QudaInvertParam *p = qq_checkInvertParam(L, 1);
  printQudaInvertParam(p);
  return 0;
}

static int
qq_ip_copy(lua_State *L)
{
  QudaInvertParam *p = qq_checkInvertParam(L, 1);
  QudaInvertParam *v = lua_newuserdata(L, sizeof (QudaInvertParam));

  luaL_getmetatable(L, mtnInvertParam);
  lua_setmetatable(L, -2);
  *v = *p;

  return 1;
}


static struct luaL_Reg mtInvertParam[] = {
  { "__tostring", qq_ip_fmt },
  { "__index",    qq_ip_get },
  { "__newindex", qq_ip_put },
  { "copy",       qq_ip_copy },
  { "print",      qq_ip_print },
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

static int
qq_initQudaMemory(lua_State *L)
{
  initQudaMemory();
  return 0;
}

static int
qq_endQuda(lua_State *L)
{
  endQuda();
  return 0;
}

static int
qq_freeGaugeQuda(lua_State *L)
{
  freeGaugeQuda();
  return 0;
}

static int
qq_freeCloverQuda(lua_State *L)
{
  freeCloverQuda();
  return 0;
}

static int
qq_qChargeCuda(lua_State *L)
{
  double q = qChargeCuda();
  lua_pushnumber(L, q);
  return 1;
}

static int
qq_openMagma(lua_State *L)
{
  openMagma();
  return 0;
}

static int
qq_closeMagma(lua_State *L)
{
  closeMagma();
  return 0;
}

static int
qq_initQudaDevice(lua_State *L)
{
  int dev = lua_gettop(L) > 0? luaL_checkint(L, 1): -1;
  initQudaDevice(dev);

  return 0;
}

static int
qq_initQuda(lua_State *L)
{
  int dev = lua_gettop(L) > 0? luaL_checkint(L, 1): -1;
  initQuda(dev);

  return 0;
}

static int
qq_plaqQuda(lua_State *L)
{
  double plaq[3];
  plaqQuda(plaq);
  lua_pushnumber(L, plaq[0]);
  lua_pushnumber(L, plaq[1]);
  lua_pushnumber(L, plaq[2]);

  return 3;
}

static int
qq_setVerbosityQuda(lua_State *L)
{
  QudaVerbosity verb = qlua2quda_QudaVerbosity(L, "verbosity", "setQudaVerbosity", luaL_checkstring(L, 1));
  const char *prefix = lua_gettop(L) >= 2? luaL_checkstring(L, 2): "QUDA> ";
  setVerbosityQuda(verb, prefix, stdout);

  return 0;
}

static int
qq_initCommsGridQuda(lua_State *L)
{
  mLattice *S = qlua_checkLattice(L, 1);
  initCommsGridQuda(S->rank, S->net, qlua_comm_map, S);
  return 0;
}

static int
quda_index(const int x[], const int lo[], const int hi[])
{
  int d, k, p, v;
  for (d = QUDA_DIM, v = 1, k = 0, p = 0; d--;) {
    k *= hi[d] - lo[d];
    k += x[d] - lo[d];
    p += x[d];
    v *= hi[d] - lo[d];
  }
  k /= 2;
  if (p & 1) k += v / 2;
  return k;
}

static void
get_gauge_field(QUDA_REAL **q, QDP_D3_ColorMatrix **U, lua_State *L, int idx, int d)
{
  QLA_D3_ColorMatrix *Ux;
  mLattice *S = NULL;
  int lo[QUDA_DIM];
  int hi[QUDA_DIM];
  int subvol;
  int i;

  lua_pushnumber(L, d + 1); /* lua indexing */
  lua_gettable(L, idx);
  *U = qlua_checkLatColMat3(L, -1, NULL, QUDA_Nc)->ptr;
  S = qlua_ObjLattice(L, -1);
  qlua_assert(S->rank == QUDA_DIM, "expected rank 4 lattice");
  qlua_sublattice(lo, hi, S->node, S);
  for (i = 0, subvol = 1; i < QUDA_DIM; i++) {
    subvol *= hi[i] - lo[i];
  }
  
  *q = qlua_malloc(L, subvol * 2 * QUDA_Nc * QUDA_Nc * sizeof (QUDA_REAL));
  Ux = QDP_D3_expose_M(*U);
  for (i = 0; i < subvol; i++) {
    int a, b, ci, x[QUDA_DIM];
    QUDA_REAL *ptr;
    QDP_get_coords_L(S->lat, x, S->node, i);
    ci = quda_index(x, lo, hi);
    ptr = (*q) + ci * 2 * QUDA_Nc * QUDA_Nc;
    for (a = 0; a < QUDA_Nc; a++) {
      for (b = 0; b < QUDA_Nc; b++) {
	ptr[2 * QUDA_Nc * a + 2 * b] = QLA_real(QLA_D3_elem_M(Ux[i], a, b));
	ptr[2 * QUDA_Nc * a + 2 * b + 1] = QLA_imag(QLA_D3_elem_M(Ux[i], a, b));
      }
    }
  }
  QDP_D3_reset_M(*U);
  lua_pop(L, 2);
}

static void
free_gauge_field(QUDA_REAL **q, QDP_D3_ColorMatrix **U, lua_State *L, int idx, int d)
{
  qlua_free(L, *q);
}


static void
get_fermion_field(QUDA_REAL **q, QDP_D3_DiracFermion *f, lua_State *L, mLattice *S)
{
  QLA_D3_DiracFermion *fx;
  int lo[QUDA_DIM];
  int hi[QUDA_DIM];
  int subvol;
  int i;

  qlua_assert(S->rank == QUDA_DIM, "expected rank 4 lattice");
  qlua_sublattice(lo, hi, S->node, S);
  for (i = 0, subvol = 1; i < QUDA_DIM; i++) {
    subvol *= hi[i] - lo[i];
  }
  *q = qlua_malloc(L, subvol * 2 * QUDA_Nc * QUDA_Ns * sizeof (QUDA_REAL));
  fx = QDP_D3_expose_D(f);
  for (i = 0; i < subvol; i++) {
    int c, d, ci, x[QUDA_DIM];
    QUDA_REAL *ptr;
    QDP_get_coords_L(S->lat, x, S->node, i);
    ci = quda_index(x, lo, hi);
    ptr = (*q) + ci * 2 * QUDA_Ns * QUDA_Nc;
    for (c = 0; c < QUDA_Nc; c++) {
      for (d = 0; d < QUDA_Ns; d++) {
	ptr[2 * QUDA_Ns * c + 2 * d] = QLA_real(QLA_D3_elem_D(fx[i], c, d));
	ptr[2 * QUDA_Ns * c + 2 * d + 1] = QLA_imag(QLA_D3_elem_D(fx[i], c, d));
      }
    }
  }
  QDP_D3_reset_D(f);
}

static void
put_fermion_field(QDP_D3_DiracFermion *f, QUDA_REAL *q, lua_State *L, mLattice *S)
{
  QLA_D3_DiracFermion *fx;
  int lo[QUDA_DIM];
  int hi[QUDA_DIM];
  int subvol;
  int i;

  qlua_assert(S->rank == QUDA_DIM, "expected rank 4 lattice");
  qlua_sublattice(lo, hi, S->node, S);
  for (i = 0, subvol = 1; i < QUDA_DIM; i++) {
    subvol *= hi[i] - lo[i];
  }
  fx = QDP_D3_expose_D(f);
  for (i = 0; i < subvol; i++) {
    int c, d, ci, x[QUDA_DIM];
    QUDA_REAL *ptr;
    QDP_get_coords_L(S->lat, x, S->node, i);
    ci = quda_index(x, lo, hi);
    ptr = q + ci * 2 * QUDA_Ns * QUDA_Nc;
    for (c = 0; c < QUDA_Nc; c++) {
      for (d = 0; d < QUDA_Ns; d++) {
	double v_re = ptr[2 * QUDA_Ns * c + 2 * d];
	double v_im = ptr[2 * QUDA_Ns * c + 2 * d + 1];
	QLA_c_eq_r_plus_ir(QLA_D3_elem_D(fx[i], c, d), v_re, v_im);
      }
    }
  }
  QDP_D3_reset_D(f);
}

static void
free_fermion_field(QUDA_REAL **q, QDP_D3_DiracFermion *f, lua_State *L, mLattice *S)
{
  qlua_free(L, *q);
}

static int
qq_loadGaugeQuda(lua_State *L)
{
  int i;
  QudaGaugeParam *p = qq_checkGaugeParam(L, 2);
  QDP_D3_ColorMatrix *U[QUDA_DIM];
  QUDA_REAL *qu[QUDA_DIM];

  luaL_checktype(L, 1, LUA_TTABLE);
  CALL_QDP(L);
  for (i = 0; i < QUDA_DIM; i++)
    get_gauge_field(&qu[i], &U[i], L, 1, i);
  loadGaugeQuda(qu, p);
  for (i = 0; i < QUDA_DIM; i++)
    free_gauge_field(&qu[i], &U[i], L, 1, i);
  
  return 0;
}

static int
qq_loadCloverQuda(lua_State *L)
{
  QudaInvertParam *p = qq_checkInvertParam(L, 1);
  CALL_QDP(L);
  loadCloverQuda(NULL, NULL, p);
  return 0;
}

static int
qq_invertQuda(lua_State *L)
{
  QDP_D3_DiracFermion *rhs = qlua_checkLatDirFerm3(L, 1, NULL, 3)->ptr;
  QudaInvertParam *p = qq_checkInvertParam(L, 2);
  mLattice *S = qlua_ObjLattice(L, 1);
  QDP_D3_DiracFermion *sol = qlua_newZeroLatDirFerm3(L, lua_gettop(L), 3)->ptr;
  QUDA_REAL *q_rhs;
  QUDA_REAL *q_sol;
  
  CALL_QDP(L);
  get_fermion_field(&q_rhs, rhs, L, S);
  get_fermion_field(&q_sol, sol, L, S);
  
  invertQuda(q_sol, q_rhs, p);
  put_fermion_field(sol, q_sol, L, S);
  free_fermion_field(&q_rhs, rhs, L, S);
  free_fermion_field(&q_sol, sol, L, S);
  return 1;
}

static int
qq_performAPEnStep(lua_State *L)
{
  int nSteps = luaL_checkint(L, 1);
  double alpha = luaL_checknumber(L, 2);
  performAPEnStep(nSteps, alpha);

  return 0;
}


static struct luaL_Reg fquda[] = {
  /* QUDA structures */
  {"GaugeParam",              qq_gauge_param           },
  {"InvertParam",             qq_invert_param          },
  /* QUDA functions */
  {"initCommsGridQuda",       qq_initCommsGridQuda     },
  {"initQuda",                qq_initQuda              },
  {"initQudaDevice",          qq_initQudaDevice        },
  {"initQudaMemory",          qq_initQudaMemory        },
  {"endQuda",                 qq_endQuda               },
  {"freeCloverQuda",          qq_freeCloverQuda        },
  {"freeGaugeQuda",           qq_freeGaugeQuda         },
  {"invertQuda",              qq_invertQuda            },
  {"loadCloverQuda",          qq_loadCloverQuda        },
  {"loadGaugeQuda",           qq_loadGaugeQuda         },
  {"performAPEnStep",         qq_performAPEnStep       },
  {"plaqQuda",                qq_plaqQuda              },
  {"qChargeCuda",             qq_qChargeCuda           },
  {"setVerbosityQuda",        qq_setVerbosityQuda      },
  {"openMagma",               qq_openMagma             },
  {"closeMagma",              qq_closeMagma            },
  { NULL,                     NULL                     }
};

int
init_quda(lua_State *L)
{
  luaL_register(L, qudalib, fquda);
  qlua_metatable(L, mtnGaugeParam, mtGaugeParam, qQudaGaugeParam);
  qlua_metatable(L, mtnInvertParam, mtInvertParam, qQudaInvertParam);
  return 0;
}

void
fini_quda(void)
{
  /* free resources of GPU ?? */
}
