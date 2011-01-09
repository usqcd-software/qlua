static const char QopqdpCloverName[] = "lattice.QopqdpClover";

typedef struct {
  QOP_D3_FermionLinksWilson *flw;
  QOP_F3_FermionLinksWilson *fflw;
  double preres;
  double flops, seconds;
  int max_iter;
  int restart_iter;
  int max_restarts;
  int cgtype;
  int float_inner;
  int iterations, restarts;
} mClover;

static mClover *qlua_newClover(lua_State *L, int Sidx);

static mClover *
qlua_checkClover(lua_State *L, int idx, mLattice *S, int live)
{
  mClover *c = qlua_checkLatticeType(L, idx, qQopqdpClover, QopqdpCloverName);

  if (S) {
    mLattice *S1 = qlua_ObjLattice(L, idx);
    if (S1->id != S->id)
      luaL_error(L, "%s on a wrong lattice", QopqdpCloverName);
    lua_pop(L, 1);
  }

  //if (live && (c->state == 0 || c->gauge == 0))
  //luaL_error(L, "using closed qcd.Clover");

  return c;
}

// D(clover, in, kappa)
static int
q_clover_D(lua_State *L)
{
  mClover *c = qlua_checkClover(L, 1, NULL, 1);
  double kappa = luaL_checknumber(L, 3);
  mLattice *S = qlua_ObjLattice(L, 1);
  int Sidx = lua_gettop(L);

  switch (qlua_qtype(L, 2)) {
  case qLatDirFerm3: {
    mLatDirFerm3 *in = qlua_checkLatDirFerm3(L, 2, S, 3);
    mLatDirFerm3 *out = qlua_newLatDirFerm3(L, Sidx, 3);

    CALL_QDP(L);
    QOP_D3_wilson_dslash_qdp(NULL, c->flw, kappa, 1, out->ptr, in->ptr,
			     QOP_EVENODD, QOP_EVENODD);

    return 1;
  }
  default:
    break;
  }
  return luaL_error(L, "bad arguments in QOPQDP Clover Dslash");
}

// Dx(clover, in, kappa)
static int
q_clover_Dx(lua_State *L)
{
  mClover *c = qlua_checkClover(L, 1, NULL, 1);
  double kappa = luaL_checknumber(L, 3);
  mLattice *S = qlua_ObjLattice(L, 1);
  int Sidx = lua_gettop(L);

  switch (qlua_qtype(L, 2)) {
  case qLatDirFerm3: {
    mLatDirFerm3 *in = qlua_checkLatDirFerm3(L, 2, S, 3);
    mLatDirFerm3 *out = qlua_newLatDirFerm3(L, Sidx, 3);

    CALL_QDP(L);
    QOP_D3_wilson_dslash_qdp(NULL, c->flw, kappa, -1, out->ptr, in->ptr,
			     QOP_EVENODD, QOP_EVENODD);

    return 1;
  }
  default:
    break;
  }
  return luaL_error(L, "bad arguments in QOPQDP Clover Dslash");
}

static void
cloverInvert(QDP_D3_DiracFermion *prop, mClover *c,
	     QDP_D3_DiracFermion *source, double kappa, double res)
{
  QLA_D_Real rsq, rsqstop, norm2in, rsqmin;
  QDP_D3_DiracFermion *Dprop, *r;
  double flops=0, secs=0;
  int iters=0, restarts=0;
  QDP_Lattice *lat = QDP_D3_get_lattice_D(source);
  QDP_Subset all = QDP_all_L(lat);

  QOP_invert_arg_t inv_arg = QOP_INVERT_ARG_DEFAULT;
  QOP_resid_arg_t res_arg = QOP_RESID_ARG_DEFAULT;
  QOP_info_t info;

  inv_arg.max_iter = c->max_iter;
  inv_arg.restart = c->restart_iter;
  inv_arg.max_restarts = c->max_restarts;
  inv_arg.evenodd = QOP_EVENODD;

  QOP_opt_t opt;
  opt.tag = "cg";
  opt.value = c->cgtype;
  QOP_wilson_invert_set_opts(&opt, 1);

  QDP_D3_r_eq_norm2_D(&norm2in, source, all);
  rsqstop = res*res*norm2in;

  r = QDP_D3_create_D_L(lat);
  Dprop = QDP_D3_create_D_L(lat);
  while(1) {
    QOP_D3_wilson_dslash_qdp(NULL, c->flw, kappa, 1, Dprop, prop,
			     QOP_EVENODD, QOP_EVENODD);
    QDP_D3_D_eq_D_minus_D(r, source, Dprop, all);
    QDP_D3_r_eq_norm2_D(&rsq, r, all);
    if(iters>0)
      printf0("iters = %4i  secs = %8f  mflops = %.0f rsq = %g\n", iters,
              secs, 1e-6*flops/secs, rsq/norm2in);
    //printf("%g\n", rsq/norm2in);
    if(rsq<rsqstop) break;
    if(iters>0) restarts++;

    rsqmin = c->preres*c->preres;
    if(rsqstop/rsq>rsqmin) rsqmin = rsqstop/rsq;
    rsqmin *= 0.99;
    res_arg.rsqmin = rsqmin;

    if (c->float_inner) {
      QDP_F3_DiracFermion *fout, *fin;
      QOP_F3_DiracFermion *qopout, *qopin;
      fout = QDP_F3_create_D_L(lat);
      fin = QDP_F3_create_D_L(lat);
      QDP_FD3_D_eq_D(fin, r, all);
      QDP_F3_D_eq_zero(fout, all);
      qopin = QOP_F3_create_D_from_qdp(fin);
      qopout = QOP_F3_create_D_from_qdp(fout);
      QOP_F3_wilson_invert(&info, c->fflw, &inv_arg, &res_arg, kappa,
			   qopout, qopin);
      QOP_F3_extract_D_to_qdp(fout, qopout);
      QDP_DF3_D_eq_D(Dprop, fout, all);
      QDP_F3_destroy_D(fout);
      QDP_F3_destroy_D(fin);
      QOP_F3_destroy_D(qopout);
      QOP_F3_destroy_D(qopin);
    } else {
      QOP_D3_DiracFermion *qopout, *qopin;
      QDP_D3_D_eq_zero(Dprop, all);
      qopin = QOP_D3_create_D_from_qdp(r);
      qopout = QOP_D3_create_D_from_qdp(Dprop);
      QOP_D3_wilson_invert(&info, c->flw, &inv_arg, &res_arg, kappa,
			   qopout, qopin);
      QOP_D3_extract_D_to_qdp(Dprop, qopout);
      QOP_D3_destroy_D(qopout);
      QOP_D3_destroy_D(qopin);
    }
    iters += res_arg.final_iter;
    secs += info.final_sec;
    flops += info.final_flop;

    QDP_D3_D_peq_D(prop, Dprop, all);
  }
  QDP_D3_destroy_D(r);
  QDP_D3_destroy_D(Dprop);
  c->flops = flops;
  c->iterations = iters;
  c->seconds = secs;
  c->restarts = restarts;
}

// solve(clover, in, kappa, res)
static int
q_clover_solve(lua_State *L)
{
  mClover *c = qlua_checkClover(L, 1, NULL, 1);
  double kappa = luaL_checknumber(L, 3);
  double res = luaL_checknumber(L, 4);
  mLattice *S = qlua_ObjLattice(L, 1);
  int Sidx = lua_gettop(L);

  switch (qlua_qtype(L, 2)) {
  case qLatDirFerm3: {
    mLatDirFerm3 *in = qlua_checkLatDirFerm3(L, 2, S, 3);
    mLatDirFerm3 *out = qlua_newLatDirFerm3(L, Sidx, 3);

    CALL_QDP(L);
    QDP_D3_D_eq_zero(out->ptr, QDP_all_L(S->lat));
    cloverInvert(out->ptr, c, in->ptr, kappa, res);

    return 1;
  }
  default:
    break;
  }
  return luaL_error(L, "bad arguments in QOPQDP Clover Dslash");
}

static int
q_clover_defaults(lua_State *L)
{
  static const char *cgtype_key = "cgtype";
  static const char *float_key = "float_inner";
  static const char *preres_key = "preres";
  static const char *max_iter_key = "max_iter";
  static const char *restart_iter_key = "restart_iter";
  static const char *max_restarts_key = "max_restarts";
  mClover *c = qlua_checkClover(L, 1, NULL, 1);

  switch (lua_gettop(L)) {
  case 1: {
    lua_createtable(L, 0, 1);
    add_default(L, cgtype_key, c->cgtype);
    add_default(L, float_key, c->float_inner);
    add_default(L, preres_key, c->preres);
    add_default(L, max_iter_key, c->max_iter);
    add_default(L, restart_iter_key, c->restart_iter);
    add_default(L, max_restarts_key, c->max_restarts);
    return 1;
  }
  case 2: {
    if (qlua_qtype(L, 2) == qTable) {

      int cg = (int) set_default(L, cgtype_key, c->cgtype);
      if (cg == 0 || cg == 1 || cg == 3)
	c->cgtype = cg;
      else
	return luaL_error(L, "bad default cgtype %d", cg);

      int fi = (int) set_default(L, float_key, c->float_inner);
      if (fi == 0 || fi == 1)
	c->float_inner = fi;
      else
	return luaL_error(L, "bad default float_inner %d", fi);

      double preres = set_default(L, preres_key, c->preres);
      if (preres >= 0)
	c->preres = preres;
      else
	return luaL_error(L, "bad default preres %d", preres);

      int max_iter = (int) set_default(L, max_iter_key, c->max_iter);
      if (max_iter >= 0)
	c->max_iter = max_iter;
      else
	return luaL_error(L, "bad default max_iter %d", max_iter);

      int restart_iter = (int)set_default(L,restart_iter_key,c->restart_iter);
      if (restart_iter >= 0)
	c->restart_iter = restart_iter;
      else
	return luaL_error(L, "bad default restart_iter %d", restart_iter);

      int max_restarts = (int)set_default(L,max_restarts_key,c->max_restarts);
      if (max_restarts >= 0)
	c->max_restarts = max_restarts;
      else
	return luaL_error(L, "bad default max_restarts %d", max_restarts);

      return 1;
    }
    break;
  }
  default:
    break;
  }
  return luaL_error(L, "bad parameters for QopqdpClover:defaults()");
}

static int
q_clover_stats(lua_State *L)
{
  static const char *flops_key = "flops";
  static const char *seconds_key = "seconds";
  static const char *iterations_key = "iterations";
  static const char *restarts_key = "restarts";
  mClover *c = qlua_checkClover(L, 1, NULL, 1);

  lua_createtable(L, 0, 1);
  add_default(L, flops_key, c->flops);
  add_default(L, seconds_key, c->seconds);
  add_default(L, iterations_key, c->iterations);
  add_default(L, restarts_key, c->restarts);

  return 1;
}

static void
set_clover(mClover *wi, QDP_D3_ColorMatrix *gauge[], double aniso,
	   double clov_s, double clov_t, QOP_bc_t *bc, int ndim)
{
  QDP_Lattice *lat = QDP_D3_get_lattice_M(gauge[0]);
  Init(lat);

  QOP_opt_t opt[4];
  opt[0].tag = "st";
  opt[0].value = 3;
  opt[1].tag = "ns";
  opt[1].value = 8;
  opt[2].tag = "nm";
  opt[2].value = 8;
  opt[3].tag = "cg";
  opt[3].value = wi->cgtype;
  QOP_wilson_invert_set_opts(opt, 4);
  QDP_set_block_size(64);

  QOP_info_t info;
  QOP_wilson_coeffs_t coeffs;
  coeffs.aniso = 1/aniso;
  coeffs.clov_s = clov_s;
  coeffs.clov_t = clov_t;
  printf0("loading QOP links... ");
  QOP_GaugeField *qg = QOP_create_G_from_qdp(gauge);
  QOP_D3_rephase_G(qg, bc, NULL);
  wi->flw = QOP_wilson_create_L_from_G(&info, &coeffs, qg);
  QOP_destroy_G(qg);
  if(!wi->flw) {
    fprintf(stderr, "error can't allocate wilson links\n");
    QDP_abort(1);
  }
  printf0("done\n");

  if (wi->float_inner) {
    printf0("loading QOP float links... ");
    QDP_F3_ColorMatrix **fgauge;
    int i, ndim;
    ndim = QDP_ndim();
    fgauge = (QDP_F3_ColorMatrix **) malloc(ndim*sizeof(QDP_F3_ColorMatrix *));
    for(i=0; i<ndim; i++) {
      fgauge[i] = QDP_F3_create_M();
      QDP_FD3_M_eq_M(fgauge[i], gauge[i], QDP_all);
    }
    QOP_F3_GaugeField *qgf = QOP_F3_create_G_from_qdp(fgauge);
    QOP_F3_rephase_G(qgf, bc, NULL);
    wi->fflw = QOP_F3_wilson_create_L_from_G(&info, &coeffs, qgf);
    QOP_F3_destroy_G(qgf);
    for(i=0; i<ndim; i++) QDP_F3_destroy_M(fgauge[i]);
    free(fgauge);
    printf0("done\n");
  }
}

/*
 * qcd.qopqdp.Clover(
 *  U,            -- 1, {U0,U1,U2,U3}, a table of color matrices
 *  c_sw,         -- 2, double, the clover term, or
 *  c_s,c_t,aniso -- 2,3,4, double, the anisotropic clover terms and anisotropy
 *  boundary)     -- last, {r/c, ...}, a table of boundary phases
 */
static int
q_clover(lua_State *L)
{
  int nargs = lua_gettop(L);
  if(nargs<1 || nargs>5)
    return luaL_error(L, "invalid number of arguments %i", nargs);

  luaL_checktype(L, 1, LUA_TTABLE);
  lua_pushnumber(L, 1);
  lua_gettable(L, 1);
  qlua_checkLatColMat3(L, -1, NULL, 3);
  mLattice *S = qlua_ObjLattice(L, -1);
  int ndim = S->rank;

  if (ndim != 4)
    return luaL_error(L, "clover is not implemented for #L=%d", ndim);
  if (QDP_Ns != 4)
    return luaL_error(L, "clover does not support Ns=%d", QDP_Ns);

  QDP_D3_ColorMatrix *gauge[ndim];
  double aniso=1, clov_s=0, clov_t=0;
  QOP_Complex bcphase[ndim];
  QOP_bc_t bc;
  bc.phase = bcphase;
  int getbc=0;

  // get gauge field
  for (int i = 0; i < ndim; i++) {
    lua_pushnumber(L, i + 1); /* [sic] lua indexing */
    lua_gettable(L, 1);
    gauge[i] = qlua_checkLatColMat3(L, -1, S, 3)->ptr;
    lua_pop(L, 1);
  }

  // get clover & anisotropy
  if(nargs>=2) {
    if(lua_type(L, 2)==LUA_TNUMBER) {
      clov_s = lua_tonumber(L, 2);
      clov_t = clov_s; // assume isotropic for now
    } else {
      if(nargs>2)
	return luaL_error(L, "invalid extra arguments %d>2", nargs);
      getbc = 1;
    }
  }
  if(nargs>=3) {
    if(lua_type(L, 3)==LUA_TNUMBER) {
      clov_t = lua_tonumber(L, 3);
      if(nargs<4 || lua_type(L, 4)!=LUA_TNUMBER)
	return luaL_error(L, "missing argument 4>%d", nargs);
      aniso = lua_tonumber(L, 4);
    } else {
      if(nargs>3)
	return luaL_error(L, "invalid extra arguments %d>3", nargs);
      getbc = 1;
    }
  }
  if(nargs>=5) {
    getbc = 1;
  }

  // get BC
  for(int i=0; i<ndim; i++) { bcphase[i].re = 1; bcphase[i].im = 0; }
  bcphase[ndim-1].re = -1; // default antiperiodic
  if(getbc) {
    luaL_checktype(L, nargs, LUA_TTABLE);
    for(int i=0; i<ndim; i++) {
      lua_pushnumber(L, i + 1);
      lua_gettable(L, nargs);
      switch (qlua_qtype(L, -1)) {
      case qReal:
	bcphase[i].re = lua_tonumber(L, -1);
	bcphase[i].im = 0;
	break;
      case qComplex:
	bcphase[i].re = QLA_real(*qlua_checkComplex(L, -1));
	bcphase[i].im = QLA_imag(*qlua_checkComplex(L, -1));
	break;
      default:
	luaL_error(L, "bad clover boundary condition type");
      }
      lua_pop(L, 1);
    }
  }

  int Sidx = lua_gettop(L);
  mClover *c = qlua_newClover(L, Sidx);
  printf0("created clover\n");
  set_clover(c, gauge, aniso, clov_s, clov_t, &bc, ndim);
  printf0("set clover\n");

  return 1;
}

static struct luaL_Reg mtClover[] = {
  //{ "__tostring",   q_CL_fmt },
  //{ "__gc",         q_CL_gc },
  //{ "close",        q_CL_close },
  { "D",            q_clover_D },
  { "Dx",           q_clover_Dx },
  { "solve",        q_clover_solve },
  { "defaults",     q_clover_defaults },
  { "stats",        q_clover_stats },
  //{ "mixed_solver", q_CL_make_mixed_solver },
  //{ "eig_deflator", q_CL_make_deflator },
    { NULL, NULL }
};

static mClover *
qlua_newClover(lua_State *L, int Sidx)
{
  mClover *c = lua_newuserdata(L, sizeof (mClover));

  c->flw = NULL;
  c->fflw = NULL;
  c->preres = 5e-6;
  c->max_iter = 2000;
  c->restart_iter = 2000;
  c->max_restarts = 5;
  c->cgtype = 1;
  c->float_inner = 1;
  c->flops = 0;
  c->iterations = 0;
  c->seconds = 0;
  c->restarts = 0;

  qlua_createLatticeTable(L, Sidx, mtClover, qQopqdpClover, QopqdpCloverName);
  lua_setmetatable(L, -2);

  return c;
}
