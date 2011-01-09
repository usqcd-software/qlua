static const char QopqdpAsqtadName[] = "lattice.QopqdpAsqtad";

typedef struct {
  QOP_D3_FermionLinksAsqtad *fla;
  QOP_F3_FermionLinksAsqtad *ffla;
  double preres;
  double flops, seconds;
  int max_iter;
  int restart_iter;
  int max_restarts;
  int cgtype;
  int float_inner;
  int iterations, restarts;
} mAsqtad;

static mAsqtad *qlua_newAsqtad(lua_State *L, int Sidx);

static mAsqtad *
qlua_checkAsqtad(lua_State *L, int idx, mLattice *S, int live)
{
  mAsqtad *c = qlua_checkLatticeType(L, idx, qQopqdpAsqtad, QopqdpAsqtadName);

  if (S) {
    mLattice *S1 = qlua_ObjLattice(L, idx);
    if (S1->id != S->id)
      luaL_error(L, "%s on a wrong lattice", QopqdpAsqtadName);
    lua_pop(L, 1);
  }

  //if (live && (c->state == 0 || c->gauge == 0))
  //luaL_error(L, "using closed qcd.Asqtad");

  return c;
}

// D(asqtad, in, mass)
static int
q_asqtad_D(lua_State *L)
{
  mAsqtad *c = qlua_checkAsqtad(L, 1, NULL, 1);
  double mass = luaL_checknumber(L, 3);
  mLattice *S = qlua_ObjLattice(L, 1);
  int Sidx = lua_gettop(L);

  switch (qlua_qtype(L, 2)) {
  case qLatColVec3: {
    mLatColVec3 *in = qlua_checkLatColVec3(L, 2, S, 3);
    mLatColVec3 *out = qlua_newLatColVec3(L, Sidx, 3);

    CALL_QDP(L);
    QOP_D3_asqtad_dslash_qdp(NULL, c->fla, mass, out->ptr, in->ptr,
			     QOP_EVENODD, QOP_EVENODD);

    return 1;
  }
  default:
    break;
  }
  return luaL_error(L, "bad arguments in QOPQDP Asqtad Dslash");
}

// Dx(asqtad, in, mass)
static int
q_asqtad_Dx(lua_State *L)
{
  mAsqtad *c = qlua_checkAsqtad(L, 1, NULL, 1);
  double mass = luaL_checknumber(L, 3);
  mLattice *S = qlua_ObjLattice(L, 1);
  int Sidx = lua_gettop(L);

  switch (qlua_qtype(L, 2)) {
  case qLatColVec3: {
    mLatColVec3 *in = qlua_checkLatColVec3(L, 2, S, 3);
    mLatColVec3 *out = qlua_newLatColVec3(L, Sidx, 3);

    CALL_QDP(L);
    QOP_D3_asqtad_dslash_qdp(NULL, c->fla, -mass, out->ptr, in->ptr,
			     QOP_EVENODD, QOP_EVENODD);
    QDP_D3_V_eqm_V(out->ptr, out->ptr, S->all);

    return 1;
  }
  default:
    break;
  }
  return luaL_error(L, "bad arguments in QOPQDP Asqtad Dslash");
}

static void
asqtadInvert(QDP_D3_ColorVector *prop, mAsqtad *c,
	     QDP_D3_ColorVector *source, double mass, double res)
{
  QLA_D_Real rsq, rsqstop, norm2in, rsqmin;
  QDP_D3_ColorVector *Dprop, *r;
  double flops=0, secs=0;
  int iters=0, restarts=0;
  QDP_Lattice *lat = QDP_D3_get_lattice_V(source);
  QDP_Subset all = QDP_all_L(lat);

  QOP_invert_arg_t inv_arg = QOP_INVERT_ARG_DEFAULT;
  QOP_resid_arg_t res_arg = QOP_RESID_ARG_DEFAULT;
  QOP_info_t info;

  inv_arg.max_iter = c->max_iter;
  inv_arg.restart = c->restart_iter;
  //inv_arg.max_restarts = c->max_restarts;
  inv_arg.max_restarts = 0;
  inv_arg.evenodd = QOP_EVENODD;

  QOP_opt_t opt;
  opt.tag = "cg";
  opt.value = c->cgtype;
  QOP_asqtad_invert_set_opts(&opt, 1);

  QDP_D3_r_eq_norm2_V(&norm2in, source, all);
  rsqstop = res*res*norm2in;

  r = QDP_D3_create_V_L(lat);
  Dprop = QDP_D3_create_V_L(lat);
  while(1) {
    QOP_D3_asqtad_dslash_qdp(NULL, c->fla, mass, Dprop, prop,
			     QOP_EVENODD, QOP_EVENODD);
    QDP_D3_V_eq_V_minus_V(r, source, Dprop, all);
    QDP_D3_r_eq_norm2_V(&rsq, r, all);
    if(iters>0) {
      //printf0("iters = %4i  secs = %8f  mflops = %.0f rsq = %g\n", iters,
      //secs, 1e-6*flops/secs, rsq/norm2in);
      printf0("iters = %4i  secs = %8f  mflops = %.0f rsq = %g\n",
	      res_arg.final_iter, info.final_sec,
	      1e-6*info.final_flop/info.final_sec, rsq/norm2in);
    }
    //printf("%g\n", rsq/norm2in);
    if(rsq<rsqstop) break;
    if(restarts>c->max_restarts) break;
    if(iters>0) restarts++;

    rsqmin = c->preres*c->preres;
    if(rsqstop/rsq>rsqmin) rsqmin = rsqstop/rsq;
    rsqmin *= 0.99;
    res_arg.rsqmin = rsqmin;

    if (c->float_inner) {
      QDP_F3_ColorVector *fout, *fin;
      QOP_F3_ColorVector *qopout, *qopin;
      fout = QDP_F3_create_V_L(lat);
      fin = QDP_F3_create_V_L(lat);
      QDP_FD3_V_eq_V(fin, r, all);
      QDP_F3_V_eq_zero(fout, all);
      qopin = QOP_F3_create_V_from_qdp(fin);
      qopout = QOP_F3_create_V_from_qdp(fout);
      QOP_F3_asqtad_invert(&info, c->ffla, &inv_arg, &res_arg, mass,
			   qopout, qopin);
      QOP_F3_extract_V_to_qdp(fout, qopout);
      QDP_DF3_V_eq_V(Dprop, fout, all);
      QDP_F3_destroy_V(fout);
      QDP_F3_destroy_V(fin);
      QOP_F3_destroy_V(qopout);
      QOP_F3_destroy_V(qopin);
    } else {
      QOP_D3_ColorVector *qopout, *qopin;
      QDP_D3_V_eq_zero(Dprop, all);
      qopin = QOP_D3_create_V_from_qdp(r);
      qopout = QOP_D3_create_V_from_qdp(Dprop);
      QOP_D3_asqtad_invert(&info, c->fla, &inv_arg, &res_arg, mass,
			   qopout, qopin);
      QOP_D3_extract_V_to_qdp(Dprop, qopout);
      QOP_D3_destroy_V(qopout);
      QOP_D3_destroy_V(qopin);
    }
    iters += res_arg.final_iter;
    secs += info.final_sec;
    flops += info.final_flop;

    QDP_D3_V_peq_V(prop, Dprop, all);
  }
  printf0("titers = %4i  secs = %8f  mflops = %.0f rsq = %g\n", iters,
	  secs, 1e-6*flops/secs, rsq/norm2in);
  QDP_D3_destroy_V(r);
  QDP_D3_destroy_V(Dprop);
  c->flops = flops;
  c->iterations = iters;
  c->seconds = secs;
  c->restarts = restarts;
}

// solve(asqtad, in, mass, res)
static int
q_asqtad_solve(lua_State *L)
{
  mAsqtad *c = qlua_checkAsqtad(L, 1, NULL, 1);
  double mass = luaL_checknumber(L, 3);
  double res = luaL_checknumber(L, 4);
  mLattice *S = qlua_ObjLattice(L, 1);
  int Sidx = lua_gettop(L);

  switch (qlua_qtype(L, 2)) {
  case qLatColVec3: {
    mLatColVec3 *in = qlua_checkLatColVec3(L, 2, S, 3);
    mLatColVec3 *out = qlua_newLatColVec3(L, Sidx, 3);

    CALL_QDP(L);
    QDP_D3_V_eq_zero(out->ptr, S->all);
    asqtadInvert(out->ptr, c, in->ptr, mass, res);

    return 1;
  }
  default:
    break;
  }
  return luaL_error(L, "bad arguments in QOPQDP Asqtad Dslash");
}

static int
q_asqtad_defaults(lua_State *L)
{
  static const char *cgtype_key = "cgtype";
  static const char *float_key = "float_inner";
  static const char *preres_key = "preres";
  static const char *max_iter_key = "max_iter";
  static const char *restart_iter_key = "restart_iter";
  static const char *max_restarts_key = "max_restarts";
  mAsqtad *c = qlua_checkAsqtad(L, 1, NULL, 1);

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
  return luaL_error(L, "bad parameters for QopqdpAsqtad:defaults()");
}

static int
q_asqtad_stats(lua_State *L)
{
  static const char *flops_key = "flops";
  static const char *seconds_key = "seconds";
  static const char *iterations_key = "iterations";
  static const char *restarts_key = "restarts";
  mAsqtad *c = qlua_checkAsqtad(L, 1, NULL, 1);

  lua_createtable(L, 0, 1);
  add_default(L, flops_key, c->flops);
  add_default(L, seconds_key, c->seconds);
  add_default(L, iterations_key, c->iterations);
  add_default(L, restarts_key, c->restarts);

  return 1;
}

static void
set_asqtad(mAsqtad *asq, QDP_D3_ColorMatrix *gauge[],
	   QOP_asqtad_coeffs_t *coeffs, QOP_bc_t *bc, int ndim)
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
  opt[3].value = asq->cgtype;
  QOP_asqtad_invert_set_opts(opt, 4);
  QDP_set_block_size(64);

  int signmask[4] = {0,1,3,7};
  QOP_staggered_sign_t ss;
  ss.signmask = signmask;
  printf0("loading QOP links... ");
  QOP_GaugeField *qg = QOP_create_G_from_qdp(gauge);
  QOP_D3_rephase_G(qg, bc, &ss);

  QOP_info_t info;
  asq->fla = QOP_asqtad_create_L_from_G(&info, coeffs, qg);
  QOP_destroy_G(qg);
  if(!asq->fla) {
    fprintf(stderr, "error can't allocate asqtad links\n");
    QDP_abort(1);
  }
  printf0("done\n");

  if (asq->float_inner) {
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
    QOP_F3_rephase_G(qgf, bc, &ss);
    asq->ffla = QOP_F3_asqtad_create_L_from_G(&info, coeffs, qgf);
    QOP_F3_destroy_G(qgf);
    for(i=0; i<ndim; i++) QDP_F3_destroy_M(fgauge[i]);
    free(fgauge);
    printf0("done\n");
  }
}

/*
 * qcd.qopqdp.Asqtad(
 *  U,            -- 1, {U0,U1,U2,U3}, a table of color matrices
 *  u0,           -- 2, double, tadpole factor, and/or
 *  coeffs,       -- 2 or 3, table of asqtad coefficients
 *  boundary)     -- last, {r/c, ...}, a table of boundary phases
 * coeffs = { one_link, three_staple, five_staple, seven_staple, lepage, naik }
 */
static int
q_asqtad(lua_State *L)
{
  int nargs = lua_gettop(L);
  if(nargs<1 || nargs>4)
    return luaL_error(L, "invalid number of arguments %i", nargs);

  luaL_checktype(L, 1, LUA_TTABLE);
  lua_pushnumber(L, 1);
  lua_gettable(L, 1);
  qlua_checkLatColMat3(L, -1, NULL, 3);
  mLattice *S = qlua_ObjLattice(L, -1);
  int ndim = S->rank;

  if (ndim != 4)
    return luaL_error(L, "asqtad is not implemented for #L=%d", ndim);
  if (QDP_Ns != 4)
    return luaL_error(L, "asqtad does not support Ns=%d", QDP_Ns);

  QDP_D3_ColorMatrix *gauge[ndim];
  double u0=1;
  QOP_Complex bcphase[ndim];
  QOP_bc_t bc;
  bc.phase = bcphase;
  int coeffindex=2, bcindex=2;

  // get gauge field
  for (int i = 0; i < ndim; i++) {
    lua_pushnumber(L, i + 1); /* [sic] lua indexing */
    lua_gettable(L, 1);
    gauge[i] = qlua_checkLatColMat3(L, -1, S, 3)->ptr;
    lua_pop(L, 1);
  }

  // get u0
  if(nargs>=2) {
    if(lua_type(L, 2)==LUA_TNUMBER) {
      u0 = lua_tonumber(L, 2);
      coeffindex = 3;
      bcindex = 3;
    }
  }

  // create default coefficients
  QOP_asqtad_coeffs_t coeffs;
  {
    double u2, u4;
    u2 = 1.0/(u0*u0);
    u4 = u2*u2;
    coeffs.one_link = 5.0/8.0;
    coeffs.three_staple = -u2/16.0;
    coeffs.five_staple = u4/64.0;
    coeffs.seven_staple = -(u4*u2)/384.0;
    coeffs.lepage = -u4/16.0;
    coeffs.naik = -u2/24.0;
  }

  if(nargs>=coeffindex) {
    if(lua_type(L, coeffindex)==LUA_TTABLE) {
      lua_pushnumber(L, 1);
      lua_gettable(L, coeffindex);
      if(lua_type(L, -1)==LUA_TNIL) { // not bc table
#define getcoeff(s) lua_pop(L, 1); lua_pushstring(L, #s); \
	lua_gettable(L, coeffindex);					\
	if(lua_type(L, -1)==LUA_TNUMBER) coeffs.s = lua_tonumber(L, -1);
	getcoeff(one_link);
	getcoeff(three_staple);
	getcoeff(five_staple);
	getcoeff(seven_staple);
	getcoeff(lepage);
	getcoeff(naik);
	bcindex = coeffindex + 1;
      }
      lua_pop(L, 1);
    }
  }

  printf0("  u0 = %g\n", u0);
#define printcoeff(s) printf("  " #s " = %g\n", coeffs.s);
  printcoeff(one_link);
  printcoeff(three_staple);
  printcoeff(five_staple);
  printcoeff(seven_staple);
  printcoeff(lepage);
  printcoeff(naik);

  // get BC
  for(int i=0; i<ndim; i++) { bcphase[i].re = 1; bcphase[i].im = 0; }
  bcphase[ndim-1].re = -1; // default antiperiodic
  if(nargs>=bcindex) {
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
	luaL_error(L, "bad asqtad boundary condition type");
      }
      lua_pop(L, 1);
    }
  }

  int Sidx = lua_gettop(L);
  mAsqtad *c = qlua_newAsqtad(L, Sidx);
  printf0("created asqtad\n");
  set_asqtad(c, gauge, &coeffs, &bc, ndim);
  printf0("set asqtad\n");

  return 1;
}

static struct luaL_Reg mtAsqtad[] = {
  //{ "__tostring",   q_CL_fmt },
  //{ "__gc",         q_CL_gc },
  //{ "close",        q_CL_close },
  { "D",            q_asqtad_D },
  { "Dx",           q_asqtad_Dx },
  { "solve",        q_asqtad_solve },
  { "defaults",     q_asqtad_defaults },
  { "stats",        q_asqtad_stats },
  //{ "mixed_solver", q_CL_make_mixed_solver },
  //{ "eig_deflator", q_CL_make_deflator },
    { NULL, NULL }
};

static mAsqtad *
qlua_newAsqtad(lua_State *L, int Sidx)
{
  mAsqtad *c = lua_newuserdata(L, sizeof (mAsqtad));

  c->fla = NULL;
  c->ffla = NULL;
  c->preres = 5e-6;
  c->max_iter = 2000;
  c->restart_iter = 2000;
  c->max_restarts = 5;
  c->cgtype = 0;
  c->float_inner = 1;
  c->flops = 0;
  c->iterations = 0;
  c->seconds = 0;
  c->restarts = 0;

  qlua_createLatticeTable(L, Sidx, mtAsqtad, qQopqdpAsqtad, QopqdpAsqtadName);
  lua_setmetatable(L, -2);

  return c;
}
