// Colored vectors
static void
Qs(wpack_colvec)(lua_State *L, void *data, SHA256_Context *ctx, void *src, int nc, WriteSize wsize)
{
#if QNc == 'N'
  typedef QLA_DN_ColorVector(nc, DataType);
#else
  typedef Qx(QLA_D,_ColorVector) DataType;
#endif
  DataType *ptr = (DataType *)src;
  int c;

  switch (wsize) {
  case WS_Double: {
    machine_complex_double *dst = (machine_complex_double *)data;
    for (c = 0; c < nc; c++, dst++) {
      QLA_D_Complex zz;
      Qx(QLA_D,_C_eq_elem_V)(QNC(nc) &zz, ptr, c);
      dst->re = QLA_real(zz);
      dst->im = QLA_imag(zz);
    }
    sha256_sum_add_doubles(ctx, data, nc * 2);
  } break;
  case WS_Float: {
    machine_complex_float *dst = (machine_complex_float *)data;
    for (c = 0; c < nc; c++, dst++) {
      QLA_D_Complex zz;
      Qx(QLA_D,_C_eq_elem_V)(QNC(nc) &zz, ptr, c);
      dst->re = QLA_real(zz);
      dst->im = QLA_imag(zz);
    }
    sha256_sum_add_floats(ctx, data, nc * 2);
  } break;
  default:
    luaL_error(L, "unknown wsize in wpack_colvecN()");
    break;
  }
}

static void
Qs(w_colvec)(lua_State *L, mHdf5File *b, mLattice *S,
             struct wopts_s *opts, struct laddr_s *laddr,
             SHA256_Sum *sum, void **data, hid_t *filetype, hid_t *memtype,
             const char **kind)
{
  Qs(mSeqColVec) *m = Qs(qlua_checkSeqColVec)(L, 3, -1);
  int nc = QC(m);
  int dsize;

  switch (opts->wsize) {
  case WS_Double:
    dsize = sizeof (machine_complex_double);
    break;
  case WS_Float:
    dsize = sizeof (machine_complex_float);
    break;
  default:
    luaL_error(L, "Unknown precision in w_colvec()");
  }
  dsize = dsize * nc;
  *data = qlua_malloc(L, dsize);
  SHA256_Context *ctx = sha256_create(L);
  Qs(wpack_colvec)(L, *data, ctx, m->ptr, nc, opts->wsize);
  sha256_sum(sum, ctx);
  sha256_destroy(ctx);
  *kind = "ColorVector";
  *filetype = get_colvec_type(L, b, nc, opts->wsize, 1);
  *memtype  = get_colvec_type(L, b, nc, opts->wsize, 0);
}

static void
Qs(w_latcolvec)(lua_State *L, mHdf5File *b, mLattice *S,
             struct wopts_s *opts, struct laddr_s *laddr,
             SHA256_Sum *sum, void **data, hid_t *filetype, hid_t *memtype,
             const char **kind)
{
  /* XXX */
}

static void
Qs(w_colmat)(lua_State *L, mHdf5File *b, mLattice *S,
             struct wopts_s *opts, struct laddr_s *laddr,
             SHA256_Sum *sum, void **data, hid_t *filetype, hid_t *memtype,
             const char **kind)
{
  /* XXX */
}

static void
Qs(w_latcolmat)(lua_State *L, mHdf5File *b, mLattice *S,
             struct wopts_s *opts, struct laddr_s *laddr,
             SHA256_Sum *sum, void **data, hid_t *filetype, hid_t *memtype,
             const char **kind)
{
  /* XXX */
}

static void
Qs(w_dirferm)(lua_State *L, mHdf5File *b, mLattice *S,
             struct wopts_s *opts, struct laddr_s *laddr,
             SHA256_Sum *sum, void **data, hid_t *filetype, hid_t *memtype,
             const char **kind)
{
  /* XXX */
}

static void
Qs(w_latdirferm)(lua_State *L, mHdf5File *b, mLattice *S,
             struct wopts_s *opts, struct laddr_s *laddr,
             SHA256_Sum *sum, void **data, hid_t *filetype, hid_t *memtype,
             const char **kind)
{
  /* XXX */
}

static void
Qs(w_dirprop)(lua_State *L, mHdf5File *b, mLattice *S,
             struct wopts_s *opts, struct laddr_s *laddr,
             SHA256_Sum *sum, void **data, hid_t *filetype, hid_t *memtype,
             const char **kind)
{
  /* XXX */
}

static void
Qs(w_latdirprop)(lua_State *L, mHdf5File *b, mLattice *S,
             struct wopts_s *opts, struct laddr_s *laddr,
             SHA256_Sum *sum, void **data, hid_t *filetype, hid_t *memtype,
             const char **kind)
{
  /* XXX */
}


#undef QNc
#undef Qcolors
#undef Qs
#undef Qx
#undef QC
#undef QNC
