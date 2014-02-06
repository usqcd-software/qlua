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
  *data = qlua_malloc(L, dsize * nc);
  SHA256_Context *ctx = sha256_create(L);
  Qs(wpack_colvec)(L, *data, ctx, m->ptr, nc, opts->wsize);
  sha256_sum(sum, ctx);
  sha256_destroy(ctx);
  *kind = knColorVector;
  *filetype = get_colvec_type(L, b, nc, opts->wsize, 1);
  *memtype  = get_colvec_type(L, b, nc, opts->wsize, 0);
}

static void
Qs(w_latcolvec)(lua_State *L, mHdf5File *b, mLattice *S,
             struct wopts_s *opts, struct laddr_s *laddr,
             SHA256_Sum *sum, void **data, hid_t *filetype, hid_t *memtype,
             const char **kind)
{
  Qs(mLatColVec) *m = Qs(qlua_checkLatColVec)(L, 3, S, -1);
  int nc = QC(m);
  int *local_x = qlua_malloc(L, laddr->rank * sizeof (int));
  SHA256_Context *ctx = sha256_create(L);
  int volume = laddr->volume;
  int rank = laddr->rank;
  int dsize, i;
#if QNc == 'N'
  typedef QLA_DN_ColorVector(nc, Vtype);
#else
  typedef Qx(QLA_D,_ColorVector) Vtype;
#endif

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
  *data = qlua_malloc(L, dsize * nc * volume);
  CALL_QDP(L);
  Vtype *locked = Qx(QDP_D,_expose_V)(m->ptr);
  char *out_mem = *data;
  for (i = 0; i < volume; i++) {
    qdp2hdf5_addr(local_x, i, laddr);
    QLUA_ASSERT(QDP_node_number_L(S->lat, local_x) == QDP_this_node);
    Vtype *ptr = &locked[QDP_index_L(S->lat, local_x)];
    sha256_reset(ctx);
    sha256_sum_add_ints(ctx, &rank, 1);
    sha256_sum_add_ints(ctx, local_x, rank);
    Qs(wpack_colvec)(L, out_mem, ctx, ptr, nc, opts->wsize);
    out_mem += dsize * nc;
    SHA256_Sum l_sum;
    sha256_sum(&l_sum, ctx);
    local_combine_checksums(sum, &l_sum);
  }
  Qx(QDP_D, _reset_V)(m->ptr);
  sha256_destroy(ctx);
  qlua_free(L, local_x);
  *kind = knLatticeColorVector;
  *filetype = get_colvec_type(L, b, nc, opts->wsize, 1);
  *memtype  = get_colvec_type(L, b, nc, opts->wsize, 0);
}

// Color Matrix
static void
Qs(wpack_colmat)(lua_State *L, void *data, SHA256_Context *ctx, void *src, int nc, WriteSize wsize)
{
#if QNc == 'N'
  typedef QLA_DN_ColorMatrix(nc, DataType);
#else
  typedef Qx(QLA_D,_ColorMatrix) DataType;
#endif
  DataType *ptr = (DataType *)src;
  int ci, cj;

  switch (wsize) {
  case WS_Double: {
    machine_complex_double *dst = (machine_complex_double *)data;
    for (ci = 0; ci < nc; ci++) {
      for (cj = 0; cj < nc; cj++, dst++) {
        QLA_D_Complex zz;
        Qx(QLA_D,_C_eq_elem_M)(QNC(nc) &zz, ptr, ci,cj);
        dst->re = QLA_real(zz);
        dst->im = QLA_imag(zz);
      }
    }
    sha256_sum_add_doubles(ctx, data, nc * nc * 2);
  } break;
  case WS_Float: {
    machine_complex_float *dst = (machine_complex_float *)data;
    for (ci = 0; ci < nc; ci++) {
      for (cj = 0; cj < nc; cj++, dst++) {
        QLA_D_Complex zz;
        Qx(QLA_D,_C_eq_elem_M)(QNC(nc) &zz, ptr, ci, cj);
        dst->re = QLA_real(zz);
        dst->im = QLA_imag(zz);
      }
    }
    sha256_sum_add_floats(ctx, data, nc * nc * 2);
  } break;
  default:
    luaL_error(L, "unknown wsize in wpack_colvecN()");
    break;
  }
}

static void
Qs(w_colmat)(lua_State *L, mHdf5File *b, mLattice *S,
             struct wopts_s *opts, struct laddr_s *laddr,
             SHA256_Sum *sum, void **data, hid_t *filetype, hid_t *memtype,
             const char **kind)
{
  Qs(mSeqColMat) *m = Qs(qlua_checkSeqColMat)(L, 3, -1);
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
    luaL_error(L, "Unknown precision in w_colmat()");
  }
  *data = qlua_malloc(L, dsize * nc * nc);
  SHA256_Context *ctx = sha256_create(L);
  Qs(wpack_colmat)(L, *data, ctx, m->ptr, nc, opts->wsize);
  sha256_sum(sum, ctx);
  sha256_destroy(ctx);
  *kind = knColorMatrix;
  *filetype = get_colmat_type(L, b, nc, opts->wsize, 1);
  *memtype  = get_colmat_type(L, b, nc, opts->wsize, 0);
}

static void
Qs(w_latcolmat)(lua_State *L, mHdf5File *b, mLattice *S,
             struct wopts_s *opts, struct laddr_s *laddr,
             SHA256_Sum *sum, void **data, hid_t *filetype, hid_t *memtype,
             const char **kind)
{
  Qs(mLatColMat) *m = Qs(qlua_checkLatColMat)(L, 3, S, -1);
  int nc = QC(m);
  int *local_x = qlua_malloc(L, laddr->rank * sizeof (int));
  SHA256_Context *ctx = sha256_create(L);
  int volume = laddr->volume;
  int rank = laddr->rank;
  int dsize, i;
#if QNc == 'N'
  typedef QLA_DN_ColorMatrix(nc, Vtype);
#else
  typedef Qx(QLA_D,_ColorMatrix) Vtype;
#endif

  switch (opts->wsize) {
  case WS_Double:
    dsize = sizeof (machine_complex_double);
    break;
  case WS_Float:
    dsize = sizeof (machine_complex_float);
    break;
  default:
    luaL_error(L, "Unknown precision in w_colmat()");
  }
  *data = qlua_malloc(L, dsize * nc * nc * volume);
  CALL_QDP(L);
  Vtype *locked = Qx(QDP_D,_expose_M)(m->ptr);
  char *out_mem = *data;
  for (i = 0; i < volume; i++) {
    qdp2hdf5_addr(local_x, i, laddr);
    QLUA_ASSERT(QDP_node_number_L(S->lat, local_x) == QDP_this_node);
    Vtype *ptr = &locked[QDP_index_L(S->lat, local_x)];
    sha256_reset(ctx);
    sha256_sum_add_ints(ctx, &rank, 1);
    sha256_sum_add_ints(ctx, local_x, rank);
    Qs(wpack_colmat)(L, out_mem, ctx, ptr, nc, opts->wsize);
    out_mem += dsize * nc * nc;
    SHA256_Sum l_sum;
    sha256_sum(&l_sum, ctx);
    local_combine_checksums(sum, &l_sum);
  }
  Qx(QDP_D, _reset_M)(m->ptr);
  sha256_destroy(ctx);
  qlua_free(L, local_x);
  *kind = knLatticeColorMatrix;
  *filetype = get_colmat_type(L, b, nc, opts->wsize, 1);
  *memtype  = get_colmat_type(L, b, nc, opts->wsize, 0);
}

// Dirac Fermions
static void
Qs(wpack_dirferm)(lua_State *L, void *data, SHA256_Context *ctx, void *src, int nc, WriteSize wsize)
{
#if QNc == 'N'
  typedef QLA_DN_DiracFermion(nc, DataType);
#else
  typedef Qx(QLA_D,_DiracFermion) DataType;
#endif
  DataType *ptr = (DataType *)src;
  int c, d;

  switch (wsize) {
  case WS_Double: {
    machine_complex_double *dst = (machine_complex_double *)data;
    for (d = 0; d < QDP_Ns; d++) {
      for (c = 0; c < nc; c++, dst++) {
        QLA_D_Complex zz;
        Qx(QLA_D,_C_eq_elem_D)(QNC(nc) &zz, ptr, c, d);
        dst->re = QLA_real(zz);
        dst->im = QLA_imag(zz);
      }
    }
    sha256_sum_add_doubles(ctx, data, nc * QDP_Ns * 2);
  } break;
  case WS_Float: {
    machine_complex_float *dst = (machine_complex_float *)data;
    for (d = 0; d < QDP_Ns; d++) {
      for (c = 0; c < nc; c++, dst++) {
        QLA_D_Complex zz;
        Qx(QLA_D,_C_eq_elem_D)(QNC(nc) &zz, ptr, c, d);
        dst->re = QLA_real(zz);
        dst->im = QLA_imag(zz);
      }
    }
    sha256_sum_add_floats(ctx, data, nc * QDP_Ns * 2);
  } break;
  default:
    luaL_error(L, "unknown wsize in wpack_dirfermN()");
    break;
  }
}

static void
Qs(w_dirferm)(lua_State *L, mHdf5File *b, mLattice *S,
             struct wopts_s *opts, struct laddr_s *laddr,
             SHA256_Sum *sum, void **data, hid_t *filetype, hid_t *memtype,
             const char **kind)
{
  Qs(mSeqDirFerm) *m = Qs(qlua_checkSeqDirFerm)(L, 3, -1);
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
    luaL_error(L, "Unknown precision in w_dirferm()");
  }
  *data = qlua_malloc(L, dsize * nc * QDP_Ns);
  SHA256_Context *ctx = sha256_create(L);
  Qs(wpack_dirferm)(L, *data, ctx, m->ptr, nc, opts->wsize);
  sha256_sum(sum, ctx);
  sha256_destroy(ctx);
  *kind = knDiracFermion;
  *filetype = get_dirferm_type(L, b, nc, opts->wsize, 1);
  *memtype  = get_dirferm_type(L, b, nc, opts->wsize, 0);
}

static void
Qs(w_latdirferm)(lua_State *L, mHdf5File *b, mLattice *S,
             struct wopts_s *opts, struct laddr_s *laddr,
             SHA256_Sum *sum, void **data, hid_t *filetype, hid_t *memtype,
             const char **kind)
{
  Qs(mLatDirFerm) *m = Qs(qlua_checkLatDirFerm)(L, 3, S, -1);
  int nc = QC(m);
  int *local_x = qlua_malloc(L, laddr->rank * sizeof (int));
  SHA256_Context *ctx = sha256_create(L);
  int volume = laddr->volume;
  int rank = laddr->rank;
  int dsize, i;
#if QNc == 'N'
  typedef QLA_DN_DiracFermion(nc, Vtype);
#else
  typedef Qx(QLA_D,_DiracFermion) Vtype;
#endif

  switch (opts->wsize) {
  case WS_Double:
    dsize = sizeof (machine_complex_double);
    break;
  case WS_Float:
    dsize = sizeof (machine_complex_float);
    break;
  default:
    luaL_error(L, "Unknown precision in w_dirferm()");
  }
  *data = qlua_malloc(L, dsize * nc * QDP_Ns * volume);
  CALL_QDP(L);
  Vtype *locked = Qx(QDP_D,_expose_D)(m->ptr);
  char *out_mem = *data;
  for (i = 0; i < volume; i++) {
    qdp2hdf5_addr(local_x, i, laddr);
    QLUA_ASSERT(QDP_node_number_L(S->lat, local_x) == QDP_this_node);
    Vtype *ptr = &locked[QDP_index_L(S->lat, local_x)];
    sha256_reset(ctx);
    sha256_sum_add_ints(ctx, &rank, 1);
    sha256_sum_add_ints(ctx, local_x, rank);
    Qs(wpack_dirferm)(L, out_mem, ctx, ptr, nc, opts->wsize);
    out_mem += dsize * nc * QDP_Ns;
    SHA256_Sum l_sum;
    sha256_sum(&l_sum, ctx);
    local_combine_checksums(sum, &l_sum);
  }
  Qx(QDP_D, _reset_D)(m->ptr);
  sha256_destroy(ctx);
  qlua_free(L, local_x);
  *kind = knLatticeDiracFermion;
  *filetype = get_dirferm_type(L, b, nc, opts->wsize, 1);
  *memtype  = get_dirferm_type(L, b, nc, opts->wsize, 0);
}

// Dirac Propagators
static void
Qs(wpack_dirprop)(lua_State *L, void *data, SHA256_Context *ctx, void *src, int nc, WriteSize wsize)
{
#if QNc == 'N'
  typedef QLA_DN_DiracPropagator(nc, DataType);
#else
  typedef Qx(QLA_D,_DiracPropagator) DataType;
#endif
  DataType *ptr = (DataType *)src;
  int ci, cj, di, dj;

  switch (wsize) {
  case WS_Double: {
    machine_complex_double *dst = (machine_complex_double *)data;
    for (di = 0; di < QDP_Ns; di++) {
      for (dj = 0; dj < QDP_Ns; dj++) {
        for (ci = 0; ci < nc; ci++) {
          for (cj = 0; cj < nc; cj++, dst++) {
            QLA_D_Complex zz;
            Qx(QLA_D,_C_eq_elem_P)(QNC(nc) &zz, ptr, ci, di, cj, dj);
            dst->re = QLA_real(zz);
            dst->im = QLA_imag(zz);
          }
        }
      }
    }
    sha256_sum_add_doubles(ctx, data, nc * nc * QDP_Ns * QDP_Ns * 2);
  } break;
  case WS_Float: {
    machine_complex_float *dst = (machine_complex_float *)data;
    for (di = 0; di < QDP_Ns; di++) {
      for (dj = 0; dj < QDP_Ns; dj++) {
        for (ci = 0; ci < nc; ci++) {
          for (cj = 0; cj < nc; cj++, dst++) {
            QLA_D_Complex zz;
            Qx(QLA_D,_C_eq_elem_P)(QNC(nc) &zz, ptr, ci, di, cj, dj);
            dst->re = QLA_real(zz);
            dst->im = QLA_imag(zz);
          }
        }
      }
    }
    sha256_sum_add_floats(ctx, data, nc * nc * QDP_Ns * QDP_Ns * 2);
  } break;
  default:
    luaL_error(L, "unknown wsize in wpack_dirpropN()");
    break;
  }
}

static void
Qs(w_dirprop)(lua_State *L, mHdf5File *b, mLattice *S,
             struct wopts_s *opts, struct laddr_s *laddr,
             SHA256_Sum *sum, void **data, hid_t *filetype, hid_t *memtype,
             const char **kind)
{
  Qs(mSeqDirProp) *m = Qs(qlua_checkSeqDirProp)(L, 3, -1);
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
    luaL_error(L, "Unknown precision in w_dirprop()");
  }
  *data = qlua_malloc(L, dsize * nc * nc * QDP_Ns * QDP_Ns);
  SHA256_Context *ctx = sha256_create(L);
  Qs(wpack_dirprop)(L, *data, ctx, m->ptr, nc, opts->wsize);
  sha256_sum(sum, ctx);
  sha256_destroy(ctx);
  *kind = knDiracPropagator;
  *filetype = get_dirprop_type(L, b, nc, opts->wsize, 1);
  *memtype  = get_dirprop_type(L, b, nc, opts->wsize, 0);
}

static void
Qs(w_latdirprop)(lua_State *L, mHdf5File *b, mLattice *S,
             struct wopts_s *opts, struct laddr_s *laddr,
             SHA256_Sum *sum, void **data, hid_t *filetype, hid_t *memtype,
             const char **kind)
{
  Qs(mLatDirProp) *m = Qs(qlua_checkLatDirProp)(L, 3, S, -1);
  int nc = QC(m);
  int *local_x = qlua_malloc(L, laddr->rank * sizeof (int));
  SHA256_Context *ctx = sha256_create(L);
  int volume = laddr->volume;
  int rank = laddr->rank;
  int dsize, i;
#if QNc == 'N'
  typedef QLA_DN_DiracPropagator(nc, Vtype);
#else
  typedef Qx(QLA_D,_DiracPropagator) Vtype;
#endif

  switch (opts->wsize) {
  case WS_Double:
    dsize = sizeof (machine_complex_double);
    break;
  case WS_Float:
    dsize = sizeof (machine_complex_float);
    break;
  default:
    luaL_error(L, "Unknown precision in w_dirprop()");
  }
  *data = qlua_malloc(L, dsize * nc * nc * QDP_Ns * QDP_Ns * volume);
  CALL_QDP(L);
  Vtype *locked = Qx(QDP_D,_expose_P)(m->ptr);
  char *out_mem = *data;
  for (i = 0; i < volume; i++) {
    qdp2hdf5_addr(local_x, i, laddr);
    QLUA_ASSERT(QDP_node_number_L(S->lat, local_x) == QDP_this_node);
    Vtype *ptr = &locked[QDP_index_L(S->lat, local_x)];
    sha256_reset(ctx);
    sha256_sum_add_ints(ctx, &rank, 1);
    sha256_sum_add_ints(ctx, local_x, rank);
    Qs(wpack_dirprop)(L, out_mem, ctx, ptr, nc, opts->wsize);
    out_mem += dsize * nc * nc * QDP_Ns * QDP_Ns;
    SHA256_Sum l_sum;
    sha256_sum(&l_sum, ctx);
    local_combine_checksums(sum, &l_sum);
  }
  Qx(QDP_D, _reset_P)(m->ptr);
  sha256_destroy(ctx);
  qlua_free(L, local_x);
  *kind = knLatticeDiracPropagator;
  *filetype = get_dirprop_type(L, b, nc, opts->wsize, 1);
  *memtype  = get_dirprop_type(L, b, nc, opts->wsize, 0);
}


#undef QNc
#undef Qcolors
#undef Qs
#undef Qx
#undef QC
#undef QNC
