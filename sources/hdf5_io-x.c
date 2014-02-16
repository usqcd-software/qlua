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
    luaL_error(L, "unknown wsize in wpack_colvecX()");
    break;
  }
}

static void
Qs(unpack_colvecD)(void *dst, int nc, const machine_complex_double *src, SHA256_Context *ctx)
{
#if QNc == 'N'
  typedef QLA_DN_ColorVector(nc, Vtype);
#else
  typedef Qx(QLA_D,_ColorVector) Vtype;
#endif
  Vtype *ptr = dst;
  int c;
  sha256_sum_add_doubles(ctx, (double *)src, 2 * nc);
  for (c = 0; c < nc; c++) {
    QLA_D_Complex zz;
    QLA_real(zz) = src->re;
    QLA_imag(zz) = src->im;
    Qx(QLA_D,_V_eq_elem_C)(QNC(nc) ptr, &zz, c);
    src++;
  }
}

static void
Qs(unpack_colvecF)(void *dst, int nc, const machine_complex_float *src, SHA256_Context *ctx)
{
#if QNc == 'N'
  typedef QLA_DN_ColorVector(nc, Vtype);
#else
  typedef Qx(QLA_D,_ColorVector) Vtype;
#endif
  Vtype *ptr = dst;
  int c;
  sha256_sum_add_floats(ctx, (float *)src, 2 * nc);
  for (c = 0; c < nc; c++) {
    QLA_D_Complex zz;
    QLA_real(zz) = src->re;
    QLA_imag(zz) = src->im;
    Qx(QLA_D,_V_eq_elem_C)(QNC(nc) ptr, &zz, c);
    src++;
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
    dsize = 0;
    luaL_error(L, "Unknown precision in w_colvec()");
    break;
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
Qs(r_colvec)(lua_State *L, struct ropts_s *ropts,
             SHA256_Context *ctx,
             int nc, WriteSize wsize, void *data)
{
  Qs(mSeqColVec) *m = Qs(qlua_newSeqColVec)(L, nc);
  switch (wsize) {
  case WS_Double:
    Qs(unpack_colvecD)(m->ptr, nc, data, ctx);
    break;
  case WS_Float:
    Qs(unpack_colvecF)(m->ptr, nc, data, ctx);
    break;
  default:
    QLUA_ABORT("unknown precision in r_colvecX()");
  }
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
    dsize = 0;
    luaL_error(L, "Unknown precision in w_latcolvec()");
    break;
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

static void
Qs(r_latcolvec)(lua_State *L, struct ropts_s *ropts,
                struct laddr_s *laddr, int *local_x,
                SHA256_Context *ctx, SHA256_Sum *sum,
                int nc, WriteSize wsize, void *data)
{
#if QNc == 'N'
  typedef QLA_DN_ColorVector(nc, Vtype);
#else
  typedef Qx(QLA_D,_ColorVector) Vtype;
#endif
  Qs(mLatColVec) *m = Qs(qlua_newLatColVec)(L, ropts->Sidx, nc);
  Vtype *dst = Qx(QDP_D, _expose_V)(m->ptr);
  int volume = laddr->volume;
  switch (wsize) {
  case WS_Double: {
    machine_complex_double *src = data;
    SHA256_Sum l_sum;
    int i;
    for (i = 0; i < volume; i++) {
      qdp2hdf5_addr(local_x, i, laddr);
      sha256_reset(ctx);
      sha256_sum_add_ints(ctx, &laddr->rank, 1);
      sha256_sum_add_ints(ctx, local_x, laddr->rank);
      Qs(unpack_colvecD)(&dst[QDP_index_L(ropts->S->lat, local_x)], nc, src, ctx);
      sha256_sum(&l_sum, ctx);
      local_combine_checksums(sum, &l_sum);
      src += nc;
    }
  } break;
  case WS_Float: {
    machine_complex_float *src = data;
    SHA256_Sum l_sum;
    int i;
    for (i = 0; i < volume; i++) {
      qdp2hdf5_addr(local_x, i, laddr);
      sha256_reset(ctx);
      sha256_sum_add_ints(ctx, &laddr->rank, 1);
      sha256_sum_add_ints(ctx, local_x, laddr->rank);
      Qs(unpack_colvecF)(&dst[QDP_index_L(ropts->S->lat, local_x)], nc, src, ctx);
      sha256_sum(&l_sum, ctx);
      local_combine_checksums(sum, &l_sum);
      src += nc;
    }
  } break;
  default:
    QLUA_ABORT("unknown precision in r_latcolvecX()");
  }
  Qx(QDP_D, _reset_V)(m->ptr);
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
    sha256_sum_add_doubles(ctx, data, 2 * nc * nc);
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
    sha256_sum_add_floats(ctx, data, 2 * nc * nc);
  } break;
  default:
    luaL_error(L, "unknown wsize in wpack_colmatX()");
    break;
  }
}

static void
Qs(unpack_colmatD)(void *dst, int nc, const machine_complex_double *src, SHA256_Context *ctx)
{
#if QNc == 'N'
  typedef QLA_DN_ColorMatrix(nc, Vtype);
#else
  typedef Qx(QLA_D,_ColorMatrix) Vtype;
#endif
  Vtype *ptr = dst;
  int ci, cj;
  sha256_sum_add_doubles(ctx, (double *)src, 2 * nc * nc);
  for (ci = 0; ci < nc; ci++) {
    for (cj = 0; cj < nc; cj++) {
      QLA_D_Complex zz;
      QLA_real(zz) = src->re;
      QLA_imag(zz) = src->im;
      Qx(QLA_D,_M_eq_elem_C)(QNC(nc) ptr, &zz, ci, cj);
      src++;
    }
  }
}

static void
Qs(unpack_colmatF)(void *dst, int nc, const machine_complex_float *src, SHA256_Context *ctx)
{
#if QNc == 'N'
  typedef QLA_DN_ColorMatrix(nc, Vtype);
#else
  typedef Qx(QLA_D,_ColorMatrix) Vtype;
#endif
  Vtype *ptr = dst;
  int ci, cj;
  sha256_sum_add_floats(ctx, (float *)src, 2 * nc * nc);
  for (ci = 0; ci < nc; ci++) {
    for (cj = 0; cj < nc; cj++) {
      QLA_D_Complex zz;
      QLA_real(zz) = src->re;
      QLA_imag(zz) = src->im;
      Qx(QLA_D,_M_eq_elem_C)(QNC(nc) ptr, &zz, ci, cj);
      src++;
    }
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
    dsize = 0;
    luaL_error(L, "Unknown precision in w_colmat()");
    break;
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
Qs(r_colmat)(lua_State *L, struct ropts_s *ropts,
             SHA256_Context *ctx,
             int nc, WriteSize wsize, void *data)
{
  Qs(mSeqColMat) *m = Qs(qlua_newSeqColMat)(L, nc);
  switch (wsize) {
  case WS_Double:
    Qs(unpack_colmatD)(m->ptr, nc, data, ctx);
    break;
  case WS_Float:
    Qs(unpack_colmatF)(m->ptr, nc, data, ctx);
    break;
  default:
    QLUA_ABORT("unknown precision in r_colmatX()");
  }
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
    dsize = 0;
    luaL_error(L, "Unknown precision in w_latcolmat()");
    break;
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

static void
Qs(r_latcolmat)(lua_State *L, struct ropts_s *ropts,
                struct laddr_s *laddr, int *local_x,
                SHA256_Context *ctx, SHA256_Sum *sum,
                int nc, WriteSize wsize, void *data)
{
#if QNc == 'N'
  typedef QLA_DN_ColorMatrix(nc, Vtype);
#else
  typedef Qx(QLA_D,_ColorMatrix) Vtype;
#endif
  Qs(mLatColMat) *m = Qs(qlua_newLatColMat)(L, ropts->Sidx, nc);
  Vtype *dst = Qx(QDP_D, _expose_M)(m->ptr);
  int volume = laddr->volume;
  switch (wsize) {
  case WS_Double: {
    machine_complex_double *src = data;
    SHA256_Sum l_sum;
    int i;
    for (i = 0; i < volume; i++) {
      qdp2hdf5_addr(local_x, i, laddr);
      sha256_reset(ctx);
      sha256_sum_add_ints(ctx, &laddr->rank, 1);
      sha256_sum_add_ints(ctx, local_x, laddr->rank);
      Qs(unpack_colmatD)(&dst[QDP_index_L(ropts->S->lat, local_x)], nc, src, ctx);
      sha256_sum(&l_sum, ctx);
      local_combine_checksums(sum, &l_sum);
      src += nc * nc;
    }
  } break;
  case WS_Float: {
    machine_complex_float *src = data;
    SHA256_Sum l_sum;
    int i;
    for (i = 0; i < volume; i++) {
      qdp2hdf5_addr(local_x, i, laddr);
      sha256_reset(ctx);
      sha256_sum_add_ints(ctx, &laddr->rank, 1);
      sha256_sum_add_ints(ctx, local_x, laddr->rank);
      Qs(unpack_colmatF)(&dst[QDP_index_L(ropts->S->lat, local_x)], nc, src, ctx);
      sha256_sum(&l_sum, ctx);
      local_combine_checksums(sum, &l_sum);
      src += nc * nc;
    }
  } break;
  default:
    QLUA_ABORT("unknown precision in r_latcolmatX()");
  }
  Qx(QDP_D, _reset_M)(m->ptr);
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
Qs(unpack_dirfermD)(void *dst, int nc, const machine_complex_double *src, SHA256_Context *ctx)
{
#if QNc == 'N'
  typedef QLA_DN_DiracFermion(nc, Vtype);
#else
  typedef Qx(QLA_D,_DiracFermion) Vtype;
#endif
  Vtype *ptr = dst;
  int d, c;
  sha256_sum_add_doubles(ctx, (double *)src, 2 * QDP_Ns * nc);
  for (d = 0; d < QDP_Ns; d++) {
    for (c = 0; c < nc; c++) {
      QLA_D_Complex zz;
      QLA_real(zz) = src->re;
      QLA_imag(zz) = src->im;
      Qx(QLA_D,_D_eq_elem_C)(QNC(nc) ptr, &zz, c, d);
      src++;
    }
  }
}

static void
Qs(unpack_dirfermF)(void *dst, int nc, const machine_complex_float *src, SHA256_Context *ctx)
{
#if QNc == 'N'
  typedef QLA_DN_DiracFermion(nc, Vtype);
#else
  typedef Qx(QLA_D,_DiracFermion) Vtype;
#endif
  Vtype *ptr = dst;
  int d, c;
  sha256_sum_add_floats(ctx, (float *)src, 2 * QDP_Ns * nc);
  for (d = 0; d < QDP_Ns; d++) {
    for (c = 0; c < nc; c++) {
      QLA_D_Complex zz;
      QLA_real(zz) = src->re;
      QLA_imag(zz) = src->im;
      Qx(QLA_D,_D_eq_elem_C)(QNC(nc) ptr, &zz, c, d);
      src++;
    }
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
    dsize = 0;
    luaL_error(L, "Unknown precision in w_dirferm()");
    break;
  }
  *data = qlua_malloc(L, QDP_Ns * dsize * nc);
  SHA256_Context *ctx = sha256_create(L);
  Qs(wpack_dirferm)(L, *data, ctx, m->ptr, nc, opts->wsize);
  sha256_sum(sum, ctx);
  sha256_destroy(ctx);
  *kind = knDiracFermion;
  *filetype = get_dirferm_type(L, b, nc, opts->wsize, 1);
  *memtype  = get_dirferm_type(L, b, nc, opts->wsize, 0);
}

static void
Qs(r_dirferm)(lua_State *L, struct ropts_s *ropts,
              SHA256_Context *ctx,
              int nc, WriteSize wsize, void *data)
{
  Qs(mSeqDirFerm) *m = Qs(qlua_newSeqDirFerm)(L, nc);
  switch (wsize) {
  case WS_Double:
    Qs(unpack_dirfermD)(m->ptr, nc, data, ctx);
    break;
  case WS_Float:
    Qs(unpack_dirfermF)(m->ptr, nc, data, ctx);
    break;
  default:
    QLUA_ABORT("unknown precision in r_dirfermX()");
  }
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
    dsize = 0;
    luaL_error(L, "Unknown precision in w_dirferm()");
    break;
  }
  *data = qlua_malloc(L, QDP_Ns * dsize * nc * volume);
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

static void
Qs(r_latdirferm)(lua_State *L, struct ropts_s *ropts,
                 struct laddr_s *laddr, int *local_x,
                 SHA256_Context *ctx, SHA256_Sum *sum,
                 int nc, WriteSize wsize, void *data)
{
#if QNc == 'N'
  typedef QLA_DN_DiracFermion(nc, Vtype);
#else
  typedef Qx(QLA_D,_DiracFermion) Vtype;
#endif
  Qs(mLatDirFerm) *m = Qs(qlua_newLatDirFerm)(L, ropts->Sidx, nc);
  Vtype *dst = Qx(QDP_D, _expose_D)(m->ptr);
  int volume = laddr->volume;
  switch (wsize) {
  case WS_Double: {
    machine_complex_double *src = data;
    SHA256_Sum l_sum;
    int i;
    for (i = 0; i < volume; i++) {
      qdp2hdf5_addr(local_x, i, laddr);
      sha256_reset(ctx);
      sha256_sum_add_ints(ctx, &laddr->rank, 1);
      sha256_sum_add_ints(ctx, local_x, laddr->rank);
      Qs(unpack_dirfermD)(&dst[QDP_index_L(ropts->S->lat, local_x)], nc, src, ctx);
      sha256_sum(&l_sum, ctx);
      local_combine_checksums(sum, &l_sum);
      src += QDP_Ns * nc;
    }
  } break;
  case WS_Float: {
    machine_complex_float *src = data;
    SHA256_Sum l_sum;
    int i;
    for (i = 0; i < volume; i++) {
      qdp2hdf5_addr(local_x, i, laddr);
      sha256_reset(ctx);
      sha256_sum_add_ints(ctx, &laddr->rank, 1);
      sha256_sum_add_ints(ctx, local_x, laddr->rank);
      Qs(unpack_dirfermF)(&dst[QDP_index_L(ropts->S->lat, local_x)], nc, src, ctx);
      sha256_sum(&l_sum, ctx);
      local_combine_checksums(sum, &l_sum);
      src += QDP_Ns * nc;
    }
  } break;
  default:
    QLUA_ABORT("unknown precision in r_latdirfermX()");
  }
  Qx(QDP_D, _reset_D)(m->ptr);
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
    sha256_sum_add_doubles(ctx, data, 2 * QDP_Ns * QDP_Ns * nc * nc);
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
    sha256_sum_add_floats(ctx, data, 2 * QDP_Ns * QDP_Ns * nc * nc);
  } break;
  default:
    luaL_error(L, "unknown wsize in wpack_dirpropN()");
    break;
  }
}

static void
Qs(unpack_dirpropD)(void *dst, int nc, const machine_complex_double *src, SHA256_Context *ctx)
{
#if QNc == 'N'
  typedef QLA_DN_DiracPropagator(nc, Vtype);
#else
  typedef Qx(QLA_D,_DiracPropagator) Vtype;
#endif
  Vtype *ptr = dst;
  int di, dj, ci, cj;
  sha256_sum_add_doubles(ctx, (double *)src, 2 * QDP_Ns * QDP_Ns * nc * nc);
  for (di = 0; di < QDP_Ns; di++) {
    for (dj = 0; dj < QDP_Ns; dj++) {
      for (ci = 0; ci < nc; ci++) {
        for (cj = 0; cj < nc; cj++) {
          QLA_D_Complex zz;
          QLA_real(zz) = src->re;
          QLA_imag(zz) = src->im;
          Qx(QLA_D,_P_eq_elem_C)(QNC(nc) ptr, &zz, ci, di, cj, dj);
          src++;
        }
      }
    }
  }
}

static void
Qs(unpack_dirpropF)(void *dst, int nc, const machine_complex_float *src, SHA256_Context *ctx)
{
#if QNc == 'N'
  typedef QLA_DN_DiracPropagator(nc, Vtype);
#else
  typedef Qx(QLA_D,_DiracPropagator) Vtype;
#endif
  Vtype *ptr = dst;
  int di, dj, ci, cj;
  sha256_sum_add_floats(ctx, (float *)src, 2 * QDP_Ns * QDP_Ns * nc * nc);
  for (di = 0; di < QDP_Ns; di++) {
    for (dj = 0; dj < QDP_Ns; dj++) {
      for (ci = 0; ci < nc; ci++) {
        for (cj = 0; cj < nc; cj++) {
          QLA_D_Complex zz;
          QLA_real(zz) = src->re;
          QLA_imag(zz) = src->im;
          Qx(QLA_D,_P_eq_elem_C)(QNC(nc) ptr, &zz, ci, di, cj, dj);
          src++;
        }
      }
    }
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
    dsize = 0;
    luaL_error(L, "Unknown precision in w_dirprop()");
    break;
  }
  *data = qlua_malloc(L, QDP_Ns * QDP_Ns * dsize * nc * nc);
  SHA256_Context *ctx = sha256_create(L);
  Qs(wpack_dirprop)(L, *data, ctx, m->ptr, nc, opts->wsize);
  sha256_sum(sum, ctx);
  sha256_destroy(ctx);
  *kind = knDiracPropagator;
  *filetype = get_dirprop_type(L, b, nc, opts->wsize, 1);
  *memtype  = get_dirprop_type(L, b, nc, opts->wsize, 0);
}

static void
Qs(r_dirprop)(lua_State *L, struct ropts_s *ropts,
             SHA256_Context *ctx,
             int nc, WriteSize wsize, void *data)
{
  Qs(mSeqDirProp) *m = Qs(qlua_newSeqDirProp)(L, nc);
  switch (wsize) {
  case WS_Double:
    Qs(unpack_dirpropD)(m->ptr, nc, data, ctx);
    break;
  case WS_Float:
    Qs(unpack_dirpropF)(m->ptr, nc, data, ctx);
    break;
  default:
    QLUA_ABORT("unknown precision in r_dirpropX()");
  }
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
    dsize = 0;
    luaL_error(L, "Unknown precision in w_dirprop()");
    break;
  }
  *data = qlua_malloc(L, QDP_Ns * QDP_Ns * dsize * nc * nc * volume);
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

static void
Qs(r_latdirprop)(lua_State *L, struct ropts_s *ropts,
                struct laddr_s *laddr, int *local_x,
                SHA256_Context *ctx, SHA256_Sum *sum,
                int nc, WriteSize wsize, void *data)
{
#if QNc == 'N'
  typedef QLA_DN_DiracPropagator(nc, Vtype);
#else
  typedef Qx(QLA_D,_DiracPropagator) Vtype;
#endif
  Qs(mLatDirProp) *m = Qs(qlua_newLatDirProp)(L, ropts->Sidx, nc);
  Vtype *dst = Qx(QDP_D, _expose_P)(m->ptr);
  int volume = laddr->volume;
  switch (wsize) {
  case WS_Double: {
    machine_complex_double *src = data;
    SHA256_Sum l_sum;
    int i;
    for (i = 0; i < volume; i++) {
      qdp2hdf5_addr(local_x, i, laddr);
      sha256_reset(ctx);
      sha256_sum_add_ints(ctx, &laddr->rank, 1);
      sha256_sum_add_ints(ctx, local_x, laddr->rank);
      Qs(unpack_dirpropD)(&dst[QDP_index_L(ropts->S->lat, local_x)], nc, src, ctx);
      sha256_sum(&l_sum, ctx);
      local_combine_checksums(sum, &l_sum);
      src += QDP_Ns * QDP_Ns * nc * nc;
    }
  } break;
  case WS_Float: {
    machine_complex_float *src = data;
    SHA256_Sum l_sum;
    int i;
    for (i = 0; i < volume; i++) {
      qdp2hdf5_addr(local_x, i, laddr);
      sha256_reset(ctx);
      sha256_sum_add_ints(ctx, &laddr->rank, 1);
      sha256_sum_add_ints(ctx, local_x, laddr->rank);
      Qs(unpack_dirpropF)(&dst[QDP_index_L(ropts->S->lat, local_x)], nc, src, ctx);
      sha256_sum(&l_sum, ctx);
      local_combine_checksums(sum, &l_sum);
      src += QDP_Ns * QDP_Ns * nc * nc;
    }
  } break;
  default:
    QLUA_ABORT("unknown precision in r_latdirpropX()");
  }
  Qx(QDP_D, _reset_P)(m->ptr);
}

#undef QNc
#undef Qcolors
#undef Qs
#undef Qx
#undef QC
#undef QNC
