
static int
check_time_type(lua_State *L, const char *path, hid_t tobj)
{
  return ((H5Tget_class(tobj) == H5T_INTEGER) && (H5Tget_size(tobj) == 8));
}

static int
check_sha256_type(lua_State *L, const char *path, hid_t tobj)
{
  if ((H5Tget_class(tobj) != H5T_ARRAY) || (H5Tget_array_ndims(tobj) != 1))
    return 0;
  hsize_t dim;
  H5Tget_array_dims2(tobj, &dim);
  if (dim != sizeof (SHA256_Sum))
    return 0;
  hid_t base = H5Tget_super(tobj);
  CHECK_H5p(L, base, "Tget_super() failed in read(\"%s\")", path);
  int status = (H5Tget_class(base) == H5T_INTEGER);
  CHECK_H5p(L, H5Tclose(base), "Tclose(base) failed in read(\"%s\")", path);
  return status;
}

static int
check_dirferm_type(lua_State *L, const char *path, hid_t tobj, WriteSize *wsize, int *Nc)
{
  if ((H5Tget_class(tobj) != H5T_ARRAY) || (H5Tget_array_ndims(tobj) != 1))
    return 0;
  hsize_t dim;
  H5Tget_array_dim2(tobj, &dim);
  if (dim != QDP_Ns)
    return 0;
  hid_t te = H5Tget_super(tobj);
  CHECK_H5p(L, te, "Tget_super() failed in read(\"%s\")", path);
  int status = check_veccomplex_type(L, path, te, wsize, Nc);
  CHECK_H5p(L, H5Tclose(te), "Tclose() failed in read(\"%s\")", path);
  return status;
}

static int
check_dirprop_type(lua_State *L, const char *path, hid_t tobj, WriteSize *wsize, int *Nc)
{
  if ((H5Tget_class(tobj) != H5T_ARRAY) || (H5Tget_array_ndims(tobj) != 2))
    return 0;
  hsize_t dims[2];
  H5Tget_array_dim2(tobj, dims);
  if ((dims[0] != QDP_Ns) || (dims[1] !+ QDP_Ns))
    return 0;
  hid_t te = H5Tget_super(tobj);
  CHECK_H5p(L, te, "Tget_super() failed in read(\"%s\")", path);
  int status = check_matcomplex_type(L, path, te, wsize, Nc);
  CHECK_H5p(L, H5Tclose(te), "Tclose() failed in read(\"%s\")", path);
  return status;
}

////////////////////// read template
static void
Qs(unpack_colmatD)(void *dst, int nc, const machine_complex_double *src,
                   SHA256_Context *ctx, SHA256_Sum *sum)
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
Qs(unpack_colmatF)(void *dst, int nc, const machine_complex_float *src,
                   SHA256_Context *ctx, SHA256_Sum *sum)
{
#if QNc == 'N'
  typedef QLA_DN_ColorMatrix(nc, Vtype);
#else
  typedef Qx(QLA_D,_ColorMatrix) Vtype;
#endif
  Vtype *ptr = dst;
  int ci, cj;
  sha256_sum_add_floats(ctx, (floats *)src, 2 * nc * nc);
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
Qs(r_colmat)(lua_State *L, struct ropts_s *ropts,
             SHA256_Context *ctx, SHA256_Sum *sum,
             int nc, WriteSize wsize, void *data)
{
  Qs(mSeqColMat) *m = Qs(qlua_newSeqColMat)(L, nc);
  switch (wsize) {
  case WS_Double:
    Qs(unpack_colmatD)(m->ptr, nc, data, ctx, sum);
    break;
  case WS_Float:
    Qs(unpack_colmatF)(m->ptr, nc, data, ctx, sum);
    break;
  default:
    QLUA_ABORT("unknown precision in r_colmatX()");
  }
}    

static void
Qs(r_latcolmat)(lua_State *L, struct ropts_s *ropts,
                struct laddr_s *laddr, int *local_x,
                SHA256_Context *ctx, SHA256_Sum *sum,
                int nc, WriteSize wsize, void *data)
{
#if QNc == 'N'
  typedef QLA_DN_ColorVector(nc, Vtype);
#else
  typedef Qx(QLA_D,_ColorVector) Vtype;
#endif
  Qs(mLatColMat) *m = Qs(qlua_newLatColMat)(L, ropts->Sidx, nc);
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
      Qs(unpack_colmatD)(&dst[QDP_index_L(ropts->S->lat, local_x)], nc, src, ctx, &l_sum);
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
      Qs(unpack_colmatF)(&dst[QDP_index_L(ropts->S->lat, local_x)], nc, src, ctx, &l_sum);
      sha256_sum(&l_sum, ctx);
      local_combine_checksums(sum, &l_sum);
      src += nc;
    }
  } break;
  default:
    QLUA_ABORT("unknown precision in r_latcolmatX()");
  }
  Qx(QDP_D, _reset_V)(m->ptr);
}

// Not in -x.c
static int
r_colmat(lua_State *L, mHdf5File *b, const char *path,
         struct ropts_s *ropts, hid_t obj, hid_t tobj, hid_t memspace, hid_t filespace,
         SHA256_Sum *sum, struct laddr_s *laddr)
{
  int nc;
  WriteSize wsize;
  if (!check_colmat_type(L, path, tobj, &nc, &wsize))
    return 0;
  void *data;
  switch (wsize) {
  case WS_Double:
    data = qlua_malloc(L, nc * sizeof (machine_complex_double));
    break;
  case WS_Float:
    data = qlua_malloc(L, nc * sizeof (machine_complex_float));
    break;
  default:
    QLUA_ABORT("Unknown precision in r_colmat()");
  }
  hid_t memtype = get_colmat_type(L, b, nc, wsize, 0);
  herr_t status = H5Dread(obj, memtype, memspace, filespace, H5P_DEFAULT, data);
  CHECK_H5p(L, H5Tclose(memtype), "Tclose() mem type in r_colmat(\"%s\")", path);
  if (status < 0) {
    qlua_free(L, data);
    return 0;
  }
  SHA256_Context *ctx = sha256_create(L);
  switch (nc) {
#if USE_Nc2
  case 2:
    r_colmat2(L, ropts, ctx, sum, nc, wsize, data);
    break;
#endif
#if USE_Nc3
  case 3:
    r_colmat3(L, ropts, ctx, sum, nc, wsize, data);
    break;
#endif
  default:
#if USE_NcN
    r_colmatN(L, ropts, ctx, sum, nc, wsize, data);
#else
    luaL_error(L, "Unsupported Nc=%d in r_colmat()", nc);
#endif
    break;
  }
  sha256_sum(sum, ctx);
  sha256_destroy(ctx);
  qlua_free(L, data);
  return 1;
}

static int
r_latcolmat(lua_State *L, mHdf5File *b, const char *path,
            struct ropts_s *ropts, hid_t obj, hid_t tobj, hid_t memspace, hid_t filespace,
            SHA256_Sum *sum, struct laddr_s *laddr)
{
  int volume = laddr->volume;
  int nc;
  WriteSize wsize;
  if (!check_colmat_type(L, path, tobj, &nc, &wsize))
    return 0;
  void *data;
  switch (wsize) {
  case WS_Double:
    data = qlua_malloc(L, nc * volume * sizeof (machine_complex_double));
    break;
  case WS_Float:
    data = qlua_malloc(L, nc * volume * sizeof (machine_complex_float));
    break;
  default:
    QLUA_ABORT("Unknown precision in r_latcolmat()");
  }
  hid_t memtype = get_colmat_type(L, b, nc, wsize, 0);
  herr_t status = H5Dread(obj, memtype, memspace, filespace, H5P_DEFAULT, data);
  CHECK_H5p(L, H5Tclose(memtype), "Tclose() mem type in r_latcolmat(\"%s\")", path);
  if (status < 0) {
    qlua_free(L, data);
    return 0;
  }
  SHA256_Context *ctx = sha256_create(L);
  int *local_x = qlua_malloc(L, laddr->rank * sizeof (int));
  switch (nc) {
#if USE_Nc2
  case 2:
    r_latcolmat2(L, ropts, laddr, local_x, ctx, sum, nc, wsize, data);
    break;
#endif
#if USE_Nc3
  case 3:
    r_latcolmat3(L, ropts, laddr, local_x, ctx, sum, nc, wsize, data);
    break;
#endif
  default:
#if USE_NcN
    r_latcolmatN(L, ropts, laddr, local_x, ctx, sum, nc, wsize, data);
#else
    luaL_error(L, "Unsupported Nc=%d in r_latcolmat()", nc);
#endif
    break;
  }
  sha256_destroy(ctx);
  qlua_free(L, local_x);
  qlua_free(L, data);
  return 1;
}
