
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
check_colvec_type(lua_State *L, const char *path, hid_t tobj, WriteSize *wsize, int *Nc)
{
  return check_veccomplex_type(L, path, tobj, wsize, Nc);
}

static int
check_colmat_type(lua_State *L, const char *path, hid_t tobj, WriteSize *wsize, int *Nc)
{
  int la, lb;
  if (check_matcomplex_type(L, path, tobj, wsize, *la, *lb) == 0)
    return 0;
  if (la != lb)
    return 0;
  *Nc = la;
  return 1;
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
static int
r_latcomplex(lua_State *L, mHdf5File *b, const char *path,
             struct ropts_s *ropts, hid_t obj, hid_t tobj, hid_t memspace, hid_t filespace,
             SHA256_Sum *sum, struct laddr_s *laddr)
{
  WriteSize wsize;
  if (!check_complex_type(L, path, tobj, &wsize))
    return 0;
  hid_t memtype = get_complex_type(L, b, wsize, 0);
  int volume = laddr->volume;
  int rank = ropts->S->rank;
  mLatComplex *m = qlua_newLatComplex(L, ropts->Sidx);
  QLA_D_Complex *locked = QDP_expose_C(m->ptr);
  SHA256_Context *ctx = sha256_create(L);
  int *local_x = qlua_malloc(L, ropts->S->rank * sizeof (int));
  herr_t status;
  switch (wsize) {
  case WS_Double: {
    machine_complex_double *ptr = qlua_malloc(L, volume * sizeof (machine_complex_double));
    status = H5Dread(obj, memtype, memspace, filespace, H5P_DEFAULT, ptr);
    int i;
    for (i = 0; i < volume; i++) {
      qdp2hdf5_addr(local_x, i, laddr);
      QLUA_ASSERT(QDP_node_number_L(ropts->S->lat, local_x) == QDP_this_node);
      QLA_D_Complex zz;
      QLA_real(zz) = ptr[i].re;
      QLA_imag(zz) = ptr[i].im;
      QLA_elem_C(locked[QDP_index_L(ropts->S->lat, local_x)]) = zz;
      sha256_reset(ctx);
      sha256_sum_add_ints(ctx, &rank, 1);
      sha256_sum_add_ints(ctx, local_x, rank);
      sha256_sum_add_doubles(ctx, (double *)&ptr[i], 2);
      SHA256_Sum l_sum;
      sha256_sum(&l_sum, ctx);
      local_combine_checksums(sum, &l_sum);
    }
    qlua_free(L, ptr);
  } break;
  case WS_Float: {
    machine_complex_float *ptr = qlua_malloc(L, volume * sizeof (machine_complex_float));
    status = H5Dread(obj, memtype, memspace, filespace, H5P_DEFAULT, ptr);
    int i;
    for (i = 0; i < volume; i++) {
      qdp2hdf5_addr(local_x, i, laddr);
      QLUA_ASSERT(QDP_node_number_L(ropts->S->lat, local_x) == QDP_this_node);
      QLA_D_Complex zz;
      QLA_real(zz) = ptr[i].re;
      QLA_imag(zz) = ptr[i].im;
      QLA_elem_C(locked[QDP_index_L(ropts->S->lat, local_x)]) = zz;
      sha256_reset(ctx);
      sha256_sum_add_ints(ctx, &rank, 1);
      sha256_sum_add_ints(ctx, local_x, rank);
      sha256_sum_add_floats(ctx, (float *)&ptr[i], 2);
      SHA256_Sum l_sum;
      sha256_sum(&l_sum, ctx);
      local_combine_checksums(sum, &l_sum);
    }
    qlua_free(L, ptr);
  } break;
  default:
    QLUA_ABORT("Unknown precision in r_latcomplex()");
  }
  CHECK_H5(L, H5Tclose(memtype), "Tclose() memory type");
  QDP_reset_C(m->ptr);
  sha256_destroy(ctx);
  qlua_free(L, local_x);
  if (status < 0) {
    lua_pop(L, 1);
    return 0;
  }
  return 1;
}

