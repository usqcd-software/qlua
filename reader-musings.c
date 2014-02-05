
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
static void
w_latint(lua_State *L, mHdf5File *b, mLattice *S,
         struct wopts_s *opts, struct laddr_s *laddr,
         SHA256_Sum *sum, void **data, hid_t *filetype, hid_t *memtype,
         const char **kind)
{
  int *local_x = qlua_malloc(L, laddr->rank * sizeof (int));
  SHA256_Context *ctx = sha256_create(L);
  int volume = laddr->volume;
  int rank = laddr->rank;
  int i;

  int *ptr = qlua_malloc(L, volume * sizeof (int));
  CALL_QDP(L);
  mLatInt *m = qlua_checkLatInt(L, 3, S);
  QLA_Int *locked = QDP_expose_I(m->ptr);
  for (i = 0; i < volume; i++) {
    qdp2hdf5_addr(local_x, i, laddr);
    QLUA_ASSERT(QDP_node_number_L(S->lat, local_x) == QDP_this_node);
    ptr[i] = QLA_elem_I(locked[QDP_index_L(S->lat, local_x)]);
    sha256_reset(ctx);
    sha256_sum_add_ints(ctx, &rank, 1);
    sha256_sum_add_ints(ctx, local_x, rank);
    sha256_sum_add_ints(ctx, &ptr[i], 1);
    SHA256_Sum l_sum;
    sha256_sum(&l_sum, ctx);
    local_combine_checksums(sum, &l_sum);
  }
  QDP_reset_I(m->ptr);
  *data = ptr;

  sha256_destroy(ctx);
  qlua_free(L, local_x);
  *kind = knLatticeInt;
  *filetype = get_int_type(L, b, 1);
  *memtype  = get_int_type(L, b, 0);
}

static int
r_latint(lua_State *L, mHdf5File *b, const char *path,
         struct ropts_s *ropts, hid_t obj, hid_t tobj, hid_t memspace, hid_t filespace,
         SHA256_Sum *sum, struct laddr_s *laddr)
{
  if (!check_int_type(L, path, tobj))
    return 0;
  int *ptr = qlua_malloc(L, volume * sizeof (int));
  // XXX may be wrong: what steps are needed for independent read?
  if (H5Dread(obj, memtype, H5S_ALL, dspace, H5P_DEFAULT, ptr) < 0) {
    qlua_free(L, ptr);
    return 0;
  mLatInt *m = qlua_newLatInt(L, ropts->Sidx);
  QLA_Int *locked = QDP_expose_I(m->ptr);
  SHA256_Context *ctx = sha256_create(L);
  int *local_x = qlua_malloc(L, ropts->rank * sizeof (int));
  int i;
  for (i = 0; i < laddr->volume; i++) {
    qdp2hdf5_addr(local_x, i, laddr);
    QLUA_ASSERT(QDP_node_number_L(S->lat, local_x) == QDP_this_node);
    QLA_elem_I(locked[QDP_index_L(S->lat, local_x)]) = ptr[i];
    sha256_reset(ctx);
    sha256_sum_add_ints(ctx, &rank, 1);
    sha256_sum_add_ints(ctx, local_x, rank);
    sha256_sum_add_ints(ctx, &ptr[i], 1);
    SHA256_Sum l_sum;
    sha256_sum(&l_sum, ctx);
    local_combine_checksums(sum, &l_sum);
  }
  QDP_reset_I(m->ptr);
  sha256_destroy(ctx);
  qlua_free(L, local_x);
  qlua_free(L, ptr);
  return 1;
}

static int
r_string(lua_State *L, mLattice *S, mHdf5File *b, const char *path,
         struct ropts_s *ropts, hid_t obj, hid_t tobj, SHA256_Sum *sum)
{
  int len;
  if (!check_string_type(L, path, tobj, &len))
    return 0;
  hid_t memtype = get_string_type(L, b, len, 0);
  char *buffer = qlua_malloc(L, len + 1);
  herr_t status = H5Dread(obj, memtype, H5S_ALL, H5S_ALL, H5P_DEFAULT, buffer);
  CHECK_H5p(L, H5Tclose(memtype), "Tclose() failed in read(\"%s\")", path);
  if (status < 0) {
    qlua_free(L, buffer);
    return 0;
  }
  buffer[len] = 0;
  sha256_sum_string(sum, buffer, len);
  lua_pushstring(L, buffer);
  return 1;
}

