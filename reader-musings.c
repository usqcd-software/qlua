
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

