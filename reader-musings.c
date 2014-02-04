#define CHECK_H5p(L, expr, message, path) do { if ((expr) < 0) luaL_error(L, message, path); } while (0)

static int
check_time_type(lua_State *L, const char *path, hid_t tobj)
{
  return ((H5Tget_class(tobj) == H5T_INTEGER) && (H5Tget_size(tobj) == 8));
}

static int
check_int_type(lua_State *L, const char *path, hid_t tobj)
{
  return ((H5Tget_class(tobj) == H5T_INTEGER) && (H5Tget_size(tobj) == 4));
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
check_string_type(lua_State *L, const char *path, hid_t tobj, int *len)
{
  if (H5Tget_class(tobj) != H5T_STRING)
    return 0;
  *len = H5Tget_size(tobj);
  return 1;
}

static int
check_real_type(lua_State *L, const char *path, hid_t tobj, WriteSize *wsize)
{
  if (H5Tget_class(tobj) != H5T_STRING)
    return 0;
  switch (H5Tget_size(tobj)) {
  case 8: *wsize = WS_Double; return 1;
  case 4: *wsize = WS_Float; return 1;
  }
  return 0;
}

static int
check_complex_type(lua_State *L, const char *path, hid_t tobj, WriteSize *wsize)
{
  if ((H5Tget_class(tobj) != H5T_COMPOUND) || (H5Tget_nmembers(tobj) != 2))
    return 0;
  WriteSize s[2];
  int has_r = 0;
  int has_i = 0;
  int i;
  for (i = 0; i < 2; i++) {
    char *name = H5Tget_member_name(tobj, i);
    if (strcmp(name, "r") == 0)
      has_r = 1;
    if (strcmp(name, "i") == 0)
      has_i = 1;
    free(name);
    hid_t et = H5T_get_member_type(tobj, i);
    CHECK_H5p(L, et, "Tget_member_type() failed in read(\"%s\")", path);
    int status = check_real_type(L, et, &s[i]);
    CHECK_H5p(L, H5Tclose(et), "Tclose() failed in read(\"%s\")", path);
    if (status == 0)
      return 0;
  }
  if ((has_r == 0) || (has_i == 0) || (s[0] != s[1]))
    return 0;
  *wsize = s[0];
  return 1;
}

static int
check_vecint_type(lua_State *L, const char *path, hid_t tobj, int *len)
{
  if ((H5Tget_class(tobj) != H5T_ARRAY) || (H5Tget_array_ndims(tobj) != 1))
    return 0;
  hsize_t dim;
  H5Tget_array_dim2(tobj, &dim);
  *len = dim;
  hid_t te = H5Tget_super(tobj);
  CHECK_H5p(L, te, "Tget_super() failed in read(\"%s\")", path);
  int status = check_int_type(L, path, te);
  CHECK_H5p(L, H5Tclose(te), "Tclose() failed in read(\"%s\")", path);
  return status;
}

static int
check_vecreal_type(lua_State *L, const char *path, hid_t tobj, WriteSize *wsize, int *len)
{
  if ((H5Tget_class(tobj) != H5T_ARRAY) || (H5Tget_array_ndims(tobj) != 1))
    return 0;
  hsize_t dim;
  H5Tget_array_dim2(tobj, &dim);
  *len = dim;
  hid_t te = H5Tget_super(tobj);
  CHECK_H5p(L, te, "Tget_super() failed in read(\"%s\")", path);
  int status = check_real_type(L, path, te, wsize);
  CHECK_H5p(L, H5Tclose(te), "Tclose() failed in read(\"%s\")", path);
  return status;
}

static int
check_veccomplex_type(lua_State *L, const char *path, hid_t tobj, WriteSize *wsize, int *len)
{
  if ((H5Tget_class(tobj) != H5T_ARRAY) || (H5Tget_array_ndims(tobj) != 1))
    return 0;
  hsize_t dim;
  H5Tget_array_dim2(tobj, &dim);
  *len = dim;
  hid_t te = H5Tget_super(tobj);
  CHECK_H5p(L, te, "Tget_super() failed in read(\"%s\")", path);
  int status = check_complex_type(L, path, te, wsize);
  CHECK_H5p(L, H5Tclose(te), "Tclose() failed in read(\"%s\")", path);
  return status;
}

static int
check_matreal_type(lua_State *L, const char *path, hid_t tobj, WriteSize *wsize, int *l_len, int *r_len)
{
  if ((H5Tget_class(tobj) != H5T_ARRAY) || (H5Tget_array_ndims(tobj) != 2))
    return 0;
  hsize_t dims[2];
  H5Tget_array_dim2(tobj, dims);
  *l_len = dim[0];
  *r_len = dim[1];
  hid_t te = H5Tget_super(tobj);
  CHECK_H5p(L, te, "Tget_super() failed in read(\"%s\")", path);
  int status = check_real_type(L, path, te, wsize);
  CHECK_H5p(L, H5Tclose(te), "Tclose() failed in read(\"%s\")", path);
  return status;
}

static int
check_matcomplex_type(lua_State *L, const char *path, hid_t tobj, WriteSize *wsize, int *l_len, int *r_len)
{
  if ((H5Tget_class(tobj) != H5T_ARRAY) || (H5Tget_array_ndims(tobj) != 2))
    return 0;
  hsize_t dims[2];
  H5Tget_array_dim2(tobj, dims);
  *l_len = dim[0];
  *r_len = dim[1];
  hid_t te = H5Tget_super(tobj);
  CHECK_H5p(L, te, "Tget_super() failed in read(\"%s\")", path);
  int status = check_complex_type(L, path, te, wsize);
  CHECK_H5p(L, H5Tclose(te), "Tclose() failed in read(\"%s\")", path);
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

typedef int (*InUnpacker_H5)(lua_State *L, mLattice *S, mHdf5File *b, const char *path,
                             struct ropts_s *ropts, hid_t obj, hid_t tobj, SHA256_Sum *sum);

static int
r_string(lua_State *L, mLattice *S, mHdf5File *b, const char *path,
         struct ropts_s *ropts, hid_t obj, hid_t tobj, SHA256_Sum *sum)
{
  int len;
  if (!check_string_type(L, path, tobj, &len))
    return 0;
  hid_t memtype = get_string_type(L, b, len, 0);
  char *buffer = qmalloc(L, len + 1);
  herr_t status = H5Dread(obj, memtype, H5S_ALL, H5S_ALL, H5P_DEFAULT, buffer);
  CHECK_H5p(L, H5Tclose(memtype), "Tclose() failed in read(\"%s\")", path);
  if (status < 0) {
    qlua_free(buffer);
    return 0;
  }
  buffer[len] = 0;
  sha256_sum_string(sum, buffer, len);
  lua_pushstring(L, buffer);
  return 1;
}

/* all other r_*()'s */

static struct {
  KindOfH kind;
  int is_parallel;
  InUnpacker_H5 unpacker;
} qortable[] = {
  { kString,                  0,  r_string        },
#if 0 /* XXXXXXXXXXXXXXXXXXXXXXXXXX */
  { kReal,                    0,  r_real          },
  { kComplex,                 0,  r_complex       },
  { kMatrixReal,              0,  r_matreal       },
  { kMatrixComplex,           0,  r_matcomplex    },
  { kVectorInt,               0,  r_vecint        },
  { kVectorReal,              0,  r_vecreal       },
  { kVectorComplex,           0,  r_veccomplex    },
  { kColorVector,             0,  r_colvec        },
  { kColorMatrix,             0,  r_colmat        },
  { kDiracFermion,            0,  r_dirferm       },
  { kDiracPropagator,         0,  r_dirprop       },
  { kLatticeInt,              1,  r_latint        },
  { kLatticeReal,             1,  r_latreal       },
  { kLatticeComplex,          1,  r_latcomplex    },
  { kLatticeColorVector,      1,  r_latcolvec     },
  { kLatticeColorMatrix,      1,  r_latcolmat     },
  { kLatticeDiracFermion,     1,  r_latdirferm    },
  { kLatticeDiracPropagator,  1,  r_latdirprop    },
#endif /* XXXXXXXXXXXXXXXXXXXXXXXXXX */
  { kNoKind,                  0,  NULL            }
};

struct ropts_s {
  KindOfH kind;
  mLattice *S;
  int forced_p;
  int check_p;
};

static struct ropts_s
process_ropts(lua_State *L)
{
  struct ropts_s ropts;
  ropts.kind = kNoKind;
  ropts.S = NULL;
  ropts.forced_p = 0;
  ropts.check_p = 1;
  if (qlua_checkopt_table(L, 3)) {
    const char *check = qlua_tabkey_stringopt(L, 3, "sha256", "check");
    if (strcmp(check, "ignore") == 0)
      opts.check_p = 0;
    else if (strcmp(check, "check") != 0)
      luaL_error(L, "unknown check value \"%s\"", check);
    const char *forced = qlua_tabkey_stringopt(L, 3, "kind", NULL);
    if (forced) {
      ropts.forced_p = 1;
      ropts.kind = name2kind(forced);
    }
    if (qlua_tabpushopt_key(L, 3, "lattice")) {
      ropts.S = qlua_checkLattice(L, -1);
      lua_pop(L, 1);
    }
  }
  return ropts;
}

typedef int (*TheReader)(lua_State *L, mHdf5File *b, const char *path,
                         struct ropts_s *ropts, hid_t obj, hid_t dtype, hid_t dspace,
                         InUnpacker_H5 unpacker);

static int
read_seq(lua_State *L, mHdf5File *b, const char *path,
         struct ropts_s *ropts, hid_t obj, hid_t dtype, hid_t dspace,
         InUnpacker_H5 unpacker, SHA256_Sum *sum)
{
  if (!H5Sis_simple(dspace) || (H5get_simple_extent_ndims(dspace) != 0))
    return 0;
  return (*unpacker)(L, NULL, b, path, ropts, obj, dtype, sum);
}

static int
read_lat(lua_State *L, mHdf5File *b, const char *path,
         struct ropts_s *ropts, hid_t obj, hid_t dtype, hid_t dspace,
         InUnpacker_H5 unpacker, SHA256_Sum *sum)
{
  if (!H5Sis_simple(dspace) || (ropts->S == NULL))
    return 0;
  int rank = H5Sget_simple_extent_ndims(dspace);
  if (rank != ropts->S->rank)
    return 0;
  hsize_t *dims = qmalloc(L, rank * sizeof (hsize_t));
  H5Sget_simple_extent_dims(dspace, NULL, dims);
  int i;
  for (i = 0; i < rank; i++) {
    if (dims[i] != ropts->S->dim[i]) {
      qlua_free(L, dims);
      return 0;
    }
  }
  qlua_free(dims);
  int status = (*unpacker)(L, ropts->S, b, path, ropts, obj, dtype, sum);
  combine_checksums(sum, 1);
  return status;
}

/* If the reader managed to get the data, it pushes
 *   (data, "OK")           if sha256 matches
 *   (data, "missing")      if there is no sha256 in the file
 *   (data, "mismatched")   if sha256 do not agree
 * the last two cases can only happen if !ropts->check_p
 * if ropts->check_p, then the last two cases result in an error report.
 * The reader returns 1 if data was read, otherwise it returns 0 (e.g., no kind or space mismatch.)
 */
static int
try_reader(lua_State *L, mHdf5File *b, const char *path,
           struct ropts_s *ropts, hid_t obj, KindOfH kind)
{
  int i;
  for (i = 0; kortable[i].unpacker; i++) {
    if (kortable[i].kind == kind)
      break;
  }
  if (kortable[i].unpacker == NULL)
    return 0;
  hid_t dspace = H5Dget_space(obj);
  CHECK_H5p(L, dset, "Dget_space() failed in read(\"%s\")", path);
  hid_t dtype = H5Dget_type(obj);
  CHECK_H5p(L, dtype, "Dget_type() failed in read(\"%s\")", path);
  TheReader reader = kortable[i].is_parallel? read_lat: read_seq;
  SHA256_Sum r_sum;
  int status = (*reader)(L, b, path, ropts, obj, dtype, dspace, kortable[i].unpacker, &r_sum);
  CHECK_H5p(L, H5Tclose(dtype), "Tclose() failed in read(\"%s\")", path);
  CHECK_H5p(L, H5Tclose(dspace), "Sclose() failed in read(\"%s\")", path);
  if (!status)
    return 0;
  SHA256_Sum f_sum;
  if (!read_sha256(L, b, obj, &f_sum)) {
    if (ropts->check_p)
      luaL_error(L, "missing checksum in read(\"%s\")", path);
    lua_pushstring(L, "missing");
  } else {
    if (!sha256_cmp(&r_sum, &f_sum)) {
      if (ropts->check_p)
        luaL_error(L, "mismatched checksum in read(\"%s\")", path);
      lua_pushstring(L, "mismatched");
    } else {
      lua_pushstring(L, "OK");
    }
  }
  return 1;
}

static void
qhdf5_read(lua_State *L)
{
  mHdf5File *b = qlua_checkHdf5File(L, 1);
  const char *path = luaL_checkstring(L, 2);
  struct ropts_s ropts = process_ropts(L);

  check_file(L, b);
  qlua_Hdf5_enter(L);
  hid_t obj = H5Dopen(path[0] == '/'? b->file: b->cwd, path, H5P_DEFAULT);
  CHECK_H5p(L, obj, "no object for read(\"%s\")", path);
  int status;
  if (ropts->forced_p) {
    status = try_reader(L, b, path, &ropts, obj, ropts->kind);
  } else {
    KindOfH kind = get_h5_kind(L, b, obj);
    status = try_reader(L, b, path, &ropts, obj, kind);
  }
  CHECK_H5p(L, H5Dclose(obj), "Dclose() failed in read(\"%s\")", path);
  if (!status)
    luaL_error(L, "no suitable reader for read(\"%s\")", path);
  return 2;
}
