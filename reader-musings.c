#define CHECK_H5p(L, expr, message, path) ....

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
  H5Tget_array_dims2(tobj, &dim); /* XXX check error? */
  if (dim != sizeof (SHA256_Sum))
    return 0;
  hid_t base = H5Tget_super(tobj);
  CHECK_H5p(L, base, "Tget_super() failed on \"%s\"", path);
  int status = (H5Tget_class(base) == H5T_INTEGER);
  CHECK_H5p(L, H5Tclose(base), "Tclose(base) failed on \"%s\"", path);
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
    CHECK_H5p(L, et, "Tget_member_type() on \"%s\"", path);
    int status = check_real_type(L, et, &s[i]);
    CHECK_H5p(L, H5Tclose(et), "Tclose() on \"%s\"", path);
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
  H5Tget_array_dim2(tobj, &dim); /* XXX */
  *len = dim;
  hid_t te = H5Tget_super(tobj);
  CHECK_H5p(L, te, "Tget_super() on \"%s\"", path);
  int status = check_int_type(L, path, te);
  CHECK_H5p(L, H5Tclose(te), "Tclose() on \"%s\"", path);
  return status;
}

static int
check_vecreal_type(lua_State *L, const char *path, hid_t tobj, WriteSize *wsize, int *len)
{
  if ((H5Tget_class(tobj) != H5T_ARRAY) || (H5Tget_array_ndims(tobj) != 1))
    return 0;
  hsize_t dim;
  H5Tget_array_dim2(tobj, &dim); /* XXX */
  *len = dim;
  hid_t te = H5Tget_super(tobj);
  CHECK_H5p(L, te, "Tget_super() on \"%s\"", path);
  int status = check_real_type(L, path, te, wsize);
  CHECK_H5p(L, H5Tclose(te), "Tclose() on \"%s\"", path);
  return status;
}

static int
check_veccomplex_type(lua_State *L, const char *path, hid_t tobj, WriteSize *wsize, int *len)
{
  if ((H5Tget_class(tobj) != H5T_ARRAY) || (H5Tget_array_ndims(tobj) != 1))
    return 0;
  hsize_t dim;
  H5Tget_array_dim2(tobj, &dim); /* XXX status? */
  *len = dim;
  hid_t te = H5Tget_super(tobj);
  CHECK_H5p(L, te, "Tget_super() on \"%s\"", path);
  int status = check_complex_type(L, path, te, wsize);
  CHECK_H5p(L, H5Tclose(te), "Tclose() on \"%s\"", path);
  return status;
}

static int
check_matreal_type(lua_State *L, const char *path, hid_t tobj, WriteSize *wsize, int *l_len, int *r_len)
{
  if ((H5Tget_class(tobj) != H5T_ARRAY) || (H5Tget_array_ndims(tobj) != 2))
    return 0;
  hsize_t dims[2];
  H5Tget_array_dim2(tobj, dims); /* XXX */
  *l_len = dim[0];
  *r_len = dim[1];
  hid_t te = H5Tget_super(tobj);
  CHECK_H5p(L, te, "Tget_super() on \"%s\"", path);
  int status = check_real_type(L, path, te, wsize);
  CHECK_H5p(L, H5Tclose(te), "Tclose() on \"%s\"", path);
  return status;
}

static int
check_matcomplex_type(lua_State *L, const char *path, hid_t tobj, WriteSize *wsize, int *l_len, int *r_len)
{
  if ((H5Tget_class(tobj) != H5T_ARRAY) || (H5Tget_array_ndims(tobj) != 2))
    return 0;
  hsize_t dims[2];
  H5Tget_array_dim2(tobj, dims); /* XXX */
  *l_len = dim[0];
  *r_len = dim[1];
  hid_t te = H5Tget_super(tobj);
  CHECK_H5p(L, te, "Tget_super() on \"%s\"", path);
  int status = check_complex_type(L, path, te, wsize);
  CHECK_H5p(L, H5Tclose(te), "Tclose() on \"%s\"", path);
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
  H5Tget_array_dim2(tobj, &dim); /* XXX status? */
  if (dim != QDP_Ns)
    return 0;
  hid_t te = H5Tget_super(tobj);
  CHECK_H5p(L, te, "Tget_super() on \"%s\"", path);
  int status = check_veccomplex_type(L, path, te, wsize, Nc);
  CHECK_H5p(L, H5Tclose(te), "Tclose() on \"%s\"", path);
  return status;
}

static int
check_dirprop_type(lua_State *L, const char *path, hid_t tobj, WriteSize *wsize, int *Nc)
{
  if ((H5Tget_class(tobj) != H5T_ARRAY) || (H5Tget_array_ndims(tobj) != 2))
    return 0;
  hsize_t dims[2];
  H5Tget_array_dim2(tobj, dims); /* XXX status? */
  if ((dims[0] != QDP_Ns) || (dims[1] !+ QDP_Ns))
    return 0;
  hid_t te = H5Tget_super(tobj);
  CHECK_H5p(L, te, "Tget_super() on \"%s\"", path);
  int status = check_matcomplex_type(L, path, te, wsize, Nc);
  CHECK_H5p(L, H5Tclose(te), "Tclose() on \"%s\"", path);
  return status;
}

typedef int (*InUnpacker_H5)(lua_State *L, mLattice *S, mHdf5File *b, const char *path,
                             struct ropts_s *ropts, hid_t obj, hid_t tobj, SHA256_Sum *sum);

static int
r_string(lua_State *L, mLattice *S, mHdf5File *b, const char *path,
         struct ropts_s *ropts, hid_t obj, hid_t tobj, SHA256_Sum *sum)
{
  /* XXXX check that obj type is STRING, get length, return 0 on failure after freeing resources */
  /* XXX allocate membuf of length bytes for data */
  /* XXX read data, free resources on failure, return 0 */
  /* XXX push the string to L stack */
  /* XXX compute local checksum, free resources, return 1 */
}

/* all other r_*()'s */

static struct {
  KindOfH kind;
  int is_parallel;
  InUnpacker_H5 unpacker;
} qortable[] = {
  { kString,                  0,  r_string        },
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
  { kNoKind,                  0,  NULL            }
};

struct ropts_s {
  KindOfH kind;
  mLattice *S;
};

static struct ropts_s
process_ropts(lua_State *L)
{
  struct ropts_s ropts;
  memset(ropts, 0, sizeof (struct ropts_s));
  /* XXX ropts: S = ropts.lattice*/
  /* XXX ropts: kind = ropts.forcedKind if present, ropts.kind otherwise or kNoKind */
  /* XXX ropts: forced_p = ropts.forcedKind is present */
  /* XXX ropts: check_p = ropts.sha256 != "ignore" */
  return ropts;
}

static
read_common(...);
  /* XXX if checksum mismatch and ropts allow it, push (nil, value), return 1 */
  /* XXX otherwise if checksum mismatch, report error */
  /* XXX push (value, value), return 1 */

typedef int (*TheReader)(lua_State *L, mHdf5File *b, const char *path,
                         struct ropts_s *ropts, hid_t obj, hid_t dtype, hid_t dspace,
                         InUnpacker_H5 unpacker);

/* The reader closes dtype and dspace */
static int
read_seq(lua_State *L, mHdf5File *b, const char *path,
         struct ropts_s *ropts, hid_t obj, hid_t dtype, hid_t dspace,
         InUnpacker_H5 unpacker)
{
  /* XXX check that dataset is SCALAR */
  /* XXX try unpacker, return 0 if failed */
  /* XXX check that everyone got the same checksum */
  /* XXX free resources */
  /* XXX return read_common(...) */
}

static int
read_par(lua_State *L, mHdf5File *b, const char *path,
         struct ropts_s *ropts, hid_t obj, hid_t dtype, hid_t dspace,
         InUnpacker_H5 unpacker)
{
  /* XXX Check that data set rank and dimensions match ropts->S */
  /* XXX try unpacker, return 0 on failure */
  /* XXX compute global checksum */
  /* XXX free resources */
  /* XXX return read_common(...) */
}


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
  CHECK_H5p(L, dset, "Dget_space() in reader(\"%s\")", path);
  hid_t dtype = H5Dget_type(obj);
  CHECK_H5p(L, dtype, "Dget_type() in reader(\"%s\")", path);
  TheReader reader = kortable[i].is_parallel? read_par: read_seq;
  return (*reader)(L, b, path, ropts, obj, dtype, dspace, kortable[i].unpacker);
}

static void
qhdf5_read(lua_State *L)
{
  mHdf5File *b = qlua_checkHdf5File(L, 1);
  const char *path = luaL_checkstring(L, 2);
  struct ropts_s = process_ropts(L);

  check_file(L, b);
  qlua_Hdf5_enter(L);
  hid_t obj = H5Dopen(path[0] == '/'? b->file: b->cwd, path, H5P_DEFAULT);
  CHECK_H5p(L, obj, "hdf5 read(): no object at \"%s\"", path);

  /* XXX is ropts.forced_p: try_reader() and report error on failure, return on success */
  KindOfH kind = get_h5_kind(L, b, obj);
  /* XXX try natice kind, return on success */
  /* XXX if there is ropts.kind, try it, report error on failure */
  /* XXX free resources, return 2 */
}
