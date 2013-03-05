#include "qlua.h"                                                    /* DEPS */
#include "qvector.h"                                                 /* DEPS */
#include "qmatrix.h"                                                 /* DEPS */
#include "hdf5_io.h"                                                 /* DEPS */
#include "sha256.h"                                                  /* DEPS */
#include <string.h>
#include <complex.h>
#include <sys/stat.h>
#include "hdf5.h"
#include "qmp.h"

const char hdf5_io[] = "hdf5";
static const char mtnReader[] = "qcd.hdf5.mtReader";
static const char mtnWriter[] = "qcd.hdf5.mtRriter";

/* type names for H5 types */
enum {
  qhCharType,
  qhIntType,
  qhRealType,
  qhComplexType,
  qhTypeCount
};

/* sizes of disk images for qh*Type */
const int felem_size[] = { 1, 8, 8, 16, 0};
const char cs_attr_name[] = "sha256";

typedef struct {
  double re;
  double im;
} qlua_machine_complex;

struct mHdf5Writer_s {
  hid_t file;
  hid_t cwd;
  hid_t mt[qhTypeCount];
  hid_t ft[qhTypeCount];
  int master;
};

struct mHdf5Reader_s {
  hid_t file;
  hid_t cwd;
  hid_t mt[qhTypeCount];
  hid_t ft[qhTypeCount];
};

/* XXX all native H5 types should be replaced with fixed size */
static hid_t
get_qh_machine_type(lua_State *L, hid_t file, hid_t xf[], int tid)
{
  hid_t v = xf[tid];
  herr_t ec;

  if (v > 0)
    return v;
  if (xf[tid] > 0)
    return xf[tid];
  switch (tid) {
  case qhCharType: v = H5Tcopy(H5T_C_S1); break;
  case qhIntType: v = H5Tcopy(H5T_NATIVE_INT); break;
  case qhRealType: v = H5Tcopy(H5T_NATIVE_DOUBLE); break;
  case qhComplexType:
    v = H5Tcreate(H5T_COMPOUND, sizeof(qlua_machine_complex));
    if (v < 0) luaL_error(L, "HDF5: complex create failed");
    ec = H5Tinsert(v, "re", HOFFSET(qlua_machine_complex, re), H5T_NATIVE_DOUBLE);    
    if (ec < 0) luaL_error(L, "HDF5: complex.re insert failed");
    ec = H5Tinsert(v, "im", HOFFSET(qlua_machine_complex, im), H5T_NATIVE_DOUBLE);    
    if (ec < 0) luaL_error(L, "HDF5: complex.im insert failed");
    break;
  default:
    luaL_error(L, "HDF5: qh Type Name out of range");
  }
  if (v < 0)
    luaL_error(L, "HDF5: qh Type creation failed");
  xf[tid] = v;
  return v;
}

static hid_t
get_qh_file_type(lua_State *L, hid_t file, hid_t xf[], int tid)
{
  hid_t v = xf[tid];
  herr_t ec;

  if (v > 0)
    return v;
  switch (tid) {
  case qhCharType: v = H5Tcopy(H5T_C_S1); break;
  case qhIntType: v = H5Tcopy(H5T_STD_I64BE); break;
  case qhRealType: v = H5Tcopy(H5T_IEEE_F64BE); break;
  case qhComplexType:
    v = H5Tcreate(H5T_COMPOUND, sizeof(qlua_machine_complex));
    if (v < 0) luaL_error(L, "HDF5: complex create failed");
    ec = H5Tinsert(v, "re", HOFFSET(qlua_machine_complex, re), H5T_IEEE_F64BE);
    if (ec < 0) luaL_error(L, "HDF5: complex.re insert failed");
    ec = H5Tinsert(v, "im", HOFFSET(qlua_machine_complex, im), H5T_IEEE_F64BE);    
    if (ec < 0) luaL_error(L, "HDF5: complex.im insert failed");
    if (file > 0) {
      ec = H5Tcommit(file, ".complex", v, H5P_DEFAULT,  H5P_DEFAULT,  H5P_DEFAULT);
      if (ec < 0) luaL_error(L, "HDF5: complex commit failed");
    }
    break;
  default:
    luaL_error(L, "HDF5: qh Type Name out of range");
    break;
  }
  if (v < 0)
    luaL_error(L, "HDF5: qh Type creation failed");
  xf[tid] = v;
  return v;
}

/* helpers */
static void
check_reader(lua_State *L, mHdf5Reader *r)
{
  if (r->file < 0)
    luaL_error(L, "closed hdf5 reader");
}

static void
check_writer(lua_State *L, mHdf5Writer *r)
{
  if (QDP_this_node == qlua_master_node) {
    if (r->file < 0)
      luaL_error(L, "closed hdf5 writer");
  }
}

void
qlua_Hdf5_enter(lua_State *L)
{
  lua_gc(L, LUA_GCCOLLECT, 0);
}

void
qlua_Hdf5_leave(void)
{
}

static hid_t
make_hdf5_path(lua_State *L, hid_t file, hid_t cwd, char *dpath)
{
  char *dname = dpath[0] == '/' ? &dpath[1] : dpath;
  char *buf = dname;
  hid_t dhandle = dpath[0] == '/' ? H5Gopen2(file, "/", H5P_DEFAULT) : cwd;
  int isused = dpath[0] != '/';
  hid_t dnext;
    
  while (buf) {
    strsep(&buf, "/");
    dnext = H5Gopen2(dhandle, dname, H5P_DEFAULT);
    if (dnext < 0) {
      dnext = H5Gcreate2(dhandle, dname, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
      if (dnext < 0)
        luaL_error(L, "mkpath failed");
    }
    if (isused == 0)
      H5Gclose(dhandle);
    dhandle = dnext;
    isused = 0;
    dname = buf;
  }
  return dhandle;
}

static hid_t
change_hdf5_path(lua_State *L, hid_t file, hid_t cwd, const char *p)
{
  char *dpath = qlua_strdup(L, p);
  char *dname = p[0] == '/' ? &dpath[1] : dpath;
  char *buf = dname;
  hid_t dhandle = p[0] == '/' ? H5Gopen2(file, "/", H5P_DEFAULT) : cwd;
  int isused = p[0] != '/';
  hid_t dnext;
    
  while (buf) {
    strsep(&buf, "/");
    if (dname[0] == 0)
      break;
    dnext = H5Gopen2(dhandle, dname, H5P_DEFAULT);
    if (dnext < 0)
      luaL_error(L, "chpath failed");
    if (isused == 0)
      H5Gclose(dhandle);
    dhandle = dnext;
    isused = 0;
    dname = buf;
  }
  qlua_free(L, dpath);
  return dhandle;
}

/* allocation */
static mHdf5Reader *
qlua_newHdf5Reader(lua_State *L)
{
  mHdf5Reader *h = lua_newuserdata(L, sizeof (mHdf5Reader));
  int i;

  h->file = -1;
  h->cwd = -1;
  for (i = 0; i < qhTypeCount; i++) {
    h->mt[i] = -1;
    h->ft[i] = -1;
  }
  luaL_getmetatable(L, mtnReader);
  lua_setmetatable(L, -2);
  return h;
}

static mHdf5Writer *
qlua_newHdf5Writer(lua_State *L)
{
  mHdf5Writer *h = lua_newuserdata(L, sizeof (mHdf5Writer));
  int i;
  
  h->master = 0;
  h->file = -1;
  h->cwd = -1;
  for (i = 0; i < qhTypeCount; i++) {
    h->mt[i]= -1;
    h->ft[i] = -1;
  }
  luaL_getmetatable(L, mtnWriter);
  lua_setmetatable(L, -2);
  
  return h;
}

/* type checking */
mHdf5Reader *
qlua_checkHdf5Reader(lua_State *L, int idx)
{
  void *v = luaL_checkudata(L, idx, mtnReader);

  luaL_argcheck(L, v != 0, idx, "qcd.hdf5.Reader expected");
  return v;
}

mHdf5Writer *
qlua_checkHdf5Writer(lua_State *L, int idx)
{
  void *v = luaL_checkudata(L, idx, mtnWriter);
  
  luaL_argcheck(L, v != 0, idx, "qcd.hdf5.Writer expected");
  return v;
}

/* conversion to string */
static int
qhdf5_r_fmt(lua_State *L)
{
  mHdf5Reader *b = qlua_checkHdf5Reader(L, 1);
  char fmt[72];
    
  if (b->file > 0)
    sprintf(fmt, "hdf5.Reader(0x%llx)", (long long)b->file);
  else
    sprintf(fmt, "hdf5.Reader(closed)");

  lua_pushstring(L, fmt);
  return 1;
}

static int
qhdf5_w_fmt(lua_State *L)
{
    mHdf5Writer *b = qlua_checkHdf5Writer(L, 1);
    char fmt[72];

    if (b->master) {
      if (b->file > 0)
        sprintf(fmt, "hdf5.Writer(0x%llx)", (long long int)b->file);
      else
        sprintf(fmt, "hdf5.Writer(closed)");
    } else {
      sprintf(fmt, "hdf5.Writer(slave)");
    }
    
    lua_pushstring(L, fmt);
    return 1;
}

/* garbage collection */
static int
qhdf5_r_gc(lua_State *L)
{
  mHdf5Reader *b = qlua_checkHdf5Reader(L, 1);
  int i;

  qlua_Hdf5_enter(L);
  for (i = 0; i < qhTypeCount; i++) {
    if (b->mt[i] > 0)
      H5Tclose(b->mt[i]);
    b->mt[i] = -1;
    if (b->ft[i] > 0)
      H5Tclose(b->ft[i]);
    b->ft[i] = -1;
  }
  if (b->file > 0)
    H5Fclose(b->file);
  b->file = -1;
  qlua_Hdf5_leave();
  return 0;
}

static int
qhdf5_w_gc(lua_State *L)
{
  mHdf5Writer *b = qlua_checkHdf5Writer(L, 1);

  qlua_Hdf5_enter(L);
  if (b->master) {
    int i;

    for (i = 0; i < qhTypeCount; i++) {
      if (b->mt[i] > 0)
        H5Tclose(b->mt[i]);
      b->mt[i] = -1;
      if (b->ft[i] > 0)
        H5Tclose(b->ft[i]);
      b->ft[i] = -1;
    }
    b->cwd = -1;
    if (b->file > 0)
      H5Fclose(b->file);
    b->file = -1;
  }
  qlua_Hdf5_leave();
    
  return 0;
}

/* closing */
static int
qhdf5_r_close(lua_State *L)
{
  mHdf5Reader *b = qlua_checkHdf5Reader(L, 1);
  int i;

  check_reader(L, b);
  qlua_Hdf5_enter(L);
  for (i = 0; i < qhTypeCount; i++) {
    if (b->mt[i] > 0)
      H5Tclose(b->mt[i]);
    b->mt[i] = -1;
    if (b->ft[i] > 0)
      H5Tclose(b->ft[i]);
    b->ft[i] = -1;
  }
  H5Fclose(b->file);
  b->file = -1;
  qlua_Hdf5_leave();
  return 0;
}

static int
qhdf5_w_close(lua_State *L)
{
  mHdf5Writer *b = qlua_checkHdf5Writer(L, 1);
  int st = 0;

  check_writer(L, b);

  qlua_Hdf5_enter(L);
  if (b->master) {
    int i;

    for (i = 0; i < qhTypeCount; i++) {
      if (b->mt[i] > 0)
        H5Tclose(b->mt[i]);
      b->mt[i] = -1;
      if (b->ft[i] > 0)
        H5Tclose(b->ft[i]);
      b->ft[i] = -1;
    }
    if (b->cwd > 0)
      H5Gclose(b->cwd);
    b->cwd = -1;
    st = H5Fclose(b->file);
    b->file = -1;
    if (st < 0)
      luaL_error(L, "hdf5 writer close error");
  }
  qlua_Hdf5_leave();

  lua_pushnil(L);
  return 1;
}

static const char *
hdf5_vm(hid_t dh, const char *vname, const char *mname, const char *unknown)
{
  hid_t ds = H5Dget_space(dh);
  const char *result = unknown;
  switch (H5Sget_simple_extent_ndims(ds)) {
  case 1: result = vname; break;
  case 2: result = mname; break;
  default: break;
  }
  H5Sclose(ds);
  return result;
}

static int
qhdf5_kind(lua_State *L, hid_t file, hid_t cwd, const char *path)
{
  hid_t obj = H5Oopen(path[0] == '/'? file: cwd, path, H5P_DEFAULT);
  const char *kind = "unknown";

  switch (H5Iget_type(obj)) {
  case H5I_FILE: kind = "hdf5.file"; break;
  case H5I_GROUP: kind = "hdf5.group"; break;
  case H5I_DATATYPE: kind = "hdf5.type"; break;
  case H5I_DATASPACE: kind = "hdf5.space"; break;
  case H5I_ATTR: kind = "hdf5.attr"; break;
  case H5I_DATASET: {
    hid_t dh = H5Dopen2(path[0] == '/'? file: cwd, path, H5P_DEFAULT);
    hid_t dt = H5Dget_type(dh);
    switch (H5Tget_class(dt)) {
    case H5T_STRING: kind = hdf5_vm(dh, "string", kind, kind); break;
    case H5T_INTEGER: kind = hdf5_vm(dh, "vector.int", kind, kind); break;
    case H5T_FLOAT: kind = hdf5_vm(dh, "vector.real", "matrix.real", kind); break;
      /* all compounds are considered complex below */
    case H5T_COMPOUND: kind = hdf5_vm(dh, "vector.complex", "matrix.complex", kind); break;
    default: kind = "hdf5.data"; break;
    }
    H5Tclose(dt);
    H5Dclose(dh);
  } break;
  default: break;
  }
  H5Oclose(obj);
  lua_pushstring(L, kind);
  return 1;
}

static int
qhdf5_r_kind(lua_State *L)
{
  mHdf5Reader *b = qlua_checkHdf5Reader(L, 1);
  const char *p = luaL_checkstring(L, 2);
  
  check_reader(L, b);
  return qhdf5_kind(L, b->file, b->cwd, p);
}

static int
qhdf5_w_kind(lua_State *L)
{
  mHdf5Writer *b = qlua_checkHdf5Writer(L, 1);
  const char *p = luaL_checkstring(L, 2);
  
  check_writer(L, b);
  return qhdf5_kind(L, b->file, b->cwd, p);
}

static herr_t
qh_get_list(hid_t loc_id, const char *name, const H5L_info_t *linfo, void *opdata)
{
  lua_State *L = opdata;
  int k = luaL_checkinteger(L, -1);

  lua_pushstring(L, name);
  lua_rawseti(L, -3, k);
  lua_pop(L, 1);
  lua_pushinteger(L, k + 1);

  return 0;
}

static int
qhdf5_r_list(lua_State *L)
{
  mHdf5Reader *b = qlua_checkHdf5Reader(L, 1);
  const char *p = luaL_checkstring(L, 2);
  hid_t dh;
  herr_t ec;
  H5G_info_t gi;

  check_reader(L, b);
  qlua_Hdf5_enter(L);
  dh = H5Gopen(p[0] == '/'? b->file: b->cwd, p, H5P_DEFAULT);
  if (dh < 0) luaL_error(L, "qcd.hdf5.Reader:list() failed");
  ec = H5Gget_info(dh, &gi);
  if (ec < 0) luaL_error(L, "HDF5: Gget_info failed");
  lua_createtable(L, gi.nlinks, 0);
  lua_pushinteger(L, 1);
  ec = H5Literate(dh, H5_INDEX_NAME, H5_ITER_INC, NULL, qh_get_list, L);
  if (ec < 0) luaL_error(L, "HDF5: Literate failed");
  H5Gclose(dh);
  lua_pop(L, 1);
  return 1;
}

/* XXX rewrite it */
static void
verify_checksum(lua_State *L, hid_t xf[], hid_t vh, hsize_t flen, int qhtype)
{
  if (H5Aexists(vh, cs_attr_name) >= 0) {
    int size = flen * felem_size[qhtype];
    unsigned char *data = qlua_malloc(L, size);
    SHA256_Context *c = sha256_create(L);
    SHA256_Sum r_sum;
    SHA256_Sum f_sum;
    hid_t ft = get_qh_file_type(L, -1, xf, qhtype);
    hid_t ah = H5Aopen(vh, cs_attr_name, H5P_DEFAULT);
    hid_t at = H5Aget_type(ah);
    herr_t ec = H5Aread(ah, at, f_sum.v);

    if (ec < 0) {
      lua_pushstring(L, "ill-formed sha256");
      lua_insert(L, -2);
      lua_pushnil(L);
      lua_insert(L, -3);
      qlua_free(L, data);
      sha256_destroy(c);
      return;
    }

    H5Tclose(at);
    H5Aclose(ah);
    ec = H5Dread(vh, ft, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);
    if (ec < 0) luaL_error(L, "read error (sha256)");
    sha256_update(c, data, size);
    sha256_sum(&r_sum, c);
    sha256_destroy(c);
    qlua_free(L, data);

    if (sha256_cmp(&f_sum, &r_sum) == 0) {
      lua_pushstring(L, "OK");
      lua_pushnil(L);
    } else {
      lua_pushstring(L, "sha256 mismatch");
      lua_insert(L, -2);
      lua_pushnil(L);
      lua_insert(L, -3);
    }
  } else {
    lua_pushstring(L, "no sha256");
    lua_insert(L, -2);
    lua_pushnil(L);
    lua_insert(L, -3);
  }
}

static int
qhdf5_r_read(lua_State *L)
{
  mHdf5Reader *b = qlua_checkHdf5Reader(L, 1);
  const char *p = luaL_checkstring(L, 2);
  hid_t dh;
  hid_t dt;
  hid_t ds;
  herr_t ec;
  int qhtype = 0;
  hsize_t total_len = 0;
  void *data_ptr = NULL;

  check_reader(L, b);
  qlua_Hdf5_enter(L);
  dh = H5Dopen2(p[0] == '/'? b->file: b->cwd, p, H5P_DEFAULT);
  if (dh < 0) luaL_error(L, "qcd.hdf5.Reader:read() failed");
  dt = H5Dget_type(dh);
  ds = H5Dget_space(dh);
  switch (H5Tget_class(dt)) {
  case H5T_STRING: {
    hid_t mt = get_qh_machine_type(L, b->file, b->mt, qhCharType);
    hsize_t len;
    char *data;
    
    if (H5Sget_simple_extent_ndims(ds) != 1)
      luaL_error(L, "read(): string is not 1d");
    ec = H5Sget_simple_extent_dims(ds, &len, NULL);
    total_len = len;
    qhtype = qhCharType;
    if (ec < 0) luaL_error(L, "HDF5: read() get_dims failed");
    data_ptr = qlua_malloc(L, len + 1);
    data = (char *)data_ptr;
    ec = H5Dread(dh, mt, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);
    if (ec < 0) luaL_error(L, "read() failed");
    data[len] = 0;
    lua_pushstring(L, data);
  } break;
  case H5T_INTEGER: {
    hid_t mt = get_qh_machine_type(L, b->file, b->mt, qhIntType);
    hsize_t len;
    mVecInt *v = NULL;
    int *data;
    int i;
    
    if (H5Sget_simple_extent_ndims(ds) != 1)
      luaL_error(L, "read(): VecInt is not 1d");
    ec = H5Sget_simple_extent_dims(ds, &len, NULL);
    if (ec < 0) luaL_error(L, "HDF5: read() get_dims failed");
    total_len = len;
    qhtype = qhIntType;
    v = qlua_newVecInt(L, len);
    data_ptr = qlua_malloc(L, len * sizeof (int));
    data = (int *)data_ptr;
    ec = H5Dread(dh, mt, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);
    if (ec < 0) luaL_error(L, "read() failed");
    for (i = 0; i < len; i++)
      v->val[i] = data[i];
  } break;
  case H5T_FLOAT: {
    hid_t mt = get_qh_machine_type(L, b->file, b->mt, qhRealType);
    double *data;
    
    switch (H5Sget_simple_extent_ndims(ds)) {
    case 1: {
      hsize_t len;
      mVecReal *v = NULL;
      int i;

      ec = H5Sget_simple_extent_dims(ds, &len, NULL);
      if (ec < 0) luaL_error(L, "HDF5: read() get_dims failed");
      total_len = len;
      qhtype = qhRealType;
      v = qlua_newVecReal(L, len);
      data_ptr = qlua_malloc(L, len * sizeof (double));
      data = (double *)data_ptr;
      ec = H5Dread(dh, mt, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);
      if (ec < 0) luaL_error(L, "read() failed");
      for (i = 0; i < len; i++)
        v->val[i] = data[i];
    } break;
    case 2: {
      hsize_t len[2];
      mMatReal *v = NULL;
      int i, j;

      ec = H5Sget_simple_extent_dims(ds, &len[0], NULL);
      if (ec < 0) luaL_error(L, "HDF5: read() get_dims failed");
      total_len = len[0] * len[1];
      qhtype = qhRealType;
      v = qlua_newMatReal(L, len[0], len[1]);
      data_ptr = qlua_malloc(L, len[0] * len[1] * sizeof (double));
      data = (double *)data_ptr;
      ec = H5Dread(dh, mt, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);
      if (ec < 0) luaL_error(L, "read() failed");
      for (i = 0; i < len[0]; i++) {
        for (j = 0; j < len[1]; j++) {
          gsl_matrix_set(v->m, i, j, data[i * len[1] + j]);
        }
      }
    } break;
    default:
      luaL_error(L, "read(): float, but not VecReal nor MatReal");
    }
  } break;
  case H5T_COMPOUND: {
    hid_t mt = get_qh_machine_type(L, b->file, b->mt, qhComplexType);
    qlua_machine_complex *data;
    
    switch (H5Sget_simple_extent_ndims(ds)) {
    case 1: {
      hsize_t len;
      mVecComplex *v = NULL;
      int i;

      ec = H5Sget_simple_extent_dims(ds, &len, NULL);
      if (ec < 0) luaL_error(L, "HDF5: read() get_dims failed");
      total_len = len;
      qhtype = qhComplexType;
      v = qlua_newVecComplex(L, len);
      data_ptr = qlua_malloc(L, len * sizeof (qlua_machine_complex));
      data = (qlua_machine_complex *)data_ptr;
      ec = H5Dread(dh, mt, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);
      if (ec < 0) luaL_error(L, "read() failed");
      for (i = 0; i < len; i++) {
        QLA_real(v->val[i]) = data[i].re;
        QLA_imag(v->val[i]) = data[i].im;
      }
    } break;
    case 2: {
      hsize_t len[2];
      mMatComplex *v = NULL;
      int i, j;

      ec = H5Sget_simple_extent_dims(ds, &len[0], NULL);
      if (ec < 0) luaL_error(L, "HDF5: read() get_dims failed");
      total_len = len[0] * len[1];
      qhtype = qhComplexType;
      v = qlua_newMatComplex(L, len[0], len[1]);
      data_ptr = qlua_malloc(L, len[0] * len[1] * sizeof (qlua_machine_complex));
      data = (qlua_machine_complex *)data_ptr;
      ec = H5Dread(dh, mt, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);
      if (ec < 0) luaL_error(L, "read() failed");
      for (i = 0; i < len[0]; i++) {
        for (j = 0; j < len[1]; j++) {
          gsl_complex zz;
          GSL_REAL(zz) = data[i * len[1] + j].re;
          GSL_IMAG(zz) = data[i * len[1] + j].im;
          gsl_matrix_complex_set(v->m, i, j, zz);
        }
      }
    } break;
    default:
      luaL_error(L, "read(): float, but not VecReal nor MatReal");
    }
  } break;

  default:
    luaL_error(L, "read(): unsupported data type");
  }
  /* XXX compute sha256 on data_ptr, total_len, qhtype. Data may be reshuffled now */
  verify_checksum(L, b->ft, dh, total_len, qhtype);
  qlua_free(L, data_ptr);
  H5Sclose(ds);
  H5Tclose(dt);
  H5Gclose(dh);
  return 3; /* data, "OK", nil | nil, syndrome, data */
}

/* SHA256 checksum using HDF5 conversions */
/* XXX make it completely in-memory */
static void
create_checksum(lua_State *L, hid_t vh, hid_t ft, hsize_t flen, int qhtype)
{
  int size = flen * felem_size[qhtype];
  unsigned char *data = qlua_malloc(L, size);
  SHA256_Context *c = sha256_create(L);
  SHA256_Sum sum;
  hid_t tp = H5Tcopy(H5T_STD_U8BE);
  hsize_t alen = sizeof(sum.v);
  hid_t ds = H5Screate_simple(1, &alen, NULL);
  hid_t ah = H5Acreate2(vh, cs_attr_name, tp, ds, H5P_DEFAULT, H5P_DEFAULT);
  herr_t ec = H5Dread(vh, ft, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);

  if (ec < 0) luaL_error(L, "write error (readback)");
  sha256_update(c, data, size);
  sha256_sum(&sum, c);
  sha256_destroy(c);
  ec = H5Awrite(ah, tp, sum.v);
  if (ec < 0) luaL_error(L, "write error (sha256)");
  H5Aclose(ah);
  H5Sclose(ds);
  H5Tclose(tp);
}

static int
qhdf5_w_write(lua_State *L)
{
  mHdf5Writer *b = qlua_checkHdf5Writer(L, 1);
  const char *p = luaL_checkstring(L, 2);

  check_writer(L, b);

  qlua_Hdf5_enter(L);
  if (b->master) {
    char *dpath = qlua_strdup(L, p);
    hid_t wdir = b->cwd;
    hid_t dataset = -1;
    hid_t dtype_file = -1;
    hid_t dtype_machine = -1;
    hid_t dataspace = -1;
    void *data = NULL;
    char *ename = strrchr(dpath, '/');
    hsize_t full_len = 0;
    int qh_type = 0;
    herr_t status;

    if (ename == NULL) {
      ename = dpath;
    } else {
      ename[0] = 0;
      ename = ename + 1;
      wdir = make_hdf5_path(L, b->file, b->cwd, dpath);
    }
    switch (qlua_qtype(L, 3)) {
    case qString: {
      const char *str = luaL_checkstring(L, 3);
      hsize_t len = strlen(str);
      data = qlua_strdup(L, str);
      full_len = len;
      qh_type = qhCharType;
      dataspace = H5Screate_simple(1, &len, NULL);
    } break;
    case qVecInt: {
      mVecInt *v = qlua_checkVecInt(L, 3);
      hsize_t len = v->size;
      int *ptr; /* XXX should be fixed size */
      int i;
      ptr = qlua_malloc(L, len * sizeof (int));
      data = ptr;
      full_len = len;
      qh_type = qhIntType;
      dataspace = H5Screate_simple(1, &len, NULL);
      for (i = 0; i < len; i++)
        ptr[i] = v->val[i];
    } break;
    case qVecReal: {
      mVecReal *v = qlua_checkVecReal(L, 3);
      hsize_t len = v->size;
      double *ptr;
      int i;
      ptr = qlua_malloc(L, len * sizeof (double));
      data = ptr;
      full_len = len;
      qh_type = qhRealType;
      dataspace = H5Screate_simple(1, &len, NULL);
      for (i = 0; i < len; i++)
        ptr[i] = v->val[i];
    } break;
    case qVecComplex: {
      mVecComplex *v = qlua_checkVecComplex(L, 3);
      hsize_t len = v->size;
      qlua_machine_complex *ptr;
      int i;
      ptr = qlua_malloc(L, v->size * sizeof (qlua_machine_complex));
      data = ptr;
      full_len = len;
      qh_type = qhComplexType;
      dataspace = H5Screate_simple(1, &len, NULL);
      for (i = 0; i < v->size; i++) {
        ptr[i].re = QLA_real(v->val[i]);
        ptr[i].im = QLA_imag(v->val[i]);
      }
    } break;
    case qMatReal: {
      mMatReal *v = qlua_checkMatReal(L, 3);
      hsize_t len[2];
      double *ptr;
      int i, j;
      ptr = qlua_malloc(L, v->l_size * v->r_size * sizeof (double));
      len[0] = v->l_size;
      len[1] = v->r_size;
      data = ptr;
      full_len = len[0] * len[1];
      qh_type = qhRealType;
      dataspace = H5Screate_simple(2, len, NULL);
      for (i = 0; i < v->l_size; i++)
        for (j = 0; j < v->r_size; j++)
          *ptr++ = gsl_matrix_get(v->m, i, j);
    } break;
    case qMatComplex: {
      mMatComplex *v = qlua_checkMatComplex(L, 3);
      hsize_t len[2];
      qlua_machine_complex *ptr;
      int i, j;
      ptr = qlua_malloc(L, v->l_size * v->r_size * sizeof (qlua_machine_complex));
      len[0] = v->l_size;
      len[1] = v->r_size;
      data = ptr;
      full_len = len[0] * len[1];
      qh_type = qhComplexType;
      dataspace = H5Screate_simple(2, len, NULL);
      for (i = 0; i < v->l_size; i++) {
        for (j = 0; j < v->r_size; j++) {
          gsl_complex zz = gsl_matrix_complex_get(v->m, i, j);
          qlua_machine_complex zv;
          zv.re = GSL_REAL(zz);
          zv.im = GSL_IMAG(zz);
          *ptr++ = zv;
        }
      }
    } break;
    default:
      luaL_error(L, "Unsupported data type");
      break;
    }
    dtype_file = get_qh_file_type(L, b->file, b->ft, qh_type);
    dtype_machine = get_qh_machine_type(L, b->file, b->mt, qh_type);
    dataset = H5Dcreate2(wdir, ename, dtype_file, dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    status = H5Dwrite(dataset, dtype_machine, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);
    if (status < 0) luaL_error(L, "write error");

    /* XXX compute and write sha256 on data, full_len, qh_type. May reorder data */
    create_checksum(L, dataset, dtype_file, full_len, qh_type);
    H5Dclose(dataset);
    qlua_free(L, data);
    qlua_free(L, dpath);
    if (wdir != b->cwd)
      H5Gclose(wdir);
  }
  qlua_Hdf5_leave();
  return 0;
}

static int
qhdf5_r_chpath(lua_State *L)
{
  mHdf5Reader *b = qlua_checkHdf5Reader(L, 1);
  const char *p = luaL_checkstring(L, 2);
  hid_t dh;

  check_reader(L, b);

  qlua_Hdf5_enter(L);
  dh = change_hdf5_path(L, b->file, b->cwd, p);
  if (dh != b->cwd)
      H5Gclose(b->cwd);
  b->cwd = dh;
  qlua_Hdf5_leave();
  return 0;
}

static int
qhdf5_w_chpath(lua_State *L)
{
  mHdf5Writer *b = qlua_checkHdf5Writer(L, 1);
  const char *p = luaL_checkstring(L, 2);
  
  check_writer(L, b);
  
  qlua_Hdf5_enter(L);
  if (b->master) {
    hid_t dir = change_hdf5_path(L, b->file, b->cwd, p);
    if (dir != b->cwd)
      H5Gclose(b->cwd);
    b->cwd = dir;
  }
  qlua_Hdf5_leave();
  return 0;
}

static int
qhdf5_w_mkpath(lua_State *L)
{
  mHdf5Writer *b = qlua_checkHdf5Writer(L, 1);
  const char *p = luaL_checkstring(L, 2);

  check_writer(L, b);

  qlua_Hdf5_enter(L);
  if (b->master) {
    char *dpath = qlua_strdup(L, p);
    hid_t dh = make_hdf5_path(L, b->file, b->cwd, dpath);
    qlua_free(L, dpath);
    if (dh != b->cwd)
      H5Gclose(dh);
  }
  qlua_Hdf5_leave();
  return 0;
}

static int
q_hdf5_reader(lua_State *L)
{
  const char *name = luaL_checkstring(L, 1);
  mHdf5Reader *r = qlua_newHdf5Reader(L);
  int status;

  qlua_Hdf5_enter(L);
  r->file = H5Fopen(name, H5F_ACC_RDONLY, H5P_DEFAULT);
  status = (r->file < 0);
  if (status == 0) {
    r->cwd = H5Gopen2(r->file, "/", H5P_DEFAULT);
    status = (r->cwd < 0);
  }
  qlua_Hdf5_leave();
  QMP_sum_int(&status);
  
  if (status == 0)
    return 1;

  return luaL_error(L, "qcd.hdf5.Reader() failed");
}

static int
q_hdf5_writer(lua_State *L)
{
  const char *name = luaL_checkstring(L, 1);
  mHdf5Writer *w = qlua_newHdf5Writer(L);
  int status = 0;

  qlua_Hdf5_enter(L);
  if (QDP_this_node == qlua_master_node) {
    struct stat st;
    w->master = 1;
    /* here a race condition is possible ... */
    if (stat(name, &st) == 0) {
      /* ... hope the file has not been removed */
      w->file = H5Fopen(name, H5F_ACC_RDWR, H5P_DEFAULT);
    } else {
      /* ... if the file was created since the stat() above, truncate it! */
      w->file = H5Fcreate(name, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    }
    if (w->file < 0) {
      w->file = -1;
    } else {
      w->cwd = H5Gopen2(w->file, "/", H5P_DEFAULT);
      status = 1;
    }
  }
  qlua_Hdf5_leave();
  QMP_sum_int(&status);

  if (status)
        return 1;

  return luaL_error(L, "qcd.hdf5.Writer() failed");
}

/* metatables */
static const struct luaL_Reg mtReader[] = {
  { "__tostring",       qhdf5_r_fmt },
  { "__gc",             qhdf5_r_gc },
  { "close",            qhdf5_r_close },
  { "read",             qhdf5_r_read },
  { "chpath",           qhdf5_r_chpath },
  { "kind",             qhdf5_r_kind },
  { "list",             qhdf5_r_list },
  { NULL,               NULL}
};

static const struct luaL_Reg mtWriter[] = {
  { "__tostring",       qhdf5_w_fmt },
  { "__gc",             qhdf5_w_gc },
  { "close",            qhdf5_w_close },
  { "write",            qhdf5_w_write },
  { "chpath",           qhdf5_w_chpath },
  { "mkpath",           qhdf5_w_mkpath },
  { "kind",             qhdf5_w_kind },
  { NULL,               NULL}
};

/* names and routines for qcd.qdpc table */
static const struct luaL_Reg fHDF5io[] = {
  { "Reader",   q_hdf5_reader },
  { "Writer",   q_hdf5_writer },
  { NULL,       NULL }
};

int
init_hdf5_io(lua_State *L)
{
  H5Eset_auto2(H5E_DEFAULT, NULL, NULL);
  lua_getglobal(L, qcdlib);
  lua_newtable(L);
  luaL_register(L, NULL, fHDF5io);
  lua_setfield(L, -2, hdf5_io);
  lua_pop(L, 1);
  qlua_metatable(L, mtnReader, mtReader, qHdf5Reader);
  qlua_metatable(L, mtnWriter, mtWriter, qHdf5Writer);
  
  return 0;
}

int
fini_hdf5_io(lua_State *L)
{
  return 0;
}
