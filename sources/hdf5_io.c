#include "qlua.h"                                                    /* DEPS */
#include "qcomplex.h"                                                /* DEPS */
#include "qvector.h"                                                 /* DEPS */
#include "qmatrix.h"                                                 /* DEPS */
#include "hdf5_io.h"                                                 /* DEPS */
#include "sha256.h"                                                  /* DEPS */
#include "qlayout.h"                                                 /* DEPS */
#include "lattice.h"                                                 /* DEPS */
#include "latint.h"                                                  /* DEPS */
#include <string.h>
#include <complex.h>
#include <sys/stat.h>
#include "hdf5.h"
#include "qmp.h"

const char hdf5_io[] = "hdf5";
static const char mtnFile[] = "qcd.hdf5.mtFile";
const char csum_attr_name[] = ".sha256";
const char kind_attr_name[] = ".kind";
const char htn_complex_double[] = ".complexDouble";
const char htn_complex_float[] = ".complexFloat";

#define FAILED_H5CALL -1
#define CHECK_H5(L,expr,message) do { if (expr < 0) luaL_error(L, "HDF5 error %s", message); } while (0)

typedef struct {
  double re;
  double im;
} machine_complex_double;

typedef struct {
  float re;
  float im;
} machine_complex_float;

/* mappings to HDF5 types */
typedef struct htype_s {
  char *name;
  hid_t mtype;
  hid_t ftype;
  struct htype_s *next;
} HType;

struct mHdf5File_s {
  hid_t file;
  hid_t cwd;
  int master;
  int writer;
  HType *htype;
};

typedef struct QObjTable_s {
  QLUA_Type qtype;
  int is_par;
  int (*writer)(lua_State *L, mHdf5File *b, const char *path);
} QObjTable;

static QObjTable qotable[];

/* common helpers */
static void
local_combine_checksums(void *a, /* const */ void *b)
{
  SHA256_Sum *p = a;
  SHA256_Sum *q = b;
  int i;
  for (i = 0; i < sizeof (SHA256_Sum); i++)
    p->v[i] = p->v[i] ^ q->v[i];
}

static SHA256_Sum *
combine_checksums(SHA256_Sum *s, int is_writer)
{
  if (!is_writer)
    memset(s, 0, sizeof (SHA256_Sum));

  /* MPI: may require changes to use a communicator */
  QMP_binary_reduction(s, sizeof (SHA256_Sum), local_combine_checksums);
  return s;
}

static void
qlua_Hdf5_enter(lua_State *L)
{
  lua_gc(L, LUA_GCCOLLECT, 0);
}

static void
qlua_Hdf5_leave(void)
{
}

static void
free_types(lua_State *L, HType *ht)
{
  HType *next;

  for (; ht; ht = next) {
    next = ht->next;
    if (ht->mtype >= 0)
      H5Tclose(ht->mtype);
    ht->mtype = -1;
    if (ht->ftype >= 0)
      H5Tclose(ht->ftype);
    ht->ftype = -1;
    qlua_free(L, ht->name);
    qlua_free(L, ht);
  }
}

#if 0 /* XXXXXX complex hdf5 type constructors */
static HType *
lookup_htype(lua_State *L, mHdf5File *b, const char *name)
{
  HType *p;

  for (p = b->htype; p; p = p->next) {
    if (strcmp(p->name, name) == 0)
      break;
  }
  if (!p) {
    p = qlua_malloc(L, sizeof (HType));
    p->name = qlua_strdup(L, name);
    p->next = b->htype;
    p->ftype = -1;
    p->mtype = -1;
    b->htype = p;
  }
  return p;
}

static hid_t
construct_ftype_complex(lua_State *L, mHdf5File *b,
                        const char *name, size_t tsize, size_t r_offset, size_t i_offset)
{
  HType *p = lookup_htype(L, b, name);

  if (p->ftype < 0) {
    hid_t v = H5Tcreate(H5T_COMPOUND, tsize);
    CHECK_H5(L, v, "file complex type create failed");
    CHECK_H5(L, H5Tinsert(v, "r", r_offset, H5T_IEEE_F64BE), "complex.r insert failed");
    CHECK_H5(L, H5Tinsert(v, "i", i_offset, H5T_IEEE_F64BE), "complex.i insert failed");
    if (b->writer)
      CHECK_H5(L, H5Tcommit(b->file, name, v, H5P_DEFAULT,  H5P_DEFAULT,  H5P_DEFAULT), "complex commit failed");
    p->ftype = v;
  }
  return H5Tcopy(p->ftype);
}

static hid_t
construct_mtype_complex(lua_State *L, mHdf5File *b,
                        const char *name, size_t tsize, size_t r_offset, size_t i_offset)
{
  HType *p = lookup_htype(L, b, name);

  if (p->mtype < 0) {
    hid_t v = H5Tcreate(H5T_COMPOUND, tsize);
    CHECK_H5(L, v, "machine complex type create failed");
    CHECK_H5(L, H5Tinsert(v, "r", r_offset, H5T_NATIVE_DOUBLE), "complex.r insert failed");
    CHECK_H5(L, H5Tinsert(v, "i", i_offset, H5T_NATIVE_DOUBLE), "complex.i insert failed");
    p->mtype = v;
  }
  return H5Tcopy(p->mtype);
}

static hid_t
construct_ftype_complex_double(lua_State *L, mHdf5File *b)
{
  return construct_ftype_complex(L, b, htn_complex_double, 2 * 8, 0, 8);
}

static hid_t
construct_mtype_complex_double(lua_State *L, mHdf5File *b)
{
  return construct_mtype_complex(L, b, htn_complex_double,
                                 sizeof (machine_complex_double),
                                 HOFFSET(machine_complex_double, re),
                                 HOFFSET(machine_complex_double, im));
}
#endif /* XXXXXXXX */

/* writer */
static void
check_file(lua_State *L, mHdf5File *b)
{
  if (b->file < 0)
    luaL_error(L, "closed hdf5 file");
}

static void
check_writer(lua_State *L, mHdf5File *w)
{
  if (! w->writer)
    luaL_error(L, "expecting hdf5 writer");
  if (w->file < 0)
    luaL_error(L, "closed hdf5 writer");
}

static hid_t
make_hdf5_path(lua_State *L, hid_t file, hid_t cwd, char *dpath)
{
  char *dname = dpath[0] == '/' ? &dpath[1] : dpath;
  char *buf = dname;
  hid_t dhandle = dpath[0] == '/' ? H5Gopen2(file, "/", H5P_DEFAULT) : cwd;
  int isused = dpath[0] != '/';

  while (buf) {
    strsep(&buf, "/");
    hid_t dnext = H5Gopen2(dhandle, dname, H5P_DEFAULT);
    if (dnext < 0) {
      dnext = H5Gcreate2(dhandle, dname, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
      CHECK_H5(L, dnext, "mkpath failed");
    }
    if (isused == 0)
      H5Gclose(dhandle);
    isused = 0;
    dhandle = dnext;
    dname = buf;
  }
  return dhandle;
}

static mHdf5File *
qlua_newHdf5Writer(lua_State *L)
{
  mHdf5File *h = lua_newuserdata(L, sizeof (mHdf5File));

  h->master = (QDP_this_node == qlua_master_node);
  h->writer = 1;
  h->file = -1;
  h->cwd = -1;
  h->htype = NULL;

  luaL_getmetatable(L, mtnFile);
  lua_setmetatable(L, -2);

  return h;
}

static mHdf5File *
qlua_checkHdf5File(lua_State *L, int idx)
{
  void *v = luaL_checkudata(L, idx, mtnFile);

  luaL_argcheck(L, v != 0 , idx, "qcd.hdf5.File expected");

  return v;
}

mHdf5File *
qlua_checkHdf5Writer(lua_State *L, int idx)
{
  mHdf5File *v = qlua_checkHdf5File(L, idx);

  luaL_argcheck(L, v->writer , idx, "qcd.hdf5.Writer expected");

  return v;
}

static int
qhdf5_fmt(lua_State *L)
{
    mHdf5File *b = qlua_checkHdf5File(L, 1);
    char fmt[72];

    if (b->file >= 0)
      sprintf(fmt, "hdf5.%s(0x%llx)", b->writer? "Writer": "Reader", (long long int)b->file);
    else
      sprintf(fmt, "hdf5.%s(closed)", b->writer? "Writer": "Reader");

    lua_pushstring(L, fmt);
    return 1;
}

static int
do_close(lua_State *L, mHdf5File *b)
{
  int status = 0;

  qlua_Hdf5_enter(L);

  free_types(L, b->htype);
  b->htype = NULL;
  if (b->cwd >= 0)
    H5Gclose(b->cwd);
  b->cwd = -1;
  if (b->file >= 0)
    status = H5Fclose(b->file);
  b->file = -1;

  qlua_Hdf5_leave();
  return status;
}

static int
qhdf5_gc(lua_State *L)
{
  mHdf5File *b = qlua_checkHdf5File(L, 1);

  do_close(L, b);
  return 0;
}

static int
qhdf5_close(lua_State *L)
{
  mHdf5File *b = qlua_checkHdf5File(L, 1);

  check_file(L, b);
  // CHECK_H5(L, do_w_close(L, b), "writer close error");
  do_close(L, b);
  lua_pushnil(L);
  return 1;
}

static int
qhdf5_chpath(lua_State *L)
{
  mHdf5File *b = qlua_checkHdf5File(L, 1);
  const char *p = luaL_checkstring(L, 2);

  check_file(L, b);
  qlua_Hdf5_enter(L);

  char *dpath = qlua_strdup(L, p);
  char *dname = p[0] == '/' ? &dpath[1] : dpath;
  char *buf = dname;
  hid_t dhandle = p[0] == '/' ? H5Gopen2(b->file, "/", H5P_DEFAULT) : b->cwd;
  int isused = p[0] != '/';

  while (buf) {
    strsep(&buf, "/");
    if (dname[0] == 0)
      break;
    hid_t dnext = H5Gopen2(dhandle, dname, H5P_DEFAULT);
    CHECK_H5(L, dnext, "chpath failed");
    if (isused == 0)
      H5Gclose(dhandle);
    isused = 0;
    dhandle = dnext;
    dname = buf;
  }
  qlua_free(L, dpath);
  if (dhandle != b->cwd)
    H5Gclose(b->cwd);
  b->cwd = dhandle;

  qlua_Hdf5_leave();
  return 0;
}

static int
qhdf5_mkpath(lua_State *L)
{
  mHdf5File *b = qlua_checkHdf5Writer(L, 1);
  const char *p = luaL_checkstring(L, 2);

  check_writer(L, b);
  qlua_Hdf5_enter(L);

  char *dpath = qlua_strdup(L, p);
  hid_t dh = make_hdf5_path(L, b->file, b->cwd, dpath);
  qlua_free(L, dpath);
  if (dh != b->cwd)
    H5Gclose(dh);

  qlua_Hdf5_leave();
  return 0;
}


static int
qhdf5_kind(lua_State *L)
{
  // mHdf5File *b = qlua_checkHdf5Writer(L, 1);
  // const char *p = luaL_checkstring(L, 2);

  /* XXX kind */
  return 0;
}

static int
qhdf5_exists(lua_State *L)
{
  // mHdf5File *b = qlua_checkHdf5Writer(L, 1);
  // const char *p = luaL_checkstring(L, 2);

  /* XXX exists */
  return 0;
}

static void
write_attrs(lua_State *L, mHdf5File *b, hid_t dset, const SHA256_Sum *sum, const char *kind)
{
  hsize_t klen = strlen(kind);
  hid_t ktype = H5Tcopy(H5T_C_S1);
  CHECK_H5(L, H5Tset_size(ktype, klen), "set kind size");
  hid_t kds = H5Screate(H5S_SCALAR);
  hid_t kattr = H5Acreate2(dset, kind_attr_name, ktype, kds, H5P_DEFAULT, H5P_DEFAULT);
  CHECK_H5(L, H5Awrite(kattr, ktype, kind), "write kind attr");
  H5Sclose(kds);
  H5Aclose(kattr);
  H5Tclose(ktype);

  hsize_t slen = sizeof(sum->v);
  hid_t stype = H5Tcopy(H5T_STD_U8BE);
  hid_t sds = H5Screate_simple(1, &slen, NULL);
  hid_t sattr = H5Acreate2(dset, csum_attr_name, stype, sds, H5P_DEFAULT, H5P_DEFAULT);
  CHECK_H5(L, H5Awrite(sattr, stype, sum->v), "write sum attr");
  H5Sclose(sds);
  H5Aclose(sattr);
  H5Tclose(stype);

}

/* XXXX writers */
#if 0 /* XXXXXXXX scalar uncolored objects */
static void
w_scalar(lua_State *L, mHdf5File *b, const char *path, const char *kind,
         hid_t ftype, hid_t mtype, hid_t dataspace, const void *data, const SHA256_Sum *sum)
{
  char *dpath = qlua_strdup(L, path);
  char *ename = strrchr(dpath, '/');
  hid_t wdir = b->cwd;
  if (ename == NULL) {
    ename = dpath;
  } else {
    ename[0] = 0;
    ename = ename + 1;
    wdir = make_hdf5_path(L, b->file, b->cwd, dpath);
  }
  hid_t dataset = H5Dcreate2(wdir, ename, ftype, dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  hid_t plist = H5Pcreate(H5P_DATASET_XFER);
  if (b->master) {
    /* Only on the master */
    CHECK_H5(L, H5Dwrite(dataset, mtype, H5S_ALL, H5S_ALL, plist, data), "write string");
  }
  /* Everyone must write attributes, which must be identical or else. */
  write_attrs(L, b, dataset, sum, kind);

  CHECK_H5(L, H5Pclose(plist), "Pclose() plist");
  CHECK_H5(L, H5Dclose(dataset), "Dclose() dataset");
  CHECK_H5(L, H5Sclose(dataspace), "Sclose() dataspace");
  if (wdir != b->cwd)
    CHECK_H5(L, H5Gclose(wdir), "Gclose() write dir");
  CHECK_H5(L, H5Tclose(ftype), "Tclose() ftype");
  CHECK_H5(L, H5Tclose(mtype), "Tclose() mtype");
  qlua_free(L, dpath);
}

static int
w_string(lua_State *L, mHdf5File *b, const char *path)
{
  const char *str = luaL_checkstring(L, 3);
  size_t len = strlen(str);
  SHA256_Sum sum;

  sha256_sum_string(&sum, str, len);
  check_writer(L, b);

  qlua_Hdf5_enter(L);

  hid_t dataspace = H5Screate(H5S_SCALAR);
  hid_t ftype = H5Tcopy(H5T_C_S1);
  hid_t mtype = H5Tcopy(H5T_C_S1);
  CHECK_H5(L, H5Tset_size(ftype, len), "set ftype size");
  CHECK_H5(L, H5Tset_size(mtype, len), "set mtype size");
  w_scalar(L, b, path, "String", ftype, mtype, dataspace, str, global_combine_checksums(&sum, b->master));

  qlua_Hdf5_leave();

  return 0;
}

static int
w_real(lua_State *L, mHdf5File *b, const char *path)
{
  double val = luaL_checknumber(L, 3);
  SHA256_Context *ctx = sha256_create(L);
  SHA256_Sum sum;

  sha256_sum_add_doubles(ctx, &val, 1);
  sha256_sum(&sum, ctx);
  sha256_destroy(ctx);
  check_writer(L, b);

  qlua_Hdf5_enter(L);

  hid_t dataspace = H5Screate(H5S_SCALAR);
  hid_t ftype = H5Tcopy(H5T_IEEE_F64BE);
  hid_t mtype = H5Tcopy(H5T_NATIVE_DOUBLE);
  w_scalar(L, b, path, "Real", ftype, mtype, dataspace, &val, global_combine_checksums(&sum, b->master));

  qlua_Hdf5_leave();

  return 0;
}

static int
w_complex(lua_State *L, mHdf5File *b, const char *path)
{
  QLA_D_Complex *v = qlua_checkComplex(L, 3);
  SHA256_Context *ctx = sha256_create(L);
  machine_complex_double cv;
  SHA256_Sum sum;

  cv.re = QLA_real(*v);
  cv.im = QLA_imag(*v);
  sha256_sum_add_doubles(ctx, (double *)&cv, 2);
  sha256_sum(&sum, ctx);
  sha256_destroy(ctx);
  check_writer(L, b);

  qlua_Hdf5_enter(L);

  hid_t dataspace = H5Screate(H5S_SCALAR);
  hid_t ftype = construct_ftype_complex_double(L, b);
  hid_t mtype = construct_mtype_complex_double(L, b);
  w_scalar(L, b, path, "Complex", ftype, mtype, dataspace, &cv, global_combine_checksums(&sum, b->master));

  qlua_Hdf5_leave();

  return 0;
}

static int
w_vecint(lua_State *L, mHdf5File *b, const char *path)
{
  mVecInt *v = qlua_checkVecInt(L, 3);
  SHA256_Context *ctx = sha256_create(L);
  SHA256_Sum sum;
  hsize_t len;

  sha256_sum_add_ints(ctx, v->val, v->size);
  sha256_sum(&sum, ctx);
  sha256_destroy(ctx);
  check_writer(L, b);

  qlua_Hdf5_enter(L);
  len = v->size;
  hid_t dataspace = H5Screate_simple(1, &len, NULL);
  hid_t ftype = H5Tcopy(H5T_STD_I64BE);
  hid_t mtype = H5Tcopy(H5T_NATIVE_INT);
  w_scalar(L, b, path, "VectorInt", ftype, mtype, dataspace, v->val, global_combine_checksums(&sum, b->master));

  qlua_Hdf5_leave();

  return 0;
}

static int
w_vecreal(lua_State *L, mHdf5File *b, const char *path)
{
  mVecReal *v = qlua_checkVecReal(L, 3);
  SHA256_Context *ctx = sha256_create(L);
  SHA256_Sum sum;
  hsize_t len;

  sha256_sum_add_doubles(ctx, v->val, v->size);
  sha256_sum(&sum, ctx);
  sha256_destroy(ctx);
  check_writer(L, b);

  qlua_Hdf5_enter(L);
  len = v->size;
  hid_t dataspace = H5Screate_simple(1, &len, NULL);
  hid_t ftype = H5Tcopy(H5T_IEEE_F64BE);
  hid_t mtype = H5Tcopy(H5T_NATIVE_DOUBLE);
  w_scalar(L, b, path, "VectorReal", ftype, mtype, dataspace, v->val, global_combine_checksums(&sum, b->master));

  qlua_Hdf5_leave();

  return 0;
}

static int
w_veccomplex(lua_State *L, mHdf5File *b, const char *path)
{
  mVecComplex *v = qlua_checkVecComplex(L, 3);
  machine_complex_double *cv = qlua_malloc(L, v->size * sizeof (machine_complex_double));
  SHA256_Context *ctx = sha256_create(L);
  SHA256_Sum sum;
  hsize_t len;
  int i;

  for (i = 0; i < v->size; i++) {
    cv[i].re = QLA_real(v->val[i]);
    cv[i].im = QLA_imag(v->val[i]);
  }

  sha256_sum_add_doubles(ctx, (double *)cv, 2 * v->size);
  sha256_sum(&sum, ctx);
  sha256_destroy(ctx);
  check_writer(L, b);

  qlua_Hdf5_enter(L);
  len = v->size;
  hid_t dataspace = H5Screate_simple(1, &len, NULL);
  hid_t ftype = construct_ftype_complex_double(L, b);
  hid_t mtype = construct_mtype_complex_double(L, b);
  w_scalar(L, b, path, "VectorComplex", ftype, mtype, dataspace, v->val, global_combine_checksums(&sum, b->master));
  qlua_free(L, cv);

  qlua_Hdf5_leave();

  return 0;
}

static int
w_matreal(lua_State *L, mHdf5File *b, const char *path)
{
  mMatReal *v = qlua_checkMatReal(L, 3);
  SHA256_Context *ctx = sha256_create(L);
  double *cv = qlua_malloc(L, v->l_size * v->r_size * sizeof (double));
  SHA256_Sum sum;
  hsize_t len[2];
  int i, j;
  double *ptr;

  for (ptr = cv, i = 0; i < v->l_size; i++) {
    for (j = 0; j < v->r_size; j++, ptr++) {
      *ptr = gsl_matrix_get(v->m, i, j);
    }
  }
  sha256_sum_add_doubles(ctx, cv, v->l_size * v->r_size);
  sha256_sum(&sum, ctx);
  sha256_destroy(ctx);
  check_writer(L, b);

  qlua_Hdf5_enter(L);
  len[0] = v->l_size;
  len[1] = v->r_size;
  hid_t dataspace = H5Screate_simple(2, len, NULL);
  hid_t ftype = H5Tcopy(H5T_IEEE_F64BE);
  hid_t mtype = H5Tcopy(H5T_NATIVE_DOUBLE);
  w_scalar(L, b, path, "MatrixReal", ftype, mtype, dataspace, cv, global_combine_checksums(&sum, b->master));
  qlua_free(L, cv);
  qlua_Hdf5_leave();

  return 0;
}

static int
w_matcomplex(lua_State *L, mHdf5File *b, const char *path)
{
  mMatComplex *v = qlua_checkMatComplex(L, 3);
  machine_complex_double *cv = qlua_malloc(L, v->l_size * v->r_size * sizeof (machine_complex_double));
  SHA256_Context *ctx = sha256_create(L);
  SHA256_Sum sum;
  hsize_t len[2];
  int i, j;
  machine_complex_double *ptr;

  for (ptr = cv, i = 0; i < v->l_size; i++) {
    for (j = 0; j < v->r_size; j++, ptr++) {
      gsl_complex zz = gsl_matrix_complex_get(v->m, i, j);
      ptr->re = GSL_REAL(zz);
      ptr->im = GSL_IMAG(zz);
    }
  }

  sha256_sum_add_doubles(ctx, (double *)cv, 2 * v->l_size * v->r_size);
  sha256_sum(&sum, ctx);
  sha256_destroy(ctx);
  check_writer(L, b);

  qlua_Hdf5_enter(L);
  len[0] = v->l_size;
  len[1] = v->r_size;
  hid_t dataspace = H5Screate_simple(2, len, NULL);
  hid_t ftype = construct_ftype_complex_double(L, b);
  hid_t mtype = construct_mtype_complex_double(L, b);
  w_scalar(L, b, path, "MatrixComplex", ftype, mtype, dataspace, cv, global_combine_checksums(&sum, b->master));
  qlua_free(L, cv);

  qlua_Hdf5_leave();

  return 0;
}
#endif /* XXXXXXXXXXXXXXXXX scalar uncolored objects */

static int
w_latint(lua_State *L, mHdf5File *b, const char *path)
{
  mLatInt *m = qlua_checkLatInt(L, 3, NULL);
  // int has_opts = qlua_checkopt_table(L, 4);
  mLattice *S = qlua_ObjLattice(L, 3);
  int rank = S->rank;
  int *iptr = qlua_malloc(L, 3 * rank * sizeof (int));
  int *local_lo = iptr;
  int *local_hi = iptr + rank;
  int *local_x = iptr + 2 * rank;
  hsize_t *hptr = qlua_malloc(L, 5 * rank * sizeof (hsize_t));
  hsize_t *offset = hptr;
  hsize_t *stride = hptr + rank;
  hsize_t *count = hptr + 2 * rank;
  hsize_t *block = hptr + 3 * rank;
  hsize_t *latdim = hptr + 4 * rank;
  SHA256_Context *ctx = sha256_create(L);
  SHA256_Sum g_sum, l_sum;
  int volume, i, j, xl;
  int *data;

  qlua_sublattice(local_lo, local_hi, QDP_this_node, S);

  for (volume = 1, j = 0; j < rank; j++) {
    volume *= local_hi[j] - local_lo[j];
    offset[j] = local_lo[j];
    stride[j] = 1;
    count[j] = 1;
    block[j] = local_hi[j] - local_lo[j];
    latdim[j] = S->dim[j];
  }
  data = qlua_malloc(L, volume * sizeof (int));

  CALL_QDP(L);
  QLA_Int *locked = QDP_expose_I(m->ptr);
  sha256_sum_clear(&g_sum);
  for (i = 0; i < volume; i++) {
    for (xl = i, j = rank; j--;) {
      local_x[j] = xl % block[j] + local_lo[j];
      xl = xl / block[j];
    }
    QLUA_ASSERT(QDP_node_number_L(S->lat, local_x) == QDP_this_node);
    data[i] = QLA_elem_I(locked[QDP_index_L(S->lat, local_x)]);
    sha256_reset(ctx);
    sha256_sum_add_ints(ctx, &rank, 1);
    sha256_sum_add_ints(ctx, local_x, rank);
    sha256_sum_add_ints(ctx, &data[i], 1);
    sha256_sum(&l_sum, ctx);
    local_combine_checksums(&g_sum, &l_sum);
  }
  QDP_reset_I(m->ptr);

  combine_checksums(&g_sum, 1);

  qlua_Hdf5_enter(L);

  char *dpath = qlua_strdup(L, path);
  char *ename = strrchr(dpath, '/');
  hid_t wdir = b->cwd;
  if (ename == NULL) {
    ename = dpath;
  } else {
    ename[0] = 0;
    ename = ename + 1;
    wdir = make_hdf5_path(L, b->file, b->cwd, dpath);
  }

  hid_t ftype = H5Tcopy(H5T_STD_I64BE);
  hid_t mtype = H5Tcopy(H5T_NATIVE_INT);
  hid_t filespace = H5Screate_simple(rank, latdim, NULL);
  hid_t dataset = H5Dcreate2(wdir, ename, ftype, filespace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  CHECK_H5(L, H5Tclose(ftype), "Tclose() file type");
  CHECK_H5(L, H5Sclose(filespace), "Sclose() file space");
  if (wdir != b->cwd)
    CHECK_H5(L, H5Gclose(wdir), "Gclose() write dir");
  hid_t memspace = H5Screate_simple(rank, block, NULL);
  filespace =  H5Dget_space(dataset);
  H5Sselect_hyperslab(filespace, H5S_SELECT_SET, offset, stride, count, block);
  hid_t plist = H5Pcreate(H5P_DATASET_XFER);

  CHECK_H5(L, H5Dwrite(dataset, mtype, memspace, filespace, plist, data), "Dwrite() data");
  CHECK_H5(L, H5Pclose(plist), "Pclose() xfer plist");
  CHECK_H5(L, H5Sclose(filespace), "Sclose() file space");
  CHECK_H5(L, H5Sclose(memspace), "Sclose() mem space");
  CHECK_H5(L, H5Tclose(mtype), "Tclose() mem type");

  write_attrs(L, b, dataset, &g_sum, "LatticeInt");
  CHECK_H5(L, H5Dclose(dataset), "Dclose() data set");

  qlua_Hdf5_leave();

  sha256_destroy(ctx);
  qlua_free(L, data);
  qlua_free(L, iptr);
  qlua_free(L, hptr);

  return 0;
}

static int
qhdf5_write(lua_State *L)
{
  mHdf5File *b = qlua_checkHdf5Writer(L, 1);
  const char *p = luaL_checkstring(L, 2);
  QLUA_Type kind = qlua_qtype(L, 3);
  int count;
  int i;

  check_writer(L, b);
  qlua_Hdf5_enter(L);
  for (i = 0; qotable[i].qtype != qNoType; i++) {
    if (qotable[i].qtype == kind)
      break;
  }
  if (qotable[i].writer == NULL)
    luaL_error(L, "unwritable data");
  count = qotable[i].writer(L, b, p);

  qlua_Hdf5_leave();
  return count;
}

static int
q_hdf5_writer(lua_State *L)
{
  const char *name = luaL_checkstring(L, 1);
  mHdf5File *w = qlua_newHdf5Writer(L);
  MPI_Comm comm = MPI_COMM_WORLD;
  MPI_Info info = MPI_INFO_NULL;
  int status = 0;
  struct stat st;
  hid_t acc_tpl1;

  qlua_Hdf5_enter(L);

  acc_tpl1 = H5Pcreate(H5P_FILE_ACCESS);
  CHECK_H5(L, acc_tpl1, "Pcreate() failed");
  CHECK_H5(L, H5Pset_fapl_mpio(acc_tpl1, comm, info), "Pset_fapl_mpio() failed");

  /* possible fs race condition here */
  status = stat(name, &st);
 /* all nodes must agree on the file existance */
  QMP_sum_int(&status);
  if (status == 0) {
    /* open existing file */
    w->file = H5Fopen(name, H5F_ACC_RDWR, acc_tpl1);
  } else {
    /* create a new file. If it was created after stat() returned, truncate it */
    w->file = H5Fcreate(name, H5F_ACC_TRUNC, H5P_DEFAULT, acc_tpl1);
  }
  CHECK_H5(L, w->file, "qcd.hdf5.Writer failed");
  CHECK_H5(L, H5Pclose(acc_tpl1), "Pclose(template) failed");
  w->cwd = H5Gopen2(w->file, "/", H5P_DEFAULT);
  CHECK_H5(L, w->cwd, "Gopen2(\"/\") failed");

  qlua_Hdf5_leave();

  return 1;
}

/* setup */
static QObjTable qotable[] = {
#if 0 /* XXXXXX scalar uncolored objects */
  { qString,                 0,  w_string        },
  { qReal,                   0,  w_real          },
  { qComplex,                0,  w_complex       },
  { qVecInt,                 0,  w_vecint        },
  { qVecReal,                0,  w_vecreal       },
  { qVecComplex,             0,  w_veccomplex    },
  { qMatReal,                0,  w_matreal       },
  { qMatComplex,             0,  w_matcomplex    },
#endif /* XXXXX */
  { qLatInt,                 1,  w_latint        },
#if 0 /* XXXXX qlua object dispatch table */
  { qLatReal,                1,  w_latreal       },
  { qLatComplex,             1,  w_latcomplex    },
  { qSeqColVec2,             0,  w_colvec2       },
  { qSeqColMat2,             0,  w_colmat2       },
  { qSeqDirFerm2,            0,  w_dirferm2      },
  { qSeqDirProp2,            0,  w_dirprop2      },
  { qLatColVec2,             1,  w_latcolvec2    },
  { qLatColMat2,             1,  w_latcolmat2    },
  { qLatDirFerm2,            1,  w_latdirferm2   },
  { qLatDirProp2,            1,  w_latdirprop2   },
  { qSeqColVec3,             0,  w_colvec3       },
  { qSeqColMat3,             0,  w_colmat3       },
  { qSeqDirFerm3,            0,  w_dirferm3      },
  { qSeqDirProp3,            0,  w_dirprop3      },
  { qLatColVec3,             1,  w_latcolvec3    },
  { qLatColMat3,             1,  w_latcolmat3    },
  { qLatDirFerm3,            1,  w_latdirferm3   },
  { qLatDirProp3,            1,  w_latdirprop3   },
  { qSeqColVecN,             0,  w_colvecN       },
  { qSeqColMatN,             0,  w_colmatN       },
  { qSeqDirFermN,            0,  w_dirfermN      },
  { qSeqDirPropN,            0,  w_dirpropN      },
  { qLatColVecN,             1,  w_latcolvecN    },
  { qLatColMatN,             1,  w_latcolmatN    },
  { qLatDirFermN,            1,  w_latdirfermN   },
  { qLatDirPropN,            1,  w_latdirpropN   },
#endif /* XXXXX qlua object dispatch table */
  { qNoType,                 0,  NULL            }
};

/* metatables */
static const struct luaL_Reg mtFile[] = {
  { "__tostring",       qhdf5_fmt     },
  { "__gc",             qhdf5_gc      },
  { "close",            qhdf5_close   },
  { "write",            qhdf5_write   },
  { "chpath",           qhdf5_chpath  },
  { "mkpath",           qhdf5_mkpath  },
  { "kind",             qhdf5_kind    },
  { "exists",           qhdf5_exists  },
#if 0 /* XXXXXXXXXXXXXXXXXXXXXXXXXX */
  { "read",             qhdf5_real    },
  { "cwd",              qhdf5_cwd     },  /* XXX ? */
  { "list",             qhdf5_list    },
  { "dims",             qhdf5_dims    },
#endif /* XXXXXXXXXXXXXXXXXXXXXXXXXX */
  { NULL,               NULL          }
};

/* names and routines for qcd.hdf5 table */
static const struct luaL_Reg fHDF5io[] = {
  //{ "Reader",   q_hdf5_reader },
  { "Writer",   q_hdf5_writer  },
  { NULL,       NULL           }
};

int
init_hdf5_io(lua_State *L)
{
  H5open();
  H5Eset_auto2(H5E_DEFAULT, NULL, NULL);
  lua_getglobal(L, qcdlib);
  lua_newtable(L);
  luaL_register(L, NULL, fHDF5io);
  lua_setfield(L, -2, hdf5_io);
  lua_pop(L, 1);
  qlua_metatable(L, mtnFile, mtFile, qHdf5File);

  return 0;
}

int
fini_hdf5_io(lua_State *L)
{
  H5close();
  return 0;
}
