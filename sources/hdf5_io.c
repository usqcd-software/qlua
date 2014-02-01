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

static const char hdf5_io[] = "hdf5";
static const char mtnFile[] = "qcd.hdf5.mtFile";
static const char csum_attr_name[] = ".sha256";
static const char kind_attr_name[] = ".kind";
static const char htn_complex_double[]     = ".complexDouble";
static const char htn_complex_float[]      = ".complexFloat";
static const char htn_vector_int_fmt[]     = ".vectorInt%d";
static const char htn_vector_real_fmt[]    = ".vectorReal%s%d";
static const char htn_vector_complex_fmt[] = ".vectorComplex%s%d";

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

struct laddr_s {
  int rank;  /* lattice rank */
  int volume; /* local volume */
  int *low;  /* local low site */
  int *high; /* local high site */
};

typedef enum {
  WS_Double,
  WS_Float
} WriteSize;

struct wopts_s {
  WriteSize wsize;
};

typedef void (*OutPacker_H5)(lua_State *L, mHdf5File *b, mLattice *S,
                             struct wopts_s *opts, struct laddr_s *laddr,
                             SHA256_Sum *sum, void **data, hid_t *filetype, hid_t *memtype,
                             const char **kind);

typedef struct QObjTable_s {
  QLUA_Type qtype;
  int is_parallel;
  OutPacker_H5 repack;
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

static const char *
wsize_name(lua_State *L, WriteSize wsize)
{
  switch (wsize) {
  case WS_Double: return "Double";
  case WS_Float: return "Float";
  default:
    luaL_error(L, "Unknown write size %d", (int)wsize);
    break;
  }
  return NULL;
}

static hid_t
wsize_real_mtype(lua_State *L, WriteSize wsize)
{
  hid_t t = -1;

  switch (wsize) {
  case WS_Double: t = H5T_NATIVE_DOUBLE; break;
  case WS_Float: t = H5T_NATIVE_FLOAT; break;
  default: break;
  }
  QLUA_ASSERT(t != -1);
  t = H5Tcopy(t);
  CHECK_H5(L, t, "Tcopy(real mem type)");
  return t;
}

static hid_t
wsize_real_ftype(lua_State *L, WriteSize wsize)
{
  hid_t t = -1;

  switch (wsize) {
  case WS_Double: t = H5T_IEEE_F64BE; break;
  case WS_Float: t = H5T_IEEE_F32BE; break;
  default: break;
  }
  QLUA_ASSERT(t != -1);
  t = H5Tcopy(t);
  CHECK_H5(L, t, "Tcopy(real file type)");
  return t;
}

static hid_t
construct_ftype_complex(lua_State *L, mHdf5File *b,
                        const char *name, hid_t ftype, size_t tsize, size_t r_offset, size_t i_offset)
{
  HType *p = lookup_htype(L, b, name);

  if (p->ftype < 0) {
    hid_t v = H5Tcreate(H5T_COMPOUND, tsize);
    CHECK_H5(L, v, "file complex type create failed");
    CHECK_H5(L, H5Tinsert(v, "r", r_offset, ftype), "complex.r insert failed");
    CHECK_H5(L, H5Tinsert(v, "i", i_offset, ftype), "complex.i insert failed");
    if (b->writer)
      CHECK_H5(L, H5Tcommit(b->file, name, v, H5P_DEFAULT,  H5P_DEFAULT,  H5P_DEFAULT), "complex commit failed");
    p->ftype = v;
  }
  return H5Tcopy(p->ftype);
}

static hid_t
construct_mtype_complex(lua_State *L, mHdf5File *b,
                        const char *name, hid_t ftype, size_t tsize, size_t r_offset, size_t i_offset)
{
  HType *p = lookup_htype(L, b, name);

  if (p->mtype < 0) {
    hid_t v = H5Tcreate(H5T_COMPOUND, tsize);
    CHECK_H5(L, v, "machine complex type create failed");
    CHECK_H5(L, H5Tinsert(v, "r", r_offset, ftype), "complex.r insert failed");
    CHECK_H5(L, H5Tinsert(v, "i", i_offset, ftype), "complex.i insert failed");
    p->mtype = v;
  }
  return H5Tcopy(p->mtype);
}

static hid_t
construct_ftype_complex_double(lua_State *L, mHdf5File *b)
{
  return construct_ftype_complex(L, b, htn_complex_double, H5T_IEEE_F64BE, 2 * 8, 0, 8);
}

static hid_t
construct_mtype_complex_double(lua_State *L, mHdf5File *b)
{
  return construct_mtype_complex(L, b, htn_complex_double,
                                 H5T_NATIVE_DOUBLE,
                                 sizeof (machine_complex_double),
                                 HOFFSET(machine_complex_double, re),
                                 HOFFSET(machine_complex_double, im));
}

static hid_t
construct_ftype_complex_float(lua_State *L, mHdf5File *b)
{
  return construct_ftype_complex(L, b, htn_complex_float, H5T_IEEE_F32BE, 2 * 4, 0, 4);
}

static hid_t
construct_mtype_complex_float(lua_State *L, mHdf5File *b)
{
  return construct_mtype_complex(L, b, htn_complex_float,
                                 H5T_NATIVE_FLOAT,
                                 sizeof (machine_complex_float),
                                 HOFFSET(machine_complex_float, re),
                                 HOFFSET(machine_complex_float, im));
}

static hid_t
construct_ftype_vecint(lua_State *L, mHdf5File *b, int size)
{
  char fmt[128];
  sprintf(fmt, htn_vector_int_fmt, size);
  HType *p = lookup_htype(L, b, fmt);
  if (p->ftype < 0) {
    hsize_t hsize = size;
    hid_t v = H5Tarray_create2(H5T_STD_I32BE, 1, &hsize);
    CHECK_H5(L, v, "file vecint type create failed");
    if (b->writer)
      CHECK_H5(L, H5Tcommit(b->file, fmt, v, H5P_DEFAULT,  H5P_DEFAULT,  H5P_DEFAULT), "vecint commit failed");
    p->ftype = v;
  }
  return H5Tcopy(p->ftype);
}

static hid_t
construct_mtype_vecint(lua_State *L, mHdf5File *b, int size)
{
  char fmt[128];
  sprintf(fmt, htn_vector_int_fmt, size);
  HType *p = lookup_htype(L, b, fmt);
  if (p->mtype < 0) {
    hsize_t hsize = size;
    hid_t v = H5Tarray_create2(H5T_NATIVE_INT, 1, &hsize);
    CHECK_H5(L, v, "machine vecint type create failed");
    p->mtype = v;
  }
  return H5Tcopy(p->mtype);
}

static hid_t
construct_ftype_vecreal(lua_State *L, mHdf5File *b, int size, WriteSize wsize)
{
#if 0 /* XXXXXXXXXXX */
  char fmt[128];
  sprintf(fmt, htn_vector_int_fmt, size);
  HType *p = lookup_htype(L, b, fmt);
  if (p->ftype < 0) {
    hsize_t hsize = size;
    hid_t v = H5Tarray_create2(H5T_STD_I32BE, 1, &hsize);
    CHECK_H5(L, v, "file vecint type create failed");
    if (b->writer)
      CHECK_H5(L, H5Tcommit(b->file, fmt, v, H5P_DEFAULT,  H5P_DEFAULT,  H5P_DEFAULT), "vecint commit failed");
    p->ftype = v;
  }
  return H5Tcopy(p->ftype);
#endif /* XXXXXXXXXXX */
  return -1;
}

static hid_t
construct_mtype_vecreal(lua_State *L, mHdf5File *b, int size, WriteSize wsize)
{
#if 0 /* XXXX */
  char fmt[128];
  sprintf(fmt, htn_vector_int_fmt, size);
  HType *p = lookup_htype(L, b, fmt);
  if (p->mtype < 0) {
    hsize_t hsize = size;
    hid_t v = H5Tarray_create2(H5T_NATIVE_INT, 1, &hsize);
    CHECK_H5(L, v, "machine vecint type create failed");
    p->mtype = v;
  }
  return H5Tcopy(p->mtype);
#endif /* XXXX */
  return -1;
}

static hid_t
construct_ftype_veccomplex(lua_State *L, mHdf5File *b, int size, WriteSize wsize)
{
#if 0 /* XXXXXXXXXXX */
  char fmt[128];
  sprintf(fmt, htn_vector_int_fmt, size);
  HType *p = lookup_htype(L, b, fmt);
  if (p->ftype < 0) {
    hsize_t hsize = size;
    hid_t v = H5Tarray_create2(H5T_STD_I32BE, 1, &hsize);
    CHECK_H5(L, v, "file vecint type create failed");
    if (b->writer)
      CHECK_H5(L, H5Tcommit(b->file, fmt, v, H5P_DEFAULT,  H5P_DEFAULT,  H5P_DEFAULT), "vecint commit failed");
    p->ftype = v;
  }
  return H5Tcopy(p->ftype);
#endif /* XXXXXXXXXXX */
  return -1;
}

static hid_t
construct_mtype_veccomplex(lua_State *L, mHdf5File *b, int size, WriteSize wsize)
{
#if 0 /* XXXX */
  char fmt[128];
  sprintf(fmt, htn_vector_int_fmt, size);
  HType *p = lookup_htype(L, b, fmt);
  if (p->mtype < 0) {
    hsize_t hsize = size;
    hid_t v = H5Tarray_create2(H5T_NATIVE_INT, 1, &hsize);
    CHECK_H5(L, v, "machine vecint type create failed");
    p->mtype = v;
  }
  return H5Tcopy(p->mtype);
#endif /* XXXX */
  return -1;
}

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

  while (buf && buf[0]) {
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

static void
write_attrs(lua_State *L, mHdf5File *b, hid_t dset, const SHA256_Sum *sum, const char *kind)
{
  hsize_t klen = strlen(kind);
  hid_t ktype = H5Tcopy(H5T_C_S1);
  CHECK_H5(L, ktype, "Tcopy(C_S1) in write_attrs()");
  CHECK_H5(L, H5Tset_size(ktype, klen), "set kind size");
  hid_t kds = H5Screate(H5S_SCALAR);
  CHECK_H5(L, kds, "Screate(SCALAR) in write_attrs()");
  hid_t kattr = H5Acreate2(dset, kind_attr_name, ktype, kds, H5P_DEFAULT, H5P_DEFAULT);
  CHECK_H5(L, kattr, "Acreate(kind) in write_attrs()");
  CHECK_H5(L, H5Awrite(kattr, ktype, kind), "Awrite(kind) in write_attrs()");
  CHECK_H5(L, H5Sclose(kds), "Sclose(kind) in write_attrs()");
  CHECK_H5(L, H5Aclose(kattr), "Aclose(kind) in write_attrs()");
  CHECK_H5(L, H5Tclose(ktype), "Tclose(kind) in write_attrs()");

  hsize_t slen = sizeof(sum->v);
  hid_t stype = H5Tcopy(H5T_STD_U8BE);
  CHECK_H5(L, stype, "Tcopy(U8BE) in write_attrs()");
  hid_t sds = H5Screate_simple(1, &slen, NULL);
  CHECK_H5(L, sds, "Screate(sum) in write_attrs()");
  hid_t sattr = H5Acreate2(dset, csum_attr_name, stype, sds, H5P_DEFAULT, H5P_DEFAULT);
  CHECK_H5(L, sattr, "Acreate(sum) in write_attrs()");
  CHECK_H5(L, H5Awrite(sattr, stype, sum->v), "Awrite(sum) in write_attrs()");
  CHECK_H5(L, H5Sclose(sds), "Sclose(sum) in write_attrs()");
  CHECK_H5(L, H5Aclose(sattr), "Aclose(sum) in write_attrs()");
  CHECK_H5(L, H5Tclose(stype), "Tclose(sum) in write_attrs()");

}

/* XXXX writers */
#if 0 /* XXXXXXXX scalar uncolored objects */

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

static void
qdp2hdf5_addr(int *local_x, int xl, struct laddr_s *laddr)
{
  int j;
  int rank = laddr->rank;

  for (j = rank; j--;) {
    int block_j = laddr->high[j] - laddr->low[j];
    local_x[j] = xl % block_j + laddr->low[j];
    xl = xl / block_j;
  }
}

static void
process_wopts(lua_State *L, struct wopts_s *opts) /* XXX */
{
  opts->wsize = WS_Double;
  if (qlua_checkopt_table(L, 4)) {
    const char *prec = qlua_tabkey_stringopt(L, 4, "precision", "double");
    if (strcmp(prec, "double") == 0)
      opts->wsize = WS_Double;
    else if (strcmp(prec, "float") == 0)
      opts->wsize = WS_Float;
    else
      luaL_error(L, "Unknown precision value \"%s\"", prec);
  }
}

static void
w_string(lua_State *L, mHdf5File *b, mLattice *S,
         struct wopts_s *opts, struct laddr_s *laddr,
         SHA256_Sum *sum, void **data, hid_t *filetype, hid_t *memtype,
         const char **kind)
{
  const char *str = luaL_checkstring(L, 3);
  size_t len = strlen(str) + 1;
  sha256_sum_string(sum, str, len);
  *kind = "String";
  *data = qlua_strdup(L, str);
  *filetype = H5Tcopy(H5T_C_S1);
  CHECK_H5(L, *filetype, "Tcopy(filetype) in w_string()");
  CHECK_H5(L, H5Tset_size(*filetype, len), "Sset_size(filetype) in w_string()");
  *memtype = H5Tcopy(H5T_C_S1);
  CHECK_H5(L, *memtype, "Tcopy(memtype) in w_string()");
  CHECK_H5(L, H5Tset_size(*memtype, len), "Sset_size(memtype) in w_string()");
}

static void
w_real(lua_State *L, mHdf5File *b, mLattice *S,
       struct wopts_s *opts, struct laddr_s *laddr,
       SHA256_Sum *sum, void **data, hid_t *filetype, hid_t *memtype,
       const char **kind)
{
  double val = luaL_checknumber(L, 3);
  SHA256_Context *ctx = sha256_create(L);

  switch (opts->wsize) {
  case WS_Double:
    *data = qlua_malloc(L, sizeof (double));
    *(double *)(*data) = val;
    sha256_sum_add_doubles(ctx, *data, 1);
    *filetype = H5Tcopy(H5T_IEEE_F64BE);
    CHECK_H5(L, *filetype, "Tcopy(filetype) in w_real()");
    *memtype = H5Tcopy(H5T_NATIVE_DOUBLE);
    CHECK_H5(L, *memtype, "Tcopy(memtype) in w_real()");
    break;
  case WS_Float:
    *data = qlua_malloc(L, sizeof (float));
    *(float *)(*data) = val;
    sha256_sum_add_floats(ctx, *data, 1);
    *filetype = H5Tcopy(H5T_IEEE_F32BE);
    CHECK_H5(L, *filetype, "Tcopy(filetype) in w_real()");
    *memtype = H5Tcopy(H5T_NATIVE_FLOAT);
    CHECK_H5(L, *memtype, "Tcopy(memtype) in w_real()");
    break;
  default:
    luaL_error(L, "Unknown precision in w_real()");
  }
  sha256_sum(sum, ctx);
  sha256_destroy(ctx);
  *kind = "Real";
}

static void
w_complex(lua_State *L, mHdf5File *b, mLattice *S,
          struct wopts_s *opts, struct laddr_s *laddr,
          SHA256_Sum *sum, void **data, hid_t *filetype, hid_t *memtype,
          const char **kind)
{
  QLA_D_Complex *v = qlua_checkComplex(L, 3);
  SHA256_Context *ctx = sha256_create(L);

  switch (opts->wsize) {
  case WS_Double: {
    machine_complex_double *cv = qlua_malloc(L, sizeof (machine_complex_double));
    cv->re = QLA_real(*v);
    cv->im = QLA_imag(*v);
    *data = cv;
    sha256_sum_add_doubles(ctx, *data, 2);
    *filetype = construct_ftype_complex_double(L, b);
    *memtype = construct_mtype_complex_double(L, b);
  } break;
  case WS_Float: {
    machine_complex_float *cv = qlua_malloc(L, sizeof (machine_complex_float));
    cv->re = QLA_real(*v);
    cv->im = QLA_imag(*v);
    *data = cv;
    sha256_sum_add_floats(ctx, *data, 2);
    *filetype = construct_ftype_complex_float(L, b);
    *memtype = construct_mtype_complex_float(L, b);
  } break;
  default:
    luaL_error(L, "Unknown precision in w_complex()");
  }
  sha256_sum(sum, ctx);
  sha256_destroy(ctx);
  *kind = "Complex";
}

static void
w_vecint(lua_State *L, mHdf5File *b, mLattice *S,
         struct wopts_s *opts, struct laddr_s *laddr,
         SHA256_Sum *sum, void **data, hid_t *filetype, hid_t *memtype,
         const char **kind)
{
  mVecInt *v = qlua_checkVecInt(L, 3);
  SHA256_Context *ctx = sha256_create(L);
  int *vec = qlua_malloc(L, v->size * sizeof (int));
  int i;

  for (i = 0; i < v->size; i++)
    vec[i] = v->val[i];
  sha256_sum_add_ints(ctx, v->val, v->size);
  sha256_sum(sum, ctx);
  sha256_destroy(ctx);
  *data = vec;
  *filetype = construct_ftype_vecint(L, b, v->size);
  *memtype = construct_mtype_vecint(L, b, v->size);
  *kind = "VectorInt";
}

static void
w_vecreal(lua_State *L, mHdf5File *b, mLattice *S,
          struct wopts_s *opts, struct laddr_s *laddr,
          SHA256_Sum *sum, void **data, hid_t *filetype, hid_t *memtype,
          const char **kind)
{
  mVecReal *v = qlua_checkVecReal(L, 3);
  SHA256_Context *ctx = sha256_create(L);
  int i;

  switch (opts->wsize) {
  case WS_Double: {
    double *vec= qlua_malloc(L, v->size * sizeof (double));
    for (i = 0; i < v->size; i++)
      vec[i] = v->val[i];
    *data = vec;
    sha256_sum_add_doubles(ctx, *data, v->size);
  } break;
  case WS_Float: {
    float *vec= qlua_malloc(L, v->size * sizeof (float));
    for (i = 0; i < v->size; i++)
      vec[i] = v->val[i];
    *data = vec;
    sha256_sum_add_floats(ctx, *data, v->size);
  } break;
  default:
    luaL_error(L, "Unknown precision in w_vecreal()");
  }
  sha256_sum(sum, ctx);
  sha256_destroy(ctx);
  *kind = "VectorReal";
  *filetype = construct_ftype_vecreal(L, b, v->size, opts->wsize);
  *memtype =  construct_mtype_vecreal(L, b, v->size, opts->wsize);
}

static void
w_veccomplex(lua_State *L, mHdf5File *b, mLattice *S,
             struct wopts_s *opts, struct laddr_s *laddr,
             SHA256_Sum *sum, void **data, hid_t *filetype, hid_t *memtype,
             const char **kind)
{
  mVecComplex *v = qlua_checkVecComplex(L, 3);
  SHA256_Context *ctx = sha256_create(L);
  int i;

  switch (opts->wsize) {
  case WS_Double: {
    machine_complex_double *vec= qlua_malloc(L, v->size * sizeof (machine_complex_double));
    for (i = 0; i < v->size; i++) {
      vec[i].re = QLA_real(v->val[i]);
      vec[i].im = QLA_imag(v->val[i]);
    }
    *data = vec;
    sha256_sum_add_doubles(ctx, *data, 2 * v->size);
  } break;
  case WS_Float: {
    machine_complex_float *vec= qlua_malloc(L, v->size * sizeof (machine_complex_float));
    for (i = 0; i < v->size; i++) {
      vec[i].re = QLA_real(v->val[i]);
      vec[i].im = QLA_imag(v->val[i]);
    }
    *data = vec;
    sha256_sum_add_floats(ctx, *data, 2 * v->size);
  } break;
  default:
    luaL_error(L, "Unknown precision in w_veccomplex()");
  }
  sha256_sum(sum, ctx);
  sha256_destroy(ctx);
  *kind = "VectorReal";
  *filetype = construct_ftype_veccomplex(L, b, v->size, opts->wsize);
  *memtype =  construct_mtype_veccomplex(L, b, v->size, opts->wsize);
}

static void
w_latint(lua_State *L, mHdf5File *b, mLattice *S,
         struct wopts_s *opts, struct laddr_s *laddr,
         SHA256_Sum *sum, void **data, hid_t *filetype, hid_t *memtype,
         const char **kind)
{
  int *local_x = qlua_malloc(L, laddr->rank * sizeof (int));
  mLatInt *m = qlua_checkLatInt(L, 3, S);
  SHA256_Context *ctx = sha256_create(L);
  int volume = laddr->volume;
  int rank = laddr->rank;
  int i;

  int *ptr = qlua_malloc(L, volume * sizeof (int));

  CALL_QDP(L);

  *data = ptr;
  *kind = "LatticeInt";
  *filetype =  H5Tcopy(H5T_STD_I32BE);
  CHECK_H5(L, *filetype, "Tcopy(i32be) in w_latint()");
  *memtype = H5Tcopy(H5T_NATIVE_INT);
  CHECK_H5(L, *filetype, "Tcopy(int) in w_latint()");
  
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

  sha256_destroy(ctx);
  qlua_free(L, local_x);
}

static int
write_lat(lua_State *L, mHdf5File *b, const char *path, OutPacker_H5 repack)
{
  struct wopts_s wopts;
  process_wopts(L, &wopts);
  mLattice *S = qlua_ObjLattice(L, 3);
  struct laddr_s laddr;
  laddr.rank = S->rank;
  laddr.low = qlua_malloc(L, laddr.rank * sizeof (int));
  laddr.high = qlua_malloc(L, laddr.rank * sizeof (int));
  hsize_t *offset = qlua_malloc(L, laddr.rank * sizeof (hsize_t));
  hsize_t *stride = qlua_malloc(L, laddr.rank * sizeof (hsize_t));
  hsize_t *count = qlua_malloc(L, laddr.rank * sizeof (hsize_t));
  hsize_t *block = qlua_malloc(L, laddr.rank * sizeof (hsize_t));
  hsize_t *hlatdim = qlua_malloc(L, laddr.rank * sizeof (hsize_t));
  
  qlua_sublattice(laddr.low, laddr.high, QDP_this_node, S);
  int volume, j;
  for (volume = 1, j = 0; j < laddr.rank; j++) {
    int extend_j = laddr.high[j] - laddr.low[j];
    volume *= extend_j;
    offset[j] = laddr.low[j];
    stride[j] = 1;
    count[j] = 1;
    block[j] = extend_j;
    hlatdim[j] = S->dim[j];
  }
  laddr.volume = volume;

  hid_t filetype, memtype;
  void *data;
  SHA256_Sum sum;
  const char *kind;
  (*repack)(L, b, S, &wopts, &laddr, &sum, &data, &filetype, &memtype, &kind);
  qlua_free(L, laddr.low);
  qlua_free(L, laddr.high);
  combine_checksums(&sum, 1);

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
  hid_t filespace = H5Screate_simple(S->rank, hlatdim, NULL);
  qlua_free(L, hlatdim);
  CHECK_H5(L, filespace, "Screate_simple() file space");
  hid_t dataset = H5Dcreate2(wdir, ename, filetype, filespace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  qlua_free(L, dpath);
  CHECK_H5(L, dataset, "Dcreate() data set");
  CHECK_H5(L, H5Tclose(filetype), "Tclose() file type");
  CHECK_H5(L, H5Sclose(filespace), "Sclose() file space");
  if (wdir != b->cwd)
    CHECK_H5(L, H5Gclose(wdir), "Gclose() write dir");
  hid_t memspace = H5Screate_simple(S->rank, block, NULL);
  CHECK_H5(L, memspace, "Screate_simple() mem space");
  filespace =  H5Dget_space(dataset);
  CHECK_H5(L, H5Sselect_hyperslab(filespace, H5S_SELECT_SET, offset, stride, count, block), "Sselect_hyperslab()");
  qlua_free(L, block);
  qlua_free(L, offset);
  qlua_free(L, stride);
  qlua_free(L, count);
  hid_t plist = H5Pcreate(H5P_DATASET_XFER);
  CHECK_H5(L, plist, "Pcreate() xfter plist");
  CHECK_H5(L, H5Dwrite(dataset, memtype, memspace, filespace, plist, data), "Dwrite() data");
  CHECK_H5(L, H5Pclose(plist), "Pclose() xfer plist");
  CHECK_H5(L, H5Sclose(filespace), "Sclose() file space");
  CHECK_H5(L, H5Sclose(memspace), "Sclose() mem space");
  CHECK_H5(L, H5Tclose(memtype), "Tclose() mem type");
  qlua_free(L, data);

  write_attrs(L, b, dataset, &sum, kind);
  CHECK_H5(L, H5Dclose(dataset), "Dclose() data set");

  qlua_Hdf5_leave();

  return 0;
}

static int
write_seq(lua_State *L, mHdf5File *b, const char *path, OutPacker_H5 repack)
{
  struct wopts_s wopts;
  process_wopts(L, &wopts);

  hid_t filetype, memtype;
  void *data;
  SHA256_Sum sum;
  const char *kind;
  (*repack)(L, b, NULL, &wopts, NULL, &sum, &data, &filetype, &memtype, &kind);
  combine_checksums(&sum, b->master);

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
  hid_t filespace = H5Screate(H5S_SCALAR);
  CHECK_H5(L, filespace, "Screate(SCALAR) in write_seq()");
  hid_t dataset = H5Dcreate2(wdir, ename, filetype, filespace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  qlua_free(L, dpath);
  CHECK_H5(L, dataset, "Dcreate() data set");
  CHECK_H5(L, H5Tclose(filetype), "Tclose() file type");
  CHECK_H5(L, H5Sclose(filespace), "Sclose() file space");
  if (wdir != b->cwd)
    CHECK_H5(L, H5Gclose(wdir), "Gclose() write dir");
  hid_t plist = H5Pcreate(H5P_DATASET_XFER);
  CHECK_H5(L, plist, "Pcreate() xfter plist");
  if (b->master)
    CHECK_H5(L, H5Dwrite(dataset, memtype, H5S_ALL, H5S_ALL, plist, data), "Dwrite() data");
  CHECK_H5(L, H5Pclose(plist), "Pclose() xfer plist");
  CHECK_H5(L, H5Tclose(memtype), "Tclose() mem type");
  qlua_free(L, data);

  write_attrs(L, b, dataset, &sum, kind);
  CHECK_H5(L, H5Dclose(dataset), "Dclose() data set");

  qlua_Hdf5_leave();

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
  for (i = 0; qotable[i].qtype != qNoType; i++) {
    if (qotable[i].qtype == kind)
      break;
  }
  if (qotable[i].repack == NULL)
    luaL_error(L, "unwritable data");
  count = (qotable[i].is_parallel? write_lat: write_seq)(L, b, p, qotable[i].repack);
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
  { qString,                 0,  w_string        },
  { qReal,                   0,  w_real          },
  { qComplex,                0,  w_complex       },
  { qVecInt,                 0,  w_vecint        },
  { qVecReal,                0,  w_vecreal       },
  { qVecComplex,             0,  w_veccomplex    },
#if 0 /* XXXXXX scalar uncolored objects */
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
#if 0 /* XXXXXXXXXXXXXXXXXXXXXXXXXX */
  { "read",             qhdf5_real    },
  { "delete",           qhdf5_delete  },
  { "list",             qhdf5_list    },
  { "dims",             qhdf5_dims    },
  { "kind",             qhdf5_kind    },
  { "exists",           qhdf5_exists  },
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
