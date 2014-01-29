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
const char csum_attr_name[] = ".sha256";
const char kind_attr_name[] = ".kind";

#define FAILED_H5CALL -1
#define CHECK_H5(L,expr,message) do { if (expr < 0) luaL_error(L, "HDF5 error %s", message); } while (0)

/* mappings to HDF5 types */
typedef struct htype_s {
  char *name;
  hid_t mtype;
  hid_t ftype;
  struct htype_s *next;
} HType;

struct mHdf5Writer_s {
  hid_t file;
  hid_t cwd;
  int master;
  HType *htype;
};

struct mHdf5Reader_s {
  hid_t file;
  hid_t cwd;
  HType *htype;
};

typedef struct QObjTable_s {
  const char *name;
  H5I_type_t  ftype;
  int (*writer)(lua_State *L, mHdf5Writer *b, const char *path, struct QObjTable_s *qot);
  int (*reader)(lua_State *L, mHdf5Reader *b, const char *path, struct QObjTable_s *qot);
} QObjTable;

static QObjTable qotable[];

/* common helpers */

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
    CHECK_H5(L, dnext, "chpath failed");
    if (isused == 0)
      H5Gclose(dhandle);
    dhandle = dnext;
    isused = 0;
    dname = buf;
  }
  qlua_free(L, dpath);
  return dhandle;
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

static int h_kind(lua_State *L, hid_t file, hid_t cwd, const char *path) { /* XXXXX */ return 0; }

static int h_exists(lua_State *L, hid_t file, hid_t cwd, const char *path) { /* XXXXX */ return 0; }

/* writer */

static void
check_writer(lua_State *L, mHdf5Writer *w)
{
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
  hid_t dnext;
    
  while (buf) {
    strsep(&buf, "/");
    dnext = H5Gopen2(dhandle, dname, H5P_DEFAULT);
    if (dnext < 0) {
      dnext = H5Gcreate2(dhandle, dname, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
      CHECK_H5(L, dnext, "mkpath failed");
    }
    if (isused == 0)
      H5Gclose(dhandle);
    dhandle = dnext;
    isused = 0;
    dname = buf;
  }
  return dhandle;
}

static mHdf5Writer *
qlua_newHdf5Writer(lua_State *L)
{
  mHdf5Writer *h = lua_newuserdata(L, sizeof (mHdf5Writer));
  
  h->master = (QDP_this_node == qlua_master_node);
  h->file = -1;
  h->cwd = -1;
  h->htype = NULL;

  luaL_getmetatable(L, mtnWriter);
  lua_setmetatable(L, -2);
  
  return h;
}

mHdf5Writer *
qlua_checkHdf5Writer(lua_State *L, int idx)
{
  void *v = luaL_checkudata(L, idx, mtnWriter);
  
  luaL_argcheck(L, v != 0, idx, "qcd.hdf5.Writer expected");
  return v;
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

static int
do_w_close(lua_State *L, mHdf5Writer *b)
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
qhdf5_w_gc(lua_State *L)
{
  mHdf5Writer *b = qlua_checkHdf5Writer(L, 1);

  do_w_close(L, b);
  return 0;
}

static int
qhdf5_w_close(lua_State *L)
{
  mHdf5Writer *b = qlua_checkHdf5Writer(L, 1);

  check_writer(L, b);
  // CHECK_H5(L, do_w_close(L, b), "writer close error");
  do_w_close(L, b);
  lua_pushnil(L);
  return 1;
}

static int
qhdf5_w_chpath(lua_State *L)
{
  mHdf5Writer *b = qlua_checkHdf5Writer(L, 1);
  const char *p = luaL_checkstring(L, 2);
  
  check_writer(L, b);
  qlua_Hdf5_enter(L);

  hid_t dir = change_hdf5_path(L, b->file, b->cwd, p);
  if (dir != b->cwd)
    H5Gclose(b->cwd);
  b->cwd = dir;

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

  char *dpath = qlua_strdup(L, p);
  hid_t dh = make_hdf5_path(L, b->file, b->cwd, dpath);
  qlua_free(L, dpath);
  if (dh != b->cwd)
    H5Gclose(dh);

  qlua_Hdf5_leave();
  return 0;
}


static int
qhdf5_w_kind(lua_State *L)
{
  mHdf5Writer *b = qlua_checkHdf5Writer(L, 1);
  const char *p = luaL_checkstring(L, 2);

  return h_kind(L, b->file, b->cwd, p);
}

static int
qhdf5_w_exists(lua_State *L)
{
  mHdf5Writer *b = qlua_checkHdf5Writer(L, 1);
  const char *p = luaL_checkstring(L, 2);

  return h_exists(L, b->file, b->cwd, p);
}

static void
write_attrs(lua_State *L, mHdf5Writer *b, hid_t dset, const SHA256_Sum *sum, const char *kind)
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

  SHA256_Sum gsum = *sum;
  QMP_broadcast(&gsum, sizeof (gsum));
  hsize_t slen = sizeof(gsum.v);
  hid_t stype = H5Tcopy(H5T_STD_U8BE);
  hid_t sds = H5Screate_simple(1, &slen, NULL);
  hid_t sattr = H5Acreate2(dset, csum_attr_name, stype, sds, H5P_DEFAULT, H5P_DEFAULT);
  CHECK_H5(L, H5Awrite(sattr, stype, gsum.v), "write sum attr");
  H5Sclose(sds);
  H5Aclose(sattr);
  H5Tclose(stype);

}

/* XXXX writers */
static int
w_string(lua_State *L, mHdf5Writer *b, const char *path, struct QObjTable_s *qot)
{
  const char *str = luaL_checkstring(L, 3);
  size_t len = strlen(str);
  SHA256_Sum sum;

  sha256_sum_string(&sum, str, len);
  check_writer(L, b);

  qlua_Hdf5_enter(L);

  hid_t ftype = H5Tcopy(H5T_C_S1);
  hid_t mtype = H5Tcopy(H5T_C_S1);
  CHECK_H5(L, H5Tset_size(ftype, len), "set ftype size");
  CHECK_H5(L, H5Tset_size(mtype, len), "set mtype size");
  
  char *dpath = qlua_strdup(L, path);
  hid_t wdir = b->cwd;
  char *ename = strrchr(dpath, '/');
  if (ename == NULL) {
    ename = dpath;
  } else {
    ename[0] = 0;
    ename = ename + 1;
    wdir = make_hdf5_path(L, b->file, b->cwd, dpath);
  }
  hid_t dataspace = H5Screate(H5S_SCALAR);
  hid_t dataset = H5Dcreate2(wdir, ename, ftype, dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  hid_t plist = H5Pcreate(H5P_DATASET_XFER);
  if (b->master) {
    /* Only on the master */
    CHECK_H5(L, H5Dwrite(dataset, mtype, H5S_ALL, H5S_ALL, plist, str), "write string");
  }
  /* Everyone must write attributes, which must be identical or else. */
  write_attrs(L, b, dataset, &sum, qot->name);

  CHECK_H5(L, H5Pclose(plist), "Pclose() plist");
  CHECK_H5(L, H5Dclose(dataset), "Dclose() dataset");
  CHECK_H5(L, H5Sclose(dataspace), "Sclose() dataspace");
  if (wdir != b->cwd)
    CHECK_H5(L, H5Gclose(wdir), "Gclose() write dir");
  CHECK_H5(L, H5Tclose(ftype), "Tclose() ftype");
  CHECK_H5(L, H5Tclose(mtype), "Tclose() mtype");
  qlua_free(L, dpath);
  qlua_Hdf5_enter(L);

  return 0;
}

static int
w_real(lua_State *L, mHdf5Writer *b, const char *path, struct QObjTable_s *qot)
{
  double val = luaL_checknumber(L, 3);
  SHA256_Context *ctx = sha256_create(L);
  SHA256_Sum sum;

  sha256_sum_add_doubles(ctx, &val, 1);
  sha256_sum(&sum, ctx);
  sha256_destroy(ctx);
  check_writer(L, b);

  qlua_Hdf5_enter(L);

  hid_t ftype = H5Tcopy(H5T_IEEE_F64BE);
  hid_t mtype = H5Tcopy(H5T_NATIVE_DOUBLE);
  
  char *dpath = qlua_strdup(L, path);
  hid_t wdir = b->cwd;
  char *ename = strrchr(dpath, '/');
  if (ename == NULL) {
    ename = dpath;
  } else {
    ename[0] = 0;
    ename = ename + 1;
    wdir = make_hdf5_path(L, b->file, b->cwd, dpath);
  }
  hid_t dataspace = H5Screate(H5S_SCALAR);
  hid_t dataset = H5Dcreate2(wdir, ename, ftype, dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  hid_t plist = H5Pcreate(H5P_DATASET_XFER);
  if (b->master) {
    /* Only on the master */
    CHECK_H5(L, H5Dwrite(dataset, mtype, H5S_ALL, H5S_ALL, plist, &val), "write string");
  }
  /* Everyone must write attributes, which must be identical or else. */
  write_attrs(L, b, dataset, &sum, qot->name);

  CHECK_H5(L, H5Pclose(plist), "Pclose() plist");
  CHECK_H5(L, H5Dclose(dataset), "Dclose() dataset");
  CHECK_H5(L, H5Sclose(dataspace), "Sclose() dataspace");
  if (wdir != b->cwd)
    CHECK_H5(L, H5Gclose(wdir), "Gclose() write dir");
  CHECK_H5(L, H5Tclose(ftype), "Tclose() ftype");
  CHECK_H5(L, H5Tclose(mtype), "Tclose() mtype");
  qlua_free(L, dpath);
  qlua_Hdf5_enter(L);

  return 0;
}

static const char *
qh5_data_kind(lua_State *L, int idx)
{
  switch (qlua_qtype(L, idx)) {
  case qString: return "String";
  case qReal: return "Real";
  case qComplex: return "Complex";
  case qMatReal: return "MatrixReal";
  case qMatComplex: return "MatrixComplex";
  case qVecInt: return "VectorInt";
  case qVecReal: return "VectorReal";
  case qVecComplex: return "VectorComplex";
  case qSeqColVec2: case qSeqColVec3: case qSeqColVecN: return "ColorVector";
  case qSeqColMat2: case qSeqColMat3: case qSeqColMatN: return "ColorMatrix";
  case qSeqDirFerm2: case qSeqDirFerm3: case qSeqDirFermN: return "DiracFermion";
  case qSeqDirProp2: case qSeqDirProp3: case qSeqDirPropN: return "DiracPropagator";
  case qLatInt: return "LatticeInt";
  case qLatReal: return "LatticeReal";
  case qLatComplex: return "LatticeComplex";
  case qLatColVec2: case qLatColVec3: case qLatColVecN: return "LatticeColorVector";
  case qLatColMat2: case qLatColMat3: case qLatColMatN: return "LatticeColorMatrix";
  case qLatDirFerm2: case qLatDirFerm3: case qLatDirFermN: return "LatticeDiracFermion";
  case qLatDirProp2: case qLatDirProp3: case qLatDirPropN: return "LatticeDiracPropagator";
  case qTable: {
    const char *kind = NULL;
    lua_pushinteger(L, 1);
    lua_gettable(L, idx);
    switch (qlua_qtype(L, -1)) {
    case qLatColMat2: case qLatColMat3: case qLatColMatN: kind = "Gauge"; break;
    case qLatDirFerm2: case qLatDirFerm3: case qLatDirFermN: kind = "DomainWallFermion"; break;
    case qLatDirProp2: case qLatDirProp3: case qLatDirPropN: kind = "DomainWallPropagator"; break;
    default:
      break;
    }
    lua_pop(L, 1);
    if (kind)
      return kind;
  } break;
  default:
    break;
  }
  return "unknown";
}

static int
qhdf5_w_write(lua_State *L)
{
  mHdf5Writer *b = qlua_checkHdf5Writer(L, 1);
  const char *p = luaL_checkstring(L, 2);
  const char *kind = qh5_data_kind(L, 3);
  int count;
  int i;

  check_writer(L, b);
  qlua_Hdf5_enter(L);
  for (i = 0; qotable[i].name; i++) {
    if (strcmp(qotable[i].name, kind) == 0)
      break;
  }
  if (qotable[i].writer == NULL)
    luaL_error(L, "unwritable data of kind %s", kind);
  count = qotable[i].writer(L, b, p, &qotable[i]);

  qlua_Hdf5_leave();
  return count;
}

static int
q_hdf5_writer(lua_State *L)
{
  const char *name = luaL_checkstring(L, 1);
  mHdf5Writer *w = qlua_newHdf5Writer(L);
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

/* reader */

static void
check_reader(lua_State *L, mHdf5Reader *r)
{
  if (r->file < 0)
    luaL_error(L, "closed hdf5 reader");
}

/* allocation */
static mHdf5Reader *
qlua_newHdf5Reader(lua_State *L)
{
  mHdf5Reader *h = lua_newuserdata(L, sizeof (mHdf5Reader));

  h->file = -1;
  h->cwd = -1;
  h->htype = NULL;

  luaL_getmetatable(L, mtnReader);
  lua_setmetatable(L, -2);
  return h;
}

mHdf5Reader *
qlua_checkHdf5Reader(lua_State *L, int idx)
{
  void *v = luaL_checkudata(L, idx, mtnReader);

  luaL_argcheck(L, v != 0, idx, "qcd.hdf5.Reader expected");
  return v;
}

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
do_r_close(lua_State *L, mHdf5Reader *b)
{
  qlua_Hdf5_enter(L);

  free_types(L, b->htype);
  b->htype = NULL;
  if (b->file > 0)
    H5Fclose(b->file);
  b->file = -1;

  qlua_Hdf5_leave();
  return 0;
}

static int
qhdf5_r_gc(lua_State *L)
{
  mHdf5Reader *b = qlua_checkHdf5Reader(L, 1);

  return do_r_close(L, b);
}

static int
qhdf5_r_close(lua_State *L)
{
  mHdf5Reader *b = qlua_checkHdf5Reader(L, 1);

  check_reader(L, b);
  return do_r_close(L, b);
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
  CHECK_H5(L, dh, "qcd.hdf5.Reader:list() failed");
  ec = H5Gget_info(dh, &gi);
  CHECK_H5(L, ec, "Gget_info() failed");
  lua_createtable(L, gi.nlinks, 0);
  lua_pushinteger(L, 1);
  ec = H5Literate(dh, H5_INDEX_NAME, H5_ITER_INC, NULL, qh_get_list, L);
  CHECK_H5(L, ec, "Literate() failed");
  H5Gclose(dh);
  lua_pop(L, 1);
  return 1;
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
qhdf5_r_kind(lua_State *L)
{
  mHdf5Reader *b = qlua_checkHdf5Reader(L, 1);
  const char *p = luaL_checkstring(L, 2);

  return h_kind(L, b->file, b->cwd, p);
}

static int
qhdf5_r_exists(lua_State *L)
{
  mHdf5Reader *b = qlua_checkHdf5Reader(L, 1);
  const char *p = luaL_checkstring(L, 2);

  return h_exists(L, b->file, b->cwd, p);
}

/* XXXX readers */
static int r_string(lua_State *L, mHdf5Reader *b, const char *path, struct QObjTable_s *qot) { /* XXXX */ return 0; }
static int r_real(lua_State *L, mHdf5Reader *b, const char *path, struct QObjTable_s *qot) { /* XXXX */ return 0; }

static int qhdf5_r_read(lua_State *L) { /* XXXXX */ return 0; }

static int
q_hdf5_reader(lua_State *L)
{
  const char *name = luaL_checkstring(L, 1);
  mHdf5Reader *r = qlua_newHdf5Reader(L);
  MPI_Comm comm = MPI_COMM_WORLD;
  MPI_Info info = MPI_INFO_NULL;
  hid_t acc_tpl1;

  qlua_Hdf5_enter(L);

  acc_tpl1 = H5Pcreate(H5P_FILE_ACCESS);
  CHECK_H5(L, acc_tpl1, "Pcreate() failed");
  CHECK_H5(L, H5Pset_fapl_mpio(acc_tpl1, comm, info), "Pset_fapl_mpio() failed");
  r->file = H5Fopen(name, H5F_ACC_RDONLY, acc_tpl1);
  CHECK_H5(L, r->file, "Fopen() failed");
  r->cwd = H5Gopen2(r->file, "/", H5P_DEFAULT);
  CHECK_H5(L, r->cwd, "Gopen(\"/\") failed");

  qlua_Hdf5_leave();
  
  return 1;
}

/* setup */
static QObjTable qotable[] = {
  { "String",                 H5I_DATASET,  w_string,      r_string     },
  { "Real",                   H5I_DATASET,  w_real,        r_real       },
#if 0 /* XXXXX qlua object dispatch table */
  { "Complex",                H5I_DATASET,  w_complex,     r_complex    },
  { "MatrixReal",             H5I_DATASET,  w_matreal,     r_matreal    },
  { "MatrixComplex",          H5I_DATASET,  w_matcomplex,  r_matcomplex },
  { "VectorInt",              H5I_DATASET,  w_varint,      r_varint     },
  { "VectorReal",             H5I_DATASET,  w_vecreal,     r_vecreal    },
  { "VectorComplex",          H5I_DATASET,  w_veccomplex,  r_veccomplex },
  { "ColorVector",            H5I_DATASET,  w_colvec,      r_colvec     },
  { "ColorMatrix",            H5I_DATASET,  w_colmat,      r_colmat     },
  { "DiracFermion",           H5I_DATASET,  w_dirferm,     r_dirferm    },
  { "DiracPropagator",        H5I_DATASET,  w_dirprop,     r_dirprop    },
  { "LatticeInt",             H5I_DATASET,  w_latint,      r_latint     },
  { "LatticeReal",            H5I_DATASET,  w_latreal,     r_latreal    },
  { "LatticeComplex",         H5I_DATASET,  w_latcomplex,  r_latcomplex },
  { "LatticeColorVector",     H5I_DATASET,  w_latcolvec,   r_latcolvec  },
  { "LatticeColorMatrix",     H5I_DATASET,  w_latcolmat,   r_latcolmat  },
  { "LatticeDiracFermion",    H5I_DATASET,  w_latdirferm,  r_latdirferm },
  { "LatticeDiracPropagator", H5I_DATASET,  w_latdirprop,  r_latdirprop },
  { "Gauge",                  H5I_GROUP,    w_gauge,       r_gauge      },
  { "DomainWallFermion",      H5I_GROUP,    w_dwferm,      r_dwferm     },
  { "DomainWallPropagator",   H5I_GROUP,    w_dwprop,      r_dwprop     },
#endif /* XXXXX qlua object dispatch table */
  { NULL,                     H5I_BADID,    NULL,          NULL         }
};


/* metatables */
static const struct luaL_Reg mtReader[] = {
  { "__tostring",       qhdf5_r_fmt },
  { "__gc",             qhdf5_r_gc },
  { "close",            qhdf5_r_close },
  { "read",             qhdf5_r_read },
  { "chpath",           qhdf5_r_chpath },
  { "kind",             qhdf5_r_kind },
  { "exists",           qhdf5_r_exists },
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
  { "exists",           qhdf5_w_exists },
  { NULL,               NULL}
};

/* names and routines for qcd.hdf5 table */
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
