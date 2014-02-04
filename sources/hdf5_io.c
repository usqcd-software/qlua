#include "modules.h"                                                 /* DEPS */
#include "qlua.h"                                                    /* DEPS */
#include "qcomplex.h"                                                /* DEPS */
#include "qvector.h"                                                 /* DEPS */
#include "qmatrix.h"                                                 /* DEPS */
#include "hdf5_io.h"                                                 /* DEPS */
#include "sha256.h"                                                  /* DEPS */
#include "qlayout.h"                                                 /* DEPS */
#include "lattice.h"                                                 /* DEPS */
#include "seqcolvec.h"                                               /* DEPS */
#include "seqcolmat.h"                                               /* DEPS */
#include "seqdirferm.h"                                              /* DEPS */
#include "seqdirprop.h"                                              /* DEPS */
#include "latint.h"                                                  /* DEPS */
#include "latreal.h"                                                 /* DEPS */
#include "latcomplex.h"                                              /* DEPS */
#include "latcolvec.h"                                               /* DEPS */
#include "latcolmat.h"                                               /* DEPS */
#include "latdirferm.h"                                              /* DEPS */
#include "latdirprop.h"                                              /* DEPS */
#include <string.h>
#include <complex.h>
#include <sys/stat.h>
#include "hdf5.h"
#include "qmp.h"

static const char hdf5_io[] = "hdf5";
static const char mtnFile[] = "qcd.hdf5.mtFile";

static const char csum_attr_name[] = ".sha256";
static const char kind_attr_name[] = ".kind";
static const char time_attr_name[] = ".time";

static const char htn_complex_fmt[]         = ".Complex%s";
static const char htn_vector_fmt[]          = ".Vector%s%s%d";
static const char htn_matrix_fmt[]          = ".Matrix%s%s%dx%d";
static const char htn_color_fmt[]           = ".Color%s%s%d";
static const char htn_dirac_fmt[]           = ".Dirac%s%s%d";

typedef enum {
  sScalar,
  sLattice,
  sOther
} SpaceOfH;

typedef enum {
  kString,
  kReal,
  kComplex,
  kMatrixReal,
  kMatrixComplex,
  kVectorInt,
  kVectorReal,
  kVectorComplex,
  kColorVector,
  kColorMatrix,
  kDiracFermion,
  kDiracPropagator,
  kLatticeInt,
  kLatticeReal,
  kLatticeComplex,
  kLatticeColorVector,
  kLatticeColorMatrix,
  kLatticeDiracFermion,
  kLatticeDiracPropagator,
  kGroup,
  kDataSpace,
  kDataSet,
  kDataType,
  kAttribute,
  kFile,
  kUnknown,
  kNoKind
} KindOfH;

static const char knString[]                  = "String";
static const char knReal[]                    = "Real";
static const char knComplex[]                 = "Complex";
static const char knMatrixReal[]              = "MatrixReal";
static const char knMatrixComplex[]           = "MatrixComplex";
static const char knVectorInt[]               = "VectorInt";
static const char knVectorReal[]              = "VectorReal";
static const char knVectorComplex[]           = "VectorComplex";
static const char knColorVector[]             = "ColorVector";
static const char knColorMatrix[]             = "ColorMatrix";
static const char knDiracFermion[]            = "DiracFermion";
static const char knDiracPropagator[]         = "DiracPropagator";
static const char knLatticeInt[]              = "LatticeInt";
static const char knLatticeReal[]             = "LatticeReal";
static const char knLatticeComplex[]          = "LatticeComplex";
static const char knLatticeColorVector[]      = "LatticeColorVector";
static const char knLatticeColorMatrix[]      = "LatticeColorMatrix";
static const char knLatticeDiracFermion[]     = "LatticeDiracFermion";
static const char knLatticeDiracPropagator[]  = "LatticeDiracPropagator";
static const char knUnknown[]                 = "Unknown";

static const struct {
  KindOfH kind;
  const char *name;
} knTable[] = {
  { kString,                  knString                  },
  { kReal,                    knReal                    },
  { kComplex,                 knComplex                 },
  { kMatrixReal,              knMatrixReal              },
  { kMatrixComplex,           knMatrixComplex           },
  { kVectorInt,               knVectorInt               },
  { kVectorReal,              knVectorReal              },
  { kVectorComplex,           knVectorComplex           },
  { kColorVector,             knColorVector             },
  { kColorMatrix,             knColorMatrix             },
  { kDiracFermion,            knDiracFermion            },
  { kDiracPropagator,         knDiracPropagator         },
  { kLatticeInt,              knLatticeInt              },
  { kLatticeReal,             knLatticeReal             },
  { kLatticeComplex,          knLatticeComplex          },
  { kLatticeColorVector,      knLatticeColorVector      },
  { kLatticeColorMatrix,      knLatticeColorMatrix      },
  { kLatticeDiracFermion,     knLatticeDiracFermion     },
  { kLatticeDiracPropagator,  knLatticeDiracPropagator  },
  { kGroup,                   "Group"                   },
  { kDataSpace,               "DataSpace"               },
  { kDataSet,                 "DataSet"                 },
  { kDataType,                "DataType"                },
  { kAttribute,               "Attribute"               },
  { kFile,                    "File"                    },
  { kUnknown,                 knUnknown                 },
  { kNoKind,                  NULL                      }
};

#define CHECK_H5(L,expr,message) do { if (expr < 0) luaL_error(L, "HDF5 error %s", message); } while (0)

typedef struct {
  double re;
  double im;
} machine_complex_double;

typedef struct {
  float re;
  float im;
} machine_complex_float;

struct mHdf5File_s {
  hid_t file;
  hid_t cwd;
  int master;
  int writer;
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

/* common helpers */
static const char *
kind2name(KindOfH k)
{
  int i;
  for (i = 0; knTable[i].name; i++) {
    if (k == knTable[i].kind)
      return knTable[i].name;
  }
  QLUA_ABORT("Unknown hdf5 kind");
  return knUnknown;
}

static KindOfH
name2kind(const char *name)
{
  int i;
  for (i = 0; knTable[i].name; i++) {
    if (strcmp(knTable[i].name, name) == 0)
      return knTable[i].kind;
  }
  QLUA_ABORT("Unknown hdf5 kind");
  return kUnknown;
}

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

static hid_t
read_htype(lua_State *L, mHdf5File *hf, const char *name, int ftype_p)
{
  if (ftype_p)
    return H5Topen2(hf->file, name, H5P_DEFAULT);
  return -1;
}

static hid_t
write_htype(lua_State *L, mHdf5File *hf, const char *name, int ftype_p, hid_t handle)
{
  if (!ftype_p || !hf->writer)
    return handle;
  CHECK_H5(L, H5Tcommit2(hf->file, name, handle, H5P_DEFAULT,  H5P_DEFAULT,  H5P_DEFAULT), "Tcommit() in write_type");
  CHECK_H5(L, H5Tclose(handle), "Tclose() in write_type");
  handle = H5Topen2(hf->file, name, H5P_DEFAULT);
  CHECK_H5(L, handle, "Topen() in write_type");
  return handle;
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
get_time_type(lua_State *L, mHdf5File *b, int ftype_p)
{
  hid_t v = ftype_p? H5Tcopy(H5T_STD_I64BE): H5Tcopy(H5T_NATIVE_LLONG);
  CHECK_H5(L, v, "Tcopy() in time_type");
  return v;
}

static hid_t
get_sha256_type(lua_State *L, mHdf5File *b, int ftype_p)
{
  hid_t ct = H5Tcopy(H5T_STD_U8BE);
  CHECK_H5(L, ct, "sha256_type(): copying u8 type");
  hsize_t vsize = sizeof (SHA256_Sum);
  hid_t v = H5Tarray_create(ct, 1, &vsize);
  CHECK_H5(L, v, "sha256_type(): set size");
  CHECK_H5(L, H5Tclose(ct), "sha256_type(): close u8 type");
  return v;
}

static hid_t
get_string_type(lua_State *L, mHdf5File *b, int size, int ftype_p)
{
  hid_t v = H5Tcopy(H5T_C_S1);
  CHECK_H5(L, v, "string_type(): copying char type");
  CHECK_H5(L, H5Tset_size(v, size), "string_type(): set size");
  return v;
}

////////
static hid_t
get_int_type(lua_State *L, mHdf5File *b, int ftype_p)
{
  hid_t v;
  if (ftype_p) {
    v = H5T_STD_I32BE;
  } else {
    v = H5T_NATIVE_INT;
  }
  v = H5Tcopy(v);
  CHECK_H5(L, v, "int_type(): copying int type");
  return v;
}

static hid_t
get_real_type(lua_State *L, mHdf5File *b, WriteSize wsize, int ftype_p)
{
  hid_t v = -1;
  if (ftype_p) {
    switch (wsize) {
    case WS_Double: v = H5T_IEEE_F64BE; break;
    case WS_Float: v = H5T_IEEE_F32BE; break;
    }
  } else {
    switch (wsize) {
    case WS_Double: v = H5T_NATIVE_DOUBLE; break;
    case WS_Float: v = H5T_NATIVE_FLOAT; break;
    }
  }
  CHECK_H5(L, v, "real_type(): unknown wsize");
  v = H5Tcopy(v);
  CHECK_H5(L, v, "real_type(): copying real type");
  return v;
}

static hsize_t
get_complex_size(lua_State *L, WriteSize wsize, int ftype_p)
{
  if (ftype_p) {
    switch (wsize) {
    case WS_Double: return 16;
    case WS_Float: return 8;
    default:
      break;
    }
  } else {
    switch (wsize) {
    case WS_Double: return sizeof (machine_complex_double);
    case WS_Float: return sizeof (machine_complex_float);
    default:
      break;
    }
  }
  CHECK_H5(L, -1, "complex_size(): unknown wsize");
  return -1;
}

static hsize_t
get_real_offset(lua_State *L, WriteSize wsize, int ftype_p)
{
  if (ftype_p) {
    switch (wsize) {
    case WS_Double: return 0;
    case WS_Float: return 0;
    default:
      break;
    }
  } else {
    switch (wsize) {
    case WS_Double: return HOFFSET(machine_complex_double, re);
    case WS_Float: return HOFFSET(machine_complex_float, re);
    default:
      break;
    }
  }
  CHECK_H5(L, -1, "real_offset(): unknown wsize");
  return -1;
}

static hsize_t
get_imag_offset(lua_State *L, WriteSize wsize, int ftype_p)
{
  if (ftype_p) {
    switch (wsize) {
    case WS_Double: return 8;
    case WS_Float: return 4;
    default:
      break;
    }
  } else {
    switch (wsize) {
    case WS_Double: return HOFFSET(machine_complex_double, im);
    case WS_Float: return HOFFSET(machine_complex_float, im);
    default:
      break;
    }
  }
  CHECK_H5(L, -1, "imag_offset(): unknown wsize");
  return -1;
}

static hid_t
get_complex_type(lua_State *L, mHdf5File *hf, WriteSize wsize, int ftype_p)
{
  char tname[72];
  sprintf(tname, htn_complex_fmt, wsize_name(L, wsize));
  hid_t v = read_htype(L, hf, tname, ftype_p);
  if (v < 0) {
    hid_t ct = get_real_type(L, hf, wsize, ftype_p);
    v = H5Tcreate(H5T_COMPOUND, get_complex_size(L, wsize, ftype_p));
    CHECK_H5(L, v, "complex_type(): creating compound");
    CHECK_H5(L, H5Tinsert(v, "r", get_real_offset(L, wsize, ftype_p), ct), "complex_type(): insert real");
    CHECK_H5(L, H5Tinsert(v, "i", get_imag_offset(L, wsize, ftype_p), ct), "complex_type(): insert imag");
    CHECK_H5(L, H5Tclose(ct), "complex_type(): close real type");
    v = write_htype(L, hf, tname, ftype_p, v);
  }
  return v;
}

static hid_t
get_vecint_type(lua_State *L, mHdf5File *hf, int n, int ftype_p)
{
  char tname[72];
  sprintf(tname, htn_vector_fmt, "Int", "", n);
  hid_t v = read_htype(L, hf, tname, ftype_p);
  if (v < 0) {
    hid_t ct = get_int_type(L, hf, ftype_p);
    hsize_t vsize = n;
    v = H5Tarray_create(ct, 1, &vsize);
    CHECK_H5(L, v, "vecint_type(): creating vecint type");
    CHECK_H5(L, H5Tclose(ct), "vecint_type(): closing int type");
    v = write_htype(L, hf, tname, ftype_p, v);
  }
  return v;
}

static hid_t
get_vecreal_type(lua_State *L, mHdf5File *hf, int n, WriteSize wsize, int ftype_p)
{
  char tname[72];
  sprintf(tname, htn_vector_fmt, "Real", wsize_name(L, wsize), n);
  hid_t v = read_htype(L, hf, tname, ftype_p);
  if (v < 0) {
    hid_t ct = get_real_type(L, hf, wsize, ftype_p);
    hsize_t vsize = n;
    v = H5Tarray_create(ct, 1, &vsize);
    CHECK_H5(L, v, "vecreal_type(): creating vecreal type");
    CHECK_H5(L, H5Tclose(ct), "vecreal_type(): closing real type");
    v = write_htype(L, hf, tname, ftype_p, v);
  }
  return v;
}

static hid_t
get_matreal_type(lua_State *L, mHdf5File *hf, int n, int m, WriteSize wsize, int ftype_p)
{
  char tname[72];
  sprintf(tname, htn_matrix_fmt, "Real", wsize_name(L, wsize), n, m);
  hid_t v = read_htype(L, hf, tname, ftype_p);
  if (v < 0) {
    hid_t ct = get_real_type(L, hf, wsize, ftype_p);
    hsize_t msize[2];
    msize[0] = n; msize[1] = m;
    v = H5Tarray_create(ct, 2, msize);
    CHECK_H5(L, v, "matreal_type(): creating matreal type");
    CHECK_H5(L, H5Tclose(ct), "matreal_type(): closing real type");
    v = write_htype(L, hf, tname, ftype_p, v);
  }
  return v;
}

static hid_t
get_veccomplex_type(lua_State *L, mHdf5File *hf, int n, WriteSize wsize, int ftype_p)
{
  char tname[72];
  sprintf(tname, htn_vector_fmt, "Complex", wsize_name(L, wsize), n);
  hid_t v = read_htype(L, hf, tname, ftype_p);
  if (v < 0) {
    hid_t ct = get_complex_type(L, hf, wsize, ftype_p);
    hsize_t vsize = n;
    v = H5Tarray_create(ct, 1, &vsize);
    CHECK_H5(L, v, "veccomplex_type(): creating veccomplex type");
    CHECK_H5(L, H5Tclose(ct), "veccomplex_type(): closing complex type");
    v = write_htype(L, hf, tname, ftype_p, v);
  }
  return v;
}

static hid_t
get_matcomplex_type(lua_State *L, mHdf5File *hf, int n, int m, WriteSize wsize, int ftype_p)
{
  char tname[72];
  sprintf(tname, htn_matrix_fmt, "Complex", wsize_name(L, wsize), n, m);
  hid_t v = read_htype(L, hf, tname, ftype_p);
  if (v < 0) {
    hid_t ct = get_complex_type(L, hf, wsize, ftype_p);
    hsize_t msize[2];
    msize[0] = n; msize[1] = m;
    v = H5Tarray_create(ct, 2, msize);
    CHECK_H5(L, v, "matcomplex_type(): creating matcomplex type");
    CHECK_H5(L, H5Tclose(ct), "matcomplex_type(): closing complex type");
    v = write_htype(L, hf, tname, ftype_p, v);
  }
  return v;
}

static hid_t
get_colvec_type(lua_State *L, mHdf5File *hf, int nc, WriteSize wsize, int ftype_p)
{
  char tname[72];
  sprintf(tname, htn_color_fmt, "Vector", wsize_name(L, wsize), nc);
  hid_t v = read_htype(L, hf, tname, ftype_p);
  if (v < 0) {
    hid_t ct = get_complex_type(L, hf, wsize, ftype_p);
    hsize_t vsize = nc;
    v = H5Tarray_create(ct, 1, &vsize);
    CHECK_H5(L, v, "colvec_type(): creating colvec type");
    CHECK_H5(L, H5Tclose(ct), "colvec_type(): closing complex type");
    v = write_htype(L, hf, tname, ftype_p, v);
  }
  return v;
}

static hid_t
get_colmat_type(lua_State *L, mHdf5File *hf, int nc, WriteSize wsize, int ftype_p)
{
  char tname[72];
  sprintf(tname, htn_color_fmt, "Matrix", wsize_name(L, wsize), nc);
  hid_t v = read_htype(L, hf, tname, ftype_p);
  if (v < 0) {
    hid_t ct = get_complex_type(L, hf, wsize, ftype_p);
    hsize_t msize[2];
    msize[0] = nc; msize[1] = nc;
    v = H5Tarray_create(ct, 2, msize);
    CHECK_H5(L, v, "colmat_type(): creating colmat type");
    CHECK_H5(L, H5Tclose(ct), "colmat_type(): closing complex type");
    v = write_htype(L, hf, tname, ftype_p, v);
  }
  return v;
}

static hid_t
get_dirferm_type(lua_State *L, mHdf5File *hf, int nc, WriteSize wsize, int ftype_p)
{
  char tname[72];
  sprintf(tname, htn_dirac_fmt, "Fermion", wsize_name(L, wsize), nc);
  hid_t v = read_htype(L, hf, tname, ftype_p);
  if (v < 0) {
    hid_t ct = get_colvec_type(L, hf, nc, wsize, ftype_p);
    hsize_t vsize = QDP_Ns;
    v = H5Tarray_create(ct, 1, &vsize);
    CHECK_H5(L, v, "dirferm_type(): creating dirferm type");
    CHECK_H5(L, H5Tclose(ct), "dirferm_type(): closing colvec type");
    v = write_htype(L, hf, tname, ftype_p, v);
  }
  return v;
}

static hid_t
get_dirprop_type(lua_State *L, mHdf5File *hf, int nc, WriteSize wsize, int ftype_p)
{
  char tname[72];
  sprintf(tname, htn_dirac_fmt, "Propagator", wsize_name(L, wsize), nc);
  hid_t v = read_htype(L, hf, tname, ftype_p);
  if (v < 0) {
    hid_t ct = get_colmat_type(L, hf, nc, wsize, ftype_p);
    hsize_t msize[] = {QDP_Ns, QDP_Ns};
    v = H5Tarray_create(ct, 2, msize);
    CHECK_H5(L, v, "dirprop_type(): creating dirprop type");
    CHECK_H5(L, H5Tclose(ct), "dirprop_type(): closing colmat type");
    v = write_htype(L, hf, tname, ftype_p, v);
  }
  return v;
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

  while (buf && buf[0]) {
    strsep(&buf, "/");
    hid_t dnext = H5Gopen2(dhandle, dname, H5P_DEFAULT);
    if (dnext < 0) {
      dnext = H5Gcreate2(dhandle, dname, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
      CHECK_H5(L, dnext, "mkpath failed");
    }
    if (dhandle != cwd)
      CHECK_H5(L, H5Gclose(dhandle), "Gclose() failed");
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

  luaL_getmetatable(L, mtnFile);
  lua_setmetatable(L, -2);

  return h;
}

static mHdf5File *
qlua_newHdf5Reader(lua_State *L)
{
  mHdf5File *h = lua_newuserdata(L, sizeof (mHdf5File));

  h->master = (QDP_this_node == qlua_master_node);
  h->writer = 0;
  h->file = -1;
  h->cwd = -1;

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

static herr_t
get_list(hid_t loc_id, const char *name, const H5L_info_t *linfo, void *opdata)
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
qhdf5_list(lua_State *L)
{
  mHdf5File *b = qlua_checkHdf5File(L, 1);
  const char *p = luaL_checkstring(L, 2);

  check_file(L, b);
  qlua_Hdf5_enter(L);
  hid_t dh = H5Gopen(p[0] == '/'? b->file: b->cwd, p, H5P_DEFAULT);
  CHECK_H5(L, dh, "qcd.hdf5.Reader:list() failed");
  H5G_info_t gi;
  CHECK_H5(L, H5Gget_info(dh, &gi), "Gget_info() failed");
  lua_createtable(L, gi.nlinks, 0);
  lua_pushinteger(L, 1);
  CHECK_H5(L, H5Literate(dh, H5_INDEX_NAME, H5_ITER_INC, NULL, get_list, L), "Literate() failed");
  CHECK_H5(L, H5Gclose(dh), "Gclose() failed");
  lua_pop(L, 1);
  qlua_Hdf5_leave();
  return 1;
}

static KindOfH
get_h5_kind(lua_State *L, mHdf5File *b, hid_t obj)
{
  switch (H5Iget_type(obj)) {
  case H5I_FILE: return kFile;
  case H5I_GROUP: return kGroup;
  case H5I_DATATYPE: return kDataType;
  case H5I_DATASPACE: return kDataSpace;
  case H5I_ATTR: return kAttribute;
  case H5I_DATASET: {
    hid_t a = H5Aopen(obj, kind_attr_name, H5P_DEFAULT);
    KindOfH kind = kDataSet;
    if (a >= 0) {
      int len = H5Aget_storage_size(a);
      if (len > 0) {
        char *val = qlua_malloc(L, len + 1);
        hid_t mtype = get_string_type(L, b, len, 0);
        CHECK_H5(L, H5Aread(a, mtype, val), "Aread() failed");
        CHECK_H5(L, H5Tclose(mtype), "Tclose() failed");
        val[len] = 0;
        kind = name2kind(val);
        qlua_free(L, val);
      }
      CHECK_H5(L, H5Aclose(a), "Aclose() failed");
    }
    return kind;
  }
  default:
    break;
  }
  return kUnknown;
}

static int
get_h5_time(lua_State *L, mHdf5File *b, hid_t obj, long long *pt)
{
  hid_t a = H5Aopen(obj, time_attr_name, H5P_DEFAULT);
  if (a >= 0) {
    hid_t mtype = get_time_type(L, b, 0);
    CHECK_H5(L, H5Aread(a, mtype, pt), "Aread() failed");
    CHECK_H5(L, H5Tclose(mtype), "Tclose() failed");
    return 1;
  }
  return 0;
}

static SpaceOfH
get_h5_space(lua_State *L, mHdf5File *b, hid_t obj, int *rank, hsize_t **dims)
{
  SpaceOfH s = sOther;
  hid_t space = H5Dget_space(obj);
  if (space >= 0) {
    if (H5Sis_simple(space)) {
      int nd = H5Sget_simple_extent_ndims(space);
      s = sScalar;
      if (nd > 0) {
        *dims = qlua_malloc(L, nd * sizeof (hsize_t));
        *rank = nd;
        H5Sget_simple_extent_dims(space, NULL, dims);
        s = sLattice;
      }
    }
    CHECK_H5(L, H5Sclose(space), "Sclose() failed");
  }
  return s;
}

static int
get_h5_type_name(lua_State *L, mHdf5File *b, hid_t obj, char **ptr)
{
  hid_t tid = H5Dget_type(obj);
  int present = 0;
  if (tid >= 0) {
    if (H5Tequal(tid, H5T_IEEE_F64BE) || H5Tequal(tid, H5T_IEEE_F32BE)) {
      *ptr = qlua_strdup(L, "real");
      present = 1;
    } else {
      ssize_t len = H5Iget_name(tid, NULL, 0);
      if (len > 0) {
        *ptr = qlua_malloc(L, len + 3);
        H5Iget_name(tid, *ptr, len + 2);
        (*ptr)[len + 2] = 0;
        present = 1;
      }
    }
    CHECK_H5(L, H5Tclose(tid), "Tclose() failed");
  }
  return present;
}

static int
qhdf5_stat(lua_State *L)
{
  mHdf5File *b = qlua_checkHdf5File(L, 1);
  const char *path = luaL_checkstring(L, 2);

  check_file(L, b);
  qlua_Hdf5_enter(L);
  hid_t obj = H5Oopen(path[0] == '/'? b->file: b->cwd, path, H5P_DEFAULT);
  if (obj < 0) {
    lua_pushnil(L);
    return 1;
  }
  KindOfH kind = get_h5_kind(L, b, obj);
  lua_createtable(L, 0, 4);
  lua_pushstring(L, kind2name(kind));
  lua_setfield(L, -2, "kind");
  long long t;
  if (get_h5_time(L, b, obj, &t)) {
    lua_pushnumber (L, t);
    lua_setfield(L, -2, "time");
  }
  int rank = 0;
  hsize_t *dims = NULL;
  SpaceOfH space = get_h5_space(L, b, obj, &rank, &dims);
  switch (space) {
  case sScalar:
    lua_pushstring(L, "scalar");
    lua_setfield(L, -2, "geometry");
    break;
  case sLattice: {
    int i;
    lua_createtable(L, 0, rank);
    for (i = 0; i < rank; i++) {
      lua_pushnumber(L, i+1);
      lua_pushnumber(L, dims[i]);
      lua_settable(L, -3);
    }
    lua_setfield(L, -2, "geometry");
    qlua_free(L, dims);
  } break;
  case sOther:
    break;
  default:
    QLUA_ABORT("unkown SpaceOfH value");
  }
  char *tname = NULL;
  if (get_h5_type_name(L, b, obj, &tname)) {
    lua_pushstring(L, tname);
    lua_setfield(L, -2, "type");
    qlua_free(L, tname);
  }
  CHECK_H5(L, H5Oclose(obj), "Oclose() failed");
  qlua_Hdf5_leave();
  return 1;
}

static int
qhdf5_cwd(lua_State *L)
{
  mHdf5File *b = qlua_checkHdf5File(L, 1);
  int res = 0;
  check_file(L, b);
  qlua_Hdf5_enter(L);
  int len = H5Iget_name(b->cwd, NULL, 0);
  if (len > 0) {
    char *path = qlua_malloc(L, len + 1);
    H5Iget_name(b->cwd, path, len + 1);
    lua_pushstring(L, path);
    res = 1;
  }
  qlua_Hdf5_leave();
  return res;
}

static int
qhdf5_remove(lua_State *L)
{
  mHdf5File *b = qlua_checkHdf5Writer(L, 1);
  const char *path = luaL_checkstring(L, 2);

  check_writer(L, b);
  qlua_Hdf5_enter(L);
  char *dpath = qlua_strdup(L, path);
  char *ename = strrchr(dpath, '/');
  hid_t wdir = b->cwd;
  if (ename == NULL) {
    ename = dpath;
  } else {
    ename[0] = 0;
    ename = ename + 1;
    char *dname = path[0] == '/'? &dpath[1]: dpath;
    char *buf = dname;
    hid_t dhandle = path[0] == '/'? H5Gopen2(b->file, "/", H5P_DEFAULT): b->cwd;
    while (buf && buf[0]) {
      strsep(&buf, "/");
      hid_t dnext = H5Gopen2(dhandle, dname, H5P_DEFAULT);
      if (dhandle != b->cwd)
        CHECK_H5(L, H5Gclose(dhandle), "Gclose() failed");
      if (dnext >= 0) {
        dhandle = dnext;
        dname = buf;
        continue;
      }
      qlua_free(L, dpath);
      return 0;
    }
    wdir = dhandle;
  }
  H5Ldelete(wdir, ename, H5P_DEFAULT);
  if (wdir != b->cwd)
    CHECK_H5(L, H5Gclose(wdir), "Gclose() failed after delete");
  qlua_free(L, dpath);
  qlua_Hdf5_leave();
  return 0;
}

static int
qhdf5_flush(lua_State *L)
{
  mHdf5File *b = qlua_checkHdf5Writer(L, 1);

  check_writer(L, b);
  qlua_Hdf5_enter(L);
  CHECK_H5(L, H5Fflush(b->file, H5F_SCOPE_GLOBAL), "Fflush() failed");
  qlua_Hdf5_leave();
  return 0;
}

static void
write_attrs(lua_State *L, mHdf5File *b, hid_t dset, const SHA256_Sum *sum, const char *kind)
{
  hid_t dspace = H5Screate(H5S_SCALAR);
  CHECK_H5(L, dspace, "Screate() in write_attrs()");

  hsize_t klen = strlen(kind);
  hid_t ktype = get_string_type(L, b, klen + 1, 0);
  hid_t kattr = H5Acreate2(dset, kind_attr_name, ktype, dspace, H5P_DEFAULT, H5P_DEFAULT);
  CHECK_H5(L, kattr, "Acreate(kind) in write_attrs()");
  CHECK_H5(L, H5Awrite(kattr, ktype, kind), "Awrite(kind) in write_attrs()");
  CHECK_H5(L, H5Aclose(kattr), "Aclose(kind) in write_attrs()");
  CHECK_H5(L, H5Tclose(ktype), "Tclose(kind) in write_attrs()");

  hid_t stype = get_sha256_type(L, b, 1);
  hid_t sattr = H5Acreate2(dset, csum_attr_name, stype, dspace, H5P_DEFAULT, H5P_DEFAULT);
  CHECK_H5(L, sattr, "Acreate(sum) in write_attrs()");
  CHECK_H5(L, H5Awrite(sattr, stype, sum->v), "Awrite(sum) in write_attrs()");
  CHECK_H5(L, H5Aclose(sattr), "Aclose(sum) in write_attrs()");
  CHECK_H5(L, H5Tclose(stype), "Tclose(sum) in write_attrs()");

  hid_t tftype = get_time_type(L, b, 1);
  hid_t tmtype = get_time_type(L, b, 0);
  long long now = (long long)(1e6 * qlua_timeofday());
  hid_t tattr = H5Acreate2(dset, time_attr_name, tftype, dspace, H5P_DEFAULT, H5P_DEFAULT);
  CHECK_H5(L, tattr, "Acreate(time) in write_attrs()");
  CHECK_H5(L, H5Awrite(tattr, tmtype, &now), "Awrite(time) in write_attrs()");
  CHECK_H5(L, H5Aclose(tattr), "Aclose(time) in write_attrs()");
  CHECK_H5(L, H5Tclose(tftype), "Tclose(time file type) in write_attrs()");
  CHECK_H5(L, H5Tclose(tmtype), "Tclose(time mem type) in write_attrs()");

  CHECK_H5(L, H5Sclose(dspace), "Sclose() in write_attrs()");
}

#if USE_Nc2
#define QNc  '2'
#define Qcolors "2"
#define Qs(a)   a ## 2
#define Qx(a,b)  a ## 2 ## b
#define QC(x)    2
#define QNC(x)
#include "hdf5_io-x.c"                                              /* DEPS */
#endif

#if USE_Nc3
#define QNc  '3'
#define Qcolors "3"
#define Qs(a)   a ## 3
#define Qx(a,b)  a ## 3 ## b
#define QC(x)    3
#define QNC(x)
#include "hdf5_io-x.c"                                              /* DEPS */
#endif

#if USE_NcN
#define QNc  'N'
#define Qcolors "N"
#define Qs(a)   a ## N
#define Qx(a,b)  a ## N ## b
#define QC(x)    (x)->nc
#define QNC(x)   (x), 
#include "hdf5_io-x.c"                                              /* DEPS */
#endif

/* writers */

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
  *kind = knString;
  *data = qlua_strdup(L, str);
  *filetype = get_string_type(L, b, len, 1);
  *memtype  = get_string_type(L, b, len, 0);
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
    break;
  case WS_Float:
    *data = qlua_malloc(L, sizeof (float));
    *(float *)(*data) = val;
    sha256_sum_add_floats(ctx, *data, 1);
    break;
  default:
    luaL_error(L, "Unknown precision in w_real()");
  }
  sha256_sum(sum, ctx);
  sha256_destroy(ctx);
  *kind = knReal;
  *filetype = get_real_type(L, b, opts->wsize, 1);
  *memtype  = get_real_type(L, b, opts->wsize, 0);
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
  } break;
  case WS_Float: {
    machine_complex_float *cv = qlua_malloc(L, sizeof (machine_complex_float));
    cv->re = QLA_real(*v);
    cv->im = QLA_imag(*v);
    *data = cv;
    sha256_sum_add_floats(ctx, *data, 2);
  } break;
  default:
    luaL_error(L, "Unknown precision in w_complex()");
  }
  sha256_sum(sum, ctx);
  sha256_destroy(ctx);
  *kind = knComplex;
  *filetype = get_complex_type(L, b, opts->wsize, 1);
  *memtype  = get_complex_type(L, b, opts->wsize, 0);
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
  *filetype = get_vecint_type(L, b, v->size, 1);
  *memtype =  get_vecint_type(L, b, v->size, 0);
  *kind = knVectorInt;
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
  *filetype = get_vecreal_type(L, b, v->size, opts->wsize, 1);
  *memtype  = get_vecreal_type(L, b, v->size, opts->wsize, 0);
  *kind = knVectorReal;
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
  *filetype = get_veccomplex_type(L, b, v->size, opts->wsize, 1);
  *memtype  = get_veccomplex_type(L, b, v->size, opts->wsize, 0);
  *kind = knVectorReal;
}

static void
w_matreal(lua_State *L, mHdf5File *b, mLattice *S,
          struct wopts_s *opts, struct laddr_s *laddr,
          SHA256_Sum *sum, void **data, hid_t *filetype, hid_t *memtype,
          const char **kind)
{
  mMatReal *v = qlua_checkMatReal(L, 3);
  SHA256_Context *ctx = sha256_create(L);
  int i, j;

  switch (opts->wsize) {
  case WS_Double: {
    double *mat = qlua_malloc(L, v->l_size * v->r_size * sizeof (double));
    double *ptr;
    for (ptr = mat, i = 0; i < v->l_size; i++) {
      for (j = 0; j < v->r_size; j++) {
        *ptr++ = gsl_matrix_get(v->m, i, j);
      }
    }
    *data = mat;
    sha256_sum_add_doubles(ctx, *data, v->l_size * v->r_size);
  } break;
  case WS_Float: {
    float *mat= qlua_malloc(L, v->l_size * v->r_size * sizeof (float));
    float *ptr;
    for (ptr = mat, i = 0; i < v->l_size; i++) {
      for (j = 0; j < v->r_size; j++) {
        *ptr++ = gsl_matrix_get(v->m, i, j);
      }
    }
    *data = mat;
    sha256_sum_add_floats(ctx, *data, v->l_size * v->r_size);
  } break;
  default:
    luaL_error(L, "Unknown precision in w_vecreal()");
  }
  sha256_sum(sum, ctx);
  sha256_destroy(ctx);
  *filetype = get_matreal_type(L, b, v->l_size, v->r_size, opts->wsize, 1);
  *memtype  = get_matreal_type(L, b, v->l_size, v->r_size, opts->wsize, 0);
  *kind = knMatrixReal;
}

static void
w_matcomplex(lua_State *L, mHdf5File *b, mLattice *S,
             struct wopts_s *opts, struct laddr_s *laddr,
             SHA256_Sum *sum, void **data, hid_t *filetype, hid_t *memtype,
             const char **kind)
{
  mMatComplex *v = qlua_checkMatComplex(L, 3);
  SHA256_Context *ctx = sha256_create(L);
  int i, j;

  switch (opts->wsize) {
  case WS_Double: {
    machine_complex_double *mat = qlua_malloc(L, v->l_size * v->r_size * sizeof (machine_complex_double));
    machine_complex_double *ptr;
    for (ptr = mat, i = 0; i < v->l_size; i++) {
      for (j = 0; j < v->r_size; j++) {
        gsl_complex zz = gsl_matrix_complex_get(v->m, i, j);
        ptr->re = GSL_REAL(zz);
        ptr->im = GSL_IMAG(zz);
        ptr++;
      }
    }
    *data = mat;
    sha256_sum_add_doubles(ctx, *data, 2 * v->l_size * v->r_size);
  } break;
  case WS_Float: {
    machine_complex_float *mat = qlua_malloc(L, v->l_size * v->r_size * sizeof (machine_complex_float));
    machine_complex_float *ptr;
    for (ptr = mat, i = 0; i < v->l_size; i++) {
      for (j = 0; j < v->r_size; j++) {
        gsl_complex zz = gsl_matrix_complex_get(v->m, i, j);
        ptr->re = GSL_REAL(zz);
        ptr->im = GSL_IMAG(zz);
      }
    }
    *data = mat;
    sha256_sum_add_floats(ctx, *data, 2 * v->l_size * v->r_size);
  } break;
  default:
    luaL_error(L, "Unknown precision in w_veccomplex()");
  }
  sha256_sum(sum, ctx);
  sha256_destroy(ctx);
  *filetype = get_matcomplex_type(L, b, v->l_size, v->r_size, opts->wsize, 1);
  *memtype  = get_matcomplex_type(L, b, v->l_size, v->r_size, opts->wsize, 0);
  *kind = knMatrixComplex;
}

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

static void
w_latreal(lua_State *L, mHdf5File *b, mLattice *S,
         struct wopts_s *opts, struct laddr_s *laddr,
         SHA256_Sum *sum, void **data, hid_t *filetype, hid_t *memtype,
         const char **kind)
{
  int *local_x = qlua_malloc(L, laddr->rank * sizeof (int));
  SHA256_Context *ctx = sha256_create(L);
  int volume = laddr->volume;
  int rank = laddr->rank;
  int i;

  switch (opts->wsize) {
  case WS_Double: {
    double *ptr = qlua_malloc(L, volume * sizeof (double));
    CALL_QDP(L);
    mLatReal *m = qlua_checkLatReal(L, 3, S);
    QLA_Real *locked = QDP_expose_R(m->ptr);
    for (i = 0; i < volume; i++) {
      qdp2hdf5_addr(local_x, i, laddr);
      QLUA_ASSERT(QDP_node_number_L(S->lat, local_x) == QDP_this_node);
      ptr[i] = QLA_elem_R(locked[QDP_index_L(S->lat, local_x)]);
      sha256_reset(ctx);
      sha256_sum_add_ints(ctx, &rank, 1);
      sha256_sum_add_ints(ctx, local_x, rank);
      sha256_sum_add_doubles(ctx, &ptr[i], 1);
      SHA256_Sum l_sum;
      sha256_sum(&l_sum, ctx);
      local_combine_checksums(sum, &l_sum);
    }
    QDP_reset_R(m->ptr);
    *data = ptr;
  } break;
  case WS_Float: {
    float *ptr = qlua_malloc(L, volume * sizeof (float));
    CALL_QDP(L);
    mLatReal *m = qlua_checkLatReal(L, 3, S);
    QLA_Real *locked = QDP_expose_R(m->ptr);
    for (i = 0; i < volume; i++) {
      qdp2hdf5_addr(local_x, i, laddr);
      QLUA_ASSERT(QDP_node_number_L(S->lat, local_x) == QDP_this_node);
      ptr[i] = QLA_elem_R(locked[QDP_index_L(S->lat, local_x)]);
      sha256_reset(ctx);
      sha256_sum_add_ints(ctx, &rank, 1);
      sha256_sum_add_ints(ctx, local_x, rank);
      sha256_sum_add_floats(ctx, &ptr[i], 1);
      SHA256_Sum l_sum;
      sha256_sum(&l_sum, ctx);
      local_combine_checksums(sum, &l_sum);
    }
    QDP_reset_R(m->ptr);
    *data = ptr;
  } break;
  default:
    luaL_error(L, "Unknown precision in w_latreal()");
  }
  sha256_destroy(ctx);
  qlua_free(L, local_x);
  *kind = knLatticeReal;
  *filetype = get_real_type(L, b, opts->wsize, 1);
  *memtype  = get_real_type(L, b, opts->wsize, 0);

}

static void
w_latcomplex(lua_State *L, mHdf5File *b, mLattice *S,
             struct wopts_s *opts, struct laddr_s *laddr,
             SHA256_Sum *sum, void **data, hid_t *filetype, hid_t *memtype,
             const char **kind)
{
  int *local_x = qlua_malloc(L, laddr->rank * sizeof (int));
  SHA256_Context *ctx = sha256_create(L);
  int volume = laddr->volume;
  int rank = laddr->rank;
  int i;

  switch (opts->wsize) {
  case WS_Double: {
    machine_complex_double *ptr = qlua_malloc(L, volume * sizeof (machine_complex_double));
    CALL_QDP(L);
    mLatComplex *m = qlua_checkLatComplex(L, 3, S);
    QLA_D_Complex *locked = QDP_expose_C(m->ptr);
    for (i = 0; i < volume; i++) {
      qdp2hdf5_addr(local_x, i, laddr);
      QLUA_ASSERT(QDP_node_number_L(S->lat, local_x) == QDP_this_node);
      QLA_D_Complex zz = QLA_elem_R(locked[QDP_index_L(S->lat, local_x)]);
      ptr[i].re = QLA_real(zz);
      ptr[i].im = QLA_imag(zz);
      sha256_reset(ctx);
      sha256_sum_add_ints(ctx, &rank, 1);
      sha256_sum_add_ints(ctx, local_x, rank);
      sha256_sum_add_doubles(ctx, (double *)&ptr[i], 2);
      SHA256_Sum l_sum;
      sha256_sum(&l_sum, ctx);
      local_combine_checksums(sum, &l_sum);
    }
    QDP_reset_C(m->ptr);
    *data = ptr;
  } break;
  case WS_Float: {
    machine_complex_float *ptr = qlua_malloc(L, volume * sizeof (machine_complex_float));
    CALL_QDP(L);
    mLatComplex *m = qlua_checkLatComplex(L, 3, S);
    QLA_D_Complex *locked = QDP_expose_C(m->ptr);
    for (i = 0; i < volume; i++) {
      qdp2hdf5_addr(local_x, i, laddr);
      QLUA_ASSERT(QDP_node_number_L(S->lat, local_x) == QDP_this_node);
      QLA_D_Complex zz = QLA_elem_R(locked[QDP_index_L(S->lat, local_x)]);
      ptr[i].re = QLA_real(zz);
      ptr[i].im = QLA_imag(zz);
      sha256_reset(ctx);
      sha256_sum_add_ints(ctx, &rank, 1);
      sha256_sum_add_ints(ctx, local_x, rank);
      sha256_sum_add_floats(ctx, (float *)&ptr[i], 2);
      SHA256_Sum l_sum;
      sha256_sum(&l_sum, ctx);
      local_combine_checksums(sum, &l_sum);
    }
    QDP_reset_C(m->ptr);
    *data = ptr;
  } break;
  default:
    luaL_error(L, "Unknown precision in w_latreal()");
  }
  sha256_destroy(ctx);
  qlua_free(L, local_x);
  *kind = knLatticeReal;
  *filetype = get_complex_type(L, b, opts->wsize, 1);
  *memtype  = get_complex_type(L, b, opts->wsize, 0);

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

static struct {
  QLUA_Type qtype;
  int is_parallel;
  OutPacker_H5 repack;
} qowtable[] = {
  { qString,                 0,  w_string        },
  { qReal,                   0,  w_real          },
  { qComplex,                0,  w_complex       },
  { qVecInt,                 0,  w_vecint        },
  { qVecReal,                0,  w_vecreal       },
  { qVecComplex,             0,  w_veccomplex    },
  { qMatReal,                0,  w_matreal       },
  { qMatComplex,             0,  w_matcomplex    },
  { qLatInt,                 1,  w_latint        },
  { qLatReal,                1,  w_latreal       },
  { qLatComplex,             1,  w_latcomplex    },
#if USE_Nc2
  { qSeqColVec2,             0,  w_colvec2       },
  { qSeqColMat2,             0,  w_colmat2       },
  { qSeqDirFerm2,            0,  w_dirferm2      },
  { qSeqDirProp2,            0,  w_dirprop2      },
  { qLatColVec2,             1,  w_latcolvec2    },
  { qLatColMat2,             1,  w_latcolmat2    },
  { qLatDirFerm2,            1,  w_latdirferm2   },
  { qLatDirProp2,            1,  w_latdirprop2   },
#endif
#if USE_Nc3
  { qSeqColVec3,             0,  w_colvec3       },
  { qSeqColMat3,             0,  w_colmat3       },
  { qSeqDirFerm3,            0,  w_dirferm3      },
  { qSeqDirProp3,            0,  w_dirprop3      },
  { qLatColVec3,             1,  w_latcolvec3    },
  { qLatColMat3,             1,  w_latcolmat3    },
  { qLatDirFerm3,            1,  w_latdirferm3   },
  { qLatDirProp3,            1,  w_latdirprop3   },
#endif
#if USE_NcN
  { qSeqColVecN,             0,  w_colvecN       },
  { qSeqColMatN,             0,  w_colmatN       },
  { qSeqDirFermN,            0,  w_dirfermN      },
  { qSeqDirPropN,            0,  w_dirpropN      },
  { qLatColVecN,             1,  w_latcolvecN    },
  { qLatColMatN,             1,  w_latcolmatN    },
  { qLatDirFermN,            1,  w_latdirfermN   },
  { qLatDirPropN,            1,  w_latdirpropN   },
#endif
  { qNoType,                 0,  NULL            }
};

static int
qhdf5_write(lua_State *L)
{
  mHdf5File *b = qlua_checkHdf5Writer(L, 1);
  const char *p = luaL_checkstring(L, 2);
  QLUA_Type kind = qlua_qtype(L, 3);
  int count;
  int i;

  check_writer(L, b);
  for (i = 0; qowtable[i].qtype != qNoType; i++) {
    if (qowtable[i].qtype == kind)
      break;
  }
  if (qowtable[i].repack == NULL)
    luaL_error(L, "unwritable data");
  count = (qowtable[i].is_parallel? write_lat: write_seq)(L, b, p, qowtable[i].repack);
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

static int
q_hdf5_reader(lua_State *L)
{
  const char *name = luaL_checkstring(L, 1);
  mHdf5File *w = qlua_newHdf5Reader(L);
  MPI_Comm comm = MPI_COMM_WORLD;
  MPI_Info info = MPI_INFO_NULL;
  hid_t acc_tpl1;

  qlua_Hdf5_enter(L);

  acc_tpl1 = H5Pcreate(H5P_FILE_ACCESS);
  CHECK_H5(L, acc_tpl1, "Pcreate() failed");
  CHECK_H5(L, H5Pset_fapl_mpio(acc_tpl1, comm, info), "Pset_fapl_mpio() failed");

  w->file = H5Fopen(name, H5F_ACC_RDONLY, acc_tpl1);
  CHECK_H5(L, w->file, "qcd.hdf5.Reader failed");
  CHECK_H5(L, H5Pclose(acc_tpl1), "Pclose(template) failed");
  w->cwd = H5Gopen2(w->file, "/", H5P_DEFAULT);
  CHECK_H5(L, w->cwd, "Gopen2(\"/\") failed");

  qlua_Hdf5_leave();

  return 1;
}

/* metatables */
static const struct luaL_Reg mtFile[] = {
  { "__tostring",       qhdf5_fmt     },
  { "__gc",             qhdf5_gc      },
  { "close",            qhdf5_close   },
  { "write",            qhdf5_write   },
  { "chpath",           qhdf5_chpath  },
  { "mkpath",           qhdf5_mkpath  },
  { "list",             qhdf5_list    },
  { "flush",            qhdf5_flush   },
  { "stat",             qhdf5_stat    },
  { "cwd",              qhdf5_cwd     },
  { "remove",           qhdf5_remove  },
#if 0 /* XXXXXXXXXXXXXXXXXXXXXXXXXX */
  { "read",             qhdf5_real    },
#endif /* XXXXXXXXXXXXXXXXXXXXXXXXXX */
  { NULL,               NULL          }
};

/* names and routines for qcd.hdf5 table */
static const struct luaL_Reg fHDF5io[] = {
  { "Reader",   q_hdf5_reader },
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
