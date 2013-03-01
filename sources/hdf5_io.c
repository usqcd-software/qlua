#include "qlua.h"                                                    /* DEPS */
#include "qvector.h"                                                 /* DEPS */
#include "qmatrix.h"                                                 /* DEPS */
#include "hdf5_io.h"                                                 /* DEPS */
#include <string.h>
#include <complex.h>
#include <sys/stat.h>
#include "hdf5.h"
#include "qmp.h"

#include <assert.h> /* XXX */
#define XXX_assert(x) assert(x)
const char hdf5_io[] = "hdf5";
static const char mtnReader[] = "qcd.hdf5.mtReader";
static const char mtnWriter[] = "qcd.hdf5.mtRriter";
#if 0 /* XXX */

struct Hdf5Reader_s {}; /* XXX */
#endif /* XXX */

/* type names for H5 types */
enum {
  qhCharType,
  qhIntType,
  qhRealType,
  qhComplexType,
  qhTypeCount
};

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
    H5Tinsert(v, "re", HOFFSET(qlua_machine_complex, re), H5T_NATIVE_DOUBLE);    
    if (ec < 0) luaL_error(L, "HDF5: complex.re insert failed");
    H5Tinsert(v, "im", HOFFSET(qlua_machine_complex, im), H5T_NATIVE_DOUBLE);    
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
    ec = H5Tcommit(file, ".complex", v, H5P_DEFAULT,  H5P_DEFAULT,  H5P_DEFAULT);
    if (ec < 0) luaL_error(L, "HDF5: complex commit failed");
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

#if 0 /* XXX */
/* helpers */
static void
check_reader(lua_State *L, mHdf5Reader *r)
{
  XXX_assert("check_reader" == "implemented"); /* XXX */
    if (r->ptr == 0)
        luaL_error(L, "closed hdf5 reader");
}
#endif /* XXX */

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
    dnext = H5Gopen2(dhandle, dname, H5P_DEFAULT);
    if (dnext < 0)
      luaL_error(L, "chpath failed");
    if (isused == 0)
      H5Gclose(dhandle);
    dhandle = dnext;
    isused = 0;
    dname = buf;
  }
  free(dpath);
  H5Gclose(cwd);
  return dhandle;
}

#if 0 /* XXX */
/* allocation */
static mHdf5Reader *
qlua_newHdf5Reader(lua_State *L, struct Hdf5Reader_s *reader)
{
  XXX_assert("qlua_newHdf5Reader" == "implemented"); /* XXX */
  return 0;
    mHdf5Reader *h = lua_newuserdata(L, sizeof (mHdf5Reader));

    h->ptr = reader;
    h->dir = hdf5_reader_root(reader);
    luaL_getmetatable(L, mtnReader);
    lua_setmetatable(L, -2);

    return h;
}
#endif /* XXX */

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

#if 0 /* XXX */
/* type checking */
mHdf5Reader *
qlua_checkHdf5Reader(lua_State *L, int idx)
{
  XXX_assert("qlua_checkHdf5Reader" == "implemented"); /* XXX */
  return 0;

    void *v = luaL_checkudata(L, idx, mtnReader);

    luaL_argcheck(L, v != 0, idx, "qcd.hdf5.Reader expected");

    return v;
}
#endif /* XXX */

mHdf5Writer *
qlua_checkHdf5Writer(lua_State *L, int idx)
{
  void *v = luaL_checkudata(L, idx, mtnWriter);
  
  luaL_argcheck(L, v != 0, idx, "qcd.hdf5.Writer expected");
  return v;
}

#if 0 /* XXX */
/* conversion to string */
static int
qhdf5_r_fmt(lua_State *L)
{
  XXX_assert("qhdf5_r_fmt" == "implemented"); /* XXX */
  return 0;

    mHdf5Reader *b = qlua_checkHdf5Reader(L, 1);
    char fmt[72];
    
    if ( b->ptr )
        sprintf(fmt, "hdf5.Reader(%p)", b->ptr);
    else
        sprintf(fmt, "hdf5.Reader(closed)");

    lua_pushstring(L, fmt);

    return 1;
}
#endif /* XXX */

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

#if 0 /* XXX */
/* garbage collection */
static int
qhdf5_r_gc(lua_State *L)
{
  XXX_assert("qhdf5_r_gc" == "implemented"); /* XXX */
  return 0;

    mHdf5Reader *b = qlua_checkHdf5Reader(L, 1);


    qlua_Hdf5_enter(L);
    if (b->ptr)
        hdf5_reader_close(b->ptr);
    b->ptr = 0;
    qlua_Hdf5_leave();
    
    return 0;
}
#endif /* XXX */

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

#if 0 /* XXX */
/* closing */
static int
qhdf5_r_close(lua_State *L)
{
  XXX_assert("qhdf5_r_close" == "implemented"); /* XXX */
  return 0;

    mHdf5Reader *b = qlua_checkHdf5Reader(L, 1);

    check_reader(L, b);
    qlua_Hdf5_enter(L);
    hdf5_reader_close(b->ptr);
    b->ptr = 0;
    qlua_Hdf5_leave();
    
    return 0;
}
#endif /* XXX */

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

#if 0 /* XXX */
static int
qhdf5_r_status(lua_State *L)
{
  XXX_assert("qhdf5_r_status" == "implemented"); /* XXX */
  return 0;
    mHdf5Reader *b = qlua_checkHdf5Reader(L, 1);
    const char *s;

    check_reader(L, b);

    qlua_Hdf5_enter(L);
    s = hdf5_reader_errstr(b->ptr);
    qlua_Hdf5_leave();
        
    if (s == 0)
        lua_pushnil(L);
    else
        lua_pushstring(L, s);

    return 1;
}
#endif /* XXX */

static int
qhdf5_w_status(lua_State *L)
{
  mHdf5Writer *b = qlua_checkHdf5Writer(L, 1);

  lua_pushnil(L);
  return 1;
}

#if 0 /* XXX */
typedef struct {
    lua_State *L;
    int        k;
} qHdf5Dir;

static void
qar_get_list(struct Hdf5Node_s *n, void *arg)
{
  XXX_assert("qar_get_list" == "implemented"); /* XXX */
    qHdf5Dir *d = arg;

    lua_pushstring(d->L, hdf5_symbol_name(hdf5_node_name(n)));
    lua_rawseti(d->L, -2, d->k);
    d->k++;
}

static int
qhdf5_r_list(lua_State *L)
{
  XXX_assert("qhdf5_r_list" == "implemented"); /* XXX */
  return 0;
    mHdf5Reader *b = qlua_checkHdf5Reader(L, 1);
    struct Hdf5Node_s *r;
    const char *p = luaL_checkstring(L, 2);
    qHdf5Dir dir;
    
    check_reader(L, b);

    qlua_Hdf5_enter(L);
    r = qlua_Hdf5ReaderChPath(b, p);
    if (r == 0)
        return luaL_error(L, hdf5_reader_errstr(b->ptr));
    lua_newtable(L);
    dir.L = L;
    dir.k = 1;
    hdf5_node_foreach(r, qar_get_list, &dir);
    qlua_Hdf5_leave();

    return 1;
}

static int
qhdf5_r_read(lua_State *L)
{
  XXX_assert("qhdf5_r_read" == "implemented"); /* XXX */
  return 0;
    mHdf5Reader *b = qlua_checkHdf5Reader(L, 1);
    const char *p = luaL_checkstring(L, 2);
    struct Hdf5Node_s *n;
    uint32_t size;

    check_reader(L, b);

    qlua_Hdf5_enter(L);
    n = qlua_Hdf5ReaderChPath(b, p);
    if (n == 0)
        goto end;
    size = hdf5_node_size(n);

    switch (hdf5_node_type(n)) {
    case hdf5NodeVoid:
        lua_pushboolean(L, 1);
        qlua_Hdf5_leave();
        return 1;
    case hdf5NodeChar: {
                char *d = qlua_malloc(L, size + 1);

        if (hdf5_node_get_char(b->ptr, n, d, size) == 0) {
            d[size] = 0;
            lua_pushstring(L, d);
            qlua_Hdf5_leave();
                        qlua_free(L, d);

            return 1;
        }
                qlua_free(L, d);
        break;
    }
    case hdf5NodeInt: {
                uint32_t *d = qlua_malloc(L, size * sizeof (uint32_t));

        if (hdf5_node_get_int(b->ptr, n, d, size) == 0) {
            mVecInt *v = qlua_newVecInt(L, size);
            int i;

            for (i = 0; i < size; i++)
                v->val[i] = d[i];
            qlua_Hdf5_leave();
                        qlua_free(L, d);

            return 1;
        }
                qlua_free(L, d);
        break;
    }
    case hdf5NodeDouble: {
                double *d = qlua_malloc(L, size * sizeof (double));

        if (hdf5_node_get_double(b->ptr, n, d, size) == 0) {
            mVecReal *v = qlua_newVecReal(L, size);
            int i;

            for (i = 0; i < size; i++)
                v->val[i] = d[i];
            qlua_Hdf5_leave();
                        qlua_free(L, d);
                        
            return 1;
        }
                qlua_free(L, d);
        break;
    }
    case hdf5NodeComplex: {
                double _Complex *d = qlua_malloc(L, size * sizeof (double _Complex));

        if (hdf5_node_get_complex(b->ptr, n, d, size) == 0) {
            mVecComplex *v = qlua_newVecComplex(L, size);
            int i;

            for (i = 0; i < size; i++) {
                QLA_real(v->val[i]) = creal(d[i]);
                QLA_imag(v->val[i]) = cimag(d[i]);
            }
            qlua_Hdf5_leave();
                        qlua_free(L, d);

            return 1;
        }
                qlua_free(L, d);
        break;
    }
    default:
        goto end;
    }
    qlua_Hdf5_leave();
    return luaL_error(L, hdf5_reader_errstr(b->ptr));

end:
    qlua_Hdf5_leave();
    return luaL_error(L, "bad arguments for HDF5 read");
}
#endif /* XXX */

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
      dtype_file = get_qh_file_type(L, b->file, b->ft, qhCharType);
      dtype_machine = get_qh_machine_type(L, b->file, b->mt, qhCharType);
      dataspace = H5Screate_simple(1, &len, NULL);
    } break;
    case qVecInt: {
      mVecInt *v = qlua_checkVecInt(L, 3);
      hsize_t len = v->size;
      int *ptr;
      int i;
      ptr = qlua_malloc(L, len * sizeof (int));
      data = ptr;
      dtype_file = get_qh_file_type(L, b->file, b->ft, qhIntType);
      dtype_machine = get_qh_machine_type(L, b->file, b->mt, qhIntType);
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
      dtype_file = get_qh_file_type(L, b->file, b->ft, qhRealType);
      dtype_machine = get_qh_machine_type(L, b->file, b->mt, qhRealType);
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
      dtype_file = get_qh_file_type(L, b->file, b->ft, qhComplexType);
      dtype_machine = get_qh_machine_type(L, b->file, b->mt, qhComplexType);
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
      dtype_file = get_qh_file_type(L, b->file, b->ft, qhRealType);
      dtype_machine = get_qh_machine_type(L, b->file, b->mt, qhRealType);
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
      dtype_file = get_qh_file_type(L, b->file, b->ft, qhComplexType);
      dtype_machine = get_qh_machine_type(L, b->file, b->mt, qhComplexType);
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
    dataset = H5Dcreate2(wdir, ename, dtype_file, dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    status = H5Dwrite(dataset, dtype_machine, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);
    if (status < 0)
      luaL_error(L, "write error");

    H5Dclose(dataset);
    H5Sclose(dataspace);
    free(data);
    free(dpath);
    if (wdir != b->cwd)
      H5Gclose(wdir);
  }
  qlua_Hdf5_leave();
  return 0;
}

#if 0 /* XXX */
static int
qhdf5_r_chpath(lua_State *L)
{
  XXX_assert("qhdf5_r_chpath" == "implemented"); /* XXX */
  return 0;
    mHdf5Reader *b = qlua_checkHdf5Reader(L, 1);
    const char *p = luaL_checkstring(L, 2);
    struct Hdf5Node_s *r;

    check_reader(L, b);

    qlua_Hdf5_enter(L);
    r = qlua_Hdf5ReaderChPath(b, p);
    b->dir = r;
    qlua_Hdf5_leave();

    if (r != NULL)
        return 0;
    else
        return luaL_error(L, hdf5_reader_errstr(b->ptr));
}
#endif /* XXX */

static int
qhdf5_w_chpath(lua_State *L)
{
  mHdf5Writer *b = qlua_checkHdf5Writer(L, 1);
  const char *p = luaL_checkstring(L, 2);
  
  check_writer(L, b);
  
  qlua_Hdf5_enter(L);
  if (b->master) {
    hid_t dir = change_hdf5_path(L, b->file, b->cwd, p);
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
    free(dpath);
    if (dh != b->cwd)
      H5Gclose(dh);
  }
  qlua_Hdf5_leave();
  return 0;
}

#if 0 /* XXX */
static int
q_hdf5_reader(lua_State *L)
{
  XXX_assert("q_hdf5_reader" == "implemented"); /* XXX */
  return 0;
    const char *name = luaL_checkstring(L, 1);
    struct Hdf5Reader_s *r;
    const char *msg;

    qlua_Hdf5_enter(L);
    r = hdf5_reader(name);
    msg = hdf5_reader_errstr(r);
    if (msg == NULL) {
        qlua_newHdf5Reader(L, r);
        qlua_Hdf5_leave();

        return 1;
    } else {
        hdf5_reader_close(r);
        qlua_Hdf5_leave();

        return luaL_error(L, msg);
    }
}
#endif /* XXX */

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
#if 0 /* XXX */
    { "__tostring",       qhdf5_r_fmt},
    { "__gc",             qhdf5_r_gc},
    { "close",            qhdf5_r_close},
    { "read",             qhdf5_r_read},
    { "chpath",           qhdf5_r_chpath },
    { "status",           qhdf5_r_status },
    { "list",             qhdf5_r_list },
#endif /* XXX */
    { NULL,               NULL}
};

static const struct luaL_Reg mtWriter[] = {
  { "__tostring",       qhdf5_w_fmt },
  { "__gc",             qhdf5_w_gc },
  { "close",            qhdf5_w_close },
  { "write",            qhdf5_w_write },
  { "chpath",           qhdf5_w_chpath },
  { "mkpath",           qhdf5_w_mkpath },
  { "status",           qhdf5_w_status },
  { NULL,               NULL}
};

/* names and routines for qcd.qdpc table */
static const struct luaL_Reg fHDF5io[] = {
#if 0  /* XXX */
  { "Reader",   q_hdf5_reader},
#endif /* XXX */
  { "Writer",   q_hdf5_writer},
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
