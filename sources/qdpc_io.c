#include "modules.h"                                                 /* DEPS */
#include "qlua.h"                                                    /* DEPS */
#include "lattice.h"                                                 /* DEPS */
#include "qdpc_io.h"                                                 /* DEPS */
#include "qio_utils.h"                                               /* DEPS */
#include "latcolmat.h"                                               /* DEPS */
#include "latcolvec.h"                                               /* DEPS */
#include "latcomplex.h"                                              /* DEPS */
#include "latdirferm.h"                                              /* DEPS */
#include "latdirprop.h"                                              /* DEPS */
#include "latint.h"                                                  /* DEPS */
#include "latrandom.h"                                               /* DEPS */
#include "latreal.h"                                                 /* DEPS */
#include "qvector.h"                                                 /* DEPS */
#ifdef HAS_GSL
#include "qmatrix.h"                                                 /* DEPS */
#endif

#include "qdp_df.h"
#if USE_Nc2
#include "qdp_f2.h"
#include "qdp_df2.h"
#endif
#if USE_Nc3
#include "qdp_f3.h"
#include "qdp_df3.h"
#endif
#if USE_NcN
#include "qdp_fn.h"
#include "qdp_dfn.h"
#endif
#include <string.h>

typedef struct {
    QDP_Reader *ptr;
} mReader;

typedef struct {
    QDP_Writer *ptr;
} mWriter;

/* this magic used to identify global matrices in QIO */
static const char qio_matrix_magic[] = "480C94A9-F941-4FD5-9882-3E70B004AACF";

static const char qdp_io[] = "qdpc";
static const char ReaderName[] = "qcd.qdpc.reader";
static const char WriterName[] = "qcd.qdpc.writer";

static mReader *q_checkReader(lua_State *L, int idx, mLattice *S);
static mWriter *q_checkWriter(lua_State *L, int idx, mLattice *S);

/* check for a live QIO object */
static void check_reader(lua_State *L, mReader *b);
static void check_writer(lua_State *L, mWriter *b);


/* plain lattice types */
#define T_QTYPE       QDP_Int
#define QLUA_NAME(x)  x ## LatInt
#define X_ID(x)       x ## I
#undef X_DF
#include "qdpc_io-x.c"                                               /* DEPS */

#define T_QTYPE       QDP_RandomState
#define QLUA_NAME(x)  x ## LatRandom
#define X_ID(x)       x ## S
#undef X_DF
#include "qdpc_io-x.c"                                               /* DEPS */

#define T_QTYPE       QDP_D_Real
#define QLUA_NAME(x)  x ## LatReal
#define X_ID(x)       x ## R
#define X_ID2(a,b)    a ## R ## b
#define X_ID3(a,b)    a ## R ## b ## R
#define T_sTYPE       QDP_F_Real
#define X_DF
#include "qdpc_io-x.c"                                               /* DEPS */

#define T_QTYPE       QDP_D_Complex
#define QLUA_NAME(x)  x ## LatComplex
#define X_ID(x)       x ## C
#define X_ID2(a,b)    a ## C ## b
#define X_ID3(a,b)    a ## C ## b ## C
#define T_sTYPE       QDP_F_Complex
#define X_DF
#include "qdpc_io-x.c"                                               /* DEPS */

#if USE_Nc2 || USE_Nc3 || USE_NcN
/* colored lattice types */
#define QT(a)         a ## ColorVector
#define QTx(a,b)      a ## ColorVector ## b
#define QN(a,b)       a ## LatColVec ## b
#define QA(p)         p ## V
#define QAx(a,b)      a ## V ## b
#define QAy(a,b)      a ## V ## b ## V
#include "qdpc_io-y.c"                                               /* DEPS */

#define QT(a)         a ## ColorMatrix
#define QTx(a,b)      a ## ColorMatrix ## b
#define QN(a,b)       a ## LatColMat ## b
#define QA(p)         p ## M
#define QAx(a,b)      a ## M ## b
#define QAy(a,b)      a ## M ## b ## M
#include "qdpc_io-y.c"                                               /* DEPS */

#define QT(a)         a ## DiracFermion
#define QTx(a,b)      a ## DiracFermion ## b
#define QN(a,b)       a ## LatDirFerm ## b
#define QA(p)         p ## D
#define QAx(a,b)      a ## D ## b
#define QAy(a,b)      a ## D ## b ## D
#include "qdpc_io-y.c"                                               /* DEPS */

#define QT(a)         a ## DiracPropagator
#define QTx(a,b)      a ## DiracPropagator ## b
#define QN(a,b)       a ## LatDirProp ## b
#define QA(p)         p ## P
#define QAx(a,b)      a ## P ## b
#define QAy(a,b)      a ## P ## b ## P
#include "qdpc_io-y.c"                                               /* DEPS */
#endif /* use any colors */

/* QDPC reader */
static int
qdpc_r_fmt(lua_State *L)
{
    mReader *b = q_checkReader(L, 1, NULL);
    char fmt[72];
    
    if ( b->ptr )
        sprintf(fmt, "qdpc.Reader(%p)", b->ptr);
    else
        sprintf(fmt, "qdpc.Reader(closed)");

    lua_pushstring(L, fmt);

    return 1;
}

static int
qdpc_r_gc(lua_State *L)
{
    mReader *b = q_checkReader(L, 1, NULL);

    if (b->ptr)
        QDP_close_read(b->ptr);
    b->ptr = 0;
    
    return 0;
}

static void
check_reader(lua_State *L, mReader *b)
{
    if (b->ptr == 0)
        luaL_error(L, "closed QDP/C reader");
}

static int
qdpc_r_skip(lua_State *L)
{
    mReader *reader = q_checkReader(L, 1, NULL);
    int status;
    int rcount;

    check_reader(L, reader);

    CALL_QDP(L);
    status = QDP_next_record(reader->ptr);
    if (status == 0) {
        lua_pushboolean(L, 1);
        rcount = 1;
    } else {
        rcount = 0;
    }

    return rcount;
}

static int
qdpc_r_close(lua_State *L)
{
    mReader *b = q_checkReader(L, 1, NULL);
    
    check_reader(L, b);

    QDP_close_read(b->ptr);
    b->ptr = 0;
    
    return 0;
}

static int
qdpc_r_iv(lua_State *L)
{
    mReader *reader = q_checkReader(L, 1, NULL);
    int len = luaL_checkint(L, 2);
    mVecInt *v = qlua_newVecInt(L, len);
    QDP_String *info;
    
    check_reader(L, reader);

    CALL_QDP(L);

    info = QDP_string_create();
    if (len <= 0)
        goto error;

    if (QDP_vread_i(reader->ptr, info, v->val, v->size) != 0)
        goto error;

    lua_pushstring(L, QDP_string_ptr(info));
    QDP_string_destroy(info);
    return 2;
error:
    QDP_string_destroy(info);
    return luaL_error(L, "qdpc read error");
}

static int
qdpc_r_rv(lua_State *L)
{
    mReader *reader = q_checkReader(L, 1, NULL);
    int len = luaL_checkint(L, 2);
    mVecReal *v = qlua_newVecReal(L, len);
    QDP_String *info;
    
    check_reader(L, reader);

    CALL_QDP(L);

    info = QDP_string_create();
    if (len <= 0)
        goto error;

    if (QDP_D_vread_r(reader->ptr, info, v->val, v->size) != 0)
        goto error;

    lua_pushstring(L, QDP_string_ptr(info));
    QDP_string_destroy(info);
    return 2;
error:
    QDP_string_destroy(info);
    return luaL_error(L, "qdpc read error");
}

static int
qdpc_r_cv(lua_State *L)
{
    mReader *reader = q_checkReader(L, 1, NULL);
    int len = luaL_checkint(L, 2);
    mVecComplex *v = qlua_newVecComplex(L, len);
    QDP_String *info;
    
    check_reader(L, reader);

    CALL_QDP(L);

    info = QDP_string_create();
    if (len <= 0)
        goto error;

    if (QDP_D_vread_c(reader->ptr, info, v->val, v->size) != 0)
        goto error;

    lua_pushstring(L, QDP_string_ptr(info));
    QDP_string_destroy(info);
    return 2;
error:
    QDP_string_destroy(info);
    return luaL_error(L, "qdpc read error");
}

#ifdef HAS_GSL
static int
qdpc_r_rm(lua_State *L)
{
    mReader *reader = q_checkReader(L, 1, NULL);
    int l_size = luaL_checkint(L, 2);
    int r_size = luaL_checkint(L, 3);
    mMatReal *v = qlua_newMatReal(L, l_size, r_size);
    QLA_Int f_size[2];
    QDP_String *info;
    int i, j;

    check_reader(L, reader);

    CALL_QDP(L);

    info = QDP_string_create();
    if ((l_size <= 0) || (r_size <= 0))
        goto error;

    if (QDP_vread_i(reader->ptr, info, f_size, 2) != 0)
        goto error;
    if ((f_size[0] != l_size) ||
        (f_size[1] != r_size))
        goto error;
    if (strcmp(QDP_string_ptr(info), qio_matrix_magic) != 0)
        goto error;
    {
        QLA_D_Real d[l_size * r_size];
        if (QDP_D_vread_r(reader->ptr, info, d, l_size * r_size) != 0) {
            goto error;
        }
        /* must only agree with the writer below */
        for (j = 0; j < r_size; j++) {
            for (i = 0; i < l_size; i++) {
                gsl_matrix_set(v->m, i, j, d[i + l_size * j]);
            }
        }
        lua_pushstring(L, QDP_string_ptr(info));
        QDP_string_destroy(info);
        return 2;
    }
error:
    QDP_string_destroy(info);
    return luaL_error(L, "qdpc read error");
}

static int
qdpc_r_cm(lua_State *L)
{
    mReader *reader = q_checkReader(L, 1, NULL);
    int l_size = luaL_checkint(L, 2);
    int r_size = luaL_checkint(L, 3);
    mMatComplex *v = qlua_newMatComplex(L, l_size, r_size);
    QLA_Int f_size[2];
    QDP_String *info;
    int i, j;

    check_reader(L, reader);

    CALL_QDP(L);

    info = QDP_string_create();
    if ((l_size <= 0) || (r_size <= 0))
        goto error;

    if (QDP_vread_i(reader->ptr, info, f_size, 2) != 0)
        goto error;
    if ((f_size[0] != l_size) ||
        (f_size[1] != r_size))
        goto error;
    if (strcmp(QDP_string_ptr(info), qio_matrix_magic) != 0)
        goto error;
    {
        QLA_D_Complex d[l_size * r_size];
        if (QDP_D_vread_c(reader->ptr, info, d, l_size * r_size) != 0) {
            goto error;
        }
        /* must only agree with the writer below */
        for (j = 0; j < r_size; j++) {
            for (i = 0; i < l_size; i++) {
                QLA_D_Complex z = d[i + l_size * j];
                gsl_complex zz;
                
                GSL_SET_COMPLEX(&zz, QLA_real(z), QLA_imag(z));
                gsl_matrix_complex_set(v->m, i, j, zz);
            }
        }
        lua_pushstring(L, QDP_string_ptr(info));
        QDP_string_destroy(info);
        return 2;
    }
error:
    QDP_string_destroy(info);
    return luaL_error(L, "qdpc read error");
}
#endif

static const struct luaL_Reg mtReader[] = {
    { "__tostring",       qdpc_r_fmt               },
    { "__gc",             qdpc_r_gc                },
    { "__newindex",       qlua_nowrite             },
    { "skip",             qdpc_r_skip              },
    { "close",            qdpc_r_close             },
    { "Int",              qdpc_r_I                 },
    { "RandomState",      qdpc_r_S                 },
    { "Real",             qdpc_r_R                 },
    { "Complex",          qdpc_r_C                 },
#if USE_Nc2 || USE_Nc3 || USE_NcN
    { "ColorVector",      qdpc_r_ColorVector       },
    { "ColorVectorN",     qdpc_r_ColorVectorN      },
    { "ColorMatrix",      qdpc_r_ColorMatrix       },
    { "ColorMatrixN",     qdpc_r_ColorMatrixN      },
    { "DiracFermion",     qdpc_r_DiracFermion      },
    { "DiracFermionN",    qdpc_r_DiracFermionN     },
    { "DiracPropagator",  qdpc_r_DiracPropagator   },
    { "DiracPropagatorN", qdpc_r_DiracPropagatorN  },
#endif /* use any colors */
    { "int_vector",       qdpc_r_iv                },
    { "real_vector",      qdpc_r_rv                },
    { "complex_vector",   qdpc_r_cv                },
#ifdef HAS_GSL
    { "real_matrix",      qdpc_r_rm                },
    { "complex_matrix",   qdpc_r_cm                },
#endif
    /* "lattice" */
    /* "a-type" */
    { NULL,               NULL                     }
};

static mReader *
q_checkReader(lua_State *L, int idx, mLattice *S)
{
    void *v = qlua_checkLatticeType(L, idx, qReader, ReaderName);
    
    if (S) {
        mLattice *S1 = qlua_ObjLattice(L, idx);
        if (S1->id != S->id)
            luaL_error(L, "%s on a wrong lattice", ReaderName);
        lua_pop(L, 1);
    }

    return (mReader *)v;
}

static mReader *
q_newReader(lua_State *L, QDP_Reader *reader)
{
    mReader *h = lua_newuserdata(L, sizeof (mReader));

    h->ptr = reader;
    qlua_createLatticeTable(L, 1, mtReader, qReader, ReaderName);
    lua_setmetatable(L, -2);

    return h;
}

/* QDP/C reader constructors
 *
 * qcd.qdpc.reader(L, file_name)
 */
static int
q_qdpc_reader (lua_State *L)
{
    mLattice *S = qlua_checkLattice(L, 1);
    const char *name = luaL_checkstring(L, 2);
    QDP_String *xml = 0;
    QDP_Reader *reader = 0;
    int rcount;

    /* first collect garbage */
    CALL_QDP(L);

    /* go through the motions of opening a QDP reader */
    xml = QDP_string_create();
    reader = QDP_open_read_L(S->lat, xml, (char *)name); /* [sic] */

    /* convert QDP results to LUA */
    if (reader == 0) {
        /* reader is not opened - return no results */
        rcount = 0;
    } else {
        /* we have a reader - return (mReader, file_xml) */
        q_newReader(L, reader);
        lua_pushstring(L, QDP_string_ptr(xml));
        rcount = 2;
    }
    /* free QDP_string now */
    QDP_string_destroy(xml);

    return rcount;
}

/* QDPC writer */

static void
check_writer(lua_State *L, mWriter *b)
{
    if (b->ptr == 0)
        luaL_error(L, "closed QDP/C writer");
}

static int
qdpc_w_fmt(lua_State *L)
{
    mWriter *b = q_checkWriter(L, 1, NULL);
    char fmt[72];
    
    if (b->ptr)
        sprintf(fmt, "qdpc.Writer(%p)", b->ptr);
    else
        sprintf(fmt, "qdpc.Writer(closed)");

    lua_pushstring(L, fmt);

    return 1;
}

static int
qdpc_w_gc(lua_State *L)
{
    mWriter *b = q_checkWriter(L, 1, NULL);

    if (b->ptr)
        QDP_close_write(b->ptr);
    b->ptr = 0;
    
    return 0;
}

static int
qdpc_w_close(lua_State *L)
{
    mWriter *b = q_checkWriter(L, 1, NULL);

    check_writer(L, b);
    
    QDP_close_write(b->ptr);
    b->ptr = 0;

    return 0;
}

static int
qdpc_w_iv(lua_State *L)
{
    mWriter *writer = q_checkWriter(L, 1, NULL);
    mVecInt *v = qlua_checkVecInt(L, 2);
    const char *info = luaL_checkstring(L, 3);
    QDP_String *xml;
    
    check_writer(L, writer);

    CALL_QDP(L);

    xml = QDP_string_create();
    QDP_string_set(xml, (char *)info); /* [ sic ] */

    if (QDP_vwrite_i(writer->ptr, xml, v->val, v->size) == 0) {
        lua_pushboolean(L, 1);
        QDP_string_destroy(xml);
        return 1;
    } else {
        QDP_string_destroy(xml);
        return luaL_error(L, "qdpc write error");
    }
}

static int
qdpc_w_rv(lua_State *L)
{
    mWriter *writer = q_checkWriter(L, 1, NULL);
    mVecReal *v = qlua_checkVecReal(L, 2);
    const char *info = luaL_checkstring(L, 3);
    QDP_String *xml;
    
    check_writer(L, writer);

    CALL_QDP(L);

    xml = QDP_string_create();
    QDP_string_set(xml, (char *)info); /* [ sic ] */

    if (QDP_D_vwrite_r(writer->ptr, xml, v->val, v->size) == 0) {
        lua_pushboolean(L, 1);
        QDP_string_destroy(xml);
        return 1;
    } else {
        QDP_string_destroy(xml);
        return luaL_error(L, "qdpc write error");
    }
}

static int
qdpc_w_cv(lua_State *L)
{
    mWriter *writer = q_checkWriter(L, 1, NULL);
    mVecComplex *v = qlua_checkVecComplex(L, 2);
    const char *info = luaL_checkstring(L, 3);
    QDP_String *xml;
    
    check_writer(L, writer);

    CALL_QDP(L);

    xml = QDP_string_create();
    QDP_string_set(xml, (char *)info); /* [ sic ] */

    if (QDP_D_vwrite_c(writer->ptr, xml, v->val, v->size) == 0) {
        lua_pushboolean(L, 1);
        QDP_string_destroy(xml);
        return 1;
    } else {
        QDP_string_destroy(xml);
        return luaL_error(L, "qdpc write error");
    }
}

#ifdef HAS_GSL
static int
qdpc_w_rm(lua_State *L)
{
    mWriter *writer = q_checkWriter(L, 1, NULL);
    mMatReal *v = qlua_checkMatReal(L, 2);
    const char *info = luaL_checkstring(L, 3);
    int l_size = v->l_size;
    int r_size = v->r_size;
    QLA_Int f_size[2];
    QDP_String *xml;
    int i, j;

    check_writer(L, writer);

    f_size[0] = l_size;
    f_size[1] = r_size;
    QLA_D_Real d[l_size * r_size];
    /* must only agree with the reader above */
    for (j = 0; j < r_size; j++) {
        for (i = 0; i < l_size; i++) {
            d[i + l_size * j] = gsl_matrix_get(v->m, i, j);
        }
    }

    CALL_QDP(L);

    xml = QDP_string_create();
    QDP_string_set(xml, (char *)qio_matrix_magic); /* [ sic ] */
    if (QDP_vwrite_i(writer->ptr, xml, f_size, 2) != 0)
        goto error;
    QDP_string_set(xml, (char *)info); /* [ sic ] */
    if (QDP_D_vwrite_r(writer->ptr, xml, d, l_size * r_size) != 0)
        goto error;
    QDP_string_destroy(xml);
    lua_pushboolean(L, 1);
    return 1;
error:
    QDP_string_destroy(xml);
    return luaL_error(L, "qdpc write error");
}

static int
qdpc_w_cm(lua_State *L)
{
    mWriter *writer = q_checkWriter(L, 1, NULL);
    mMatComplex *v = qlua_checkMatComplex(L, 2);
    const char *info = luaL_checkstring(L, 3);
    int l_size = v->l_size;
    int r_size = v->r_size;
    QLA_Int f_size[2];
    QDP_String *xml;
    int i, j;

    check_writer(L, writer);

    f_size[0] = l_size;
    f_size[1] = r_size;
    QLA_D_Complex d[l_size * r_size];
    /* must only agree with the reader above */
    for (j = 0; j < r_size; j++) {
        for (i = 0; i < l_size; i++) {
            gsl_complex z = gsl_matrix_complex_get(v->m, i, j);
            QLA_real(d[i + l_size * j]) = GSL_REAL(z);
            QLA_imag(d[i + l_size * j]) = GSL_IMAG(z);
        }
    }

    CALL_QDP(L);

    xml = QDP_string_create();
    QDP_string_set(xml, (char *)qio_matrix_magic); /* [ sic ] */
    if (QDP_vwrite_i(writer->ptr, xml, f_size, 2) != 0)
        goto error;
    QDP_string_set(xml, (char *)info); /* [ sic ] */
    if (QDP_D_vwrite_c(writer->ptr, xml, d, l_size * r_size) != 0)
        goto error;
    QDP_string_destroy(xml);
    lua_pushboolean(L, 1);
    return 1;
error:
    QDP_string_destroy(xml);
    return luaL_error(L, "qdpc write error");
}
#endif

static const struct luaL_Reg mtWriter[] = {
    { "__tostring",       qdpc_w_fmt             },
    { "__gc",             qdpc_w_gc              },
    { "__newindex",       qlua_nowrite           },
    { "close",            qdpc_w_close           },
    { "Int",              qdpc_w_I               },
    { "RandomState",      qdpc_w_S               },
    { "Real",             qdpc_w_R               },
    { "Complex",          qdpc_w_C               },
#if USE_Nc2 || USE_Nc3 || USE_NcN
    { "ColorVector",      qdpc_w_ColorVector     },
    { "ColorMatrix",      qdpc_w_ColorMatrix     },
    { "DiracFermion",     qdpc_w_DiracFermion    },
    { "DiracPropagator",  qdpc_w_DiracPropagator },
#endif /* use any colors */
    { "int_vector",       qdpc_w_iv              },
    { "real_vector",      qdpc_w_rv              },
    { "complex_vector",   qdpc_w_cv              },
#ifdef HAS_GSL
    { "real_matrix",      qdpc_w_rm              },
    { "complex_matrix",   qdpc_w_cm              },
#endif
    /* "lattice" */
    /* "a-type" */
    { NULL,               NULL                   }
};

static mWriter *
q_checkWriter(lua_State *L, int idx, mLattice *S)
{
    void *v = qlua_checkLatticeType(L, idx, qWriter, WriterName);
    
    if (S) {
        mLattice *S1 = qlua_ObjLattice(L, idx);
        if (S1->id != S->id)
            luaL_error(L, "%s on a wrong lattice", WriterName);
        lua_pop(L, 1);
    }

    return (mWriter *)v;
}

static mWriter *
q_newWriter(lua_State *L, QDP_Writer *writer)
{
    mWriter *h = lua_newuserdata(L, sizeof (mWriter));

    h->ptr = writer;
    qlua_createLatticeTable(L, 1, mtWriter, qWriter, WriterName);
    lua_setmetatable(L, -2);

    return h;
}

/* QDP/C writer constructors
 *
 *  qcd.qdpc.writer(L, file_name, file_info)
 *  qcd.qdpc.writer(L, file_name, file_info, "single")
 *  qcd.qdpc.writer(L, file_name, file_info, "multi")
 */
static int
q_qdpc_writer (lua_State *L)
{
    mLattice *S = qlua_checkLattice(L, 1);
    const char *name = luaL_checkstring(L, 2);
    const char *info = luaL_checkstring(L, 3);
    int volfmt = qlua_qio_volume_format(L, 4, lua_gettop(L));
    QDP_String *xml = 0;
    QDP_Writer *writer = 0;
    int rcount;

    /* first collect garbage */
    CALL_QDP(L);
    
    /* open the QDP writer */
    xml = QDP_string_create();
    QDP_string_set(xml, (char *)info); /* [ sic ] */
    writer = QDP_open_write_L(S->lat, xml, (char *)name, volfmt); /* [ sic ] */

    /* convert QDP results into LUA values */
    if (writer == 0) {
        /* failure - zero results */
        rcount = 0;
    } else {
        /* success -- return writer */
        q_newWriter(L, writer);
        rcount = 1;
    }
    QDP_string_destroy(xml);

    return rcount;
}

/* names and routines for qcd.qdpc table */
static const struct luaL_Reg fQDPio[] = {
    { "Reader",   q_qdpc_reader},
    { "Writer",   q_qdpc_writer},
    { NULL,       NULL }
};

int
init_qdpc_io(lua_State *L)
{
    lua_getglobal(L, qcdlib);
    lua_newtable(L);
    luaL_register(L, NULL, fQDPio);
    lua_setfield(L, -2, qdp_io);
    lua_pop(L, 1);
    
    return 0;
}

int
fini_qdpc_io(lua_State *L)
{
    return 0;
}
