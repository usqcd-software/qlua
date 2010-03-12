#include "qlua.h"                                                    /* DEPS */
#include "qdpc_io.h"                                                 /* DEPS */
#include "latcolmat.h"                                               /* DEPS */
#include "latcolvec.h"                                               /* DEPS */
#include "latcomplex.h"                                              /* DEPS */
#include "latdirferm.h"                                              /* DEPS */
#include "latdirprop.h"                                              /* DEPS */
#include "latint.h"                                                  /* DEPS */
#include "latrandom.h"                                               /* DEPS */
#include "latreal.h"                                                 /* DEPS */
#include "qvector.h"                                                 /* DEPS */
#include "qmatrix.h"                                                 /* DEPS */

#include <string.h>

typedef struct {
    QDP_Reader *ptr;
} mReader;

typedef struct {
    QDP_Writer *ptr;
} mWriter;

/* this magic used to identify global matrices in QIO */
static const char qio_matrix_magic[] = "480C94A9-F941-4FD5-9882-3E70B004AACF";

const char qdp_io[] = "qdpc";
static const char mtnReader[] = "qcd.qdpc.reader";
static const char mtnWriter[] = "qcd.qdpc.writer";

/* helpers */
static void
check_reader(lua_State *L, mReader *b)
{
    if (b->ptr == 0)
        luaL_error(L, "closed QDP/C reader");
}

static void
check_writer(lua_State *L, mWriter *b)
{
    if (b->ptr == 0)
        luaL_error(L, "closed QDP/C writer");
}

/* allocation */
static mReader *
q_newReader(lua_State *L, QDP_Reader *reader)
{
    mReader *h = lua_newuserdata(L, sizeof (mReader));

    h->ptr = reader;
    luaL_getmetatable(L, mtnReader);
    lua_setmetatable(L, -2);

    return h;
}

static mWriter *
q_newWriter(lua_State *L, QDP_Writer *writer)
{
    mWriter *h = lua_newuserdata(L, sizeof (mWriter));

    h->ptr = writer;
    luaL_getmetatable(L, mtnWriter);
    lua_setmetatable(L, -2);

    return h;
}

/* type checking */
static mReader *
q_checkReader(lua_State *L, int idx)
{
    void *v = luaL_checkudata(L, idx, mtnReader);

    luaL_argcheck(L, v != 0, idx, "qcd.qdpc.Reader expected");

    return v;
}

static mWriter *
q_checkWriter(lua_State *L, int idx)
{
    void *v = luaL_checkudata(L, idx, mtnWriter);

    luaL_argcheck(L, v != 0, idx, "qcd.qdpc.Writer expected");

    return v;
}

/* conversion to string */
static int
qdpc_r_fmt(lua_State *L)
{
    mReader *b = q_checkReader(L, 1);
    char fmt[72];
    
    if ( b->ptr )
        sprintf(fmt, "qdpc.Reader(%p)", b->ptr);
    else
        sprintf(fmt, "qdpc.Reader(closed)");

    lua_pushstring(L, fmt);

    return 1;
}

static int
qdpc_w_fmt(lua_State *L)
{
    mWriter *b = q_checkWriter(L, 1);
    char fmt[72];
    
    if (b->ptr)
        sprintf(fmt, "qdpc.Writer(%p)", b->ptr);
    else
        sprintf(fmt, "qdpc.Writer(closed)");

    lua_pushstring(L, fmt);

    return 1;
}

/* garbage collection */
static int
qdpc_r_gc(lua_State *L)
{
    mReader *b = q_checkReader(L, 1);

    if (b->ptr)
        QDP_close_read(b->ptr);
    b->ptr = 0;
    
    return 0;
}

static int
qdpc_w_gc(lua_State *L)
{
    mWriter *b = q_checkWriter(L, 1);

    if (b->ptr)
        QDP_close_write(b->ptr);
    b->ptr = 0;
    
    return 0;
}

/* reader manipulations */
#if 0
/* QIO error handling broken in MPI implementation */
static int
qdpc_r_info(lua_State *L)
{
    mReader *reader = q_checkReader(L, 1);
    QDP_String *xml;
    int status;
    int rcount;
#if 1
    static int count = 0;
    count++;
    printf("r_info/%d: enter\n", count);
    fflush(stdout);
#endif

    check_reader(L, reader);

    CALL_QDP(L);
    xml = QDP_string_create();
    status = QDP_read_record_info(reader->ptr, xml);
    if (status == 0) {
        lua_pushboolean(L, 1 );
        lua_pushstring(L, QDP_string_ptr(xml));
#if 1
        printf("r_info/%d: 2 results\n", count);
        fflush(stdout);
#endif
        rcount = 2;
    } else {
#if 1
        printf("r_info/%d: no results\n", count);
        fflush(stdout);
#endif
        rcount = 0;
    }
    QDP_string_destroy(xml);

    return rcount;
}
#endif

static int
qdpc_r_skip(lua_State *L)
{
    mReader *reader = q_checkReader(L, 1);
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

/* closing */
static int
qdpc_r_close(lua_State *L)
{
    mReader *b = q_checkReader(L, 1);
    
    check_reader(L, b);

    QDP_close_read(b->ptr);
    b->ptr = 0;
    
    return 0;
}

static int
qdpc_w_close(lua_State *L)
{
    mWriter *b = q_checkWriter(L, 1);

    check_writer(L, b);
    
    QDP_close_write(b->ptr);
    b->ptr = 0;

    return 0;
}


/* lattice readers and writers are almost identical for all lattice types */
#define T_QTYPE       QDP_Int
#define QLUA_NAME(x)  x ## LatInt
#define X_ID(x)       x ## I
#include "qdpc_io_template.c"                                        /* DEPS */

#define T_QTYPE       QDP_RandomState
#define QLUA_NAME(x)  x ## LatRandom
#define X_ID(x)       x ## S
#include "qdpc_io_template.c"                                        /* DEPS */

#define T_QTYPE       QDP_Real
#define QLUA_NAME(x)  x ## LatReal
#define X_ID(x)       x ## R
#include "qdpc_io_template.c"                                        /* DEPS */

#define T_QTYPE       QDP_Complex
#define QLUA_NAME(x)  x ## LatComplex
#define X_ID(x)       x ## C
#include "qdpc_io_template.c"                                        /* DEPS */

#define T_QTYPE       QDP_ColorVector
#define QLUA_NAME(x)  x ## LatColVec
#define X_ID(x)       x ## V
#include "qdpc_io_template.c"                                        /* DEPS */

#define T_QTYPE       QDP_ColorMatrix
#define QLUA_NAME(x)  x ## LatColMat
#define X_ID(x)       x ## M
#include "qdpc_io_template.c"                                        /* DEPS */

#define T_QTYPE       QDP_DiracFermion
#define QLUA_NAME(x)  x ## LatDirFerm
#define X_ID(x)       x ## D
#include "qdpc_io_template.c"                                        /* DEPS */

#define T_QTYPE       QDP_DiracPropagator
#define QLUA_NAME(x)  x ## LatDirProp
#define X_ID(x)       x ## P
#include "qdpc_io_template.c"                                        /* DEPS */

static int
qdpc_r_iv(lua_State *L)
{
    mReader *reader = q_checkReader(L, 1);
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
qdpc_w_iv(lua_State *L)
{
    mWriter *writer = q_checkWriter(L, 1);
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
qdpc_r_rv(lua_State *L)
{
    mReader *reader = q_checkReader(L, 1);
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
qdpc_w_rv(lua_State *L)
{
    mWriter *writer = q_checkWriter(L, 1);
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
qdpc_r_cv(lua_State *L)
{
    mReader *reader = q_checkReader(L, 1);
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

static int
qdpc_w_cv(lua_State *L)
{
    mWriter *writer = q_checkWriter(L, 1);
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

static int
qdpc_r_rm(lua_State *L)
{
    mReader *reader = q_checkReader(L, 1);
    int l_size = luaL_checkint(L, 2);
    int r_size = luaL_checkint(L, 3);
    mMatReal *v = qlua_newMatReal(L, l_size, r_size);
    QLA_Int f_size[2];
    QLA_Real *d;
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
    d = qlua_malloc(L, l_size * r_size *  sizeof (QLA_Real));
    if (QDP_D_vread_r(reader->ptr, info, d, l_size * r_size) != 0) {
        qlua_free(L, d);
        goto error;
    }
    /* must only agree with the writer below */
    for (j = 0; j < r_size; j++) {
        for (i = 0; i < l_size; i++) {
            gsl_matrix_set(v->m, i, j, d[i + l_size * j]);
        }
    }
    qlua_free(L, d);
    lua_pushstring(L, QDP_string_ptr(info));
    QDP_string_destroy(info);
    return 2;
error:
    QDP_string_destroy(info);
    return luaL_error(L, "qdpc read error");
}

static int
qdpc_w_rm(lua_State *L)
{
    mWriter *writer = q_checkWriter(L, 1);
    mMatReal *v = qlua_checkMatReal(L, 2);
    const char *info = luaL_checkstring(L, 3);
    int l_size = v->l_size;
    int r_size = v->r_size;
    QLA_Int f_size[2];
    QLA_Real *d;
    QDP_String *xml;
    int i, j;

    check_writer(L, writer);

    f_size[0] = l_size;
    f_size[1] = r_size;
    d = qlua_malloc(L, l_size * r_size * sizeof (QLA_Real));
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
    qlua_free(L, d);
    lua_pushboolean(L, 1);
    return 1;
error:
    QDP_string_destroy(xml);
    qlua_free(L, d);
    return luaL_error(L, "qdpc write error");
}

static int
qdpc_r_cm(lua_State *L)
{
    mReader *reader = q_checkReader(L, 1);
    int l_size = luaL_checkint(L, 2);
    int r_size = luaL_checkint(L, 3);
    mMatComplex *v = qlua_newMatComplex(L, l_size, r_size);
    QLA_Int f_size[2];
    QLA_Complex *d;
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
    d = qlua_malloc(L, l_size * r_size *  sizeof (QLA_Complex));
    if (QDP_D_vread_c(reader->ptr, info, d, l_size * r_size) != 0) {
        qlua_free(L, d);
        goto error;
    }
    /* must only agree with the writer below */
    for (j = 0; j < r_size; j++) {
        for (i = 0; i < l_size; i++) {
            QLA_Complex z = d[i + l_size * j];
            gsl_complex zz;

            GSL_SET_COMPLEX(&zz, QLA_real(z), QLA_imag(z));
            gsl_matrix_complex_set(v->m, i, j, zz);
        }
    }
    qlua_free(L, d);
    lua_pushstring(L, QDP_string_ptr(info));
    QDP_string_destroy(info);
    return 2;
error:
    QDP_string_destroy(info);
    return luaL_error(L, "qdpc read error");
}

static int
qdpc_w_cm(lua_State *L)
{
    mWriter *writer = q_checkWriter(L, 1);
    mMatComplex *v = qlua_checkMatComplex(L, 2);
    const char *info = luaL_checkstring(L, 3);
    int l_size = v->l_size;
    int r_size = v->r_size;
    QLA_Int f_size[2];
    QLA_Complex *d;
    QDP_String *xml;
    int i, j;

    check_writer(L, writer);

    f_size[0] = l_size;
    f_size[1] = r_size;
    d = qlua_malloc(L, l_size * r_size * sizeof (QLA_Complex));
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
    qlua_free(L, d);
    lua_pushboolean(L, 1);
    return 1;
error:
    QDP_string_destroy(xml);
    qlua_free(L, d);
    return luaL_error(L, "qdpc write error");
}

/* end of readers and writers for QLUA data types */

static int
q_qdpc_reader (lua_State *L)
{
    const char *name = luaL_checkstring(L, 1);
    QDP_String *xml = 0;
    QDP_Reader *reader = 0;
    int rcount;

    /* first collect garbage */
    CALL_QDP(L);

    /* go through the motions of opening a QDP reader */
    xml = QDP_string_create();
    reader = QDP_open_read(xml, (char *)name); /* [sic] */

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

static int
q_qdpc_writer (lua_State *L)
{
    const char *name = luaL_checkstring(L, 1);
    const char *info = luaL_checkstring(L, 2);
    QDP_String *xml = 0;
    QDP_Writer *writer = 0;
    int rcount;
    int volfmt;

    switch (lua_type(L, 3)) {
    case LUA_TNONE:
    case LUA_TNIL:
        volfmt = QDP_SINGLEFILE;
        break;
    default: { /* must be LUA_TSTRING, which we check */
        const char *fmt = luaL_checkstring(L, 3);
        if (strcmp(fmt, "single") == 0)
            volfmt = QDP_SINGLEFILE;
        else if (strcmp(fmt, "multi") == 0)
            volfmt = QDP_MULTIFILE;
        else
            return luaL_error(L, "unsupported file format");
    }
    }

    /* first collect garbage */
    CALL_QDP(L);
    
    /* open the QDP writer */
    xml = QDP_string_create();
    QDP_string_set(xml, (char *)info); /* [ sic ] */
    writer = QDP_open_write(xml, (char *)name, volfmt); /* [ sic ] */

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

/* metatables */
static const struct luaL_Reg mtReader[] = {
    { "__tostring",       qdpc_r_fmt},
    { "__gc",             qdpc_r_gc},
#if 0
/* QIO error handling broken in MPI  */
    { "info",             qdpc_r_info },
#endif
    { "skip",             qdpc_r_skip },
    { "close",            qdpc_r_close},
    { "Int",              qdpc_r_I},
    { "RandomState",      qdpc_r_S},
    { "Real",             qdpc_r_R},
    { "Complex",          qdpc_r_C},
    { "ColorVector",      qdpc_r_V},
    { "ColorMatrix",      qdpc_r_M},
    { "DiracFermion",     qdpc_r_D},
    { "DiracPropagator",  qdpc_r_P},
    { "int_vector",       qdpc_r_iv },
    { "real_vector",      qdpc_r_rv },
    { "complex_vector",   qdpc_r_cv },
    { "real_matrix",      qdpc_r_rm },
    { "complex_matrix",   qdpc_r_cm },
    { NULL,               NULL}
};

static const struct luaL_Reg mtWriter[] = {
    { "__tostring",       qdpc_w_fmt},
    { "__gc",             qdpc_w_gc},
    { "close",            qdpc_w_close},
    { "Int",              qdpc_w_I},
    { "RandomState",      qdpc_w_S},
    { "Real",             qdpc_w_R},
    { "Complex",          qdpc_w_C},
    { "ColorVector",      qdpc_w_V},
    { "ColorMatrix",      qdpc_w_M},
    { "DiracFermion",     qdpc_w_D},
    { "DiracPropagator",  qdpc_w_P},
    { "int_vector",       qdpc_w_iv },
    { "real_vector",      qdpc_w_rv },
    { "complex_vector",   qdpc_w_cv },
    { "real_matrix",      qdpc_w_rm },
    { "complex_matrix",   qdpc_w_cm },
    { NULL,               NULL}
};

/* names and routines for qcd.qdpc table */
static const struct luaL_Reg fQDPio[] = {
    { "Reader",   q_qdpc_reader},
    { "Writer",   q_qdpc_writer},
    { NULL,       NULL }
};

int
init_qdpc_io(lua_State *L)
{
    int i;

    lua_getglobal(L, qcdlib);
    lua_newtable(L);
    for (i = 0; fQDPio[i].name; i++) {
        lua_pushcfunction(L, fQDPio[i].func);
        lua_setfield(L, -2, fQDPio[i].name);
    }
    lua_setfield(L, -2, qdp_io);
    qlua_metatable(L, mtnReader, mtReader);
    qlua_metatable(L, mtnWriter, mtWriter);
    
    return 0;
}

int
fini_qdpc_io(lua_State *L)
{
    return 0;
}
