#include <qlua.h>                                                    /* DEPS */
#include <qdpc_io.h>                                                 /* DEPS */
#include <latcolmat.h>                                               /* DEPS */
#include <latcolvec.h>                                               /* DEPS */
#include <latcomplex.h>                                              /* DEPS */
#include <latdirferm.h>                                              /* DEPS */
#include <latdirprop.h>                                              /* DEPS */
#include <latint.h>                                                  /* DEPS */
#include <latrandom.h>                                               /* DEPS */
#include <latreal.h>                                                 /* DEPS */

typedef struct {
    QDP_Reader *ptr;
} mReader;

typedef struct {
    QDP_Writer *ptr;
} mWriter;

const char qdp_io[] = "qdpc";
static const char mtnReader[] = "qcd.qdpc.reader";
static const char mtnWriter[] = "qcd.qdpc.writer";

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
static int
qdpc_r_info(lua_State *L)
{
    mReader *reader = q_checkReader(L, 1);
    QDP_String *xml;
    int status;
    int rcount;

    lua_gc(L, LUA_GCCOLLECT, 0);
    xml = QDP_string_create();
    status = QDP_read_record_info(reader->ptr, xml);
    if (status == 0) {
        lua_pushboolean(L, 1 );
        lua_pushstring(L, QDP_string_ptr(xml));
        rcount = 2;
    } else {
        rcount = 0;
    }
    QDP_string_destroy(xml);

    return rcount;
}

static int
qdpc_r_skip(lua_State *L)
{
    mReader *reader = q_checkReader(L, 1);
    int status;
    int rcount;

    lua_gc(L, LUA_GCCOLLECT, 0);
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
    
    if (b->ptr)
        QDP_close_read(b->ptr);
    b->ptr = 0;
    
    return 0;
}

static int
qdpc_w_close(lua_State *L)
{
    mWriter *b = q_checkWriter(L, 1);
    
    if (b->ptr)
        QDP_close_write(b->ptr);
    b->ptr = 0;

    return 0;
}


/* lattice readers and writers are almost identical for all lattice types */
#define T_QTYPE       QDP_Int
#define QLUA_NAME(x)  x ## LatInt
#define X_ID(x)       x ## I
#include <qdpc_io_template.c>                                        /* DEPS */

#define T_QTYPE       QDP_RandomState
#define QLUA_NAME(x)  x ## LatRandom
#define X_ID(x)       x ## S
#include <qdpc_io_template.c>                                        /* DEPS */

#define T_QTYPE       QDP_Real
#define QLUA_NAME(x)  x ## LatReal
#define X_ID(x)       x ## R
#include <qdpc_io_template.c>                                        /* DEPS */

#define T_QTYPE       QDP_Complex
#define QLUA_NAME(x)  x ## LatComplex
#define X_ID(x)       x ## C
#include <qdpc_io_template.c>                                        /* DEPS */

#define T_QTYPE       QDP_ColorVector
#define QLUA_NAME(x)  x ## LatColVec
#define X_ID(x)       x ## V
#include <qdpc_io_template.c>                                        /* DEPS */

#define T_QTYPE       QDP_ColorMatrix
#define QLUA_NAME(x)  x ## LatColMat
#define X_ID(x)       x ## M
#include <qdpc_io_template.c>                                        /* DEPS */

#define T_QTYPE       QDP_DiracFermion
#define QLUA_NAME(x)  x ## LatDirFerm
#define X_ID(x)       x ## D
#include <qdpc_io_template.c>                                        /* DEPS */

#define T_QTYPE       QDP_DiracPropagator
#define QLUA_NAME(x)  x ## LatDirProp
#define X_ID(x)       x ## P
#include <qdpc_io_template.c>                                        /* DEPS */


static int
q_qdpc_reader (lua_State *L)
{
    const char *name = luaL_checkstring(L, 1);
    QDP_String *xml = 0;
    QDP_Reader *reader = 0;
    int rcount;

    /* first collect garbage */
    lua_gc(L, LUA_GCCOLLECT, 0);

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

    /* first collect garbage */
    lua_gc(L, LUA_GCCOLLECT, 0);
    
    /* open the QDP writer */
    xml = QDP_string_create();
    QDP_string_set(xml, (char *)info); /* [ sic ] */
    writer = QDP_open_write(xml, (char *)name, QDP_SINGLEFILE); /* [ sic ] */

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
    { "info",             qdpc_r_info },
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
    { NULL,               NULL}
};

/* names and routines for qcd.qdpc table */
static const struct {
    char *name;
    int (*func)(lua_State *L);
} fQDPio[] = {
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
