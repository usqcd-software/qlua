#include <qlua.h>                                                    /* DEPS */
#include <qvector.h>                                                 /* DEPS */
#include <lhpc-aff.h>
#include <aff_io.h>                                                  /* DEPS */
#include <string.h>
#include <complex.h>

const char aff_io[] = "aff";
static const char mtnReader[] = "qcd.aff.reader";
static const char mtnWriter[] = "qcd.aff.writer";

/* aff replacement allocator */
static lua_State *qL = NULL; /* to pass the LUA state to aff_realloc */
void *
aff_realloc(void *ptr, size_t size)
{
    if (ptr == 0) {
        return qlua_malloc(qL, (int)size);
    } else if (size == 0) {
        qlua_free(qL, ptr);
        return NULL;
    }
    luaL_error(qL, "bad parameters to aff_realloc");
    /* can not happen */

    return NULL;
}

/* allocation */
static mAffReader *
qlua_newAffReader(lua_State *L, struct AffReader_s *reader)
{
    mAffReader *h = lua_newuserdata(L, sizeof (mAffReader));

    h->ptr = reader;
    h->dir = aff_reader_root(reader);
    luaL_getmetatable(L, mtnReader);
    lua_setmetatable(L, -2);

    return h;
}

static mAffWriter *
qlua_newAffWriter(lua_State *L, struct AffWriter_s *writer)
{
    mAffWriter *h = lua_newuserdata(L, sizeof (mAffWriter));

    h->ptr = writer;
    h->dir = aff_writer_root(writer);
    luaL_getmetatable(L, mtnWriter);
    lua_setmetatable(L, -2);

    return h;
}

/* type checking */
mAffReader *
qlua_checkAffReader(lua_State *L, int idx)
{
    void *v = luaL_checkudata(L, idx, mtnReader);

    luaL_argcheck(L, v != 0, idx, "qcd.aff.Reader expected");

    return v;
}

mAffWriter *
qlua_checkAffWriter(lua_State *L, int idx)
{
    void *v = luaL_checkudata(L, idx, mtnWriter);

    luaL_argcheck(L, v != 0, idx, "qcd.aff.Writer expected");

    return v;
}

/* conversion to string */
static int
qaff_r_fmt(lua_State *L)
{
    mAffReader *b = qlua_checkAffReader(L, 1);
    char fmt[72];
    
    if ( b->ptr )
        sprintf(fmt, "aff.Reader(%p)", b->ptr);
    else
        sprintf(fmt, "aff.Reader(closed)");

    lua_pushstring(L, fmt);

    return 1;
}

static int
qaff_w_fmt(lua_State *L)
{
    mAffWriter *b = qlua_checkAffWriter(L, 1);
    char fmt[72];
    
    if (b->ptr)
        sprintf(fmt, "aff.Writer(%p)", b->ptr);
    else
        sprintf(fmt, "aff.Writer(closed)");

    lua_pushstring(L, fmt);

    return 1;
}

/* garbage collection */
static int
qaff_r_gc(lua_State *L)
{
    mAffReader *b = qlua_checkAffReader(L, 1);

    qL = L;
    if (b->ptr)
        aff_reader_close(b->ptr);
    b->ptr = 0;
    
    return 0;
}

static int
qaff_w_gc(lua_State *L)
{
    mAffWriter *b = qlua_checkAffWriter(L, 1);

    qL = L;
    if (b->ptr)
        aff_writer_close(b->ptr);
    b->ptr = 0;
    
    return 0;
}

/* closing */
static int
qaff_r_close(lua_State *L)
{
    mAffReader *b = qlua_checkAffReader(L, 1);
    
    qL = L;
    if (b->ptr)
        aff_reader_close(b->ptr);
    b->ptr = 0;
    
    return 0;
}

static int
qaff_w_close(lua_State *L)
{
    mAffWriter *b = qlua_checkAffWriter(L, 1);
    
    qL = L;
    if (b->ptr) {
        const char *s = aff_writer_close(b->ptr);

        b->ptr = 0;
        if (s == 0)
            lua_pushnil(L);
        else
            lua_pushstring(L, s);

        return 1;
    }

    return 0;
}

static int
qaff_r_status(lua_State *L)
{
    mAffReader *b = qlua_checkAffReader(L, 1);
    const char *s;

    if (b->ptr == 0)
        return luaL_error(L, "closed reader");

    qL = L;
    s = aff_reader_errstr(b->ptr);
        
    if (s == 0)
        lua_pushnil(L);
    else
        lua_pushstring(L, s);

    return 1;
}

static int
qaff_w_status(lua_State *L)
{
    mAffWriter *b = qlua_checkAffWriter(L, 1);
    const char *s;
    
    if (b->ptr == 0)
        return luaL_error(L, "closed writer");

    qL = L;
    s = aff_writer_errstr(b->ptr);

    if (s == 0) 
        lua_pushnil(L);
    else
        lua_pushstring(L, s);
        
    return 1;
}

/* XXX */ static int qaff_r_list(lua_State *L) { return 0; }
static int
qaff_r_read(lua_State *L)
{
    mAffReader *b = qlua_checkAffReader(L, 1);
    const char *p = luaL_checkstring(L, 2);
    struct AffNode_s *r;
    struct AffNode_s *n;
    uint32_t size;

    if (b->ptr == 0)
        goto end;

    qL = L;
    r = p[0] == '/'? NULL: b->dir;
    n = aff_reader_chpath(b->ptr, r, p);
    if (n == 0)
        goto end;
    size = aff_node_size(n);

    switch (aff_node_type(n)) {
    case affNodeVoid:
        lua_pushboolean(L, 1);
        return 1;
    case affNodeChar: {
        char *d = qlua_malloc(L, size + 1);

        if (aff_node_get_char(b->ptr, n, d, size) == 0) {
            d[size] = 0;
            lua_pushstring(L, d);
            qlua_free(L, d);

            return 1;
        } else {
            qlua_free(L, d);
        }
        break;
    }
    case affNodeInt: {
        uint32_t *d = qlua_malloc(L, size * sizeof (uint32_t));

        if (aff_node_get_int(b->ptr, n, d, size) == 0) {
            tVecInt *v = qlua_newVecInt(L, size);
            int i;

            v->size = size;
            for (i = 0; i < size; i++)
                v->val[i] = d[i];

            qlua_free(L, d);

            return 1;
        } else {
            qlua_free(L, d);
        }
        break;
    }
    case affNodeDouble: {
        double *d = qlua_malloc(L, size * sizeof (double));

        if (aff_node_get_double(b->ptr, n, d, size) == 0) {
            tVecReal *v = qlua_newVecReal(L, size);
            int i;

            v->size = size;
            for (i = 0; i < size; i++)
                v->val[i] = d[i];

            qlua_free(L, d);

            return 1;
        } else {
            qlua_free(L, d);
        }
        break;
    }
    case affNodeComplex: {
        double _Complex *d = qlua_malloc(L, size * sizeof (double _Complex));

        if (aff_node_get_complex(b->ptr, n, d, size) == 0) {
            tVecComplex *v = qlua_newVecComplex(L, size);
            int i;

            v->size = size;
            for (i = 0; i < size; i++) {
                QLA_real(v->val[i]) = creal(d[i]);
                QLA_imag(v->val[i]) = cimag(d[i]);
            }

            qlua_free(L, d);

            return 1;
        } else {
            qlua_free(L, d);
        }
        break;
    }
    default:
        goto end;
    }
    return luaL_error(L, aff_reader_errstr(b->ptr));

end:
    return luaL_error(L, "bad arguments");
}

static int
qaff_w_write(lua_State *L)
{
    mAffWriter *b = qlua_checkAffWriter(L, 1);
    const char *p = luaL_checkstring(L, 2);
    struct AffNode_s *r;
    struct AffNode_s *n;

    if (b->ptr == 0)
        goto end;

    qL = L;
    r = p[0] == '/'? NULL: b->dir;
    n = aff_writer_mkpath(b->ptr, r, p);
    if (n == 0)
        goto end;

    switch (qlua_gettype(L, 3)) {
    case qString: {
        const char *str = luaL_checkstring(L, 3);

        if (aff_node_put_char(b->ptr, n, str, strlen(str)))
            return luaL_error(L, aff_writer_errstr(b->ptr));

        return 0;
    }
    case qVecInt: {
        tVecInt *v = qlua_checkVecInt(L, 3);
        int size = v->size;
        uint32_t *d = qlua_malloc(L, size * sizeof (uint32_t));
        int i;

        for (i = 0; i < size; i++)
            d[i] = v->val[i];

        if (aff_node_put_int(b->ptr, n, d, size))
            return luaL_error(L, aff_writer_errstr(b->ptr));

        qlua_free(L, d);

        return 0;
    }
    case qVecReal: {
        tVecReal *v = qlua_checkVecReal(L, 3);
        int size = v->size;
        double *d = qlua_malloc(L, size * sizeof (double));
        int i;

        for (i = 0; i < size; i++)
            d[i] = v->val[i];

        if (aff_node_put_double(b->ptr, n, d, size))
            return luaL_error(L, aff_writer_errstr(b->ptr));

        qlua_free(L, d);

        return 0;
    }
    case qVecComplex: {
        tVecComplex *v = qlua_checkVecComplex(L, 3);
        int size = v->size;
        double _Complex *d = qlua_malloc(L, size * sizeof (double _Complex));
        int i;

        for (i = 0; i < size; i++) {
            d[i] = QLA_real(v->val[i]) + I * QLA_imag(v->val[i]);
        }

        if (aff_node_put_complex(b->ptr, n, d, size))
            return luaL_error(L, aff_writer_errstr(b->ptr));

        qlua_free(L, d);

        return 0;
    }
    }
end:
    return luaL_error(L, "bad arguments");
}

static int
qaff_r_chpath(lua_State *L)
{
    mAffReader *b = qlua_checkAffReader(L, 1);
    const char *p = luaL_checkstring(L, 2);
    struct AffNode_s *r;

    if (b->ptr == 0)
        return luaL_error(L, "closed reader");

    qL = L;
    r = p[0] == '/'? NULL: b->dir;
    r = aff_reader_chpath(b->ptr, r, p);
    if (r == 0)
        return luaL_error(L, aff_reader_errstr(b->ptr));

    b->dir = r;
    return 0;
}

static int
qaff_w_chpath(lua_State *L)
{
    mAffWriter *b = qlua_checkAffWriter(L, 1);
    const char *p = luaL_checkstring(L, 2);
    struct AffNode_s *r;

    if (b->ptr == 0)
        return luaL_error(L, "closed writer");

    qL = L;
    r = p[0] == '/'? NULL: b->dir;
    r = aff_writer_mkpath(b->ptr, r, p);
    if (r == NULL)
        return luaL_error(L, aff_writer_errstr(b->ptr));

    b->dir = r;
    return 0;
}

static int
qaff_w_mkpath(lua_State *L)
{
    mAffWriter *b = qlua_checkAffWriter(L, 1);
    const char *p = luaL_checkstring(L, 2);
    struct AffNode_s *r;

    if (b->ptr == 0)
        return luaL_error(L, "closed writer");

    qL = L;
    r = p[0] == '/'? NULL: b->dir;
    if (aff_writer_mkpath(b->ptr, r, p) == 0)
        return luaL_error(L, aff_writer_errstr(b->ptr));

    return 0;
}

static int
q_aff_reader(lua_State *L)
{
    const char *name = luaL_checkstring(L, 1);
    struct AffReader_s *r;
    const char *msg;

    qL = L;
    r = aff_reader(name);
    msg = aff_reader_errstr(r);
    if (msg == NULL) {
        qlua_newAffReader(L, r);

        return 1;
    } else {
        aff_reader_close(r);

        return luaL_error(L, msg);
    }
}

static int
q_aff_writer(lua_State *L)
{
    const char *name = luaL_checkstring(L, 1);
    struct AffWriter_s *w;
    const char *msg;

    qL = L;
    w = aff_writer(name);
    msg = aff_writer_errstr(w);
    if (msg == NULL) {
        qlua_newAffWriter(L, w);

        return 1;
    } else {
        aff_writer_close(w);

        return luaL_error(L, msg);
    }
}

/* metatables */
static const struct luaL_Reg mtReader[] = {
    { "__tostring",       qaff_r_fmt},
    { "__gc",             qaff_r_gc},
    { "close",            qaff_r_close},
    { "read",             qaff_r_read},
    { "chpath",           qaff_r_chpath },
    { "status",           qaff_r_status },
    { "list",             qaff_r_list },
    { NULL,               NULL}
};

static const struct luaL_Reg mtWriter[] = {
    { "__tostring",       qaff_w_fmt},
    { "__gc",             qaff_w_gc},
    { "close",            qaff_w_close},
    { "write",            qaff_w_write},
    { "chpath",           qaff_w_chpath },
    { "mkpath",           qaff_w_mkpath },
    { "status",           qaff_w_status },
    { NULL,               NULL}
};

/* names and routines for qcd.qdpc table */
static const struct {
    char *name;
    int (*func)(lua_State *L);
} fAFFio[] = {
    { "Reader",   q_aff_reader},
    { "Writer",   q_aff_writer},
    { NULL,       NULL }
};

int
init_aff_io(lua_State *L)
{
    int i;

    lua_getglobal(L, qcdlib);
    lua_newtable(L);
    for (i = 0; fAFFio[i].name; i++) {
        lua_pushcfunction(L, fAFFio[i].func);
        lua_setfield(L, -2, fAFFio[i].name);
    }
    lua_setfield(L, -2, aff_io);
    qlua_metatable(L, mtnReader, mtReader);
    qlua_metatable(L, mtnWriter, mtWriter);
    
    return 0;
}

int
fini_aff_io(lua_State *L)
{
    return 0;
}
