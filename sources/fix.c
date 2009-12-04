#include <qlua.h>                                                    /* DEPS */
#include <fix.h>                                                     /* DEPS */
#include <modules.h>
#include <qmp.h>
#include <string.h>
#include <sys/time.h>
#include <stdio.h>

static char self[72];

static const char mtnFile[] = "qlua.file";

static char qlib_path[] = "./?.qlua;" QLUA_LIB "/?.qlua;./qlib/?.qlua";

static FILE *rf = 0;

enum {
    qf_closed,
    qf_stdout,
    qf_stderr,
    qf_dummy,
    qf_other
} qFKind;

typedef struct {
    int   kind;
    FILE *file;
} mFile;

static mFile *
qlua_newFile(lua_State *L)
{
    mFile *v = lua_newuserdata(L, sizeof (mFile));

    v->kind = qf_closed;
    v->file = 0;
    luaL_getmetatable(L, mtnFile);
    lua_setmetatable(L, -2);

    return v;
}

static mFile *
qlua_checkFile(lua_State *L, int idx)
{
    void *v = luaL_checkudata(L, idx, mtnFile);
    
    luaL_argcheck(L, v != 0, idx, "File expected");

    return v;
}

static int
qf_fmt(lua_State *L)
{
    mFile *v = qlua_checkFile(L, 1);
    
    switch (v->kind) {
    case qf_closed: lua_pushstring(L, "File[closed]"); break;
    case qf_stdout: lua_pushstring(L, "File[stdout]"); break;
    case qf_stderr: lua_pushstring(L, "File[stderr]"); break;
    case qf_other:
    case qf_dummy: {
        char fmt[72];
        sprintf(fmt, "File[%p]", v->file);
        lua_pushstring(L, fmt);
        break;
    }
    default:
        return luaL_error(L, "bad file kind %d", v->kind);
    }
    return 1;
}

static int
qf_gc(lua_State *L)
{
    mFile *v = qlua_checkFile(L, 1);

    if (v->kind == qf_other)
        fclose(v->file);
    v->file = 0;
    v->kind = qf_closed;

    return 0;
}

static int
qf_close(lua_State *L)
{
    mFile *v = qlua_checkFile(L, 1);

    switch (v->kind) {
    case qf_stdout: fflush(stdout); break;
    case qf_stderr: fflush(stderr); break;
    case qf_other:
        fclose(v->file);
        /* through */
    case qf_dummy:
        v->file = 0;
        v->kind = qf_closed;
        break;
    default:
        return luaL_error(L, "bad file kind in close (%d)", v->kind);
    }

    return 0;
}

static int
qf_flush(lua_State *L)
{
    mFile *v = qlua_checkFile(L, 1);
    
    switch (v->kind) {
    case qf_stdout: fflush(stdout); break;
    case qf_stderr: fflush(stderr); break;
    case qf_dummy: break;
    case qf_other: fflush(v->file); break;
    default:
        return luaL_error(L, "bad file kind in flush (%d)", v->kind);
    }

    return 0;
}

static int
qf_write(lua_State *L)
{
    int status;

    if (qlua_primary_node) {
        mFile *v = qlua_checkFile(L, 1);
        const char *s= luaL_checkstring(L, 2);
        int n = strlen(s);
        FILE *f;
        
        switch (v->kind) {
        case qf_stdout: f = stdout; break;
        case qf_stderr: f = stderr; break;
        case qf_other: f = v->file; break;
        default:
            return luaL_error(L, "bad kind of file in write (%d)", v->kind);
        }
        status = fwrite(s, n, 1, f) == 1;
    } else {
        status = 0;
    }
    QMP_sum_int(&status);

    if (status == 0)
        return luaL_error(L, "write failed");

    return 0;
}

static int
read_line(lua_State *L, FILE *f)
{
    luaL_Buffer b;

    luaL_buffinit(L, &b);
    for (;;) {
        size_t l;
        char *p = luaL_prepbuffer(&b);
        if (fgets(p, LUAL_BUFFERSIZE, f) == NULL) {  /* eof? */
            luaL_pushresult(&b);  /* close buffer */
            return (lua_objlen(L, -1) > 0);  /* check whether read something */
        }
        l = strlen(p);
        if (l == 0 || p[l-1] != '\n')
            luaL_addsize(&b, l);
        else {
            luaL_addsize(&b, l - 1);  /* do not include `eol' */
            luaL_pushresult(&b);  /* close buffer */
            return 1;  /* read at least an `eol' */
        }
    }
}

static int
q_readline(lua_State *L)
{
    mFile *v = (mFile *)lua_touserdata(L, lua_upvalueindex(1));
    int success;

    if (v->file == 0)
        luaL_error(L, "file is already closed");
    success = read_line(L, v->file);
    if (ferror(v->file))
        return luaL_error(L, "line reading error");
    if (success)
        return 1;
    if (lua_toboolean(L, lua_upvalueindex(2))) {
        lua_settop(L, 0);
        lua_pushvalue(L, lua_upvalueindex(1));
        qf_close(L);
    }
    return 0;
}

static int
qf_lines(lua_State *L)
{
    qlua_checkFile(L, 1);
    lua_pushboolean(L, 0);
    lua_pushcclosure(L, q_readline, 2);
    return 1;
}

static int
q_lines(lua_State *L)
{
    const char *n = luaL_checkstring(L, 1);
    mFile *v = qlua_newFile(L);

    v->file = fopen(n, "rt");
    v->kind = qf_other;
    if (v->file == 0)
        return luaL_error(L, "file open failed");

    lua_pushboolean(L, 1);
    lua_pushcclosure(L, q_readline, 2);
    return 1;
}

static int
q_file(lua_State *L)
{
    mFile *f = qlua_newFile(L);
    int v;

    if (qlua_primary_node) {
        const char *n = luaL_checkstring(L, 1);
        const char *m = luaL_checkstring(L, 2);
        f->file = fopen(n, m);
        f->kind = qf_other;

        v = (f->file != 0);
    } else {
        f->file = (FILE *)1;
        f->kind = qf_dummy;
        v = 0;
    }
    QMP_sum_int(&v);

    if (v == 0)
        return luaL_error(L, "file open failed");

    return 1;
}

static int
qlua_print(lua_State *L)
{
    int n = lua_gettop(L);
    luaL_Buffer b;
    int i;
    const char *str;

    luaL_buffinit(L, &b);
    for (i = 0, str = self; i < n; i++, str = "\t") {
        luaL_addstring(&b, str);
        lua_getglobal(L, "tostring");
        lua_pushvalue(L, i + 1);
        if (lua_pcall(L, 1, 1, 0))
            return luaL_error(L, qlua_checkstring(L, -1,
                                                  "no #%d:tostring()", i + 1));
        str = qlua_checkstring(L, -1, "#%d:tostring() is no string", i + 1);
        luaL_addstring(&b, str);
        lua_pop(L, 1);
    }
    luaL_addstring(&b, "\n");
    luaL_pushresult(&b);
    str = lua_tostring(L, -1);
    printf("%s", str);
    
    return 0;
}

static int
qlua_exit(lua_State *L)
{
    int c = lua_gettop(L) == 0? 0: luaL_checknumber(L, 1);
    
    QDP_finalize();
    exit(c);
    /* never happens */
    return 0;
}

static int
qlua_timeofday(lua_State *L)
{
    struct timeval t;
    double v;
    
    if (qlua_primary_node) {
        gettimeofday(&t, NULL);
        v = t.tv_sec + 1e-6 * t.tv_usec;
    } else {
        v = 0;
    }
    QMP_sum_double(&v);
    lua_pushnumber(L, v);

    return 1;
}

static int
qlua_random(lua_State *L)
{
    uint32_t v;

    if (fread(&v, sizeof (v), 1, rf) != 1)
        return luaL_error(L, "error reading random source");
    lua_pushnumber(L, v);

    return 1;
}

static struct luaL_Reg mtFile[] = {
    { "__tostring", qf_fmt },
    { "__gc",       qf_gc },
    { "close",      qf_close },
    { "write",      qf_write },
    { "lines",      qf_lines },
    { "flush",      qf_flush },
    { NULL,         NULL}
};

int
init_qlua_io(lua_State *L)
{
    int n = QMP_get_number_of_nodes();
    int l;

    if (n > 1) {
        sprintf(self, "[%d]:", n);
        l = strlen(self);
        sprintf(self, "[%0*d]:", l - 3, QDP_this_node);
    } else {
        self[0] = 0;
    }

    lua_pushcfunction(L, qlua_print);
    lua_setglobal(L, "print");

    qlua_metatable(L, mtnFile, mtFile);
    lua_createtable(L, 0, 3);
    qlua_newFile(L)->kind = qf_stdout;
    lua_setfield(L, -2, "stdout");
    qlua_newFile(L)->kind = qf_stderr;
    lua_setfield(L, -2, "stderr");
    lua_pushcfunction(L, q_file);
    lua_setfield(L, -2, "open");
    lua_pushcfunction(L, q_lines);
    lua_setfield(L, -2, "lines");
    lua_setglobal(L, "io");

    /* fix os.exit() */
    lua_getglobal(L, "os");
    lua_pushcfunction(L, qlua_exit);
    lua_setfield(L, -2, "exit");
    lua_pushcfunction(L, qlua_timeofday);
    lua_setfield(L, -2, "time");
    rf = fopen("/dev/urandom", "rb");
    if (rf) {
        lua_pushcfunction(L, qlua_random);
        lua_setfield(L, -2, "random");
    }
    
    lua_pop(L, 1);

    /* fix package.path */
    lua_getglobal(L, "package");
    lua_pushstring(L, qlib_path);
    lua_setfield(L, -2, "path");
    lua_pop(L, 1);

    return 0;
}

int
fini_qlua_io(lua_State *L)
{
    if (rf) {
        fclose(rf);
        rf = 0;
    }
    return 0;
}
