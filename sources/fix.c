#include <qlua.h>                                                    /* DEPS */
#include <fix.h>                                                     /* DEPS */
#include <qmp.h>
#include <string.h>

static char self[72];

static const char mtnFile[] = "qlua.file";

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
            return luaL_error(L, luaL_checkstring(L, -1));
        str = luaL_checkstring(L, -1);
        luaL_addstring(&b, str);
        lua_pop(L, 1);
    }
    luaL_addstring(&b, "\n");
    luaL_pushresult(&b);
    str = luaL_checkstring(L, -1);
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

static struct luaL_Reg mtFile[] = {
    { "__tostring", qf_fmt },
    { "__gc",       qf_gc },
    { "close",      qf_close },
    { "write",      qf_write },
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
    lua_setglobal(L, "io");

    /* fix os.exit() */
    lua_getglobal(L, "os");
    lua_pushcfunction(L, qlua_exit);
    lua_setfield(L, -2, "exit");
    
    return 0;
}

int
fini_qlua_io(lua_State *L)
{
    return 0;
}