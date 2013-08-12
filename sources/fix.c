#include "modules.h"                                                 /* DEPS */
#include "qlua.h"                                                    /* DEPS */
#include "fix.h"                                                     /* DEPS */
#include "qmp.h"
#include <string.h>
#include <sys/time.h>
#include <sys/stat.h>
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

    if (QDP_this_node == qlua_master_node) {
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
    }
    XMP_dist_int_array(qlua_master_node, 1, &status);

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
    const char *n = luaL_checkstring(L, 1);
    const char *m = luaL_checkstring(L, 2);

    if (!strcmp(m, "r") || !strcmp(m, "rb") || !strcmp(m, "rt")) {
      double v;
      /* read files are open everywhere */
      f->file = fopen(n, m);
      f->kind = qf_other;
      v = (f->file != 0);
      /* make sure all nodes manage to open the file */
      QMP_min_double(&v);
      if (v < 1)
        if (v == 0)
          return luaL_error(L, "file open failed");
    } else {
      /* other modes are write, open the file on the master node only */
      int v;
      if (QDP_this_node == qlua_master_node) {
        f->file = fopen(n, m);
        f->kind = qf_other;
        v = (f->file != 0);
      } else {
        f->file = NULL;
        f->kind = qf_dummy;
      }
      XMP_dist_int_array(qlua_master_node, 1, &v);
      if (v == 0)
        return luaL_error(L, "file open failed");
    }

    return 1;
}

static int
q_fexists(lua_State *L)
{
  const char *fname = luaL_checkstring(L, 1);
  struct stat st;

  if ((stat(fname, &st) == 0) &&
      S_ISREG(st.st_mode))
    lua_pushboolean(L, 1);
  else
    lua_pushboolean(L, 0);

  return 1;
}

static int
q_dexists(lua_State *L)
{
  const char *fname = luaL_checkstring(L, 1);
  struct stat st;

  if ((stat(fname, &st) == 0) &&
      S_ISDIR(st.st_mode))
    lua_pushboolean(L, 1);
  else
    lua_pushboolean(L, 0);

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
    
    if (QDP_this_node == qlua_master_node) {
        gettimeofday(&t, NULL);
        v = t.tv_sec + 1e-6 * t.tv_usec;
    }
    XMP_dist_double_array(qlua_master_node, 1, &v);
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
    { "__newindex", qlua_nowrite },
    { "close",      qf_close },
    { "write",      qf_write },
    { "lines",      qf_lines },
    { "flush",      qf_flush },
    { NULL,         NULL}
};

static int
qlua_getmetatable(lua_State *L)
{
    if (lua_type(L, 1) != LUA_TUSERDATA)
        return lua_getmetatable(L, 1);
#if 0 /* ??? should work, but does not ... */
    lua_createtable(L, 0, 0); /* proxy */
    lua_createtable(L, 0, 3); /* proxy metatable */
    lua_pushcfunction(L, qlua_nowrite);
    lua_setfield(L, -2, "__newindex");
    lua_pushvalue(L, -3);
    lua_setfield(L, -2, "__index");
    lua_pushvalue(L, -1);
    lua_setfield(L, -2, "__metatable");
    lua_pushvalue(L, -1);
    lua_setmetatable(L, -2);
    lua_setmetatable(L, -2);
#endif
    return 1;
}

/* This should be a simple table, but for now this will do */
static char *
q_qtypename(lua_State *L, int idx, char *def)
{
    char *t = def;
    QLUA_Type tp = qlua_qtype(L, idx);

    switch (tp) {
    case qComplex: t = "complex"; break;
    case qGamma: t = "gamma"; break;
    case qMatReal: t = "matrix.real"; break;
    case qMatComplex: t = "matrix.complex"; break;
    case qSeqColVec2: case qSeqColVec3: case qSeqColVecN: t = "color.vector"; break;
    case qSeqColMat2: case qSeqColMat3: case qSeqColMatN: t = "color.matrix"; break;
    case qSeqDirFerm2: case qSeqDirFerm3: case qSeqDirFermN: t = "dirac.fermion"; break;
    case qSeqDirProp2: case qSeqDirProp3: case qSeqDirPropN: t = "dirac.propagator"; break;
    case qLatInt: t = "lattice.int"; break;
    case qLatReal: t = "lattice.real"; break;
    case qLatComplex: t = "lattice.complex"; break;
    case qLatColVec2: case qLatColVec3: case qLatColVecN: t = "lattice.color.vector"; break;
    case qLatColMat2: case qLatColMat3: case qLatColMatN: t = "lattice.color.matrix"; break;
    case qLatDirFerm2: case qLatDirFerm3: case qLatDirFermN: t = "lattice.dirac.fermion"; break;
    case qLatDirProp2: case qLatDirProp3: case qLatDirPropN: t = "lattice.dirac.propagator"; break;
    case qLattice: t = "lattice"; break;
    case qLatMulti: t = "multiset"; break;
    case qLatSubset: t = "subset"; break;
    case qVecInt: t = "vector.int"; break;
    case qVecReal: t = "vector.real"; break;
    case qVecComplex: t = "vector.complex"; break;
    case qSeqRandom: t = "random.state"; break;
    case qLatRandom: t = "lattice.random.state"; break;
    case qReader: t = "qio.reader"; break;
    case qWriter: t = "qio.writer"; break;
    case qAffReader: t = "aff.reader"; break;
    case qAffWriter: t = "aff.writer"; break;
    case qClover: t = "clover"; break;
    case qCloverDeflator: t = "clover.deflator"; break;
    case qCloverDeflatorState: t = "clover.deflator.state"; break;
    case qMDWF: t = "mdwf"; break;
    case qMDWFDeflator: t = "mdwf.deflator"; break;
    case qMDWFDeflatorState: t = "mdwf.deflator.state"; break;
    default: break;
    }
    return t;
}

static int
q_type(lua_State *L)
{
    const char *t = "userdata";
    luaL_checkany(L, 1);
    switch (lua_type(L, 1)) {
    case LUA_TNIL: t = "nil"; break;
    case LUA_TBOOLEAN: t = "boolean"; break;
    case LUA_TLIGHTUSERDATA: t = "userdata"; break;
    case LUA_TNUMBER: t = "number"; break;
    case LUA_TSTRING: t = "string"; break;
    case LUA_TTABLE: t = "table"; break;
    case LUA_TFUNCTION: t = "function"; break;
    case LUA_TUSERDATA: t = q_qtypename(L, 1, "userdata"); break;
    case LUA_TTHREAD: t = "thread"; break;
    }
    lua_pushstring(L, t);
    return 1;
}

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

    lua_pushcfunction(L, q_type);
    lua_setglobal(L, "type");
    
    lua_pushcfunction(L, qlua_print);
    lua_setglobal(L, "print");

    qlua_metatable(L, mtnFile, mtFile, qNoType);
    lua_createtable(L, 0, 3);
    qlua_newFile(L)->kind = qf_stdout;
    lua_setfield(L, -2, "stdout");
    qlua_newFile(L)->kind = qf_stderr;
    lua_setfield(L, -2, "stderr");
    lua_pushcfunction(L, q_file);
    lua_setfield(L, -2, "open");
    lua_pushcfunction(L, q_fexists);
    lua_setfield(L, -2, "fexists");
    lua_pushcfunction(L, q_dexists);
    lua_setfield(L, -2, "dexists");
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

    /* fix package.path -- try to get QLUALIB from the environment first */
        {
                char *envlib = getenv("QLUALIB");

                lua_getglobal(L, "package");
                if (envlib != NULL) {
                        lua_pushstring(L, envlib);
                } else {
                        lua_pushstring(L, qlib_path);
                }
                lua_setfield(L, -2, "path");
                lua_pop(L, 1);
        }

    /* fix getmetatable */
    lua_pushcfunction(L, qlua_getmetatable);
    lua_setglobal(L, "getmetatable");

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
