#include <qlua.h>

const char *progname = "qlua";
const char *qcdlib = "qcd";

void
message(const char *fmt, ...)
{
    if (QDP_this_node == 0) {
        va_list va;

        va_start(va, fmt);
        vfprintf(stderr, fmt, va);
        va_end(va);
    }
}

void
report(lua_State *L, const char *fname, int status)
{
    if (QDP_this_node == 0) {
        if (status && !lua_isnil(L, -1)) {
            const char *msg = lua_tostring(L, -1);
            if (msg == NULL) msg = "(error object is not a string)";
            message("%s ERROR:: %s\n", progname, msg);
            lua_pop(L, 1);
        }
    }
}

int
main(int argc, char *argv[])
{
    int status = 1;
    int i;
    lua_State *L = NULL;

    if (QDP_initialize(&argc, &argv)) {
        fprintf(stderr, "QDP initialization failed\n");
        return 1;
    }
    L = lua_open();
    if (L == NULL) {
        message("can not create Lua state");
        goto end;
    }
    qlua_init(L);  /* open libraries */

    for (i = 1; i < argc; i++) {
        status = luaL_dofile(L, argv[i]);
        report(L, argv[i], status);
        if (status)
            break;
    }
    lua_close(L);
end:
    QDP_finalize();
    return status;
}
