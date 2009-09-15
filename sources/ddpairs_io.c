#include <qlua.h>                                                    /* DEPS */
#include <ddpairs_io.h>                                              /* DEPS */
#include <qio_utils.h>                                               /* DEPS */
#include <qio.h>

static const char ddpairs_io[] = "ddpairs";

static int
ddpairs_read(lua_State *L)
{
    return luaL_error(L, "XXX qcd.ddpairs.read() not implemented");
}

static int
ddpairs_write(lua_State *L)
{
    return luaL_error(L, "XXX qcd.ddpairs.write() not implemented");
}

static struct luaL_Reg fDDPAIRSio[] = {
    { "read",     ddpairs_read },
    { "write",    ddpairs_write },
    { NULL,       NULL }
};

int
init_ddpairs_io(lua_State *L)
{
    int i;

    lua_getglobal(L, qcdlib);
    lua_newtable(L);
    for (i = 0; fDDPAIRSio[i].name; i++) {
        lua_pushcfunction(L, fDDPAIRSio[i].func);
        lua_setfield(L, -2, fDDPAIRSio[i].name);
    }
    lua_setfield(L, -2, ddpairs_io);
    lua_pop(L, 1);
    return 0;
}

int
fini_ddpairs_io(lua_State *L)
{
    return 0;
}
