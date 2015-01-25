#include "modules.h"                                                 /* DEPS */
#include "qlua.h"                                                    /* DEPS */
#include "lattice.h"                                                 /* DEPS */
#include "ddpairs_io.h"                                              /* DEPS */
#include "qio_utils.h"                                               /* DEPS */
#include "latdirprop.h"                                              /* DEPS */
#include "qio.h"
#include <string.h>

static const char ddpairs_io[] = "ddpairs";

#if USE_Nc2
#define Qs(a)    a ## 2
#define Qx(a,b)  a ## 2 ## b
#define QC(x)    2
#define QNc     '2'
#define Qcolors  "2"
#include "ddpairs_io-x.c"                                            /* DEPS */
#endif

#if USE_Nc3
#define Qs(a)    a ## 3
#define Qx(a,b)  a ## 3 ## b
#define QC(x)    3
#define QNc     '3'
#define Qcolors  "3"
#include "ddpairs_io-x.c"                                            /* DEPS */
#endif

#if USE_NcN
#define Qs(a)    a ## N
#define Qx(a,b)  a ## N ## b
#define QC(x)    (x)->nc
#define QNc     'N'
#define Qcolors  "N"
#include "ddpairs_io-x.c"                                            /* DEPS */
#endif

/* qcd.ddpairs.read(Lattice,       -- 1, mLattice
 *                  fname)         -- 2, file name string
 */
static int
ddpairs_read(lua_State *L)
{
    mLattice *S = qlua_checkLattice(L, 1);

    switch (S->nc) {
#if USE_Nc2
    case 2:   return dd_read2(L, 0, 2);
#endif
#if USE_Nc3
    case 3:   return dd_read3(L, 0, 3);
#endif
#if USE_NcN
    default:  return dd_readN(L, 0, S->nc);
#else
    default: return luaL_error(L, "bad number of colors");
#endif
    }
}

/* qcd.ddpairs.readN(nc,           -- 1, int number of colors
 *                  Lattice,       -- 2, mLattice
 *                  fname)         -- 3, file name string
 */
static int
ddpairs_readN(lua_State *L)
{
    int nc = luaL_checkint(L, 1);

    if (nc < 1)
        return luaL_error(L, "bad number of colors in qcd.ddpairs.readN()");
    switch (nc) {
#if USE_Nc2
    case 2:   return dd_read2(L, 1, 2);
#endif
#if USE_Nc3
    case 3:   return dd_read3(L, 1, 3);
#endif
#if USE_NcN
    default:  return dd_readN(L, 1, nc);
#else
    default: return luaL_error(L, "bad number of colors");
#endif
    }
}

/* qcd.ddpairs.write(precision,    -- 1, "F" or "D", output size float/double
 *                   fname,        -- 2, file name string
 *                   file_into,    -- 3, string (xml)
 *                   Source,       -- 4, DiracPropagator (source)
 *                   source_info,  -- 5 string (xml)
 *                   time_slice,   -- 6 int >= 0, < Lattice[#Lattice - 1]
 *                   Prop,         -- 7, DiracPropagator (solution)
 *                   prop_info,    -- 8, string (xml)
 *                  [fmt])         -- 9, optional "single" or "multi"
 */
static int
ddpairs_write(lua_State *L)
{
    switch (qlua_qtype(L, 4)) { /* get nc from Source */
#if USE_Nc2
    case qLatDirProp2: return dd_write2(L, 2);
#endif
#if USE_Nc3
    case qLatDirProp3: return dd_write3(L, 3);
#endif
#if USE_NcN
    case qLatDirPropN: {
        mLatDirPropN *src = qlua_checkLatDirPropN(L, 5, NULL, -1);
        return dd_writeN(L, src->nc);
    }
#endif
    default:
        return luaL_error(L, "bad arguments for qcd.ddpairs.write()");
    }
}


static struct luaL_Reg fDDPAIRSio[] = {
    { "read",     ddpairs_read },
    { "readN",    ddpairs_readN },
    { "write",    ddpairs_write },
    { NULL,       NULL }
};

int
init_ddpairs_io(lua_State *L)
{
    lua_getglobal(L, qcdlib);
    lua_newtable(L);
    luaL_register(L, NULL, fDDPAIRSio);
    lua_setfield(L, -2, ddpairs_io);
    lua_pop(L, 1);
    return 0;
}

void
fini_ddpairs_io(void)
{
}
