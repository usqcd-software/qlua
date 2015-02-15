#include "modules.h"                                                 /* DEPS */
#include "qlua.h"                                                    /* DEPS */
#include "lattice.h"                                                 /* DEPS */
#include "qdpcc_io.h"                                                /* DEPS */
#include "qio_utils.h"                                               /* DEPS */
#include "latdirferm.h"                                              /* DEPS */
#include "latdirprop.h"                                              /* DEPS */
#include "qio.h"

static const char qdpcc_io[] = "qdpcc";

#if USE_Nc2
#define Qs(a)    a ## 2
#define Qx(a,b)  a ## 2 ## b
#define QC(x)    2
#define QNc     '2'
#define Qcolors  "2"
#include "qdpcc_io-x.c"                                              /* DEPS */
#endif

#if USE_Nc3
#define Qs(a)    a ## 3
#define Qx(a,b)  a ## 3 ## b
#define QC(x)    3
#define QNc     '3'
#define Qcolors  "3"
#include "qdpcc_io-x.c"                                              /* DEPS */
#endif

#if USE_NcN
#define Qs(a)    a ## N
#define Qx(a,b)  a ## N ## b
#define QC(x)    (x)->nc
#define QNc     'N'
#define Qcolors  "N"
#include "qdpcc_io-x.c"                                              /* DEPS */
#endif


/* qcd.qdpcc.read_prop(Lattice,    -- 1 mLattice
 *                     file_name)  -- 2 file name string
 */
static int
qdpcc_read_P(lua_State *L)
{
    mLattice *S = qlua_checkLattice(L, 1);

    switch (S->nc) {
#if USE_Nc2
    case 2:  return cc_read_P_2(L, 0, 2);
#endif
#if USE_Nc3
    case 3:  return cc_read_P_3(L, 0, 3);
#endif
#if USE_NcN
    default: return cc_read_P_N(L, 0, S->nc);
#else
    default: return luaL_error(L, "bad number of colors");
#endif
    }
}

/* qcd.qdpcc.read_propN(nc,         -- 1 number of colors
 *                      Lattice,    -- 2 mLattice
 *                      file_name)  -- 3 file name string
 */
static int
qdpcc_read_P_N(lua_State *L)
{
    int nc = luaL_checkint(L, 1);

    if (nc < 1)
        return luaL_error(L, "bad number of colors for qcd.qdpcc.read_propN()");

    switch (nc) {
#if USE_Nc2
    case 2:  return cc_read_P_2(L, 1, 2);
#endif
#if USE_Nc3
    case 3:  return cc_read_P_3(L, 1, 3);
#endif
#if USE_NcN
    default: return cc_read_P_N(L, 1, nc);
#else
    default: return luaL_error(L, "bad number of colors");
#endif
    }
}

/* qcd.qdpcc.write_prop(fmt,        -- 1 "F" for float, "D" for double
 *                      file_name,  -- 2 file name string
 *                      file_info,  -- 3 file info string (xml)
 *                      prop,       -- 4 DiracPropagator to write
 *                      prop_info,  -- 5 record info string (xml)
 *                      [volfmt])   -- 6 optional "single" or "multi"
 */
static int
qdpcc_write_P(lua_State *L)
{
    switch (qlua_qtype(L, 4)) {
#if USE_Nc2
    case qLatDirProp2:  return cc_write_P_2(L, 2);
#endif
#if USE_Nc3
    case qLatDirProp3:  return cc_write_P_3(L, 3);
#endif
#if USE_NcN
    case qLatDirPropN: {
        mLatDirPropN *prop = qlua_checkLatDirPropN(L, 4, NULL, -1);
        return cc_write_P_N(L, prop->nc);
    }
#endif
    default:
        return luaL_error(L, "bad arguments for qcd.qdpcc.write_prop()");
    }
}

static struct luaL_Reg fQDPCCio[] = {
    { "read_propN",    qdpcc_read_P_N },
    { "read_prop",     qdpcc_read_P   },
    { "write_prop",    qdpcc_write_P  },
    { NULL,       NULL }
};

int
init_qdpcc_io(lua_State *L)
{
    lua_getglobal(L, qcdlib);
    lua_newtable(L);
    luaL_register(L, NULL, fQDPCCio);
    lua_setfield(L, -2, qdpcc_io);
    lua_pop(L, 1);
    return 0;
}

void
fini_qdpcc_io(void)
{
}
