#include "modules.h"                                                 /* DEPS */
#include "qlua.h"                                                    /* DEPS */
#include "lattice.h"                                                 /* DEPS */
#include "qquda.h"                                                   /* DEPS */
// #include "qqudaxx.hpp"                                               /* DEPS */
#include "quda_milc_interface.h"
#include <string.h>

/* Interface:
 *
 ** Create QUDA context for lattice L (can only be called once!)
 *   qobj = qcd.quda {
 *    lattice = L,              -- qlua lattice object
 *    device = 0,               -- (default 0) CUDA device number
 *    verbosity = "summary"     -- (default "silent"): "silent", "verbose", "debug", "summary"
 *   };
 *
 ** 
 */

static int in_quda_p = 0;

/* qcd.quda function: extract lattice.size and lattice.network, call qudaInit(), remember the result */

static int
create_quda(lua_State *L)
{
  if (in_quda_p == 0) {
    qlua_checktable(L, 1, "qcd.quda()");
    int dev = qlua_tabkey_intopt(L, 1, "device", 0);
    const char *verbose_name = qlua_tabkey_stringopt(L, 1, "verbosity", "silent");
    int Sidx = qlua_tabkey(L, 1, "lattice");
    mLattice *S = qlua_checkLattice(L, Sidx);
    QudaInitArgs_t args;
    
    if (strcmp(verbose_name, "silent") == 0)
    args.verbosity = QUDA_SILENT;
    else if (strcmp(verbose_name, "verbose") == 0)
      args.verbosity = QUDA_VERBOSE;
    else if (strcmp(verbose_name, "debug") == 0)
      args.verbosity = QUDA_DEBUG_VERBOSE;
    else if (strcmp(verbose_name, "summary") == 0)
      args.verbosity = QUDA_SUMMARIZE;
    else
      luaL_error(L, "Unexpected verbosity: %s", verbose_name);
    
    if (S->rank != 4)
      luaL_error(L, "Expected 4-d lattice: lattice rank = %d", S->rank);
    
    args.layout.device = dev;
    args.layout.latsize = S->dim;
    args.layout.machsize = S->net;
    
    qudaInit(args);
    in_quda_p = 1;
    
    /* XXX -- create the object */
    return 0;
  }
  luaL_error(L, "qcd.quda() called second time");
  return 0;
}

static struct luaL_Reg fquda[] = {
  { "quda",   create_quda },
  { NULL,     NULL }
};

int
init_quda(lua_State *L)
{
  luaL_register(L, qcdlib, fquda);
  return 0;
}

void
fini_quda(void)
{
  if (in_quda_p)
    qudaFinalize();
  in_quda_p = 0;
}
