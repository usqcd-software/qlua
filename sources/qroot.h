#ifndef MARK_968348DF_41CD_4236_BF8E_38B4F285E1F0
#define MARK_968348DF_41CD_4236_BF8E_38B4F285E1F0

#define GSL_RANGE_CHECK_OFF
#define HAVE_INLINE
#include "gsl/gsl_errno.h"
#include "gsl/gsl_multiroots.h"

int init_root(lua_State *L);
int fini_root(lua_State *L);

#endif /* defined(MARK_968348DF_41CD_4236_BF8E_38B4F285E1F0) */
