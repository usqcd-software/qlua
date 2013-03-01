#ifndef MARK_5E85078A_F7B5_4BB6_A219_FA360F86D805
#define MARK_5E85078A_F7B5_4BB6_A219_FA360F86D805
#define GSL_RANGE_CHECK_OFF
#define HAVE_INLINE
#include "gsl/gsl_errno.h"
#include "gsl/gsl_multimin.h"

double vnorm2_gsl(gsl_vector *v, int dim);
int init_min(lua_State *L);
int fini_min(lua_State *L);

#endif /* defined(MARK_5E85078A_F7B5_4BB6_A219_FA360F86D805) */
