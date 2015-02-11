#ifndef MARK_91468699_C123_4148_86B5_3DAB4308E77E
#define MARK_91468699_C123_4148_86B5_3DAB4308E77E
#include "qla.h"

extern const char mtnComplex[];

int init_complex(lua_State *L);
void fini_complex(void);

QLA_D_Complex *qlua_checkComplex(lua_State *L, int idx);
QLA_D_Complex *qlua_newComplex(lua_State *L);

int q_c_gaussian(lua_State *L);

#endif /* !defined(MARK_91468699_C123_4148_86B5_3DAB4308E77E) */
